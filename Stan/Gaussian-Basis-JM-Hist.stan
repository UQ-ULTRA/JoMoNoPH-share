// =============================================================================
// JOINT LONGITUDINAL–SURVIVAL MODEL
// BASELINE HAZARD: GAUSSIAN BASIS (GB)
// =============================================================================
//
// OVERVIEW:
// This Stan program implements the same joint longitudinal–survival
// model structure as the BP version, but uses a Gaussian basis expansion
// for the baseline hazard instead of Bernstein polynomials.
//
// The GB model is included to assess sensitivity to the choice of
// baseline hazard representation.
//
// HOW THIS FILE IS USED:
//   - This model is called from the same R pipeline as the BP model.
//   - Selection is controlled by joint_model_type = "GB".
//   - All data inputs, priors, and initial values are constructed in R;
//     this file assumes they are already correctly formatted.
//
// BASELINE HAZARD (KEY IDEA):
//   - The log-baseline hazard is represented as a weighted sum of
//     Gaussian basis functions evaluated over time.
//   - Basis locations and scales are fixed upstream; only weights
//     are estimated here.
//   - Smoothness is controlled through prior regularization rather
//     than explicit penalties.
//
// JOINT MODEL LINKAGE:
//   - Survival depends on the *current longitudinal value* via
//     alpha_tilde, identical to the BP formulation.
//   - This allows direct comparison of baseline hazard flexibility
//     without changing the longitudinal or association structure.
//
// IMPORTANT DESIGN NOTES:
//   - This model shares parameter naming and structure with the
//     Bernstein polynomial model to allow like-for-like comparison.
//   - Differences in performance should be interpreted as driven
//     by baseline hazard representation, not other model components.
//
// WHAT A READER SHOULD NOT WORRY ABOUT:
//   - Knot placement or scaling: handled upstream.
//   - Whether BP or GB is “correct”: both are exploratory options
//     in a simulation study.
//   - PH vs AFT interpretation: this is an AFT-style joint model.
//
// =============================================================================

functions {
  // AFT joint survival log likelihood
  vector loglik_aft_jm_gauss(
    vector time,
    vector beta_surv,
    vector beta_long,
    matrix b_long,
    vector gamma,
    vector status,
    matrix X_surv,
    real alpha,
    vector Y_long_surv,
    vector mu,            // FIXED-BASIS: knots passed in
    vector sigma          // FIXED-BASIS: widths passed in
  )
  {
    int n = num_elements(status);
    int m = num_elements(gamma);
    vector[n] log_lik;
    vector[n] h0;
    vector[n] H0;
    vector[n] y;

    // current-value linkage (unchanged)
    vector[n] C1 = X_surv * beta_surv + alpha * (beta_long[1] + b_long[,1]);
    vector[n] C2 = alpha * (beta_long[2] + X_surv[,1] * beta_long[3] + b_long[,2]);

    for (i in 1:n){
      if (C2[i] != 0) {
        y[i] = exp(-C1[i]) * (1 - exp(-C2[i] * time[i])) / C2[i];
      } else {
        y[i] = exp(-C1[i]) * time[i];
      }
    }
    real eps = 1e-6;
    real tau_aft = max(y) + eps;

    // scale to [0,1] (unchanged)
    vector[n] y_alt = fmin(fmax(y ./ tau_aft, eps), 1 - eps);

    // ---- FIXED-BASIS: remove all sorting/quantile logic; use supplied mu/sigma ----

    // build density & cdf matrices (unchanged apart from using fixed mu/sigma)
    matrix[n,m] f_mat;
    matrix[n,m] F_mat;
    for (k in 1:m) {
      real sk = fmax(sigma[k], 1e-6);
      for (i in 1:n) {
        f_mat[i,k] = normal_lpdf(y_alt[i] | mu[k], sk);
        F_mat[i,k] = normal_lcdf(y_alt[i] | mu[k], sk);
      }
    }

    // f_mat = exp(f_mat) ./ tau_aft;  // f_y(y) = f_u(u)/tau
    // F_mat = exp(F_mat);             // F_u(u)
    f_mat = exp(f_mat); // test on 20251002
    F_mat = exp(F_mat).* tau_aft; // test on 20251002

    h0 = f_mat * gamma;
    H0 = F_mat * gamma;

    log_lik = ((log(h0) - (X_surv * beta_surv + alpha * Y_long_surv)) .* status) - H0;
    return log_lik;
  }
}

// Data
data {
  //// survival data sizes
  int<lower=1> K_long;
  int<lower=1> q;
  int<lower=1> m;
  
  //// smoothing parameters
  int<lower=1> r;  // order of difference (1 = first difference, 2 = second difference, etc.)
  // real<lower=0> tau_h;  // smoothing parameter (try adaptive for now)

  // prior locations and scales
  vector[K_long] beta_long_prior_mean;
  vector<lower=0>[K_long] beta_long_prior_scale;
  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;

  //// longitudinal data
  int<lower=1> N_long;
  vector[N_long] Y_long;
  matrix[N_long, K_long] X_long;
  int<lower=1> N_1_long;
  array[N_long] int<lower=1> J_1_long;

  // individual/group-level predictor values (longitudinal)
  vector[N_long] Z_1_1_long;
  vector[N_long] Z_1_2_long;

  //// survival data
  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=0>[n] time;
  matrix[n, q] X_surv;

  //// linking data
  array[n] int<lower=1> J_1_unique;
  matrix[n, K_long] X_long_surv;
  real s_long;  // scale for longitudinal process
  real s_surv;  // scale for survival linear predictor
  real<lower=0> alpha_tilde_sd;  // prior scale for alpha_tilde (clinically calibrated)

}

// ---- FIXED-BASIS: define knots once, independent of parameters ----
transformed data {
  vector[m] mu_fix;
  vector[m] sigma_fix;
  
  mu_fix[1] = 0;
  mu_fix[m] = 1;
  
  sigma_fix = rep_vector(2.0/(3.0*(m-1.0)),m);    
  
  for (k in 1:(m-2)) {
    mu_fix[k+1] = k * 1.0 / (m-1);  // currently equally spaced in (0,1)
  }
  
  // ---- R-TH DIFFERENCE PENALTY MATRIX ----
  
  // Build the r-th difference penalty matrix Delta_r
  int n_diff = m - r;  // number of r-th differences (total degree - order of difference)
  matrix[n_diff, m] Delta_r = rep_matrix(0, n_diff, m);
  
  // Fill the difference matrix
  for (i in 1:n_diff) {
    // Compute binomial coefficients for r-th difference
    for (j in 0:r) {
      int sign = (j % 2 == 0) ? 1 : -1;
      real binom_coef = 1;
      
      // Compute binomial coefficient C(r,j)
      if (j > 0) {
        for (k in 1:j) {
          binom_coef *= (r - k + 1) * 1.0 / k;
        }
      }
      
      Delta_r[i, i + j] = sign * binom_coef;
    }
  }
  
  // Compute Delta_r^T * Delta_r for the quadratic penalty
  matrix[m, m] penalty_matrix = Delta_r' * Delta_r;
  
  // Compute rank of penalty matrix (number of non-zero eigenvalues)
  int rho = n_diff;  // For r-th differences, rank is typically m-r
}

// Parameters
parameters {
  //// longitudinal
  vector[K_long] beta_long;
  real<lower=0> sigma_long;
  vector<lower=1e-3>[2] sd_1_long;
  matrix[2, N_1_long] z_1_long;
  cholesky_factor_corr[2] L_1_long;

  //// survival
  vector[q] beta_surv;
  vector<lower=0>[m] gamma;  // Gaussian-basis weights (arm = 0)
  real<lower=0> tau_h; // smoothing parameter

  //// linking
  real alpha_tilde;
  
}

// Transformed Parameters
transformed parameters {
  matrix[2, N_1_long] b_raw = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';

  // predicted Y_i(t) at event/censoring time (unchanged)
  vector[n] Y_long_surv;
  for (i in 1:n) {
    Y_long_surv[i] = beta_long[1]
                   + beta_long[2] * X_long_surv[i, 2]
                   + beta_long[3] * X_long_surv[i, 3]
                   + b_long[J_1_unique[i], 1]
                   + b_long[J_1_unique[i], 2] * time[i];
  }
  
  real alpha = alpha_tilde * (s_surv / s_long);  // scale according to relative difference between longitudinal and survival
  
  // likelihood with fixed basis
  vector[n] log_lik = loglik_aft_jm_gauss(
    time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha,
    Y_long_surv, mu_fix, sigma_fix   // <-- pass fixed knots
  );
  
  
}

// Model
model {
  //// longitudinal priors
  beta_long ~ normal(beta_long_prior_mean, beta_long_prior_scale);
  sigma_long ~ cauchy(0, 5) T[0, ];
  sd_1_long ~ cauchy(0, 5) T[0, ];
  L_1_long ~ lkj_corr_cholesky(2);
  to_vector(z_1_long) ~ std_normal();

  // longitudinal likelihood
  vector[N_long] mu_long = X_long * beta_long;
  for (i in 1:N_long) {
    mu_long[i] += b_long[J_1_long[i], 1] * Z_1_1_long[i]
                + b_long[J_1_long[i], 2] * Z_1_2_long[i];
  }
  Y_long ~ normal(mu_long, sigma_long);

  //// survival priors
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  alpha_tilde ~ normal(0, alpha_tilde_sd);
  
  // Penalized spline coefficients with r-th difference penalty
  // This implements: p(γ_h0 | τ_h) ∝ τ_h^(ρ/2) * exp(-τ_h/2 * γ_h0^T * Δ_r^T * Δ_r * γ_h0)
  tau_h ~ gamma(0.5, 0.005); // or half-normal(0, 1) on sqrt(tau_h) WILL NEED PPC ON THIS
  target += 0.5 * rho * log(tau_h);  // Normalization constant τ_h^(ρ/2)
  target += -0.5 * tau_h * quad_form(penalty_matrix, gamma);  // Quadratic penalty
  gamma ~ normal(0, 1) T[0, ];  // Keep base prior for positivity constraint

  // survival likelihood
  target += sum(log_lik);
}
