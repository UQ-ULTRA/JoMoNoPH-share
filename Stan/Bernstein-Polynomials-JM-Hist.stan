// =============================================================================
// JOINT LONGITUDINAL–SURVIVAL MODEL
// BASELINE HAZARD: BERNSTEIN POLYNOMIAL (BP)
// =============================================================================
//
// OVERVIEW:
// This Stan program implements a Bayesian joint model for:
//   - A continuous longitudinal outcome (linear mixed model)
//   - A time-to-event outcome with a flexible baseline hazard
//     represented using Bernstein polynomials
//
// This model is used as one option in a broader simulation and
// model-comparison framework (BP vs Gaussian Basis vs GP).
//
// HOW THIS FILE IS USED:
//   - This model is NOT standalone.
//   - Data are prepared in R via convert_data_from_models(), which:
//       * Fits an LME and a simple AFT survival model
//       * Constructs all design matrices and prior scales
//       * Passes the result as a fully-formed Stan data list
//   - This file is selected at runtime based on joint_model_type = "BP".
//
// BASELINE HAZARD (KEY IDEA):
//   - The baseline hazard is modeled as a Bernstein polynomial expansion
//     with m basis functions (m = k_bases in R).
//   - Bernstein weights (gamma) are constrained to be positive,
//     ensuring a valid hazard function.
//   - A second-order difference penalty encourages smoothness.
//
// JOINT MODEL LINKAGE:
//   - The survival submodel depends on the *current value* of the
//     longitudinal process via the association parameter alpha_tilde.
//   - No PH assumption is made; interpretation is AFT-style.
//
// IMPORTANT DESIGN NOTES:
//   - Time is assumed to be scaled to the study follow-up window.
//   - Priors are intentionally weakly informative and partially
//     data-adaptive (empirical priors supplied from R).
//   - This model is designed for simulation studies and comparison,
//     not as a final clinical analysis template.
//
// WHAT A READER SHOULD NOT WORRY ABOUT:
//   - The choice of k_bases: controlled externally via R.
//   - Data wrangling or centering: all handled upstream.
//   - Alternative link types (e.g. slope-based): not used here.
//
// =============================================================================


functions {
  // AFT joint survival log likelihood 
  vector loglik_aft_jm(
    vector time,
    vector beta_surv,
    vector beta_long,
    matrix b_long,
    vector gamma,
    vector status,
    matrix X_surv,
    real alpha,
    vector Y_long_surv
  )
 {
    int n = num_elements(status);
    int m = num_elements(gamma);
    vector[n] log_lik;
    vector[n] h0;
    vector[n] H0;
    vector[n] y;
    
    
    // assuming current value linkage between the longitudinal model and the AFT model
    // given the mean effect (without the random error):
    // Y*(t) = beta0 + beta1 * t + beta2 * t * arm + b0 + b1 * t
    // Equivalence within data generation is seen at around line 155 in Data_gen-I2a-nobD1.R
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
    vector[n] y_alt = fmin(fmax(y ./ tau_aft, eps), 1 - eps);

    matrix[n,m] b2;
    matrix[n,m] B2;
    
    for (k in 1:m) {
      for (i in 1:n) {
        b2[i,k] = beta_lpdf(y_alt[i] | k, (m - k + 1));
        B2[i,k] = beta_lcdf(y_alt[i] | k, (m - k + 1));
      }
    }
    
    // b2 = exp(b2) ./ tau_aft;
    // B2 = exp(B2);
    b2 = exp(b2); // test on 20251002
    B2 = exp(B2) .* tau_aft; // test on 20251002
    h0 = b2 * gamma;
    H0 = B2 * gamma;

    log_lik = ((log(h0) - (X_surv * beta_surv + alpha * Y_long_surv)) .* status) - H0;
    return log_lik;
  }
}  

// Data
data {
  //// survival data sizes (must come first because they are used for prior dims)
  int<lower=1> K_long; // number of population-level effects
  int<lower=1> q; // number of survival covariates
  int<lower=1> m; // Bernstein polynomial degree
  
  //// smoothing parameters
  int<lower=1> r;  // order of difference (1 = first difference, 2 = second difference, etc.)
  real<lower=0> tau_h;  // smoothing parameter

  // prior locations and scales
  vector[K_long] beta_long_prior_mean;
  vector<lower=0>[K_long] beta_long_prior_scale;
  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;
  
  // scale parameters for association parameter
  real s_long;  // scale for longitudinal process
  real s_surv;  // scale for survival linear predictor
  real<lower=0> alpha_tilde_sd;  // prior scale for alpha_tilde (clinically calibrated)

  //// longitudinal data
  int<lower=1> N_long; // total number of observations
  vector[N_long] Y_long; // response variable
  matrix[N_long, K_long] X_long; // population-level design matrix
  int<lower=1> N_1_long; // number of individual levels (should equate to number of individuals)
  array[N_long] int<lower=1> J_1_long; // group indicator from ID (rank, ordered)
  
  // individual/group-level predictor values (longitudinal)
  vector[N_long] Z_1_1_long; // intercept
  vector[N_long] Z_1_2_long; // time 

  //// survival data
  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=0>[n] time;
  matrix[n, q] X_surv;

  //// linking data
  array[n] int<lower=1> J_1_unique;
  matrix[n, K_long] X_long_surv;
}

// transformed data
transformed data {
  // ---- R-TH DIFFERENCE PENALTY MATRIX ----
  // Build the r-th difference penalty matrix Delta_r
  int n_diff = m - r;  // number of r-th differences
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
  //// longitudinal parameters
  vector[K_long] beta_long; // regression coefficients 
  real<lower=0> sigma_long; // dispersion parameter
  vector<lower=1e-3>[2] sd_1_long; // group/individual-level standard deviations
  matrix[2, N_1_long] z_1_long; // individual-level random effects
  cholesky_factor_corr[2] L_1_long; // cholesky factor of correlation matrix

  //// survival parameters
  
  vector[q] beta_surv;  
  vector<lower=0>[m] gamma; // BP basis weights for baseline hazard (arm = 0)

  //// linking parameter
  // association between longitudinal and survival submodels
  real alpha_tilde; 
  
}

// Transformed Parameters
transformed parameters {

  //// longitudinal transformed parameters
  matrix[2, N_1_long] b_raw = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';

  //// survival transformed parameters
  
  // Construct subject-specific predicted value Y_i(t) at event/censoring time
  // using fixed and random effects, to be used in survival model
  vector[n] Y_long_surv;
  for (i in 1:n) {
    Y_long_surv[i] = beta_long[1]
                   + beta_long[2] * X_long_surv[i, 2]
                   + beta_long[3] * X_long_surv[i, 3]
                   + b_long[J_1_unique[i], 1]
                   + b_long[J_1_unique[i], 2] * time[i];
  }
  
  real alpha = alpha_tilde * (s_surv / s_long);  // scale according to relative difference between longitudinal and survival

  // Construct log likelihood 
  vector[n] log_lik = loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv);

}

// Model
model {
  //// longitudinal priors
  beta_long ~ normal(beta_long_prior_mean, beta_long_prior_scale);
  sigma_long ~ cauchy(0, 5);
  sd_1_long ~ cauchy(0, 5);
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
  target += 0.5 * rho * log(tau_h);  // Normalization constant τ_h^(ρ/2)
  target += -0.5 * tau_h * quad_form(penalty_matrix, gamma);  // Quadratic penalty
  gamma ~ normal(0, 5) T[0, ];  // Keep base prior for positivity constraint

  // survival likelihood
  target += sum(log_lik);
}

