functions {
  // AFT joint survival log likelihood 
  vector loglik_aft(
    vector time,
    vector beta_surv,
    vector beta_long,
    matrix b_long,
    vector gamma,
    vector status,
    matrix X_surv
  )
 {
    int n = num_elements(status);
    int m = num_elements(gamma);
    vector[n] log_lik;
    vector[n] h0;
    vector[n] H0;
    vector[n] y;
    
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
  
  //// survival data
  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=0>[n] time;
  matrix[n, q] X_surv;

  //// linking data
  array[n] int<lower=1> J_1_unique;
  matrix[n, K_long] X_long_surv;
}

// Parameters
parameters {
  //// survival parameters
  vector[q] beta_surv;  
  vector<lower=0>[m] gamma; // BP basis weights for baseline hazard (arm = 0)
  
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

