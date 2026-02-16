
// Functions
functions {
  // AFT joint survival log likelihood
  vector loglik_aft(
  vector time,
  vector beta_surv,
  vector gamma,
  vector status,
  matrix X_surv,
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

    y = exp(-(X_surv*beta_surv))  .* time;
    
    real eps = 1e-6;
    real tau_aft = max(y) + eps;
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

    log_lik = ((log(h0) - (X_surv * beta_surv)) .* status) - H0;
    return log_lik;
  }
}


// Data
data {
  //// survival data sizes (must come first because they are used for prior dims)
  int<lower=1> q; // number of survival covariates
  int<lower=1> m; // number of Gaussian basis functions

  // prior locations and scales
  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;
  
  //// survival data
  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=0>[n] time;
  matrix[n, q] X_surv;
}

// Transformed Data
transformed data {
  vector[m] mu_fix;
  vector[m] sigma_fix;
  
  mu_fix[1] = 0;
  mu_fix[m] = 1;
  
  sigma_fix = rep_vector(2.0/(3.0*(m-1.0)), m);    
  
  for (k in 1:(m-2)) {
    mu_fix[k+1] = k * 1.0 / (m-1);  // equally spaced in (0,1)
  }
}

// Parameters
parameters {
  //// survival parameters
  vector[q] beta_surv;  
  vector<lower=0>[m] gamma; // BP basis weights for baseline hazard (arm = 0)
  
}

// Transformed Parameters
transformed parameters {
  // Construct log likelihood 
  vector[n] log_lik = loglik_aft(time, beta_surv, gamma, status, X_surv, mu_fix, sigma_fix);

}

// Model
model {
  //// survival priors
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  gamma ~ normal(0, 5) T[0, ];  // Keep base prior for positivity constraint

  // survival likelihood
  target += sum(log_lik);
}

