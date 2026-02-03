
// Functions
functions {
  // AFT joint survival log likelihood 
  vector loglik_aft(
    vector time,
    vector beta_surv,
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
    
    y = exp(-(X_surv*beta_surv))  .* time;
    
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

    log_lik = ((log(h0) - (X_surv * beta_surv)) .* status) - H0;
    return log_lik;
  }
  
}  

// Data
data {
  //// survival data sizes (must come first because they are used for prior dims)
  int<lower=1> q; // number of survival covariates
  int<lower=1> m; // Bernstein polynomial degree

  // prior locations and scales
  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;
  
  //// survival data
  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=0>[n] time;
  matrix[n, q] X_surv;
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
  vector[n] log_lik = loglik_aft(time, beta_surv, gamma, status, X_surv);

}

// Model
model {
  //// survival priors
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  gamma ~ normal(0, 5) T[0, ];  // Keep base prior for positivity constraint

  // survival likelihood
  target += sum(log_lik);
}

