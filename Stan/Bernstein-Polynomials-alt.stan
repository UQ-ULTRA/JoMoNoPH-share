
// Functions
functions {
  ////// bernstein polynomial functions (A-star)
  
  // H0 (need to call H0 first since h0 contingent on H0)
  vector bp_H0(vector x, vector theta){
    
    int m = num_elements(theta)-1;
    int N = num_elements(x);
    
    real eps = 1e-12;
    vector[N] x_star = fmin(1-eps, fmax(eps,x));
    
    vector[N] H0 = rep_vector(0.0, N);
    
    for (j in 0:m){
      H0 += theta[j+1]*choose(m,j).*pow(x_star, j).*pow(1 - x_star,m-j);
    }
    
    return H0; 
  }

  // log h0
  vector bp_log_h0(vector x, vector theta, real tau, vector lin_pred){

    int m = num_elements(theta)-1;
    int N = num_elements(x);
    
    real eps = 1e-12;
    vector[N] x_star = fmin(1 - eps, fmax(eps, x));
    
    vector[N] log_h0 = rep_vector(-log(tau)+ log(m), N) - lin_pred;
    vector[N] nat_sum = rep_vector(0.0, N); 
    
    // sum on the natural scale
    for (j in 1:m){
      real delta_j = theta[j+1]-theta[j]; // should be >0 due to monotonicity
      nat_sum += delta_j*choose(m-1, j-1).*pow(x_star, j-1).*exp((m-j) * log1m(x_star));
    }
    
    log_h0 += log(nat_sum+eps);
    return log_h0; 
  }
  
  // h0
  vector bp_h0(vector x, vector theta, real tau, vector lin_pred){

    int m = num_elements(theta)-1;
    int N = num_elements(x);
    
    real eps = 1e-12;
    vector[N] x_star = fmin(1 - eps, fmax(eps, x));
    
    vector[N] h0 = rep_vector((1/tau)*m,N).*exp(-lin_pred);
    vector[N] nat_sum = rep_vector(0.0, N); 
    
    // sum on the natural scale
    for (j in 1:m){
      real delta_j = theta[j+1]-theta[j]; // should be >0 due to monotonicity
      nat_sum += delta_j*choose(m-1, j-1).*pow(x_star, j-1).*exp((m-j) * log1m(x_star));
    }
    
    h0 = h0 .* (nat_sum+eps);
    return h0; 
  }

  // p(t_i, delta_i | X_i, beta)
  vector survival_loglik_bp(vector time, vector delta, vector lin_pred, vector theta){
    int N = num_elements(time);

    vector[N] log_baseline_hazard;
    vector[N] log_survival;
    
    real eps = 1e-12;
    vector[N] kappa = exp(-lin_pred).*time;
    real tau_aft = max(kappa) + eps;
    vector[N] phi = kappa/tau_aft; // scaled to [0,1]

    log_baseline_hazard = bp_log_h0(phi, theta, tau_aft, lin_pred);
    log_survival = -bp_H0(phi, theta);

    return delta.*(log_baseline_hazard) + log_survival;
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
  // bp params
  ordered[m+1] theta_raw; // raw bernstein polynomials weights (equiv to thetas)
  real<lower=0> gamma_bp; // bp scaling parameter
  
}

// Transformed Parameters
transformed parameters {
  // Construct log likelihood 
  // vector[n] log_lik = loglik_aft(time, beta_surv, gamma, status, X_surv);
  
  vector[n] lin_pred = X_surv * beta_surv;
  
    // bp weights
  vector[m+1] theta = gamma_bp*(theta_raw-theta_raw[1])/(theta_raw[m+1]-theta_raw[1] +1e-12);

}

// Model
model {
  //// survival priors
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  
  // bp
  gamma_bp ~ lognormal(0, 1);
  theta_raw ~ normal(0, 1);   // allowed because ordered supports it elementwise

  // survival likelihood
  // target += sum(log_lik);
  
  target += sum(survival_loglik_bp(time, status, lin_pred, theta));
}

