// Functions
functions {
    ////// bernstein polynomial functions (A-star)
  
  // H0 (need to call H0 first since h0 contingent on H0)
  vector bp_H0(vector x, vector theta){
    
    int m = num_elements(theta)-1;
    int N = num_elements(x);
    
    real eps = 1e-12;
    vector[N] x_star = fmin(1-eps, fmax(eps,x));
    vector[N] omx = 1 - x_star;
    
    vector[N] H0 = rep_vector(0.0, N);
    vector[N] x_pw = rep_vector(1.0, N);  // x_star^j,   starts at j=0
    vector[N] omx_pw = pow(omx, m);         // (1-x_star)^(m-j), starts at j=0
    
    for (j in 0:m){
      H0 += theta[j+1] * choose(m,j) * x_pw .* omx_pw;
      if (j < m) {
        x_pw .*= x_star;  // advance to x_star^(j+1)
        omx_pw ./= omx;     // advance to (1-x_star)^(m-j-1)
      }
    }
    
    return H0; 
  }

  // log h0
  vector bp_log_h0(vector x, vector theta, real tau, vector lin_pred){

    int m = num_elements(theta)-1;
    int N = num_elements(x);
    
    real eps = 1e-12;
    vector[N] x_star = fmin(1 - eps, fmax(eps, x));
    vector[N] omx    = 1 - x_star;
    
    vector[N] log_h0  = rep_vector(-log(tau) + log(m), N) - lin_pred;
    vector[N] nat_sum = rep_vector(0.0, N);
    vector[N] x_pw    = rep_vector(1.0, N);   // x_star^(j-1), starts at j=1
    vector[N] omx_pw  = pow(omx, m - 1);      // (1-x_star)^(m-j), starts at j=1
    
    // sum on the natural scale
    for (j in 1:m){
      real delta_j = theta[j+1] - theta[j]; // should be >0 due to monotonicity
      nat_sum += delta_j * choose(m-1, j-1) * x_pw .* omx_pw;
      if (j < m) {
        x_pw   .*= x_star;  // advance to x_star^j
        omx_pw ./= omx;     // advance to (1-x_star)^(m-j-1)
      }
    }
    
    log_h0 += log(nat_sum + eps);
    return log_h0; 
  }
  
  // h0
  vector bp_h0(vector x, vector theta, real tau, vector lin_pred){

    int m = num_elements(theta)-1;
    int N = num_elements(x);
    
    real eps = 1e-12;
    vector[N] x_star = fmin(1 - eps, fmax(eps, x));
    vector[N] omx = 1 - x_star;
    
    vector[N] h0 = rep_vector((1/tau)*m, N) .* exp(-lin_pred);
    vector[N] nat_sum = rep_vector(0.0, N);
    vector[N] x_pw = rep_vector(1.0, N);   // x_star^(j-1), starts at j=1
    vector[N] omx_pw = pow(omx, m - 1);      // (1-x_star)^(m-j), starts at j=1
    
    // sum on the natural scale
    for (j in 1:m){
      real delta_j = theta[j+1] - theta[j]; // should be >0 due to monotonicity
      nat_sum += delta_j * choose(m-1, j-1) * x_pw .* omx_pw;
      if (j < m) {
        x_pw .*= x_star;
        omx_pw ./= omx;
      }
    }
    
    h0 = h0 .* (nat_sum + eps);
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
  
  
    // p(t_i, delta_i | X_i, beta)
  vector survival_loglik_bp_jm(vector time, vector delta, vector lin_pred, vector theta){
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
  
  //// longitudinal data
  int<lower=1> K_long; // number of population-level effects
  int<lower=1> N_long; // total number of observations
  vector[N_long] Y_long; // response variable
  matrix[N_long, K_long] X_long; // population-level design matrix
  int<lower=1> N_1_long; // number of individual levels (should equate to number of individuals)
  array[N_long] int<lower=1> J_1_long; // group indicator from ID (rank, ordered)
  
  // individual/group-level predictor values (longitudinal)
  vector[N_long] Z_1_1_long; // intercept
  vector[N_long] Z_1_2_long; // time 
  
  
  //// survival data
  int<lower=1> N_surv; // number of observations
  int<lower=1> K_surv; // number of covariates
  vector[N_surv] y_surv; // observed times
  matrix[N_surv,K_surv] X_surv; // covariate matrix
  vector[N_surv] delta; // censoring indicators [1: event, 0: censored]
  
  // prior locations and scales
  vector[K_long] beta_long_prior_mean;
  vector<lower=0>[K_long] beta_long_prior_scale;
  vector[K_surv] beta_surv_prior_mean;
  vector<lower=0>[K_surv] beta_surv_prior_scale;
  
  // aft mode
  // int<lower=0> mode; // 0: weibull, 1: loglogistic, 2: bernstein polynomial, 4: gaussian process
  
  // bp specific
  int<lower=1> m; // order of bernstein polynomial
  
  
  //// linking data
  array[N_surv] int<lower=1> J_1_unique;
  matrix[N_surv, K_long] X_long_surv;
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
  
  // real alpha; // intercept
  vector[K_surv] gamma; // regression coefficients

  // bp params
  ordered[m+1] theta_raw; // raw bernstein polynomials weights
  real<lower=0> gamma_bp; // bp scaling parameter
  // vector[N] phi = exp(-lin_pred).* y/tau; // scaled times for bp
  
  
  //// linking parameter
  real alpha; // association parameter
}

// Transformed parameters
transformed parameters{
  
  //// longitudinal transformed parameters
  matrix[2, N_1_long] b_raw = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';
  
  vector[N_surv] lin_pred; // linear predictior
  vector[N_surv] Y_long_surv; // predicted y_i(t) value at event/censoring time
  
  for (i in 1:N_surv) {
    Y_long_surv[i] = beta_long[1]
                   + beta_long[2] * X_long_surv[i, 2]
                   + beta_long[3] * X_long_surv[i, 3]
                   + b_long[J_1_unique[i], 1]
                   + b_long[J_1_unique[i], 2] * y_surv[i];
  }
  
  
  // lin_pred = alpha+X*gamma;
  lin_pred = X_surv*gamma + Y_long_surv*alpha;
  
  // bp weights
  vector[m+1] theta = gamma_bp*(theta_raw-theta_raw[1])/(theta_raw[m+1]-theta_raw[1] +1e-12); 

}

// Model
model {
  
  //// longitudinal priors
  beta_long ~ normal(beta_long_prior_mean, beta_long_prior_scale);
  sd_1_long ~ cauchy(0, 5);
  sigma_long ~ cauchy(0, 5);
  L_1_long ~ lkj_corr_cholesky(2);
  to_vector(z_1_long) ~ std_normal();
  
  // longitudinal likelihood
  vector[N_long] mu_long = X_long * beta_long;
  for (i in 1:N_long) {
    mu_long[i] += b_long[J_1_long[i], 1] * Z_1_1_long[i]
                + b_long[J_1_long[i], 2] * Z_1_2_long[i];
  }
  Y_long ~ normal(mu_long, sigma_long);
  
  //// time to event priors 
  alpha ~ normal(0,1);
  gamma ~ normal(beta_surv_prior_mean, beta_surv_prior_scale); // regression coefficients 
  
  
  // bp priors
  gamma_bp ~ cauchy(0, 3);
  theta_raw ~ normal(0, 1);   // allowed because ordered supports it elementwise
  
  //// likelihood
  
  // Bernstein polynomial mode A*
  target += sum(survival_loglik_bp_jm(y_surv, delta, lin_pred, theta));

}


