
// Functions
functions {

  // aft joint survival log likelihood function
  vector loglik_aft_jm(
    vector time, // survival times
    vector beta_surv, // survival model population-level effects
    vector beta_long, // longitudinal model population-level effects
    matrix b_long, // longitudinal model group-level effects
    vector gamma, // baseline hazard coefficients
    vector status, // event indicators
    matrix X_surv, // survival model design matrix 
    real alpha, // association parameter
    vector Y_long_surv, // longitudinal model fitted values at survival times
    vector bp_pdf_coef, // bernstein basis coefficients for the pdf (binomial coefficients)
    vector bp_pdf_to_cdf_coef, // coefficients to convert from pdf to cdf basis (m / k)
    int m // degree of the bernstein polynomial for the baseline hazard
  )
 {
    int n = num_elements(status); // number of individuals (length of survival dataset)
    vector[n] log_lik; // log-likelihood
    vector[n] h0; // baseline hazard
    vector[n] H0; // cumulative baseline hazard
    vector[n] kappa; // accelerated time
    vector[n] linpred_surv; // survival linear predictor
    
    
    // assuming current value linkage between the longitudinal model and the AFT model
    // given the mean effect (without the random error):
    // Y*(t) = beta0 + beta1 * t + beta2 * t * arm + b0 + b1 * t
    vector[n] eta_surv = X_surv * beta_surv;
    vector[n] C1 = eta_surv + alpha * (beta_long[1] + b_long[,1]);
    vector[n] C2 = alpha * (beta_long[2] + X_surv[,1] * beta_long[3] + b_long[,2]);
    linpred_surv = eta_surv + alpha * Y_long_surv;
    
    for (i in 1:n){
      if (C2[i] != 0) {
        kappa[i] = exp(-C1[i]) * (-expm1(-C2[i] * time[i])) / C2[i];
      } else {
        kappa[i] = exp(-C1[i]) * time[i];
      }
    }
    
    real eps = 1e-6; // small constant
    real tau_aft = max(kappa) + eps; // maximum kappa
    vector[n] kappa_scaled = fmin(fmax(kappa ./ tau_aft, eps), 1 - eps);
    vector[n] omy = 1 - kappa_scaled;
    
    // initialise baseline and cumulative baseline hazards as zero vectors
    h0 = rep_vector(0.0, n); 
    H0 = rep_vector(0.0, n);

    vector[m] gamma_cum = cumulative_sum(gamma);
    vector[n] y_pw = rep_vector(1.0, n);
    vector[n] omy_pw = pow(omy, m - 1);
    
    
    for (k in 1:m) {
      vector[n] bern_pdf = bp_pdf_coef[k] * y_pw .* omy_pw;
      h0 += gamma[k] * bern_pdf;
      H0 += gamma_cum[k] * bp_pdf_to_cdf_coef[k] * (kappa_scaled .* bern_pdf);
      if (k < m) {
        y_pw .*= kappa_scaled;
        omy_pw ./= omy;
      }
    }
    
    h0 *= m / tau_aft;
    h0 = fmax(h0, 1e-12);

    log_lik = ((log(h0) - linpred_surv) .* status) - H0;
    return log_lik;
  }
}  

// Data
data {
  
  //// survival data sizes 
  int<lower=1> K_long; // number of population-level effects
  int<lower=1> q; // number of survival covariates
  int<lower=1> m; // Bernstein polynomial degree

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

// Transformed data
transformed data {
  
  // construction of the bernstein basis coefficients, beta/binomal coefficients that are used to normalise cdf/pdf
  // bp_pdf_coef[k] = choose(m - 1, k - 1), used for B_{k-1,m-1}(y).
  // bp_pdf_to_cdf_coef[k] = m / k, because
  // B_{k,m}(y) = (m / k) * y * B_{k-1,m-1}(y).
  vector[m] bp_pdf_coef;
  vector[m] bp_pdf_to_cdf_coef;
  for (k in 1:m) {
    bp_pdf_coef[k] = choose(m - 1, k - 1);
    bp_pdf_to_cdf_coef[k] = m * 1.0 / k;
  }
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

// Transformed parameters
transformed parameters {

  //// longitudinal transformed parameters
  matrix[2, N_1_long] b_raw = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';

  //// survival transformed parameters
  
  // Construct the subject-specific longitudinal mean at event/censoring time
  // using fixed and random effects, to be used in the survival model.
  vector[n] Y_long_surv = rep_vector(beta_long[1], n)
                        + beta_long[2] * X_long_surv[, 2]
                        + beta_long[3] * X_long_surv[, 3]
                        + b_long[J_1_unique, 1]
                        + b_long[J_1_unique, 2] .* time;
  
  real alpha = alpha_tilde * (s_surv / s_long);  // scale according to relative difference between longitudinal and survival

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
  vector[N_long] mu_long = X_long * beta_long
                         + b_long[J_1_long, 1] .* Z_1_1_long
                         + b_long[J_1_long, 2] .* Z_1_2_long;
  Y_long ~ normal(mu_long, sigma_long);

  //// survival priors
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  alpha_tilde ~ normal(0, alpha_tilde_sd);
  
  // bernstein polynomial weight
  gamma ~ normal(0, 1) T[0, ];

  // survival likelihood
  target += sum(loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv, bp_pdf_coef, bp_pdf_to_cdf_coef, m));
}

// Generated quantities
generated quantities {
  
  // participant-level joint log likelihood generation, used for LOO/WAIC statistics
  vector[n] log_lik;
  vector[n] log_lik_surv = loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv, bp_pdf_coef, bp_pdf_to_cdf_coef, m);
  vector[n] log_lik_long = rep_vector(0, n);
  vector[N_long] mu_long_ic = X_long * beta_long + b_long[J_1_long, 1] .* Z_1_1_long + b_long[J_1_long, 2] .* Z_1_2_long;

  for (j in 1:N_long) {
    log_lik_long[J_1_long[j]] += normal_lpdf(Y_long[j] | mu_long_ic[j], sigma_long);
  }

  log_lik = log_lik_surv + log_lik_long;
}
