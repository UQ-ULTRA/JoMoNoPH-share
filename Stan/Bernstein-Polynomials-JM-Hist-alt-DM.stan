// =============================================================================
// JOINT LONGITUDINAL-SURVIVAL MODEL
// BASELINE HAZARD: BERNSTEIN POLYNOMIAL (BP)
// =============================================================================
//
// This file is statistically equivalent to Bernstein-Polynomials-JM-Hist.stan,
// but constructs the Bernstein basis directly rather than using beta_lpdf() and
// beta_cdf(). The Stan data interface, parameter names, priors, likelihood, and
// generated quantities follow Bernstein-Polynomials-JM-Hist.stan.
//
// =============================================================================


functions {
  // AFT joint survival log likelihood.
  // The original BP implementation used:
  //   beta_lpdf(y | k, m-k+1)
  //   beta_cdf(y | k, m-k+1)
  // For integer Bernstein bases, these are equivalent to:
  //   pdf: m * choose(m-1, k-1) * y^(k-1) * (1-y)^(m-k)
  //   cdf: sum_{j=k}^m choose(m, j) * y^j * (1-y)^(m-j)
  vector loglik_aft_jm(
    vector time,
    vector beta_surv,
    vector beta_long,
    matrix b_long,
    vector gamma,
    vector status,
    matrix X_surv,
    real alpha,
    vector Y_long_surv,
    vector bp_pdf_coef,
    vector bp_pdf_to_cdf_coef,
    int m
  )
 {
    int n = num_elements(status);
    vector[n] log_lik;
    vector[n] h0;
    vector[n] H0;
    vector[n] y;
    vector[n] linpred_surv;
    
    
    // assuming current value linkage between the longitudinal model and the AFT model
    // given the mean effect (without the random error):
    // Y*(t) = beta0 + beta1 * t + beta2 * t * arm + b0 + b1 * t
    // Equivalence within data generation is seen at around line 155 in Data_gen-I2a-nobD1.R
    vector[n] eta_surv = X_surv * beta_surv;
    vector[n] C1 = eta_surv + alpha * (beta_long[1] + b_long[,1]);
    vector[n] C2 = alpha * (beta_long[2] + X_surv[,1] * beta_long[3] + b_long[,2]);
    linpred_surv = eta_surv + alpha * Y_long_surv;
    
    for (i in 1:n){
      if (C2[i] != 0) {
        // Algebraically equivalent to the original expression below, but more
        // stable when C2[i] * time[i] is close to zero.
        y[i] = exp(-C1[i]) * (-expm1(-C2[i] * time[i])) / C2[i];
        // Original version:
        // y[i] = exp(-C1[i]) * (1 - exp(-C2[i] * time[i])) / C2[i];
      } else {
        y[i] = exp(-C1[i]) * time[i];
      }
    }
    
    real eps = 1e-6;
    real tau_aft = max(y) + eps;
    vector[n] y_alt = fmin(fmax(y ./ tau_aft, eps), 1 - eps);
    vector[n] omy = 1 - y_alt;
    
    h0 = rep_vector(0.0, n);
    H0 = rep_vector(0.0, n);

    // Direct hazard and cumulative hazard sums. Let y = y_alt and
    // B_{a,d}(y) = choose(d, a) * y^a * (1 - y)^(d - a).
    //
    // Baseline hazard:
    //   h0(y) = (m / tau_aft) *
    //           sum_{k=1}^m gamma[k] * B_{k-1,m-1}(y)
    // where this is equivalent to:
    //   sum_k gamma[k] * beta_pdf(y | k, m-k+1) / tau_aft.
    //
    // Cumulative baseline hazard:
    //   H0(y) = sum_{j=1}^m gamma_cum[j] * B_{j,m}(y),
    //   gamma_cum[j] = sum_{k=1}^j gamma[k].
    //
    // Instead of running a second loop over the degree-m basis for H0, use
    //   B_{k,m}(y) = (m / k) * y * B_{k-1,m-1}(y).
    // The j = 0 Bernstein term is not included because its integrated-basis
    // coefficient is zero.
    {
      vector[m] gamma_cum = cumulative_sum(gamma);
      vector[n] y_pw = rep_vector(1.0, n);
      vector[n] omy_pw = pow(omy, m - 1);
      for (k in 1:m) {
        vector[n] bern_pdf = bp_pdf_coef[k] * y_pw .* omy_pw;
        h0 += gamma[k] * bern_pdf;
        H0 += gamma_cum[k] * bp_pdf_to_cdf_coef[k] * (y_alt .* bern_pdf);
        if (k < m) {
          y_pw .*= y_alt;
          omy_pw ./= omy;
        }
      }
      h0 *= m / tau_aft;
    }

    h0 = fmax(h0, 1e-12);

    log_lik = ((log(h0) - linpred_surv) .* status) - H0;
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

// Transformed Data
transformed data {
  // Bernstein basis constants. These depend only on the basis degree, so
  // compute once rather than during every HMC evaluation.
  //
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

// Transformed Parameters
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
  
  // Regularise Bernstein baseline-hazard weights to avoid extreme flexible
  // baseline shapes that can drive very deep NUTS trajectories in simulations.
  // Previous weaker prior:
  // gamma ~ normal(0, 5) T[0, ];
  gamma ~ normal(0, 1) T[0, ];

  // Survival likelihood. The longitudinal likelihood is already included above
  // through Y_long ~ normal(...); the full joint log_lik is generated only for IC.
  target += sum(loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv, bp_pdf_coef, bp_pdf_to_cdf_coef, m));
}

generated quantities {
  // Full patient-level joint log likelihood for PSIS-LOO/WAIC.
  // This is intentionally generated after sampling, not used in target, to avoid
  // double-counting the longitudinal likelihood and unnecessary HMC work.
  vector[n] log_lik;
  {
    vector[n] log_lik_surv = loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv, bp_pdf_coef, bp_pdf_to_cdf_coef, m);
    vector[n] log_lik_long = rep_vector(0, n);
    vector[N_long] mu_long_ic = X_long * beta_long
                              + b_long[J_1_long, 1] .* Z_1_1_long
                              + b_long[J_1_long, 2] .* Z_1_2_long;

    for (j in 1:N_long) {
      // J_1_long maps each repeated longitudinal observation back to its patient.
      log_lik_long[J_1_long[j]] += normal_lpdf(Y_long[j] | mu_long_ic[j], sigma_long);
    }

    log_lik = log_lik_surv + log_lik_long;
  }
}
