// =============================================================================
// JOINT LONGITUDINAL-SURVIVAL MODEL
// BASELINE HAZARD: BERNSTEIN POLYNOMIAL (BP)
// GAMMA PRIOR SENSITIVITY: half-normal SD 0.5
// =============================================================================
//
// This file is the no-shape-penalty gamma-prior sensitivity version of
// Bernstein-Polynomials-JM-Hist.stan.
//
// Difference from the current standard gamma1 model:
//   gamma ~ normal(0, 0.5) T[0, ];
//
// No normalized-gamma shape roughness penalty is used here.
// Sampler controls such as max_treedepth and adapt_delta are set in R.
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
    vector[n] linpred_surv;

    vector[n] C1 = X_surv * beta_surv + alpha * (beta_long[1] + b_long[,1]);
    vector[n] C2 = alpha * (beta_long[2] + X_surv[,1] * beta_long[3] + b_long[,2]);
    linpred_surv = X_surv * beta_surv + alpha * Y_long_surv;

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

    // Vectorised over i for each k: eliminates n*m scalar beta_lpdf calls.
    // beta_cdf does not vectorise in this Stan version, so inner loop retained.
    for (k in 1:m) {
      real lbeta_k = lbeta(k, m - k + 1);
      b2[, k] = (k - 1) * log(y_alt) + (m - k) * log1m(y_alt) - lbeta_k;
      for (i in 1:n) {
        B2[i, k] = beta_cdf(y_alt[i] | k, m - k + 1);
      }
    }

    b2 = exp(b2) ./ tau_aft;
    h0 = fmax(b2 * gamma, 1e-12);
    H0 = B2 * gamma;

    log_lik = ((log(h0) - linpred_surv) .* status) - H0;
    return log_lik;
  }
}

data {
  int<lower=1> K_long;
  int<lower=1> q;
  int<lower=1> m;

  int<lower=1> r;
  real<lower=0> tau_h;

  vector[K_long] beta_long_prior_mean;
  vector<lower=0>[K_long] beta_long_prior_scale;
  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;

  real s_long;
  real s_surv;
  real<lower=0> alpha_tilde_sd;

  int<lower=1> N_long;
  vector[N_long] Y_long;
  matrix[N_long, K_long] X_long;
  int<lower=1> N_1_long;
  array[N_long] int<lower=1> J_1_long;

  vector[N_long] Z_1_1_long;
  vector[N_long] Z_1_2_long;

  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=0>[n] time;
  matrix[n, q] X_surv;

  array[n] int<lower=1> J_1_unique;
  matrix[n, K_long] X_long_surv;
}

parameters {
  vector[K_long] beta_long;
  real<lower=0> sigma_long;
  vector<lower=1e-3>[2] sd_1_long;
  matrix[2, N_1_long] z_1_long;
  cholesky_factor_corr[2] L_1_long;

  vector[q] beta_surv;
  vector<lower=0>[m] gamma;

  real alpha_tilde;
}

transformed parameters {
  matrix[2, N_1_long] b_raw = diag_pre_multiply(sd_1_long, L_1_long) * z_1_long;
  matrix[N_1_long, 2] b_long = b_raw';

  vector[n] Y_long_surv;
  for (i in 1:n) {
    Y_long_surv[i] = beta_long[1]
                   + beta_long[2] * X_long_surv[i, 2]
                   + beta_long[3] * X_long_surv[i, 3]
                   + b_long[J_1_unique[i], 1]
                   + b_long[J_1_unique[i], 2] * time[i];
  }

  real alpha = alpha_tilde * (s_surv / s_long);

}

model {
  beta_long ~ normal(beta_long_prior_mean, beta_long_prior_scale);
  sigma_long ~ cauchy(0, 5);
  sd_1_long ~ cauchy(0, 5);
  L_1_long ~ lkj_corr_cholesky(2);
  to_vector(z_1_long) ~ std_normal();

  vector[N_long] mu_long = X_long * beta_long;
  for (i in 1:N_long) {
    mu_long[i] += b_long[J_1_long[i], 1] * Z_1_1_long[i]
                + b_long[J_1_long[i], 2] * Z_1_2_long[i];
  }
  Y_long ~ normal(mu_long, sigma_long);

  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  alpha_tilde ~ normal(0, alpha_tilde_sd);

  // Tighter half-normal prior for sensitivity testing.
  // Previous priors tested:
  // gamma ~ normal(0, 5) T[0, ];
  // gamma ~ normal(0, 1) T[0, ];
  gamma ~ normal(0, 0.5) T[0, ];

  // Survival likelihood. The longitudinal likelihood is already included above
  // through Y_long ~ normal(...); the full joint log_lik is generated only for IC.
  target += sum(loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv));
}

generated quantities {
  // Full patient-level joint log likelihood for PSIS-LOO/WAIC.
  // This is intentionally generated after sampling, not used in target, to avoid
  // double-counting the longitudinal likelihood and unnecessary HMC work.
  vector[n] log_lik;
  {
    vector[n] log_lik_surv = loglik_aft_jm(time, beta_surv, beta_long, b_long, gamma, status, X_surv, alpha, Y_long_surv);
    vector[n] log_lik_long = rep_vector(0, n);
    vector[N_long] mu_long_ic = X_long * beta_long;

    for (j in 1:N_long) {
      mu_long_ic[j] += b_long[J_1_long[j], 1] * Z_1_1_long[j]
                    + b_long[J_1_long[j], 2] * Z_1_2_long[j];
      // J_1_long maps each repeated longitudinal observation back to its patient.
      log_lik_long[J_1_long[j]] += normal_lpdf(Y_long[j] | mu_long_ic[j], sigma_long);
    }

    log_lik = log_lik_surv + log_lik_long;
  }
}
