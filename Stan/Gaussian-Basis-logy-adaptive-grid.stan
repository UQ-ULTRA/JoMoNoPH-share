functions {
  // AFT survival log likelihood with Gaussian bases on log(y)
  vector loglik_aft_logy(
    vector time,
    vector beta_surv,
    vector gamma,
    vector status,
    matrix X_surv,
    vector mu,
    vector sigma
  ) {
    int n = num_elements(status);
    int m = num_elements(gamma);
    vector[n] log_lik;
    vector[n] h0;
    vector[n] H0;
    vector[n] log_y;
    vector[n] y;

    log_y = log(time) - (X_surv * beta_surv);
    y = exp(log_y);

    matrix[n, m] f_mat;
    matrix[n, m] F_mat;
    for (k in 1:m) {
      real sk = fmax(sigma[k], 1e-6);
      for (i in 1:n) {
        f_mat[i, k] = exp(normal_lpdf(log_y[i] | mu[k], sk)) / y[i];
        F_mat[i, k] = exp(normal_lcdf(log_y[i] | mu[k], sk));
      }
    }

    h0 = f_mat * gamma;
    H0 = F_mat * gamma;
    log_lik = ((log(h0) - (X_surv * beta_surv)) .* status) - H0;
    return log_lik;
  }
}

data {
  int<lower=1> q;
  int<lower=1> m;

  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;

  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=1e-12>[n] time;
  matrix[n, q] X_surv;

  // Base grid on log(y) before adaptive shift/stretch.
  real grid_lower;
  real grid_upper;

  // Priors for mu_k = a + b * grid_k, b > 0.
  real a_prior_mean;
  real<lower=0> a_prior_scale;
  real log_b_prior_mean;
  real<lower=0> log_b_prior_scale;
}

transformed data {
  vector[m] grid_raw;

  if (grid_upper <= grid_lower) {
    reject("grid_upper must be greater than grid_lower.");
  }

  if (m == 1) {
    grid_raw[1] = 0.5 * (grid_lower + grid_upper);
  } else {
    for (k in 1:m) {
      grid_raw[k] = grid_lower + (k - 1.0) * (grid_upper - grid_lower) / (m - 1.0);
    }
  }
}

parameters {
  vector[q] beta_surv;
  vector<lower=0>[m] gamma;

  real a;
  real<lower=-5, upper=5> log_b;
}

transformed parameters {
  vector[m] mu_fix;
  vector[m] sigma_fix;
  real b = exp(log_b);
  real base_spacing;
  vector[n] log_lik;

  mu_fix = a + b * grid_raw;

  if (m == 1) {
    sigma_fix[1] = 0.5 * b * (grid_upper - grid_lower);
  } else {
    base_spacing = (grid_upper - grid_lower) / (m - 1.0);
    sigma_fix = rep_vector(2.0 * b * base_spacing / 3.0, m);
  }

  log_lik = loglik_aft_logy(
    time, beta_surv, gamma, status, X_surv, mu_fix, sigma_fix
  );
}

model {
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  gamma ~ normal(0, 5) T[0, ];

  a ~ normal(a_prior_mean, a_prior_scale);
  log_b ~ normal(log_b_prior_mean, log_b_prior_scale);

  target += sum(log_lik);
}
