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
        // d/dy Phi((log y - mu) / sigma) = Normal(log y | mu, sigma) / y
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
  int<lower=1> q; // number of survival covariates
  int<lower=1> m; // number of Gaussian basis functions

  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;

  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=1e-12>[n] time;
  matrix[n, q] X_surv;

  // Fixed knot span on the log(y) scale.
  real log_y_lower;
  real log_y_upper;
}

transformed data {
  vector[m] mu_fix;
  vector[m] sigma_fix;
  real span = log_y_upper - log_y_lower;

  if (span <= 0) {
    reject("log_y_upper must be greater than log_y_lower.");
  }

  if (m == 1) {
    mu_fix[1] = 0.5 * (log_y_lower + log_y_upper);
    sigma_fix[1] = 0.5 * span;
  } else {
    for (k in 1:m) {
      mu_fix[k] = log_y_lower + (k - 1.0) * span / (m - 1.0);
    }
    sigma_fix = rep_vector(2.0 * span / (3.0 * (m - 1.0)), m);
  }
}

parameters {
  vector[q] beta_surv;
  vector<lower=0>[m] gamma;
}

transformed parameters {
  vector[n] log_lik = loglik_aft_logy(
    time, beta_surv, gamma, status, X_surv, mu_fix, sigma_fix
  );
}

model {
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  gamma ~ normal(0, 5) T[0, ];

  target += sum(log_lik);
}
