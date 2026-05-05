functions {
  // AFT survival-only log-likelihood with Bernstein polynomials on
  // dynamically normalized log accelerated time:
  //   y = t * exp(-x beta)
  //   u = (log(y) - min(log(y))) / (max(log(y)) - min(log(y)))
  //
  // This mirrors the dynamic scaling idea in Bernstein-Polynomials.stan,
  // but uses log accelerated time so the normalizing range is less volatile.
  vector loglik_aft_logy(
    vector time,
    vector beta_surv,
    vector gamma,
    vector status,
    matrix X_surv
  ) {
    int n = num_elements(status);
    int m = num_elements(gamma);
    real eps = 1e-6;
    real log_y_lower;
    real span;
    vector[n] linpred_surv;
    vector[n] log_y;
    vector[n] u;
    vector[n] log_lik;
    vector[n] h_u;
    vector[n] H0;
    matrix[n, m] b2;
    matrix[n, m] B2;

    linpred_surv = X_surv * beta_surv;
    log_y = log(time) - linpred_surv;
    log_y_lower = min(log_y);
    span = fmax(max(log_y) - log_y_lower, eps);
    u = fmin(fmax((log_y - log_y_lower) / span, eps), 1 - eps);

    for (k in 1:m) {
      for (i in 1:n) {
        b2[i, k] = exp(beta_lpdf(u[i] | k, m - k + 1));
        B2[i, k] = beta_cdf(u[i] | k, m - k + 1);
      }
    }

    h_u = fmax(b2 * gamma, 1e-12);
    H0 = B2 * gamma;

    // Change of variables:
    //   h_Y(y) = h_U(u) / (span * y)
    // and h_T(t | x) = h_Y(t * exp(-x beta)) * exp(-x beta)
    //               = h_U(u) / (span * t).
    log_lik = ((log(h_u) - log(span) - log(time)) .* status) - H0;
    return log_lik;
  }
}

data {
  int<lower=1> q; // number of survival covariates
  int<lower=1> m; // number of Bernstein basis terms

  vector[q] beta_surv_prior_mean;
  vector<lower=0>[q] beta_surv_prior_scale;

  int<lower=1> n;
  vector<lower=0, upper=1>[n] status;
  vector<lower=1e-12>[n] time;
  matrix[n, q] X_surv;
}

parameters {
  vector[q] beta_surv;
  vector<lower=0>[m] gamma;
}

transformed parameters {
  vector[n] log_lik = loglik_aft_logy(
    time, beta_surv, gamma, status, X_surv
  );
}

model {
  beta_surv ~ normal(beta_surv_prior_mean, beta_surv_prior_scale);
  gamma ~ normal(0, 5) T[0, ];

  target += sum(log_lik);
}
