
# Load packages
library(tidyverse) # for data manipulation and visualization
library(cmdstanr) # for Bayesian modeling

# Inverse CDF functions for survival time generation
ll_inv <- function(u, shape, scale) scale*((u / (1 - u))^(1 / shape))
wb_inv <- function(u, shape, scale) scale*(-log(1 - u))^(1 / shape)

# Primary simulation function
simulate_joint_dataset <- function(D, beta_0, beta_1, beta_2, sigma_e,
                                   log_AF, alpha_AFT,
                                   loglogistic_shape, loglogistic_scale,
                                   visit, seed, n_patients, max_FU, lambda_c = -1,
                                   aft_mode = "loglogistic", link_type = "value",
                                   ...) {  # absorb unused PH/Weibull placeholders
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- n_patients
  arm <- rep(0:1, each = n / 2)
  b <- mvrnorm(n, mu = c(0, 0), Sigma = D)
  b0 <- b[, 1]; b1 <- b[, 2]
  
  # Acceleration components (vectorized over participants)
  # Centre the current-value link at beta_0 so the baseline population mean
  # does not feed into the AFT acceleration; only individual deviations (b0)
  # from that reference drive the initial hazard.
  Y_ref <- beta_0
  eta0 <- if (link_type == "value") b0 else rep(0, n)
  C1 <- log_AF * arm + alpha_AFT * eta0
  C2 <- if (link_type %in% c("value", "slope")) {
    alpha_AFT * (beta_1 + beta_2 * arm + b1)
  } else {
    rep(0, n)
  }
  
  # Inverse-sample base survival times
  U <- runif(n, 1e-6, 1 - 1e-6)
  kappa <- if (aft_mode == "loglogistic") {
    ll_inv(U, loglogistic_shape, loglogistic_scale)
  } else {
    wb_inv(U, loglogistic_shape, loglogistic_scale)  # reuse shape/scale args
  }
  
  # Survival times: T = -log(1 - C2*exp(C1)*kappa) / C2 when C2 != 0
  A <- pmin(C2 * exp(C1) * kappa, 1 - 1e-8)
  T_i <- ifelse(abs(C2) < 1e-8, exp(C1) * kappa, -log(1 - A) / C2)
  
  # Censoring
  if (lambda_c == 0) {
    T_obs <- T_i
    status <- rep(1L, n)
  } else if (lambda_c > 0) {
    C_i <- rexp(n, rate = lambda_c)
    T_obs <- pmin(T_i, C_i, max_FU)
    status <- as.integer(T_i <= C_i & T_i <= max_FU)
  } else {
    T_obs <- pmin(T_i, max_FU)
    status <- as.integer(T_i <= max_FU)
  }
  
  # Longitudinal data: expand each participant to their valid visit times
  long_data <- do.call(rbind, lapply(seq_len(n), function(i) {
    jit <- rnorm(length(visit), 0, 1)
    jit[visit == 0] <- 0
    vt <- pmax(visit + jit, 0)
    vt <- vt[vt < T_obs[i]]
    if (length(vt) == 0) return(NULL)
    y_true <- beta_0 + (beta_1 + b1[i]) * vt + beta_2 * arm[i] * vt + b0[i]
    data.frame(id = i, time = vt, Y_true = y_true,
               Y_obs = y_true + rnorm(length(vt), 0, sigma_e),
               arm = arm[i], T_obs = T_obs[i])
  }))
  
  surv_data <- data.frame(id = seq_len(n), arm = arm,
                          T_true = T_i, T_obs = T_obs, status = status,
                          b0 = b0, b1 = b1)
  
  trt_df <- data.frame(id = seq_len(n),randgrp = factor(arm, levels = c(0, 1), labels = c("Control", "Experimental")))
  long_data <- merge(long_data, trt_df, by = "id")
  surv_data <- merge(surv_data, trt_df, by = "id")
  
  list(longitudinal = long_data, survival = surv_data)
}





