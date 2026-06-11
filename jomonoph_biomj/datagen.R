
# Load packages
library(MASS) # for mvrnorm

# Inverse CDF functions for survival time generation
ll_inv <- function(u, shape, scale) scale*((u / (1 - u))^(1 / shape))
wb_inv <- function(u, shape, scale) scale*(-log(1 - u))^(1 / shape)

# Primary simulation function
simulate_joint_dataset <- function(D = matrix(c(15^2, -0.10*15*0.20, -0.10*15*0.20, 0.20^2), 2, 2), 
                                   beta_0 = 73, beta_1 = -0.04, beta_2 = 0.00, sigma_e = 12,
                                   log_AF = 0.00, alpha_AFT = 0.012,
                                   weibull_shape = 0.90, weibull_scale = 38,
                                   loglogistic_shape = 1.20, loglogistic_scale = 23,
                                   visit = c(0, 1, seq(3, 92, 3)), 
                                   seed, n_patients = 1100, max_FU = 120, lambda_c = -1,
                                   aft_mode = "loglogistic", link_type = "value",
                                   ...) {  # absorb unused PH/Weibull placeholders
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialise parameters
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
    wb_inv(U, weibull_shape, weibull_scale)  # reuse shape/scale args
  }
  
  # Survival times: T = -log(1 - C2*exp(C1)*kappa) / C2 when C2 != 0
  if (aft_mode == "weibull") {
    power_term <- (-log(1 - U))^(1 / weibull_shape)
    A <- pmin(C2 * weibull_scale * exp(C1) * power_term, 1 - 1e-8)
    T_i <- ifelse(abs(C2) < 1e-8,
                  weibull_scale * exp(C1) * power_term,
                  log(1 - A) / (-C2))
  } else {
    A <- pmin(C2 * exp(C1) * kappa, 1 - 1e-8)
    T_i <- ifelse(abs(C2) < 1e-8, exp(C1) * kappa, -log(1 - A) / C2)
  }
  
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

# Code for initialising Stan parameter values, for stability purposes in MCMC computation
make_init <- function(chains = 4, lmm_fit, surv_fit) {
  beta_long_init <- as.numeric(fixef(lmm_fit))
  sigma_long_init <- summary(lmm_fit)$sigma
  
  vc <- nlme::VarCorr(lmm_fit)
  sd_b_init <- sqrt(as.numeric(vc[c("(Intercept)", "time"), "Variance"]))
  sd_b_init <- pmax(sd_b_init, 1e-3)
  
  beta_surv_init <- as.numeric(coef(surv_fit)["arm"])
  if (!is.finite(beta_surv_init)) beta_surv_init <- 0
  
  lapply(seq_len(chains), function(chain_id) {
    rho <- -0.1
    R <- matrix(c(1, rho, rho, 1), 2, 2)
    L_init <- t(chol(R))
    
    list(
      beta_long = beta_long_init + rnorm(stan_data$K_long, 0, 0.05),
      sigma_long = abs(rnorm(1, sigma_long_init, 0.25)),
      sd_1_long = pmax(sd_b_init + rnorm(2, 0, 0.02), 1e-3),
      z_1_long = matrix(rnorm(2 * stan_data$N_1_long, 0, 0.2), nrow = 2),
      L_1_long = L_init,
      beta_surv = rnorm(stan_data$q, beta_surv_init, 0.05),
      gamma = pmax(rnorm(stan_data$m, 0.1, 0.01), 1e-4),
      alpha_tilde = rnorm(1, 0, 0.05)
    )
  })
}






