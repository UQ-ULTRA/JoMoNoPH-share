# Shared generator for all five coefficient scenarios.
# simulation_parameters(scenario): integer 1:5 -> coefficient/design list.
# simulate_joint_dataset(): covariance D, coefficients, baseline parameters,
# visit schedule, seed, sample size, follow-up and censoring -> longitudinal
# and survival data frames. Time is in months; use aft_mode Weibull/loglogistic
# and link_type value. log_HR/alpha_PH are unused historical interface arguments.
# Random draw order is preserved: random effects, visit jitter, event uniform,
# optional censoring exponential, then observation errors for each patient.

simulation_parameters <- function(scenario) {
  coefficients <- data.frame(beta_2 = c(0, 0.04, -0.04, 0.04, -0.04), log_AF = c(
    0, 0.9,
    0.9, -0.9, -0.9
  ))
  stopifnot(length(scenario) == 1L, scenario %in% 1:5)
  sigma_0 <- 15
  sigma_1 <- 0.2
  rho_01 <- -0.1
  sigma_01 <- rho_01 * sigma_0 * sigma_1
  list(
    beta_0 = 73, beta_1 = -0.04, beta_2 = coefficients$beta_2[scenario], sigma_e = 12,
    D = matrix(c(sigma_0^2, sigma_01, sigma_01, sigma_1^2), 2, 2), log_HR = 99, log_AF = coefficients$log_AF[scenario],
    alpha_AFT = 0.012, visit = c(0, 1, seq(3, 92, 3))
  )
}

simulate_joint_dataset <- function(
  D, beta_0, beta_1, beta_2, sigma_e, log_HR, alpha_PH, log_AF,
  alpha_AFT, weibull_shape, weibull_scale, loglogistic_shape, loglogistic_scale, visit, seed,
  n_patients, max_FU, lambda_c, aft_mode, link_type
) {
  stopifnot(link_type == "value")
  if (!is.null(seed)) {
    set.seed(seed)
  }
  id <- 1:n_patients
  arm <- rep(0:1, each = n_patients / 2)
  b_mat <- mvrnorm(n_patients, mu = c(0, 0), Sigma = D)
  long_data <- data.frame()
  surv_data <- data.frame()
  m_func_raw <- function(t, arm_i, b0, b1) {
    beta_0 + beta_1 * t + beta_2 * (arm_i * t) + b0 + b1 * t
  }
  for (i in 1:n_patients) {
    arm_i <- arm[i]
    b0 <- b_mat[i, 1]
    b1 <- b_mat[i, 2]
    jitter_vis <- rnorm(length(visit), mean = 0, sd = 1)
    if (aft_mode == "Weibull") {
      v <- runif(1, min = 1e-06, max = 1)
      eta0 <- if (link_type == "value") {
        beta_0 + b0
      } else {
        0
      }
      C1 <- log_AF * arm_i + alpha_AFT * eta0
      C2 <- if (link_type %in% c("value", "slope")) {
        alpha_AFT * (beta_1 + beta_2 * arm_i + b1)
      } else {
        0
      }
      power_term <- (-log(v))^(1 / weibull_shape)
      A <- C2 * weibull_scale * exp(C1) * power_term
      if (abs(C2) < 1e-08) {
        T_i <- weibull_scale * exp(C1) * power_term
      } else {
        A <- min(A, 1 - 1e-08)
        T_i <- log(1 - A) / (-C2)
      }
    } else if (aft_mode == "loglogistic") {
      v <- runif(1, min = 1e-06, max = 1)
      eta0 <- if (link_type == "value") {
        beta_0 + b0
      } else {
        0
      }
      C1 <- log_AF * arm_i + alpha_AFT * eta0
      C2 <- if (link_type %in% c("value", "slope")) {
        alpha_AFT * (beta_1 + beta_2 * arm_i + b1)
      } else {
        0
      }
      shape <- loglogistic_shape
      scale <- loglogistic_scale
      inv_v_term <- (1 / v - 1)^(1 / shape)
      if (abs(C2) > 1e-08) {
        A <- C2 * scale * exp(C1) * inv_v_term
        A <- min(A, 1 - 1e-08)
        T_i <- -log(1 - A) / C2
      } else {
        T_i <- scale * exp(C1) * inv_v_term
      }
    } else {
      stop("Use aft_mode Weibull or loglogistic for these simulations.")
    }
    if (lambda_c == 0) {
      T_obs <- T_i
      status <- 1L
    } else if (!is.null(lambda_c) && is.finite(lambda_c) && lambda_c > 0) {
      C_i <- rexp(1, rate = lambda_c)
      T_obs <- min(T_i, C_i, max_FU)
      status <- as.integer(T_i <= C_i & T_i <= max_FU)
    } else {
      T_obs <- min(T_i, max_FU)
      status <- as.integer(T_i <= max_FU)
    }
    jittered_times <- visit + jitter_vis
    jittered_times[visit == 0] <- 0
    visit_times <- pmax(jittered_times, 0)
    valid_visits <- visit_times[visit_times < T_obs]
    visits <- valid_visits
    if (length(visits) > 0) {
      true_vals <- sapply(visits, function(t) m_func_raw(t, arm_i, b0, b1))
      obs_vals <- true_vals + rnorm(length(visits), 0, sigma_e)
      long_data <- rbind(long_data, data.frame(
        id = i, time = visits, Y_true = true_vals,
        Y_obs = obs_vals, arm = arm_i, T_obs = T_obs
      ))
    }
    surv_data <- rbind(surv_data, data.frame(
      id = i, arm = arm_i, T_true = T_i, T_obs = T_obs,
      status = status, b0 = b0, b1 = b1
    ))
  }
  trt_df <- data.frame(id = id, randgrp = factor(arm, levels = c(0, 1), labels = c(
    "Control",
    "Experimental"
  )))
  long_data <- merge(long_data, trt_df, by = "id")
  surv_data <- merge(surv_data, trt_df, by = "id")
  return(list(longitudinal = long_data, survival = surv_data))
}
