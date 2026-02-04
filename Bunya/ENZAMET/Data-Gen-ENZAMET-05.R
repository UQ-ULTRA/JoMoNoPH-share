# ------------------------------------------------------------
# Joint Longitudinal and Survival Data Simulation 
# ------------------------------------------------------------
# Simulates individual-level data for a two-arm randomized trial with:
#   - Linear mixed model for a continuous longitudinal outcome
#   - Weibull survival model (PH, AFT, or closed-form AFT)
#   - Time-to-event depends on current value or slope of the longitudinal process
#
# Longitudinal Model:
#   Y_ij = β₀ + β₁·time_ij + β₂·time_ij·treatment_i + b₀i + b₁i·time_ij + ε_ij
#
# Survival Model:
#   - Weibull baseline hazard (shape and scale parameters)
#   - Survival linked to Y(t) or dY(t)/dt
#   - Supports PH (hazard), AFT (log time), or closed-form AFT (inverse transform)
#
# Arguments:
#   D:               2×2 covariance matrix for random intercept and slope
#   beta_0:          Fixed intercept
#   beta_1:          Fixed slope for time
#   beta_2:          Fixed time-by-treatment interaction
#   sigma_e:         Residual standard deviation
#
#   weibull_shape:   Weibull shape parameter (ν)
#   weibull_scale:   Weibull scale parameter (λ)

#   Define baseline loglogistic shape and scale
#   loglogistic_shape
#   loglogistic_scale
#
#   aft_mode:        Survival model type: "PH", or "Weibull" or "loglogistic"
#   link_type:       Association type: "value" = Y(t), "slope" = dY(t)/dt
#
#   log_HR:          Log hazard ratio (treatment effect under PH)
#   alpha_PH:        Association with Y(t) or dY(t)/dt under PH (log hazard ratio per unit)
#   log_AF:          Log acceleration factor (treatment effect under AFT)
#   alpha_AFT:       Association under AFT (log acceleration factor per unit)
#
#   visit:           Vector of scheduled visit times
#   seed:            Random seed
#   n_patients:      Total sample size (must be even for 1:1 allocation)
#   max_FU:          Maximum follow-up time
#   lambda_c:        Exponential censoring rate (optional; calibrated if required):
#                    <0 for administrative censoring, controlled by max_FU
#                    =0 for no censoring
#                    >0 for random censoring + administrative censoring
#
# ------------------------------------------------------------



library(tidyverse)
library(nlme)
library(survival)
library(ggplot2)
library(MASS)
library(simsurv)

# Longitudinal model
# Fixed effects
beta_0 <- 73       # Baseline HRQoL at time 0
beta_1 <- -0.04          # Change over time in control group
beta_2 <-  -0.04          # Treatment-by-time interaction (treatment worse)
sigma_e <- 12         # Residual error SD


# Random effects
sigma_0 <- 15         # Random intercept SD
sigma_1 <- 0.20      # Random slope SD
rho_01 <- -0.10       # Negative intercept-slope correlation
sigma_01 <- rho_01 * sigma_0 * sigma_1
D <- matrix(c(sigma_0^2, sigma_01, sigma_01, sigma_1^2), 2, 2)


# Survival model

# Loglogistic parameters
loglogistic_shape <-  1.75   # Unimodal hazard (peak hazard followed by decline)
loglogistic_scale <- 12.0    # Survival centered  

# Baseline and treatment effects on log-hazard scale
log_baseline_hazard <- 99          # Baseline log-hazard (intercept)
log_HR <- 99                      # Log hazard ratio for treatment (HR = exp(log_HR))

# Weibull shape parameter (ν): controls time dependence of the hazard
weibull_shape <- 0.692                 # Shape > 1 ⇒ increasing hazard over time

# Convert baseline log-hazard to hazard rate (λ)
baseline_hazard_rate <- exp(log_baseline_hazard)  # Used as λ in simsurv

# Convert to Weibull scale parameter for use in rweibull/survreg
weibull_scale <- 13.823

# Compute acceleration factor (AF) under Weibull interpretation
log_AF <- -log_HR / weibull_shape

log_AF <- -0.90   # treatment worse

# Association parameter: PH and AFT linkage
# alpha_PH  <- 99                                # log hazard ratio per unit Y(t)
# alpha_AFT <- -alpha_PH / weibull_shape            # log acceleration factor per unit Y(t)
# alpha_AFT <- 0.012
alpha_AFT <- 0 # try no association 20260113

# Under a Weibull baseline, the AFT-scale association (log time) is derived by
# dividing the PH-scale association (log hazard) by –shape, preserving consistency across scales.


# Study design
visit <- c(0, 1, seq(3, 92, 3))
jitter_time <- c(0, .25, rep(.5, length(visit)-2))


n_patients <- 1100
max_FU <- 72  # maximum follow-up


# ============================================================
# Censoring-calibration for the joint simulator
# ============================================================
calibrate_lambda_c <- function(lambda_bounds = c(1e-4, 2),
                               target_cens = 0.50,
                               B = 2000,
                               seed = 123,
                               ...) {
  
  # wrapper to simulate B patients for a given lambda_c
  est_cens <- function(lambda_c) {
    dat <- simulate_joint_dataset(
      lambda_c = lambda_c,
      n_patients = B,
      seed = seed,
      ...
    )$survival
    
    mean(dat$status == 0)   # proportion censored
  }
  
  # compute censoring at bounds
  c_low  <- est_cens(lambda_bounds[1])
  c_high <- est_cens(lambda_bounds[2])
  
  message("Censoring at bounds:")
  message("  λ_low  = ", lambda_bounds[1], " → ", round(c_low, 3))
  message("  λ_high = ", lambda_bounds[2], " → ", round(c_high,3))
  
  # Check whether target censoring is bracketed
  bracketed <- (target_cens >= min(c_low, c_high)) &&
    (target_cens <= max(c_low, c_high))
  
  if (!bracketed) {
    warning("Target censoring not bracketed. Searching over enlarged grid.")
    
    # Expand over a grid of lambdas
    grid <- exp(seq(log(1e-4), log(5), length.out = 40))
    cens_grid <- sapply(grid, est_cens)
    
    best <- grid[ which.min(abs(cens_grid - target_cens)) ]
    
    return(best)
  }
  
  # Solve via uniroot
  f <- function(lambda) est_cens(lambda) - target_cens
  
  out <- uniroot(f, interval = lambda_bounds, tol = 1e-4)
  out$root
}

# lambda_c_new <- calibrate_lambda_c(
#   target_cens = 0.50,
#   B = 2000,
#   seed = 999,
#   D = D,
#   beta_0 = beta_0,
#   beta_1 = beta_1,
#   beta_2 = beta_2,
#   sigma_e = sigma_e,
#   log_HR = log_HR,
#   alpha_PH = alpha_PH,
#   log_AF = log_AF,
#   alpha_AFT = alpha_AFT,
#   weibull_shape = weibull_shape,
#   weibull_scale = weibull_scale,
#   loglogistic_shape = loglogistic_shape,
#   loglogistic_scale = loglogistic_scale,
#   visit = visit,
#   max_FU = max_FU,
#   aft_mode = "loglogistic",
#   link_type = "value"
# )

# lambda_c <- lambda_c_new


#--- Simulation Function ---
simulate_joint_dataset <- function(D, 
                                   beta_0, beta_1, beta_2, sigma_e,
                                   log_HR, alpha_PH, log_AF, alpha_AFT, 
                                   weibull_shape, weibull_scale,
                                   loglogistic_shape, loglogistic_scale,
                                   visit, seed, n_patients, max_FU, lambda_c,  
                                   aft_mode, link_type) {
  
  if (!is.null(seed)) set.seed(seed)
  id <- 1:n_patients
  arm <- rep(0:1, each = n_patients / 2)
  b_mat <- mvrnorm(n_patients, mu = c(0, 0), Sigma = D)
  
  design <- data.frame(id = id, arm = arm, b0 = b_mat[, 1], b1 = b_mat[, 2])
  
  long_data <- data.frame()
  surv_data <- data.frame()
  
  m_func_raw <- function(t, arm_i, b0, b1) {
    beta_0 + beta_1 * t + beta_2 * (arm_i * t) + b0 + b1 * t
  }
  
  lambda_wb <- weibull_scale^(-weibull_shape)
  rho_wb <- weibull_shape
  
  
  times_fun <- function(t, x, betas) {
    beta_0 <- betas["beta_0"]
    beta_1 <- betas["beta_1"]
    beta_2 <- betas["beta_2"]
    log_HR<- betas["log_HR"]
    alpha_PH <- betas["alpha_PH"]
    lambda_wb <- betas["lambda_wb"]
    rho_wb <- betas["rho_wb"]
    
    value_t <- beta_0 + beta_1 * t + beta_2 * (x$arm * t) + x$b0 + x$b1 * t
    slope_t <- beta_1 + beta_2 * x$arm + x$b1
    eta <- if (link_type == "slope") slope_t else value_t
    
    lambda_wb * rho_wb * t^(rho_wb - 1) * exp(log_HR * x$arm + alpha_PH * eta)
  }
  
  for (i in 1:n_patients) {
    arm_i <- arm[i]
    b0 <- b_mat[i, 1]
    b1 <- b_mat[i, 2]
    jitter_vis <- rnorm(length(visit), mean = 0, sd = 1)  
    
    if (aft_mode == "Weibull") {
      v <- runif(1, min = 1e-6, max = 1)
      
      # --- Step 1: Compute C1 ---
      eta0 <- if (link_type == "value") {
        beta_0 + b0 
      } else {
        0
      }
      C1 <- log_AF * arm_i + alpha_AFT * eta0
      
      # --- Step 2: Compute C2 ---
      C2 <- if (link_type %in% c("value", "slope")) {
        alpha_AFT * (beta_1 + beta_2 * arm_i + b1)
      } else {
        0
      }
      
      # --- Step 3: Sample survival time ---
      power_term <- (-log(v))^(1 / weibull_shape)
      A <- C2 * weibull_scale* exp(C1) * power_term
      
      if (abs(C2) < 1e-8) {
        T_i <- weibull_scale* exp(C1) * power_term
      } else {
        A <- min(A, 1 - 1e-8)
        T_i <- log(1 - A) / (-C2)
      }
    }
        else if (aft_mode == "loglogistic") {
      v <- runif(1, min = 1e-6, max = 1)  # avoid v = 0 or 1 for numerical stability
      
      # --- Step 1: Compute C1 ---
      eta0 <- if (link_type == "value") {
        beta_0 + b0
      } else {
        0
      }
      C1 <- log_AF * arm_i + alpha_AFT * eta0
      
      # --- Step 2: Compute C2 ---
      C2 <- if (link_type %in% c("value", "slope")) {
        alpha_AFT * (beta_1 + beta_2 * arm_i + b1)
      } else {
        0
      }
      
      # --- Step 3: Sample survival time ---
      shape <- loglogistic_shape
      scale <- loglogistic_scale
      inv_v_term <- (1 / v - 1)^(1 / shape)
      
      if (abs(C2) > 1e-8) {
        A <- C2 * scale * exp(C1) * inv_v_term
        A <- min(A, 1 - 1e-8)  # prevent log(0) or negative time
        T_i <- -log(1 - A) / C2
      } else {
        T_i <- scale * exp(C1) * inv_v_term
      }
    }
    

    else if (aft_mode == "PH") {
      single_pt <- data.frame(id = i, arm = arm_i, b0 = b0, b1 = b1)
      surv_sim <- simsurv(hazard = times_fun, 
                          x = single_pt,
                          maxt = max_FU, 
                          betas=c(
                            beta_0 = beta_0,
                            beta_1 = beta_1,
                            beta_2 = beta_2,
                            log_HR= log_HR,
                            alpha_PH = alpha_PH,
                            lambda_wb = lambda_wb,
                            rho_wb = rho_wb
                          ))
      T_i <- surv_sim$eventtime
      status <- as.integer(T_i <= max_FU)
      T_obs <- min(T_i, max_FU)
      
    } else if (aft_mode == "AFT") {
      m_func <- function(t) {
      beta_0 + beta_1 * t + beta_2 * (arm_i * t) + b0 + b1 * t
      }
      slope_func <- function(t) {
        beta_1 + beta_2 * arm_i + b1
      }
      eps_i <- log(rweibull(1, shape = weibull_shape, scale = weibull_scale))
      root_fun <- function(t) {
        if (t <= 0) return(1e6)
        eta <- if (link_type == "slope") slope_func(t) else m_func(t)
        log(t) - (eps_i + log_AF * arm_i + alpha_AFT * eta)
      }
      T_i <- tryCatch(uniroot(root_fun, interval = c(1e-8, 1e8))$root, error = function(e) NA)
      if (is.na(T_i)) next
    } else {
      stop("Invalid aft_mode")
    }
    
    # status <- as.integer(T_i <= max_FU)
    # T_obs <- min(T_i, max_FU)

    if (lambda_c == 0) {
      # No censoring
      T_obs <- T_i
      status <- 1L
      ## --- Independent exponential censoring + administrative censoring ---
    } else if (!is.null(lambda_c) && is.finite(lambda_c) && lambda_c > 0) {
      C_i   <- rexp(1, rate = lambda_c)
      T_obs <- min(T_i, C_i, max_FU)
      status <- as.integer(T_i <= C_i & T_i <= max_FU)
    } else {
      # no random censoring, only administrative
      T_obs <- min(T_i, max_FU)
      status <- as.integer(T_i <= max_FU)
    }
    
    # Define allowable window 
    jittered_times <- visit + jitter_vis
    jittered_times[visit == 0] <- 0  # Ensure baseline time is exactly 0
    visit_times <- pmax(jittered_times, 0)
    valid_visits <- visit_times[visit_times < T_obs]
    visits <- valid_visits
    
    if (length(visits) > 0) {
      true_vals <- sapply(visits, function(t) m_func_raw(t, arm_i, b0, b1))
      obs_vals <- true_vals + rnorm(length(visits), 0, sigma_e)
      long_data <- rbind(long_data, data.frame(id = i, time = visits,
                                               Y_true = true_vals, Y_obs = obs_vals,
                                               arm = arm_i, T_obs = T_obs))
    }
    
    surv_data <- rbind(surv_data, data.frame(id = i, arm = arm_i,
                                             T_true = T_i, T_obs = T_obs, status = status,
                                             b0 = b0, b1 = b1))
  }
  
  trt_df <- data.frame(id = id,
                       randgrp = factor(arm, levels = c(0, 1), labels = c("Control", "Experimental")))
  
  long_data <- merge(long_data, trt_df, by = "id")
  surv_data <- merge(surv_data, trt_df, by = "id")

  # === Return both versions ===
  return(list(
    longitudinal = long_data,
    survival = surv_data
  ))

}


#--- Simulation ---
# sim_data <- simulate_joint_dataset(
#   D = D,
#   beta_0 = beta_0,
#   beta_1 = beta_1,
#   beta_2 = beta_2,
#   sigma_e = sigma_e,
#   log_HR = log_HR,
#   log_AF = log_AF,
#   alpha_PH = alpha_PH,
#   alpha_AFT = alpha_AFT,
#   weibull_shape= weibull_shape,
#   weibull_scale= weibull_scale,
#   loglogistic_shape = loglogistic_shape,
#   loglogistic_scale = loglogistic_scale,
#   visit = visit,
#   seed = 301,
#   n_patients = n_patients,
#   max_FU = max_FU,
#   lambda_c = lambda_c,
#   aft_mode =  "loglogistic",  # "PH" or "Weibull" or "loglogistic"
#   link_type = "value" # "value" or "slope"
# )

