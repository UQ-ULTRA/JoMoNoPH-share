# ------------------------------------------------------------
# Joint Longitudinal and Survival Data Simulation 
# ------------------------------------------------------------
# Simulates individual-level data for a two-arm randomized trial with:
#   - Linear mixed model for a continuous longitudinal outcome
#   - Weibull survival model (PH, AFT, or closed-form AFT)
#   - Time-to-event depends on current value or slope of the longitudinal process
#
# DESIGN / USAGE NOTES:
#   - This script is intentionally written to be flexible across PH and AFT
#     formulations for methodological comparison and future extensions.
#   - In practice, simulations for this project use AFT formulations
#     ("Weibull" or "loglogistic") only.
#   - The PH formulation is retained for completeness and is not used
#     in production runs; PH-specific parameters may therefore appear
#     as placeholders.
#
# Longitudinal Model:
#   Y_ij = β₀ + β₁·time_ij + β₂·time_ij·treatment_i + b₀i + b₁i·time_ij + ε_ij
#
# Survival Model:
#   - Weibull or log-logistic baseline hazard
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
#
#   loglogistic_shape: Log-logistic shape parameter
#   loglogistic_scale: Log-logistic scale parameter
#
#   aft_mode:        Survival model type: "PH", "Weibull", or "loglogistic"
#                    (PH retained for completeness; not used in practice)
#   link_type:       Association type:
#                    "value" = current Y(t)
#                    "slope" = dY(t)/dt
#
#   log_HR:          Log hazard ratio (PH; placeholder unless PH used)
#   alpha_PH:        Longitudinal association under PH (placeholder unless PH used)
#   log_AF:          Log acceleration factor (AFT treatment effect)
#   alpha_AFT:       Longitudinal association under AFT
#
#   visit:           Vector of scheduled visit times (with jittering applied)
#   seed:            Random seed
#   n_patients:      Total sample size (must be even for 1:1 allocation)
#   max_FU:          Maximum follow-up time
#   lambda_c:        Exponential censoring rate (optional; calibrated if required)
#
# ------------------------------------------------------------




library(tidyverse)
library(nlme)
library(survival)
library(ggplot2)
library(MASS)
library(simsurv)

# ============================================================
# Censoring-calibration for the joint simulator
# ============================================================
# calibrate_lambda_c <- function(lambda_bounds = c(1e-4, 2),
#                                target_cens = 0.50,
#                                B = 2000,
#                                seed = 123,
#                                ...) {
#   
#   # wrapper to simulate B patients for a given lambda_c
#   est_cens <- function(lambda_c) {
#     dat <- simulate_joint_dataset(
#       lambda_c = lambda_c,
#       n_patients = B,
#       seed = seed,
#       ...
#     )$survival
#     
#     mean(dat$status == 0)   # proportion censored
#   }
#   
#   # compute censoring at bounds
#   c_low  <- est_cens(lambda_bounds[1])
#   c_high <- est_cens(lambda_bounds[2])
#   
#   message("Censoring at bounds:")
#   message("  λ_low  = ", lambda_bounds[1], " → ", round(c_low, 3))
#   message("  λ_high = ", lambda_bounds[2], " → ", round(c_high,3))
#   
#   # Check whether target censoring is bracketed
#   bracketed <- (target_cens >= min(c_low, c_high)) &&
#     (target_cens <= max(c_low, c_high))
#   
#   if (!bracketed) {
#     warning("Target censoring not bracketed. Searching over enlarged grid.")
#     
#     # Expand over a grid of lambdas
#     grid <- exp(seq(log(1e-4), log(5), length.out = 40))
#     cens_grid <- sapply(grid, est_cens)
#     
#     best <- grid[ which.min(abs(cens_grid - target_cens)) ]
#     
#     return(best)
#   }
#   
#   # Solve via uniroot
#   f <- function(lambda) est_cens(lambda) - target_cens
#   
#   out <- uniroot(f, interval = lambda_bounds, tol = 1e-4)
#   out$root
# }

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


#### Simulation Function (load this) ####
simulate_joint_dataset <- function(D, 
                                   beta_0, beta_1, beta_2, sigma_e,
                                   log_HR, alpha_PH, log_AF, alpha_AFT, 
                                   weibull_shape, weibull_scale,
                                   loglogistic_shape, loglogistic_scale,
                                   visit, seed, n_patients, max_FU, lambda_c,  
                                   aft_mode, link_type) {
  # browser()
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
    
    ## --- Independent exponential censoring + administrative censoring ---
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












# Example usage:

# # Longitudinal model
# # Fixed effects
# beta_0 <- 73       # Baseline HRQoL at time 0
# beta_1 <- -0.04          # Change over time in control group
# beta_2 <-  0.00          # Treatment-by-time interaction (no treatment effect)
# sigma_e <- 12         # Residual error SD
# 
# 
# # Random effects
# sigma_0 <- 15         # Random intercept SD
# sigma_1 <- 0.20      # Random slope SD
# rho_01 <- -0.10       # Negative intercept-slope correlation
# sigma_01 <- rho_01 * sigma_0 * sigma_1
# D <- matrix(c(sigma_0^2, sigma_01, sigma_01, sigma_1^2), 2, 2)
# 
# 
# # Survival model
# 
# # Loglogistic parameters
# # Old Bernstein Polynomial approximation:
# # loglogistic_shape <- 1.75
# # loglogistic_scale <- 12.0
# loglogistic_shape <- 1.20
# loglogistic_scale <- 23.0
# 
# # Baseline and treatment effects on log-hazard scale
# log_baseline_hazard <- 99          # Baseline log-hazard (intercept)
# log_HR <- 99                      # Log hazard ratio for treatment (HR = exp(log_HR))
# 
# # Weibull shape parameter (ν): controls time dependence of the hazard
# # Old Bernstein Polynomial approximation:
# # weibull_shape <- 0.692
# weibull_shape <- 0.70
# 
# # Convert baseline log-hazard to hazard rate (λ)
# baseline_hazard_rate <- exp(log_baseline_hazard)  # Used as λ in simsurv
# 
# # Convert to Weibull scale parameter for use in rweibull/survreg
# # Old Bernstein Polynomial approximation:
# # weibull_scale <- 13.823
# weibull_scale <- 38.0
# 
# # Compute acceleration factor (AF) under Weibull interpretation
# log_AF <- -log_HR / weibull_shape
# 
# log_AF <- 0.00   
# 
# # Association parameter: PH and AFT linkage
# # alpha_PH  <- 99                                # log hazard ratio per unit Y(t)
# # alpha_AFT <- -alpha_PH / weibull_shape            # log acceleration factor per unit Y(t)
# # alpha_AFT <- 0.012
# alpha_AFT <- 0 # try no association 20260113
# 
# # Under a Weibull baseline, the AFT-scale association (log time) is derived by
# # dividing the PH-scale association (log hazard) by –shape, preserving consistency across scales.
# 
# 
# # Study design
# visit <- c(0, 1, seq(3, 92, 3))
# jitter_time <- c(0, .25, rep(.5, length(visit)-2))
# 
# 
# n_patients <- 1100
# max_FU <- 72  # maximum follow-up



###### above copied from Data-Gen-ENZAMET-01.R ######

###### single test of simulate_joint_dataset() function ######
# single test for joint model data generation
sim_data <- simulate_joint_dataset(
  D = matrix(c(15^2, -0.3, -0.3, 0.20^2), 2, 2),
  beta_0 = 73,
  beta_1 = -0.04,
  beta_2 = 0.00,
  sigma_e = 12,
  log_HR = 99,
  log_AF = 0.00,
  alpha_PH = 99,
  alpha_AFT = 0,
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  seed = 32,
  n_patients = 1100,
  max_FU = 72,
  lambda_c = 0.05331733, # 0.023 for around 50% censoring under Weibull, 0.013 for around 50% censoring loglogistic
  aft_mode =  "loglogistic",  # "PH" or "Weibull" or "loglogistic"
  link_type = "value" # "value" or "slope"
)

# single test for aft only data generation
sim_data <- simulate_joint_dataset(
  D = matrix(c(0, 0, 0, 0), 2, 2),
  beta_0 = 0,
  beta_1 = 0,
  beta_2 = 0,
  sigma_e = 0,
  log_HR = 0,
  log_AF = 0,
  alpha_PH = 0,
  alpha_AFT = 0,
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  seed = 32,
  n_patients = 1100,
  max_FU = 72,
  lambda_c = 0, # 0 for no censoring, <0 for administrative censoring only by max_FU
  aft_mode =  "loglogistic",  # "PH" or "Weibull" or "loglogistic"
  link_type = "value" # "value" or "slope"
)
mean(sim_data$survival$status)




############## average censoring function (slow version, ignore and skip) ###############
# create a function to use a number of repeated simulations (randomness controlled by seeds) to calculate the average censoring proportion
avg_censoring <- function(lambda_c, seeds, n_patients, n_cores = 4, ...) {
  censoring_prop <- numeric(length(seeds))
  
  for (i in seq_along(seeds)) {
    sim_data <- simulate_joint_dataset(
      seed = seeds[i],
      n_patients = n_patients,
      lambda_c = lambda_c,
      ...
    )
    
    # status == 0 means censored
    censoring_prop[i] <- mean(sim_data$survival$status == 0)
  }
  
  list(
    censoring_proportion = censoring_prop,
    avg_censoring = mean(censoring_prop)
  )
}

###### Example simulation to test avg_censoring() function (ignore and skip) ######
avg_censoring(
  lambda_c = 0.05331733,
  # generate 10 random integers as seeds
  seeds = sample(1:1e8, 10, replace = FALSE),
  n_patients = 1100,
  n_cores = 8, # adjust based on your machine
  # other parameters
  D = matrix(c(15^2, -0.3, -0.3, 0.20^2), 2, 2),
  beta_0 = 73,
  beta_1 = -0.04,
  beta_2 = 0.00,
  sigma_e = 12,
  log_HR = 99,
  log_AF = 0.00,
  alpha_PH = 99,
  alpha_AFT = 0,
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  max_FU = 72,
  aft_mode =  "loglogistic",
  link_type = "value"
)

############# average censoring function parallel version (on Macbook Pro 14-inch M2pro, load this) ###############
library(parallel)

avg_censoring_parallel <- function(lambda_c,
                                   seeds,
                                   n_patients,
                                   n_cores = 4,
                                   verbose = FALSE,
                                   ...) {
  
  censoring_prop <- mclapply(
    seeds,
    function(seed_i) {
      tryCatch(
        {
          sim_data <- simulate_joint_dataset(
            seed = seed_i,
            n_patients = n_patients,
            lambda_c = lambda_c,
            ...
          )
          
          mean(sim_data$survival$status == 0)
        },
        error = function(e) {
          if (verbose) {
            message(sprintf("Seed %d failed: %s", seed_i, e$message))
          }
          NA_real_
        }
      )
    },
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  
  censoring_prop <- unlist(censoring_prop)
  
  list(
    censoring_proportion = censoring_prop,
    avg_censoring = mean(censoring_prop, na.rm = TRUE),
    n_failed = sum(is.na(censoring_prop))
  )
}
###### Example simulation to test avg_censoring_parallel() function (ignore and skip) ###### 
avg_censoring_parallel(
  lambda_c = 0.05331733,
  seeds = sample(1:1e8, 10, replace = FALSE),
  n_patients = 1100,
  n_cores = 8,
  # other parameters
  D = matrix(c(15^2, -0.3, -0.3, 0.20^2), 2, 2),
  beta_0 = 73,
  beta_1 = -0.04,
  beta_2 = 0.00,
  sigma_e = 12,
  log_HR = 99,
  log_AF = 0.00,
  alpha_PH = 99,
  alpha_AFT = 0,
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  max_FU = 72,
  aft_mode =  "loglogistic",
  link_type = "value"
)

############# a function to find lambda_c for target censoring proportion (load this) ###############

# a function to calibrate lambda_c for target censoring proportion using bisection method
calibrate_lambda_c <- function(target_cens,
                               seeds,
                               n_patients,
                               n_cores = 4,
                               lambda_bounds = c(1e-4, 1),
                               tol = 1e-3,
                               max_iter = 20,
                               verbose = TRUE,
                               ...) {

  eval_avg <- function(lambda) {
    res <- avg_censoring_parallel(
      lambda_c = lambda,
      seeds = seeds,
      n_patients = n_patients,
      n_cores = n_cores,
      ...
    )

    if (res$n_failed > 0 && verbose) {
      message(sprintf("  %d/%d seeds failed", res$n_failed, length(seeds)))
    }

    res$avg_censoring
  }
  
  lambda_lo <- lambda_bounds[1]
  lambda_hi <- lambda_bounds[2]
  
  c_lo <- eval_avg(lambda_lo)
  c_hi <- eval_avg(lambda_hi)

  
  if (!is.finite(c_lo) || !is.finite(c_hi)) {
    stop("Non-finite censoring estimates at initial bounds. Check simulation stability.")
  }
  
  if (verbose) {
    message("Initial bracket:")
    message(sprintf("  lambda = %.6g, avg censoring = %.6f", lambda_lo, c_lo))
    message(sprintf("  lambda = %.6g, avg censoring = %.6f", lambda_hi, c_hi))
  }
  
  if (!(target_cens >= c_lo && target_cens <= c_hi)) {
    stop("Target censoring proportion not bracketed by lambda_bounds.")
  }
  
  for (iter in seq_len(max_iter)) {
    
    lambda_mid <- 0.5 * (lambda_lo + lambda_hi)
    c_mid <- eval_avg(lambda_mid)
    
    if (!is.finite(c_mid)) next
    
    if (verbose) {
      message(sprintf(
        "Iter %02d: lambda = %.6f, avg censoring = %.6f",
        iter, lambda_mid, c_mid
      ))
    }
    
    if (abs(c_mid - target_cens) < tol) {
      message("Converged.")
      return(list(
        lambda_c = lambda_mid,
        avg_censoring = c_mid,
        iterations = iter,
        converged = TRUE
      ))
    }
    
    if (c_mid < target_cens) {
      lambda_lo <- lambda_mid
    } else {
      lambda_hi <- lambda_mid
    }
  }
  
  list(
    lambda_c = lambda_mid,
    avg_censoring = c_mid,
    iterations = max_iter,
    converged = FALSE
  )
}


###### Example usage of calibrate_lambda_c() function using 100 seeds ######
# use the same parameters as what used in the actual simulation to get accurate lambda_c
res <- calibrate_lambda_c(
  target_cens = 0.10,
  seeds = 1:100, # use 1:100 seeds
  n_patients = 1100,
  n_cores = 8, # adjust based on your machine
  lambda_bounds = c(1e-4, 1),  # initial bounds for lambda_c
  tol = 0.001, # tolerance for difference between target censoring and average censoring proportions
  D = matrix(c(15^2, -0.10, -0.10, 0.20^2), 2, 2),
  beta_0 = 73, beta_1 = -0.04,
  beta_2 = 0.00, # scenario 1 0.00, scenario 2 0.04, scenario 3 -0.04, scenario 4 0.04, scenario 5 -0.04
  sigma_e = 12,
  log_HR = 99,
  alpha_PH = 99,
  log_AF = 0.00, # scenario 1 0.00, scenario 2 0.90, scenario 3 0.90, scenario 4 -0.90, scenario 5 -0.90
  alpha_AFT = 0, # try no association
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  max_FU = 72,
  aft_mode = "loglogistic", # loglogistic or Weibull
  link_type = "value"
)

res$lambda_c
res$avg_censoring

# calibration results (log-logistic, 50% censoring, 1:100 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.05331733
# scenario 2 (beta_2 = 0.04, log_AF = 0.90): lambda_c = 0.0332998
# scenario 3 (beta_2 = -0.04, log_AF = 0.90): lambda_c = 0.0332998
# scenario 4 (beta_2 = 0.04, log_AF = -0.90): lambda_c = 0.08236716
# scenario 5 (beta_2 = -0.04, log_AF = -0.90): lambda_c = 0.08236716


# calibration results (Weibull, 50% censoring, 1:100 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.08114658
# scenario 2 (beta_2 = 0.04, log_AF = 0.90): lambda_c = 0.05234087
# scenario 3 (beta_2 = -0.04, log_AF = 0.90): lambda_c = 0.05234087
# scenario 4 (beta_2 = 0.04, log_AF = -0.90): lambda_c = 0.1289934
# scenario 5 (beta_2 = -0.04, log_AF = -0.90): lambda_c = 0.1289934





# seeds = 1e6 * as.integer(as.numeric(as.factor(paste(n_patients, lambda_c, aft_mode, k_bases, scenario, max_FU, sep = "_")))) + 1:
# seeds = 1e6 * as.integer(as.numeric(as.factor(paste(c(1100), c(0.05322266, 0.03295898, 0.08300781), c("Weibull", "loglogistic"), c(5), c(1,2,4), c(72), sep = "_"))))


##### find lambda_c's (control censoring proportion) for "aft only" data generation using 500 seeds (run this) #####
res <- calibrate_lambda_c(
  target_cens = 0.50,
  # seeds = 1:100, # use 1:100 seeds
  seeds = 1:500, # use 1:500 seeds
  n_patients = 1100,
  n_cores = 8, # adjust based on your machine
  lambda_bounds = c(0, 1),  # initial bounds for lambda_c
  tol = 0.001, # tolerance for difference between target censoring and average censoring proportions
  D = matrix(c(0^2, 0, 0, 0^2), 2, 2),
  beta_0 = 0, beta_1 = 0,
  beta_2 = 0.00, # scenario 1 0.00, scenario 2 0.04, scenario 3 -0.04, scenario 4 0.04, scenario 5 -0.04
  sigma_e = 0,
  log_HR = 0,
  alpha_PH = 0,
  log_AF = -0.90, # scenario 1 0.00, scenario 2 0.90, scenario 3 0.90, scenario 4 -0.90, scenario 5 -0.90
  alpha_AFT = 0, # try no association, aft only
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  max_FU = 72,
  aft_mode = "Weibull", # loglogistic or Weibull
  link_type = "value"
)

res$lambda_c
res$avg_censoring

# calibration results (Weibull, 50% censoring, 1:100 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the old Bernstein Polynomial approximation: weibull_shape = 0.692, weibull_scale = 13.823,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.08203125
# scenario 2 (beta_2 = 0.04, log_AF = 0.90): lambda_c = 0.05175781
# scenario 3 (beta_2 = -0.04, log_AF = 0.90): lambda_c = 0.05175781
# scenario 4 (beta_2 = 0.04, log_AF = -0.90): lambda_c = 0.1298828
# scenario 5 (beta_2 = -0.04, log_AF = -0.90): lambda_c = 0.1298828

# calibration results (log-logistic, 50% censoring, 1:100 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the old Bernstein Polynomial approximation: loglogistic_shape = 1.75, loglogistic_scale = 12.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.05322266
# scenario 2 (beta_2 = 0.04, log_AF = 0.90): lambda_c = 0.03295898; (1e6*1+1):(1e6*1+100) seeds produce lambda_c = 0.03320312
# scenario 3 (beta_2 = -0.04, log_AF = 0.90): lambda_c = 0.03295898
# scenario 4 (beta_2 = 0.04, log_AF = -0.90): lambda_c = 0.08300781
# scenario 5 (beta_2 = -0.04, log_AF = -0.90): lambda_c = 0.08300781



# calibration results (Weibull, 50% censoring, 1:500 seeds, n=1100, alpha_AFT=0, all beta=0, max_FU=72):
# for the old Bernstein Polynomial approximation: weibull_shape = 0.692, weibull_scale = 13.823,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.08154297, avg_censoring = 0.5001691
# scenario 2 (beta_2 = 0.04, log_AF = 0.90): lambda_c = 0.05224609, avg_censoring = 0.5003527
# scenario 3 (beta_2 = -0.04, log_AF = 0.90): lambda_c = 0.05224609, avg_censoring = 0.5003527
# scenario 4 (beta_2 = 0.04, log_AF = -0.90): lambda_c = 0.1289062, avg_censoring = 0.4998873
# scenario 5 (beta_2 = -0.04, log_AF = -0.90): lambda_c = 0.1289062, avg_censoring = 0.4998873

# 20260427
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 0.70, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.02832031, avg_censoring = 0.5005309
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.015625, avg_censoring = 0.5001
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.04638672, avg_censoring = 0.4999345

# 20260506
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 0.40, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.04052734, avg_censoring = 0.5001927
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.02441406, avg_censoring = 0.4993745
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.06445312, avg_censoring = 0.4990636

# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 0.90, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.02563477, avg_censoring = 0.4992855
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.0135498, avg_censoring = 0.4999655
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.04223633, avg_censoring = 0.4997545

# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 1.30, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.02368164, avg_censoring = 0.5006127
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.01171875, avg_censoring = 0.5009891
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.03808594, avg_censoring = 0.4994727

# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 2.10, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.02233887, avg_censoring = 0.50036
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.009887695, avg_censoring = 0.5008618
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.03466797, avg_censoring = 0.5002764

# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 4.20, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.02062988, avg_censoring = 0.5003127
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.006469727, avg_censoring = 0.5009564
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.03149414, avg_censoring = 0.5007145


# calibration results (log-logistic, 50% censoring, 1:500 seeds, n=1100, alpha_AFT=0, all beta=0, max_FU=72): 
# for the old Bernstein Polynomial approximation: loglogistic_shape = 1.75, loglogistic_scale = 12.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.05322266, avg_censoring = 0.5001255
# scenario 2 (beta_2 = 0.04, log_AF = 0.90): lambda_c = 0.03320312, avg_censoring = 0.5007873
# scenario 3 (beta_2 = -0.04, log_AF = 0.90): lambda_c = 0.03320312, avg_censoring = 0.5007873
# scenario 4 (beta_2 = 0.04, log_AF = -0.90): lambda_c = 0.08251953, avg_censoring = 0.4998182
# scenario 5 (beta_2 = -0.04, log_AF = -0.90): lambda_c = 0.08251953, avg_censoring = 0.4998182

# 20260427
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): loglogistic_shape = 1.20, loglogistic_scale = 23.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = 0.02587891, avg_censoring = 0.5007
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = 0.01391602, avg_censoring = 0.5005055
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = 0.04150391, avg_censoring = 0.4992745



##### check administrative censoring only "aft-only" using 500 seeds (run this) #####
# lambda_c < 0 triggers administrative censoring only:
#   T_obs = min(T_i, max_FU), status = I(T_i <= max_FU)
admin_censoring <- avg_censoring_parallel(
  lambda_c = -1,
  seeds = 1:500,
  n_patients = 1100,
  n_cores = 8,
  D = matrix(c(0^2, 0, 0, 0^2), 2, 2),
  beta_0 = 0, beta_1 = 0,
  beta_2 = 0.00, # scenario 1 0.00, scenario 2 0.04, scenario 3 -0.04, scenario 4 0.04, scenario 5 -0.04
  sigma_e = 0,
  log_HR = 0,
  alpha_PH = 0,
  log_AF = -0.90, # scenario 1 0.00, scenario 2 0.90, scenario 3 0.90, scenario 4 -0.90, scenario 5 -0.90
  alpha_AFT = 0, # try no association, aft only
  # Old Bernstein Polynomial approximation:
  # weibull_shape = 0.692,
  # weibull_scale = 13.823,
  # loglogistic_shape = 1.75,
  # loglogistic_scale = 12.0,
  # Previous amended BP setting:
  # weibull_shape = 1.30,
  # weibull_shape = 2.10,
  # weibull_shape = 0.90,
  # weibull_shape = 0.40,
  # weibull_shape = 4.20,
  weibull_shape = 0.70,
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  visit = c(0, 1, seq(3, 92, 3)),
  max_FU = 72,
  aft_mode = "Weibull", # loglogistic or Weibull
  link_type = "value"
)

admin_censoring$avg_censoring
admin_censoring$n_failed

# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 0.70, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.2090255
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.32196
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.1311945

# 20260506
# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 0.90, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.1689982
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.3109982
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.09368182

# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 1.30, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.1001582
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.2949364
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.05034182

# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 2.10, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.02129636
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.2904255
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.01054182

# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 4.20, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.35774
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0


# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=84): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 2.10, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.004929091
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.2269218
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.002345455

# calibration results (Weibull, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=96):
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 0.70, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.1471982
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.2542164
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.08718545

# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 1.30, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.03498909
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.1950909
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.01736

# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): weibull_shape = 2.10, weibull_scale = 38.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.0008709091
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.1741564
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.0004345455



# calibration results (loglogistic, administrative censoring only, 1:500 seeds, n=1100, alpha_AFT=0, max_FU=72): 
# for the amended BP approximation (b2 = exp(b2) ./ tau_aft, B2 = exp(B2)): loglogistic_shape = 1.20, loglogistic_scale = 23.0,
# scenario 1 (beta_2 = 0.00, log_AF = 0.00): lambda_c = -1 (any negative value), avg_censoring = 0.2023836
# scenario 2 and 3 (beta_2 = 0.00, log_AF = 0.90): lambda_c = -1 (any negative value), avg_censoring = 0.3155127
# scenario 4 and 5 (beta_2 = 0.00, log_AF = -0.90): lambda_c = -1 (any negative value), avg_censoring = 0.1409127



##### find lambda_c's (control censoring proportion) for sAFT-JM data generation using Table 2 scenarios (load all four functions here) #####
# This section uses the rounded ENZAMET-inspired joint-model parameters in Table 2.
# In the data generation code, gamma in Table 2 is passed as log_AF.
# alpha_AFT is the longitudinal association parameter in the AFT event-time model.

jm_scenarios_table2 <- tibble::tribble(
  ~scenario, ~beta_1, ~beta_2, ~beta_1_plus_beta_2, ~log_AF,
  1L,        -0.04,   0.00,    -0.04,               0.00,
  2L,        -0.04,   0.04,     0.00,               0.90,
  3L,        -0.04,  -0.04,    -0.08,               0.90,
  4L,        -0.04,   0.04,     0.00,              -0.90,
  5L,        -0.04,  -0.04,    -0.08,              -0.90
)

calibrate_lambda_c_jm <- function(scenario_id,
                                  target_cens = 0.50,
                                  seeds = 1:500,
                                  n_patients = 1100,
                                  n_cores = 8,
                                  lambda_bounds = c(0, 1),
                                  tol = 0.001,
                                  max_iter = 20,
                                  max_FU = 72,
                                  aft_mode = "Weibull",
                                  weibull_shape = 0.70,
                                  weibull_scale = 38.0,
                                  loglogistic_shape = 1.20,
                                  loglogistic_scale = 23.0,
                                  alpha_AFT = 0.012,
                                  verbose = TRUE) {

  scenario_row <- jm_scenarios_table2[jm_scenarios_table2$scenario == scenario_id, ]
  if (nrow(scenario_row) != 1) {
    stop("scenario_id must be one of: ", paste(jm_scenarios_table2$scenario, collapse = ", "))
  }

  sigma_0 <- 15
  sigma_1 <- 0.20
  rho_01 <- -0.10
  sigma_01 <- rho_01 * sigma_0 * sigma_1
  D_jm <- matrix(c(sigma_0^2, sigma_01, sigma_01, sigma_1^2), 2, 2)

  res <- calibrate_lambda_c(
    target_cens = target_cens,
    seeds = seeds,
    n_patients = n_patients,
    n_cores = n_cores,
    lambda_bounds = lambda_bounds,
    tol = tol,
    max_iter = max_iter,
    verbose = verbose,
    D = D_jm,
    beta_0 = 73,
    beta_1 = scenario_row$beta_1,
    beta_2 = scenario_row$beta_2,
    sigma_e = 12,
    log_HR = 99,
    alpha_PH = 99,
    log_AF = scenario_row$log_AF,
    alpha_AFT = alpha_AFT,
    weibull_shape = weibull_shape,
    weibull_scale = weibull_scale,
    loglogistic_shape = loglogistic_shape,
    loglogistic_scale = loglogistic_scale,
    visit = c(0, 1, seq(3, 92, 3)),
    max_FU = max_FU,
    aft_mode = aft_mode,
    link_type = "value"
  )

  dplyr::bind_cols(
    scenario_row,
    tibble::tibble(
      aft_mode = aft_mode,
      weibull_shape = weibull_shape,
      weibull_scale = weibull_scale,
      loglogistic_shape = loglogistic_shape,
      loglogistic_scale = loglogistic_scale,
      alpha_AFT = alpha_AFT,
      max_FU = max_FU,
      target_cens = target_cens,
      lambda_c = res$lambda_c,
      avg_censoring = res$avg_censoring,
      iterations = res$iterations,
      converged = res$converged
    )
  )
}

calibrate_lambda_c_jm_all_scenarios <- function(scenarios = jm_scenarios_table2$scenario,
                                                ...) {
  dplyr::bind_rows(lapply(scenarios, function(scenario_id) {
    calibrate_lambda_c_jm(scenario_id = scenario_id, ...)
  }))
}

avg_admin_censoring_jm <- function(scenario_id,
                                   seeds = 1:500,
                                   n_patients = 1100,
                                   n_cores = 8,
                                   max_FU = 72,
                                   aft_mode = "Weibull",
                                   weibull_shape = 0.70,
                                   weibull_scale = 38.0,
                                   loglogistic_shape = 1.20,
                                   loglogistic_scale = 23.0,
                                   alpha_AFT = 0.012,
                                   verbose = TRUE) {

  scenario_row <- jm_scenarios_table2[jm_scenarios_table2$scenario == scenario_id, ]
  if (nrow(scenario_row) != 1) {
    stop("scenario_id must be one of: ", paste(jm_scenarios_table2$scenario, collapse = ", "))
  }

  sigma_0 <- 15
  sigma_1 <- 0.20
  rho_01 <- -0.10
  sigma_01 <- rho_01 * sigma_0 * sigma_1
  D_jm <- matrix(c(sigma_0^2, sigma_01, sigma_01, sigma_1^2), 2, 2)

  res <- avg_censoring_parallel(
    lambda_c = -1,
    seeds = seeds,
    n_patients = n_patients,
    n_cores = n_cores,
    verbose = verbose,
    D = D_jm,
    beta_0 = 73,
    beta_1 = scenario_row$beta_1,
    beta_2 = scenario_row$beta_2,
    sigma_e = 12,
    log_HR = 99,
    alpha_PH = 99,
    log_AF = scenario_row$log_AF,
    alpha_AFT = alpha_AFT,
    weibull_shape = weibull_shape,
    weibull_scale = weibull_scale,
    loglogistic_shape = loglogistic_shape,
    loglogistic_scale = loglogistic_scale,
    visit = c(0, 1, seq(3, 92, 3)),
    max_FU = max_FU,
    aft_mode = aft_mode,
    link_type = "value"
  )

  dplyr::bind_cols(
    scenario_row,
    tibble::tibble(
      aft_mode = aft_mode,
      weibull_shape = weibull_shape,
      weibull_scale = weibull_scale,
      loglogistic_shape = loglogistic_shape,
      loglogistic_scale = loglogistic_scale,
      alpha_AFT = alpha_AFT,
      max_FU = max_FU,
      lambda_c = -1,
      avg_censoring = res$avg_censoring,
      n_failed = res$n_failed
    )
  )
}

avg_admin_censoring_jm_all_scenarios <- function(scenarios = jm_scenarios_table2$scenario,
                                                 ...) {
  dplyr::bind_rows(lapply(scenarios, function(scenario_id) {
    avg_admin_censoring_jm(scenario_id = scenario_id, ...)
  }))
}

###### Example usage for one Table 2 sAFT-JM scenario (run this) ######
# jm_res <- calibrate_lambda_c_jm(
#   scenario_id = 1,       # Table 2 scenario: 1, 2, 3, 4, or 5
#   target_cens = 0.50,
#   seeds = 1:500,
#   n_patients = 1100,
#   n_cores = 8,
#   lambda_bounds = c(0, 1),
#   tol = 0.001,
#   max_FU = 72,
#   aft_mode = "Weibull",  # "Weibull" or "loglogistic"
#   # Candidate Weibull baselines:
#   weibull_shape = 0.70,
#   # weibull_shape = 0.90,
#   # weibull_shape = 1.30,
#   # weibull_shape = 2.10,
#   weibull_scale = 38.0,
#   loglogistic_shape = 1.20,
#   loglogistic_scale = 23.0,
#   alpha_AFT = 0.012
# )
# 
# jm_res

###### Optional: calibrate all five Table 2 sAFT-JM scenarios for the chosen baseline ######
jm_res_all <- calibrate_lambda_c_jm_all_scenarios(
  target_cens = 0.50,
  seeds = 1:500,
  n_patients = 1100,
  n_cores = 8,
  lambda_bounds = c(0, 1),
  tol = 0.001,
  max_FU = 120,
  aft_mode = "Weibull", # "Weibull" or "loglogistic"
  weibull_shape = 1.30, # 0.7, 0.9, 1.3, 2.1
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  alpha_AFT = 0.012
)

jm_res_all

# max_FU = 120, target_cens = 0.50, loglogistic(1.20, 23.0), lambda_c values of Table 2 scenario 1-5: 0.009948730 0.003906250 0.004272461 0.017028809 0.017211914, and avg_censoring 0.5003491 0.5002273 0.5009836 0.4994745 0.4994018
# max_FU = 120, target_cens = 0.50, Weibull(0.70, 38.0), lambda_c values of Table 2 scenario 1-5: 0.010864258 0.004638672 0.004882812 0.019042969 0.019287109, and avg_censoring 0.4993636 0.4996927 0.5003436 0.5008727 0.5001236
# max_FU = 120, target_cens = 0.50, Weibull(0.90, 38.0), lambda_c values of Table 2 scenario 1-5: 0.009643555 0.003479004 0.003784180 0.017211914 0.017456055, and avg_censoring 0.5000636 0.5007182 0.5002345 0.5005545 0.4991873
# max_FU = 120, target_cens = 0.50, Weibull(1.30, 38.0), lambda_c values of Table 2 scenario 1-5: 0.008666992 0.002227783 0.002593994 0.015563965 0.015716553, and avg_censoring 0.5009473 0.4996927 0.5008764 0.5009491 0.4997655
# max_FU = 120, target_cens = 0.50, Weibull(2.10, 38.0), lambda_c values of Table 2 scenario 1-5: 0.008056641 0.001098633 0.001464844 0.014099121 0.014160156, and avg_censoring 0.4996691 0.4997836 0.4991836 0.5008964 0.4992164






###### Example usage for administrative censoring only in one Table 2 sAFT-JM scenario (run this) ######
# jm_admin_res <- avg_admin_censoring_jm(
#   scenario_id = 1,       # Table 2 scenario: 1, 2, 3, 4, or 5
#   seeds = 1:500,
#   n_patients = 1100,
#   n_cores = 8,
#   max_FU = 72,
#   aft_mode = "Weibull",  # "Weibull" or "loglogistic"
#   # Candidate Weibull baselines:
#   weibull_shape = 0.70,
#   # weibull_shape = 0.90,
#   # weibull_shape = 1.30,
#   # weibull_shape = 2.10,
#   weibull_scale = 38.0,
#   loglogistic_shape = 1.20,
#   loglogistic_scale = 23.0,
#   alpha_AFT = 0.012
# )
# 
# jm_admin_res

###### Optional: administrative censoring only for all five Table 2 sAFT-JM scenarios ######
jm_admin_res_all <- avg_admin_censoring_jm_all_scenarios(
  seeds = 1:500,
  n_patients = 1100,
  n_cores = 8,
  max_FU = 120,
  aft_mode = "Weibull", # "Weibull" or "loglogistic"
  weibull_shape = 0.9, # 0.7, 0.9, 1.3, 2.1
  weibull_scale = 38.0,
  loglogistic_shape = 1.20,
  loglogistic_scale = 23.0,
  alpha_AFT = 0.012
)

jm_admin_res_all

# max_FU = 72, loglogistic(1.20, 23.0), avg_censoring of Table 2 scenario 1-5: 0.4163582 0.5477473 0.5430855 0.3079255 0.3040436
# max_FU = 72, Weibull(0.70, 38.0), avg_censoring of Table 2 scenario 1-5: 0.4234927 0.5296364 0.5268291 0.3143327 0.3108055
# max_FU = 72, Weibull(0.90, 38.0), avg_censoring of Table 2 scenario 1-5: 0.4381073 0.5665818 0.5625055 0.3019745 0.2981564
# max_FU = 72, Weibull(1.30, 38.0), avg_censoring of Table 2 scenario 1-5: 0.4686127 0.6295800 0.6250236 0.2846818 0.2798309
# max_FU = 72, Weibull(2.10, 38.0), avg_censoring of Table 2 scenario 1-5: 0.5208509 0.7135055 0.7097036 0.2775891 0.2743527

# max_FU = 96, loglogistic(1.20, 23.0), avg_censoring of Table 2 scenario 1-5: 0.3358327 0.4696655 0.4618291 0.2443091 0.2393473
# max_FU = 96, Weibull(0.70, 38.0), avg_censoring of Table 2 scenario 1-5: 0.3486964 0.4608709 0.4559527 0.2460382 0.2422836
# max_FU = 96, Weibull(0.90, 38.0), avg_censoring of Table 2 scenario 1-5: 0.3432418 0.4843345 0.4780073 0.2219727 0.2175782
# max_FU = 96, Weibull(1.30, 38.0), avg_censoring of Table 2 scenario 1-5: 0.3333727 0.5224927 0.5150164 0.1867309 0.1834345
# max_FU = 96, Weibull(2.10, 38.0), avg_censoring of Table 2 scenario 1-5: 0.3157691 0.5753436 0.5674418 0.1607618 0.1599491

# max_FU = 120, loglogistic(1.20, 23.0), avg_censoring of Table 2 scenario 1-5: 0.2782564 0.4080382 0.3983382 0.1999982 0.1952182
# max_FU = 120, Weibull(0.70, 38.0), avg_censoring of Table 2 scenario 1-5: 0.2909709 0.4065145 0.3995782 0.1975455 0.1937400
# max_FU = 120, Weibull(0.90, 38.0), avg_censoring of Table 2 scenario 1-5: 0.2701600 0.4167018 0.4077364 0.1660855 0.1617945
# max_FU = 120, Weibull(1.30, 38.0), avg_censoring of Table 2 scenario 1-5: 0.2321873 0.4332055 0.4228909 0.1240382 0.1222836
# max_FU = 120, Weibull(2.10, 38.0), avg_censoring of Table 2 scenario 1-5: 0.17286000 0.46067818 0.44663091 0.08716182 0.08703636
