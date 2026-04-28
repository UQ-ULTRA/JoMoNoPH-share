# =============================================================================
# BAYESIAN JOINT MODELLING - MODEL COMPARISON (BP/GB/GP) - SLURM VERSION
# ENZAMET SCENARIO
# =============================================================================
# 
# OVERVIEW:
# This script performs simulation and Bayesian inference for joint models,
# comparing different baseline hazard specifications.
# Designed for parallel execution on HPC clusters using SLURM array jobs.
# ADAPTED FOR ENZAMET DATA GENERATION SCHEME
#
# HOW TO READ THIS SCRIPT (SIGNPOSTS):
#   - Think of this file as an end-to-end pipeline: simulate → fit LME/AFT → fit Stan JM → save batch.
#   - The ENZAMET data generator is sourced dynamically inside run_simulation() via:
#       ~/JoMoNoPH-share/Bunya/ENZAMET/Data-Gen-ENZAMET-XX.R  (scenario = 1..7)
#     so the simulation settings (D, betas, Weibull/loglogistic params, visit schedule, etc.)
#     come from the chosen generator script for that scenario.
#   - We retain flexibility for multiple baseline hazard models (BP/GB/GB_Quantile/GP).
#     Only one is fit per run, controlled by joint_model_type.
#   - This is an HPC/SLURM “array job” script: each SLURM task runs one *batch*
#     (e.g., 100 trials) for one configuration row in config_grid, then writes one .rds.
#
# WORKFLOW (HIGH LEVEL):
# 1. Pre-processing and setup (libraries, CmdStan path, helper sources)
# 2. convert_data_from_models(): takes lme() + survreg() fits and builds Stan data list
# 3. run_simulation(): (a) source ENZAMET generator for scenario
#                    (b) simulate joint data
#                    (c) fit LME on longitudinal outcome (unless AFT-only model)
#                    (d) fit simple AFT model (survreg) as a staging model / prior anchor
#                    (e) fit chosen Stan model (BP/GB/GB_Quantile/GP/BP_aft_alt/BP_aft_logy/GB_aft/BP_aft)
#                    (f) return summaries + diagnostics + runtime + LOO/WAIC
# 4. SLURM job mapping: SLURM_ARRAY_TASK_ID → (config_id, batch_id) → trial_ids
# 5. Results saved as one RDS per batch; script skips work if output already exists
#
# MODEL OPTIONS (baseline hazard in the Stan joint model):
# - BP: Bernstein Polynomial baseline hazard (joint model)
# - GB: Gaussian Basis baseline hazard (joint model)
# - GB_Quantile: Gaussian Basis with quantile-based knots (joint model)
# - GP: Gaussian Process baseline hazard (joint model)
# - BP_aft_alt: AFT model only, Bernstein Polynomials with ordered parameters, no roughness penalty
# - BP_aft_logy: AFT model only, Bernstein Polynomials on log(y), no roughness penalty
# - GB_aft: AFT model only, Gaussian Basis, no roughness penalty
# - GB_aft_logy: AFT model only, Gaussian Basis on log(y), no roughness penalty
# - GB_aft_logy_adaptive: AFT model only, Gaussian Basis on log(y) with adaptive grid, no roughness penalty
# - BP_aft: AFT model only, Bernstein Polynomials, no roughness penalty
#
# IMPORTANT DESIGN NOTES / ASSUMPTIONS:
# - The survival “staging” model fit by survreg() uses dist = "exponential".
#   This is used to construct consistent design matrices and initialize priors/starting values;
#   it is not intended to be the final survival model of interest.
# - AFT/PH language: the ENZAMET generator supports multiple event-time formulations.
#   This HPC pipeline is primarily used with AFT-style data generation (see config_grid aft_mode).
# - Reproducibility: seeds are made globally unique via unique_seed = 1e6 * config_hash + i.
#   This prevents collisions across SLURM tasks/configurations.
# - Resource model: cpus_per_sim = chains * threads_per_chain; workers are set from
#   SLURM_CPUS_PER_TASK / cpus_per_sim, so each worker runs one Stan fit at a time.
# - Timeouts: each trial is capped by an elapsed wall-time limit (safe_elapsed).
#
# WHERE TO CHANGE THINGS (COMMON EDIT POINTS):
# - Choose baseline hazard model: joint_model_type <- "BP" / "GB" / "GB_Quantile" / "GP" / "BP_aft_alt" / "BP_aft_logy" / "GB_aft" / "GB_aft_logy" / "GB_aft_logy_adaptive" / "BP_aft"
# - Choose basis dimension: k_bases (often overridden by config_grid)
# - Choose scenarios and censoring: config_base, lambda_map, cfg_default/cfg_zero
# - Batch sizing: batch_size, sims_per_config
# - Stan sampling controls: iter_warmup, iter_sampling, max_treedepth, chains/threads
#
# OUTPUTS (PER BATCH RDS):
# - Posterior summaries for the joint model + diagnostic_summary()
# - LME fixed effects for comparison
# - Runtime breakdown (LME, Stan compile, Stan sampling, total)
# - LOO + WAIC estimates (if computation succeeds)
#
# =============================================================================


# Pre-processing
##########

rm(list = ls())

library(tidyverse)
library(survival)
library(MASS)
library(nlme)
library(future)
library(furrr)
library(simsurv)
library(cmdstanr)
library(R.utils)

# Set CmdStan path globally for main session
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")

# Load helper functions
source("~/JoMoNoPH-share/Bunya/Auxiliary/BP-functions.R")
# ENZAMET data generator will be sourced dynamically based on scenario parameter

##########


# Model parameters
##########

chains <- 4
threads_per_chain <- 1
parallel_chains <- chains
cpus_per_sim <- chains * threads_per_chain

# Number of basis functions for baseline hazard
# Will be set from config_grid during execution
k_bases = 5  # Default, will be overridden

# Margin on the log-time scale used to extend fixed support for log(y) models.
# Current value preserves the existing support rule:
#   exp(2) ≈ 7.4x wider on each side of the observed time range.
log_support_margin <- 2

# Store a compact posterior reconstruction of the AFT-only baseline
# hazard/cumulative hazard for plotting. This avoids saving full CmdStan fit
# objects; joint-model baseline reconstruction needs separate model-specific
# logic and is intentionally not handled here yet.
baseline_recon_grid_size <- 200
baseline_recon_draws_max <- 500

# Joint model selection: "BP" (Bernstein Polynomial), "GB" (Gaussian Basis), "GB_Quantile", "GP" (Gaussian Process), "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive", or "BP_aft"
joint_model_type <- "BP_aft"  # Change to "BP", "GB", "GB_Quantile", "GP", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive", or "BP_aft" as needed

if (joint_model_type %in% c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")) {
  library(aftgee) # alternative semiparametric AFT model
  library(smoothSurv) # alternative semiparametric AFT model
  library(rstpm2) # alternative semiparametric AFT model
  # library(spsurv) # alternative semiparametric AFT model
}
##########


# AFT-only baseline reconstruction helpers
##########

extract_indexed_draws <- function(draws_mat, parameter) {
  idx <- grep(paste0("^", parameter, "\\["), colnames(draws_mat))
  if (length(idx) == 0) {
    stop("No posterior draws found for parameter: ", parameter)
  }

  parameter_index <- as.integer(sub(
    paste0("^", parameter, "\\[([0-9]+)\\]$"),
    "\\1",
    colnames(draws_mat)[idx]
  ))
  idx <- idx[order(parameter_index)]
  draws_mat[, idx, drop = FALSE]
}

select_reconstruction_draws <- function(n_draws, max_draws) {
  n_keep <- min(n_draws, max_draws)
  unique(round(seq(1, n_draws, length.out = n_keep)))
}

summarise_curve_draws <- function(draw_matrix, prefix) {
  tibble::tibble(
    !!paste0(prefix, "_mean") := colMeans(draw_matrix),
    !!paste0(prefix, "_q2.5") := apply(draw_matrix, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    !!paste0(prefix, "_q97.5") := apply(draw_matrix, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
}

build_fixed_basis_knots <- function(m) {
  if (m == 1) {
    return(list(mu = 0.5, sigma = 0.5))
  }

  list(
    mu = seq(0, 1, length.out = m),
    sigma = rep(2.0 / (3.0 * (m - 1.0)), m)
  )
}

build_fixed_logy_gaussian_knots <- function(m, log_y_lower, log_y_upper) {
  span <- log_y_upper - log_y_lower
  if (span <= 0) {
    stop("log_y_upper must be greater than log_y_lower.")
  }

  if (m == 1) {
    return(list(
      mu = 0.5 * (log_y_lower + log_y_upper),
      sigma = 0.5 * span
    ))
  }

  list(
    mu = seq(log_y_lower, log_y_upper, length.out = m),
    sigma = rep(2.0 * span / (3.0 * (m - 1.0)), m)
  )
}

bernstein_pdf_matrix <- function(u, m) {
  basis <- sapply(seq_len(m), function(k) stats::dbeta(u, shape1 = k, shape2 = m - k + 1))
  if (m == 1) matrix(basis, ncol = 1) else basis
}

bernstein_cdf_matrix <- function(u, m) {
  basis <- sapply(seq_len(m), function(k) stats::pbeta(u, shape1 = k, shape2 = m - k + 1))
  if (m == 1) matrix(basis, ncol = 1) else basis
}

gaussian_pdf_matrix <- function(x, mu, sigma) {
  basis <- sapply(seq_along(mu), function(k) stats::dnorm(x, mean = mu[k], sd = sigma[k]))
  if (length(mu) == 1) matrix(basis, ncol = 1) else basis
}

gaussian_cdf_matrix <- function(x, mu, sigma) {
  basis <- sapply(seq_along(mu), function(k) stats::pnorm(x, mean = mu[k], sd = sigma[k]))
  if (length(mu) == 1) matrix(basis, ncol = 1) else basis
}

build_aft_baseline_reconstruction <- function(fit,
                                              stan_data,
                                              joint_model_type,
                                              t_max,
                                              n_grid = 200,
                                              max_draws = 500) {
  supported_models <- c("BP_aft", "GB_aft", "BP_aft_logy", "GB_aft_logy")
  # BP_aft_alt, GB_aft_logy_adaptive, and joint models need separate
  # reconstruction logic; do not silently use the wrong basis/support rule.
  if (!joint_model_type %in% supported_models) {
    return(NULL)
  }

  draws_mat <- as.matrix(posterior::as_draws_matrix(
    fit$draws(variables = c("beta_surv", "gamma"))
  ))
  keep_rows <- select_reconstruction_draws(nrow(draws_mat), max_draws)
  draws_mat <- draws_mat[keep_rows, , drop = FALSE]

  beta_draws <- extract_indexed_draws(draws_mat, "beta_surv")
  gamma_draws <- extract_indexed_draws(draws_mat, "gamma")
  m <- ncol(gamma_draws)
  n_draws <- nrow(gamma_draws)
  eps <- 1e-6

  is_logy_model <- joint_model_type %in% c("BP_aft_logy", "GB_aft_logy")
  t_grid <- seq(if (is_logy_model) eps else 0, t_max, length.out = n_grid)

  h0_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_grid)
  H0_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_grid)
  tau_aft <- rep(NA_real_, n_draws)

  if (joint_model_type %in% c("BP_aft", "GB_aft")) {
    time <- stan_data$time
    X_surv <- stan_data$X_surv

    for (draw_id in seq_len(n_draws)) {
      linpred <- as.vector(X_surv %*% beta_draws[draw_id, ])
      y <- exp(-linpred) * time
      tau_aft[draw_id] <- max(y) + eps
      u <- pmin(pmax(t_grid / tau_aft[draw_id], eps), 1 - eps)

      if (joint_model_type == "BP_aft") {
        pdf_basis <- bernstein_pdf_matrix(u, m)
        cdf_basis <- bernstein_cdf_matrix(u, m)
      } else {
        knots <- build_fixed_basis_knots(m)
        pdf_basis <- gaussian_pdf_matrix(u, knots$mu, knots$sigma)
        cdf_basis <- gaussian_cdf_matrix(u, knots$mu, knots$sigma)
      }

      h0_draws[draw_id, ] <- as.vector(pdf_basis %*% gamma_draws[draw_id, ]) / tau_aft[draw_id]
      H0_draws[draw_id, ] <- as.vector(cdf_basis %*% gamma_draws[draw_id, ])
      H0_draws[draw_id, t_grid == 0] <- 0
    }
  } else if (joint_model_type == "BP_aft_logy") {
    span <- stan_data$log_y_upper - stan_data$log_y_lower
    u <- pmin(
      pmax((log(t_grid) - stan_data$log_y_lower) / span, eps),
      1 - eps
    )
    pdf_basis <- bernstein_pdf_matrix(u, m)
    cdf_basis <- bernstein_cdf_matrix(u, m)
    h_basis <- sweep(pdf_basis, 1, span * t_grid, "/")

    h0_draws <- gamma_draws %*% t(h_basis)
    H0_draws <- gamma_draws %*% t(cdf_basis)
  } else if (joint_model_type == "GB_aft_logy") {
    knots <- build_fixed_logy_gaussian_knots(
      m,
      stan_data$log_y_lower,
      stan_data$log_y_upper
    )
    log_t_grid <- log(t_grid)
    pdf_basis <- gaussian_pdf_matrix(log_t_grid, knots$mu, knots$sigma)
    cdf_basis <- gaussian_cdf_matrix(log_t_grid, knots$mu, knots$sigma)
    h_basis <- sweep(pdf_basis, 1, t_grid, "/")

    h0_draws <- gamma_draws %*% t(h_basis)
    H0_draws <- gamma_draws %*% t(cdf_basis)
  }

  curve <- dplyr::bind_cols(
    tibble::tibble(time = t_grid),
    summarise_curve_draws(h0_draws, "h0"),
    summarise_curve_draws(H0_draws, "H0")
  )

  list(
    model = joint_model_type,
    time_scale = "baseline accelerated time",
    n_draws_used = n_draws,
    tau_aft = if (all(is.na(tau_aft))) NULL else c(
      mean = mean(tau_aft),
      q2.5 = unname(stats::quantile(tau_aft, 0.025)),
      q97.5 = unname(stats::quantile(tau_aft, 0.975))
    ),
    curve = curve
  )
}

##########


# Data conversion functions
##########

# For joint models (BP/GB/GP)
convert_data_from_models <- function(fit_long, fit_surv, surv_data, k_bases = 5) {
  if (!inherits(fit_long, "lme")) stop("'fit_long' must be an object of class 'lme'")
  if (!inherits(fit_surv, "survreg")) stop("'fit_surv' must be an object of class 'survreg'")
  
  # --- Longitudinal design ---
  mf_long_raw <- nlme::getData(fit_long)
  mf_long <- model.frame(formula(fit_long), data = mf_long_raw, na.action = na.omit)
  idx_long <- as.numeric(rownames(mf_long))
  mf_long$id <- mf_long_raw$id[idx_long]
  mf_long$time <- mf_long_raw$time[idx_long]
  mf_long$arm <- mf_long_raw$arm[idx_long]
  
  X_long <- model.matrix(~ time + time_by_arm, data = mf_long)
  
  y_long <- model.response(mf_long)
  id_long <- as.factor(mf_long$id)
  J_1_long <- as.integer(id_long)
  Y_long <- as.numeric(y_long)
  Z_1_1_long <- rep(1, length(Y_long))
  Z_1_2_long <- mf_long$time
  N_1_long <- length(unique(J_1_long))
  
  # --- Survival design ---
  mf_surv <- model.frame(fit_surv)
  idx_surv <- as.numeric(rownames(mf_surv))
  mf_surv$id <- surv_data$id[idx_surv]
  mf_surv$arm <- surv_data$arm[idx_surv]
  y_surv <- model.response(mf_surv)
  mf_surv$time <- y_surv[, "time"]
  mf_surv$status <- y_surv[, "status"]
  mf_surv$arm <- as.numeric(mf_surv$arm)
  X_surv <- model.matrix(~ 0 + arm, data = mf_surv)
  
  # Link to individuals
  J_1_unique <- match(unique(mf_surv$id), levels(id_long))
  
  # No centering 
  mf_surv$time_by_arm <- mf_surv$time * mf_surv$arm
  X_long_surv <- model.matrix(~ time + time_by_arm, data = mf_surv)
  
  # --- Empirical prior components ---
  mf_long$X_long <- X_long  
  beta_long_mean <- nlme::fixef(fit_long)
  beta_long_sd <- rep(10, length(beta_long_mean))
  beta_surv_mean <- coef(fit_surv)[-1]  # omit intercept
  beta_surv_sd <- rep(10, length(beta_surv_mean))
  
  vc_long <- VarCorr(fit_long)
  sd_b <- sqrt(as.numeric(vc_long[c("(Intercept)", "time"), "Variance"]))
  
  s_long <- summary(fit_long)$sigma
  s_surv <- summary(fit_surv)$scale
  
  # Calibrate alpha_tilde prior scale using clinically interpretable parameterization
  # Delta_y: MCID expressed in SD units of Y_long (e.g., 1.96 ~ 2 SD)
  # max_f: Maximum plausible multiplicative change in survival time (e.g., 2x)
  Delta_y <- 1.96
  max_f <- 2
  sigma_alpha <- log(max_f) / Delta_y
  alpha_tilde_sd <- sigma_alpha * (s_long / s_surv)
  
  list(
    # Longitudinal
    N_long = nrow(X_long),
    K_long = ncol(X_long),
    Y_long = Y_long,
    X_long = X_long,
    J_1_long = J_1_long,
    N_1_long = N_1_long,
    Z_1_1_long = Z_1_1_long,
    Z_1_2_long = Z_1_2_long,
    
    # Survival
    n = nrow(X_surv),
    q = ncol(X_surv),
    status = as.numeric(mf_surv$status),
    time = mf_surv$time,
    X_surv = X_surv,
    
    # Bernstein basis / Gaussian basis / GP
    m = k_bases,
    
    # Linking
    X_long_surv = X_long_surv,
    J_1_unique = J_1_unique,
    
    # Priors
    beta_long_prior_mean = beta_long_mean,
    beta_long_prior_scale = beta_long_sd,
    beta_surv_prior_mean = beta_surv_mean,
    beta_surv_prior_scale = beta_surv_sd,
    sd_b_prior_mean = sd_b,
    sd_b_prior_scale = rep(5, length(sd_b)),
    s_long = s_long,
    s_surv = s_surv,
    alpha_tilde_sd = alpha_tilde_sd,
    
    # Priors for gamma scaling only - all others left unstandardized
    location_beta_real = 0,
    scale_beta_real = 1,
    std_real = 1,
    means_real = 0,
    
    # Smoothing parameters
    r = 2,        # second-order difference penalty
    tau_h = 0.01  # smoothing parameter
  )
}

convert_data_aft_only <- function(fit_surv,
                                  surv_data,
                                  k_bases = 5) {
  
  if (!inherits(fit_surv, "survreg")) {
    stop("'fit_surv' must be an object of class 'survreg'")
  }
  
  # -------------------------------------------------
  # Survival design
  # -------------------------------------------------
  mf_surv <- model.frame(fit_surv)
  idx_surv <- as.numeric(rownames(mf_surv))
  
  mf_surv$id     <- surv_data$id[idx_surv]
  mf_surv$arm    <- as.numeric(surv_data$arm[idx_surv])
  
  y_surv <- model.response(mf_surv)
  mf_surv$time   <- y_surv[, "time"]
  mf_surv$status <- y_surv[, "status"]
  
  X_surv <- model.matrix(~ 0 + arm, data = mf_surv)
  
  # -------------------------------------------------
  # Empirical priors (survival only)
  # -------------------------------------------------
  beta_surv_mean <- coef(fit_surv)[-1]
  beta_surv_sd   <- rep(10, length(beta_surv_mean))
  
  s_surv <- summary(fit_surv)$scale
  
  # -------------------------------------------------
  # Stan data list
  # -------------------------------------------------
  list(
    n       = nrow(X_surv),
    q       = ncol(X_surv),
    time    = mf_surv$time,
    status  = as.numeric(mf_surv$status),
    X_surv  = X_surv,
    
    # Bernstein basis
    m = k_bases,
    
    # Priors
    beta_surv_prior_mean  = beta_surv_mean,
    beta_surv_prior_scale = beta_surv_sd,
    s_surv = s_surv,
    
    # Smoothing
    r     = 2,
    tau_h = 0.01
  )
}

##########


# Simulation wrapper function
##########

run_simulation <- function(i, n_patients, lambda_c, aft_mode, max_FU, config_hash, joint_model_type = "BP", k_bases = 5, scenario = 2) {
  set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")
  
  # Load the appropriate ENZAMET data generator for this scenario
  scenario_file <- sprintf("~/JoMoNoPH-share/Bunya/ENZAMET/Data-Gen-ENZAMET-%02d.R", scenario)
  source(scenario_file)
  
  # New (globally unique):
  unique_seed <- 1e6 * config_hash + i    
  
  # Generate data using ENZAMET data generator
  sim_dat <- simulate_joint_dataset(
    # D = D,
    # beta_0 = beta_0,
    # beta_1 = beta_1,
    # beta_2 = beta_2,
    # sigma_e = sigma_e,
    D = matrix(c(0^2, 0, 0, 0^2), 2, 2), # CAUTION: overwrite values for aft only scenarios
    beta_0 = 0,
    beta_1 = 0,
    beta_2 = 0,
    sigma_e = 0,
    log_HR = log_HR, # beta_surv = log_HR for PH models, should always be set in the data generator, double check the data generator files
    log_AF = log_AF, # beta_surv = log_AF for AFT models, should always be set in the data generator, double check the data generator files
    # alpha_PH = alpha_PH,
    # alpha_AFT = alpha_AFT,
    alpha_PH = 0, # CAUTION: overwrite values for aft only scenarios
    alpha_AFT = 0,
    weibull_shape = weibull_shape,
    weibull_scale = weibull_scale,
    loglogistic_shape = loglogistic_shape, 
    loglogistic_scale = loglogistic_scale,
    visit = visit,
    seed = unique_seed,
    n_patients = n_patients,
    max_FU = max_FU,
    lambda_c = lambda_c,
    aft_mode = aft_mode,  # "PH" or "Weibull" or "loglogistic"
    link_type = "value"   # "value" or "slope"
  )
  
  # ------------------------------------------------------------------------------ 
  # Fit longitudinal model on raw variables 
  # ------------------------------------------------------------------------------ 
  # Augment longitudinal data with time-by-arm interaction
  sim_dat$longitudinal <- sim_dat$longitudinal %>%
    mutate(time_by_arm = time * arm)
  
  # Track runtime for LME fit
  t_lme_start <- Sys.time()
  fit_long <- NULL
  if (!joint_model_type %in% c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")) {
    fit_long <- nlme::lme(
      fixed = Y_obs ~ time + time_by_arm,
      random = ~ time | id,
      data = sim_dat$longitudinal,
      na.action = na.omit,
      control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
    )
  }
  t_lme_elapsed <- as.numeric(difftime(Sys.time(), t_lme_start, units = "secs"))
  
  # ------------------------------------------------------------------------------ 
  # Prepare survival data with consistent variable names (rename only) 
  # ------------------------------------------------------------------------------ 
  sim_dat$survival <- sim_dat$survival %>%
    dplyr::rename(time = T_obs)
  
  
  #### fit alternative semiparametric AFT models (aftgee, aftssr, smoothSurv, rstpm2, spsurv) ####
  if (joint_model_type %in% c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")) {
    # Weibull parametric AFT model
    fit_weibull <- tryCatch({
      survival::survreg(Surv(time, status) ~ arm, data = sim_dat$survival, dist = "weibull")
    }, error = function(e) {
      message("Weibull AFT fit failed: ", e$message)
      NULL
    })
    beta_surv_weibull = if (is.null(fit_weibull)) NA_real_ else fit_weibull$coefficients[-1]
    fit_loglogistic <- tryCatch({
      survival::survreg(Surv(time, status) ~ arm, data = sim_dat$survival, dist = "loglogistic")
    }, error = function(e) {
      message("Log-logistic AFT fit failed: ", e$message)
      NULL
    })
    beta_surv_loglogistic = if (is.null(fit_loglogistic)) NA_real_ else fit_loglogistic$coefficients[-1]
    # least-square-based method (with bootstrap method for variance estimation)
    fit_aftgee <- tryCatch({
      aftgee::aftgee(Surv(time, status) ~ arm, data = sim_dat$survival)
    }, error = function(e) {
      message("aftgee fit failed: ", e$message)
      NULL
    })
    beta_surv_aftgee = if (is.null(fit_aftgee)) NA_real_ else fit_aftgee$coef.res[-1]
    # rank-based method (with induced smoothing method for variance estimation)
    fit_aftsrr <- tryCatch({
      aftgee::aftsrr(Surv(time, status) ~ arm, data = sim_dat$survival, se = "ISMB")
    }, error = function(e) {
      message("aftgee fit failed: ", e$message)
      NULL
    })
    beta_surv_aftsrr = if (is.null(fit_aftsrr)) NA_real_ else fit_aftsrr$beta
    # Komarek's likelihood-based error distribution smoothing method
    fit_smoothSurv <- tryCatch({
      smoothSurv::smoothSurvReg(Surv(time, status) ~ arm, data = sim_dat$survival)
    }, error = function(e) {
      message("smoothSurv fit failed: ", e$message)
      NULL
    })
    beta_surv_smoothSurv = if (is.null(fit_smoothSurv)) {
      NA_real_
    } else {
      fit_smoothSurv$regres[-c(1, NROW(fit_smoothSurv$regres)-1, NROW(fit_smoothSurv$regres)), 1]
    }
    # Crowther et al. (2022) flexible parametric AFT method
    fit_rstpm2 <- tryCatch({
      rstpm2::aft(Surv(time, status) ~ arm, data = sim_dat$survival, df = k_bases)
    }, error = function(e) {
      message("rstpm2 fit failed: ", e$message)
      NULL
    })
    # beta_surv_rstpm2 = summary(fit_rstpm2)@coef[-((NROW(summary(fit_rstpm2)@coef)-k_bases+1):NROW(summary(fit_rstpm2)@coef)),1]
    beta_surv_rstpm2 = if (is.null(fit_rstpm2)) {
      NA_real_
    } else {
      fit_rstpm2@coef[-c((length(fit_rstpm2@coef)-k_bases+1):length(fit_rstpm2@coef))]
    }
    
    # # spsurv semiparametric AFT with Bernstein polynomials (mle)
    # fit_spsurv_mle <- tryCatch({
    #   spsurv::spbp(Surv(time, status) ~ arm, data = sim_dat$survival, approach = "mle", model = "aft", degree = k_bases)
    # }, error = function(e) {
    #   message("spsurv (mle) fit failed: ", e$message)
    #   NULL
    # })
    # beta_surv_spsurv_mle = fit_spsurv_mle$coefficients[-((length(fit_spsurv_mle$coefficients)-k_bases+1):length(fit_spsurv_mle$coefficients))]
    # # spsurv semiparametric AFT with Bernstein polynomials (Bayesian using Stan)
    # fit_spsurv_bayes <- tryCatch({
    #   spsurv::spbp(Surv(time, status) ~ arm, data = sim_dat$survival, approach = "bayes", model = "aft", cores = cpus_per_sim, degree = k_bases, iter = 2000, chains = 4)
    # }, error = function(e) {
    #   message("spsurv (bayes) fit failed: ", e$message)
    #   NULL
    # })
    # spsurv_stan_summary = rstan::summary(fit_spsurv_bayes$stanfit)$summary
    # beta_surv_spsurv_bayes = spsurv_stan_summary[grepl("^beta\\[[0-9]+\\]$", rownames(spsurv_stan_summary)), 1]
    # rm(spsurv_stan_summary)
  }
  
  # ------------------------------------------------------------------------------ 
  # Fit survival model on raw variables 
  # ------------------------------------------------------------------------------ 
  fit_surv <- survival::survreg(
    Surv(time, status) ~ arm,
    data = sim_dat$survival,
    dist = "exponential",
    x = TRUE
  )
  
  # Prepare data for joint model (BP, GB, GB_Quantile, GP, BP_aft_alt, BP_aft_logy, GB_aft, GB_aft_logy, GB_aft_logy_adaptive, or BP_aft)
  if (joint_model_type == "BP") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist.stan"
    # stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist-DM.stan"  # 20260129
  } else if (joint_model_type == "GB") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis-JM-Hist.stan"
  } else if (joint_model_type == "GB_Quantile") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis-JM-Hist-Quantile.stan"
  } else if (joint_model_type == "BP_aft_alt") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Bernstein-Polynomials-alt.stan"
  } else if (joint_model_type == "BP_aft_logy") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Bernstein-Polynomials-logy.stan"
  } else if (joint_model_type == "GB_aft") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis.stan"
  } else if (joint_model_type == "GB_aft_logy") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis-logy.stan"
  } else if (joint_model_type == "GB_aft_logy_adaptive") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis-logy-adaptive-grid.stan"
  } else if (joint_model_type == "GP") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Process-JM-Hist.stan"
  } else if (joint_model_type == "BP_aft") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Bernstein-Polynomials.stan"
  } else {
    stop(sprintf("Unknown joint_model_type: %s. Use 'BP', 'GB', 'GB_Quantile', 'GP', 'BP_aft_alt', 'BP_aft_logy', 'GB_aft', 'GB_aft_logy', 'GB_aft_logy_adaptive', or 'BP_aft'.", joint_model_type))
  }
  
  if (joint_model_type %in% c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")) {
    stan_data_joint <- convert_data_aft_only(fit_surv = fit_surv,
                                             surv_data = sim_dat$survival,
                                             k_bases = k_bases)
  } else {
    stan_data_joint <- convert_data_from_models(fit_long = fit_long,
                                                fit_surv = fit_surv,
                                                surv_data = sim_dat$survival,
                                                k_bases = k_bases)
  }

  if (joint_model_type %in% c("BP_aft_logy", "GB_aft_logy")) {
    log_t <- log(stan_data_joint$time)
    L <- min(log_t)
    U <- max(log_t)

    stan_data_joint$log_y_lower <- L - log_support_margin
    stan_data_joint$log_y_upper <- U + log_support_margin
  }

  if (joint_model_type == "GB_aft_logy_adaptive") {
    log_t <- log(stan_data_joint$time)
    L <- min(log_t)
    U <- max(log_t)

    stan_data_joint$grid_lower <- L - log_support_margin
    stan_data_joint$grid_upper <- U + log_support_margin

    stan_data_joint$a_prior_mean <- mean(log_t) # centre at the mean instead of 0 for better scaling
    stan_data_joint$a_prior_scale <- 2
    stan_data_joint$log_b_prior_mean <- 0
    stan_data_joint$log_b_prior_scale <- 0.5
  }
  
  # Add quantile-based knots for GB_Quantile model
  if (joint_model_type == "GB_Quantile") {
    source("~/JoMoNoPH-share/Bunya/calculate_quantile_knots.R")
    knots <- calculate_quantile_knots(time = sim_dat$survival$time, m = k_bases, coverage = 0.6)
    stan_data_joint$mu_knots <- knots$mu_knots
    stan_data_joint$sigma_knots <- knots$sigma_knots
  }
  
  # Add GP hyperparameters for GP model
  if (joint_model_type == "GP") {
    # Inverse-gamma priors: IG(alpha, beta) with mean = beta/(alpha-1) for alpha > 1
    # For marginal SD: target mean ~0.5, use IG(3, 1) -> mean = 0.5
    # For length scale: target mean ~0.2, use IG(3, 0.4) -> mean = 0.2
    stan_data_joint$gp_scale_prior_alpha <- 3
    stan_data_joint$gp_scale_prior_beta <- 1
    stan_data_joint$gp_length_prior_alpha <- 3
    stan_data_joint$gp_length_prior_beta <- 0.4
  }
  
  
  # Initialisation function for joint models
  init_fun_joint <- function(sim_dat, stan_data, chains = 4, joint_model_type = "BP") {
    # --- Fallback values ---
    # These are used if the lme() or survreg() model fails to converge or throws an error.
    # They are deliberately safe, neutral guesses meant to allow Stan to initialize without biasing the sampler.
    fallback <- list(
      beta_long = rep(0, stan_data$K_long),     # Intercept, time, time_by_arm
      sigma_long = 10,                          # Reasonable residual SD
      sd_1_long = c(1, 1),                      # Conservative random effects SDs
      beta_surv = rep(0, stan_data$q),          # survival covariates
      gamma_scaled = rep(0.1, stan_data$m),     # Small positive Bernstein weights
      alpha_tilde = 0                           # No association
    )
    
    # --- Estimate longitudinal model using lme() ---
    # These provide empirical starting values for beta_long, sigma_long, and sd_1_long
    LMM <- nlme::lme(
      fixed = Y_obs ~ time + time_by_arm,
      random = ~ time | id,
      data = sim_dat$longitudinal,
      control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
    )
    
    # Fixed effects (intercept and covariate effects)
    fallback$beta_long <- as.numeric(nlme::fixef(LMM))
    
    # Residual error SD
    fallback$sigma_long <- summary(LMM)$sigma
    
    # Random intercept/slope SDs
    vc <- nlme::VarCorr(LMM)
    fallback$sd_1_long <- sqrt(as.numeric(vc[1:2, "Variance"]))
    
    # --- Estimate survival treatment effect using survreg() ---
    # This is used to initialise MCMC sampling
    AFT <- survival::survreg(Surv(time, status) ~ arm, data = sim_dat$survival)
    fallback$beta_surv <- rep(coef(AFT)[2], stan_data$q)
    
    # --- Construct per-chain initial values ---
    # Adds small jitter to allow overdispersed starting values across chains
    init_list <- lapply(1:chains, function(chain_id) {
      # Jittered SDs
      sd_jittered <- fallback$sd_1_long + abs(rnorm(2, 0, 0.05))
      sd_jittered <- pmax(sd_jittered, 1e-3)
      
      # Correlation matrix and Cholesky
      rho <- 0.2
      R <- matrix(c(1, rho, rho, 1), 2, 2)
      L_1_long_init <- t(chol(R))
      
      # Base initialization for all models
      init_vals <- list(
        z_1_long = matrix(rnorm(2 * stan_data$N_1_long), nrow = 2),
        L_1_long = L_1_long_init,
        beta_long = fallback$beta_long + rnorm(stan_data$K_long, 0, 0.1),
        sigma_long = abs(rnorm(1, fallback$sigma_long, 1)),
        sd_1_long = sd_jittered
      )
      
      # Add survival components for joint models
      init_vals$beta_surv <- rnorm(stan_data$q, fallback$beta_surv, 0.05)
      init_vals$alpha_tilde <- rnorm(1, 0, 0.1)
      
      # Model-specific parameters
      if (joint_model_type == "GP") {
        # GP-specific initialization
        init_vals$gp_length_scale <- abs(rnorm(1, 0.2, 0.05))
        init_vals$gp_marginal_sd <- abs(rnorm(1, 0.5, 0.1))
        init_vals$f_gp_raw <- rnorm(stan_data$m, 0, 0.5)
      } else {
        # BP/GB/GB_Quantile use gamma weights
        init_vals$gamma <- abs(rnorm(stan_data$m, 0.1, 0.01))
      }
      
      init_vals
    })
    
    return(init_list)
  }
  
  init_fun_aft_only <- function(sim_dat,
                                stan_data,
                                chains = 4,
                                joint_model_type = "BP") {
    
    # ------------------------------------------------------------------
    # 1. Fallback values (aft-only)
    # ------------------------------------------------------------------
    # BP_aft_alt needs m+1 gamma values (for degree m polynomial)
    n_gamma <- ifelse(joint_model_type == "BP_aft_alt", stan_data$m + 1, stan_data$m)
    
    fallback <- list(
      beta_surv = rep(0, stan_data$q),    # survival covariates
      gamma     = rep(0.1, n_gamma)       # Bernstein weights
    )
    
    # ------------------------------------------------------------------
    # 2. Empirical AFT initialisation
    # ------------------------------------------------------------------
    AFT <- tryCatch(
      survival::survreg(Surv(time, status) ~ arm, data = sim_dat$survival),
      error = function(e) NULL
    )
    
    if (!is.null(AFT)) {
      # treatment effect only (consistent with convert_data_aft_only)
      fallback$beta_surv <- rep(coef(AFT)[2], stan_data$q)
    }
    
    # ------------------------------------------------------------------
    # 3. Per-chain initial values
    # ------------------------------------------------------------------
    init_list <- lapply(seq_len(chains), function(chain_id) {
      
      init_vals <- list(
        beta_surv = rnorm(
          stan_data$q,
          mean = fallback$beta_surv,
          sd   = 0.05
        )
      )

      if (joint_model_type == "GB_aft_logy_adaptive") {
        init_vals$a <- rnorm(1, mean = stan_data$a_prior_mean, sd = 0.1)
        init_vals$log_b <- rnorm(1, mean = stan_data$log_b_prior_mean, sd = 0.1)
      }
      
      # Model-specific baseline hazard parameters
      if (joint_model_type == "GP") {
        init_vals$gp_length_scale <- abs(rnorm(1, 0.2, 0.05))
        init_vals$gp_marginal_sd  <- abs(rnorm(1, 0.5, 0.1))
        init_vals$f_gp_raw        <- rnorm(stan_data$m, 0, 0.5)
      } else if (joint_model_type == "BP_aft_alt") {
        # BP_aft_alt requires ordered (monotonically increasing) gamma values
        n_gamma <- length(fallback$gamma)
        gamma_increments <- abs(rnorm(n_gamma, 0.1, 0.02))
        init_vals$gamma <- cumsum(gamma_increments)
      } else {
        # BP / BP_aft_logy / GB / GB_Quantile / GB_aft / GB_aft_logy / GB_aft_logy_adaptive
        n_gamma <- length(fallback$gamma)
        init_vals$gamma <- pmax(
          fallback$gamma + rnorm(n_gamma, 0, 0.01),
          1e-4
        )
      }
      
      init_vals
    })
    
    return(init_list)
  }
  
  # ------
  # Generate initial values for joint model
  if (joint_model_type %in% c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")) {
    init_values_joint <- init_fun_aft_only(sim_dat = sim_dat, stan_data = stan_data_joint, chains = chains, joint_model_type = joint_model_type)
  } else {
    init_values_joint <- init_fun_joint(sim_dat = sim_dat, stan_data = stan_data_joint, chains = chains, joint_model_type = joint_model_type)
  }
  # ========== Fit Joint Model (BP, GB, GB_Quantile, or GP) ==========
  t_compile_joint_start <- Sys.time()
  mod_joint <- cmdstan_model(stan_model_file_joint,
                             cpp_options = list(stan_threads = TRUE),
                             force_recompile = FALSE)
  Sys.chmod(mod_joint$exe_file(), mode = "0755")
  if (!file.exists(mod_joint$exe_file())) {
    stop("Stan joint model compilation failed; binary not found.")
  }
  t_compile_joint_elapsed <- as.numeric(difftime(Sys.time(), t_compile_joint_start, units = "secs"))
  
  t_stan_joint_start <- Sys.time()
  fit_joint <- mod_joint$sample(
    data = stan_data_joint,
    seed = unique_seed,
    init = init_values_joint,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = 1000,
    iter_sampling = 1000,
    max_treedepth = 10)
  t_stan_joint_elapsed <- as.numeric(difftime(Sys.time(), t_stan_joint_start, units = "secs"))
  
  # Compute LOO and WAIC for joint model
  loo_joint <- tryCatch({
    fit_joint$loo(save_psis = FALSE)
  }, error = function(e) {
    message("LOO computation failed for joint model: ", e$message)
    NULL
  })
  
  waic_joint <- tryCatch({
    fit_joint$waic()
  }, error = function(e) {
    message("WAIC computation failed for joint model: ", e$message)
    NULL
  })
  
  # Collect results from joint model (with 2.5th and 97.5th percentiles)
  summ_joint <- fit_joint$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975))
  )
  summ_joint$sim_id <- i
  summ_joint$model <- joint_model_type
  summ_joint$rhat <- fit_joint$summary()$rhat # added to store rhat for all parameters in the summary table
  
  diag_joint <- fit_joint$diagnostic_summary()
  diag_joint$sim_id <- i
  diag_joint$model <- joint_model_type
  
  # LME coefficients for comparison
  coefs <- NULL
  if (!is.null(fit_long)) {
    coefs <- nlme::fixef(fit_long)
  }
  
  # Collect runtimes
  runtimes <- data.frame(
    sim_id = i,
    joint_model = joint_model_type,
    k_bases = k_bases,
    lme_time = t_lme_elapsed,
    joint_compile_time = t_compile_joint_elapsed,
    joint_stan_time = t_stan_joint_elapsed,
    total_time = t_lme_elapsed + t_compile_joint_elapsed + t_stan_joint_elapsed
  )
  
  # Extract LOO/WAIC estimates if available
  loo_joint_est <- if (!is.null(loo_joint)) {
    data.frame(
      elpd_loo = loo_joint$estimates["elpd_loo", "Estimate"],
      p_loo = loo_joint$estimates["p_loo", "Estimate"],
      looic = loo_joint$estimates["looic", "Estimate"]
    )
  } else {
    data.frame(elpd_loo = NA, p_loo = NA, looic = NA)
  }
  
  waic_joint_est <- if (!is.null(waic_joint)) {
    data.frame(
      elpd_waic = waic_joint$estimates["elpd_waic", "Estimate"],
      p_waic = waic_joint$estimates["p_waic", "Estimate"],
      waic = waic_joint$estimates["waic", "Estimate"]
    )
  } else {
    data.frame(elpd_waic = NA, p_waic = NA, waic = NA)
  }

  logy_support <- NULL
  if (all(c("log_y_lower", "log_y_upper") %in% names(stan_data_joint))) {
    logy_support <- list(
      log_y_lower = stan_data_joint$log_y_lower,
      log_y_upper = stan_data_joint$log_y_upper
    )
  }

  baseline_recon <- tryCatch({
    build_aft_baseline_reconstruction(
      fit = fit_joint,
      stan_data = stan_data_joint,
      joint_model_type = joint_model_type,
      t_max = max_FU,
      n_grid = baseline_recon_grid_size,
      max_draws = baseline_recon_draws_max
    )
  }, error = function(e) {
    message("Baseline reconstruction failed: ", e$message)
    NULL
  })
  
  if (joint_model_type %in% c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")) {
    return(list(
      censor_prop = mean(sim_dat$survival$status == 0),
      unique_seed = unique_seed,
      # fit_aftgee = fit_aftgee, # it took too much memory to save the full model objects for aftgee, aftsrr, smoothSurv, rstpm2 and spsurv
      # fit_aftsrr = fit_aftsrr, # so decide to just save the coefficients of interest (treatment effect estimates)
      # fit_smoothSurv = fit_smoothSurv,
      # fit_rstpm2 = fit_rstpm2,
      beta_sAFT = rbind(beta_surv_aftgee, beta_surv_aftsrr, beta_surv_smoothSurv, beta_surv_rstpm2, beta_surv_weibull, beta_surv_loglogistic),
      # beta_sAFT = rbind(beta_surv_aftgee, beta_surv_aftsrr, beta_surv_smoothSurv, beta_surv_rstpm2, beta_surv_spsurv_mle, beta_surv_spsurv_bayes, beta_surv_weibull, beta_surv_loglogistic),
      joint_summary = summ_joint,
      joint_diagnostics = diag_joint,
      lme_coefs = coefs,
      logy_support = logy_support,
      baseline_recon = baseline_recon,
      runtimes = runtimes,
      loo_joint = loo_joint_est,
      waic_joint = waic_joint_est))
  } else {
    return(list(
      censor_prop = mean(sim_dat$survival$status == 0),
      unique_seed = unique_seed,
      joint_summary = summ_joint,
      joint_diagnostics = diag_joint,
      lme_coefs = coefs,
      logy_support = logy_support,
      baseline_recon = baseline_recon,
      runtimes = runtimes,
      loo_joint = loo_joint_est,
      waic_joint = waic_joint_est))
  }
}

##########


# SLURM Job Mapping Logic
##########

# Read job ID from environment
sim_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
message(sprintf("SLURM_ARRAY_TASK_ID (also the sim_id): %d", sim_id))
if (sim_id == 0) stop("SLURM_ARRAY_TASK_ID not found.")

#### Define configuration grid for ENZAMET ####
# ENZAMET uses n_patients=1100, max_FU=72, loglogistic distribution
config_base <- expand.grid(
  n_patients = c(1100),
  aft_mode   = c("loglogistic"), # current options: "Weibull" and/or "loglogistic"
  k_bases    = c(5), # change from 5 to 10 on 20260128
  # scenario   = c(1, 2, 3, 4, 5, 6, 7),  # Data generation scenarios
  # scenario   = c(1, 2, 3, 4, 5),  # Data generation scenarios
  scenario   = c(1, 2, 4),  # Data generation scenarios for aft only
  max_FU     = 72,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

# lambda_map <- c(weibull = 0.023, loglogistic = 0.013) # under scenario 1, give ~50% censoring
lambda_c_table <- data.frame( # need different lambda_c values for different scenarios
  scenario = rep(1:5, each = 2),
  aft_mode = rep(c("Weibull", "loglogistic"), times = 5),
  lambda_c = c(
    # Old Bernstein Polynomial approximation:
    # scenario 1: 0.08154297, 0.05322266
    # scenario 2: 0.05224609, 0.03320312
    # scenario 3: 0.05224609, 0.03320312
    # scenario 4: 0.1289062,  0.08251953
    # scenario 5: 0.1289062,  0.08251953
    #
    # Updated 20260427 after amended Bernstein Polynomial approximation:
    # scenario 1
    0.02832031, 0.02587891, # for aft only and 50% censoring
    # scenario 2
    0.015625, 0.01391602, # for aft only and 50% censoring
    # scenario 3
    0.015625, 0.01391602, # for aft only and 50% censoring
    # scenario 4
    0.04638672, 0.04150391, # for aft only and 50% censoring
    # scenario 5
    0.04638672, 0.04150391 # for aft only and 50% censoring
  ),
  stringsAsFactors = FALSE
)

# cfg_default <- transform(config_base,
#                          lambda_c = lambda_map[tolower(aft_mode)])
cfg_default <- merge(
  config_base,
  lambda_c_table,
  by = c("scenario", "aft_mode"),
  all.x = TRUE,
  sort = FALSE
)
cfg_default$cens_label <- "50"

cfg_zero <- transform(config_base,
                      lambda_c = 0,
                      cens_label = "0") # 0 for no censoring

# Administrative censoring only at max_FU.
# lambda_c < 0 triggers T_obs = min(T_i, max_FU), status = I(T_i <= max_FU).
# Data-Gen-test-DM.R 20260427 administrative-only (aft only) check, loglogistic:
# scenario 1 avg censoring = 0.2023836
# scenario 2/3 avg censoring = 0.3155127
# scenario 4/5 avg censoring = 0.1409127
# Weibull:
# scenario 1 avg censoring = 0.2090255
# scenario 2/3 avg censoring = 0.32196
# scenario 4/5 avg censoring = 0.1311945
cfg_admin <- transform(config_base,
                       lambda_c = -1,
                       cens_label = "admin")

config_grid <- rbind(cfg_default, cfg_zero, cfg_admin)
row.names(config_grid) <- NULL

# Define batch structure (e.g. 100 sims per batch, 1000 total sims)
# batch_size <- 100
# sims_per_config <- 1000
batch_size <- 10 # try 10 first for testing
sims_per_config <- 100 # try 100 first for testing
batches_per_config <- sims_per_config / batch_size  

# Decode SLURM array index into config ID and batch ID
config_id <- ((sim_id - 1) %/% batches_per_config) + 1      
batch_id  <- ((sim_id - 1) %% batches_per_config) + 1       

# Pull scenario settings for this SLURM job
config <- config_grid[config_id, ]
n_patients <- config$n_patients
lambda_c <- config$lambda_c
aft_mode <- config$aft_mode
max_FU <- config$max_FU
k_bases <- config$k_bases
scenario <- config$scenario
cens_label <- config$cens_label

# Define trial numbers for this batch
trial_ids <- ((batch_id - 1) * batch_size + 1):(batch_id * batch_size)

# Hash the config (e.g., using integers from factors)
# config_hash <- as.integer(as.numeric(as.factor(
#   paste(n_patients, lambda_c, aft_mode, k_bases, scenario, max_FU, sep = "_")
# )))
config_hash <- config_id # use config_id as hash since it's already a unique integer identifier for the config
# print config and config_hash for debugging
message(sprintf("Config ID: %d, Config Hash: %d, n_patients:
%d, lambda_c: %.5f, aft_mode: %s, k_bases: %d, scenario: %d, max_FU: %d",
config_id, config_hash, n_patients, lambda_c, aft_mode, k_bases, scenario, max_FU))

##########
# ---- Skip work if this batch output already exists ----
out_file <- sprintf(
  "~/JoMoNoPH-share/Bunya/ENZAMET/output/%s_k%d_n%d_cens%s_%s_scen%d_FU%d_batch%03d.rds",
  joint_model_type, k_bases, n_patients, cens_label, aft_mode, scenario, max_FU, batch_id
)

if (file.exists(out_file)) {
  message(sprintf("SKIP: output already exists for this batch: %s", out_file))
  quit(save = "no", status = 0)
}
# -------------------------------------------------------


# Run all simulations in this batch
##########

message(sprintf(
  "Running batch %d (config_id %d: n=%d, censoring=%s, k=%d, scenario=%d, model=%s, ENZAMET)",
  batch_id, config_id, n_patients, cens_label, k_bases, scenario, joint_model_type
))

# Detect number of available CPUs from SLURM
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
workers <- floor(cpus / cpus_per_sim)
workers <- min(workers, length(trial_ids))

# Set up parallel backend
plan(multisession, workers = workers)

safe_elapsed <- 7200  # 2 hours wall time per trial

results <- future_map(
  trial_ids,
  ~ tryCatch({
    t0 <- Sys.time()
    # Elapsed-only timeout using base R
    ans <- tryCatch({
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
      setTimeLimit(elapsed = safe_elapsed, cpu = Inf, transient = TRUE)
      run_simulation(.x, n_patients, lambda_c, aft_mode, max_FU, config_hash, joint_model_type, k_bases, scenario)
    }, error = function(e) {
      if (grepl("reached elapsed time limit", conditionMessage(e))) NULL else stop(e)
    })
    
    cat(sprintf("Trial %d elapsed: %.1f min\n",
                .x, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    ans
  },
  error = function(e) {
    message(sprintf("Simulation %d failed or timed out: %s", .x, e$message))
    NULL
  }
  ),
  .options = furrr_options(seed = TRUE)
)

##########


# Save the results to file
##########

saveRDS(results, file = out_file)


##########
