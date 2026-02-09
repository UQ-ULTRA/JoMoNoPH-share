# =============================================================================
# BAYESIAN JOINT MODELLING - MODEL COMPARISON (BP/GB/GP)
# ENZAMET SCENARIO - SINGLE SIMULATION
# =============================================================================
# 
# OVERVIEW:
# This script performs simulation and Bayesian inference for joint models,
# comparing different baseline hazard specifications.
# ADAPTED FOR ENZAMET DATA GENERATION SCHEME
#
# WORKFLOW:
# 1. Pre-processing and setup
# 2. Data conversion functions for joint models
# 3. Simulation wrapper function with model selection
# 4. Configuration and execution with runtime tracking
#
# MODEL OPTIONS:
# - BP: Bernstein Polynomial baseline hazard (joint model)
# - GB: Gaussian Basis baseline hazard (joint model)
# - GB_Quantile: Gaussian Basis with quantile-based knots (joint model)
# - GP: Gaussian Process baseline hazard (joint model)
# - aft_only_BP_noRP: AFT model only using Bernstein Polynomials (no joint modelling)
#
# OUTPUTS:
# - Stan model summaries (BP/GB/GP)
# - LME fixed effects comparison
# - Runtime statistics for all fits
#
# =============================================================================

# Pre-processing
##########

rm(list = ls())

library(tidyverse)
library(survival)
library(MASS)
library(nlme)
library(simsurv)
library(rstan)
library(brms)
library(cmdstanr)
library(here)

# Set CmdStan path globally for main session
# set_cmdstan_path(here::here("..", ".cmdstan", "cmdstan-2.36.0"))
set_cmdstan_path("/Users/dingma/.cmdstan/cmdstan-2.36.0")

# Load helper functions
source(here::here("Bunya", "Auxiliary", "BP-functions.R"))
# ENZAMET data generator will be sourced dynamically based on scenario parameter

##########


# Model parameters
##########

chains <- 4
threads_per_chain <- 1
degree <- 5

# Joint model selection: "BP" (Bernstein Polynomial), "GB" (Gaussian Basis), "GB_Quantile", "GP" (Gaussian Process), or "aft_only_BP_noRP (BP without roughness penalty)"
joint_model_type <- "aft_only_BP_noRP"  # Change to "BP", "GB", "GB_Quantile", "GP", or "aft_only_BP_noRP" as needed

##########


# Data conversion functions
##########

# For joint models (BP/GB/GP)
convert_data_from_models <- function(fit_long, fit_surv, surv_data, degree = 5) {
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
    
    # Bernstein basis
    m = degree,
    
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
    
    # Priors for gamma scaling only — all others left unstandardized
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
                                  degree = 5) {

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
    m = degree,

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

run_simulation <- function(i, n_patients, lambda_c, aft_mode, max_FU, unique_seed, joint_model_type = "BP", scenario = 2) {
  # set_cmdstan_path(here::here("..", "cmdstan", "cmdstan-2.36.0"))
  set_cmdstan_path("/Users/dingma/.cmdstan/cmdstan-2.36.0")
  s
  # Load the appropriate ENZAMET data generator for this scenario
  scenario_file <- here::here("Bunya", "ENZAMET", sprintf("Data-Gen-ENZAMET-%02d.R", scenario))
  source(scenario_file)    
    
    # Generate data using ENZAMET data generator
    sim_dat <- simulate_joint_dataset(
      # D = D,
      # beta_0 = beta_0,
      # beta_1 = beta_1,
      # beta_2 = beta_2,
      # sigma_e = sigma_e,
      D = matrix(c(0^2, 0, 0, 0^2), 2, 2),
      beta_0 = 0,
      beta_1 = 0,
      beta_2 = 0,
      sigma_e = 0,
      # log_HR = log_HR,
      # log_AF = log_AF,
      log_HR = 0,
      log_AF = 0,
      # alpha_PH = alpha_PH,
      # alpha_AFT = alpha_AFT,
      alpha_PH = 0,
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
    if (joint_model_type != "aft_only_BP_noRP") {
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
    
    # ------------------------------------------------------------------------------ 
    # Fit survival model on raw variables 
    # ------------------------------------------------------------------------------ 
    fit_surv <- survival::survreg(
      Surv(time, status) ~ arm,
      data = sim_dat$survival,
      dist = "exponential",
      x = TRUE
    )
    
    
    # Prepare data for joint model (BP, GB, GB_Quantile, GP, or aft_only_BP_noRP)
    if (joint_model_type == "BP") {
      stan_model_file_joint <- "Bernstein-Polynomials-JM-Hist.stan"
    } else if (joint_model_type == "GB") {
      stan_model_file_joint <- "Gaussian-Basis-JM-Hist.stan"
    } else if (joint_model_type == "GB_Quantile") {
      stan_model_file_joint <- "Gaussian-Basis-JM-Hist-Quantile.stan"
    } else if (joint_model_type == "GP") {
      stan_model_file_joint <- "Gaussian-Process-JM-Hist.stan"
    } else if (joint_model_type == "aft_only_BP_noRP") {
      stan_model_file_joint <- "Bernstein-Polynomials.stan"
    } else {
      stop(sprintf("Unknown joint_model_type: %s. Use 'BP', 'GB', 'GB_Quantile', 'GP', or 'aft_only_BP_noRP'.", joint_model_type))
    }
    
    if (joint_model_type == "aft_only_BP_noRP") {
      stan_data_joint <- convert_data_aft_only(fit_surv = fit_surv,
                                               surv_data = sim_dat$survival,
                                               degree = 5)
    } else {
      stan_data_joint <- convert_data_from_models(fit_long = fit_long,
                                                  fit_surv = fit_surv,
                                                  surv_data = sim_dat$survival,
                                                  degree = 5)
    }
    
    # Add quantile-based knots for GB_Quantile model
    if (joint_model_type == "GB_Quantile") {
      source(here::here("Bunya", "calculate_quantile_knots.R"))
      knots <- calculate_quantile_knots(time = sim_dat$survival$time, m = 5, coverage = 0.6)
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
      fallback <- list(
        beta_surv = rep(0, stan_data$q),    # survival covariates
        gamma     = rep(0.1, stan_data$m)   # Bernstein weights
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

        # Model-specific baseline hazard parameters
        if (joint_model_type == "GP") {
          init_vals$gp_length_scale <- abs(rnorm(1, 0.2, 0.05))
          init_vals$gp_marginal_sd  <- abs(rnorm(1, 0.5, 0.1))
          init_vals$f_gp_raw        <- rnorm(stan_data$m, 0, 0.5)
        } else {
          # BP / GB / GB_Quantile
          init_vals$gamma <- pmax(
            fallback$gamma + rnorm(stan_data$m, 0, 0.01),
            1e-4
          )
        }

        init_vals
      })

      return(init_list)
    }
    
    # ------
    # Generate initial values for joint model
    if (joint_model_type == "aft_only_BP_noRP") {
      init_values_joint <- init_fun_aft_only(sim_dat = sim_dat, stan_data = stan_data_joint, chains = chains, joint_model_type = joint_model_type)
    } else {
      init_values_joint <- init_fun_joint(sim_dat = sim_dat, stan_data = stan_data_joint, chains = chains, joint_model_type = joint_model_type)
    }
    # ========== Fit Joint Model (BP, GB, GB_Quantile, or GP) ==========
    t_compile_joint_start <- Sys.time()
    mod_joint <- cmdstan_model(here::here("Stan", stan_model_file_joint),
                               cpp_options = list(stan_threads = TRUE),
                               force_recompile = FALSE)
    Sys.chmod(mod_joint$exe_file(), mode = "0755")
    t_compile_joint_elapsed <- as.numeric(difftime(Sys.time(), t_compile_joint_start, units = "secs"))
    
    t_stan_joint_start <- Sys.time()
    fit_joint <- mod_joint$sample(
      data = stan_data_joint,
      seed = unique_seed,
      init = init_values_joint,
      chains = chains,
      parallel_chains = chains,
      threads_per_chain = threads_per_chain,
      iter_warmup = 1000,
      iter_sampling = 1000,
      max_treedepth = 10)
    t_stan_joint_elapsed <- as.numeric(difftime(Sys.time(), t_stan_joint_start, units = "secs"))
    
    # Collect results from joint model
    summ_joint <- fit_joint$summary()
    summ_joint$sim_id <- i
    summ_joint$model <- joint_model_type
    
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
      lme_time = t_lme_elapsed,
      joint_compile_time = t_compile_joint_elapsed,
      joint_stan_time = t_stan_joint_elapsed,
      total_time = t_lme_elapsed + t_compile_joint_elapsed + t_stan_joint_elapsed
    )
    
    return(list(
      joint_summary = summ_joint,
      joint_diagnostics = diag_joint,
      lme_coefs = coefs,
      runtimes = runtimes))
}

##########


#### Configuration and execution ####
##########

# Define simulation settings for ENZAMET
n_patients <- 1100
lambda_c <- 0.01
aft_mode <- "Weibull"  # "Weibull" or "loglogistic"
max_FU <- 72
scenario <- 2  # Data generation scenario (1-7)

# Number of simulations to run
n_sims <- 1

# Set base seed
base_seed <- 51234

# Run simulations sequentially
results <- vector("list", n_sims)

for (i in 1:n_sims) {
  cat(sprintf("\n--- Running ENZAMET simulation %d of %d (Scenario %d, Joint: %s) ---\n", i, n_sims, scenario, joint_model_type))
  
  t0 <- Sys.time()
  unique_seed <- base_seed + i
  
  results[[i]] <- run_simulation(
    i = i,
    n_patients = n_patients,
    lambda_c = lambda_c,
    aft_mode = aft_mode,
    max_FU = max_FU,
    unique_seed = unique_seed,
    joint_model_type = joint_model_type,
    scenario = scenario
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Trial %d completed in %.1f min\n", i, elapsed))
}

# Aggregate results
##########

runtime_summary <- do.call(rbind, lapply(results, function(r) r$runtimes))

# Runtime table
runtime_table <- data.frame(
  Component = c("LME", sprintf("%s Compile", joint_model_type), sprintf("%s Fit", joint_model_type), 
                "Total"),
  Mean_sec = c(mean(runtime_summary$lme_time), mean(runtime_summary$joint_compile_time),
               mean(runtime_summary$joint_stan_time), mean(runtime_summary$total_time)),
  SD_sec = c(sd(runtime_summary$lme_time), sd(runtime_summary$joint_compile_time),
             sd(runtime_summary$joint_stan_time), sd(runtime_summary$total_time))
)

# Extract key parameter estimates
joint_params <- do.call(rbind, lapply(results, function(r) {
  s <- r$joint_summary
  s[s$variable %in% c("beta_long[1]", "beta_long[2]", "beta_long[3]", "alpha_tilde", "beta_surv[1]", "beta_surv[2]"), 
    c("variable", "mean", "sd", "q5", "q95", "rhat")]
}))
joint_params$model <- joint_model_type

lme_coefs <- NULL
if (joint_model_type != "aft_only_BP_noRP") {
  lme_coefs <- do.call(rbind, lapply(seq_along(results), function(i) {
    data.frame(
      sim_id = i,
      param = names(results[[i]]$lme_coefs),
      estimate = as.numeric(results[[i]]$lme_coefs),
      model = "LME"
    )
  }))
}

# Parameter estimates summary
param_summary <- aggregate(cbind(mean, sd, rhat) ~ variable + model, joint_params, mean)

cat("\n=== ENZAMET Runtime Summary (seconds) ===\n")
print(runtime_table, row.names = FALSE, digits = 2)

cat("\n=== Parameter Estimates Summary ===\n")
print(param_summary, row.names = FALSE, digits = 3)

lme_summary <- NULL
if (!is.null(lme_coefs)) {
  cat("\n=== LME Fixed Effects (Mean across simulations) ===\n")
  lme_summary <- aggregate(estimate ~ param, lme_coefs, function(x) c(mean = mean(x), sd = sd(x)))
  lme_summary <- data.frame(param = lme_summary$param, 
                            mean = lme_summary$estimate[, "mean"],
                            sd = lme_summary$estimate[, "sd"])
  print(lme_summary, row.names = FALSE, digits = 3)
}
# Save results to ENZAMET output directory
output_file <- here::here(
  "Bunya", "ENZAMET", "output",
  sprintf("%s_LMM_n%d_cens%.2f_%s_scen%d_FU%d_local.rds", joint_model_type, n_patients, lambda_c, aft_mode, scenario, max_FU)
)

output_dir <- dirname(output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(list(results = results, 
             runtime_table = runtime_table, 
             param_summary = param_summary,
             lme_summary = lme_summary), 
        file = output_file)
cat(sprintf("\nResults saved to: %s\n", output_file))

##########
