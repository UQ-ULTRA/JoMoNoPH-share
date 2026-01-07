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
#                    (c) fit LME on longitudinal outcome
#                    (d) fit simple AFT model (survreg) as a staging model / prior anchor
#                    (e) fit chosen Stan joint model (BP/GB/GB_Quantile/GP)
#                    (f) return summaries + diagnostics + runtime + LOO/WAIC
# 4. SLURM job mapping: SLURM_ARRAY_TASK_ID → (config_id, batch_id) → trial_ids
# 5. Results saved as one RDS per batch; script skips work if output already exists
#
# MODEL OPTIONS (baseline hazard in the Stan joint model):
# - BP: Bernstein Polynomial baseline hazard (joint model)
# - GB: Gaussian Basis baseline hazard (joint model)
# - GB_Quantile: Gaussian Basis with quantile-based knots (joint model)
# - GP: Gaussian Process baseline hazard (joint model)
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
# - Choose baseline hazard model: joint_model_type <- "BP" / "GB" / "GB_Quantile" / "GP"
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
k_bases <- 5  # Default, will be overridden

# Joint model selection: "BP" (Bernstein Polynomial), "GB" (Gaussian Basis), "GB_Quantile", or "GP" (Gaussian Process)
joint_model_type <- "GB"  # Change to "BP", "GB", "GB_Quantile", or "GP" as needed

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
    
    # Priors for gamma scaling only â all others left unstandardized
    location_beta_real = 0,
    scale_beta_real = 1,
    std_real = 1,
    means_real = 0,
    
    # Smoothing parameters
    r = 2,        # second-order difference penalty
    tau_h = 0.01  # smoothing parameter
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
    D = D,
    beta_0 = beta_0,
    beta_1 = beta_1,
    beta_2 = beta_2,
    sigma_e = sigma_e,
    log_HR = log_HR,
    log_AF = log_AF,
    alpha_PH = alpha_PH,
    alpha_AFT = alpha_AFT,
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
  fit_long <- nlme::lme(
    fixed = Y_obs ~ time + time_by_arm,
    random = ~ time | id,
    data = sim_dat$longitudinal,
    na.action = na.omit,
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
  )
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
  
  # Prepare data for joint model (BP, GB, GB_Quantile, or GP)
  if (joint_model_type == "BP") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist.stan"
  } else if (joint_model_type == "GB") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis-JM-Hist.stan"
  } else if (joint_model_type == "GB_Quantile") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Basis-JM-Hist-Quantile.stan"
  } else if (joint_model_type == "GP") {
    stan_model_file_joint <- "~/JoMoNoPH-share/Stan/Gaussian-Process-JM-Hist.stan"
  } else {
    stop(sprintf("Unknown joint_model_type: %s. Use 'BP', 'GB', 'GB_Quantile', or 'GP'.", joint_model_type))
  }
  
  stan_data_joint <- convert_data_from_models(fit_long, 
                                              fit_surv, 
                                              surv_data = sim_dat$survival,
                                              k_bases = k_bases)
  
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
  
  # ------
  # Generate initial values for joint model
  init_values_joint <- init_fun_joint(sim_dat = sim_dat, stan_data = stan_data_joint, chains = chains, joint_model_type = joint_model_type)
  
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
  
  diag_joint <- fit_joint$diagnostic_summary()
  diag_joint$sim_id <- i
  diag_joint$model <- joint_model_type
  
  # LME coefficients for comparison
  coefs <- nlme::fixef(fit_long)
  
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
  
  return(list(
    joint_summary = summ_joint,
    joint_diagnostics = diag_joint,
    lme_coefs = coefs,
    runtimes = runtimes,
    loo_joint = loo_joint_est,
    waic_joint = waic_joint_est))
}

##########


# SLURM Job Mapping Logic
##########

# Read job ID from environment
sim_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
if (sim_id == 0) stop("SLURM_ARRAY_TASK_ID not found.")

# Define configuration grid for ENZAMET
# ENZAMET uses n_patients=1100, max_FU=72, loglogistic distribution
config_base <- expand.grid(
  n_patients = c(1100),
  aft_mode   = c("Weibull"),
  k_bases    = c(5),
  # scenario   = c(1, 2, 3, 4, 5, 6, 7),  # Data generation scenarios
  scenario   = c(1, 2, 3, 4, 5),  # Data generation scenarios
  max_FU     = 72,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

lambda_map <- c(weibull = 0.023, loglogistic = 0.013) # under scenario 1, give ~50% censoring

cfg_default <- transform(config_base,
                         lambda_c = lambda_map[tolower(aft_mode)])

cfg_zero <- transform(config_base,
                      lambda_c = 0)

config_grid <- rbind(cfg_default, cfg_zero)
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

# Define trial numbers for this batch
trial_ids <- ((batch_id - 1) * batch_size + 1):(batch_id * batch_size)

# Hash the config (e.g., using integers from factors)
config_hash <- as.integer(as.numeric(as.factor(
  paste(n_patients, lambda_c, aft_mode, k_bases, scenario, max_FU, sep = "_")
)))

# label for censoring used in filename & messages
cens_label <- if (lambda_c == 0) "0" else "50"

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
