# ------------------------------------------------------------------------------
# Script Title: Bayesian Joint Modelling (Bernstein-Polynomials)
# Description :
#   This script performs simulation and Bayesian inference for a joint 
#   longitudinal-survival model. The steps include:
#     - Compiling a Stan model (e.g., `Bernstein-Polynomials-JM-Hist.stan`)
#     - Setting model and prior parameters
#     - Simulating synthetic joint data (longitudinal + survival)
#     - Converting data to Stan-compatible format
#     - Fitting the model using `cmdstanr`
#     - Extracting and summarizing posterior draws
#     - Relabelling and presenting results for interpretation 
# ------------------------------------------------------------------------------
rm(list = ls())


library(tidyverse)
library(survival)
library(MASS)
library(future)
library(furrr)
library(simsurv)
library(rstan)
library(brms)
library(cmdstanr)
library(R.utils)

# Set CmdStan path globally for main session
set_cmdstan_path("~/cmdstan/cmdstan-2.36.0")

# Load helper functions and data generator
source("~/JoMoNoPH/Bunya/Auxiliary/BP-functions.R")
# source("~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET.R") # ENZAMET


scenario_tbl <- tibble::tribble(
  ~scenario_id, ~scenario_name, ~data_gen_file,                                    ~out_dir,
  1,            "scenario-1",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-01.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-1",
  2,            "scenario-2",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-02.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-2",
  3,            "scenario-3",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-03.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-3",
  4,            "scenario-4",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-04.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-4",
  5,            "scenario-5",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-05.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-5",
  6,            "scenario-6",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-06.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-6",
  7,            "scenario-7",   "~/JoMoNoPH/Bunya/ENZAMET/Data-Gen-ENZAMET-07.R",  "~/JoMoNoPH/Bunya/ENZAMET/output/scenario-7"
)




# Set model parameters
chains <- 4
threads_per_chain <- 1
parallel_chains <- chains
cpus_per_sim <- chains * threads_per_chain

degree <- 5

# Batch structure is set later in the script based on SLURM logic:

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
  
  vc_long <- tryCatch(VarCorr(fit_long), error = function(e) NULL)
  sd_b <- if (!is.null(vc_long)) {
    sqrt(as.numeric(vc_long[c("(Intercept)", "time"), "Variance"]))
  } else rep(1.0, 2)
  
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
    
    # Smoothing parameters for r-th difference penalty
    r = 2,          # Second-order differences (penalize curvature)
    tau_h = 1.0,    # Smoothing parameter (moderate smoothing)
    
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
    means_real = 0
  )
}


# --- Simulation wrapper function for one iteration ---
run_simulation <- function(i, n_patients, lambda_c, aft_mode, max_FU, config_hash) {
  tryCatch({
  set_cmdstan_path("~/cmdstan/cmdstan-2.36.0")

# New (globally unique):
unique_seed <- 1e6 * config_hash + i    
      
# Generate data
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
    weibull_shape= weibull_shape,
    weibull_scale= weibull_scale,
    loglogistic_shape = loglogistic_shape, 
    loglogistic_scale = loglogistic_scale,
    visit = visit,
    seed = unique_seed,
    n_patients = n_patients,
    max_FU = max_FU,
    lambda_c = lambda_c,
    aft_mode =  aft_mode,  # "PH" or "Weibull" or "loglogistic"
    link_type = "value" # "value" or "slope"
  )

# ------------------------------------------------------------------------------ 
# Fit longitudinal model on raw variables 
# ------------------------------------------------------------------------------ 
# Augment longitudinal data with time-by-arm interaction
sim_dat$longitudinal <- sim_dat$longitudinal %>%
  mutate(time_by_arm = time * arm)

fit_long <- tryCatch({
  nlme::lme(
    fixed = Y_obs ~ time + time_by_arm,
    random = ~ time | id,
    data = sim_dat$longitudinal,
    na.action = na.omit,
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
  )
}, error = function(e) {
  message(sprintf("Skipping simulation %d: lme() failed — %s", i, e$message))
  return(NULL)
})

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


# Convert all inputs to Stan data format
stan_data <- convert_data_from_models(fit_long, 
                                      fit_surv, 
                                      surv_data = sim_dat$survival,
                                      degree = 5)



# Initialisation function
init_fun <- function(sim_dat, stan_data, chains = 4) {
  # --- Fallback values ---
    # These are used if the lme() or survreg() model fails to converge or throws an error.
    # They are deliberately safe, neutral guesses meant to allow Stan to initialize without biasing the sampler.
  fallback <- list(
    beta_long = rep(0, stan_data$K_long),     # Intercept, time, time_by_arm
    sigma_long = 10,                          # Reasonable residual SD
    sd_1_long = c(1, 1),                      # Conservative random effects SDs
    beta_surv = rep(0, stan_data$q),          # survival covariates
    gamma_scaled = rep(0.1, stan_data$m),     # Small positive Bernstein weights
    alpha_tilde = 0                                 # No association
  )
  
  # --- Estimate longitudinal model using lme() ---
    # These provide empirical starting values for beta_long, sigma_long, and sd_1_long
  LMM <- tryCatch(
    nlme::lme(
      fixed = Y_obs ~ time + time_by_arm,
      random = ~ time | id,
      data = sim_dat$longitudinal,
      control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
    ),
    error = function(e) NULL
  )
  
  if (!is.null(LMM)) {
      # Fixed effects (intercept and covariate effects)
      fallback$beta_long <- tryCatch(as.numeric(nlme::fixef(LMM)), error = function(e) fallback$beta_long)
      
      # Residual error SD
      fallback$sigma_long <- tryCatch(summary(LMM)$sigma, error = function(e) fallback$sigma_long)
      
      # Random intercept/slope SDs
      vc <- tryCatch(nlme::VarCorr(LMM), error = function(e) NULL)
      if (!is.null(vc)) {
        fallback$sd_1_long <- tryCatch(
          sqrt(as.numeric(vc[1:2, "Variance"])),  # First two rows are intercept and slope variances
          error = function(e) fallback$sd_1_long
        )
      }
    }
    
  # --- Estimate survival treatment effect using survreg() ---
    # This is used to initialise MCMC sampling)
    AFT <- tryCatch(
      survival::survreg(Surv(time, status) ~ arm, data = sim_dat$survival),
      error = function(e) NULL
    )
    
    if (!is.null(AFT)) {
    fallback$beta_surv <- rep(coef(AFT)[2], stan_data$q)
    }
    
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
      

      list(
      # Latent random effects
        z_1_long = matrix(rnorm(2 * stan_data$N_1_long), nrow = 2),
        
        L_1_long = L_1_long_init,
        
      # Longitudinal parameters
      beta_long = fallback$beta_long + rnorm(stan_data$K_long, 0, 0.1),
      sigma_long = abs(rnorm(1, fallback$sigma_long, 1)),
      sd_1_long = sd_jittered,
      
        # Survival treatment effect from survreg(), jittered
      beta_surv = rnorm(stan_data$q, fallback$beta_surv, 0.05),
        
        # Positive weights for spline-based hazard: jittered around small positive default
      gamma = abs(rnorm(stan_data$m, 0.1, 0.01)),
      
      # Association
      alpha_tilde = rnorm(1, 0, 0.1)
      )
    })
    
    return(init_list)
  }
  
extract_lmm_components <- function(fit_long) {
  stopifnot(inherits(fit_long, "lme"))
  
  # Fixed effects
  fe <- nlme::fixef(fit_long)
  
  # Residual SD
  sigma_eps <- summary(fit_long)$sigma
  
  # Random effects SDs (intercept + time slope)
  vc <- nlme::VarCorr(fit_long)
  
  # VarCorr() prints like a matrix; safest is to pull StdDev by row name
  sd_b0 <- as.numeric(vc["(Intercept)", "StdDev"])
  sd_b1 <- as.numeric(vc["time",        "StdDev"])
  
  c(
    fe,
    "sd_1_long[1]" = sd_b0,
    "sd_1_long[2]" = sd_b1,
    "sigma_long"   = sigma_eps
  )
}


  # ------
# Generate initial values for Stan chains
init_values <- init_fun(sim_dat = sim_dat, stan_data = stan_data, chains = chains)
  
  mod <- cmdstan_model("~/JoMoNoPH/Stan/Bernstein-Polynomials-JM-Hist.stan",
                       cpp_options = list(stan_threads = TRUE),
                       force_recompile = FALSE)
  
  Sys.chmod(mod$exe_file(), mode = "0755")
  if (!file.exists(mod$exe_file())) {
    stop("Stan model compilation failed; binary not found.")
  }
  
  # --- Fit model ---
  fit <- mod$sample(
    data = stan_data,
    seed = unique_seed,
    init = init_values,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = 1000, # 100
    iter_sampling = 1000, # 100
    max_treedepth = 10)
  
    vars <- c(
      "beta_long[1]", "beta_long[2]", "beta_long[3]",
      "sd_1_long[1]", "sd_1_long[2]", "sigma_long",
      "alpha", "beta_surv[1]",
      paste0("gamma[", 1:stan_data$m, "]")
    )
  
    summ <- fit$summary(variables = vars) %>%
      dplyr::select(variable, mean, sd, rhat, ess_bulk, ess_tail)
    
    summ$sim_id <- i
  

    CrI <- posterior::summarise_draws(
      fit$draws(variables = vars),
      "mean", "median",
      ~quantile(.x, probs = c(0.025, 0.975))
    )
    CrI$sim_id <- i
     
    diag <- fit$diagnostic_summary()
    diag$sim_id <- i
      
    coefs <- tryCatch(
      {
        out <- extract_lmm_components(fit_long)
        out <- as.list(out)
        out$sim_id <- i
        out
      },
      error = function(e) {
        message(sprintf("Skipping LMM coefs for sim %d: %s", i, e$message))
        list(sim_id = i)
      }
    )
    
    
  return(list(
    CrI = CrI,
    summary = summ,
    diagnostics = diag,
    fit_long = coefs))
  
  
  }, error = function(e) {
    message(sprintf("Simulation %d failed — %s", i, e$message))
    return(NULL)
  })
}


# === GLOBAL SETTINGS FOR CONFIG GRID + BATCHING (define ONCE) =================



# =============================================================================
# SLURM job mapping + batching
# Goal: 1 SLURM array task runs exactly ONE scenario × ONE config × ONE batch.
# Within that batch, we parallelise over trial_ids using future_map().
# =============================================================================

# --- Build config grid ONCE (common across scenarios) -------------------------
# Each row is one "config" (e.g., censoring mechanism / lambda_c setting).
config_base <- expand.grid(
  n_patients = c(1100), #110
  aft_mode   = c("Weibull", "loglogistic"),
  max_FU     = 60,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Map AFT baseline type to the exponential censoring rate used in that setting
lambda_map <- c(weibull = 0.023, loglogistic = 0.013)

# Default censoring (your calibrated lambda_c)
cfg_default <- transform(config_base,
                         lambda_c = lambda_map[tolower(aft_mode)])

# Additional config with no independent exponential censoring
cfg_zero <- transform(config_base,
                      lambda_c = 0)

# Final grid of configs (rows are configs)
config_grid <- rbind(cfg_default, cfg_zero)
row.names(config_grid) <- NULL


# --- Batching within a config ------------------------------------------------
batch_size <- 10
sims_per_config <- 100
batches_per_config <- sims_per_config / batch_size
stopifnot(batches_per_config == as.integer(batches_per_config))  # safety check

# === DERIVED TASK COUNTS (DEFINE ONCE) ========================================
n_scenarios <- nrow(scenario_tbl)
n_configs   <- nrow(config_grid)

tasks_per_scenario <- n_configs * batches_per_config
total_tasks <- n_scenarios * tasks_per_scenario



# --- Decode SLURM array index ------------------------------------------------
# SLURM_ARRAY_TASK_ID enumerates all scenario × config × batch combinations.
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
if (task_id == 0) stop("SLURM_ARRAY_TASK_ID not found.")

array_offset <- as.integer(Sys.getenv("ARRAY_OFFSET", unset = "0"))
global_task_id <- task_id + array_offset

if (global_task_id > total_tasks) {
  stop(sprintf("global_task_id=%d exceeds total_tasks=%d.", global_task_id, total_tasks))
}

scenario_index <- ((global_task_id - 1) %/% tasks_per_scenario) + 1
task_within_scenario <- ((global_task_id - 1) %% tasks_per_scenario) + 1

config_id <- ((task_within_scenario - 1) %/% batches_per_config) + 1
batch_id  <- ((task_within_scenario - 1) %% batches_per_config) + 1



# --- Load scenario-specific DGM parameters -----------------------------------
# IMPORTANT: Only one scenario's parameters are sourced per SLURM task.
scenario_id   <- scenario_tbl$scenario_id[scenario_index]
scenario_name <- scenario_tbl$scenario_name[scenario_index]
out_dir       <- scenario_tbl$out_dir[scenario_index]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf(
  "\n====================\nRunning %s (task %d/%d)\nconfig_id=%d/%d, batch_id=%d/%d\n====================\n",
  scenario_name, global_task_id, total_tasks,
  config_id, n_configs,
  batch_id, batches_per_config
))


# Remove previous scenario params (defensive, in case of re-use)
rm(list = c("D","beta_0","beta_1","beta_2","sigma_e","log_HR","log_AF",
            "alpha_PH","alpha_AFT","weibull_shape","weibull_scale",
            "loglogistic_shape","loglogistic_scale","visit",
            "n_patients","lambda_c","max_FU"),
   envir = .GlobalEnv)

# Source the scenario-specific data generating mechanism file
source(scenario_tbl$data_gen_file[scenario_index])


config <- config_grid[config_id, ]
aft_mode <- config$aft_mode

if (scenario_id == 7) {
  if (!exists("n_patients") || !exists("lambda_c") || !exists("max_FU")) {
    stop("Scenario 7 requires n_patients, lambda_c, and max_FU to be defined in the DGM file.")
  }
  message(sprintf(
    "Scenario 7 override: using DGM n_patients=%s, lambda_c=%s, max_FU=%s",
    n_patients, lambda_c, max_FU
  ))
} else {
  n_patients <- config$n_patients
  lambda_c   <- config$lambda_c
  max_FU     <- config$max_FU
}


# Trial IDs for this batch (e.g., batch 1 -> 1:10, batch 2 -> 11:20, ...)
trial_ids <- ((batch_id - 1) * batch_size + 1):(batch_id * batch_size)

# Hash ensures unique seeds across scenario/config/batch while remaining reproducible
config_hash <- as.integer(as.numeric(as.factor(
  paste(scenario_id, n_patients, lambda_c, aft_mode, max_FU, sep = "_")
)))

message(sprintf(
  "Running batch %d (config_id %d): n=%d, lambda_c=%.4f, aft_mode=%s, max_FU=%d",
  batch_id, config_id, n_patients, lambda_c, aft_mode, max_FU
))


# --- Parallel execution within the batch -------------------------------------
# Each Stan fit uses cpus_per_sim cores (typically = chains × threads_per_chain).
# We therefore allocate workers = floor(total_cpus / cpus_per_sim).
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
workers <- floor(cpus / cpus_per_sim)
workers <- max(workers, 1)
workers <- min(workers, length(trial_ids))

plan(multisession, workers = workers)

# NOTE: This timeout is PER TRIAL (not per batch). Change if you want batch-level.
safe_elapsed <- 10800  # 3 hours elapsed limit per trial


# Path where this batch would be written
out_file <- file.path(
  out_dir,
  sprintf(
    "BP_%s_n%d_cens%.4f_%s_FU%d_batch%03d.rds",
    scenario_name, n_patients, lambda_c, aft_mode, max_FU, batch_id
  )
)

# Skip if already completed
if (file.exists(out_file)) {
  message("Batch already exists, skipping: ", out_file)
  quit(save = "no", status = 0)
}




results <- future_map(
  trial_ids,
  ~ tryCatch({
    t0 <- Sys.time()
    
    ans <- tryCatch({
      # Ensure the time limit resets even if errors occur
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
      
      # Elapsed-only timeout for the whole run_simulation() call
      setTimeLimit(elapsed = safe_elapsed, cpu = Inf, transient = TRUE)
      
      run_simulation(.x, n_patients, lambda_c, aft_mode, max_FU, config_hash)
    }, error = function(e) {
      # Convert elapsed timeout into a NULL result; rethrow other errors
      if (grepl("reached elapsed time limit", conditionMessage(e))) NULL else stop(e)
    })
    
    cat(sprintf("Trial %d elapsed: %.1f min\n",
                .x, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    ans
  },
  error = function(e) {
    message(sprintf("Simulation %d failed or timed out: %s", .x, e$message))
    NULL
  }),
  .options = furrr_options(seed = TRUE)
)

# --- Save output for this scenario/config/batch -------------------------------
saveRDS(results, file = file.path(
  out_dir,
  sprintf("BP_%s_n%d_cens%.4f_%s_FU%d_batch%03d.rds",
          scenario_name, n_patients, lambda_c, aft_mode, max_FU, batch_id)
))
