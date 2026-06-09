

# Pre-processing
##########

rm(list = ls())

library(tidyverse)
library(survival)
library(MASS)
library(nlme)
library(simsurv)
library(cmdstanr)
library(JMbayes2)
library(here)

# Auto-detect CmdStan: prefer pinned version, fall back to latest installed
# Uses cmdstan_default_install_path() so the path resolves correctly on both
# Windows (USERPROFILE/.cmdstan) and Linux/Mac (~/.cmdstan)
.find_cmdstan <- function() {
  base      <- cmdstanr::cmdstan_default_install_path()
  preferred <- file.path(base, "cmdstan-2.37.0")
  if (dir.exists(preferred)) return(preferred)
  avail <- sort(list.dirs(base, recursive = FALSE, full.names = TRUE))
  if (length(avail) > 0) return(avail[length(avail)])
  stop("No CmdStan installation found in: ", base)
}
set_cmdstan_path(.find_cmdstan())

source(here("Bunya", "Auxiliary", "BP-functions.R"))

##########



# Settings
##########

# Want to look over three 'problematic' scenarios
# 1) LL 1.20 shape 23.0 scale, scenario 1
# 2) WB 0.90 shape 38.0 scale, scenario 3
# 3) WB 0.90 shape 38.0 scale, scenario 5


# Edit these to change the scenario, sample size, and sampler controls


scenario <- 1 # ENZAMET data-gen scenario (1, 2, 4, ...)
n_patients <- 1100 # smaller than ENZAMET full (1100) for local speed
max_FU <- 120
lambda_c <- 0.009948730 # scenario-1 censoring rate; set 0 for no censoring
k_bases <- 5  # Bernstein basis dimension
chains <- 4
iter_warmup <- 1000
iter_sampling <- 1000

# config_hash
config_hash <- 1

# Construct unique_seed from config_hash
unique_seed <- 1e6 * config_hash + 1  
seed <- unique_seed

# lambda_c = c(
#   0.009948730, 0.003906250, 0.004272461, 0.017028809, 0.017211914, # loglogistic(1.20, 23.0), 50% censoring
#   0.009643555, 0.003479004, 0.003784180, 0.017211914, 0.017456055, # Weibull(0.90, 38.0) # administrative censoring
#   0.008666992, 0.002227783, 0.002593994, 0.015563965, 0.015716553  # Weibull(1.30, 38.0) # administrative censoring
# )

##########



# Helper functions
##########

### Data conversion: lme() + survreg() fits -> Stan data list

convert_data_from_models <- function(fit_long, fit_surv, surv_data, k_bases = 5) {
  mf_long_raw <- nlme::getData(fit_long)
  mf_long <- model.frame(formula(fit_long), data = mf_long_raw, na.action = na.omit)
  idx_long <- as.numeric(rownames(mf_long))
  mf_long$id <- mf_long_raw$id[idx_long]
  mf_long$time <- mf_long_raw$time[idx_long]
  mf_long$arm <- mf_long_raw$arm[idx_long]
  
  X_long   <- model.matrix(~ time + time_by_arm, data = mf_long)
  Y_long   <- as.numeric(model.response(mf_long))
  id_long  <- as.factor(mf_long$id)
  J_1_long <- as.integer(id_long)
  N_1_long <- length(unique(J_1_long))
  
  mf_surv <- model.frame(fit_surv)
  idx_surv <- as.numeric(rownames(mf_surv))
  mf_surv$id <- surv_data$id[idx_surv]
  mf_surv$arm <- surv_data$arm[idx_surv]
  y_surv <- model.response(mf_surv)
  mf_surv$time <- y_surv[, "time"]
  mf_surv$status <- y_surv[, "status"]
  mf_surv$arm <- as.numeric(mf_surv$arm)
  X_surv <- model.matrix(~ 0 + arm, data = mf_surv)
  J_1_unique <- match(unique(mf_surv$id), levels(id_long))
  
  mf_surv$time_by_arm <- mf_surv$time * mf_surv$arm
  X_long_surv <- model.matrix(~ time + time_by_arm, data = mf_surv)
  
  beta_long_mean <- nlme::fixef(fit_long)
  beta_surv_mean <- coef(fit_surv)[-1]
  vc_long <- VarCorr(fit_long)
  sd_b <- sqrt(as.numeric(vc_long[c("(Intercept)", "time"), "Variance"]))
  s_long <- summary(fit_long)$sigma
  s_surv <- summary(fit_surv)$scale
  alpha_tilde_sd <- (log(2) / 1.96) * (s_long / s_surv)
  
  list(
    N_long = nrow(X_long), K_long = ncol(X_long),
    Y_long = Y_long, X_long = X_long,
    J_1_long = J_1_long, N_1_long = N_1_long,
    Z_1_1_long = rep(1, length(Y_long)),
    Z_1_2_long = mf_long$time,
    n = nrow(X_surv), q = ncol(X_surv),
    status = as.numeric(mf_surv$status),
    time   = mf_surv$time,
    X_surv = X_surv,
    m = k_bases,
    X_long_surv = X_long_surv,
    J_1_unique  = J_1_unique,
    beta_long_prior_mean  = beta_long_mean,
    beta_long_prior_scale = rep(10, length(beta_long_mean)),
    beta_surv_prior_mean  = beta_surv_mean,
    beta_surv_prior_scale = rep(10, length(beta_surv_mean)),
    sd_b_prior_mean  = sd_b,
    sd_b_prior_scale = rep(5, length(sd_b)),
    s_long = s_long, s_surv = s_surv,
    alpha_tilde_sd = alpha_tilde_sd,
    location_beta_real = 0, scale_beta_real = 1,
    std_real = 1, means_real = 0
  )
}


### Chain initialisation: overdispersed starting values across chains

make_init <- function(sim_dat, stan_data, chains) {
  lmm <- nlme::lme(
    fixed = Y_obs ~ time + time_by_arm, random = ~ time | id,
    data = sim_dat$longitudinal,
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
  )
  aft <- survival::survreg(Surv(time, status) ~ arm, data = sim_dat$survival)
  b_fe <- as.numeric(nlme::fixef(lmm))
  sig <- summary(lmm)$sigma
  vc <- nlme::VarCorr(lmm)
  sds <- sqrt(as.numeric(vc[1:2, "Variance"]))
  b_sv <- rep(coef(aft)[2], stan_data$q)
  rho <- 0.2
  L <- t(chol(matrix(c(1, rho, rho, 1), 2, 2)))
  
  lapply(seq_len(chains), function(.) {
    list(
      z_1_long = matrix(rnorm(2 * stan_data$N_1_long), nrow = 2),
      L_1_long = L,
      beta_long = b_fe + rnorm(stan_data$K_long, 0, 0.1),
      sigma_long = abs(rnorm(1, sig, 1)),
      sd_1_long = pmax(sds + abs(rnorm(2, 0, 0.05)), 1e-3),
      beta_surv = rnorm(stan_data$q, b_sv, 0.05),
      alpha_tilde = rnorm(1, 0, 0.1),
      gamma = abs(rnorm(stan_data$m, 0.1, 0.01))
    )
  })
}


### Stan model fit wrapper: compile, sample, return fit + elapsed time

fit_stan_model <- function(stan_file, stan_data, init_vals, seed,
                           chains, iter_warmup, iter_sampling) {
  mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE),
                       force_recompile = FALSE)
  t0  <- Sys.time()
  fit <- mod$sample(
    data              = stan_data,
    seed              = seed,
    init              = init_vals,
    chains            = chains,
    parallel_chains   = chains,
    threads_per_chain = 1,
    iter_warmup       = iter_warmup,
    iter_sampling     = iter_sampling,
    max_treedepth     = 10,
    refresh           = 100
  )
  list(fit = fit, elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

##########



# Simulate data
##########

set.seed(seed)
source(here("Bunya", "ENZAMET", sprintf("Data-Gen-ENZAMET-%02d.R", scenario)))

sim_dat <- simulate_joint_dataset(
  D = D, beta_0 = beta_0, beta_1 = beta_1, beta_2 = beta_2, sigma_e = sigma_e,
  log_HR = log_HR, log_AF = log_AF,
  alpha_PH = alpha_PH, alpha_AFT = alpha_AFT,
  weibull_shape = weibull_shape, weibull_scale = weibull_scale,
  loglogistic_shape = loglogistic_shape, loglogistic_scale = loglogistic_scale,
  visit = visit, seed = seed,
  n_patients = n_patients, max_FU = max_FU, lambda_c = lambda_c,
  aft_mode = "loglogistic", link_type = "value"
)

# Save this data to a CSV file
write.csv(sim_dat$longitudinal, file = here("Bunya", "ENZAMET", sprintf("simulated_longitudinal_data.csv")), row.names = FALSE)
write.csv(sim_dat$survival, file = here("Bunya", "ENZAMET", sprintf("simulated_survival_data.csv")), row.names = FALSE)



# Center longitudinal outcome and keep the sample mean for back-transformation later.
Y_obs_sample_mean <- mean(sim_dat$longitudinal$Y_obs, na.rm = TRUE)
sim_dat$longitudinal <- sim_dat$longitudinal %>%
  mutate(
    Y_obs = Y_obs - Y_obs_sample_mean,
    time_by_arm = time * arm
  )
sim_dat$Y_obs_sample_mean <- Y_obs_sample_mean

sim_dat$survival <- sim_dat$survival %>% dplyr::rename(time = T_obs)

message(sprintf("Censoring proportion: %.1f%%", 100*mean(sim_dat$survival$status == 0)))

##########



# Model fitting
##########

### Model 1: LMM (random intercept + slope)

t0 <- Sys.time()
fit_lmm <- nlme::lme(
  fixed     = Y_obs ~ time + time_by_arm,
  random    = ~ time | id,
  data      = sim_dat$longitudinal,
  na.action = na.omit,
  control   = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
)
lmm_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf("LMM fitted in %.1f s", lmm_time))
print(summary(fit_lmm)$tTable)


### Staging survival model (prior anchor only; not a model of interest)
# Used solely to anchor the Stan priors and build the shared data list

fit_surv_stage <- survival::survreg(
  Surv(time, status) ~ arm, data = sim_dat$survival,
  dist = "exponential", x = TRUE
)

stan_data <- convert_data_from_models(fit_lmm, fit_surv_stage, sim_dat$survival, k_bases)
init_vals <- make_init(sim_dat, stan_data, chains)

### Model 2: Cox PH joint model (JMbayes2)

fit_cox <- survival::coxph(Surv(time, status) ~ arm, data = sim_dat$survival, x = TRUE)

t0 <- Sys.time()
fit_jmb2 <- suppressMessages(
  JMbayes2::jm(fit_cox, list(fit_lmm), time_var = "time",
               n_iter = 3500L, n_burnin = 500L, n_chains = 3L)
)
jmb2_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf("JMbayes2 fitted in %.1f s", jmb2_time))
print(summary(fit_jmb2))

### Model 3: BP (Bernstein Polynomial baseline hazard)

# res_bp <- fit_stan_model(
#   stan_file = here("Stan", "Bernstein-Polynomials-JM-Hist.stan"),
#   stan_data = stan_data,
#   init_vals = init_vals,
#   seed = seed,
#   chains = chains,
#   iter_warmup = iter_warmup,
#   iter_sampling = iter_sampling
# )
# message(sprintf("BP fitted in %.1f s", res_bp$elapsed))
# 
# summ_bp <- res_bp$fit$summary(NULL, posterior::default_summary_measures(),
#                               quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)))
# summ_bp$model <- "BP"
# print(res_bp$fit$diagnostic_summary())


### Model 4: BP_Penalty (BP with second-order difference penalty)

res_bpp <- fit_stan_model(
  stan_file     = here("Stan", "Bernstein-Polynomials-JM-Hist-Penalty.stan"),
  stan_data     = stan_data,
  init_vals     = init_vals,
  seed          = seed,
  chains        = chains,
  iter_warmup   = iter_warmup,
  iter_sampling = iter_sampling
)
message(sprintf("BP_Penalty fitted in %.1f s", res_bpp$elapsed))

summ_bpp       <- res_bpp$fit$summary(NULL, posterior::default_summary_measures(),
                                      quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)))
summ_bpp$model <- "BP_Penalty"
print(res_bpp$fit$diagnostic_summary())

### Model 5: 

res_bp_2 <- fit_stan_model(
  stan_file = here("Stan", "Bernstein-Polynomials-JM-Hist.stan"),
  stan_data = stan_data,
  init_vals = init_vals,
  seed = seed,
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling
)
message(sprintf("BP fitted in %.1f s", res_bp$elapsed))

summ_bp <- res_bp$fit$summary(NULL, posterior::default_summary_measures(),
                              quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)))
summ_bp$model <- "BP"
print(res_bp$fit$diagnostic_summary())


##########



# Results
##########

params_of_interest <- c("alpha_tilde", "beta_surv[1]",
                        "beta_long[1]", "beta_long[2]", "beta_long[3]")

### LMM fixed effects
cat("\nLMM fixed effects\n")
print(round(nlme::fixef(fit_lmm), 4))

### BP posterior summary (key parameters)
cat("\nBP posterior summary\n")
print(summ_bp %>%
        filter(variable %in% params_of_interest) %>%
        select(variable, mean, median, sd, q2.5, q97.5, rhat, ess_bulk))

### BP_Penalty posterior summary (key parameters)
# cat("\nBP_Penalty posterior summary\n")
# print(summ_bpp %>%
#         filter(variable %in% params_of_interest) %>%
#         select(variable, mean, median, sd, q2.5, q97.5, rhat, ess_bulk))

### JMbayes2 survival and association summaries
s <- summary(fit_jmb2)
cat("\nJMbayes2 survival summary\n");    print(s$Survival)
cat("\nJMbayes2 association summary\n"); print(s$Associations)

### Runtimes
print(data.frame(
  model   = c("LMM", "BP (Stan)", "BP_Penalty (Stan)", "JMbayes2 Cox PH"),
  seconds = round(c(lmm_time, res_bp$elapsed, res_bpp$elapsed, jmb2_time), 1)
))


##########
