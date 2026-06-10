
# Load packages
library(cmdstanr) # for use of Stan
library(lme4) # for LMM
library(nlme) # for LMM
library(survival) # for survival analysis
library(here) # for file directories

### REMOVE THIS ONCE MOVED INTO https://github.com/UQ-ULTRA/JoMoNoPH_BiomJ 
setwd(here("jomonoph_biomj"))
### 

# Set seed for reproducibility
global_seed <- 4159125
set.seed(global_seed)

# Choice to use pre-compiled CSV files or generate data for a specific seed
simulate_own_data <- TRUE # FALSE = use the pre-complied CSV files, TRUE = simulate data for a specific seed (global_seed)

## Load data
if (!simulate_own_data) {
  longitudinal_data <- read.csv("simulated_longitudinal_data.csv")
  survival_data <- read.csv("simulated_survival_data.csv")
} else {
  # Simulate data given a specific seed (global_seed)
  source("datagen.R")
  
  
}


# Load Stan model
saftjm_model <- cmdstan_model("saft_jm.stan")

# Fit longitudinal and survival submodels separately
lmm_fit <- lme(fixed = Y_obs ~ time+time:arm, random = ~ time|id, data = longitudinal_data) # Fit with LMM (linear mixed model)
surv_fit <- survreg(Surv(T_obs, status) ~ arm, data = survival_data, dist = "weibull") # Fit with AFT (survival only component)


### Construct data necessary for Stan

# ID levels for each participant, used for matching longitudinal and time-to-event observations
id_levels <- sort(unique(longitudinal_data$id))
J_1_long <- match(longitudinal_data$id, id_levels)
J_1_unique <- match(survival_data$id, id_levels)

# Covariate lists for longitudinal, survival and longitudinal-survival (used for inserting current value into time-to-event model) submodels
X_long <- cbind(Intercept = 1, time = longitudinal_data$time, time_arm = longitudinal_data$time*longitudinal_data$arm)
X_surv <- matrix(survival_data$arm, ncol = 1)
X_long_surv <- cbind(Intercept = 1, time = survival_data$T_obs, time_arm = survival_data$T_obs*survival_data$arm)


# Stan data to be inputted to the Stan sAFT model
stan_data <- list(
  K_long = ncol(X_long), # number of longitudinal population level covariates
  q = ncol(X_surv), # number of survival covariates
  m = 5, # degree of bernstein polynomial
  s_long = sd(longitudinal_data$Y_obs), # sample standard deviation of longitudinal outcome
  s_surv = sd(survival_data$T_obs), # sample standard deviation of survival scale
  alpha_tilde_sd = 1,
  N_long = nrow(longitudinal_data),
  Y_long = longitudinal_data$Y_obs,
  X_long = X_long,
  N_1_long = length(id_levels),
  J_1_long = J_1_long,
  Z_1_1_long = rep(1, nrow(longitudinal_data)),
  Z_1_2_long = longitudinal_data$time,
  n = nrow(survival_data),
  status = survival_data$status,
  time = survival_data$T_obs,
  X_surv = X_surv,
  J_1_unique = J_1_unique,
  X_long_surv = X_long_surv
)

# Fit Stan model
sAFT_fit <- saftjm_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  seed = global_seed,
  iter_warmup = 1000,
  iter_sampling = 1000
)


# Comparison of results
result_comparison <- data.frame(
  Parameter = c("beta_long_time", "beta_long_time_arm", "beta_surv_arm", "alpha_tilde"),
  LMM = c(fixef(lmm_fit)["time"], fixef(lmm_fit)["time:arm"], NA, NA),
  AFT = c(NA, NA, surv_fit$coefficients["arm"], NA),
  sAFT = c(sAFT_fit$summary(c("beta_long[2]", "beta_long[3]"))$mean, sAFT_fit$summary("beta_surv[1]")$mean, sAFT_fit$summary("alpha_tilde")$mean)
)
result_comparison


print(summary(lmm_fit)$tTable)
print(summary(surv_fit))
print(sAFT_fit$summary(c("beta_long", "beta_surv", "alpha_tilde")))


