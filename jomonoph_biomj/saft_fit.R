
# Load packages
library(cmdstanr) # for use of Stan
library(lme4) # for LMM
library(nlme) # for LMM
library(survival) # for survival analysis
library(here) # for file directories

### REMOVE THIS ONCE MOVED INTO https://github.com/UQ-ULTRA/JoMoNoPH_BiomJ 
setwd(here("jomonoph_biomj"))
### 

# Load data
longitudinal_data <- read.csv("simulated_longitudinal_data.csv")
survival_data <- read.csv("simulated_survival_data.csv")

# Load Stan model
saftjm_model <- cmdstan_model("saft_jm.stan")

# Fit with LMM (linear mixed model)
lmm_fit <- lme(fixed = Y_obs ~ time+time:arm, random = ~ time|id, data = longitudinal_data)

# Fit with AFT (survival only component)
surv_fit <- survreg(Surv(T_obs, status) ~ arm, data = survival_data, dist = "weibull")

# Construct data necessary for Stan
id_levels <- sort(unique(longitudinal_data$id))
J_1_long <- match(longitudinal_data$id, id_levels)
J_1_unique <- match(survival_data$id, id_levels)

X_long <- cbind(Intercept = 1, time = longitudinal_data$time, time_arm = longitudinal_data$time * longitudinal_data$arm)
X_surv <- matrix(survival_data$arm, ncol = 1)
X_long_surv <- cbind(Intercept = 1, time = survival_data$T_obs, time_arm = survival_data$T_obs * survival_data$arm)


# Stan data to be inputted to the Stan sAFT model
stan_data <- list(
  K_long = ncol(X_long), # number of fixed effects in the longitudinal submodel
  q = ncol(X_surv), # number of fixed effects in the survival submodel
  m = 5, # number of random effects in the longitudinal submodel (intercept and slope)
  beta_long_prior_mean = rep(0, ncol(X_long)),
  beta_long_prior_scale = rep(10, ncol(X_long)),
  beta_surv_prior_mean = rep(0, ncol(X_surv)),
  beta_surv_prior_scale = rep(10, ncol(X_surv)),
  s_long = sd(longitudinal_data$Y_obs),
  s_surv = 1,
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
  iter_warmup = 1000,
  iter_sampling = 1000
)

