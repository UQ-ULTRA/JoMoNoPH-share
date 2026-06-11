
# Load packages
library(cmdstanr) # for use of Stan
library(lme4) # for LMM
library(nlme) # for LMM
library(survival) # for survival analysis
library(here) # for file directories
library(tidyverse) # for data manipulation and plotting

source("datagen.R")  # Source file for data generation and MCMC initialisation functions

# Set seed for reproducibility
global_seed <- 854098 # Change to whichever seed one wants simulation for (854098)
set.seed(global_seed)

# Choice to use pre-compiled CSV files or generate data for a specific seed
simulate_own_data <- FALSE # FALSE = use the pre-complied CSV files, TRUE = simulate data for a specific seed (global_seed)

#### Scenarios used within manuscript

## By default, the simulate_joint_dataset function simulates treatment 'scenario 4' (beta_2 = 0.04, gamma (log_AF)= -0.9) under the loglogistic setting (aft_mode="loglogistic") and administrative censoring only (lambda_c=0)
## To simulate over different scenarios, true data generation distributions and censoring proportions as described within the main manuscript, input the following parameters into 'simulate_joint_dataset'

# Treatment effect scenarios: Scenario 1 (beta_2=0, log_AF=0), Scenario 2 (beta_2=0.04, log_AF=0.9), Scenario 3 (beta_2=-0.04, log_AF=0.9), Scenario 4 (beta_2=0.04, log_AF=-0.9), Scenario 5 (beta_2=-0.04, log_AF=-0.9), 
# Baseline hazard distributions: Log-logistic (aft_mode="loglogistic"), Weibull with shape=0.90 (aft_mode="weibull",weibull_shape=0.9), Weibull with shape=1.30 (aft_mode="weibull",weibull_shape=1.3), Weibull with shape=2.10 (aft_mode="weibull",weibull_shape=2.1)

# Censoring proportion (if using 50% censoring) is treatment scenario and baseline hazard distribution dependent. If only using administrative censoring, use 'lambda_c=-1', otherwise, use the following 'lambda_c = XX' values for your specific treatment + baseline hazard scenario

# Log-logistic: 0.009948730 (scen 1), 0.003906250 (scen 2), 0.004272461 (scen 3), 0.017028809 (scen 4), 0.017211914 (scen 5) 
# Weibull (shape 0.9): 0.009643555 (scen 1), 0.003479004 (scen 2), 0.003784180 (scen 3), 0.017211914, (scen 4) 0.017456055 (scen 5) 
# Weibull (shape 1.3): 0.008666992 (scen 1), 0.002227783 (scen 2), 0.002593994 (scen 3), 0.015563965 (scen 4), 0.015716553 (scen 5)
# Weibull (shape 2.1): 0.008056641 (scen 1), 0.001098633 (scen 2), 0.001464844 (scen 3), 0.014099121 (scen 4), 0.014160156 (scen 5)




## Load data
if (!simulate_own_data) {
  # Use pre-generated CSV files. These correspond to data generation for scenario 4 under global_seed=1
  longitudinal_data <- read.csv("simulated_longitudinal_data.csv")
  survival_data <- read.csv("simulated_survival_data.csv")
} else {
  # Simulate data given a specific seed (global_seed)
  sim_data <- simulate_joint_dataset(seed = global_seed, beta_2 = 0.04, log_AF = -0.90, aft_mode = "loglogistic", lambda_c = -1) # Default settings use 'scenario 4' for the loglogistic setting with administrative censoring only 
  longitudinal_data <- sim_data$longitudinal
  survival_data <- sim_data$survival
}

Y_obs_mean <- mean(longitudinal_data$Y_obs)
longitudinal_data$Y_obs_centred <- longitudinal_data$Y_obs - Y_obs_mean

# Load Stan model
saftjm_model <- cmdstan_model("saft_jm.stan")

# Fit longitudinal and survival submodels separately
lmm_fit <- lme(fixed = Y_obs_centred ~ time+time:arm, random = ~ time|id, data = longitudinal_data) # Fit with LMM (linear mixed model)
surv_fit <- survreg(Surv(T_obs, status) ~ arm, data = survival_data, dist = "exponential") # Fit with AFT (survival only component)


### Construct data necessary for Stan

# ID levels for each participant, used for matching longitudinal and time-to-event observations
id_levels <- sort(unique(longitudinal_data$id))
J_1_long <- match(longitudinal_data$id, id_levels)
J_1_unique <- match(survival_data$id, id_levels)

# Covariate lists for longitudinal, survival and longitudinal-survival (used for inserting current value into time-to-event model) submodels
X_long <- cbind(Intercept = 1, time = longitudinal_data$time, time_arm = longitudinal_data$time*longitudinal_data$arm)
X_surv <- matrix(survival_data$arm, ncol = 1)
X_long_surv <- cbind(Intercept = 1, time = survival_data$T_obs, time_arm = survival_data$T_obs*survival_data$arm)

# Scale for longitudinal and survival outcomes
s_long <- sd(longitudinal_data$Y_obs)
s_surv <- sd(survival_data$T_obs)

# Stan data to be inputted to the Stan sAFT model
stan_data <- list(
  K_long = ncol(X_long), # number of longitudinal population level covariates
  q = ncol(X_surv), # number of survival covariates
  m = 5, # degree of bernstein polynomial
  s_long = s_long, # sample standard deviation of longitudinal outcome
  s_surv = s_surv, # sample standard deviation of survival scale
  alpha_sd = (log(2) / 1.96) * (s_long / s_surv), # standard deviation of the prior for the association parameter alpha, scaled by the ratio of the standard deviations of the longitudinal and survival outcomes
  N_long = nrow(longitudinal_data), # number of longitudinal observations
  Y_long = longitudinal_data$Y_obs_centred, # longitudinal outcome (centred)
  X_long = X_long, # longitudinal covariate matrix
  N_1_long = length(id_levels), # number of unique participants in the longitudinal data (used for indexing random effects)
  J_1_long = J_1_long, # indexing variable for matching longitudinal observations to participants (used for indexing random effects)
  Z_1_1_long = rep(1, nrow(longitudinal_data)),  # random intercept design matrix (currently only including random intercepts, but can be extended to include random slopes by adding additional columns and modifying the model code accordingly)
  Z_1_2_long = longitudinal_data$time, # random slope design matrix (currently only including random slopes for time, but can be extended to include additional random effects by adding additional columns and modifying the model code accordingly)
  n = nrow(survival_data), # number of survival observations
  status = survival_data$status, # event indicator for survival data
  time = survival_data$T_obs, # observed time for survival data
  X_surv = X_surv, # survival covariate matrix
  J_1_unique = J_1_unique, # indexing variable for matching survival observations to participants (used for indexing random effects)
  X_long_surv = X_long_surv # longitudinal covariates mapped to the survival observation times (used for including the current value of the longitudinal outcome in the survival submodel
)

# Initialisation values for MCMC fitting
init_values <- make_init(chains=4, lmm_fit, surv_fit)

# Fit Stan model
sAFT_JM_fit <- saftjm_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  seed = global_seed,
  init = init_values,
  iter_warmup = 1000,
  iter_sampling = 1000
)



### Formatting of results

# Extract linear mixed model estimate and standard error
lmm_est <- fixed.effects(lmm_fit) 
lmm_se <- sqrt(diag(vcov(lmm_fit))) 

# Extract sAFT-JM model estimate and standard error 
sAFT_JM_vars <- c("beta_long[1]", "beta_long[2]", "beta_long[3]", "gamma[1]", "alpha") # Parameters of interest
sAFT_JM_draws <- posterior::as_draws_df(sAFT_JM_fit$draws(variables = sAFT_JM_vars)) # Extract posterior draws 
sAFT_JM_est <- sapply(sAFT_JM_vars, function(x) mean(sAFT_JM_draws[[x]])) 
sAFT_JM_se <- sapply(sAFT_JM_vars, function(x) sd(sAFT_JM_draws[[x]])) 
sAFT_JM_l95 <- sapply(sAFT_JM_vars, function(x) quantile(sAFT_JM_draws[[x]], 0.025)) 
sAFT_JM_u95 <- sapply(sAFT_JM_vars, function(x) quantile(sAFT_JM_draws[[x]], 0.975)) 

# Construct dataframe of results
result_comparison <- data.frame(Parameter = c("beta_long_intercept", "beta_long_time", "beta_long_time_arm", "gamma", "alpha"), 
                                LMM = c(lmm_est["(Intercept)"] + Y_obs_mean, lmm_est["time"], lmm_est["time:arm"], NA, NA), 
                                LMM_SE = c(lmm_se["(Intercept)"], lmm_se["time"], lmm_se["time:arm"], NA, NA), 
                                LMM_L95 = c(lmm_est["(Intercept)"] - 1.96*lmm_se["(Intercept)"] + Y_obs_mean, lmm_est["time"] - 1.96*lmm_se["time"], lmm_est["time:arm"] - 1.96*lmm_se["time:arm"], NA, NA), 
                                LMM_U95 = c(lmm_est["(Intercept)"] + 1.96*lmm_se["(Intercept)"] + Y_obs_mean, lmm_est["time"] + 1.96*lmm_se["time"], lmm_est["time:arm"] + 1.96*lmm_se["time:arm"], NA, NA), 
                                sAFT_JM = c(sAFT_JM_est["beta_long[1]"] + Y_obs_mean, sAFT_JM_est["beta_long[2]"], sAFT_JM_est["beta_long[3]"], sAFT_JM_est["gamma[1]"], sAFT_JM_est["alpha"]), 
                                sAFT_JM_SE = c(sAFT_JM_se["beta_long[1]"], sAFT_JM_se["beta_long[2]"], sAFT_JM_se["beta_long[3]"], sAFT_JM_se["gamma[1]"], sAFT_JM_se["alpha"]), 
                                sAFT_JM_L95 = c(sAFT_JM_l95["beta_long[1].2.5%"] + Y_obs_mean, sAFT_JM_l95["beta_long[2].2.5%"], sAFT_JM_l95["beta_long[3].2.5%"], sAFT_JM_l95["gamma[1].2.5%"], sAFT_JM_l95["alpha.2.5%"]), 
                                sAFT_JM_U95 = c(sAFT_JM_u95["beta_long[1].97.5%"] + Y_obs_mean, sAFT_JM_u95["beta_long[2].97.5%"], sAFT_JM_u95["beta_long[3].97.5%"], sAFT_JM_u95["gamma[1].97.5%"], sAFT_JM_u95["alpha.97.5%"])) 


# Format result comparison function (SE; 95% CI)
format_result <- function(est, se, l95, u95) {
  ifelse(is.na(est),"", paste0(round(est, 4)," ( ",round(l95, 4),", ",round(u95, 4),")"))
}

# Construct clean result comparison
result_comparison_clean <- data.frame(
  Parameter = result_comparison$Parameter,
  LMM = format_result(result_comparison$LMM, result_comparison$LMM_SE, result_comparison$LMM_L95, result_comparison$LMM_U95),
  sAFT_JM = format_result(result_comparison$sAFT_JM, result_comparison$sAFT_JM_SE, result_comparison$sAFT_JM_L95, result_comparison$sAFT_JM_U95)
)
result_comparison_clean


