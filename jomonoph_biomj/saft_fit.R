
# Load packages
library(cmdstanr) # for use of Stan
library(lme4) # for LMM
library(nlme) # for LMM
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

# Fit Stan model


