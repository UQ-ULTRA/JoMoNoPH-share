
# Load packages
library(cmdstanr) # for use of Stan
library(survreg) # for AFT model
library(lme4) # for LMM

# Load data
data <- read.csv("joint_model_data.csv")


# Load Stan model
model <- cmdstan_model("saft_jm.stan")

# Fit with LMM (linear mixed model)


# Fit with AFT (survival only component)

# Fit Stan model


