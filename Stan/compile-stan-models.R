# library(tidyverse)
# library(survival)
# library(MASS)
# library(future)
# library(furrr)
# library(simsurv)
# library(rstan)
# library(brms)
library(cmdstanr)

# Set CmdStan path globally for main session
set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")


# compile and save the models
# cmdstan_model("~/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist.stan",
#               cpp_options = list(stan_threads = TRUE),
#               force_recompile = TRUE)
# 
# cmdstan_model("~/JoMoNoPH-share/Stan/Gaussian-Basis-JM-Hist.stan",
#               cpp_options = list(stan_threads = TRUE),
#               force_recompile = TRUE)

cmdstan_model("~/JoMoNoPH-share/Stan/Bernstein-Polynomials.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)

cmdstan_model("~/JoMoNoPH-share/Stan/Gaussian-Basis.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)

cmdstan_model("~/JoMoNoPH-share/Stan/Gaussian-Basis-logy.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)

cmdstan_model("~/JoMoNoPH-share/Stan/Gaussian-Basis-logy-adaptive-grid.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)

cmdstan_model("~/JoMoNoPH-share/Stan/Bernstein-Polynomials-alt.stan",
              cpp_options = list(stan_threads = TRUE),
              force_recompile = TRUE)
