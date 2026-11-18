# ------------------------------------------------------------------------------
# Case-study model fit, for audit on simulated mock data
#
# Fits the LMM and the sAFT joint model to mock-data.rds and writes the table of
# results. The mock data are simulated and are not from any real trial, so these
# results say nothing about the trial reported in the manuscript.
#
# Evaluated by run_checks.R with task = "case-study" from the package root.
# Reads this folder's mock-data.rds and the shared model under production/Stan.
# Writes under output/case-study/: mock_table1_update.csv, mock_table1_update_rows.tex,
#   mock_table1_update.rds, mock_fit_objects.rds
#
# Steps:
#   1. Read the data
#   2. Fit the LMM and an exponential AFT model (staging fits)
#   3. Build the Stan data
#   4. Generate initial values
#   5. Fit the sAFT joint model in Stan
#   6. Build the results table
#   7. Save the results
# ------------------------------------------------------------------------------

library(tidyverse)
library(survival)
library(nlme)
library(cmdstanr)

seed   <- 10235   # seeds both the initial values and the Stan sampler
chains <- 4
degree <- 5       # Bernstein polynomial degree (baseline weights theta_1 to theta_5)

# Compiled on the first run only; later runs reuse the executable
mod <- cmdstan_model(file.path("production", "Stan", "Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"),
                     cpp_options = list(stan_threads = TRUE),
                     force_recompile = FALSE)

# ==============================================================================
# 1. Read the data
# ==============================================================================
# QoL:  one row per HRQoL assessment: id, time (months), D6 score, arm
# cPFS: one row per patient: id, time (months), status (1 = event), arm

dat  <- readRDS(file.path("case-study", "source", "mock-data.rds"))
QoL  <- dat$QoL
cPFS <- dat$cPFS

# Time-by-arm interaction. There is no main effect of arm: the randomised arms
# share an intercept and differ only in slope.
QoL$time_by_arm <- QoL$time * QoL$arm

# Centre D6 to avoid a very large intercept in Stan. The centring constant is
# added back to the intercept in step 6.
longitudinal_response_center <- mean(QoL$D6, na.rm = TRUE)
QoL$D6 <- QoL$D6 - longitudinal_response_center

# Survival times must be positive: events at time 0 are moved to 1e-6 months
# (all times shift by the same negligible amount)
cPFS$time <- pmax(cPFS$time, 0) + 1e-6

# ==============================================================================
# 2. Staging fits
# ==============================================================================
# The LMM is the LMM column of the results table. Both fits also supply the
# Stan data, initial values and scale factors.

# LMM: random intercept and slope per patient, fitted by REML
fit_long <- lme(
  fixed = D6 ~ time + time_by_arm,
  random = ~ time | id,
  data = QoL,
  na.action = na.omit,
  control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
)

# Exponential AFT model for the treatment effect. Used only for an initial
# value and a scale factor; the final survival model is the joint model.
fit_surv <- survreg(Surv(time, status) ~ arm, data = cPFS,
                    dist = "exponential", x = TRUE)

# ==============================================================================
# 3. Build the Stan data
# ==============================================================================

# --- Longitudinal part: records with a non-missing D6 ---
mf_long_raw <- getData(fit_long)
rownames(mf_long_raw) <- NULL  # row names must be positions
mf_long <- model.frame(formula(fit_long), data = mf_long_raw, na.action = na.omit)
idx_long <- as.numeric(rownames(mf_long))
mf_long$id <- mf_long_raw$id[idx_long]

X_long  <- model.matrix(~ time + time_by_arm, data = mf_long)  # fixed effects
id_long <- as.factor(mf_long$id)

# --- Survival part: one row per patient ---
mf_surv <- model.frame(fit_surv)
idx_surv <- as.numeric(rownames(mf_surv))
mf_surv$id <- cPFS$id[idx_surv]
y_surv <- model.response(mf_surv)
mf_surv$time   <- y_surv[, "time"]
mf_surv$status <- y_surv[, "status"]
mf_surv$arm    <- as.numeric(cPFS$arm[idx_surv])
mf_surv$time_by_arm <- mf_surv$time * mf_surv$arm
X_surv <- model.matrix(~ 0 + arm, data = mf_surv)

# Position of each survival patient's random effects in the longitudinal model.
# Every patient needs at least one D6 value.
J_1_unique <- match(unique(mf_surv$id), levels(id_long))
if (anyNA(J_1_unique)) stop("Some patients have no D6 values.")

# --- Scale factors and priors ---
s_long <- summary(fit_long)$sigma   # LMM residual SD
s_surv <- summary(fit_surv)$scale   # exponential AFT scale
sd_b <- sqrt(as.numeric(VarCorr(fit_long)[c("(Intercept)", "time"), "Variance"]))

# Treatment effect from the exponential fit, used as an initial value only
beta_surv_init_mean <- coef(fit_surv)[-1]
beta_surv_init_mean[!is.finite(beta_surv_init_mean)] <- 0

# The raw alpha prior has SD log(2)/1.96, hence approximately 95% prior
# probability between -log(2) and log(2). Stan samples alpha_tilde, with
# alpha = alpha_tilde * s_surv /
# s_long, so its prior SD is rescaled accordingly.
alpha_prior_sd <- log(2) / 1.96
alpha_tilde_sd <- alpha_prior_sd * (s_long / s_surv)

# Names must match the data block of the .stan file
stan_data <- list(
  # Longitudinal
  N_long = nrow(X_long),
  K_long = ncol(X_long),
  Y_long = as.numeric(model.response(mf_long)),
  X_long = X_long,
  J_1_long = as.integer(id_long),              # patient index of each record
  N_1_long = nlevels(id_long),
  Z_1_1_long = rep(1, nrow(X_long)),           # random-intercept design
  Z_1_2_long = mf_long$time,                   # random-slope design

  # Survival
  n = nrow(X_surv),
  q = ncol(X_surv),
  status = as.numeric(mf_surv$status),
  time = mf_surv$time,
  X_surv = X_surv,

  # Bernstein baseline
  m = degree,
  r = 2,          # required legacy data field; no penalty is used here
  tau_h = 0.01,   # required by the data block; unused by this model

  # Linking: longitudinal design at each patient's survival time
  X_long_surv = model.matrix(~ time + time_by_arm, data = mf_surv),
  J_1_unique = J_1_unique,

  # Priors: Normal(0, 10) for the regression coefficients
  beta_long_prior_mean = rep(0, ncol(X_long)),
  beta_long_prior_scale = rep(10, ncol(X_long)),
  beta_surv_prior_mean = rep(0, length(beta_surv_init_mean)),
  beta_surv_prior_scale = rep(10, length(beta_surv_init_mean)),
  beta_surv_init_mean = beta_surv_init_mean,
  sd_b_prior_mean = sd_b,
  sd_b_prior_scale = rep(5, length(sd_b)),
  s_long = s_long,
  s_surv = s_surv,
  alpha_tilde_sd = alpha_tilde_sd
)

# ==============================================================================
# 4. Initial values
# ==============================================================================
# One set per chain: the staging estimates plus random jitter, so the chains
# start from different points

set.seed(seed)

beta_long_start <- as.numeric(fixef(fit_long))
L_start <- t(chol(matrix(c(1, 0.2, 0.2, 1), 2, 2)))  # random-effects correlation 0.2

init_values <- lapply(seq_len(chains), function(chain) {
  sd_start <- pmax(sd_b + abs(rnorm(2, 0, 0.05)), 1e-3)
  list(
    z_1_long    = matrix(rnorm(2 * stan_data$N_1_long), nrow = 2),
    L_1_long    = L_start,
    beta_long   = beta_long_start + rnorm(stan_data$K_long, 0, 0.1),
    sigma_long  = abs(rnorm(1, s_long, 1)),
    sd_1_long   = sd_start,
    beta_surv   = rnorm(stan_data$q, as.numeric(beta_surv_init_mean), 0.05),
    gamma       = abs(rnorm(stan_data$m, 0.1, 0.01)),  # baseline weights
    alpha_tilde = rnorm(1, 0, 0.1)
  )
})

# ==============================================================================
# 5. Fit the sAFT joint model
# ==============================================================================

fit <- mod$sample(
  data = stan_data,
  seed = seed,
  init = init_values,
  chains = chains,
  parallel_chains = chains,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  max_treedepth = 10)

summ <- fit$summary()             # includes R-hat and ESS for every parameter
diag <- fit$diagnostic_summary()  # divergences and treedepth hits by chain

# ==============================================================================
# 6. Build the results table
# ==============================================================================

# "estimate (lower, upper)" to 3 decimals
format_est_ci <- function(mean, lower, upper) {
  sprintf("%.3f (%.3f, %.3f)", mean, lower, upper)
}

# sAFT-JM column: posterior mean and 95% equal-tailed credible interval
saft <- fit$draws(c("beta_long", "sd_1_long", "sigma_long",
                    "beta_surv", "alpha", "gamma")) %>%
  posterior::summarise_draws("mean", ~quantile(.x, probs = c(0.025, 0.975))) %>%
  transmute(variable, mean, lower = `2.5%`, upper = `97.5%`)

# LMM column: Wald intervals for the fixed effects; intervals() for the SDs
fixed_est <- fixef(fit_long)
fixed_se  <- sqrt(diag(vcov(fit_long)))
lmm_int   <- intervals(fit_long)
re_int    <- lmm_int$reStruct$id[c("sd((Intercept))", "sd(time)"), ]

lmm <- tibble(
  variable = c("beta_long[1]", "beta_long[2]", "beta_long[3]",
               "sd_1_long[1]", "sd_1_long[2]", "sigma_long"),
  mean  = c(fixed_est, sd_b, s_long),
  lower = c(fixed_est - qnorm(0.975) * fixed_se, re_int[, "lower"],
            lmm_int$sigma["lower"]),
  upper = c(fixed_est + qnorm(0.975) * fixed_se, re_int[, "upper"],
            lmm_int$sigma["upper"])
)

# Put the intercept back on the original D6 scale
uncentre <- function(tab) {
  row <- tab$variable == "beta_long[1]"
  tab[row, c("mean", "lower", "upper")] <-
    tab[row, c("mean", "lower", "upper")] + longitudinal_response_center
  tab
}

# Rows in table order, with manuscript labels. beta_surv[1] is gamma in Table 1;
# gamma[1] to gamma[5] are the Bernstein baseline weights theta_1 to theta_5.
labels <- c(
  "beta_long[1]" = "$\\beta_0$",
  "beta_long[2]" = "$\\beta_1$",
  "beta_long[3]" = "$\\beta_2$",
  "sd_1_long[1]" = "$\\sigma_{b0}$",
  "sd_1_long[2]" = "$\\sigma_{b1}$",
  "sigma_long"   = "$\\sigma_\\epsilon$",
  "beta_surv[1]" = "$\\gamma$",
  "alpha"        = "$\\alpha$",
  "gamma[1]"     = "$\\theta_1$",
  "gamma[2]"     = "$\\theta_2$",
  "gamma[3]"     = "$\\theta_3$",
  "gamma[4]"     = "$\\theta_4$",
  "gamma[5]"     = "$\\theta_5$"
)

table1_update <- tibble(variable = names(labels)) %>%
  left_join(uncentre(saft) %>%
              transmute(variable, sAFT_JM = format_est_ci(mean, lower, upper)),
            by = "variable") %>%
  left_join(uncentre(lmm) %>%
              transmute(variable, LMM = format_est_ci(mean, lower, upper)),
            by = "variable") %>%
  mutate(parameter = unname(labels[variable]),
         LMM = replace_na(LMM, "")) %>%   # the LMM has no survival parameters
  dplyr::select(parameter, variable, sAFT_JM, LMM)

print(table1_update)

# ==============================================================================
# 7. Save the results
# ==============================================================================

prior_info <- list(longitudinal_response_center = longitudinal_response_center,
                   s_long = s_long, s_surv = s_surv,
                   alpha_prior_sd = alpha_prior_sd, alpha_tilde_sd = alpha_tilde_sd)

output_dir <- file.path("output", "case-study")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(table1_update, file.path(output_dir, "mock_table1_update.csv"), row.names = FALSE)
writeLines(sprintf("%s & %s & %s \\\\", table1_update$parameter,
                   table1_update$sAFT_JM, table1_update$LMM),
           file.path(output_dir, "mock_table1_update_rows.tex"))
saveRDS(list(table1_update = table1_update, prior_info = prior_info),
        file.path(output_dir, "mock_table1_update.rds"))

# The audit runner omits the R6 fit pointer when saving this bundle; numeric
# summaries and diagnostics remain portable, and full draws are in Stan CSVs.
saveRDS(list(fit = fit, summ = summ, diag = diag,
             table1_update = table1_update, prior_info = prior_info),
        file.path(output_dir, "mock_fit_objects.rds"))
