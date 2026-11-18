# ------------------------------------------------------------------------------
# Script Title: Bayesian Joint Modelling 
# Description :
#   This script performs simulation and Bayesian inference for a joint 
#   longitudinal-survival model. The steps include:
#     - Compiling a Stan model (e.g., `Bernstein-Polynomials-JM.stan`)
#     - Setting model and prior parameters
#     - Simulating synthetic joint data (longitudinal + survival)
#     - Converting data to Stan-compatible format
#     - Fitting the model using `cmdstanr`
#     - Extracting and summarizing posterior draws
#     - Relabelling and presenting results for interpretation
# ------------------------------------------------------------------------------
rm(list = ls())

# Read the list back in
out_obj <- readRDS("ENZAMET-data.rds")

# Unpack into the global environment
list2env(out_obj, envir = .GlobalEnv)




library(tidyverse)
library(survival)
library(MASS)
library(rstan)
library(brms)
library(cmdstanr)
library(here)
library(JMbayes2)

# Compile the Stan model
# mod <- cmdstan_model(here("Gaussian-Basis-JM-Hist.stan"),
#                      cpp_options = list(stan_threads = TRUE),
#                      force_recompile = FALSE)


mod <- cmdstan_model(here("Stan", "Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"),
                     cpp_options = list(stan_threads = TRUE),
                     force_recompile = TRUE)




# The current alt-DM Stan model computes the Bernstein basis internally, so the
# legacy BP helper file is not needed for preparing stan_data or init values.
# source(here("Bunya", "Auxiliary", "BP-functions.R"))


QoL <- QoL %>% 
  rename(id = SUBJECTNUMBERSTR_0s,
         # time = study_months,
         arm = arm_num) %>% 
  mutate(time_by_arm = time * arm)

# Match the current ENZAMET simulation strategy: center the longitudinal
# response before fitting to avoid a very large Stan intercept. Keep the raw
# value for later back-transformation of beta_long[1] if required.
longitudinal_response_center <- mean(QoL$D6, na.rm = TRUE)
if (!is.finite(longitudinal_response_center)) {
  stop("Observed D6 mean is not finite; cannot center the longitudinal response.")
}
QoL <- QoL %>%
  mutate(
    D6_raw = D6,
    D6 = D6 - longitudinal_response_center
  )


fit_long <-   model <- lme(
  fixed = D6~ time + time_by_arm,
  #fixed = D6 ~ ns(study_months, 3) * arm_num,
  random = ~ time | id,
  data = QoL,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(msverbose = FALSE, msMaxIter = 500, msMaxEval = 500, opt = "optim")
)

summary(fit_long)  



cPFS <- cPFS %>% 
  rename(id = SUBJECTNUMBERSTR_0s,
         # time = study_months,
         arm = arm_num)

# Optional Cox-JM comparator fit retained for reference. This runs before the
# proposed sAFT(BP)-JM fit and uses fit_surv_cox so it does not overwrite the
# later exponential survreg() staging fit used for Stan data/init values.
if (TRUE) {
  fit_surv_cox <- coxph(
    Surv(time, status) ~ arm,
    data = cPFS,
    ties = "breslow",
    x = TRUE
  )

  model_jm <- jm(
    fit_surv_cox,
    fit_long,
    time_var = "time",
    control = list(n_iter = 5000, n_burnin = 1000, n_chains = 4, cores = 4)
  )

  print(model_jm)
}



# Set model parameters
chains <- 4
threads_per_chain <- 1
parallel_chains <- chains
degree <- 5

alpha_prior_calibration <- function() {
  # Clinically interpretable prior for the unstandardised AFT association
  # alpha: by default, a 1.96-unit longitudinal-response change is expected
  # to change survival time by at most a factor of 2 on the 95% prior scale.
  # Override with SIM_ALPHA_DELTA_Y / SIM_ALPHA_MAX_F only for sensitivity runs.
  delta_y <- as.numeric(Sys.getenv("SIM_ALPHA_DELTA_Y", unset = "1.96"))
  max_f <- as.numeric(Sys.getenv("SIM_ALPHA_MAX_F", unset = "2"))
  if (!is.finite(delta_y) || delta_y <= 0) {
    stop("SIM_ALPHA_DELTA_Y must be a positive finite number.")
  }
  if (!is.finite(max_f) || max_f < 1) {
    stop("SIM_ALPHA_MAX_F must be finite and >= 1.")
  }
  list(
    delta_y = delta_y,
    max_f = max_f,
    sigma_alpha = log(max_f) / delta_y
  )
}


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
  if (anyNA(J_1_unique)) {
    missing_ids <- unique(mf_surv$id)[is.na(J_1_unique)]
    stop(
      "Some survival IDs are absent from the fitted longitudinal model: ",
      paste(utils::head(missing_ids, 10), collapse = ", "),
      if (length(missing_ids) > 10) " ..." else "",
      ". The current Stan joint model requires every survival subject to have ",
      "a longitudinal-model random-effect level."
    )
  }
  
  # The response has already been centered before fit_long is fitted. Do not
  # center survival times or covariates here.
  mf_surv$time_by_arm <- mf_surv$time * mf_surv$arm
  
  X_long_surv <- model.matrix(~ time + time_by_arm, data = mf_surv)
  
  
  # --- Empirical prior components ---
  
  mf_long$X_long <- X_long  
  
  # Coefficient priors follow the current simulation settings: zero-centered,
  # weakly informative Normal(0, 10). The lme() and exponential survreg()
  # estimates are used for initial values only, not as prior means.
  # Previous empirical prior centering:
  # beta_long_mean <- nlme::fixef(fit_long)
  beta_long_mean <- rep(0, ncol(X_long))
  beta_long_sd <- rep(10, length(beta_long_mean))

  # Previous empirical prior centering:
  # beta_surv_mean <- coef(fit_surv)[-1]  # omit intercept
  beta_surv_init_mean <- coef(fit_surv)[-1]  # omit intercept
  beta_surv_init_mean[!is.finite(beta_surv_init_mean)] <- 0
  beta_surv_mean <- rep(0, length(beta_surv_init_mean))
  beta_surv_sd <- rep(10, length(beta_surv_mean))
  
  vc_long <- tryCatch(VarCorr(fit_long), error = function(e) NULL)
  sd_b <- if (!is.null(vc_long)) {
    sqrt(as.numeric(vc_long[c("(Intercept)", "time"), "Variance"]))
  } else rep(1.0, 2)

  s_long <- summary(fit_long)$sigma
  s_surv <- summary(fit_surv)$scale

  # Stan samples alpha_tilde and transforms internally via
  # alpha = alpha_tilde * (s_surv / s_long). Therefore alpha_tilde_sd is scaled
  # so that the implied prior SD for raw alpha equals log(max_f) / delta_y.
  alpha_cal <- alpha_prior_calibration()
  alpha_tilde_sd <- alpha_cal$sigma_alpha * (s_long / s_surv)
  
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
    tau_h = 0.01,   # Kept for Stan data compatibility; unused by this alt-DM model
    
    # Linking
    X_long_surv = X_long_surv,
    J_1_unique = J_1_unique,
    
    # Priors
    beta_long_prior_mean = beta_long_mean,
    beta_long_prior_scale = beta_long_sd,
    beta_surv_prior_mean = beta_surv_mean,
    beta_surv_prior_scale = beta_surv_sd,
    beta_surv_init_mean = beta_surv_init_mean,
    sd_b_prior_mean = sd_b,
    sd_b_prior_scale = rep(5, length(sd_b)),
    s_long = s_long,
    s_surv = s_surv,
    alpha_tilde_sd = alpha_tilde_sd
  )
}
# ------------------------------------------------------------------------------ 
# Fit longitudinal model on raw variables 
# ------------------------------------------------------------------------------ 

fit_long <- tryCatch({
  nlme::lme(
    fixed = D6 ~ time + time_by_arm,
    random = ~ time | id,
    data = QoL,
    na.action = na.omit,
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
  )
}, error = function(e) {
  message(sprintf("lme() failed: %s", e$message))
  return(NULL)
})


# ------------------------------------------------------------------------------ 
# Fit survival model on raw variables 
# ------------------------------------------------------------------------------ 

cPFS$time <- pmax(cPFS$time, 0) + 1e-6

fit_surv <- survival::survreg(
  Surv(time, status) ~ arm,
  data = cPFS,
  # This staging model is only used to build Stan data and initial values.
  # The final survival model is the flexible BP AFT-JM in Stan.
  dist = "exponential",
  x = TRUE
)

coxph(Surv(time, status) ~ arm, 
      data =cPFS, 
      ties="breslow",
      x = TRUE)


# Convert all inputs to Stan data format
stan_data <- convert_data_from_models(fit_long, 
                                      fit_surv, 
                                      surv_data = cPFS,
                                      degree = 5)



safe_surv_beta_init <- function(stan_data, fallback_value = 0) {
  # beta_surv initial values come from beta_surv_init_mean, which is built from
  # survreg(..., dist = "exponential"). The prior remains centered at zero.
  beta <- tryCatch(as.numeric(stan_data$beta_surv_init_mean), error = function(e) numeric(0))
  if (length(beta) != stan_data$q) beta <- rep(fallback_value, stan_data$q)
  beta[!is.finite(beta)] <- fallback_value
  beta
}

# Initialisation function
init_fun <- function(stan_data, fit_long = NULL, chains = 4) {
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
  
  # Use the already-fitted lme() object for empirical starting values. These
  # are initial values only; coefficient priors are centered at zero in stan_data.
  if (inherits(fit_long, "lme")) {
    # Fixed effects (intercept and covariate effects)
    fallback$beta_long <- tryCatch(as.numeric(nlme::fixef(fit_long)), error = function(e) fallback$beta_long)
    
    # Residual error SD
    fallback$sigma_long <- tryCatch(summary(fit_long)$sigma, error = function(e) fallback$sigma_long)
    
    # Random intercept/slope SDs
    vc <- tryCatch(nlme::VarCorr(fit_long), error = function(e) NULL)
    if (!is.null(vc)) {
      fallback$sd_1_long <- tryCatch(
        sqrt(as.numeric(vc[1:2, "Variance"])),  # First two rows are intercept and slope variances
        error = function(e) fallback$sd_1_long
      )
    }
  }
  
  # Exponential survreg() staging estimate, used only for beta_surv init.
  fallback$beta_surv <- safe_surv_beta_init(stan_data, fallback_value = 0)
  
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
      alpha_tilde = rnorm(1, fallback$alpha_tilde, 0.1)
    )
  })
  
  return(init_list)
}

# ------

init_values <- init_fun(stan_data = stan_data, fit_long = fit_long, chains = chains)

i <- 235
# --- Fit model ---
fit <- mod$sample(
  data = stan_data,
  seed = 10000 + i,
  init = init_values,
  chains = chains,
  parallel_chains = parallel_chains,
  threads_per_chain = threads_per_chain,
  iter_warmup = 1000,
  iter_sampling = 1000,
  max_treedepth = 10)

summ <- fit$summary()

sum <- posterior::summarise_draws(fit,
                "mean", "median",
                ~quantile(.x, probs = c(0.025, 0.975)))



diag <- fit$diagnostic_summary()

summary<-summ %>% 
  filter(grepl("alpha|beta|gamma|sigma", variable))

format_est_ci <- function(mean, lower, upper, digits = 3) {
  ifelse(
    is.na(mean) | is.na(lower) | is.na(upper),
    "",
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
      mean, lower, upper
    )
  )
}

summarise_saft_for_table <- function(fit, longitudinal_response_center) {
  # These posterior summaries are table-ready. beta_long[1] is shifted back
  # from the centered D6 scale to the original D6 scale; all other parameters
  # are unchanged by response centering.
  draws <- fit$draws(
    variables = c("beta_long", "sd_1_long", "sigma_long", "beta_surv", "alpha", "gamma")
  )
  tab <- posterior::summarise_draws(
    draws,
    "mean",
    ~quantile(.x, probs = c(0.025, 0.975))
  ) %>%
    as.data.frame() %>%
    transmute(
      variable,
      mean = mean,
      lower = `2.5%`,
      upper = `97.5%`
    )

  beta0_row <- tab$variable == "beta_long[1]"
  tab[beta0_row, c("mean", "lower", "upper")] <-
    tab[beta0_row, c("mean", "lower", "upper")] + longitudinal_response_center
  tab
}

summarise_lmm_for_table <- function(fit_long, longitudinal_response_center) {
  # LMM estimates are also on the centered D6 scale because fit_long was fitted
  # after centering. Shift only the intercept back to the original D6 scale.
  fixed_est <- nlme::fixef(fit_long)
  fixed_vcov <- stats::vcov(fit_long)
  fixed_se <- sqrt(diag(fixed_vcov))
  fixed_names <- c("(Intercept)", "time", "time_by_arm")
  fixed_tab <- data.frame(
    variable = c("beta_long[1]", "beta_long[2]", "beta_long[3]"),
    mean = as.numeric(fixed_est[fixed_names]),
    lower = as.numeric(fixed_est[fixed_names] - stats::qnorm(0.975) * fixed_se[fixed_names]),
    upper = as.numeric(fixed_est[fixed_names] + stats::qnorm(0.975) * fixed_se[fixed_names]),
    stringsAsFactors = FALSE
  )
  fixed_tab[fixed_tab$variable == "beta_long[1]", c("mean", "lower", "upper")] <-
    fixed_tab[fixed_tab$variable == "beta_long[1]", c("mean", "lower", "upper")] +
    longitudinal_response_center

  vc <- nlme::VarCorr(fit_long)
  sd_est <- sqrt(as.numeric(vc[c("(Intercept)", "time"), "Variance"]))
  sigma_est <- summary(fit_long)$sigma

  int <- tryCatch(nlme::intervals(fit_long), error = function(e) NULL)
  random_lower <- random_upper <- rep(NA_real_, 2)
  sigma_lower <- sigma_upper <- NA_real_
  if (!is.null(int)) {
    re_mat <- NULL
    if (!is.null(int$reStruct)) {
      re_mat <- if (is.list(int$reStruct)) int$reStruct[[1]] else int$reStruct
    }
    if (!is.null(re_mat)) {
      re_rows <- grep("sd", rownames(re_mat), ignore.case = TRUE)
      if (length(re_rows) < 2) {
        re_rows <- setdiff(seq_len(nrow(re_mat)), grep("cor", rownames(re_mat), ignore.case = TRUE))
      }
      if (length(re_rows) >= 2) {
        random_lower <- as.numeric(re_mat[re_rows[1:2], "lower"])
        random_upper <- as.numeric(re_mat[re_rows[1:2], "upper"])
      }
    }
    if (!is.null(int$sigma)) {
      sigma_lower <- as.numeric(int$sigma["lower"])
      sigma_upper <- as.numeric(int$sigma["upper"])
    }
  }

  random_tab <- data.frame(
    variable = c("sd_1_long[1]", "sd_1_long[2]", "sigma_long"),
    mean = c(sd_est, sigma_est),
    lower = c(random_lower, sigma_lower),
    upper = c(random_upper, sigma_upper),
    stringsAsFactors = FALSE
  )

  bind_rows(fixed_tab, random_tab)
}

param_order <- c(
  "beta_long[1]", "beta_long[2]", "beta_long[3]",
  "sd_1_long[1]", "sd_1_long[2]", "sigma_long",
  "beta_surv[1]", "alpha",
  "gamma[1]", "gamma[2]", "gamma[3]", "gamma[4]", "gamma[5]"
)
# Manuscript notation differs from the Stan object names:
# - beta_surv[1] is the AFT treatment coefficient shown as gamma in Table 1.
# - gamma[1:5] are the Bernstein-polynomial baseline hazard weights, shown as
#   theta_1,...,theta_5 in Table 1.
param_labels <- c(
  "beta_long[1]" = "$\\beta_0$",
  "beta_long[2]" = "$\\beta_1$",
  "beta_long[3]" = "$\\beta_2$",
  "sd_1_long[1]" = "$\\sigma_{b0}$",
  "sd_1_long[2]" = "$\\sigma_{b1}$",
  "sigma_long" = "$\\sigma_\\epsilon$",
  "beta_surv[1]" = "$\\gamma$",
  "alpha" = "$\\alpha$",
  "gamma[1]" = "$\\theta_1$",
  "gamma[2]" = "$\\theta_2$",
  "gamma[3]" = "$\\theta_3$",
  "gamma[4]" = "$\\theta_4$",
  "gamma[5]" = "$\\theta_5$"
)

saft_table <- summarise_saft_for_table(fit, longitudinal_response_center) %>%
  transmute(variable, sAFT_JM = format_est_ci(mean, lower, upper))
lmm_table <- summarise_lmm_for_table(fit_long, longitudinal_response_center) %>%
  transmute(variable, LMM = format_est_ci(mean, lower, upper))

table1_update <- tibble(variable = param_order) %>%
  left_join(saft_table, by = "variable") %>%
  left_join(lmm_table, by = "variable") %>%
  mutate(
    parameter = unname(param_labels[variable]),
    sAFT_JM = tidyr::replace_na(sAFT_JM, ""),
    LMM = tidyr::replace_na(LMM, "")
  ) %>%
  select(parameter, variable, sAFT_JM, LMM)

table1_latex_rows <- table1_update %>%
  transmute(row = sprintf("%s & %s & %s \\\\", parameter, sAFT_JM, LMM)) %>%
  pull(row)

alpha_prior_info <- list(
  longitudinal_response_center = longitudinal_response_center,
  s_long = stan_data$s_long,
  s_surv = stan_data$s_surv,
  alpha_tilde_sd = stan_data$alpha_tilde_sd,
  alpha_prior_delta_y = alpha_prior_calibration()$delta_y,
  alpha_prior_max_f = alpha_prior_calibration()$max_f,
  alpha_prior_sd_raw = stan_data$alpha_tilde_sd * stan_data$s_surv / stan_data$s_long
)

write.csv(table1_update, "ENZAMET_table1_update.csv", row.names = FALSE)
writeLines(table1_latex_rows, "ENZAMET_table1_update_rows.tex")
saveRDS(
  list(table1_update = table1_update, alpha_prior_info = alpha_prior_info),
  file = "ENZAMET_table1_update.rds"
)
message(
  sprintf(
    "Alpha prior raw SD = %.8f; expected default log(2)/1.96 = %.8f",
    alpha_prior_info$alpha_prior_sd_raw,
    log(2) / 1.96
  )
)

# Put everything you care about into one list object
out_obj <- list(
  fit           = fit,
  summ          = summ,
  diag          = diag,
  summary       = summary,
  table1_update = table1_update,
  alpha_prior_info = alpha_prior_info,
  longitudinal_response_center = longitudinal_response_center
)

# Write to RDS
saveRDS(out_obj, file = "BP_objects.rds")

# ------------------------------------------------------------------------------
# Legacy comparison/table code below is retained for reference, but is not needed
# for the single proposed sAFT(BP)-JM ENZAMET fit. Leaving it active can make an
# otherwise successful Stan run fail later if GB_objects.rds or brms inputs are
# not present in the working directory.
# ------------------------------------------------------------------------------
if (FALSE) {

BP_objects <- readRDS("BP_objects.rds")
GB_objects <- readRDS("GB_objects.rds")


library(dplyr)
library(stringr)



library(brms)

# same structure as your lme:
# fixed = D6 ~ study_months + study_months:arm_num
# random = ~ study_months | SUBJECTNUMBERSTR_0s

bayes_long <- brm(
  formula = D6 ~ study_months * arm_num + (study_months | SUBJECTNUMBERSTR_0s),
  data    = out_obj$QoL,
  family  = gaussian(),
  prior   = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 2),  class = "b"),          # regression coefficients
    prior(cauchy(0, 2),  class = "sd"),         # random effect SDs
    prior(cauchy(0, 2),  class = "sigma")       # residual SD
  ),
  cores   = 4,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)


smry_long <- posterior_summary(
  bayes_long,
  variable = c(
    "b_Intercept",
    "b_study_months",
    "b_arm_num",
    "b_study_months:arm_num",
    "sd_SUBJECTNUMBERSTR_0s__Intercept",
    "sd_SUBJECTNUMBERSTR_0s__study_months",
    "sigma"
  ),
  probs = c(0.025, 0.975)
)

tab_long <- smry_long %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  transmute(
    variable = dplyr::recode(
      param,
      # fixed effects
      "b_Intercept"               = "beta_long[1]",
      "b_study_months"            = "beta_long[2]",
      "b_arm_num"                 = "beta_long[3]",          # or map to [4] if you prefer
      "b_study_months:arm_num"    = "beta_long[4]",          # optional: add in param_map / order
      # random effects SDs
      "sd_SUBJECTNUMBERSTR_0s__Intercept"     = "sd_1_long[1]",
      "sd_SUBJECTNUMBERSTR_0s__study_months"  = "sd_1_long[2]",
      # residual SD
      "sigma"                     = "sigma_long"
    ),
    mean_long = Estimate,
    q25_long   = Q2.5,
    q975_long  = Q97.5
  )



## 1. Extract just what we need from each model -----------------------------

tab_BP <- posterior::summarise_draws(BP_objects$fit,
                             "mean", "median",
                             ~quantile(.x, probs = c(0.025, 0.975))) %>% 
  dplyr::select(variable,
         mean_BP = mean,
         q25_BP   = `2.5%`,
         q975_BP  = `97.5%`)

tab_GB <- posterior::summarise_draws(GB_objects$fit,
                                     "mean", "median",
                                     ~quantile(.x, probs = c(0.025, 0.975))) %>% 
  dplyr::select(variable,
                mean_GB = mean,
                q25_GB   = `2.5%`,
                q975_GB  = `97.5%`)




## Join BP + GB by parameter name
tab <- full_join(tab_BP, tab_GB, by = "variable") %>%
  full_join(tab_long, by = "variable")

## 2. Map variable names to LaTeX symbols ----------------------------------

param_map <- c(
  "beta_long[1]" = "$\\beta_{0}$",
  "beta_long[2]" = "$\\beta_{1}$",
  "beta_long[3]" = "$\\beta_{2}$",
  # add these if they exist in your summaries:
  "sd_1_long[1]"     = "$\\sigma_{b0}$",
  "sd_1_long[2]"     = "$\\sigma_{b1}$",
  "sigma_long"   = "$\\sigma_{\\varepsilon}$",
  "beta_surv[1]"     = "$\\gamma_{1}$",
  "alpha"  = "$\\alpha$",
  "gamma[1]" = "$\\beta_{u_1}$",
  "gamma[2]" = "$\\beta_{u_2}$",
  "gamma[3]" = "$\\beta_{u_3}$",
  "gamma[4]" = "$\\beta_{u_4}$",
  "gamma[5]" = "$\\beta_{u_5}$"
)

## Order of rows in the table
param_order <- c(
  "beta_long[1]",
  "beta_long[2]",
  "beta_long[3]",
  "sd_1_long[1]",
  "sd_1_long[2]",
  "sigma_long",
  "beta_surv[1]",
  "alpha",
  "gamma[1]", "gamma[2]", "gamma[3]", "gamma[4]", "gamma[5]"
)

## 3. Helper to format mean (q5, q95) --------------------------------------

fmt <- function(m, l, u) {
  ifelse(
    is.na(m),
    "",
    sprintf("%.3f (%.3f, %.3f)", m, l, u)
  )
}

## 4. Build the table rows --------------------------------------------------

tab_latex <- tab %>%
  filter(variable %in% param_order) %>%
  mutate(
    latex_name = param_map[variable],
    BP   = fmt(mean_BP, q25_BP, q975_BP),
    GB   = fmt(mean_GB, q25_GB, q975_GB),
    Long = fmt(mean_long, q25_long, q975_long)   # new column
  ) %>%
  arrange(factor(variable, levels = param_order))

latex_rows <- sapply(seq_len(nrow(tab_latex)), function(i) {
  r <- tab_latex[i, ]
  end <- if (i == nrow(tab_latex)) "\\\\" else "\\\\[4pt]"
  sprintf("%s  & %s  & %s  & %s  %s",
          r$latex_name, r$BP, r$GB, r$Long, end)
})

cat(paste(latex_rows, collapse = "\n"))
}
