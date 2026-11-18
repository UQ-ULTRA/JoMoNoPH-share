# sAFT-JM simulation functions: local BP joint-model path only.
# run_simulation(): one original replicate i/config_hash, with explicit baseline
# and censoring inputs; returns coefficient summaries, diagnostics and runtimes.
# The seed is 1e6 * config_hash + i. No MCMC controls are shortened here.
# convert_data_from_models(): lme/survreg fits + survival rows + basis count ->
# Stan data, prior parameters and survival initial-value anchors.
# compute_waic_from_log_lik(): fitted Stan object and variable name -> WAIC.
# alpha_prior_calibration(): scoped prior options -> contrast and alpha SD.

compute_waic_from_log_lik <- function(fit, log_lik_var = "log_lik") {
  log_lik_draws <- fit$draws(variables = log_lik_var, format = "matrix")
  if (is.null(log_lik_draws) || ncol(log_lik_draws) == 0) {
    stop("No draws found for Stan log-likelihood variable: ", log_lik_var)
  }
  loo::waic(log_lik_draws)
}

alpha_prior_calibration <- function() {
  delta_y <- as.numeric(Sys.getenv("SIM_ALPHA_DELTA_Y", unset = "1.96"))
  max_f <- as.numeric(Sys.getenv("SIM_ALPHA_MAX_F", unset = "2"))
  if (!is.finite(delta_y) || delta_y <= 0) {
    stop("SIM_ALPHA_DELTA_Y must be a positive finite number.")
  }
  if (!is.finite(max_f) || max_f < 1) {
    stop("SIM_ALPHA_MAX_F must be finite and >= 1.")
  }
  list(delta_y = delta_y, max_f = max_f, sigma_alpha = log(max_f) / delta_y)
}

chains <- 4

threads_per_chain <- 1

parallel_chains <- chains

k_bases <- 5

joint_model_type <- "BP"

convert_data_from_models <- function(fit_long, fit_surv, surv_data, k_bases = 5) {
  if (!inherits(fit_long, "lme")) {
    stop("'fit_long' must be an object of class 'lme'")
  }
  if (!inherits(fit_surv, "survreg")) {
    stop("'fit_surv' must be an object of class 'survreg'")
  }
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
  mf_long$X_long <- X_long
  beta_long_mean <- rep(0, ncol(X_long))
  beta_long_sd <- rep(10, length(beta_long_mean))
  beta_surv_init_mean <- coef(fit_surv)[-1]
  beta_surv_init_mean[!is.finite(beta_surv_init_mean)] <- 0
  beta_surv_mean <- rep(0, length(beta_surv_init_mean))
  beta_surv_sd <- rep(10, length(beta_surv_mean))
  vc_long <- VarCorr(fit_long)
  sd_b <- sqrt(as.numeric(vc_long[c("(Intercept)", "time"), "Variance"]))
  s_long <- summary(fit_long)$sigma
  s_surv <- summary(fit_surv)$scale
  alpha_cal <- alpha_prior_calibration()
  sigma_alpha <- alpha_cal$sigma_alpha
  alpha_tilde_sd <- sigma_alpha * (s_long / s_surv)
  list(
    N_long = nrow(X_long), K_long = ncol(X_long), Y_long = Y_long, X_long = X_long, J_1_long = J_1_long,
    N_1_long = N_1_long, Z_1_1_long = Z_1_1_long, Z_1_2_long = Z_1_2_long, n = nrow(X_surv),
    q = ncol(X_surv), status = as.numeric(mf_surv$status), time = mf_surv$time, X_surv = X_surv,
    m = k_bases, X_long_surv = X_long_surv, J_1_unique = J_1_unique, beta_long_prior_mean = beta_long_mean,
    beta_long_prior_scale = beta_long_sd, beta_surv_prior_mean = beta_surv_mean, beta_surv_prior_scale = beta_surv_sd,
    beta_surv_init_mean = beta_surv_init_mean, sd_b_prior_mean = sd_b, sd_b_prior_scale = rep(
      5,
      length(sd_b)
    ), s_long = s_long, s_surv = s_surv, alpha_tilde_sd = alpha_tilde_sd,
    location_beta_real = 0, scale_beta_real = 1, std_real = 1, means_real = 0, r = 2, tau_h = 0.01
  )
}

run_simulation <- function(
  i, n_patients, lambda_c, aft_mode, max_FU, config_hash, joint_model_type = "BP",
  k_bases = 5, scenario = 2, baseline_label = aft_mode, cfg_weibull_shape, cfg_weibull_scale,
  cfg_loglogistic_shape, cfg_loglogistic_scale
) {
  force(i)
  force(n_patients)
  force(lambda_c)
  force(aft_mode)
  force(max_FU)
  force(config_hash)
  force(joint_model_type)
  if (!identical(joint_model_type, "BP")) {
    stop("This local module implements the BP joint model only.")
  }
  force(k_bases)
  force(scenario)
  force(baseline_label)
  force(cfg_weibull_shape)
  force(cfg_weibull_scale)
  force(cfg_loglogistic_shape)
  force(cfg_loglogistic_scale)
  list2env(simulation_parameters(scenario), envir = environment())
  # Data generation uses the original configuration/replicate seed rule.
  unique_seed <- 1e+06 * config_hash + i
  print(unique_seed)
  sim_dat <- simulate_joint_dataset(
    D = D, beta_0 = beta_0, beta_1 = beta_1, beta_2 = beta_2,
    sigma_e = sigma_e, log_HR = log_HR, log_AF = log_AF, alpha_PH = 99, alpha_AFT = alpha_AFT,
    weibull_shape = cfg_weibull_shape, weibull_scale = cfg_weibull_scale, loglogistic_shape = cfg_loglogistic_shape,
    loglogistic_scale = cfg_loglogistic_scale, visit = visit, seed = unique_seed, n_patients = n_patients,
    max_FU = max_FU, lambda_c = lambda_c, aft_mode = aft_mode, link_type = "value"
  )
  sim_dat$longitudinal <- sim_dat$longitudinal %>% mutate(time_by_arm = time * arm)
  longitudinal_response_center <- 0
  longitudinal_response_center <- mean(sim_dat$longitudinal$Y_obs, na.rm = TRUE)
  if (!is.finite(longitudinal_response_center)) {
    stop("Observed longitudinal response mean is not finite; cannot center Y_obs.")
  }
  sim_dat$longitudinal$Y_obs <- sim_dat$longitudinal$Y_obs - longitudinal_response_center
  message(sprintf(
    "Centered longitudinal response by subtracting observed mean %.6f before fitting.",
    longitudinal_response_center
  ))
  t_lme_start <- Sys.time()
  fit_long <- NULL
  fit_long <- nlme::lme(
    fixed = Y_obs ~ time + time_by_arm, random = ~ time | id, data = sim_dat$longitudinal,
    na.action = na.omit, control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
  )
  t_lme_elapsed <- as.numeric(difftime(Sys.time(), t_lme_start, units = "secs"))
  sim_dat$survival <- sim_dat$survival %>% dplyr::rename(time = T_obs)
  # Preliminary exponential AFT fit supplies design/scaling and initial values.
  fit_surv <- survival::survreg(Surv(time, status) ~ arm,
    data = sim_dat$survival, dist = "exponential",
    x = TRUE
  )
  stan_model_file_joint <- stan_model_path
  stan_model_file_override <- Sys.getenv("SIM_STAN_MODEL_FILE", unset = "")
  if (nzchar(stan_model_file_override)) {
    stan_model_file_joint <- stan_model_file_override
  }
  stan_model_file_joint <- path.expand(stan_model_file_joint)
  message("Stan joint model file: ", stan_model_file_joint)
  stan_data_joint <- convert_data_from_models(
    fit_long = fit_long, fit_surv = fit_surv, surv_data = sim_dat$survival,
    k_bases = k_bases
  )
  safe_surv_beta_init <- function(stan_data, fallback_value = 0) {
    beta <- tryCatch(as.numeric(stan_data$beta_surv_init_mean), error = function(e) numeric(0))
    if (length(beta) != stan_data$q) {
      beta <- rep(fallback_value, stan_data$q)
    }
    beta[!is.finite(beta)] <- fallback_value
    beta
  }
  init_fun_joint <- function(sim_dat, stan_data, chains = 4, joint_model_type = "BP") {
    fallback <- list(beta_long = rep(0, stan_data$K_long), sigma_long = 10, sd_1_long = c(
      1,
      1
    ), beta_surv = rep(0, stan_data$q), gamma_scaled = rep(0.1, stan_data$m), alpha_tilde = 0)
    LMM <- nlme::lme(
      fixed = Y_obs ~ time + time_by_arm, random = ~ time | id, data = sim_dat$longitudinal,
      control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
    )
    fallback$beta_long <- as.numeric(nlme::fixef(LMM))
    fallback$sigma_long <- summary(LMM)$sigma
    vc <- nlme::VarCorr(LMM)
    fallback$sd_1_long <- sqrt(as.numeric(vc[1:2, "Variance"]))
    fallback$beta_surv <- safe_surv_beta_init(stan_data, fallback_value = 0)
    init_list <- lapply(1:chains, function(chain_id) {
      sd_jittered <- fallback$sd_1_long + abs(rnorm(2, 0, 0.05))
      sd_jittered <- pmax(sd_jittered, 0.001)
      rho <- 0.2
      R <- matrix(c(1, rho, rho, 1), 2, 2)
      L_1_long_init <- t(chol(R))
      init_vals <- list(
        z_1_long = matrix(rnorm(2 * stan_data$N_1_long), nrow = 2), L_1_long = L_1_long_init,
        beta_long = fallback$beta_long + rnorm(stan_data$K_long, 0, 0.1), sigma_long = abs(rnorm(
          1,
          fallback$sigma_long, 1
        )), sd_1_long = sd_jittered
      )
      init_vals$beta_surv <- rnorm(stan_data$q, fallback$beta_surv, 0.05)
      init_vals$alpha_tilde <- rnorm(1, 0, 0.1)
      init_vals$gamma <- abs(rnorm(stan_data$m, 0.1, 0.01))
      init_vals
    })
    return(init_list)
  }
  # Preserve per-chain jitter and its random-number draw order.
  init_values_joint <- init_fun_joint(
    sim_dat = sim_dat, stan_data = stan_data_joint, chains = chains,
    joint_model_type = joint_model_type
  )
  compact_init_values <- function(init_values) {
    lapply(seq_along(init_values), function(chain_id) {
      init_chain <- init_values[[chain_id]]
      keep_names <- intersect(c(
        "beta_long", "beta_surv", "alpha_tilde", "gamma", "sigma_long",
        "sd_1_long", "L_1_long", "gp_length_scale", "gp_marginal_sd", "f_gp_raw", "a",
        "log_b"
      ), names(init_chain))
      kept <- init_chain[keep_names]
      kept$chain <- chain_id
      kept
    })
  }
  init_info <- list(
    note = "Jittered Stan initial values by chain; z_1_long omitted to avoid storing a large random-effects matrix.",
    values = compact_init_values(init_values_joint)
  )
  # The recording wrapper saves exact inputs, starts and persistent draws.
  t_compile_joint_start <- Sys.time()
  mod_joint <- cmdstan_model(stan_model_file_joint,
    cpp_options = list(stan_threads = TRUE),
    force_recompile = FALSE
  )
  Sys.chmod(mod_joint$exe_file(), mode = "0755")
  if (!file.exists(mod_joint$exe_file())) {
    stop("Stan joint model compilation failed; binary not found.")
  }
  t_compile_joint_elapsed <- as.numeric(difftime(Sys.time(), t_compile_joint_start, units = "secs"))
  t_stan_joint_start <- Sys.time()
  fit_joint <- mod_joint$sample(
    data = stan_data_joint, seed = unique_seed, init = init_values_joint,
    chains = chains, parallel_chains = parallel_chains, threads_per_chain = threads_per_chain,
    iter_warmup = 1000, iter_sampling = 1000, max_treedepth = 10
  )
  t_stan_joint_elapsed <- as.numeric(difftime(Sys.time(), t_stan_joint_start, units = "secs"))
  loo_joint <- tryCatch(
    {
      fit_joint$loo(save_psis = FALSE)
    },
    error = function(e) {
      message("LOO computation failed for joint model: ", e$message)
    }
  )
  waic_joint <- tryCatch(
    {
      compute_waic_from_log_lik(fit_joint, log_lik_var = "log_lik")
    },
    error = function(e) {
      message("WAIC computation failed for joint model: ", e$message)
    }
  )
  summ_joint <- fit_joint$summary(
    variables = NULL, posterior::default_summary_measures(),
    quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975))
  )
  summ_joint$sim_id <- i
  summ_joint$model <- joint_model_type
  summ_joint$rhat <- fit_joint$summary()$rhat
  if (is.finite(longitudinal_response_center) && longitudinal_response_center != 0) {
    beta0_row <- summ_joint$variable == "beta_long[1]"
    shift_cols <- intersect(c("mean", "median", "q2.5", "q97.5"), names(summ_joint))
    if (any(beta0_row) && length(shift_cols) > 0) {
      summ_joint[beta0_row, shift_cols] <- lapply(
        summ_joint[beta0_row, shift_cols, drop = FALSE],
        function(x) x + longitudinal_response_center
      )
    }
    attr(summ_joint, "longitudinal_response_center") <- longitudinal_response_center
  }
  diag_joint <- fit_joint$diagnostic_summary()
  diag_joint$sim_id <- i
  diag_joint$model <- joint_model_type
  coefs <- NULL
  if (!is.null(fit_long)) {
    coefs <- nlme::fixef(fit_long)
    if (is.finite(longitudinal_response_center) && longitudinal_response_center != 0 &&
      "(Intercept)" %in% names(coefs)) {
      coefs["(Intercept)"] <- coefs["(Intercept)"] + longitudinal_response_center
      attr(coefs, "longitudinal_response_center") <- longitudinal_response_center
    }
  }
  runtimes <- data.frame(
    sim_id = i, joint_model = joint_model_type, k_bases = k_bases, lme_time = t_lme_elapsed,
    joint_compile_time = t_compile_joint_elapsed, joint_stan_time = t_stan_joint_elapsed,
    total_time = t_lme_elapsed + t_compile_joint_elapsed + t_stan_joint_elapsed
  )
  loo_joint_est <- if (!is.null(loo_joint)) {
    data.frame(elpd_loo = loo_joint$estimates["elpd_loo", "Estimate"], p_loo = loo_joint$estimates[
      "p_loo",
      "Estimate"
    ], looic = loo_joint$estimates["looic", "Estimate"])
  } else {
    data.frame(elpd_loo = NA, p_loo = NA, looic = NA)
  }
  waic_joint_est <- if (!is.null(waic_joint)) {
    data.frame(elpd_waic = waic_joint$estimates["elpd_waic", "Estimate"], p_waic = waic_joint$estimates[
      "p_waic",
      "Estimate"
    ], waic = waic_joint$estimates["waic", "Estimate"])
  } else {
    data.frame(elpd_waic = NA, p_waic = NA, waic = NA)
  }
  stan_data_field <- function(name, default = NULL) {
    if (name %in% names(stan_data_joint)) {
      stan_data_joint[[name]]
    } else {
      default
    }
  }
  as_numeric_unnamed <- function(x) {
    if (is.null(x)) {
      return(numeric(0))
    }
    unname(as.numeric(x))
  }
  s_long_saved <- as_numeric_unnamed(stan_data_field("s_long", NA_real_))[1]
  s_surv_saved <- as_numeric_unnamed(stan_data_field("s_surv", NA_real_))[1]
  alpha_tilde_sd_saved <- as_numeric_unnamed(stan_data_field("alpha_tilde_sd", NA_real_))[1]
  alpha_prior_sd_raw <- if (is.finite(s_long_saved) && is.finite(s_surv_saved) && is.finite(alpha_tilde_sd_saved) &&
    s_long_saved != 0) {
    alpha_tilde_sd_saved * s_surv_saved / s_long_saved
  } else {
    NA_real_
  }
  beta_long_prior_mean_centered <- as_numeric_unnamed(stan_data_field("beta_long_prior_mean"))
  beta_long_prior_mean_raw <- beta_long_prior_mean_centered
  if (length(beta_long_prior_mean_raw) >= 1 && is.finite(longitudinal_response_center) &&
    longitudinal_response_center != 0) {
    beta_long_prior_mean_raw[1] <- beta_long_prior_mean_raw[1] + longitudinal_response_center
  }
  prior_info <- list(
    longitudinal_response_center = longitudinal_response_center, alpha_prior_delta_y = alpha_prior_calibration()$delta_y,
    alpha_prior_max_f = alpha_prior_calibration()$max_f, alpha_prior_sd_target = alpha_prior_calibration()$sigma_alpha,
    s_long = s_long_saved, s_surv = s_surv_saved, alpha_tilde_prior_mean = 0, alpha_tilde_sd = alpha_tilde_sd_saved,
    alpha_prior_mean_raw = 0, alpha_prior_sd_raw = alpha_prior_sd_raw, beta_long_prior_mean = beta_long_prior_mean_centered,
    beta_long_prior_mean_raw = beta_long_prior_mean_raw, beta_long_prior_scale = as_numeric_unnamed(stan_data_field("beta_long_prior_scale")),
    beta_surv_prior_mean = as_numeric_unnamed(stan_data_field("beta_surv_prior_mean")),
    beta_surv_prior_scale = as_numeric_unnamed(stan_data_field("beta_surv_prior_scale")),
    beta_surv_init_mean = as_numeric_unnamed(stan_data_field("beta_surv_init_mean")), sd_b_prior_mean = as_numeric_unnamed(stan_data_field("sd_b_prior_mean")),
    sd_b_prior_scale = as_numeric_unnamed(stan_data_field("sd_b_prior_scale"))
  )
  logy_support <- NULL
  if (all(c("log_y_lower", "log_y_upper") %in% names(stan_data_joint))) {
    logy_support <- list(log_y_lower = stan_data_joint$log_y_lower, log_y_upper = stan_data_joint$log_y_upper)
  }
  baseline_recon <- NULL
  return(list(
    censor_prop = mean(sim_dat$survival$status == 0), unique_seed = unique_seed,
    longitudinal_response_center = longitudinal_response_center, baseline_label = baseline_label,
    aft_mode = aft_mode, weibull_shape = cfg_weibull_shape, weibull_scale = cfg_weibull_scale,
    loglogistic_shape = cfg_loglogistic_shape, loglogistic_scale = cfg_loglogistic_scale,
    joint_summary = summ_joint, joint_diagnostics = diag_joint, lme_coefs = coefs, logy_support = logy_support,
    baseline_recon = baseline_recon, runtimes = runtimes, prior_info = prior_info, init_info = init_info,
    loo_joint = loo_joint_est, waic_joint = waic_joint_est
  ))
}
