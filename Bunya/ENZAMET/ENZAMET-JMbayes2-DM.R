#!/usr/bin/env Rscript

# =============================================================================
# JMBAYES2 JOINT MODEL COMPARISON - ENZAMET DM SIMULATIONS - SLURM VERSION
# =============================================================================
#
# This script mirrors the ENZAMET-DM.R simulation/SLURM layout, but replaces the
# Stan sAFT-JM fit with JMbayes2 comparison fits. It intentionally keeps the
# same data-generation grid and seed convention so outputs can be paired with
# the completed Stan simulation results.
#
# Default grid:
#   4 baseline distributions x 5 scenarios x 2 censoring settings = 40 configs.
#
# Current intended use:
#   set SIM_BASELINE_FILTER=LL120_230 to run only the 10 log-logistic configs.
#
# Main controls:
#   SIM_JMB2_MODELS       Cox,Weibull,loglogistic, or comma-separated/all
#   SIM_JMB2_SAVE_MODE    summary, draws, or full. The summary mode keeps
#                         compact posterior summaries, including derived
#                         random-effect SD summaries from D posterior draws.
#   SIM_BASELINE_FILTER   e.g. LL120_230
#   SIM_BATCH_SIZE        replicates per batch
#   SIMS_PER_CONFIG       replicates per config
#
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(MASS)
  library(nlme)
  library(future)
  library(furrr)
  library(simsurv)
})

if (!requireNamespace("JMbayes2", quietly = TRUE)) {
  stop("Package 'JMbayes2' is required for ENZAMET-JMbayes2-DM.R.")
}

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file_arg) > 0) {
  sub("^--file=", "", script_file_arg[[1]])
} else {
  file.path("Bunya", "ENZAMET", "ENZAMET-JMbayes2-DM.R")
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < 1) {
    stop(name, " must be a positive integer.")
  }
  value
}

env_nonnegative_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < 0) {
    stop(name, " must be a non-negative integer.")
  }
  value
}

parse_jmb2_models <- function(value) {
  value <- trimws(value)
  if (!nzchar(value) || tolower(value) == "all") {
    return(c("Cox", "Weibull", "loglogistic"))
  }

  pieces <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  pieces <- pieces[nzchar(pieces)]
  aliases <- c(
    cox = "Cox",
    coxph = "Cox",
    ph = "Cox",
    weibull = "Weibull",
    weib = "Weibull",
    loglogistic = "loglogistic",
    log_logistic = "loglogistic",
    llogis = "loglogistic",
    ll = "loglogistic"
  )

  out <- unname(aliases[tolower(pieces)])
  if (any(is.na(out))) {
    stop(
      "Unknown SIM_JMB2_MODELS value(s): ",
      paste(pieces[is.na(out)], collapse = ", "),
      ". Use Cox, Weibull, loglogistic, or all."
    )
  }
  unique(out)
}

selected_jmb2_models <- parse_jmb2_models(Sys.getenv("SIM_JMB2_MODELS", unset = "all"))
selected_model_label <- if (setequal(selected_jmb2_models, c("Cox", "Weibull", "loglogistic"))) {
  "All"
} else {
  paste(selected_jmb2_models, collapse = "_")
}
output_model_prefix <- paste0("JMbayes2_", selected_model_label)

# Keep JMbayes2's default iteration/burn-in/thinning schedule for posterior
# summaries so the comparator is not disadvantaged by a shorter MCMC run. The
# production launcher overrides only n_chains to 4, matching the Stan sAFT-JM
# diagnostic convention while otherwise retaining JMbayes2's MCMC defaults.
jmb2_n_iter <- env_int("SIM_JMB2_N_ITER", 3500)
jmb2_n_burnin <- env_int("SIM_JMB2_N_BURNIN", 500)
jmb2_n_thin <- env_int("SIM_JMB2_N_THIN", 1)
jmb2_n_chains <- env_int("SIM_JMB2_N_CHAINS", 3)
jmb2_cores <- env_int("SIM_JMB2_CORES", jmb2_n_chains)
jmb2_retry_cores <- env_nonnegative_int("SIM_JMB2_RETRY_CORES", 0)
jmb2_save_mode <- tolower(Sys.getenv("SIM_JMB2_SAVE_MODE", unset = "summary"))
if (!jmb2_save_mode %in% c("summary", "draws", "full")) {
  stop("SIM_JMB2_SAVE_MODE must be one of: summary, draws, full.")
}
jmb2_draws_pattern <- Sys.getenv(
  "SIM_JMB2_DRAWS_PATTERN",
  unset = "beta|alpha|sigma|D|tau|scale|shape"
)
jmb2_full_fit_dir <- path.expand(Sys.getenv("SIM_JMB2_FULL_FIT_DIR", unset = ""))

cpus_per_sim <- max(1L, jmb2_cores)

standardise_summary_table <- function(x, component, model, sim_id) {
  if (is.null(x)) return(tibble::tibble())

  if (is.list(x) && !is.data.frame(x) && !is.matrix(x)) {
    return(purrr::imap_dfr(x, function(xx, nm) {
      standardise_summary_table(xx, paste(component, nm, sep = "::"), model, sim_id)
    }))
  }

  df <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(tibble::tibble())

  if (!"parameter" %in% names(df)) {
    rn <- rownames(df)
    if (is.null(rn) || any(!nzchar(rn))) rn <- as.character(seq_len(nrow(df)))
    df$parameter <- rn
  }

  clean_names <- gsub("[^a-z0-9]+", "", tolower(names(df)))
  get_col <- function(patterns) {
    hit <- which(clean_names %in% patterns)
    if (length(hit) > 0) return(names(df)[hit[[1]]])
    hit <- which(vapply(patterns, function(p) any(grepl(p, clean_names)), logical(1)))
    if (length(hit) > 0) {
      col_hit <- grep(patterns[[hit[[1]]]], clean_names)
      return(names(df)[col_hit[[1]]])
    }
    NA_character_
  }

  mean_col <- get_col(c("mean", "posteriormean", "estimate", "value"))
  sd_col <- get_col(c("sd", "stdev", "stddev", "stderr", "standarderror", "se"))
  median_col <- get_col(c("median", "med"))
  q025_col <- get_col(c("25", "025", "x25", "q25", "ci25", "lower", "lcl", "x0025", "q0025"))
  q975_col <- get_col(c("975", "x975", "q975", "ci975", "upper", "ucl", "x0975", "q0975"))
  rhat_col <- get_col(c("rhat", "rhats"))
  ess_col <- get_col(c("ess", "essbulk", "neff"))

  as_num_col <- function(col) {
    if (is.na(col) || !col %in% names(df)) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(df[[col]]))
  }

  tibble::tibble(
    model = model,
    component = component,
    parameter = as.character(df$parameter),
    mean = as_num_col(mean_col),
    median = as_num_col(median_col),
    sd = as_num_col(sd_col),
    q2.5 = as_num_col(q025_col),
    q97.5 = as_num_col(q975_col),
    rhat = as_num_col(rhat_col),
    ess = as_num_col(ess_col),
    sim_id = sim_id
  )
}

as_mcmc_chain_matrices <- function(x) {
  if (is.null(x)) return(list())

  if (inherits(x, "mcmc.list")) {
    return(lapply(x, function(xx) as.matrix(xx)))
  }
  if (inherits(x, "mcmc")) {
    return(list(as.matrix(x)))
  }
  if (is.matrix(x) || is.data.frame(x)) {
    return(list(as.matrix(x)))
  }
  if (is.array(x) && length(dim(x)) == 2) {
    return(list(as.matrix(x)))
  }
  if (is.list(x)) {
    out <- purrr::map(x, as_mcmc_chain_matrices)
    return(purrr::compact(unlist(out, recursive = FALSE)))
  }

  list()
}

find_D_diag_column <- function(D_mat, diag_id) {
  nms <- colnames(D_mat)
  if (!is.null(nms)) {
    pattern <- sprintf("^D\\[\\s*%d\\s*,\\s*%d\\s*\\]$", diag_id, diag_id)
    hit <- grep(pattern, nms)
    if (length(hit) > 0) return(hit[[1]])
  }

  p <- ncol(D_mat)
  q <- (sqrt(8 * p + 1) - 1) / 2
  if (!is.finite(q) || abs(q - round(q)) > sqrt(.Machine$double.eps)) {
    return(NA_integer_)
  }
  q <- as.integer(round(q))
  if (diag_id > q) return(NA_integer_)

  idx <- lower.tri(matrix(NA_real_, q, q), diag = TRUE)
  diag_positions <- which(row(idx)[idx] == col(idx)[idx])
  diag_positions[[diag_id]]
}

summarise_draws <- function(draws, chain_draws = NULL) {
  draws <- suppressWarnings(as.numeric(draws))
  draws <- draws[is.finite(draws)]
  if (length(draws) == 0) {
    return(list(
      mean = NA_real_,
      median = NA_real_,
      sd = NA_real_,
      q2.5 = NA_real_,
      q97.5 = NA_real_,
      rhat = NA_real_,
      ess = NA_real_
    ))
  }

  rhat <- NA_real_
  ess <- NA_real_
  if (requireNamespace("coda", quietly = TRUE) &&
      length(chain_draws) > 0 &&
      length(chain_draws) >= 2) {
    chain_draws <- lapply(chain_draws, function(x) {
      x <- suppressWarnings(as.numeric(x))
      x[is.finite(x)]
    })
    chain_draws <- chain_draws[vapply(chain_draws, length, integer(1)) > 1]
    if (length(chain_draws) >= 2) {
      mcmc_list <- coda::mcmc.list(lapply(chain_draws, coda::mcmc))
      rhat <- tryCatch({
        as.numeric(coda::gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1, "Point est."])
      }, error = function(e) NA_real_)
      ess <- tryCatch({
        sum(as.numeric(coda::effectiveSize(mcmc_list)), na.rm = TRUE)
      }, error = function(e) NA_real_)
    }
  }

  list(
    mean = mean(draws),
    median = stats::median(draws),
    sd = stats::sd(draws),
    q2.5 = unname(stats::quantile(draws, 0.025, na.rm = TRUE)),
    q97.5 = unname(stats::quantile(draws, 0.975, na.rm = TRUE)),
    rhat = rhat,
    ess = ess
  )
}

jmb2_random_sd_summary <- function(fit, model, sim_id) {
  D_chains <- as_mcmc_chain_matrices(fit$mcmc$D)
  if (length(D_chains) == 0 || !all(vapply(D_chains, ncol, integer(1)) >= 1)) {
    return(tibble::tibble())
  }

  purrr::map_dfr(seq_len(2), function(diag_id) {
    col_id <- find_D_diag_column(D_chains[[1]], diag_id)
    if (is.na(col_id)) return(NULL)

    chain_draws <- lapply(D_chains, function(chain) {
      if (ncol(chain) < col_id) return(numeric())
      sqrt(pmax(suppressWarnings(as.numeric(chain[, col_id])), 0))
    })
    draws <- unlist(chain_draws, use.names = FALSE)
    s <- summarise_draws(draws, chain_draws)

    tibble::tibble(
      model = model,
      component = "DerivedRandomEffects",
      parameter = sprintf("sqrt(D[%d,%d])", diag_id, diag_id),
      variable = sprintf("sd_1_long[%d]", diag_id),
      mean = s$mean,
      median = s$median,
      sd = s$sd,
      q2.5 = s$q2.5,
      q97.5 = s$q97.5,
      rhat = s$rhat,
      ess = s$ess,
      sim_id = sim_id,
      source = paste0(model, "-JM")
    )
  })
}

map_parameter_to_common <- function(component, parameter) {
  comp <- tolower(component)
  par <- as.character(parameter)
  par_clean <- tolower(gsub("[^a-z0-9_]+", "", par))

  dplyr::case_when(
    grepl("longitudinal|outcome", comp) & par %in% c("(Intercept)", "Intercept", "intercept") ~ "beta_long[1]",
    grepl("longitudinal|outcome", comp) & par_clean %in% c("time") ~ "beta_long[2]",
    grepl("longitudinal|outcome", comp) & par_clean %in% c("time_by_arm", "timebyarm", "timearm") ~ "beta_long[3]",
    grepl("longitudinal|outcome", comp) & grepl("sigma", par_clean) ~ "sigma_long",
    grepl("survival", comp) & par_clean %in% c("arm", "arm1") ~ "beta_surv[1]",
    grepl("association|alpha", comp) | grepl("alpha|value", par_clean) ~ "alpha",
    TRUE ~ NA_character_
  )
}

tidy_jmb2_summary <- function(fit, model, sim_id, longitudinal_response_center = 0) {
  s <- summary(fit)
  random_sd_tidy <- jmb2_random_sd_summary(fit, model, sim_id)
  raw <- purrr::imap_dfr(s, function(x, nm) {
    standardise_summary_table(x, nm, model, sim_id)
  })

  if (nrow(raw) == 0) {
    return(list(raw_summary = s, tidy_summary = random_sd_tidy))
  }

  tidy <- raw %>%
    mutate(
      variable = map_parameter_to_common(component, parameter),
      source = paste0(model, "-JM")
    ) %>%
    filter(!is.na(variable))

  if (nrow(tidy) > 0 &&
      is.finite(longitudinal_response_center) &&
      longitudinal_response_center != 0) {
    beta0_row <- tidy$variable == "beta_long[1]"
    shift_cols <- intersect(c("mean", "median", "q2.5", "q97.5"), names(tidy))
    if (any(beta0_row) && length(shift_cols) > 0) {
      tidy[beta0_row, shift_cols] <- lapply(
        tidy[beta0_row, shift_cols, drop = FALSE],
        function(x) x + longitudinal_response_center
      )
    }
  }

  tidy <- dplyr::bind_rows(tidy, random_sd_tidy) %>%
    dplyr::distinct(.data$source, .data$variable, .data$parameter, .data$component, .keep_all = TRUE)

  list(raw_summary = s, tidy_summary = tidy)
}

tidy_lme_summary <- function(fit_long, sim_id, longitudinal_response_center = 0) {
  ttab <- as.data.frame(summary(fit_long)$tTable)
  ttab$parameter <- rownames(ttab)

  interval_df <- tryCatch({
    ints <- nlme::intervals(fit_long, which = "fixed")$fixed
    out <- as.data.frame(ints)
    out$parameter <- rownames(ints)
    names(out) <- sub("^lower$", "q2.5", names(out))
    names(out) <- sub("^est\\.$", "mean_interval", names(out))
    names(out) <- sub("^upper$", "q97.5", names(out))
    out
  }, error = function(e) NULL)

  fixed <- ttab %>%
    transmute(
      model = "LMM",
      component = "Longitudinal",
      parameter = parameter,
      variable = dplyr::case_when(
        parameter == "(Intercept)" ~ "beta_long[1]",
        parameter == "time" ~ "beta_long[2]",
        parameter == "time_by_arm" ~ "beta_long[3]",
        TRUE ~ NA_character_
      ),
      mean = .data$Value,
      median = NA_real_,
      sd = .data$Std.Error,
      sim_id = sim_id,
      source = "LMM"
    ) %>%
    filter(!is.na(variable))

  if (!is.null(interval_df)) {
    fixed <- fixed %>%
      dplyr::left_join(
        interval_df %>% dplyr::select(parameter, q2.5, q97.5),
        by = "parameter"
      )
  } else {
    fixed$q2.5 <- NA_real_
    fixed$q97.5 <- NA_real_
  }

  if (is.finite(longitudinal_response_center) &&
      longitudinal_response_center != 0) {
    beta0_row <- fixed$variable == "beta_long[1]"
    fixed[beta0_row, c("mean", "q2.5", "q97.5")] <- lapply(
      fixed[beta0_row, c("mean", "q2.5", "q97.5"), drop = FALSE],
      function(x) x + longitudinal_response_center
    )
  }

  vc <- tryCatch(nlme::VarCorr(fit_long), error = function(e) NULL)
  random_sd <- tibble::tibble()
  if (!is.null(vc)) {
    sd_vals <- suppressWarnings(sqrt(as.numeric(vc[c("(Intercept)", "time"), "Variance"])))
    random_sd <- tibble::tibble(
      model = "LMM",
      component = "Longitudinal",
      parameter = c("sd_(Intercept)", "sd_time"),
      variable = c("sd_1_long[1]", "sd_1_long[2]"),
      mean = sd_vals,
      median = NA_real_,
      sd = NA_real_,
      q2.5 = NA_real_,
      q97.5 = NA_real_,
      sim_id = sim_id,
      source = "LMM"
    )
  }

  sigma_row <- tibble::tibble(
    model = "LMM",
    component = "Longitudinal",
    parameter = "sigma",
    variable = "sigma_long",
    mean = as.numeric(summary(fit_long)$sigma),
    median = NA_real_,
    sd = NA_real_,
    q2.5 = NA_real_,
    q97.5 = NA_real_,
    sim_id = sim_id,
    source = "LMM"
  )

  dplyr::bind_rows(fixed, random_sd, sigma_row) %>%
    dplyr::select(model, source, component, parameter, variable, mean, median, sd, q2.5, q97.5, sim_id)
}

empty_initial_survreg_tidy <- function() {
  tibble::tibble(
    model = character(),
    source = character(),
    component = character(),
    parameter = character(),
    variable = character(),
    mean = numeric(),
    median = numeric(),
    sd = numeric(),
    q2.5 = numeric(),
    q97.5 = numeric(),
    sim_id = integer()
  )
}

tidy_initial_survreg <- function(fit_surv, model, sim_id) {
  if (!inherits(fit_surv, "survreg")) {
    return(empty_initial_survreg_tidy())
  }

  coefs <- stats::coef(fit_surv)
  vc <- tryCatch(stats::vcov(fit_surv), error = function(e) NULL)

  se_for <- function(name) {
    if (is.null(vc) || is.null(dimnames(vc)) || !name %in% rownames(vc)) {
      return(NA_real_)
    }
    out <- suppressWarnings(sqrt(vc[name, name]))
    ifelse(is.finite(out), out, NA_real_)
  }

  normal_interval <- function(est, se) {
    if (!is.finite(est) || !is.finite(se)) return(c(NA_real_, NA_real_))
    est + c(-1, 1) * stats::qnorm(0.975) * se
  }

  lognormal_interval <- function(log_est, se_log) {
    if (!is.finite(log_est) || !is.finite(se_log)) return(c(NA_real_, NA_real_))
    exp(log_est + c(-1, 1) * stats::qnorm(0.975) * se_log)
  }

  inverse_lognormal_interval <- function(log_est, se_log) {
    if (!is.finite(log_est) || !is.finite(se_log)) return(c(NA_real_, NA_real_))
    sort(exp(-(log_est + c(-1, 1) * stats::qnorm(0.975) * se_log)))
  }

  rows <- list()
  source <- paste0(model, "-survreg-initial")

  if ("arm" %in% names(coefs)) {
    est <- unname(coefs["arm"])
    se <- se_for("arm")
    ci <- normal_interval(est, se)
    rows$arm <- tibble::tibble(
      model = model,
      source = source,
      component = "InitialSurvreg",
      parameter = "arm",
      variable = "beta_surv[1]",
      mean = est,
      median = NA_real_,
      sd = se,
      q2.5 = ci[[1]],
      q97.5 = ci[[2]],
      sim_id = sim_id
    )
  }

  if ("(Intercept)" %in% names(coefs)) {
    log_scale_est <- unname(coefs["(Intercept)"])
    se_log_scale_est <- se_for("(Intercept)")
    scale_ci <- lognormal_interval(log_scale_est, se_log_scale_est)
    rows$baseline_scale <- tibble::tibble(
      model = model,
      source = source,
      component = "InitialSurvreg",
      parameter = "exp((Intercept))",
      variable = "baseline_scale",
      mean = exp(log_scale_est),
      median = NA_real_,
      sd = exp(log_scale_est) * se_log_scale_est,
      q2.5 = scale_ci[[1]],
      q97.5 = scale_ci[[2]],
      sim_id = sim_id
    )
  }

  log_scale_param <- log(fit_surv$scale)
  log_scale_name <- character()
  if (!is.null(vc) && !is.null(rownames(vc))) {
    log_scale_name <- setdiff(rownames(vc), names(coefs))
  }
  se_log_scale_param <- if (length(log_scale_name) == 1) {
    se_for(log_scale_name[[1]])
  } else {
    NA_real_
  }
  shape_ci <- inverse_lognormal_interval(log_scale_param, se_log_scale_param)
  rows$baseline_shape <- tibble::tibble(
    model = model,
    source = source,
    component = "InitialSurvreg",
    parameter = "1 / survreg_scale",
    variable = "baseline_shape",
    mean = 1 / fit_surv$scale,
    median = NA_real_,
    sd = (1 / fit_surv$scale) * se_log_scale_param,
    q2.5 = shape_ci[[1]],
    q97.5 = shape_ci[[2]],
    sim_id = sim_id
  )

  dplyr::bind_rows(rows)
}

compact_jmb2_draws <- function(fit, pattern) {
  mcmc_obj <- fit$mcmc
  if (is.null(mcmc_obj)) return(NULL)
  if (is.list(mcmc_obj) && !is.null(names(mcmc_obj))) {
    keep <- grepl(pattern, names(mcmc_obj), ignore.case = TRUE)
    return(mcmc_obj[keep])
  }
  mcmc_obj
}

save_full_fit_if_requested <- function(fit, model, unique_seed, sim_id) {
  if (!nzchar(jmb2_full_fit_dir)) return(fit)

  if (!dir.exists(jmb2_full_fit_dir)) {
    dir.create(jmb2_full_fit_dir, recursive = TRUE, showWarnings = FALSE)
  }
  path <- file.path(
    jmb2_full_fit_dir,
    sprintf("JMbayes2_%s_seed%s_sim%s.rds", model, unique_seed, sim_id)
  )
  saveRDS(fit, path)
  list(full_fit_path = path)
}

fit_one_jmb2_model <- function(model, fit_long, surv_data, unique_seed, sim_id,
                               longitudinal_response_center = 0) {
  status <- "ok"
  error_message <- NA_character_
  fit_surv <- NULL
  fit_jm <- NULL
  retry_used <- FALSE
  retry_error_message <- NA_character_

  t_surv_start <- Sys.time()
  fit_surv <- tryCatch({
    if (model == "Cox") {
      survival::coxph(
        Surv(time, status) ~ arm,
        data = surv_data,
        x = TRUE,
        model = TRUE
      )
    } else if (model == "Weibull") {
      survival::survreg(
        Surv(time, status) ~ arm,
        data = surv_data,
        dist = "weibull",
        x = TRUE,
        y = TRUE,
        model = TRUE
      )
    } else if (model == "loglogistic") {
      survival::survreg(
        Surv(time, status) ~ arm,
        data = surv_data,
        dist = "loglogistic",
        x = TRUE,
        y = TRUE,
        model = TRUE
      )
    } else {
      stop("Unknown JMbayes2 comparison model: ", model)
    }
  }, error = function(e) {
    status <<- "survival_fit_failed"
    error_message <<- conditionMessage(e)
    NULL
  })
  t_surv_elapsed <- as.numeric(difftime(Sys.time(), t_surv_start, units = "secs"))

  if (is.null(fit_surv)) {
    return(list(
      model = model,
      status = status,
      error_message = error_message,
      tidy_summary = tibble::tibble(),
      raw_summary = NULL,
      initial_survreg_tidy = empty_initial_survreg_tidy(),
      object_size_mb = NA_real_,
      saved_fit = NULL,
      selected_draws = NULL,
      runtime = tibble::tibble(
        model = model,
        survival_fit_time = t_surv_elapsed,
        jmb2_time = NA_real_,
        total_model_time = t_surv_elapsed,
        object_size_mb = NA_real_,
        status = status
      )
    ))
  }
  initial_survreg_tidy <- tidy_initial_survreg(fit_surv, model, sim_id)

  fit_jm_once <- function(cores) {
    suppressMessages(
      JMbayes2::jm(
        fit_surv,
        list(fit_long),
        time_var = "time",
        data_Surv = surv_data,
        id_var = "id",
        control = list(
          n_iter = jmb2_n_iter,
          n_burnin = jmb2_n_burnin,
          n_thin = jmb2_n_thin,
          n_chains = jmb2_n_chains,
          cores = cores,
          seed = unique_seed + match(model, c("Cox", "Weibull", "loglogistic")),
          save_random_effects = FALSE,
          save_logLik_contributions = FALSE
        )
      )
    )
  }

  t_jmb2_start <- Sys.time()
  fit_jm <- tryCatch({
    fit_jm_once(jmb2_cores)
  }, error = function(e) {
    status <<- "jm_fit_failed"
    error_message <<- conditionMessage(e)
    NULL
  })

  if (is.null(fit_jm) &&
      jmb2_retry_cores > 0L &&
      jmb2_retry_cores != jmb2_cores) {
    message(sprintf(
      "JMbayes2 %s fit failed with cores=%d; retrying with cores=%d. First error: %s",
      model, jmb2_cores, jmb2_retry_cores, error_message
    ))
    retry_used <- TRUE
    first_error_message <- error_message
    fit_jm <- tryCatch({
      status <<- "ok"
      error_message <<- NA_character_
      fit_jm_once(jmb2_retry_cores)
    }, error = function(e) {
      status <<- "jm_fit_failed"
      retry_error_message <<- conditionMessage(e)
      error_message <<- paste(
        "initial:", first_error_message,
        "| retry:", retry_error_message
      )
      NULL
    })
  }
  t_jmb2_elapsed <- as.numeric(difftime(Sys.time(), t_jmb2_start, units = "secs"))

  if (is.null(fit_jm)) {
    return(list(
      model = model,
      status = status,
      error_message = error_message,
      tidy_summary = tibble::tibble(),
      raw_summary = NULL,
      initial_survreg_tidy = initial_survreg_tidy,
      object_size_mb = NA_real_,
      saved_fit = NULL,
      selected_draws = NULL,
      retry_used = retry_used,
      retry_cores = if (retry_used) jmb2_retry_cores else NA_integer_,
      retry_error_message = retry_error_message,
      runtime = tibble::tibble(
        model = model,
        survival_fit_time = t_surv_elapsed,
        jmb2_time = t_jmb2_elapsed,
        total_model_time = t_surv_elapsed + t_jmb2_elapsed,
        object_size_mb = NA_real_,
        status = status
      )
    ))
  }

  object_size_mb <- as.numeric(utils::object.size(fit_jm)) / 1024^2
  summary_bundle <- tidy_jmb2_summary(
    fit_jm,
    model = model,
    sim_id = sim_id,
    longitudinal_response_center = longitudinal_response_center
  )

  saved_fit <- NULL
  selected_draws <- NULL
  if (jmb2_save_mode == "draws") {
    selected_draws <- compact_jmb2_draws(fit_jm, jmb2_draws_pattern)
  } else if (jmb2_save_mode == "full") {
    saved_fit <- save_full_fit_if_requested(fit_jm, model, unique_seed, sim_id)
  }

  list(
    model = model,
    status = status,
    error_message = error_message,
    tidy_summary = summary_bundle$tidy_summary,
    raw_summary = summary_bundle$raw_summary,
    initial_survreg_tidy = initial_survreg_tidy,
    object_size_mb = object_size_mb,
    saved_fit = saved_fit,
    selected_draws = selected_draws,
    retry_used = retry_used,
    retry_cores = if (retry_used) jmb2_retry_cores else NA_integer_,
    retry_error_message = retry_error_message,
    runtime = tibble::tibble(
      model = model,
      survival_fit_time = t_surv_elapsed,
      jmb2_time = t_jmb2_elapsed,
      total_model_time = t_surv_elapsed + t_jmb2_elapsed,
      object_size_mb = object_size_mb,
      status = status
    )
  )
}

run_simulation <- function(i, n_patients, lambda_c, aft_mode, max_FU, config_hash,
                           k_bases = 5, scenario = 2, baseline_label = aft_mode,
                           cfg_weibull_shape,
                           cfg_weibull_scale,
                           cfg_loglogistic_shape,
                           cfg_loglogistic_scale) {
  force(i)
  force(n_patients)
  force(lambda_c)
  force(aft_mode)
  force(max_FU)
  force(config_hash)
  force(k_bases)
  force(scenario)
  force(baseline_label)
  force(cfg_weibull_shape)
  force(cfg_weibull_scale)
  force(cfg_loglogistic_shape)
  force(cfg_loglogistic_scale)

  scenario_file <- file.path(
    repo_root,
    "Bunya", "ENZAMET",
    sprintf("Data-Gen-ENZAMET-%02d.R", scenario)
  )
  source(scenario_file, local = TRUE)

  unique_seed <- 1e6 * config_hash + i
  message("unique_seed: ", unique_seed)

  sim_dat <- simulate_joint_dataset(
    D = D,
    beta_0 = beta_0,
    beta_1 = beta_1,
    beta_2 = beta_2,
    sigma_e = sigma_e,
    log_HR = log_HR,
    log_AF = log_AF,
    alpha_PH = 99,
    alpha_AFT = alpha_AFT,
    weibull_shape = cfg_weibull_shape,
    weibull_scale = cfg_weibull_scale,
    loglogistic_shape = cfg_loglogistic_shape,
    loglogistic_scale = cfg_loglogistic_scale,
    visit = visit,
    seed = unique_seed,
    n_patients = n_patients,
    max_FU = max_FU,
    lambda_c = lambda_c,
    aft_mode = aft_mode,
    link_type = "value"
  )

  sim_dat$longitudinal <- sim_dat$longitudinal %>%
    mutate(time_by_arm = time * arm)

  longitudinal_response_center <- mean(sim_dat$longitudinal$Y_obs, na.rm = TRUE)
  if (!is.finite(longitudinal_response_center)) {
    stop("Observed longitudinal response mean is not finite; cannot center Y_obs.")
  }
  sim_dat$longitudinal$Y_obs <- sim_dat$longitudinal$Y_obs -
    longitudinal_response_center
  message(sprintf(
    "Centered longitudinal response by subtracting observed mean %.6f before fitting.",
    longitudinal_response_center
  ))

  t_lme_start <- Sys.time()
  fit_long <- nlme::lme(
    fixed = Y_obs ~ time + time_by_arm,
    random = ~ time | id,
    data = sim_dat$longitudinal,
    na.action = na.omit,
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
  )
  t_lme_elapsed <- as.numeric(difftime(Sys.time(), t_lme_start, units = "secs"))

  sim_dat$survival <- sim_dat$survival %>%
    dplyr::rename(time = T_obs)

  lmm_tidy <- tidy_lme_summary(
    fit_long = fit_long,
    sim_id = i,
    longitudinal_response_center = longitudinal_response_center
  )

  model_results <- purrr::map(
    selected_jmb2_models,
    fit_one_jmb2_model,
    fit_long = fit_long,
    surv_data = sim_dat$survival,
    unique_seed = unique_seed,
    sim_id = i,
    longitudinal_response_center = longitudinal_response_center
  )
  names(model_results) <- selected_jmb2_models

  jmb2_tidy <- purrr::map_dfr(model_results, "tidy_summary")
  comparison_tidy <- bind_rows(lmm_tidy, jmb2_tidy)

  runtimes <- bind_rows(
    tibble::tibble(
      model = "LMM",
      survival_fit_time = NA_real_,
      jmb2_time = NA_real_,
      total_model_time = t_lme_elapsed,
      object_size_mb = as.numeric(utils::object.size(fit_long)) / 1024^2,
      status = "ok"
    ),
    purrr::map_dfr(model_results, "runtime")
  ) %>%
    mutate(
      sim_id = i,
      unique_seed = unique_seed,
      lme_time = t_lme_elapsed,
      selected_model_label = selected_model_label
    )

  model_status <- purrr::imap_dfr(model_results, function(x, model_name) {
    tibble::tibble(
      model = model_name,
      status = x$status,
      error_message = x$error_message,
      retry_used = isTRUE(x$retry_used),
      retry_cores = if (is.null(x$retry_cores)) NA_integer_ else x$retry_cores,
      retry_error_message = if (is.null(x$retry_error_message)) NA_character_ else x$retry_error_message,
      object_size_mb = x$object_size_mb
    )
  })

  list(
    censor_prop = mean(sim_dat$survival$status == 0),
    unique_seed = unique_seed,
    longitudinal_response_center = longitudinal_response_center,
    baseline_label = baseline_label,
    aft_mode = aft_mode,
    weibull_shape = cfg_weibull_shape,
    weibull_scale = cfg_weibull_scale,
    loglogistic_shape = cfg_loglogistic_shape,
    loglogistic_scale = cfg_loglogistic_scale,
    k_bases = k_bases,
    selected_jmb2_models = selected_jmb2_models,
    jmb2_save_mode = jmb2_save_mode,
    comparison_tidy = comparison_tidy,
    lme_coefs = setNames(
      lmm_tidy$mean[match(c("beta_long[1]", "beta_long[2]", "beta_long[3]"), lmm_tidy$variable)],
      c("(Intercept)", "time", "time_by_arm")
    ),
    jmb2_results = model_results,
    model_status = model_status,
    runtimes = runtimes
  )
}

# =============================================================================
# SLURM job mapping and configuration grid
# =============================================================================

array_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
array_task_offset <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_OFFSET", unset = "0"))
sim_id <- array_task_id + array_task_offset
message(sprintf(
  "SLURM_ARRAY_TASK_ID: %d, SLURM_ARRAY_TASK_OFFSET: %d, global sim_id: %d",
  array_task_id, array_task_offset, sim_id
))
if (array_task_id == 0) stop("SLURM_ARRAY_TASK_ID not found.")

baseline_grid <- data.frame(
  baseline_label = c("LL120_230", "W090_380", "W130_380", "W210_380"),
  aft_mode = c("loglogistic", "Weibull", "Weibull", "Weibull"),
  weibull_shape = c(NA_real_, 0.90, 1.30, 2.10),
  weibull_scale = c(NA_real_, 38.0, 38.0, 38.0),
  loglogistic_shape = c(1.20, NA_real_, NA_real_, NA_real_),
  loglogistic_scale = c(23.0, NA_real_, NA_real_, NA_real_),
  stringsAsFactors = FALSE
)

config_base <- merge(
  expand.grid(
    n_patients = c(1100),
    k_bases = c(5),
    scenario = c(1, 2, 3, 4, 5),
    max_FU = 120,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ),
  baseline_grid,
  by = NULL
)

lambda_c_table <- data.frame(
  baseline_label = rep(baseline_grid$baseline_label, each = 5),
  scenario = rep(1:5, times = nrow(baseline_grid)),
  lambda_c = c(
    0.009948730, 0.003906250, 0.004272461, 0.017028809, 0.017211914,
    0.009643555, 0.003479004, 0.003784180, 0.017211914, 0.017456055,
    0.008666992, 0.002227783, 0.002593994, 0.015563965, 0.015716553,
    0.008056641, 0.001098633, 0.001464844, 0.014099121, 0.014160156
  ),
  stringsAsFactors = FALSE
)

cfg_default <- merge(
  config_base,
  lambda_c_table,
  by = c("baseline_label", "scenario"),
  all.x = TRUE,
  sort = FALSE
)
cfg_default$cens_label <- "50"

cfg_admin <- transform(
  config_base,
  lambda_c = -1,
  cens_label = "admin"
)

config_grid <- rbind(cfg_default, cfg_admin)
config_grid <- config_grid[order(
  match(config_grid$baseline_label, baseline_grid$baseline_label),
  config_grid$scenario,
  match(config_grid$cens_label, c("50", "admin"))
), ]
row.names(config_grid) <- NULL
config_grid$full_config_id <- seq_len(nrow(config_grid))

baseline_filter <- Sys.getenv("SIM_BASELINE_FILTER", unset = "")
if (nzchar(baseline_filter)) {
  keep_baselines <- trimws(strsplit(baseline_filter, ",", fixed = TRUE)[[1]])
  keep_baselines <- keep_baselines[nzchar(keep_baselines)]
  if (length(keep_baselines) == 0) {
    stop("SIM_BASELINE_FILTER did not contain any valid labels.")
  }
  config_grid <- config_grid[config_grid$baseline_label %in% keep_baselines, , drop = FALSE]
  if (nrow(config_grid) == 0) {
    stop(sprintf(
      "SIM_BASELINE_FILTER matched no configs. Requested: %s; available: %s",
      paste(keep_baselines, collapse = ", "),
      paste(unique(baseline_grid$baseline_label), collapse = ", ")
    ))
  }
  row.names(config_grid) <- NULL
  message(sprintf(
    "Applied SIM_BASELINE_FILTER=%s; retained %d configs.",
    paste(keep_baselines, collapse = ","),
    nrow(config_grid)
  ))
}

batch_size <- env_int("SIM_BATCH_SIZE", 1)
sims_per_config <- env_int("SIMS_PER_CONFIG", 10)
if (sims_per_config %% batch_size != 0) {
  stop(sprintf(
    "SIMS_PER_CONFIG (%d) must be divisible by SIM_BATCH_SIZE (%d).",
    sims_per_config, batch_size
  ))
}
batches_per_config <- sims_per_config %/% batch_size

config_id <- ((sim_id - 1) %/% batches_per_config) + 1
batch_id <- ((sim_id - 1) %% batches_per_config) + 1
total_array_tasks <- nrow(config_grid) * batches_per_config
message(sprintf(
  "Configured simulation: %d configs x %d batches/config = %d array tasks; batch_size=%d, sims_per_config=%d",
  nrow(config_grid), batches_per_config, total_array_tasks, batch_size, sims_per_config
))
if (sim_id > total_array_tasks) {
  stop(sprintf("SLURM_ARRAY_TASK_ID %d exceeds configured array size %d.", sim_id, total_array_tasks))
}

config <- config_grid[config_id, ]
n_patients <- config$n_patients
lambda_c <- config$lambda_c
aft_mode <- config$aft_mode
baseline_label <- config$baseline_label
weibull_shape <- config$weibull_shape
weibull_scale <- config$weibull_scale
loglogistic_shape <- config$loglogistic_shape
loglogistic_scale <- config$loglogistic_scale
max_FU <- config$max_FU
k_bases <- config$k_bases
scenario <- config$scenario
cens_label <- config$cens_label
trial_ids <- ((batch_id - 1) * batch_size + 1):(batch_id * batch_size)
config_hash <- config$full_config_id

message(sprintf(
  "Config ID: %d, Config Hash: %d, n_patients: %d, lambda_c: %.6f, baseline: %s, aft_mode: %s, k_bases: %d, scenario: %d, max_FU: %d, models: %s",
  config_id, config_hash, n_patients, lambda_c, baseline_label, aft_mode,
  k_bases, scenario, max_FU, paste(selected_jmb2_models, collapse = ",")
))

output_dir <- Sys.getenv("SIM_OUTPUT_DIR", unset = "")
if (!nzchar(output_dir)) {
  output_dir <- file.path(repo_root, "Bunya", "ENZAMET", "output_jmbayes2")
}
output_dir <- path.expand(output_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
message("Batch output directory: ", output_dir)

out_file <- sprintf(
  "%s_k%d_n%d_cens%s_%s_scen%d_FU%d_batch%03d.rds",
  output_model_prefix, k_bases, n_patients, cens_label, baseline_label,
  scenario, max_FU, batch_id
)
out_file <- file.path(output_dir, out_file)

batch_output_complete <- function(path, expected_n) {
  if (!file.exists(path)) return(FALSE)
  existing <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(existing) || length(existing) != expected_n) return(FALSE)
  if (any(vapply(existing, is.null, logical(1)))) return(FALSE)
  has_marker <- vapply(existing, function(x) {
    is.list(x) &&
      "comparison_tidy" %in% names(x) &&
      "longitudinal_response_center" %in% names(x)
  }, logical(1))
  has_initial_survreg <- vapply(existing, function(x) {
    if (!is.list(x) || is.null(x$jmb2_results)) return(FALSE)
    aft_results <- x$jmb2_results[intersect(c("Weibull", "loglogistic"), names(x$jmb2_results))]
    if (length(aft_results) == 0) return(TRUE)
    all(vapply(aft_results, function(res) {
      is.list(res) && "initial_survreg_tidy" %in% names(res)
    }, logical(1)))
  }, logical(1))
  has_jmb2_random_sd <- vapply(existing, function(x) {
    if (!is.list(x) || is.null(x$jmb2_results)) return(FALSE)
    all(vapply(x$jmb2_results, function(res) {
      if (!is.list(res) || !identical(res$status, "ok")) return(FALSE)
      tidy <- res$tidy_summary
      is.data.frame(tidy) &&
        all(c("sd_1_long[1]", "sd_1_long[2]") %in% tidy$variable)
    }, logical(1)))
  }, logical(1))
  has_ok_status <- vapply(existing, function(x) {
    if (!is.list(x) || is.null(x$model_status) || nrow(x$model_status) == 0) {
      return(FALSE)
    }
    all(!is.na(x$model_status$status) & x$model_status$status == "ok")
  }, logical(1))
  all(has_marker & has_initial_survreg & has_jmb2_random_sd & has_ok_status)
}

if (batch_output_complete(out_file, length(trial_ids))) {
  message(sprintf("SKIP: complete output already exists for this batch: %s", out_file))
  quit(save = "no", status = 0)
}
if (file.exists(out_file)) {
  message(sprintf("RERUN: existing output is incomplete and will be overwritten: %s", out_file))
}

message(sprintf(
  "Running batch %d (config_id %d: n=%d, censoring=%s, baseline=%s, k=%d, scenario=%d, prefix=%s)",
  batch_id, config_id, n_patients, cens_label, baseline_label,
  k_bases, scenario, output_model_prefix
))

cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
workers <- floor(cpus / cpus_per_sim)
workers <- max(1L, min(workers, length(trial_ids)))
message(sprintf(
  "Parallel setup: cpus=%d, cpus_per_sim=%d, workers=%d, JMbayes2 cores/model=%d",
  cpus, cpus_per_sim, workers, jmb2_cores
))
future::plan(future::multisession, workers = workers)

safe_elapsed <- env_int("SIM_SAFE_ELAPSED", 18000)

results <- furrr::future_map(
  trial_ids,
  ~ tryCatch({
    t0 <- Sys.time()
    ans <- tryCatch({
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
      setTimeLimit(elapsed = safe_elapsed, cpu = Inf, transient = TRUE)
      run_simulation(
        .x, n_patients, lambda_c, aft_mode, max_FU, config_hash,
        k_bases, scenario, baseline_label,
        weibull_shape, weibull_scale, loglogistic_shape, loglogistic_scale
      )
    }, error = function(e) {
      if (grepl("reached elapsed time limit", conditionMessage(e))) NULL else stop(e)
    })
    message(sprintf(
      "Trial %d elapsed: %.1f min",
      .x, as.numeric(difftime(Sys.time(), t0, units = "mins"))
    ))
    ans
  }, error = function(e) {
    message(sprintf("Simulation %d failed or timed out: %s", .x, e$message))
    NULL
  }),
  .options = furrr::furrr_options(seed = TRUE)
)

saveRDS(results, file = out_file)
message("Saved batch output: ", out_file)
