#!/usr/bin/env Rscript

# =============================================================================
# SUMMARISE JMBAYES2 ENZAMET DM COMPARISON RUNS
# =============================================================================
#
# Reads batch RDS files produced by:
#   Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R
#
# The output is aligned with the proposed-method simulation summaries, with one
# row per method/parameter/config. The main intended comparison is longitudinal
# fixed effects across:
#   LMM, Cox-JM, Weibull-JM, loglogistic-JM.
#
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(here)
})

env_optional_int <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(NA_integer_)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) {
    stop(sprintf("Environment variable %s must be an integer if set.", name))
  }
  out
}

env_bool <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  value <- tolower(trimws(value))
  if (value %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  stop(sprintf("Environment variable %s must be boolean-like if set.", name))
}

analysis_output_dir <- Sys.getenv("SIM_OUTPUT_DIR", unset = "")
if (!nzchar(analysis_output_dir)) {
  analysis_output_dir <- here::here("Bunya", "ENZAMET", "output_20260602_ENZAMET-JMbayes2-DM-all40-1000rep-randSDsummary")
}
analysis_output_dir <- path.expand(analysis_output_dir)

base_output_dir <- dirname(analysis_output_dir)
summary_dir <- Sys.getenv("SIM_SUMMARY_DIR", unset = "")
if (!nzchar(summary_dir)) {
  summary_dir <- file.path(
    base_output_dir,
    paste0("summary_tables_", sub("^output_", "", basename(analysis_output_dir)))
  )
}
summary_dir <- path.expand(summary_dir)
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

figure_dir <- Sys.getenv("SIM_FIGURE_DIR", unset = "")
if (!nzchar(figure_dir)) {
  figure_dir <- file.path(
    base_output_dir,
    paste0("figures_", sub("^output_", "", basename(analysis_output_dir)))
  )
}
figure_dir <- path.expand(figure_dir)

variables_to_analyse <- trimws(strsplit(
  Sys.getenv(
    "SIM_JMB2_SUMMARY_VARIABLES",
    unset = "beta_long[1],beta_long[2],beta_long[3],sigma_long,sd_1_long[1],sd_1_long[2],alpha,beta_surv[1]"
  ),
  ",",
  fixed = TRUE
)[[1]])
variables_to_analyse <- variables_to_analyse[nzchar(variables_to_analyse)]

methods_to_analyse <- trimws(strsplit(
  Sys.getenv(
    "SIM_JMB2_SUMMARY_METHODS",
    unset = "LMM,Cox-JM,Weibull-JM,loglogistic-JM"
  ),
  ",",
  fixed = TRUE
)[[1]])
methods_to_analyse <- methods_to_analyse[nzchar(methods_to_analyse)]

distributions_to_analyse <- trimws(strsplit(
  Sys.getenv("SIM_JMB2_DISTS", unset = "Weibull,loglogistic"),
  ",",
  fixed = TRUE
)[[1]])
distributions_to_analyse <- distributions_to_analyse[nzchar(distributions_to_analyse)]

configured_batch_size <- env_optional_int("SIM_BATCH_SIZE")
configured_sims_per_config <- env_optional_int("SIMS_PER_CONFIG")
jmb2_models_for_status <- c("Cox", "Weibull", "loglogistic")

# The updated JMbayes2 driver already saves compact per-model tidy summaries.
# Reusing them avoids reparsing each stored raw JMbayes2 summary for every
# replicate. Set FALSE for legacy outputs if you need to rebuild from raw_summary.
use_saved_tidy_only <- env_bool("SIM_JMB2_USE_SAVED_TIDY_ONLY", TRUE)

message("Using JMbayes2 output directory: ", analysis_output_dir)
message("Writing summaries to: ", summary_dir)
message("Writing figures to: ", figure_dir)
message("Using saved tidy summaries only: ", use_saved_tidy_only)

parse_batch_number <- function(path) {
  m <- stringr::str_match(basename(path), "_batch([0-9]{3})\\.rds$")
  suppressWarnings(as.integer(m[, 2]))
}

parse_baseline_label <- function(label) {
  label <- as.character(label)

  if (tolower(label) == "weibull") {
    return(list(
      baseline_label = label,
      dist = "Weibull",
      dist_short = "W",
      baseline_shape = NA_real_,
      baseline_scale = NA_real_
    ))
  }

  if (tolower(label) == "loglogistic") {
    return(list(
      baseline_label = label,
      dist = "loglogistic",
      dist_short = "LL",
      baseline_shape = NA_real_,
      baseline_scale = NA_real_
    ))
  }

  w <- stringr::str_match(label, "^W([0-9]+)_([0-9]+)$")
  if (!is.na(w[1, 1])) {
    shape <- as.numeric(w[1, 2]) / 100
    scale <- as.numeric(w[1, 3]) / 10
    return(list(
      baseline_label = label,
      dist = "Weibull",
      dist_short = sprintf("W%.2f", shape),
      baseline_shape = shape,
      baseline_scale = scale
    ))
  }

  ll <- stringr::str_match(label, "^LL([0-9]+)_([0-9]+)$")
  if (!is.na(ll[1, 1])) {
    shape <- as.numeric(ll[1, 2]) / 100
    scale <- as.numeric(ll[1, 3]) / 10
    return(list(
      baseline_label = label,
      dist = "loglogistic",
      dist_short = sprintf("LL%.2f", shape),
      baseline_shape = shape,
      baseline_scale = scale
    ))
  }

  list(
    baseline_label = label,
    dist = label,
    dist_short = label,
    baseline_shape = NA_real_,
    baseline_scale = NA_real_
  )
}

parse_config <- function(config_str) {
  parts <- stringr::str_match(
    config_str,
    "^(JMbayes2(?:_[A-Za-z0-9]+)*)_k(\\d+)_n(\\d+)_cens([A-Za-z0-9]+)_([A-Za-z0-9_]+)_scen(\\d+)_?(FU\\d+)?$"
  )
  if (is.na(parts[1, 1])) return(NULL)
  baseline <- parse_baseline_label(parts[6])

  list(
    bf = parts[2],
    k_bases = as.numeric(parts[3]),
    n = as.numeric(parts[4]),
    cens = parts[5],
    baseline_label = baseline$baseline_label,
    dist = baseline$dist,
    dist_short = baseline$dist_short,
    baseline_shape = baseline$baseline_shape,
    baseline_scale = baseline$baseline_scale,
    scenario = as.numeric(parts[7]),
    fu = parts[8]
  )
}

truth_values_for_scenario <- function(scenario) {
  stopifnot(scenario %in% 1:5)
  c(
    "beta_long[1]" = 73,
    "beta_long[2]" = -0.04,
    "beta_long[3]" = c(0.00, 0.04, -0.04, 0.04, -0.04)[scenario],
    "sd_1_long[1]" = 15,
    "sd_1_long[2]" = 0.20,
    "sigma_long" = 12,
    "alpha" = 0.012,
    "beta_surv[1]" = c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario]
  )
}

jmb2_survival_sources <- c("Cox-JM", "Weibull-JM", "loglogistic-JM")
jmb2_survival_variables <- c("alpha", "beta_surv[1]")

truth_for_jmbayes2_summary <- function(variable, source, parsed, scenario) {
  scenario_truth <- truth_values_for_scenario(scenario)
  base_truth <- unname(scenario_truth[variable])
  if (is.na(base_truth)) return(NA_real_)

  is_jmb2_survival <- source %in% jmb2_survival_sources &&
    variable %in% jmb2_survival_variables

  if (!is_jmb2_survival) return(base_truth)

  # The data generators set alpha and beta_surv on the AFT log-time scale.
  # For Weibull-generated data, the corresponding PH-scale coefficient is
  # -shape * AFT coefficient. JMbayes2 survival summaries are compared on this
  # hazard-side scale. For log-logistic-generated data there is no single exact
  # PH truth, so survival-side truth is left undefined.
  if (tolower(parsed$dist) == "weibull" && is.finite(parsed$baseline_shape)) {
    return(-parsed$baseline_shape * base_truth)
  }

  if (tolower(parsed$dist) == "loglogistic") {
    return(NA_real_)
  }

  NA_real_
}

format_cens_short <- function(cens) {
  cens <- as.character(cens)
  dplyr::case_when(
    cens == "admin" ~ "admin",
    stringr::str_detect(cens, "^\\d+$") ~ paste0(cens, "%"),
    TRUE ~ cens
  )
}

shorten_config <- function(cfg) {
  parsed <- parse_config(cfg)
  if (is.null(parsed)) return(cfg)
  paste0(
    parsed$bf,
    " k=", parsed$k_bases,
    ", c=", format_cens_short(parsed$cens),
    ", d=", parsed$dist_short,
    ", S", parsed$scenario
  )
}

empty_tidy <- function() {
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
    rhat = numeric(),
    ess = numeric(),
    sim_id = integer()
  )
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

empty_initial_survreg_summary <- function() {
  tibble::tibble(
    source = character(),
    model = character(),
    variable = character(),
    truth = numeric(),
    mean = numeric(),
    median = numeric(),
    empirical_sd = numeric(),
    mean_post_sd = numeric(),
    sd = numeric(),
    q2.5_est = numeric(),
    q97.5_est = numeric(),
    bias = numeric(),
    rmse = numeric(),
    mean_ci_width = numeric(),
    coverage = numeric(),
    mcse_mean = numeric(),
    mcse_coverage = numeric(),
    n_sims = integer(),
    uncertainty_type = character(),
    interval_type = character()
  )
}

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

normalise_saved_jmb2_tidy <- function(x, model_name, sim_id) {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
    return(tibble::tibble())
  }
  if (!"variable" %in% names(x)) {
    return(tibble::tibble())
  }

  if (!"model" %in% names(x)) x$model <- model_name
  if (!"source" %in% names(x)) x$source <- paste0(model_name, "-JM")
  if (!"sim_id" %in% names(x)) x$sim_id <- sim_id
  if (!"rhat" %in% names(x)) x$rhat <- NA_real_
  if (!"ess" %in% names(x)) x$ess <- NA_real_

  x %>%
    dplyr::filter(!is.na(.data$variable)) %>%
    dplyr::mutate(
      model = dplyr::coalesce(as.character(.data$model), model_name),
      source = dplyr::coalesce(as.character(.data$source), paste0(model_name, "-JM")),
      sim_id = dplyr::coalesce(as.integer(.data$sim_id), as.integer(sim_id))
    )
}

tidy_jmb2_from_raw_summary <- function(sim, sim_id) {
  if (is.null(sim$jmb2_results) || length(sim$jmb2_results) == 0) {
    return(tibble::tibble())
  }

  purrr::imap_dfr(sim$jmb2_results, function(model_result, model_name) {
    saved_tidy <- normalise_saved_jmb2_tidy(
      model_result$tidy_summary,
      model_name = model_name,
      sim_id = sim_id
    )
    if (use_saved_tidy_only && nrow(saved_tidy) > 0) {
      return(saved_tidy)
    }

    raw_summary <- model_result$raw_summary
    if (is.null(raw_summary)) return(saved_tidy)

    raw <- purrr::imap_dfr(raw_summary, function(x, nm) {
      standardise_summary_table(x, nm, model_name, sim_id)
    })
    if (nrow(raw) == 0) return(saved_tidy)

    tidy <- raw %>%
      dplyr::mutate(
        variable = map_parameter_to_common(component, parameter),
        source = paste0(model_name, "-JM")
      ) %>%
      dplyr::filter(!is.na(variable))

    if (nrow(tidy) > 0 &&
        is.finite(sim$longitudinal_response_center) &&
        sim$longitudinal_response_center != 0) {
      beta0_row <- tidy$variable == "beta_long[1]"
      shift_cols <- intersect(c("mean", "median", "q2.5", "q97.5"), names(tidy))
      if (any(beta0_row) && length(shift_cols) > 0) {
        tidy[beta0_row, shift_cols] <- lapply(
          tidy[beta0_row, shift_cols, drop = FALSE],
          function(x) x + sim$longitudinal_response_center
        )
      }
    }

    dplyr::bind_rows(tidy, saved_tidy) %>%
      dplyr::distinct(.data$source, .data$variable, .data$parameter, .data$component, .keep_all = TRUE)
  })
}

repair_comparison_tidy <- function(sim, sim_id) {
  x <- sim$comparison_tidy
  if (is.null(x)) x <- empty_tidy()
  if (!"sim_id" %in% names(x)) x$sim_id <- sim_id

  lmm_rows <- x %>%
    dplyr::filter(.data$source == "LMM")
  jmb2_rows <- tidy_jmb2_from_raw_summary(sim, sim_id)

  repaired <- dplyr::bind_rows(lmm_rows, jmb2_rows)
  if (nrow(repaired) > 0) {
    repaired <- repaired %>%
      dplyr::distinct(.data$source, .data$variable, .data$parameter, .data$component, .keep_all = TRUE)
  }
  repaired
}

tidy_initial_survreg_from_sim <- function(sim, sim_id) {
  if (is.null(sim$jmb2_results) || length(sim$jmb2_results) == 0) {
    return(empty_initial_survreg_tidy())
  }

  out <- purrr::imap_dfr(sim$jmb2_results, function(model_result, model_name) {
    x <- model_result$initial_survreg_tidy
    if (is.null(x) || nrow(x) == 0) return(NULL)
    if (!"sim_id" %in% names(x)) x$sim_id <- sim_id
    x
  })

  if (nrow(out) == 0) return(empty_initial_survreg_tidy())
  out
}

truth_for_initial_survreg <- function(variable, model, parsed, scenario) {
  model_matches_data <- tolower(model) == tolower(parsed$dist)
  scenario_truth <- truth_values_for_scenario(scenario)

  # Initial survreg rows are survival-only AFT fits. For beta_surv[1], keep the
  # raw data-generating AFT coefficient as a reference even when the fitted
  # survival family is misspecified. Baseline shape/scale truths are only exact
  # for the matching generating family.
  dplyr::case_when(
    variable == "baseline_shape" & model_matches_data ~ parsed$baseline_shape,
    variable == "baseline_scale" & model_matches_data ~ parsed$baseline_scale,
    variable == "beta_surv[1]" ~ unname(scenario_truth["beta_surv[1]"]),
    TRUE ~ NA_real_
  )
}

summarise_initial_survreg <- function(per_sim, parsed) {
  if (nrow(per_sim) == 0) return(empty_initial_survreg_summary())

  per_sim %>%
    dplyr::mutate(
      point_estimate = .data$mean,
      uncertainty = .data$sd,
      truth = purrr::pmap_dbl(
        list(.data$variable, .data$model),
        ~ truth_for_initial_survreg(..1, ..2, parsed, parsed$scenario)
      ),
      error = point_estimate - truth,
      covered = dplyr::if_else(
        !is.na(truth),
        q2.5 <= truth & truth <= q97.5,
        NA
      ),
      ci_width = q97.5 - q2.5
    ) %>%
    dplyr::group_by(source, model, variable) %>%
    dplyr::summarise(
      truth = dplyr::first(stats::na.omit(truth), default = NA_real_),
      n_sims = sum(!is.na(point_estimate)),
      mean = mean(.data$point_estimate, na.rm = TRUE),
      median = median(.data$point_estimate, na.rm = TRUE),
      empirical_sd = sd(.data$point_estimate, na.rm = TRUE),
      mean_post_sd = mean(.data$uncertainty, na.rm = TRUE),
      sd = mean_post_sd,
      q2.5_est = stats::quantile(.data$point_estimate, 0.025, na.rm = TRUE),
      q97.5_est = stats::quantile(.data$point_estimate, 0.975, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      coverage = mean(covered, na.rm = TRUE),
      mcse_mean = empirical_sd / sqrt(n_sims),
      mcse_coverage = sqrt(coverage * (1 - coverage) / n_sims),
      uncertainty_type = "frequentist_se",
      interval_type = "confidence",
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ replace(.x, is.nan(.x), NA_real_))
    ) %>%
    dplyr::select(
      source, model, variable, truth, mean, median, empirical_sd,
      mean_post_sd, sd, q2.5_est, q97.5_est, bias, rmse,
      mean_ci_width, coverage, mcse_mean, mcse_coverage, n_sims,
      uncertainty_type, interval_type
    ) %>%
    dplyr::arrange(source, variable)
}

process_jmbayes2_config <- function(config, output_dir) {
  parsed <- parse_config(config)
  if (is.null(parsed)) return(NULL)

  pattern <- paste0(
    config,
    "(?:_.*)?_batch\\d{3}\\.rds$"
  )
  rds_files <- list.files(output_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(rds_files) == 0) {
    message("No matching files found for: ", pattern)
    return(NULL)
  }
  rds_files <- sort(rds_files)

  batch_contents <- purrr::map(rds_files, readRDS)
  batch_sizes <- vapply(batch_contents, length, integer(1))
  nominal_batch_size <- if (!is.na(configured_batch_size)) {
    configured_batch_size
  } else if (length(batch_sizes) == 0 || max(batch_sizes, na.rm = TRUE) == 0) {
    NA_integer_
  } else {
    max(batch_sizes, na.rm = TRUE)
  }

  observed_batch_numbers <- parse_batch_number(rds_files)
  observed_batch_numbers <- observed_batch_numbers[!is.na(observed_batch_numbers)]
  expected_batches_per_config <- if (!is.na(configured_sims_per_config) &&
      !is.na(nominal_batch_size)) {
    if (configured_sims_per_config %% nominal_batch_size != 0) {
      stop(sprintf(
        "SIMS_PER_CONFIG (%d) must be divisible by SIM_BATCH_SIZE/nominal batch size (%d).",
        configured_sims_per_config, nominal_batch_size
      ))
    }
    configured_sims_per_config %/% nominal_batch_size
  } else if (!is.na(observed_max_batch_number) && observed_max_batch_number > 0) {
    observed_max_batch_number
  } else if (length(observed_batch_numbers) > 0) {
    max(observed_batch_numbers, na.rm = TRUE)
  } else {
    NA_integer_
  }

  batch_records <- purrr::imap_dfr(batch_contents, function(batch_sims, batch_idx) {
    batch_file <- rds_files[[batch_idx]]
    batch_number <- parse_batch_number(batch_file)
    if (length(batch_sims) == 0) {
      return(tibble::tibble(
        batch_file = basename(batch_file),
        batch_number = batch_number,
        position_in_batch = NA_integer_,
        replicate_in_config = NA_integer_,
        is_null = TRUE,
        sim = list(NULL)
      ))
    }

    purrr::map_dfr(seq_along(batch_sims), function(pos) {
      replicate_in_config <- if (is.na(batch_number) || is.na(nominal_batch_size)) {
        NA_integer_
      } else {
        (batch_number - 1L) * nominal_batch_size + pos
      }
      tibble::tibble(
        batch_file = basename(batch_file),
        batch_number = batch_number,
        position_in_batch = pos,
        replicate_in_config = replicate_in_config,
        is_null = is.null(batch_sims[[pos]]),
        sim = list(batch_sims[[pos]])
      )
    })
  })

  missing_batch_numbers <- integer()
  if (!is.na(expected_batches_per_config)) {
    missing_batch_numbers <- setdiff(seq_len(expected_batches_per_config), observed_batch_numbers)
  }
  missing_batch_records <- purrr::map_dfr(missing_batch_numbers, function(batch_number) {
    if (is.na(nominal_batch_size)) return(tibble::tibble())
    purrr::map_dfr(seq_len(nominal_batch_size), function(pos) {
      tibble::tibble(
        batch_file = sprintf("%s_batch%03d.rds", config, batch_number),
        batch_number = batch_number,
        position_in_batch = pos,
        replicate_in_config = (batch_number - 1L) * nominal_batch_size + pos,
        is_null = TRUE,
        failure_type = "missing_batch_file",
        sim = list(NULL)
      )
    })
  })

  if (nrow(missing_batch_records) > 0) {
    if (!"failure_type" %in% names(batch_records)) {
      batch_records$failure_type <- NA_character_
    }
    batch_records <- dplyr::bind_rows(batch_records, missing_batch_records) %>%
      dplyr::arrange(.data$batch_number, .data$position_in_batch)
  } else if (!"failure_type" %in% names(batch_records)) {
    batch_records$failure_type <- NA_character_
  }

  sims <- batch_records %>%
    filter(!is_null) %>%
    pull(sim)

  failed_replicates <- batch_records %>%
    filter(is_null) %>%
    transmute(
      batch_file,
      batch_number,
      position_in_batch,
      replicate_in_config,
      failure_type = dplyr::coalesce(.data$failure_type, "null_saved_replicate")
    )

  if (length(sims) == 0) {
    return(list(
      config = config,
      parsed = parsed,
      per_sim = empty_tidy(),
      summary = tibble::tibble(),
      initial_survreg_per_sim = empty_initial_survreg_tidy(),
      initial_survreg_summary = empty_initial_survreg_summary(),
      runtime_per_sim = tibble::tibble(),
      runtime_summary = tibble::tibble(),
      status_per_sim = tibble::tibble(),
      diagnostic_summary = tibble::tibble(
        n_batch_files = length(rds_files),
        n_raw_sims = nrow(batch_records),
        n_failed_sims = nrow(failed_replicates),
        n_missing_batch_files = length(missing_batch_numbers),
        n_missing_batch_replicates = sum(failed_replicates$failure_type == "missing_batch_file"),
        n_null_saved_replicates = sum(failed_replicates$failure_type == "null_saved_replicate"),
        n_valid_sims = 0L,
        mean_censor_prop = NA_real_,
        mean_object_size_mb = NA_real_,
        max_object_size_mb = NA_real_
      ),
      failed_replicates = failed_replicates
    ))
  }

  per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
    x <- repair_comparison_tidy(sim, idx)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x %>%
      dplyr::mutate(
        replicate_index = idx,
        unique_seed = sim$unique_seed,
        censor_prop = sim$censor_prop,
        longitudinal_response_center = sim$longitudinal_response_center
      )
  }) %>%
    filter(
      source %in% methods_to_analyse,
      variable %in% variables_to_analyse
    ) %>%
    mutate(
      uncertainty_type = dplyr::case_when(
        source == "LMM" & !is.na(sd) ~ "frequentist_se",
        source == "LMM" ~ "none",
        !is.na(sd) ~ "posterior_sd",
        TRUE ~ "none"
      ),
      interval_type = dplyr::case_when(
        source == "LMM" & !is.na(q2.5) & !is.na(q97.5) ~ "confidence",
        source == "LMM" ~ "none",
        !is.na(q2.5) & !is.na(q97.5) ~ "credible",
        TRUE ~ "none"
      )
    )

  initial_survreg_per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
    x <- tidy_initial_survreg_from_sim(sim, idx)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x %>%
      dplyr::mutate(
        replicate_index = idx,
        unique_seed = sim$unique_seed,
        censor_prop = sim$censor_prop
      )
  })
  if (nrow(initial_survreg_per_sim) == 0) {
    initial_survreg_per_sim <- empty_initial_survreg_tidy()
  }

  initial_survreg_summary <- summarise_initial_survreg(initial_survreg_per_sim, parsed)

  summary <- per_sim %>%
    dplyr::mutate(
      posterior_mean = .data$mean,
      posterior_sd = .data$sd,
      truth = purrr::pmap_dbl(
        list(.data$variable, .data$source),
        ~ truth_for_jmbayes2_summary(..1, ..2, parsed, parsed$scenario)
      ),
      error = posterior_mean - truth,
      covered = dplyr::if_else(
        !is.na(truth),
        q2.5 <= truth & truth <= q97.5,
        NA
      ),
      ci_width = q97.5 - q2.5
    ) %>%
    dplyr::group_by(source, variable) %>%
    dplyr::summarise(
      truth = dplyr::first(stats::na.omit(truth), default = NA_real_),
      n_sims = sum(!is.na(posterior_mean)),
      mean = mean(.data$posterior_mean, na.rm = TRUE),
      median = median(.data$posterior_mean, na.rm = TRUE),
      empirical_sd = sd(.data$posterior_mean, na.rm = TRUE),
      mean_post_sd = mean(.data$posterior_sd, na.rm = TRUE),
      sd = mean_post_sd,
      q2.5_est = stats::quantile(.data$posterior_mean, 0.025, na.rm = TRUE),
      q97.5_est = stats::quantile(.data$posterior_mean, 0.975, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      coverage = mean(covered, na.rm = TRUE),
      mcse_mean = empirical_sd / sqrt(n_sims),
      mcse_coverage = sqrt(coverage * (1 - coverage) / n_sims),
      uncertainty_type = dplyr::first(.data$uncertainty_type),
      interval_type = dplyr::first(.data$interval_type),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ replace(.x, is.nan(.x), NA_real_))
    ) %>%
    dplyr::select(
      source, variable, truth, mean, median, empirical_sd, mean_post_sd, sd,
      q2.5_est, q97.5_est, bias, rmse, mean_ci_width, coverage,
      mcse_mean, mcse_coverage, n_sims, uncertainty_type, interval_type
    ) %>%
    dplyr::arrange(source, variable)

  runtime_per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
    rt <- sim$runtimes
    if (is.null(rt) || nrow(rt) == 0) return(NULL)
    rt %>%
      mutate(
        replicate_index = idx,
        unique_seed = sim$unique_seed,
        censor_prop = sim$censor_prop
      )
  })

  runtime_summary <- runtime_per_sim %>%
    group_by(model) %>%
    summarise(
      n_sims_with_runtime = sum(!is.na(total_model_time)),
      mean_total_time_min = mean(total_model_time / 60, na.rm = TRUE),
      median_total_time_min = median(total_model_time / 60, na.rm = TRUE),
      sd_total_time_min = sd(total_model_time / 60, na.rm = TRUE),
      min_total_time_min = min(total_model_time / 60, na.rm = TRUE),
      max_total_time_min = max(total_model_time / 60, na.rm = TRUE),
      q90_total_time_min = stats::quantile(total_model_time / 60, 0.90, na.rm = TRUE),
      q95_total_time_min = stats::quantile(total_model_time / 60, 0.95, na.rm = TRUE),
      mean_object_size_mb = mean(object_size_mb, na.rm = TRUE),
      max_object_size_mb = max(object_size_mb, na.rm = TRUE),
      .groups = "drop"
    )

  status_per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
    st <- sim$model_status
    if (is.null(st) || nrow(st) == 0) return(NULL)
    st %>%
      mutate(
        replicate_index = idx,
        unique_seed = sim$unique_seed,
        censor_prop = sim$censor_prop
      )
  })

  model_failure_counts <- status_per_sim %>%
    dplyr::filter(is.na(.data$status) | .data$status != "ok") %>%
    dplyr::count(.data$model, name = "n_failed")
  model_failure_msg <- if (nrow(model_failure_counts) == 0) {
    "none"
  } else {
    paste(
      paste0(model_failure_counts$model, "=", model_failure_counts$n_failed),
      collapse = ", "
    )
  }

  diagnostic_summary <- tibble::tibble(
    n_batch_files = length(rds_files),
    n_raw_sims = nrow(batch_records),
    n_failed_sims = nrow(failed_replicates),
    n_missing_batch_files = length(missing_batch_numbers),
    n_missing_batch_replicates = sum(failed_replicates$failure_type == "missing_batch_file"),
    n_null_saved_replicates = sum(failed_replicates$failure_type == "null_saved_replicate"),
    n_valid_sims = length(sims),
    mean_censor_prop = mean(vapply(sims, function(x) x$censor_prop, numeric(1)), na.rm = TRUE),
    mean_object_size_mb = mean(status_per_sim$object_size_mb, na.rm = TRUE),
    max_object_size_mb = max(status_per_sim$object_size_mb, na.rm = TRUE)
  )

  message(
    "Processed ", config,
    ": valid N=", length(sims),
    ", failed/null replicates=", nrow(failed_replicates),
    ", model failures=", model_failure_msg
  )

  list(
    config = config,
    parsed = parsed,
    per_sim = per_sim,
    summary = summary,
    initial_survreg_per_sim = initial_survreg_per_sim,
    initial_survreg_summary = initial_survreg_summary,
    runtime_per_sim = runtime_per_sim,
    runtime_summary = runtime_summary,
    status_per_sim = status_per_sim,
    diagnostic_summary = diagnostic_summary,
    failed_replicates = failed_replicates
  )
}

stopifnot(dir.exists(analysis_output_dir))
filenames <- list.files(analysis_output_dir, pattern = "\\.[rR][dD][sS]$")
observed_batch_numbers_all <- parse_batch_number(filenames)
observed_batch_numbers_all <- observed_batch_numbers_all[!is.na(observed_batch_numbers_all)]
observed_max_batch_number <- if (length(observed_batch_numbers_all) > 0) {
  max(observed_batch_numbers_all)
} else {
  NA_integer_
}
configs <- substr(
  tools::file_path_sans_ext(filenames),
  1,
  nchar(tools::file_path_sans_ext(filenames)) - 9
) %>%
  unique()

configs <- configs[
  purrr::map_lgl(configs, function(cfg) {
    parsed <- parse_config(cfg)
    !is.null(parsed) && parsed$dist %in% distributions_to_analyse
  })
]

if (length(configs) == 0) {
  stop("No JMbayes2 config batch files found under: ", analysis_output_dir)
}

print(configs)

results_by_config <- purrr::compact(
  setNames(configs, configs) %>%
    purrr::map(process_jmbayes2_config, output_dir = analysis_output_dir)
)

config_info <- tibble::tibble(config = names(results_by_config)) %>%
  mutate(parsed = purrr::map(config, parse_config)) %>%
  tidyr::unnest_wider(parsed)

analysis_dist_title <- paste(sort(unique(config_info$dist)), collapse = "_")
analysis_dist_file <- tolower(gsub("[^A-Za-z0-9]+", "_", analysis_dist_title))
if (!nzchar(analysis_dist_file)) analysis_dist_file <- "all"

all_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$per_sim, config = .y, cfg_short = shorten_config(.y))
)

all_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$summary, config = .y, cfg_short = shorten_config(.y))
)

all_initial_survreg_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$initial_survreg_per_sim, config = .y, cfg_short = shorten_config(.y))
)

all_initial_survreg_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$initial_survreg_summary, config = .y, cfg_short = shorten_config(.y))
)

all_runtime_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$runtime_per_sim, config = .y, cfg_short = shorten_config(.y))
)

all_runtime_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$runtime_summary, config = .y, cfg_short = shorten_config(.y))
)

all_status_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$status_per_sim, config = .y, cfg_short = shorten_config(.y))
)

all_diagnostic_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$diagnostic_summary, config = .y, cfg_short = shorten_config(.y))
)

all_failed_replicates <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$failed_replicates, config = .y, cfg_short = shorten_config(.y), .before = 1)
)

config_status_grid <- config_info %>%
  dplyr::mutate(cfg_short = purrr::map_chr(.data$config, shorten_config)) %>%
  dplyr::select(config, cfg_short) %>%
  tidyr::crossing(model = jmb2_models_for_status)

status_source <- all_status_per_sim
if (!"retry_used" %in% names(status_source)) {
  status_source$retry_used <- FALSE
}

status_counts <- status_source %>%
  dplyr::mutate(
    model = as.character(.data$model),
    status = dplyr::coalesce(as.character(.data$status), "missing_status"),
    retry_used = !is.na(.data$retry_used) & .data$retry_used
  ) %>%
  dplyr::group_by(.data$config, .data$cfg_short, .data$model) %>%
  dplyr::summarise(
    n_status_rows = dplyr::n(),
    n_model_ok = sum(.data$status == "ok"),
    n_model_failed = sum(.data$status != "ok"),
    n_retry_used = sum(.data$retry_used, na.rm = TRUE),
    status_breakdown = paste(
      paste0(names(table(.data$status)), ":", as.integer(table(.data$status))),
      collapse = "; "
    ),
    .groups = "drop"
  )

failed_null_counts <- all_failed_replicates %>%
  dplyr::group_by(.data$config, .data$cfg_short) %>%
  dplyr::summarise(
    n_failed_null_replicates = dplyr::n(),
    failed_null_breakdown = paste(
      paste0(names(table(.data$failure_type)), ":", as.integer(table(.data$failure_type))),
      collapse = "; "
    ),
    .groups = "drop"
  )

diagnostic_counts <- all_diagnostic_summary %>%
  dplyr::select(
    config, cfg_short, n_batch_files, n_raw_sims, n_valid_sims,
    n_failed_sims, n_missing_batch_files, n_missing_batch_replicates,
    n_null_saved_replicates
  )

model_failure_count_cols <- c(
  "n_batch_files", "n_raw_sims", "n_valid_sims", "n_failed_sims",
  "n_missing_batch_files", "n_missing_batch_replicates",
  "n_null_saved_replicates", "n_status_rows", "n_model_ok",
  "n_model_failed", "n_retry_used", "n_failed_null_replicates"
)

all_model_failure_summary <- config_status_grid %>%
  dplyr::left_join(diagnostic_counts, by = c("config", "cfg_short")) %>%
  dplyr::left_join(status_counts, by = c("config", "cfg_short", "model")) %>%
  dplyr::left_join(failed_null_counts, by = c("config", "cfg_short")) %>%
  dplyr::mutate(
    dplyr::across(dplyr::any_of(model_failure_count_cols), ~ dplyr::coalesce(as.integer(.x), 0L)),
    n_missing_model_status = pmax(.data$n_valid_sims - .data$n_status_rows, 0L),
    n_expected_replicates = .data$n_valid_sims + .data$n_failed_null_replicates,
    n_failed_or_null = .data$n_model_failed +
      .data$n_missing_model_status +
      .data$n_failed_null_replicates,
    failure_or_null_rate = dplyr::if_else(
      .data$n_expected_replicates > 0L,
      .data$n_failed_or_null / .data$n_expected_replicates,
      NA_real_
    ),
    status_breakdown = dplyr::coalesce(.data$status_breakdown, "none"),
    failed_null_breakdown = dplyr::coalesce(.data$failed_null_breakdown, "none")
  ) %>%
  dplyr::transmute(
    .data$config,
    .data$cfg_short,
    .data$model,
    n_expected_replicates = .data$n_expected_replicates,
    n_valid_replicates = .data$n_valid_sims,
    n_failed_null_replicates = .data$n_failed_null_replicates,
    n_missing_batch_replicates = .data$n_missing_batch_replicates,
    n_null_saved_replicates = .data$n_null_saved_replicates,
    n_status_rows = .data$n_status_rows,
    n_model_ok = .data$n_model_ok,
    n_model_failed = .data$n_model_failed,
    n_missing_model_status = .data$n_missing_model_status,
    n_failed_or_null = .data$n_failed_or_null,
    failure_or_null_rate = .data$failure_or_null_rate,
    n_retry_used = .data$n_retry_used,
    status_breakdown = .data$status_breakdown,
    failed_null_breakdown = .data$failed_null_breakdown,
    n_batch_files = .data$n_batch_files,
    n_missing_batch_files = .data$n_missing_batch_files
  ) %>%
  dplyr::arrange(.data$config, .data$model)

write_outputs <- function(x, stem) {
  csv_path <- file.path(summary_dir, paste0(stem, "_", analysis_dist_file, ".csv"))
  rds_path <- file.path(summary_dir, paste0(stem, "_", analysis_dist_file, ".rds"))
  readr::write_csv(x, csv_path)
  saveRDS(x, rds_path)
  message("Wrote: ", csv_path)
  message("Wrote: ", rds_path)
}

generate_jmbayes2_figures <- function(all_per_sim, config_info, fig_dir) {
  key_vars <- intersect(
    c("alpha", "beta_long[1]", "beta_long[2]", "beta_long[3]", "beta_surv[1]"),
    variables_to_analyse
  )
  if (nrow(all_per_sim) == 0 || length(key_vars) == 0) {
    message("Skipping JMbayes2 figures because no per-simulation plotting rows are available.")
    return(invisible(NULL))
  }

  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  perf_dir <- file.path(fig_dir, "performance")
  if (!dir.exists(perf_dir)) dir.create(perf_dir, recursive = TRUE)

  method_order <- c("LMM", "Cox-JM", "Weibull-JM", "loglogistic-JM")
  param_labels <- c(
    "alpha" = "alpha",
    "beta_long[1]" = "beta[0]",
    "beta_long[2]" = "beta[1]",
    "beta_long[3]" = "beta[2]",
    "beta_surv[1]" = "gamma"
  )

  config_meta <- config_info %>%
    dplyr::select(
      "config", "dist", "dist_short", "baseline_shape",
      "scenario", "cens", "k_bases"
    )

  plot_data <- all_per_sim %>%
    dplyr::filter(.data$variable %in% key_vars) %>%
    dplyr::left_join(config_meta, by = "config") %>%
    dplyr::mutate(
      truth = purrr::pmap_dbl(
        list(.data$variable, .data$source, .data$dist, .data$baseline_shape, .data$scenario),
        function(variable, source, dist, baseline_shape, scenario) {
          truth_for_jmbayes2_summary(
            variable,
            source,
            list(dist = dist, baseline_shape = baseline_shape),
            scenario
          )
        }
      ),
      Method = factor(.data$source, levels = method_order),
      Parameter = factor(.data$variable, levels = key_vars, labels = unname(param_labels[key_vars])),
      Baseline = factor(.data$dist_short),
      Scenario = factor(.data$scenario),
      Cens = factor(format_cens_short(.data$cens)),
      k = factor(.data$k_bases)
    ) %>%
    dplyr::filter(!is.na(.data$Method), !is.na(.data$Parameter))

  if (nrow(plot_data) == 0) {
    message("Skipping JMbayes2 figures because plotting data are empty after filtering.")
    return(invisible(NULL))
  }

  facet_labeller <- ggplot2::labeller(Parameter = ggplot2::label_parsed)

  reference_lines <- plot_data %>%
    dplyr::filter(!is.na(.data$truth)) %>%
    dplyr::distinct(.data$Parameter, .data$truth)

  fig_main <- plot_data %>%
    dplyr::filter(is.finite(.data$mean)) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$cfg_short, y = .data$mean, fill = .data$Method)) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.size = 0.25, position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_hline(
      data = reference_lines,
      ggplot2::aes(yintercept = .data$truth),
      inherit.aes = FALSE,
      color = "red",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::facet_wrap(~ Parameter, scales = "free_y", labeller = facet_labeller) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      x = "Configuration",
      y = "Estimate",
      fill = "Method",
      title = paste0("JMbayes2/LMM estimates - ENZAMET (", analysis_dist_title, ", scenarios 1-5)")
    )

  main_path <- file.path(fig_dir, paste0("fig_ENZAMET_JMbayes2_", analysis_dist_file, ".pdf"))
  ggplot2::ggsave(main_path, fig_main, width = 18, height = 10, units = "in")
  message("Wrote: ", main_path)

  alpha_data <- plot_data %>%
    dplyr::filter(
      .data$variable == "alpha",
      .data$source %in% jmb2_survival_sources,
      is.finite(.data$mean)
    )

  if (nrow(alpha_data) > 0) {
    alpha_truth_lines <- alpha_data %>%
      dplyr::filter(!is.na(.data$truth)) %>%
      dplyr::distinct(.data$Cens, .data$Baseline, .data$truth)

    p_alpha_est <- ggplot2::ggplot(
      alpha_data,
      ggplot2::aes(x = .data$Scenario, y = .data$mean, fill = .data$Method)
    ) +
      ggplot2::geom_boxplot(outlier.size = 0.25, position = ggplot2::position_dodge(width = 0.8)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey40", linetype = "dotted") +
      ggplot2::geom_hline(
        data = alpha_truth_lines,
        ggplot2::aes(yintercept = .data$truth),
        inherit.aes = FALSE,
        color = "red",
        linetype = "dashed",
        linewidth = 0.5
      ) +
      ggplot2::facet_grid(Cens ~ Baseline, scales = "free_y") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        x = "Scenario",
        y = "Alpha estimate",
        fill = "Method",
        title = paste0("JMbayes2 alpha estimates across simulations (", analysis_dist_title, ", scenarios 1-5)"),
        subtitle = "Red dashed lines show converted PH-scale truth where exact; log-logistic generated data have no exact PH truth."
      )

    alpha_est_path <- file.path(perf_dir, paste0("alpha_estimates_jmbayes2_", analysis_dist_file, "_scen1_5.png"))
    ggplot2::ggsave(alpha_est_path, p_alpha_est, width = 14, height = 8, units = "in", dpi = 300)
    message("Wrote: ", alpha_est_path)

    alpha_bias_data <- alpha_data %>%
      dplyr::filter(!is.na(.data$truth)) %>%
      dplyr::mutate(alpha_bias = .data$mean - .data$truth)

    if (nrow(alpha_bias_data) > 0) {
      p_alpha_bias <- ggplot2::ggplot(
        alpha_bias_data,
        ggplot2::aes(x = .data$Scenario, y = .data$alpha_bias, fill = .data$Method)
      ) +
        ggplot2::geom_boxplot(outlier.size = 0.25, position = ggplot2::position_dodge(width = 0.8)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::facet_grid(Cens ~ Baseline, scales = "free_y") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(
          x = "Scenario",
          y = "Alpha bias (estimate - converted truth)",
          fill = "Method",
          title = paste0("JMbayes2 alpha bias across simulations (", analysis_dist_title, ", scenarios 1-5)"),
          subtitle = "Bias is shown only where an exact PH-scale truth is available."
        )

      alpha_bias_path <- file.path(perf_dir, paste0("alpha_bias_jmbayes2_", analysis_dist_file, "_scen1_5.png"))
      ggplot2::ggsave(alpha_bias_path, p_alpha_bias, width = 14, height = 8, units = "in", dpi = 300)
      message("Wrote: ", alpha_bias_path)
    } else {
      message("Skipping alpha-bias figure because no exact PH-scale alpha truth is available.")
    }
  }

  perf_sim <- plot_data %>%
    dplyr::filter(!is.na(.data$truth), is.finite(.data$mean)) %>%
    dplyr::mutate(
      bias = .data$mean - .data$truth,
      abs_err = abs(.data$bias),
      sq_error = .data$bias^2,
      covered = dplyr::if_else(
        !is.na(.data$q2.5) & !is.na(.data$q97.5),
        as.integer(.data$q2.5 <= .data$truth & .data$truth <= .data$q97.5),
        NA_integer_
      )
    )

  if (nrow(perf_sim) == 0) {
    message("Skipping performance figures because all plotting truths are undefined.")
    return(invisible(NULL))
  }

  p_bias <- ggplot2::ggplot(perf_sim, ggplot2::aes(x = .data$Scenario, y = .data$bias, fill = .data$Method)) +
    ggplot2::geom_boxplot(outlier.size = 0.25, position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(Parameter ~ Cens + Baseline, scales = "free_y", labeller = facet_labeller) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "Scenario",
      y = "Bias (estimate - truth)",
      fill = "Method",
      title = paste0("Bias distributions across simulations (", analysis_dist_title, ", scenarios 1-5)")
    )

  p_abserr <- ggplot2::ggplot(perf_sim, ggplot2::aes(x = .data$Scenario, y = .data$abs_err, fill = .data$Method)) +
    ggplot2::geom_boxplot(outlier.size = 0.25, position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::facet_grid(Parameter ~ Cens + Baseline, scales = "free_y", labeller = facet_labeller) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "Scenario",
      y = "Absolute error |estimate - truth|",
      fill = "Method",
      title = paste0("Absolute error distributions across simulations (", analysis_dist_title, ", scenarios 1-5)")
    )

  rmse_sum <- perf_sim %>%
    dplyr::group_by(.data$Parameter, .data$Method, .data$Baseline, .data$Scenario, .data$Cens, .data$k) %>%
    dplyr::summarise(rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)), .groups = "drop")

  p_rmse <- ggplot2::ggplot(rmse_sum, ggplot2::aes(x = .data$Scenario, y = .data$rmse, color = .data$Method, group = .data$Method)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.35)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35), alpha = 0.6) +
    ggplot2::facet_grid(Parameter ~ Cens + Baseline, scales = "free_y", labeller = facet_labeller) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "Scenario",
      y = "RMSE across simulations",
      color = "Method",
      title = paste0("RMSE by scenario (", analysis_dist_title, ", scenarios 1-5)")
    )

  cov_sum <- perf_sim %>%
    dplyr::filter(!is.na(.data$covered)) %>%
    dplyr::group_by(.data$Parameter, .data$Method, .data$Baseline, .data$Scenario, .data$Cens, .data$k) %>%
    dplyr::summarise(
      coverage = mean(.data$covered, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      se = sqrt(.data$coverage * (1 - .data$coverage) / .data$n),
      lo = pmax(0, .data$coverage - 1.96 * .data$se),
      hi = pmin(1, .data$coverage + 1.96 * .data$se)
    )

  p_cov <- ggplot2::ggplot(cov_sum, ggplot2::aes(x = .data$Scenario, y = .data$coverage, color = .data$Method, group = .data$Method)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.35)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      width = 0.15,
      position = ggplot2::position_dodge(width = 0.35)
    ) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed") +
    ggplot2::facet_grid(Parameter ~ Cens + Baseline, labeller = facet_labeller) +
    ggplot2::coord_cartesian(ylim = c(0.85, 1.00)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = "Scenario",
      y = "Coverage (mean +/- 1.96 SE)",
      color = "Method",
      title = paste0("Coverage across simulations (", analysis_dist_title, ", scenarios 1-5)")
    )

  perf_paths <- c(
    bias = file.path(perf_dir, paste0("bias_distributions_", analysis_dist_file, "_scen1_5.png")),
    abs_error = file.path(perf_dir, paste0("abs_error_distributions_", analysis_dist_file, "_scen1_5.png")),
    rmse = file.path(perf_dir, paste0("rmse_summary_", analysis_dist_file, "_scen1_5.png")),
    coverage = file.path(perf_dir, paste0("coverage_", analysis_dist_file, "_scen1_5.png"))
  )

  ggplot2::ggsave(perf_paths["bias"], p_bias, width = 18, height = 10, units = "in", dpi = 300)
  ggplot2::ggsave(perf_paths["abs_error"], p_abserr, width = 18, height = 10, units = "in", dpi = 300)
  ggplot2::ggsave(perf_paths["rmse"], p_rmse, width = 18, height = 10, units = "in", dpi = 300)
  if (nrow(cov_sum) > 0) {
    ggplot2::ggsave(perf_paths["coverage"], p_cov, width = 18, height = 10, units = "in", dpi = 300)
  }
  message("Wrote JMbayes2 performance figures to: ", perf_dir)
  invisible(NULL)
}

write_outputs(all_per_sim, "jmbayes2_comparison_per_sim")
write_outputs(all_summary, "jmbayes2_comparison_summary")
if (nrow(all_initial_survreg_per_sim) == 0) {
  message(
    "No initial survreg rows found. These are available only for outputs ",
    "generated by the updated ENZAMET-JMbayes2-DM.R driver."
  )
}
write_outputs(all_initial_survreg_per_sim, "jmbayes2_initial_survreg_per_sim")
write_outputs(all_initial_survreg_summary, "jmbayes2_initial_survreg_summary")
write_outputs(all_runtime_per_sim, "jmbayes2_runtime_per_sim")
write_outputs(all_runtime_summary, "jmbayes2_runtime_summary")
write_outputs(all_status_per_sim, "jmbayes2_status_per_sim")
write_outputs(all_diagnostic_summary, "jmbayes2_diagnostic_summary")
write_outputs(all_failed_replicates, "jmbayes2_failed_null_replicates")
write_outputs(all_model_failure_summary, "jmbayes2_model_failure_summary")

longitudinal_fixed <- all_summary %>%
  filter(variable %in% c("beta_long[1]", "beta_long[2]", "beta_long[3]"))
write_outputs(longitudinal_fixed, "jmbayes2_longitudinal_fixed_summary")

generate_jmbayes2_figures(all_per_sim, config_info, figure_dir)

message("JMbayes2 ENZAMET summary complete.")
