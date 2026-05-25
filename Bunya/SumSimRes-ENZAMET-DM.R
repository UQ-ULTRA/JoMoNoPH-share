# =============================================================================
# SIMULATION RESULTS SUMMARY AND VISUALIZATION - ENZAMET
# =============================================================================
# 
# OVERVIEW:
# This script aggregates and visualizes Bayesian joint modelling simulation 
# results for ENZAMET scenarios (7 different data generating mechanisms).
# Compares Stan posterior estimates across scenarios.
#
# WORKFLOW:
# 1. File preprocessing and configuration parsing (including scenario)
# 2. Helper function definitions
# 3. Main processing function for simulation aggregation
# 4. Cross-configuration analysis and summary tables
# 5. Data preparation for visualization
# 6. Figure generation (Weibull/loglogistic hazards)
#
# OUTPUTS:
# - Combined summary tables (Stan estimates by scenario)
# - Alternative AFT method summaries for AFT-only runs
# - Boxplots of posterior means by configuration and scenario
# - R-hat convergence diagnostics
# - Baseline hazard reconstruction plots
#
# =============================================================================

# Pre-processing
##########

# Load relevant libraries
library(tidyverse)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(conflicted)
library(here)

rm(list=ls())

# =============================================================================
# CONFIGURATION: Select which models to plot for baseline hazards
# =============================================================================
# Options: "BP", "GB", "GB_Quantile", "GP", "BP_aft", "BP_aft_alt",
#          "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive"
# Example: models_to_plot <- c("BP", "GP")  # for BP and GP
#          models_to_plot <- c("BP", "GB")  # for BP and GB
# models_to_plot <- c("BP")  # <-- CHANGE THIS to select models
models_to_plot <- c("BP")
distributions_to_analyse <- c("Weibull", "loglogistic") # use c("Weibull", "loglogistic") to analyse both
distribution_pattern <- paste(distributions_to_analyse, collapse = "|")
# =============================================================================

aft_only_models <- c("BP_aft", "BP_aft_alt", "BP_aft_logy", "GB_aft", "GB_aft_logy", "GB_aft_logy_adaptive")
# Cumulative baseline hazard plots currently use the AFT-only baseline_recon
# object saved by ENZAMET-DM.R. Joint-model and legacy gamma-only
# reconstructions are intentionally not used in this section.
baseline_supported_models <- c("BP_aft", "BP_aft_logy", "GB_aft", "GB_aft_logy")
has_aft_only_models <- any(models_to_plot %in% aft_only_models)
all_models_are_aft_only <- length(models_to_plot) > 0 && all(models_to_plot %in% aft_only_models)

if (has_aft_only_models && !all_models_are_aft_only) {
  stop("models_to_plot cannot mix AFT-only and non-AFT models because truth/reference values differ.")
}

##########


# File operations and setup
##########

bunya <- FALSE
analysis_output_dir <- here::here("Bunya", "ENZAMET", "output_20260523_test_100rep_W210_centered_td10_gamma5_altDM_zeroprior_initinfo") # Set manually to pin one run; empty means use newest output_* folder.

find_latest_output_dir <- function(base_dir) {
  output_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  output_dirs <- output_dirs[stringr::str_detect(
    basename(output_dirs),
    "^output_[0-9]{8}_[A-Za-z0-9._-]+$"
  )]

  if (length(output_dirs) == 0) {
    stop("No output_yyyymmdd_jobID folders found under: ", base_dir)
  }

  output_dirs[which.max(file.info(output_dirs)$mtime)]
}

if (bunya==TRUE) {
  base_output_dir <- here::here("JoMoNoPH-share", "Bunya", "ENZAMET")
  fig_dir <- here::here("JoMoNoPH-share", "Bunya", "ENZAMET", "figures")
  summary_dir <- here::here("JoMoNoPH-share", "Bunya", "ENZAMET", "summary_tables")
} else {
  base_output_dir <- here::here("Bunya", "ENZAMET")
  fig_dir <- here::here("Bunya", "ENZAMET", "figures")
  summary_dir <- here::here("Bunya", "ENZAMET", "summary_tables")
}

out_dir <- if (nzchar(analysis_output_dir)) {
  analysis_output_dir
} else if (nzchar(Sys.getenv("SIM_OUTPUT_DIR", unset = ""))) {
  Sys.getenv("SIM_OUTPUT_DIR")
} else {
  find_latest_output_dir(base_output_dir)
}
out_dir <- path.expand(out_dir)
message("Using simulation output directory: ", out_dir)

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

# out_dir <- "C:/Users/uqamar43/OneDrive - The University of Queensland/02 shared HERA ULTRA/04 Stat Methods/08 JoMoNoPH/Biometrical-Journal/Sim-Results"


##########


# Helper functions
##########

format_cens_short <- function(cens) {
  cens <- as.character(cens)
  dplyr::case_when(
    cens == "admin" ~ "admin",
    stringr::str_detect(cens, "^\\d+$") ~ paste0(cens, "%"),
    TRUE ~ cens
  )
}

format_cens_label <- function(cens) {
  cens <- as.character(cens)
  dplyr::case_when(
    cens == "admin" ~ "Administrative censoring only",
    stringr::str_detect(cens, "^\\d+$") ~ paste0(cens, "% censoring"),
    TRUE ~ cens
  )
}

format_cens_table <- function(cens) {
  cens <- as.character(cens)
  dplyr::case_when(
    cens == "admin" ~ "Admin only",
    stringr::str_detect(cens, "^\\d+$") ~ paste0(cens, "%"),
    TRUE ~ cens
  )
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

truth_values_for_scenario <- function(scenario, aft_only = FALSE) {
  stopifnot(scenario %in% 1:5)

  if (aft_only) {
    return(c(
      "beta_long[1]" = 0,
      "beta_long[2]" = 0,
      "beta_long[3]" = c(0, 0, 0, 0, 0)[scenario],
      "sd_1_long[1]" = 0,
      "sd_1_long[2]" = 0,
      "sigma_long" = 0,
      "alpha" = 0,
      "beta_surv[1]" = c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario]
    ))
  }

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

summary_variable_pattern <- "alpha|beta_|sigma_long|sd_1_long|gamma\\[|gp_length_scale|gp_marginal_sd|f_gp"

# --- Helper: make short, readable labels for facets/plots
shorten_config <- function(cfg) {
  # Updated pattern to handle BP, GB, GB_Quantile, GP and scenario
  parts <- stringr::str_match(
    cfg,
    "^([A-Za-z_]+)_k(\\d+)_n(\\d+)_cens([A-Za-z0-9]+)_([A-Za-z0-9_]+)_scen(\\d+)_?(FU\\d+)?$"
  )
  if (is.na(parts[1, 1])) return(cfg)
  
  bf   <- parts[2]
  k    <- parts[3]
  n    <- parts[4]
  cens <- parts[5]
  baseline <- parse_baseline_label(parts[6])
  dist_short <- baseline$dist_short
  scen <- parts[7]
  
  paste0(bf, " k=", k, ", c=", format_cens_short(cens), ", d=", dist_short, ", S", scen)
}


# Function to parse a config string 
parse_config <- function(config_str) {
  # Updated pattern to handle BP, GB, GB_Quantile, GP and scenario
  parts <- stringr::str_match(config_str, "^([A-Za-z_]+)_k(\\d+)_n(\\d+)_cens([A-Za-z0-9]+)_([A-Za-z0-9_]+)_scen(\\d+)_?(FU\\d+)?$")
  
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

parse_batch_number <- function(path) {
  batch_match <- stringr::str_match(basename(path), "_batch(\\d+)\\.rds$")
  if (is.na(batch_match[1, 1])) return(NA_integer_)
  as.integer(batch_match[1, 2])
}

summarise_runtime_table <- function(runtime_df, config_label = NULL) {
  expected_cols <- c(
    "total_time_min", "joint_stan_time_min", "lme_time_min",
    "joint_compile_time_min"
  )
  if (is.null(runtime_df) || nrow(runtime_df) == 0) {
    runtime_df <- tibble::tibble(
      total_time_min = numeric(),
      joint_stan_time_min = numeric(),
      lme_time_min = numeric(),
      joint_compile_time_min = numeric()
    )
  } else {
    for (col in expected_cols) {
      if (!col %in% names(runtime_df)) runtime_df[[col]] <- NA_real_
    }
  }

  stat_or_na <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    fun(x)
  }

  out <- tibble::tibble(
    n_sims_with_runtime = if (is.null(runtime_df) || nrow(runtime_df) == 0) {
      0L
    } else {
      sum(is.finite(runtime_df$total_time_min))
    },
    mean_total_time_min = stat_or_na(runtime_df$total_time_min, mean),
    median_total_time_min = stat_or_na(runtime_df$total_time_min, stats::median),
    sd_total_time_min = stat_or_na(runtime_df$total_time_min, stats::sd),
    min_total_time_min = stat_or_na(runtime_df$total_time_min, min),
    max_total_time_min = stat_or_na(runtime_df$total_time_min, max),
    q90_total_time_min = stat_or_na(
      runtime_df$total_time_min,
      function(x) stats::quantile(x, 0.90, names = FALSE)
    ),
    q95_total_time_min = stat_or_na(
      runtime_df$total_time_min,
      function(x) stats::quantile(x, 0.95, names = FALSE)
    ),
    mean_joint_stan_time_min = stat_or_na(runtime_df$joint_stan_time_min, mean),
    max_joint_stan_time_min = stat_or_na(runtime_df$joint_stan_time_min, max),
    mean_lme_time_min = stat_or_na(runtime_df$lme_time_min, mean),
    mean_joint_compile_time_min = stat_or_na(runtime_df$joint_compile_time_min, mean)
  )

  if (!is.null(config_label)) {
    out <- dplyr::bind_cols(tibble::tibble(config = config_label), out)
  }
  out
}

empty_prior_info_per_sim <- function() {
  tibble::tibble(
    sim_id = integer(),
    unique_seed = numeric(),
    longitudinal_response_center = numeric(),
    s_long = numeric(),
    s_surv = numeric(),
    alpha_tilde_prior_mean = numeric(),
    alpha_tilde_sd = numeric(),
    alpha_prior_mean_raw = numeric(),
    alpha_prior_sd_raw = numeric(),
    beta_long_prior_mean_1 = numeric(),
    beta_long_prior_mean_2 = numeric(),
    beta_long_prior_mean_3 = numeric(),
    beta_long_prior_mean_raw_1 = numeric(),
    beta_long_prior_mean_raw_2 = numeric(),
    beta_long_prior_mean_raw_3 = numeric(),
    beta_long_prior_scale_1 = numeric(),
    beta_long_prior_scale_2 = numeric(),
    beta_long_prior_scale_3 = numeric(),
    beta_surv_prior_mean_1 = numeric(),
    beta_surv_prior_scale_1 = numeric(),
    sd_b_prior_mean_1 = numeric(),
    sd_b_prior_mean_2 = numeric(),
    sd_b_prior_scale_1 = numeric(),
    sd_b_prior_scale_2 = numeric()
  )
}

safe_prior_scalar <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0) return(NA_real_)
  x[1]
}

prior_vector_cols <- function(x, prefix, n) {
  vals <- rep(NA_real_, n)
  x <- suppressWarnings(as.numeric(x))
  if (length(x) > 0) {
    n_fill <- min(length(x), n)
    vals[seq_len(n_fill)] <- x[seq_len(n_fill)]
  }
  tibble::as_tibble_row(stats::setNames(as.list(vals), paste0(prefix, seq_len(n))))
}

extract_prior_info_row <- function(sim, idx) {
  prior_info <- sim$prior_info
  if (is.null(prior_info)) return(NULL)

  dplyr::bind_cols(
    tibble::tibble(
      sim_id = idx,
      unique_seed = sim$unique_seed,
      longitudinal_response_center = safe_prior_scalar(prior_info$longitudinal_response_center),
      s_long = safe_prior_scalar(prior_info$s_long),
      s_surv = safe_prior_scalar(prior_info$s_surv),
      alpha_tilde_prior_mean = safe_prior_scalar(prior_info$alpha_tilde_prior_mean),
      alpha_tilde_sd = safe_prior_scalar(prior_info$alpha_tilde_sd),
      alpha_prior_mean_raw = safe_prior_scalar(prior_info$alpha_prior_mean_raw),
      alpha_prior_sd_raw = safe_prior_scalar(prior_info$alpha_prior_sd_raw)
    ),
    prior_vector_cols(prior_info$beta_long_prior_mean, "beta_long_prior_mean_", 3),
    prior_vector_cols(prior_info$beta_long_prior_mean_raw, "beta_long_prior_mean_raw_", 3),
    prior_vector_cols(prior_info$beta_long_prior_scale, "beta_long_prior_scale_", 3),
    prior_vector_cols(prior_info$beta_surv_prior_mean, "beta_surv_prior_mean_", 1),
    prior_vector_cols(prior_info$beta_surv_prior_scale, "beta_surv_prior_scale_", 1),
    prior_vector_cols(prior_info$sd_b_prior_mean, "sd_b_prior_mean_", 2),
    prior_vector_cols(prior_info$sd_b_prior_scale, "sd_b_prior_scale_", 2)
  )
}

summarise_prior_info_table <- function(prior_info_df, config_label = NULL) {
  if (is.null(prior_info_df) || nrow(prior_info_df) == 0) {
    out <- tibble::tibble(n_sims_with_prior_info = 0L)
    if (!is.null(config_label)) {
      out <- dplyr::bind_cols(tibble::tibble(config = config_label), out)
    }
    return(out)
  }

  stat_or_na <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    fun(x)
  }

  numeric_cols <- setdiff(names(prior_info_df), c("sim_id", "unique_seed"))
  summary_cols <- purrr::map(
    numeric_cols,
    function(col) {
      x <- prior_info_df[[col]]
      tibble::tibble(
        !!paste0(col, "_mean") := stat_or_na(x, mean),
        !!paste0(col, "_median") := stat_or_na(x, stats::median),
        !!paste0(col, "_sd") := stat_or_na(x, stats::sd),
        !!paste0(col, "_min") := stat_or_na(x, min),
        !!paste0(col, "_max") := stat_or_na(x, max)
      )
    }
  )

  out <- dplyr::bind_cols(
    tibble::tibble(n_sims_with_prior_info = nrow(prior_info_df)),
    !!!summary_cols
  )
  if (!is.null(config_label)) {
    out <- dplyr::bind_cols(tibble::tibble(config = config_label), out)
  }
  out
}

parse_short_config <- function(lbl) {
  m <- stringr::str_match(lbl, "^([A-Za-z_]+)\\s+k=(\\d+),\\s+c=([A-Za-z0-9]+)%?(?:,\\s+d=([^,]+))?,\\s+S(\\d+)$")
  if (is.na(m[1,1])) {
    return(tibble::tibble(
      bf = NA_character_,
      k_bases = NA_real_,
      cens = NA_character_,
      baseline_label = NA_character_,
      dist = NA_character_,
      dist_short = NA_character_,
      scenario = NA_real_
    ))
  }
  dist_token <- m[1,5]
  dist <- dplyr::case_when(
    stringr::str_detect(dist_token, "^W") ~ "Weibull",
    stringr::str_detect(dist_token, "^LL") ~ "loglogistic",
    TRUE ~ dist_token
  )
  tibble::tibble(
    bf       = m[1,2],
    k_bases  = as.numeric(m[1,3]),
    cens     = m[1,4],
    baseline_label = dist_token,
    dist     = dist,
    dist_short = dist_token,
    scenario = as.numeric(m[1,6])
  )
}

##########


# Main processing function
##########

# This function reads, processes, and aggregates simulation results
# from Bayesian joint modelling runs

# --- Main: process all sims for a single config
process_simulation_results <- function(bf, k_bases, n, cens, baseline_label, dist, scenario, fu,
                                       output_dir) {
  
  # Build filename pattern for this config
  fu_part <- if (!is.null(fu) && !is.na(fu) && nzchar(fu)) {
    paste0("_", fu)
  } else {
    ""
  }
  
  pattern <- paste0(
    bf, "_k", k_bases, "_n", n, "_cens", cens, "_", baseline_label,
    "_scen", scenario,
    fu_part,
    "(?:_.*)?_batch\\d{3}\\.rds$"
  )
  
  
  rds_files <- list.files(path = output_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  
  if (length(rds_files) == 0) {
    message("No matching files found for: ", pattern)
    return(NULL)
  }
  rds_files <- sort(rds_files)

  # browser()
  # ---- Read all batches while preserving batch/replicate positions.
  # Failed/time-limited replicates are saved as NULL entries; if we compact first,
  # their exact batch and within-batch position are lost.
  batch_contents <- purrr::map(rds_files, readRDS)
  batch_sizes <- vapply(batch_contents, length, integer(1))
  nominal_batch_size <- if (length(batch_sizes) == 0 || max(batch_sizes, na.rm = TRUE) == 0) {
    NA_integer_
  } else {
    max(batch_sizes, na.rm = TRUE)
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
        is_null = is.null(batch_sims),
        sim = list(NULL)
      ))
    }

    purrr::map_dfr(seq_along(batch_sims), function(pos) {
      sim <- batch_sims[[pos]]
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
        is_null = is.null(sim),
        sim = list(sim)
      )
    })
  })

  n_raw_sims <- nrow(batch_records)
  n_failed_sims <- sum(batch_records$is_null, na.rm = TRUE)
  failed_replicates <- batch_records %>%
    dplyr::filter(is_null) %>%
    dplyr::transmute(
      batch_file,
      batch_number,
      position_in_batch,
      replicate_in_config,
      failure_type = "null_saved_replicate"
    )
  sims <- batch_records %>%
    dplyr::filter(!is_null) %>%
    dplyr::pull(sim)
  
  if (length(sims) == 0) {
    message("All simulations were NULL/failed for: ", pattern)
    empty_summary <- tibble::tibble(
      variable = character(), truth = numeric(), mean = numeric(), median = numeric(),
      empirical_sd = numeric(), mean_post_sd = numeric(), sd = numeric(),
      bias = numeric(), rmse = numeric(), mean_ci_width = numeric(),
      coverage = numeric(), mcse_mean = numeric(), mcse_coverage = numeric(),
      n_sims = integer(), source = character()
    )
    empty_filtered_data <- tibble::tibble(
      variable = character(), mean = numeric(), median = numeric(), sd = numeric(),
      q2.5 = numeric(), q97.5 = numeric(), rhat = numeric(), sim_id = integer()
    )
    empty_diag_per_sim <- tibble::tibble(
      sim_id = integer(), censor_prop = numeric(), unique_seed = numeric(),
      max_rhat = numeric(), mean_rhat = numeric(), median_rhat = numeric(),
      sd_rhat = numeric(), n_rhat_gt_1_01 = numeric(), n_rhat_gt_1_05 = numeric()
    )
    empty_runtime_per_sim <- tibble::tibble(
      sim_id = integer(), unique_seed = numeric(), lme_time_sec = numeric(),
      joint_compile_time_sec = numeric(), joint_stan_time_sec = numeric(),
      total_time_sec = numeric(), lme_time_min = numeric(),
      joint_compile_time_min = numeric(), joint_stan_time_min = numeric(),
      total_time_min = numeric()
    )
    empty_prior_info_per_sim_tbl <- empty_prior_info_per_sim()
    return(list(
      fig           = NULL,
      table         = NULL,
      filtered_data = empty_filtered_data,
      summary       = empty_summary,
      diag_per_sim  = empty_diag_per_sim,
      diag_summary  = tibble::tibble(
        mean_censor_prop = NA_real_,
        sd_censor_prop   = NA_real_,
        mean_max_rhat    = NA_real_,
        max_max_rhat     = NA_real_,
        mean_mean_rhat   = NA_real_,
        mean_median_rhat = NA_real_,
        mean_n_rhat_gt_1_01 = NA_real_,
        mean_n_rhat_gt_1_05 = NA_real_,
        n_sims = 0L,
        n_batch_files = length(rds_files),
        n_raw_sims = n_raw_sims,
        n_failed_sims = n_failed_sims,
        n_valid_sims = 0L
      ),
      runtime_per_sim = empty_runtime_per_sim,
      runtime_summary = tibble::tibble(
        n_sims_with_runtime = 0L,
        mean_total_time_min = NA_real_,
        median_total_time_min = NA_real_,
        sd_total_time_min = NA_real_,
        min_total_time_min = NA_real_,
        max_total_time_min = NA_real_,
        q90_total_time_min = NA_real_,
        q95_total_time_min = NA_real_,
        mean_joint_stan_time_min = NA_real_,
        max_joint_stan_time_min = NA_real_,
        mean_lme_time_min = NA_real_,
        mean_joint_compile_time_min = NA_real_
      ),
      prior_info_per_sim = empty_prior_info_per_sim_tbl,
      prior_info_summary = summarise_prior_info_table(empty_prior_info_per_sim_tbl),
      failed_replicates = failed_replicates,
      beta_saft_raw     = NULL,
      beta_saft_summary = NULL
    ))
  }
  
  # ---- Aggregate Stan summaries (keep parameter names + sim id)
  all_summaries_list <- sims %>%
    purrr::imap(function(sim, idx) {
      s <- sim$joint_summary
      
      # Skip if summary is NULL or empty
      if (is.null(s)) return(NULL)
      
      s_df <- as.data.frame(s)
      if (nrow(s_df) == 0) return(NULL)
      s_df <- s_df %>%
        dplyr::filter(stringr::str_detect(variable, summary_variable_pattern))
      if (nrow(s_df) == 0) return(NULL)
      
      # joint_summary already has a variable column, don't add it again
      if (!"sim_id" %in% colnames(s_df)) {
        s_df$sim_id <- idx
      }
      s_df
    }) %>%
    purrr::compact()
  
  # If nothing valid, bail out for this config
  if (length(all_summaries_list) == 0) {
    message("No valid Stan summaries for: ", pattern)
    return(NULL)
  }
  
  all_summaries <- dplyr::bind_rows(all_summaries_list)
  
  beta_saft_per_sim <- NULL
  beta_saft_summary <- NULL
  
  # ---- Extract semiparametric / parametric AFT beta estimates per sim
  if (bf %in% aft_only_models) {
    beta_saft_per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
      
      if (is.null(sim$beta_sAFT)) return(NULL)
      
      beta_mat <- sim$beta_sAFT
      beta_mat <- as.matrix(beta_mat)
      method_names <- c("aftgee", "aftsrr", "smoothSurv", "rstpm2", "weibull", "loglogistic")
      if (nrow(beta_mat) <= length(method_names)) {
        rownames(beta_mat) <- method_names[seq_len(nrow(beta_mat))]
      }
      
      # Ensure matrix with method names
      beta_df <- as.data.frame(beta_mat)
      beta_df$method <- rownames(beta_mat)
      
      beta_df %>%
        pivot_longer(
          cols = -method,
          names_to = "parameter",
          values_to = "estimate"
        ) %>%
        mutate(
          parameter = dplyr::if_else(parameter == "arm", "beta_surv[1]", parameter),
          sim_id = idx
        )
    })
    # Attach true values (scenario-dependent)
    true_beta_aft <- c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario]
    # enrich the long table
    if (nrow(beta_saft_per_sim) > 0) {
      beta_saft_per_sim <- beta_saft_per_sim %>%
        mutate(truth = true_beta_aft, bias = estimate - truth)
      # Aggregate across simulation replicates
      beta_saft_summary <- beta_saft_per_sim %>%
        group_by(method, parameter) %>%
        summarise(
          Mean     = mean(estimate, na.rm = TRUE),
          Median   = median(estimate, na.rm = TRUE),
          Q2.5     = stats::quantile(estimate, 0.025, na.rm = TRUE),
          Q97.5    = stats::quantile(estimate, 0.975, na.rm = TRUE),
          SD       = sd(estimate, na.rm = TRUE),
          Bias     = mean(bias, na.rm = TRUE),
          RMSE     = sqrt(mean(bias^2, na.rm = TRUE)),
          n_sims   = sum(!is.na(estimate)),
          .groups  = "drop"
        ) %>%
        mutate(Source = "AFT (frequentist)")
    }
  }
  
  
  # ---- Extract censoring proportion and R-hat diagnostics per sim
  diag_per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
    
    # Skip failed sims
    if (is.null(sim$joint_summary)) return(NULL)
    
    rhat <- sim$joint_summary$rhat
    rhat <- rhat[is.finite(rhat)]
    
    tibble::tibble(
      sim_id          = idx,
      censor_prop     = sim$censor_prop,
      unique_seed     = sim$unique_seed,
      max_rhat        = max(rhat, na.rm = TRUE),
      mean_rhat       = mean(rhat, na.rm = TRUE),
      median_rhat     = median(rhat, na.rm = TRUE),
      sd_rhat         = sd(rhat, na.rm = TRUE),
      n_rhat_gt_1_01  = sum(rhat > 1.01, na.rm = TRUE),
      n_rhat_gt_1_05  = sum(rhat > 1.05, na.rm = TRUE)
    )
  })
  # ---- Aggregate diagnostics across replicates
  diag_summary <- diag_per_sim %>%
    summarise(
      mean_censor_prop = mean(censor_prop),
      sd_censor_prop   = sd(censor_prop),
      
      mean_max_rhat    = mean(max_rhat),
      max_max_rhat     = max(max_rhat),
      
      mean_mean_rhat   = mean(mean_rhat),
      mean_median_rhat = mean(median_rhat),
      
      mean_n_rhat_gt_1_01 = mean(n_rhat_gt_1_01),
      mean_n_rhat_gt_1_05 = mean(n_rhat_gt_1_05),
      
      n_sims = n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_batch_files = length(rds_files),
      n_raw_sims = n_raw_sims,
      n_failed_sims = n_failed_sims,
      n_valid_sims = length(sims)
    )

  # ---- Extract runtime information saved by ENZAMET-DM.R.
  # Times are stored in seconds in each replicate; summaries below report minutes
  # because those are easier to interpret for SLURM planning.
  runtime_per_sim <- purrr::imap_dfr(sims, function(sim, idx) {
    rt <- sim$runtimes
    if (is.null(rt) || nrow(rt) == 0) return(NULL)

    lme_time <- if ("lme_time" %in% names(rt)) as.numeric(rt$lme_time[1]) else NA_real_
    joint_compile_time <- if ("joint_compile_time" %in% names(rt)) as.numeric(rt$joint_compile_time[1]) else NA_real_
    joint_stan_time <- if ("joint_stan_time" %in% names(rt)) as.numeric(rt$joint_stan_time[1]) else NA_real_
    total_time <- if ("total_time" %in% names(rt)) as.numeric(rt$total_time[1]) else NA_real_

    tibble::tibble(
      sim_id = idx,
      unique_seed = sim$unique_seed,
      lme_time_sec = lme_time,
      joint_compile_time_sec = joint_compile_time,
      joint_stan_time_sec = joint_stan_time,
      total_time_sec = total_time,
      lme_time_min = lme_time / 60,
      joint_compile_time_min = joint_compile_time / 60,
      joint_stan_time_min = joint_stan_time / 60,
      total_time_min = total_time / 60
    )
  })

  runtime_summary <- summarise_runtime_table(runtime_per_sim)

  # ---- Extract empirical prior/calibration information saved by ENZAMET-DM.R.
  # This is mainly for reviewer-facing sensitivity checks of alpha and staging
  # priors. These fields document the priors used for fitting; they are not
  # posterior estimates.
  prior_info_per_sim <- purrr::imap_dfr(sims, extract_prior_info_row)
  if (nrow(prior_info_per_sim) == 0) {
    prior_info_per_sim <- empty_prior_info_per_sim()
  }
  prior_info_summary <- summarise_prior_info_table(prior_info_per_sim)
  
  
  # Filter parameters of interest (including GP-specific parameters)
  filtered_data <- all_summaries
  
  message("Total simulations found for ", bf, "_n", n, "_cens", cens, "_", baseline_label, "_scen", scenario, "_", fu,
          ": valid N=", length(unique(filtered_data$sim_id)), ", failed/null=", n_failed_sims)
  
  # Quick histogram per parameter (diagnostic)
  fig_hist <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = mean)) +
    ggplot2::geom_histogram(bins = 30, color = "black", fill = "steelblue") +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::labs(title = "Distribution of Posterior Means (per parameter)")
  
  # Define true parameter values for ENZAMET
  # Note: beta_long[3] and beta_surv[1] vary by scenario
  # Common values across all scenarios:
  # beta_0 = 73, beta_1 = -0.04, sigma_e = 12
  # sigma_0 = 15, sigma_1 = 0.20, alpha_AFT = 0.012
  # Scenario-specific:
  # Scenario 1: beta_2 = 0.00, log_AF = 0.00
  # Scenario 2: beta_2 = 0.04, log_AF = 0.90
  # Scenario 3: beta_2 = -0.04, log_AF = 0.90
  # Scenario 4: beta_2 = 0.04, log_AF = -0.90
  # Scenario 5: beta_2 = -0.04, log_AF = -0.90
  #### TRUE VALUES (scenarios 1–5 only), from supplied true_params ####
  true_values <- truth_values_for_scenario(scenario, aft_only = all_models_are_aft_only)
  
  # Compute coverage: check if true value is within [q2.5, q97.5] for each sim
  coverage_data <- filtered_data %>%
    dplyr::filter(variable %in% names(true_values)) %>%
    dplyr::mutate(
      true_value = true_values[variable],
      covered = (q2.5 <= true_value) & (true_value <= q97.5)
    ) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      coverage = mean(covered, na.rm = TRUE),
      n_sims = sum(!is.na(covered)),
      .groups = "drop"
    )
  
  # Simulation summary across replicates.
  # CP is based on each replicate's posterior 95% credible interval.
  post_summary <- filtered_data %>%
    dplyr::mutate(
      posterior_mean = .data$mean,
      posterior_sd = .data$sd,
      truth = unname(true_values[variable]),
      error = posterior_mean - truth,
      covered = dplyr::if_else(
        !is.na(truth),
        q2.5 <= truth & truth <= q97.5,
        NA
      ),
      ci_width = q97.5 - q2.5
    ) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      truth = dplyr::first(stats::na.omit(truth), default = NA_real_),
      n_sims = sum(!is.na(.data$posterior_mean)),
      mean = mean(.data$posterior_mean, na.rm = TRUE),
      median = median(.data$posterior_mean, na.rm = TRUE),
      empirical_sd = sd(.data$posterior_mean, na.rm = TRUE),
      mean_post_sd = mean(.data$posterior_sd, na.rm = TRUE),
      sd = mean_post_sd,  # legacy alias for older downstream code
      q2.5_est = stats::quantile(.data$posterior_mean, 0.025, na.rm = TRUE),
      q97.5_est = stats::quantile(.data$posterior_mean, 0.975, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      coverage = mean(covered, na.rm = TRUE),
      mcse_mean = empirical_sd / sqrt(n_sims),
      mcse_coverage = sqrt(coverage * (1 - coverage) / n_sims),
      .groups = "drop"
    ) %>%
    dplyr::mutate(source = "Stan posterior")
  
  # ENZAMET files don't have LMM - skip LMM processing
  
  # Combined summary is just Stan posterior
  combined_summary <- post_summary %>%
    dplyr::select(
      variable, truth, mean, median, empirical_sd, mean_post_sd, sd,
      bias, rmse, mean_ci_width, coverage, mcse_mean, mcse_coverage, n_sims, source
    ) %>%
    dplyr::arrange(source, variable)
  
  # Flextable: rounded values + simple header
  table <- combined_summary %>%
    dplyr::mutate(
      mean = round(mean, 3),
      empirical_sd = round(empirical_sd, 3),
      mean_post_sd = round(mean_post_sd, 3),
      bias = round(bias, 3),
      rmse = round(rmse, 3),
      mean_ci_width = round(mean_ci_width, 3),
      coverage = round(coverage, 3),
      mcse_mean = round(mcse_mean, 3),
      mcse_coverage = round(mcse_coverage, 3)
    ) %>%
    flextable::flextable() %>%
    flextable::align(j = "source", align = "left")
  
  list(
    fig           = fig_hist,        # diagnostic histograms
    table         = table,           # Stan summary table
    filtered_data = filtered_data,   # long Stan posterior means per sim
    summary       = combined_summary, # aggregated table (Stan only)
    diag_per_sim  = diag_per_sim,    # NEW
    diag_summary  = diag_summary,    # NEW
    runtime_per_sim = runtime_per_sim,
    runtime_summary = runtime_summary,
    prior_info_per_sim = prior_info_per_sim,
    prior_info_summary = prior_info_summary,
    failed_replicates = failed_replicates,
    beta_saft_raw    = beta_saft_per_sim,     # NEW
    beta_saft_summary = beta_saft_summary     # NEW
  )
}

##########


# Configuration parsing and execution
##########

# ===== Run across all configs =====

# Get all unique config stubs from files (strip last 9 chars before .rds)
stopifnot(dir.exists(out_dir))
filenames <- list.files(out_dir, pattern = "\\.[rR][dD][sS]$")
configs <- substr(tools::file_path_sans_ext(filenames),
                  1, nchar(tools::file_path_sans_ext(filenames)) - 9) %>%
  unique()
configs <- configs[
  purrr::map_lgl(configs, function(cfg) {
    cfg_parsed <- parse_config(cfg)
    !is.null(cfg_parsed) &&
      cfg_parsed$bf %in% models_to_plot &&
      cfg_parsed$dist %in% distributions_to_analyse
  })
]

# # Get all unique config stubs from BP files only (strip last 9 chars before .rds)
# filenames <- list.files(
#   out_dir,
#   pattern = "^BP_.*_loglogistic_scen[1-5]_.*\\.rds$",
#   full.names = FALSE
# )
# 
# configs <- substr(
#   tools::file_path_sans_ext(filenames),
#   1, nchar(tools::file_path_sans_ext(filenames)) - 9
# ) |>
#   unique()


print(configs)

# Parse, run, and keep only successful configs
results_by_config <- purrr::compact(
  setNames(configs, configs) %>%
    purrr::map(function(cfg) {
      cfg_parsed <- parse_config(cfg)
      if (is.null(cfg_parsed)) return(NULL)
      process_simulation_results(
        bf      = cfg_parsed$bf,
        k_bases = cfg_parsed$k_bases,
        n       = cfg_parsed$n,
        cens    = cfg_parsed$cens,
        baseline_label = cfg_parsed$baseline_label,
        dist    = cfg_parsed$dist,
        scenario = cfg_parsed$scenario,
        fu      = cfg_parsed$fu,
        output_dir = out_dir
      )
    })
)

config_info_full <- if (length(results_by_config) > 0) {
  tibble::tibble(config = names(results_by_config)) %>%
    dplyr::mutate(parsed = purrr::map(config, parse_config)) %>%
    tidyr::unnest_wider(parsed)
} else {
  tibble::tibble(config = character(), bf = character(), k_bases = numeric(),
                 n = numeric(), cens = character(), baseline_label = character(),
                 dist = character(), dist_short = character(),
                 baseline_shape = numeric(), baseline_scale = numeric(),
                 scenario = numeric(), fu = character())
}
config_info <- config_info_full

analysis_dist_title <- paste(sort(unique(config_info_full$dist)), collapse = "/")
analysis_dist_file <- tolower(gsub("[^A-Za-z0-9]+", "_", analysis_dist_title))
if (!nzchar(analysis_dist_file)) {
  analysis_dist_file <- "all"
  analysis_dist_title <- "All distributions"
}


# ---- Combine diagnostic summaries across configs
all_diag_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$diag_summary, config = .y)
)
all_diag_summary
# mean_n_rhat_gt_1_01 column: mean number of r_hat's greater than 1.01
# mean_n_rhat_gt_1_05 column: mean number of r_hat's greater than 1.05

readr::write_csv(
  all_diag_summary,
  file.path(summary_dir, paste0("jm_diagnostic_summary_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_diag_summary,
  file.path(summary_dir, paste0("jm_diagnostic_summary_", analysis_dist_file, ".rds"))
)

# diagnostic plot: boxplot of max R-hat per config
# long form, useful for boxplots of max_rhat or realised censoring
all_diag_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$diag_per_sim, config = .y)
)

readr::write_csv(
  all_diag_per_sim,
  file.path(summary_dir, paste0("jm_diagnostic_per_sim_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_diag_per_sim,
  file.path(summary_dir, paste0("jm_diagnostic_per_sim_", analysis_dist_file, ".rds"))
)

# ---- Failed/null replicate locations for rescue checks
all_failed_replicates <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$failed_replicates, config = .y, .before = 1)
)

readr::write_csv(
  all_failed_replicates,
  file.path(summary_dir, paste0("jm_failed_null_replicates_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_failed_replicates,
  file.path(summary_dir, paste0("jm_failed_null_replicates_", analysis_dist_file, ".rds"))
)

# ---- Runtime summaries across configurations
all_runtime_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$runtime_per_sim, config = .y)
)

all_runtime_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$runtime_summary, config = .y)
)

overall_runtime_summary <- summarise_runtime_table(all_runtime_per_sim, config_label = "OVERALL")

all_runtime_summary_with_overall <- dplyr::bind_rows(
  all_runtime_summary,
  overall_runtime_summary
)

readr::write_csv(
  all_runtime_per_sim,
  file.path(summary_dir, paste0("jm_runtime_per_sim_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_runtime_per_sim,
  file.path(summary_dir, paste0("jm_runtime_per_sim_", analysis_dist_file, ".rds"))
)
readr::write_csv(
  all_runtime_summary_with_overall,
  file.path(summary_dir, paste0("jm_runtime_summary_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_runtime_summary_with_overall,
  file.path(summary_dir, paste0("jm_runtime_summary_", analysis_dist_file, ".rds"))
)

# ---- Empirical prior/calibration information across configurations
all_prior_info_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$prior_info_per_sim, config = .y, .before = 1)
)

all_prior_info_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$prior_info_summary, config = .y)
)

overall_prior_info_input <- if ("config" %in% names(all_prior_info_per_sim)) {
  dplyr::select(all_prior_info_per_sim, -config)
} else {
  all_prior_info_per_sim
}
overall_prior_info_summary <- summarise_prior_info_table(
  overall_prior_info_input,
  config_label = "OVERALL"
)

all_prior_info_summary_with_overall <- dplyr::bind_rows(
  all_prior_info_summary,
  overall_prior_info_summary
)

readr::write_csv(
  all_prior_info_per_sim,
  file.path(summary_dir, paste0("jm_prior_info_per_sim_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_prior_info_per_sim,
  file.path(summary_dir, paste0("jm_prior_info_per_sim_", analysis_dist_file, ".rds"))
)
readr::write_csv(
  all_prior_info_summary_with_overall,
  file.path(summary_dir, paste0("jm_prior_info_summary_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_prior_info_summary_with_overall,
  file.path(summary_dir, paste0("jm_prior_info_summary_", analysis_dist_file, ".rds"))
)

shorten_config_1 <- function(config) {
  purrr::map_chr(config, function(cfg) {
    parsed <- parse_config(cfg)
    if (is.null(parsed)) return(cfg)
    paste(parsed$bf, paste0("k", parsed$k_bases), paste0("n", parsed$n),
          parsed$dist_short, parsed$fu, sep = ", ")
  })
}
all_diag_per_sim <- all_diag_per_sim %>%
  dplyr::mutate(cfg_short = shorten_config_1(config))
all_diag_per_sim <- all_diag_per_sim %>%
  dplyr::mutate(
    cens = stringr::str_extract(config, "cens(?:\\d+|admin)"),
    scen = stringr::str_extract(config, "scen\\d+")
  )
all_diag_per_sim <- all_diag_per_sim %>%
  dplyr::mutate(
    cens = format_cens_label(stringr::str_remove(cens, "^cens")),
    scen = factor(scen)
  )
ggplot(all_diag_per_sim, aes(x = cfg_short, y = max_rhat)) +
  geom_boxplot() +
  # coord_flip() +
  facet_grid(cens ~ scen) +
  labs(
    y = "Maximum R-hat per replicate",
    x = "Model configuration"
  )
# diagnostic plot: histogram of censoring proportions across sims, faceted by config
ggplot(all_diag_per_sim, aes(x = censor_prop)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ shorten_config_1(config))


# Build a long table of summaries with the config label attached
all_summary_tables <- purrr::map2_dfr(
  results_by_config, names(results_by_config),
  ~ dplyr::mutate(.x$summary, config = .y),
  .id = NULL
)
print(all_summary_tables, n = nrow(all_summary_tables))

readr::write_csv(
  all_summary_tables,
  file.path(summary_dir, paste0("jm_parameter_summary_long_", analysis_dist_file, ".csv"))
)
saveRDS(
  all_summary_tables,
  file.path(summary_dir, paste0("jm_parameter_summary_long_", analysis_dist_file, ".rds"))
)

# Wider view (mean/sd columns per config)
if (nrow(all_summary_tables) > 0 && any(c("mean", "sd") %in% colnames(all_summary_tables))) {
  summary_wide <- all_summary_tables %>%
    dplyr::filter(config %in% config_info_full$config) %>%
    tidyr::pivot_longer(cols = dplyr::any_of(c("mean", "sd")), names_to = "stat", values_to = "value") %>%
    tidyr::unite("config_stat", config, stat) %>%
    tidyr::pivot_wider(names_from = config_stat, values_from = value)
} else {
  summary_wide <- tibble::tibble()
  message("No summary data available to pivot")
}

# Display a compact flextable (optional)
if (nrow(summary_wide) > 0) {
  summary_wide %>%
    dplyr::arrange(variable) %>%
    flextable::flextable()
}

# Base reference lines (scenario-independent parameters)
if (all_models_are_aft_only) {
  reference_lines <- tibble::tibble(
    variable = c("beta_long[1]","beta_long[2]","beta_long[3]",
                 "sd_1_long[1]","sd_1_long[2]","sigma_long",
                 "beta_surv[1]","alpha"),
    reference = c(0, 0, 0, 0, 0, 0, NA, 0)
  )
} else {
  reference_lines <- tibble::tibble(
    variable = c("beta_long[1]","beta_long[2]","beta_long[3]",
                 "sd_1_long[1]","sd_1_long[2]","sigma_long",
                 "beta_surv[1]","alpha"),
    reference = c(73, -0.04, NA, 15, 0.20, 12, NA, 0.012)
  )
}

key_vars <- c("alpha", "beta_long[1]", "beta_long[2]", "beta_long[3]", "beta_surv[1]")

if (nrow(all_summary_tables) > 0) {
  # Attach true values
  # Use full config metadata to determine scenario and distribution.
  config_info <- config_info_full   # bf, n, cens, dist, scenario, fu
  
  results_tab <- all_summary_tables %>%
    dplyr::left_join(config_info, by = "config") %>%
    left_join(reference_lines, by = "variable") %>%
    mutate(
      reference = dplyr::coalesce(truth, reference),
      rel_bias = if_else(reference != 0, bias / reference, NA_real_)
    )
  
  results_tab2 <- results_tab
  # publication-ready table
  nice_tab_wide <- results_tab2 %>%
    dplyr::filter(
      variable %in% key_vars,
      source == "Stan posterior"
    ) %>%
    mutate(
      Parameter = variable,
      Basis     = bf,
      k         = as.integer(k_bases),
      n         = as.integer(n),
      Cens      = format_cens_table(cens),
      Dist      = dist,
      Baseline  = baseline_label,
      Shape     = baseline_shape,
      Scale     = baseline_scale,
      Scenario  = as.integer(scenario),
      Method    = "JM (Stan)"
    ) %>%
    transmute(
      Parameter,
      Basis,
      k,
      n,
      Cens,
      Dist,
      Baseline,
      Shape,
      Scale,
      Scenario,
      Method,
      Truth    = reference,
      Mean     = mean,
      Bias     = bias,
      EmpiricalSD = empirical_sd,
      MeanPostSD = mean_post_sd,
      RelBias  = rel_bias,
      RMSE     = rmse,
      MeanCIWidth = mean_ci_width,
      CP       = coverage,
      MCSE_Mean = mcse_mean,
      MCSE_CP = mcse_coverage,
      n_sims
    ) %>%
    arrange(Parameter, Basis, k, n, Cens, Dist, Scenario, Method)
  nice_tab_wide
  
  nice_ft <- nice_tab_wide %>%
    mutate(across(c(Truth, Mean, Bias, EmpiricalSD, MeanPostSD, RelBias, RMSE, MeanCIWidth, CP, MCSE_Mean, MCSE_CP), ~ round(.x, 3))) %>%
    flextable::flextable() %>%
    flextable::autofit()
  
  print(nice_ft)
  readr::write_csv(
    nice_tab_wide,
    file.path(summary_dir, paste0("jm_simulation_stats_", analysis_dist_file, ".csv"))
  )
  saveRDS(
    nice_tab_wide,
    file.path(summary_dir, paste0("jm_simulation_stats_", analysis_dist_file, ".rds"))
  )
  
  # Create gamma summary tables for basis weights (BP/GB/GB_Quantile only)
  # GP model uses different parameters (gp_length_scale, gp_marginal_sd, f_gp)
  gamma_tab_wide <- results_tab2 %>%
    dplyr::filter(
      stringr::str_detect(variable, "^gamma\\["),
      source == "Stan posterior"
    ) %>%
    mutate(
      Parameter = variable,
      Basis     = bf,
      k         = as.integer(k_bases),
      n         = as.integer(n),
      Cens      = format_cens_table(cens),
      Dist      = dist,
      Baseline  = baseline_label,
      Scenario  = as.integer(scenario)
    ) %>%
    transmute(
      Parameter,
      Basis,
      k,
      n,
      Cens,
      Dist,
      Baseline,
      Scenario,
      Mean     = mean,
      SD       = sd,
      Coverage = coverage
    ) %>%
    arrange(Basis, k, n, Cens, Dist, Scenario, Parameter)
  
  # Create GP-specific parameter summary table
  gp_tab_wide <- results_tab2 %>%
    dplyr::filter(
      stringr::str_detect(variable, "^(gp_length_scale|gp_marginal_sd)$"),
      source == "Stan posterior"
    ) %>%
    mutate(
      Parameter = variable,
      Basis     = bf,
      k         = as.integer(k_bases),
      n         = as.integer(n),
      Cens      = format_cens_table(cens),
      Dist      = dist,
      Baseline  = baseline_label,
      Scenario  = as.integer(scenario)
    ) %>%
    transmute(
      Parameter,
      Basis,
      k,
      n,
      Cens,
      Dist,
      Baseline,
      Scenario,
      Mean     = mean,
      SD       = sd,
      Coverage = coverage
    ) %>%
    arrange(Basis, k, n, Cens, Dist, Scenario, Parameter)
  
  # Display gamma table only if there are gamma parameters (BP/GB/GB_Quantile)
  if (nrow(gamma_tab_wide) > 0) {
    gamma_tab_with_scenario <- gamma_tab_wide %>%
      mutate(
        across(c(Mean, SD, Coverage), ~ round(.x, 3)),
        Scenario_label = paste(Basis, k, n, Cens, Dist, Scenario, sep = "_")
      )
    
    # Create color mapping for unique scenarios
    unique_scenarios <- unique(gamma_tab_with_scenario$Scenario_label)
    scenario_colors <- rep(c("#F0F0F0", "#FFFFFF", "#E8F4F8", "#FFF9E6"), 
                           length.out = length(unique_scenarios))
    names(scenario_colors) <- unique_scenarios
    
    # Build flextable with scenario coloring
    gamma_ft <- gamma_tab_with_scenario %>%
      flextable::flextable() %>%
      flextable::autofit()
    
    # Apply background colors row by row
    for (i in seq_len(nrow(gamma_tab_with_scenario))) {
      scenario_id <- gamma_tab_with_scenario$Scenario_label[i]
      gamma_ft <- flextable::bg(gamma_ft, i = i, bg = scenario_colors[scenario_id])
    }
    
    # Add borders between different scenarios
    for (i in seq(2, nrow(gamma_tab_with_scenario))) {
      if (gamma_tab_with_scenario$Scenario_label[i] != gamma_tab_with_scenario$Scenario_label[i-1]) {
        gamma_ft <- flextable::border(gamma_ft, i = i, 
                                      border.top = officer::fp_border(color = "gray30", width = 2))
      }
    }
    
    # Remove the Scenario_label helper column
    gamma_ft <- flextable::delete_columns(gamma_ft, j = "Scenario_label")
    
    message("\n=== GAMMA BASIS WEIGHTS SUMMARY ===")
    print(gamma_ft)
  } else {
    message("\n=== No gamma parameters found (possibly GP models only) ===")
  }
  
  # Display GP hyperparameters table if there are GP parameters
  if (nrow(gp_tab_wide) > 0) {
    gp_tab_with_scenario <- gp_tab_wide %>%
      mutate(
        across(c(Mean, SD, Coverage), ~ round(.x, 3)),
        Scenario_label = paste(Basis, k, n, Cens, Dist, Scenario, sep = "_")
      )
    
    # Create color mapping for unique scenarios
    unique_scenarios_gp <- unique(gp_tab_with_scenario$Scenario_label)
    scenario_colors_gp <- rep(c("#F0F0F0", "#FFFFFF", "#E8F4F8", "#FFF9E6"), 
                              length.out = length(unique_scenarios_gp))
    names(scenario_colors_gp) <- unique_scenarios_gp
    
    # Build flextable with scenario coloring
    gp_ft <- gp_tab_with_scenario %>%
      flextable::flextable() %>%
      flextable::autofit()
    
    # Apply background colors row by row
    for (i in seq_len(nrow(gp_tab_with_scenario))) {
      scenario_id <- gp_tab_with_scenario$Scenario_label[i]
      gp_ft <- flextable::bg(gp_ft, i = i, bg = scenario_colors_gp[scenario_id])
    }
    
    # Add borders between different scenarios
    for (i in seq(2, nrow(gp_tab_with_scenario))) {
      if (gp_tab_with_scenario$Scenario_label[i] != gp_tab_with_scenario$Scenario_label[i-1]) {
        gp_ft <- flextable::border(gp_ft, i = i, 
                                   border.top = officer::fp_border(color = "gray30", width = 2))
      }
    }
    
    # Remove the Scenario_label helper column
    gp_ft <- flextable::delete_columns(gp_ft, j = "Scenario_label")
    
    message("\n=== GP HYPERPARAMETERS SUMMARY ===")
    print(gp_ft)
  }
  
} # end if (nrow(all_summary_tables) > 0)

# ---- Combine sAFT beta_surv across configurations
all_beta_saft_raw <- tibble::tibble()
all_beta_saft_summary <- tibble::tibble()
aft_tab <- tibble::tibble()

if (has_aft_only_models) {
  all_beta_saft_raw <- purrr::map2(
    results_by_config,
    names(results_by_config),
    ~ {
      if (is.null(.x$beta_saft_raw) || nrow(.x$beta_saft_raw) == 0) {
        return(NULL)
      }
      dplyr::mutate(.x$beta_saft_raw, config = .y)
    }
  ) %>%
    purrr::compact() %>%
    dplyr::bind_rows()

  if (nrow(all_beta_saft_raw) > 0) {
    all_beta_saft_raw <- all_beta_saft_raw %>%
      dplyr::left_join(config_info, by = "config")

    readr::write_csv(
      all_beta_saft_raw,
      file.path(summary_dir, paste0("alternative_aft_beta_surv_per_sim_", analysis_dist_file, ".csv"))
    )
  }

  all_beta_saft_summary <- purrr::map2(
    results_by_config,
    names(results_by_config),
    ~ {
      if (is.null(.x$beta_saft_summary) || nrow(.x$beta_saft_summary) == 0) {
        return(NULL)
      }
      mutate(.x$beta_saft_summary, config = .y)
    }
  ) %>%
    purrr::compact() %>%
    dplyr::bind_rows()
  
  if (nrow(all_beta_saft_summary) > 0) {
    all_beta_saft_summary <- all_beta_saft_summary %>%
      dplyr::left_join(config_info, by = "config")
  
    # produce a publication-ready table of sAFT beta_surv estimates by scenario
    aft_tab <- all_beta_saft_summary %>%
      mutate(
        Basis    = bf,
        k        = as.integer(k_bases),
        n        = as.integer(n),
        Cens     = format_cens_table(cens),
        Dist     = dist,
        Scenario = as.integer(scenario),
        Method   = method
      ) %>%
      transmute(
        Parameter = parameter,
        Basis,
        k,
        n,
        Cens,
        Dist,
        Scenario,
        Method,
        Truth = c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario],
        Mean,
        Median,
        Q2.5,
        Q97.5,
        Bias,
        SD,
        RMSE,
        n_sims
      ) %>%
      arrange(Parameter, Basis, k, n, Cens, Dist, Scenario, Method)

    readr::write_csv(
      aft_tab,
      file.path(summary_dir, paste0("alternative_aft_beta_surv_summary_", analysis_dist_file, ".csv"))
    )
    saveRDS(
      aft_tab,
      file.path(summary_dir, paste0("alternative_aft_beta_surv_summary_", analysis_dist_file, ".rds"))
    )
  
    aft_ft <- aft_tab %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      flextable::flextable() %>%
      flextable::autofit()
  
    print(aft_ft)
    message(
      "Saved alternative AFT beta_surv summaries to: ",
      file.path(summary_dir, paste0("alternative_aft_beta_surv_summary_", analysis_dist_file, ".csv"))
    )
  } else {
    message("No beta_sAFT summaries found for selected AFT-only models.")
  }
}

##########


# Data preparation for visualization
##########

### Build plotting datasets

# All Stan posterior means with short labels
all_filtered_data <- purrr::map2(
  results_by_config, names(results_by_config),
  ~ dplyr::mutate(.x$filtered_data, config = shorten_config(.y))
) %>%
  purrr::list_rbind()

if (nrow(all_filtered_data) > 0) {
  
  desired_order <- tibble::tibble(config = unique(all_filtered_data$config)) %>%
    dplyr::mutate(parsed = purrr::map(config, parse_short_config)) %>%
    tidyr::unnest(parsed) %>%
    dplyr::arrange(dist, cens, scenario, k_bases, bf) %>%
    dplyr::pull(config)
  
  all_filtered_data$config <- factor(all_filtered_data$config, levels = desired_order)

  saveRDS(
    all_filtered_data,
    file.path(summary_dir, paste0("jm_posterior_summaries_per_sim_", analysis_dist_file, ".rds"))
  )
  
  # Reference lines for key parameters
  # Note: beta_long[3] and beta_surv[1] will be added dynamically per scenario in plots
  if (all_models_are_aft_only) {
    reference_lines_plot <- tibble::tibble(
      variable = c("beta_long[1]","beta_long[2]",
                   "sd_1_long[1]","sd_1_long[2]","sigma_long",
                   "alpha"),
      reference = c(0, 0, 0, 0, 0, 0)
    )
  } else {
    reference_lines_plot <- tibble::tibble(
      variable = c("beta_long[1]","beta_long[2]",
                   "sd_1_long[1]","sd_1_long[2]","sigma_long",
                   "alpha"),
      reference = c(73, -0.04, 15, 0.20, 12, 0.012)
    )
  }
  
  # Scenario-specific reference values for faceting
  if (all_models_are_aft_only) {
    scenario_refs <- tibble::tibble(
      variable = rep(c("beta_long[3]", "beta_surv[1]"), each = 7),
      scenario = rep(1:7, 2),
      reference = c(0, 0, 0, 0, 0,0, 0,  # beta_long[3] for scenarios 1-7
                    0.00, 0.90, 0.90, -0.90, -0.90, 2.70, -2.70)   # beta_surv[1] for scenarios 1-7
    )
  } else {
    scenario_refs <- tibble::tibble(
      variable = rep(c("beta_long[3]", "beta_surv[1]"), each = 7),
      scenario = rep(1:7, 2),
      # reference = c(0.00, 0.04, -0.04, 0.08, -0.08, 0.12, -0.12,  # beta_long[3] for scenarios 1-7
      #               0.00, 0.90, -0.90, 1.80, -1.80, 2.70, -2.70)   # beta_surv[1] for scenarios 1-7
      reference = c(0.00, 0.04, -0.04, 0.04, -0.04, 0.12, -0.12,  # beta_long[3] for scenarios 1-7
                    0.00, 0.90, 0.90, -0.90, -0.90, 2.70, -2.70)   # beta_surv[1] for scenarios 1-7
    )
  }
  
  # The combined posterior-mean figure uses one facet per parameter and all
  # scenarios on the x-axis. For scenario-varying parameters, draw all distinct
  # true values observed in the plotted configs as red dashed horizontal lines.
  plotted_scenarios <- tibble::tibble(config = unique(all_filtered_data$config)) %>%
    dplyr::mutate(parsed = purrr::map(config, parse_short_config)) %>%
    tidyr::unnest(parsed) %>%
    dplyr::filter(!is.na(scenario)) %>%
    dplyr::distinct(scenario) %>%
    dplyr::pull(scenario)
  
  scenario_reference_lines_plot <- scenario_refs %>%
    dplyr::filter(scenario %in% plotted_scenarios) %>%
    dplyr::distinct(variable, reference)
  
  # Shared palette - for BP/GB combinations
  # Generate 84 distinct colors (7 scenarios × 2 censoring × 2 basis × 3 k values)
  custom_palette <- grDevices::hcl.colors(max(length(desired_order), 1), palette = "Dynamic")
  names(custom_palette) <- desired_order
  
  ##########
  
  
  # Figure generation
  ##########
  
  ### Stan posterior means
  
  # Desired facet order
  var_order <- reference_lines$variable
  
  
  #### Distribution-specific figures ####
  
  message("Generating ", analysis_dist_title, " figures for ENZAMET scenarios...")
  message("Total rows in all_filtered_data: ", nrow(all_filtered_data))
  message("Unique configs:\n", paste(unique(all_filtered_data$config), collapse = "\n"))
  
  fig_1 <- all_filtered_data %>%
    dplyr::filter(!grepl("gamma\\[", variable), !grepl("f_gp\\[", variable), !grepl("^gp_", variable)) %>%
    dplyr::filter(variable %in% var_order) %>%  # Only include variables in reference lines
    dplyr::mutate(variable = factor(variable, levels = var_order)) %>%  # enforce order
    ggplot2::ggplot(ggplot2::aes(x = config, y = mean, fill = config)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::facet_wrap(
      ~ variable,
      scales = "free_y",
      drop = FALSE   # keep empty facets if any levels missing
    ) +
    ggplot2::scale_fill_manual(values = custom_palette) +
    # Add reference lines for non-scenario-specific parameters
    ggplot2::geom_hline(
      data = reference_lines_plot %>%
        dplyr::mutate(variable = factor(variable, levels = var_order)),
      ggplot2::aes(yintercept = reference),
      color = "red", linetype = "dashed", linewidth = 0.6
    ) +
    # Add multiple reference lines for scenario-specific true values:
    # beta_long[3] and beta_surv[1] vary by scenario, so a single line would
    # be misleading in the combined ENZAMET figure.
    ggplot2::geom_hline(
      data = scenario_reference_lines_plot %>%
        dplyr::mutate(variable = factor(variable, levels = var_order)),
      ggplot2::aes(yintercept = reference),
      color = "red", linetype = "dashed", linewidth = 0.6
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position.inside = c(0.95, 0.05),
      legend.justification = c("right", "bottom")
    ) +
    ggplot2::labs(
      x = "Configuration",
      y = "Posterior Mean",
      title = paste0("Posterior Means - ENZAMET (", analysis_dist_title, ", 5 Scenarios, 0%/50%/Admin Censoring)")
    )
  
  print(fig_1)
  
  # Save figures to figures directory
  ##########
  
  message("\n=== SAVING FIGURES ===")
  
  # Create figures directory if it doesn't exist
  fig_dir <- if (bunya==TRUE) {
    here::here("JoMoNoPH-share", "Bunya", "ENZAMET", "figures")
  } else {
    here::here("Bunya", "ENZAMET", "figures")
  }
  
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
    message("Created directory: ", fig_dir)
  }
  
  # Save ENZAMET posterior mean plot
  fig_filename <- paste0("fig_ENZAMET_", analysis_dist_file, ".pdf")
  ggplot2::ggsave(filename = file.path(fig_dir, fig_filename),
                  plot = fig_1, width = 18, height = 10, units = "in")
  
  message("✓ Saved ENZAMET figure to: ", fig_dir, "\n  - ", fig_filename)
} else {
  message("No data available for figure generation")
}

##########

#### Baseline hazard plots ####
##########
# Cumulative baseline hazard plots are AFT-only at this stage.
# They require the compact baseline_recon object saved by ENZAMET-DM.R.
# The previous gamma-only reconstruction is intentionally not used because, in
# the AFT model, tau_aft depends on beta_surv and changes across posterior
# draws/replicates.

if (length(results_by_config) == 0) {
  message("Skipping baseline hazard plots because no configuration metadata was built.")
} else {
# Select configurations to plot - dynamically discover from available files
# Get all unique configs from output directory
all_rds_files <- list.files(out_dir, pattern = "\\.[rR][dD][sS]$", full.names = FALSE)
all_configs <- substr(tools::file_path_sans_ext(all_rds_files),
                      1, nchar(tools::file_path_sans_ext(all_rds_files)) - 9) %>%
  unique()

# Filter for selected models with k_bases to plot baseline hazards
# ENZAMET: n=1100, numeric/admin censoring labels, Weibull/loglogistic, FU72, scenarios 1-5 ONLY
# Use the models_to_plot vector defined at the top of the script
baseline_models_to_plot <- intersect(models_to_plot, baseline_supported_models)
configs_to_plot <- character(0)

if (length(baseline_models_to_plot) > 0) {
  models_pattern <- paste0("^(", paste(baseline_models_to_plot, collapse = "|"), ")_k[579]_n1100_cens(\\d+|admin)_(", distribution_pattern, ")_scen[1-5]_FU72$")
  configs_to_plot <- all_configs[grepl(models_pattern, all_configs)]
}

if (length(baseline_models_to_plot) == 0) {
  message("Skipping cumulative baseline hazard plots because selected models are not supported by the current AFT-only baseline_recon plotting section.")
} else if (length(configs_to_plot) == 0) {
  message("No matching configs found for baseline hazard plots.\n\nAvailable configs:\n", paste(head(all_configs, 10), collapse = "\n"))
  message("Looking for models:\n", paste(baseline_models_to_plot, collapse = "\n"))
}

message("\n=== BASELINE HAZARD ESTIMATION ===")
message("Generating AFT-only cumulative baseline hazard plots from saved baseline_recon objects...")
message("Focusing on scenarios 1-5 only")

# List to store all cumulative hazard plots for compilation
all_cumhaz_plots <- list()
baseline_recon_status <- tibble::tibble()

for (config in configs_to_plot) {
  
  message("\nProcessing configuration: ", config)
  
  # Parse config to extract k_bases, scenario, and basis type
  config_parsed <- parse_config(config)
  if (is.null(config_parsed)) {
    message("Could not parse configuration: ", config)
    next
  }
  k_bases <- config_parsed$k_bases
  scenario <- config_parsed$scenario
  bf_type <- config_parsed$bf  # "GB", "BP", "GB_Quantile", or "GP"
  message("  k_bases: ", k_bases)
  message("  Scenario: ", scenario)
  message("  Model type: ", bf_type)
  
  # Find matching files
  pattern <- paste0(config, "(?:_.*)?_batch\\d{3}\\.rds$")
  rds_files <- list.files(path = out_dir, pattern = pattern, full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("No files found for baseline hazard plot: ", config)
    next
  }
  
  # Read all batches so the cumulative baseline hazard summary uses all
  # available simulation replicates for this configuration.
  sims <- purrr::compact(unlist(lapply(rds_files, readRDS), recursive = FALSE))
  n_sims_to_plot <- min(5, length(sims))
  
  if (n_sims_to_plot == 0) {
    message("No simulations in file: ", rds_files[1])
    next
  }
  
  # Extract model type from summaries
  model_type <- sims[[1]]$joint_summary$model[1]
  message("  Model type from data: ", model_type)
  message("  Found ", length(sims), " simulations; averaging all saved baseline curves and showing up to ", n_sims_to_plot, " replicate curves")

  n_with_baseline_recon <- sum(vapply(
    sims,
    function(sim) !is.null(sim$baseline_recon) && !is.null(sim$baseline_recon$curve),
    logical(1)
  ))

  baseline_recon_status <- dplyr::bind_rows(
    baseline_recon_status,
    tibble::tibble(
      config = config,
      n_batch_files = length(rds_files),
      n_sims = length(sims),
      n_with_baseline_recon = n_with_baseline_recon,
      can_plot_cumhaz = n_with_baseline_recon > 0
    )
  )

  if (!model_type %in% baseline_supported_models) {
    message("Skipping cumulative baseline hazard plot for unsupported AFT-only baseline_recon model type: ", model_type)
    next
  }

  basis_name <- dplyr::case_when(
    model_type == "BP_aft" ~ "BP AFT",
    model_type == "BP_aft_logy" ~ "BP log(y) AFT",
    model_type == "GB_aft" ~ "Gaussian basis AFT",
    model_type == "GB_aft_logy" ~ "Gaussian basis log(y) AFT",
    TRUE ~ model_type
  )
  
  # Source data generation file for true parameter values
  # ENZAMET data gen files are scenario-specific
  data_gen_file <- sprintf("Data-Gen-ENZAMET-%02d.R", scenario)
  data_gen_path <- if (bunya==TRUE) {
    here::here("JoMoNoPH-share", "Bunya", "ENZAMET", data_gen_file)
  } else {
    here::here("Bunya", "ENZAMET", data_gen_file)
  }
  
  if (file.exists(data_gen_path)) {
    source(data_gen_path)
  } else {
    message("Could not find data generation file: ", data_gen_path)
    next
  }

  baseline_curve_data <- purrr::imap_dfr(sims, function(sim, sim_idx) {
    baseline_recon <- sim$baseline_recon
    if (is.null(baseline_recon) || is.null(baseline_recon$curve)) {
      return(NULL)
    }

    sim_id <- if (!is.null(sim$joint_summary) && "sim_id" %in% names(sim$joint_summary)) {
      unique(sim$joint_summary$sim_id)[1]
    } else if (!is.null(sim$unique_seed)) {
      sim$unique_seed
    } else {
      sim_idx
    }

    baseline_recon$curve %>%
      dplyr::mutate(sim_id = factor(sim_id))
  })

  if (nrow(baseline_curve_data) == 0) {
    message("No saved baseline_recon object found for: ", config)
    message("Re-run ENZAMET-DM.R with the updated baseline reconstruction output to draw scale-correct cumulative baseline hazard plots.")
    next
  }

  estimated_cumhaz_summary <- baseline_curve_data %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      H0_est = mean(H0_mean, na.rm = TRUE),
      H0_q2.5 = stats::quantile(H0_mean, 0.025, na.rm = TRUE),
      H0_q97.5 = stats::quantile(H0_mean, 0.975, na.rm = TRUE),
      n_sims = dplyr::n(),
      .groups = "drop"
    )

  replicate_ids_to_plot <- unique(baseline_curve_data$sim_id)[
    seq_len(min(5, dplyr::n_distinct(baseline_curve_data$sim_id)))
  ]
  replicate_curves_to_plot <- baseline_curve_data %>%
    dplyr::filter(sim_id %in% replicate_ids_to_plot)

  t_grid <- sort(unique(baseline_curve_data$time))
  if (tolower(config_parsed$dist) == "loglogistic") {
    true_cumhaz <- log1p((t_grid / loglogistic_scale)^loglogistic_shape)
  } else if (tolower(config_parsed$dist) == "weibull") {
    true_cumhaz <- (t_grid / weibull_scale)^weibull_shape
  } else {
    message("Unsupported true baseline distribution for cumulative hazard plot: ", config_parsed$dist)
    next
  }

  true_cumhaz_data <- tibble::tibble(
    time = t_grid,
    H0 = true_cumhaz
  )

  p_cumhaz <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = estimated_cumhaz_summary,
      ggplot2::aes(x = time, ymin = H0_q2.5, ymax = H0_q97.5),
      fill = "#7AA6C2",
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = replicate_curves_to_plot,
      ggplot2::aes(x = time, y = H0_mean, group = sim_id),
      color = "#7A7A7A",
      linewidth = 0.45,
      linetype = "dashed",
      alpha = 0.6
    ) +
    ggplot2::geom_line(
      data = estimated_cumhaz_summary,
      ggplot2::aes(x = time, y = H0_est),
      color = "#1F78B4",
      linewidth = 1.1
    ) +
    ggplot2::geom_line(
      data = true_cumhaz_data,
      ggplot2::aes(x = time, y = H0),
      color = "black",
      linewidth = 1.1
    ) +
    ggplot2::labs(
      title = paste0("Cumulative Baseline Hazard: ENZAMET Scenario ", scenario, " (", basis_name, ", k=", k_bases, ")"),
      subtitle = paste0(
        "n=", config_parsed$n, ", ", format_cens_label(config_parsed$cens),
        ". Black = true; blue = mean over ", max(estimated_cumhaz_summary$n_sims),
        " replicates; ribbon = 2.5-97.5% across replicate posterior means."
      ),
      x = "Time t (months; baseline X = 0)",
      y = "Cumulative Baseline Hazard H_0(t)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )

  print(p_cumhaz)

  fig_name_cumhaz <- paste0("fig_cumhaz_", gsub("_FU72", "", config))
  ggplot2::ggsave(
    filename = file.path(fig_dir, paste0(fig_name_cumhaz, ".pdf")),
    plot = p_cumhaz,
    width = 10,
    height = 6,
    units = "in"
  )

  all_cumhaz_plots[[config]] <- p_cumhaz
  message("  ✓ Cumulative baseline hazard plot created and saved from saved baseline_recon: ", fig_name_cumhaz)
}

if (nrow(baseline_recon_status) > 0) {
  baseline_status_file <- file.path(summary_dir, paste0("baseline_recon_status_", analysis_dist_file, ".csv"))
  readr::write_csv(baseline_recon_status, baseline_status_file)
  message("Saved baseline reconstruction status to: ", baseline_status_file)
}

# Compile all cumulative hazard plots into a single PDF
if (length(all_cumhaz_plots) > 0) {
  message("\n=== Compiling all cumulative baseline hazard plots into single PDF ===")
  
  # Create combined plot using patchwork
  # Arrange plots in a grid (2 columns for better readability)
  combined_cumhaz <- patchwork::wrap_plots(all_cumhaz_plots, ncol = 2)
  
  # Save combined PDF
  combined_filename <- file.path(fig_dir, paste0("fig_cumhaz_", analysis_dist_file, "_ALL_scenarios_1_5.pdf"))
  ggplot2::ggsave(
    filename = combined_filename,
    plot = combined_cumhaz,
    width = 20,
    height = 6 * ceiling(length(all_cumhaz_plots) / 2),
    units = "in",
    limitsize = FALSE
  )
  
  message("  ✓ Combined cumulative baseline hazard plot saved: ", basename(combined_filename))
  message("  ✓ Total plots included: ", length(all_cumhaz_plots))
} else {
  message("\n=== No cumulative baseline hazard plots generated ===")
}

message("\n=== BASELINE HAZARD ESTIMATION COMPLETE ===\n")
}

##########

# ============================================================
# PERFORMANCE PLOTS (Bias / Abs Error / RMSE / Coverage) — scenarios 1–5
# WITH manuscript-friendly parameter labels:
#   beta_long[1] -> beta_0
#   beta_long[2] -> beta_1
#   beta_long[3] -> beta_2
#   beta_surv[1] -> gamma_1   (as per your supplied truth table)
#   alpha        -> alpha
#
# REQUIRES in environment:
#   - all_filtered_data with columns:
#       variable, mean, q2.5, q97.5, sim_id, config
#     where config is the SHORT label like: "BP k=5, c=0%, S1"
# ============================================================

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("arrange", "dplyr")
conflicted::conflict_prefer("summarise", "dplyr")

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
key_vars <- c("alpha", "beta_long[1]", "beta_long[2]", "beta_long[3]", "beta_surv[1]")

# Mapping to manuscript labels (parsed expressions for nice facet strips)
param_labels <- c(
  "alpha"        = "alpha",
  "beta_long[1]" = "beta[0]",
  "beta_long[2]" = "beta[1]",
  "beta_long[3]" = "beta[2]",
  "beta_surv[1]" = "gamma[1]"
)

# True values (scenarios 1–5) using exactly what you supplied
if (all_models_are_aft_only) {
  truth_15 <- tibble::tribble(
    ~scenario, ~variable,        ~truth,
    1,         "beta_long[1]",   0,
    2,         "beta_long[1]",   0,
    3,         "beta_long[1]",   0,
    4,         "beta_long[1]",   0,
    5,         "beta_long[1]",   0,
    
    1,         "beta_long[2]",   0,
    2,         "beta_long[2]",   0,
    3,         "beta_long[2]",   0,
    4,         "beta_long[2]",   0,
    5,         "beta_long[2]",   0,
    
    1,         "beta_long[3]",   0.00,
    2,         "beta_long[3]",   0.04,
    3,         "beta_long[3]",   -0.04,
    4,         "beta_long[3]",   0.04,
    5,         "beta_long[3]",   -0.04,
    
    1,         "alpha",          0.000,
    2,         "alpha",          0.000,
    3,         "alpha",          0.000,
    4,         "alpha",          0.000,
    5,         "alpha",          0.000,
    
    # you used gamma_1 as the survival covariate effect
    1,         "beta_surv[1]",   0.00,
    2,         "beta_surv[1]",   0.90,
    3,         "beta_surv[1]",   0.90,
    4,         "beta_surv[1]",   -0.90,
    5,         "beta_surv[1]",   -0.90
  )
} else {
  truth_15 <- tibble::tribble(
    ~scenario, ~variable,        ~truth,
    1,         "beta_long[1]",   73,
    2,         "beta_long[1]",   73,
    3,         "beta_long[1]",   73,
    4,         "beta_long[1]",   73,
    5,         "beta_long[1]",   73,
    
    1,         "beta_long[2]",   -0.04,
    2,         "beta_long[2]",   -0.04,
    3,         "beta_long[2]",   -0.04,
    4,         "beta_long[2]",   -0.04,
    5,         "beta_long[2]",   -0.04,
    
    1,         "beta_long[3]",   0.00,
    2,         "beta_long[3]",   0.04,
    3,         "beta_long[3]",   -0.04,
    4,         "beta_long[3]",   0.04,
    5,         "beta_long[3]",   -0.04,
    
    1,         "alpha",          0.012,
    2,         "alpha",          0.012,
    3,         "alpha",          0.012,
    4,         "alpha",          0.012,
    5,         "alpha",          0.012,
    
    # you used gamma_1 as the survival covariate effect
    1,         "beta_surv[1]",   0.00,
    2,         "beta_surv[1]",   0.90,
    3,         "beta_surv[1]",   0.90,
    4,         "beta_surv[1]",   -0.90,
    5,         "beta_surv[1]",   -0.90
  )
}

config_info <- tibble(config = unique(all_filtered_data$config)) %>%
  mutate(parsed = map(config, parse_short_config)) %>%
  tidyr::unnest(parsed)

# ------------------------------------------------------------
# Build per-simulation performance table (creates distributions)
# ------------------------------------------------------------
perf_sim <- all_filtered_data %>%
  dplyr::filter(variable %in% key_vars) %>%
  left_join(config_info, by = "config") %>%
  dplyr::filter(!is.na(scenario), scenario %in% 1:5) %>%
  left_join(truth_15, by = c("scenario", "variable")) %>%
  dplyr::filter(!is.na(truth)) %>%
  mutate(
    err      = mean - truth,
    bias     = err,
    abs_err  = abs(err),   # distribution across sims (absolute error)
    sq_error = err^2,
    covered  = as.integer(q2.5 <= truth & truth <= q97.5),
    
    # Use manuscript-friendly labels for facets
    Parameter = factor(variable, levels = key_vars, labels = unname(param_labels[key_vars])),
    
    Dist      = factor(dist),
    Baseline  = factor(dist_short),
    Scenario  = factor(scenario),
    Cens      = factor(format_cens_table(cens)),
    k         = factor(k_bases)
  )

# Sanity checks
counts <- perf_sim %>%
  count(Parameter, Dist, Baseline, Scenario, Cens, k) %>%
  arrange(Parameter, Dist, Baseline, Cens, Scenario, k)
print(counts)

if (nrow(perf_sim) == 0) {
  stop("perf_sim is empty. Check all_filtered_data$config format, key_vars names, and that scenarios 1–5 exist.")
}

# Optional diagnostic: check beta_0 (intercept) centering
beta0_check <- perf_sim %>%
  dplyr::filter(as.character(Parameter) == "beta[0]") %>%
  group_by(Dist, Baseline, Cens, Scenario, k) %>%
  summarise(
    mean_of_means   = mean(mean),
    sd_of_means     = sd(mean),
    median_abs_err  = median(abs_err),
    .groups = "drop"
  ) %>%
  arrange(Dist, Baseline, Cens, Scenario, k)
print(beta0_check)

# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------

facet_labeller <- ggplot2::labeller(
  Parameter = label_parsed
)

# Bias distribution
p_bias <- ggplot(perf_sim, aes(x = Scenario, y = bias)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Cens + Baseline ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Bias (posterior mean − truth)",
    title = paste0("Bias distributions across simulations (", analysis_dist_title, ", scenarios 1–5), by censoring")
  )

# Absolute error distribution (NOT RMSE)
p_abserr <- ggplot(perf_sim, aes(x = Scenario, y = abs_err)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_grid(Cens + Baseline ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Absolute error |posterior mean − truth|",
    title = paste0("Absolute error distributions across simulations (", analysis_dist_title, ", scenarios 1–5), by censoring")
  )

# True RMSE is a group summary (one value per Scenario × Cens × Parameter × k)
rmse_sum <- perf_sim %>%
  group_by(Parameter, Baseline, Scenario, Cens, k) %>%
  summarise(
    rmse = sqrt(mean(sq_error)),
    .groups = "drop"
  )

p_rmse <- ggplot(rmse_sum, aes(x = Scenario, y = rmse, group = k)) +
  geom_point(position = position_dodge(width = 0.35)) +
  facet_grid(Cens + Baseline ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "RMSE across simulations",
    title = paste0("RMSE (summary) by scenario (", analysis_dist_title, ", scenarios 1–5), by censoring")
  )

# Coverage summary + CI (binomial normal approx)
cov_sum <- perf_sim %>%
  group_by(Parameter, Baseline, Scenario, Cens, k) %>%
  summarise(
    coverage = mean(covered),
    n = n(),
    se = sqrt(coverage * (1 - coverage) / n), # ?
    lo = pmax(0, coverage - 1.96 * se), # ?
    hi = pmin(1, coverage + 1.96 * se), # ?
    .groups = "drop"
  )

p_cov <- ggplot(cov_sum, aes(x = Scenario, y = coverage, group = k)) +
  geom_point(position = position_dodge(width = 0.35)) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15,
    position = position_dodge(width = 0.35)
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_grid(Cens + Baseline ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Coverage (mean ± 1.96×SE)",
    title = paste0("Coverage across simulations (", analysis_dist_title, ", scenarios 1–5), by censoring")
  )

# ------------------------------------------------------------
# SAVE as PNGs
# ------------------------------------------------------------
out_dir <- here("Bunya/ENZAMET/figures/performance")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, paste0("bias_distributions_", analysis_dist_file, "_scen1_5.png")),     p_bias,   width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, paste0("abs_error_distributions_", analysis_dist_file, "_scen1_5.png")), p_abserr, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, paste0("rmse_summary_", analysis_dist_file, "_scen1_5.png")),           p_rmse,   width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, paste0("coverage_", analysis_dist_file, "_scen1_5.png")),               p_cov,    width = 14, height = 8, dpi = 300)

message("Saved plots to: ", out_dir)

# ============================================================
# The short config labels now include d=W/d=LL so Weibull and log-logistic
# summaries are not accidentally merged when both are analysed together.
# ============================================================
