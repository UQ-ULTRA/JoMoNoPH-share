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
# 6. Figure generation (loglogistic hazards)
#
# OUTPUTS:
# - Combined summary tables (Stan estimates by scenario)
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
models_to_plot <- c("BP_aft")
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

if (bunya==TRUE) {
  out_dir <- here::here("JoMoNoPH-share", "Bunya", "ENZAMET", "output_20260428_24052266")
} else {
  out_dir <- here::here("Bunya", "ENZAMET", "output_20260428_24052266")
}

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

# --- Helper: make short, readable labels for facets/plots
shorten_config <- function(cfg) {
  # Updated pattern to handle BP, GB, GB_Quantile, GP and scenario
  parts <- stringr::str_match(
    cfg,
    "^([A-Za-z_]+)_k(\\d+)_n(\\d+)_cens([A-Za-z0-9]+)_([A-Za-z]+)_scen(\\d+)_?(FU\\d+)?$"
  )
  if (is.na(parts[1, 1])) return(cfg)
  
  bf   <- parts[2]
  k    <- parts[3]
  n    <- parts[4]
  cens <- parts[5]
  dist_letter <- toupper(substr(parts[6], 1, 1))  # "L" for loglogistic
  scen <- parts[7]
  
  paste0(bf, " k=", k, ", c=", format_cens_short(cens), ", S", scen)
}


# Function to parse a config string 
parse_config <- function(config_str) {
  # Updated pattern to handle BP, GB, GB_Quantile, GP and scenario
  parts <- stringr::str_match(config_str, "^([A-Za-z_]+)_k(\\d+)_n(\\d+)_cens([A-Za-z0-9]+)_([A-Za-z]+)_scen(\\d+)_?(FU\\d+)?$")
  
  if (is.na(parts[1, 1])) return(NULL)
  
  list(
    bf = parts[2],
    k_bases = as.numeric(parts[3]),
    n = as.numeric(parts[4]),
    cens = parts[5],
    dist = parts[6],
    scenario = as.numeric(parts[7]),
    fu = parts[8]
  )
}

parse_short_config <- function(lbl) {
  m <- stringr::str_match(lbl, "^([A-Za-z_]+)\\s+k=(\\d+),\\s+c=([A-Za-z0-9]+)%?,\\s+S(\\d+)$")
  if (is.na(m[1,1])) {
    return(tibble::tibble(
      bf = NA_character_,
      k_bases = NA_real_,
      cens = NA_character_,
      scenario = NA_real_
    ))
  }
  tibble::tibble(
    bf       = m[1,2],
    k_bases  = as.numeric(m[1,3]),
    cens     = m[1,4],
    scenario = as.numeric(m[1,5])
  )
}

##########


# Main processing function
##########

# This function reads, processes, and aggregates simulation results
# from Bayesian joint modelling runs

# --- Main: process all sims for a single config
process_simulation_results <- function(bf, k_bases, n, cens, dist, scenario, fu,
                                       output_dir) {
  
  # Build filename pattern for this config
  fu_part <- if (!is.null(fu) && !is.na(fu) && nzchar(fu)) {
    paste0("_", fu)
  } else {
    ""
  }
  
  pattern <- paste0(
    bf, "_k", k_bases, "_n", n, "_cens", cens, "_", dist,
    "_scen", scenario,
    fu_part,
    "(?:_.*)?_batch\\d{3}\\.rds$"
  )
  
  
  rds_files <- list.files(path = output_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  
  if (length(rds_files) == 0) {
    message("No matching files found for: ", pattern)
    return(NULL)
  }
  # browser()
  # ---- Read all batches; flatten to one list of sims; drop NULLs (timeouts/failures)
  sims <- rds_files %>%
    purrr::map(readRDS) %>%
    purrr::flatten() %>%
    purrr::compact()
  
  if (length(sims) == 0) {
    message("All simulations were NULL/failed for: ", pattern)
    return(NULL)
  }
  
  # ---- Aggregate Stan summaries (keep parameter names + sim id)
  all_summaries_list <- sims %>%
    purrr::imap(function(sim, idx) {
      s <- sim$joint_summary
      
      # Skip if summary is NULL or empty
      if (is.null(s)) return(NULL)
      
      s_df <- as.data.frame(s)
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
          SD       = sd(estimate, na.rm = TRUE),
          Bias     = mean(bias, na.rm = TRUE),
          RMSE     = sqrt(mean(bias^2, na.rm = TRUE)),
          n_sims   = n(),
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
      
      n_sims = n()
    )
  
  
  # Filter parameters of interest (including GP-specific parameters)
  filtered_data <- all_summaries %>%
    dplyr::filter(stringr::str_detect(variable, "rho|alpha|beta_|sigma_long|sd_1_long|gamma|sigma|gp_length_scale|gp_marginal_sd|f_gp"))
  
  message("Total simulations found for ", bf, "_n", n, "_cens", cens, "_", dist, "_scen", scenario, "_", fu,
          ": N=", length(unique(filtered_data$sim_id)))
  
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
  stopifnot(scenario %in% 1:5)
  
  if (all_models_are_aft_only) {
    true_values <- c(
      "beta_long[1]"  = 0,     # beta_0
      "beta_long[2]"  = 0,  # beta_1
      "beta_long[3]"  = c(0, 0, 0, 0, 0)[scenario],  # beta_2 (scen 1–5)
      
      "sd_1_long[1]"  = 0,     # sigma_b0
      "sd_1_long[2]"  = 0,    # sigma_b1
      "sigma_long"    = 0,     # sigma_epsilon
      
      # "alpha"         = 0.012,  # alpha (scen 1–5)
      "alpha"         = 0,      # try no association 20260113
      "beta_surv[1]"  = c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario]   # gamma_1 (scen 1–5)
    )
  } else {
    true_values <- c(
      "beta_long[1]"  = 73,     # beta_0
      "beta_long[2]"  = -0.04,  # beta_1
      "beta_long[3]"  = c(0.00, 0.04, -0.04, 0.04, -0.04)[scenario],  # beta_2 (scen 1–5)
      
      "sd_1_long[1]"  = 15,     # sigma_b0
      "sd_1_long[2]"  = 0.2,    # sigma_b1
      "sigma_long"    = 12,     # sigma_epsilon
      
      # "alpha"         = 0.012,  # alpha (scen 1–5)
      "alpha"         = 0,      # try no association 20260113
      "beta_surv[1]"  = c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario]   # gamma_1 (scen 1–5)
    ) 
  }
  
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
  
  # Posterior summary across sims
  post_summary <- filtered_data %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      mean      = mean(mean, na.rm = TRUE),
      median    = median(mean, na.rm = TRUE),
      sd        = mean(sd, na.rm = TRUE),  # Average posterior SD across sims
      .groups = "drop"
    ) %>%
    dplyr::mutate(source = "Stan posterior") %>%
    dplyr::left_join(coverage_data, by = "variable")
  
  # ENZAMET files don't have LMM - skip LMM processing
  
  # Combined summary is just Stan posterior
  combined_summary <- post_summary %>%
    dplyr::select(variable, mean, sd, source, coverage, n_sims) %>%
    dplyr::arrange(source, variable)
  
  # Flextable: rounded values + simple header
  table <- combined_summary %>%
    dplyr::mutate(
      mean = round(mean, 3),
      sd   = round(sd, 3),
      coverage = round(coverage, 3)
    ) %>%
    flextable::flextable() %>%
    flextable::add_header_row(
      values   = c("Variable", "Mean", "SD", "Source", "Coverage", "N Sims"),
      colwidths = c(1, 1, 1, 1, 1, 1)
    ) %>%
    flextable::align(j = "source", align = "left")
  
  list(
    fig           = fig_hist,        # diagnostic histograms
    table         = table,           # Stan summary table
    filtered_data = filtered_data,   # long Stan posterior means per sim
    summary       = combined_summary, # aggregated table (Stan only)
    diag_per_sim  = diag_per_sim,    # NEW
    diag_summary  = diag_summary,    # NEW
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
        dist    = cfg_parsed$dist,
        scenario = cfg_parsed$scenario,
        fu      = cfg_parsed$fu,
        output_dir = out_dir
      )
    })
)


# ---- Combine diagnostic summaries across configs
all_diag_summary <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$diag_summary, config = .y)
)
all_diag_summary
# diagnostic plot: boxplot of max R-hat per config
# long form, useful for boxplots of max_rhat or realised censoring
all_diag_per_sim <- purrr::map2_dfr(
  results_by_config,
  names(results_by_config),
  ~ dplyr::mutate(.x$diag_per_sim, config = .y)
)

shorten_config_1 <- function(config) {
  
  tibble::tibble(config = config) %>%
    dplyr::mutate(
      model = stringr::str_replace(config, "_k\\d+.*$", ""),
      k     = stringr::str_extract(config, "k\\d+"),
      n     = stringr::str_extract(config, "n\\d+"),
      dist  = stringr::str_extract(config, "Weibull|loglogistic"),
      FU    = stringr::str_extract(config, "FU\\d+")
    ) %>%
    dplyr::mutate(
      label = paste(model, k, n, dist, FU, sep = ", ")
    ) %>%
    dplyr::pull(label)
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

# Wider view (mean/sd columns per config)
if (nrow(all_summary_tables) > 0 && any(c("mean", "sd") %in% colnames(all_summary_tables))) {
  summary_wide <- all_summary_tables %>%
    # Filter to only include BP and loglogistic
    dplyr::filter(stringr::str_detect(config, "^BP_") & stringr::str_detect(config, "_loglogistic_")) %>%
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
    # reference = c(73, -0.04, NA, 15, 0.20, 12, NA, 0.012)
    reference = c(73, -0.04, NA, 15, 0.20, 12, NA, 0) # try no association 20260113
  )
}

key_vars <- c("alpha", "beta_long[1]", "beta_long[2]", "beta_long[3]", "beta_surv[1]")

if (nrow(all_summary_tables) > 0) {
  # Attach true values
  # First get config info to determine scenario
  config_info <- tibble::tibble(config = names(results_by_config)) %>%
    dplyr::mutate(parsed = purrr::map(config, parse_config)) %>%
    tidyr::unnest_wider(parsed)   # bf, n, cens, dist, scenario, fu
  
  results_tab <- all_summary_tables %>%
    dplyr::left_join(config_info, by = "config") %>%
    left_join(reference_lines, by = "variable") %>%
    mutate(
      # Update scenario-specific true values
      reference = case_when(
        variable == "beta_long[3]" & scenario == 1 ~ 0.00,
        variable == "beta_long[3]" & scenario == 2 ~ 0.04,
        variable == "beta_long[3]" & scenario == 3 ~ -0.04,
        # variable == "beta_long[3]" & scenario == 4 ~ 0.08,
        variable == "beta_long[3]" & scenario == 4 ~ 0.04,
        # variable == "beta_long[3]" & scenario == 5 ~ -0.08,
        variable == "beta_long[3]" & scenario == 5 ~ -0.04,
        variable == "beta_long[3]" & scenario == 6 ~ 0.12,
        variable == "beta_long[3]" & scenario == 7 ~ -0.12,
        variable == "beta_surv[1]" & scenario == 1 ~ 0.00,
        variable == "beta_surv[1]" & scenario == 2 ~ 0.90,
        # variable == "beta_surv[1]" & scenario == 3 ~ -0.90,
        variable == "beta_surv[1]" & scenario == 3 ~ 0.90,
        # variable == "beta_surv[1]" & scenario == 4 ~ 1.80,
        variable == "beta_surv[1]" & scenario == 4 ~ -0.90,
        # variable == "beta_surv[1]" & scenario == 5 ~ -1.80,
        variable == "beta_surv[1]" & scenario == 5 ~ -0.90,
        variable == "beta_surv[1]" & scenario == 6 ~ 2.70,
        variable == "beta_surv[1]" & scenario == 7 ~ -2.70,
        TRUE ~ reference
      ),
      bias      = mean - reference,
      rel_bias  = if_else(reference != 0, bias / reference, NA_real_),
      rmse      = sqrt(sd^2 + bias^2)   # RMSE^2 = Var + Bias^2
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
      Scenario,
      Method,
      Truth    = reference,
      Mean     = mean,
      SD       = sd,
      RelBias  = rel_bias,
      RMSE     = rmse,
      Coverage = coverage
    ) %>%
    arrange(Parameter, Basis, k, n, Cens, Dist, Scenario, Method)
  
  nice_ft <- nice_tab_wide %>%
    mutate(across(c(Truth, Mean, SD, RelBias, RMSE, Coverage), ~ round(.x, 3))) %>%
    flextable::flextable() %>%
    flextable::autofit()
  
  print(nice_ft)
  
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
      Scenario  = as.integer(scenario)
    ) %>%
    transmute(
      Parameter,
      Basis,
      k,
      n,
      Cens,
      Dist,
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
      Scenario  = as.integer(scenario)
    ) %>%
    transmute(
      Parameter,
      Basis,
      k,
      n,
      Cens,
      Dist,
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
if (has_aft_only_models) {
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
      left_join(config_info, by = "config")
  
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
        Bias,
        SD,
        RMSE,
        n_sims
      ) %>%
      arrange(Parameter, Basis, k, n, Cens, Dist, Scenario, Method)
  
    aft_ft <- aft_tab %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      flextable::flextable() %>%
      flextable::autofit()
  
    print(aft_ft)
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
    dplyr::arrange(cens, scenario, k_bases, bf) %>%
    dplyr::pull(config)
  
  all_filtered_data$config <- factor(all_filtered_data$config, levels = desired_order)
  
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
      # reference = c(73, -0.04, 15, 0.20, 12, 0.012)
      reference = c(73, -0.04, 15, 0.20, 12, 0) # try no association 20260113
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
  
  
  #### Loglogistic distribution figures (ENZAMET only uses loglogistic) ####
  
  message("Generating loglogistic figures for ENZAMET scenarios...")
  message("Total rows in all_filtered_data: ", nrow(all_filtered_data))
  message("Unique configs:\n", paste(unique(all_filtered_data$config), collapse = "\n"))
  
  fig_1L <- all_filtered_data %>%
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
      # title = "Posterior Means - ENZAMET (Loglogistic Hazard, 7 Scenarios, 0%/50% Censoring)"
      title = "Posterior Means - ENZAMET (5 Scenarios, 0%/50% Censoring)"
    )
  
  print(fig_1L)
  
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
  
  # Save ENZAMET loglogistic plot
  ggplot2::ggsave(filename = file.path(fig_dir, "fig_ENZAMET_loglogistic.pdf"),
                  plot = fig_1L, width = 18, height = 10, units = "in")
  
  message("✓ Saved ENZAMET figure to: ", fig_dir, "\n  - fig_ENZAMET_loglogistic.pdf")
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
# ENZAMET: n=1100, numeric/admin censoring labels, loglogistic, FU72, scenarios 1-5 ONLY
# Use the models_to_plot vector defined at the top of the script
baseline_models_to_plot <- intersect(models_to_plot, baseline_supported_models)
configs_to_plot <- character(0)

if (length(baseline_models_to_plot) > 0) {
  models_pattern <- paste0("^(", paste(baseline_models_to_plot, collapse = "|"), ")_k[579]_n1100_cens(\\d+|admin)_loglogistic_scen[1-5]_FU72$")
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

# Compile all cumulative hazard plots into a single PDF
if (length(all_cumhaz_plots) > 0) {
  message("\n=== Compiling all cumulative baseline hazard plots into single PDF ===")
  
  # Create combined plot using patchwork
  # Arrange plots in a grid (2 columns for better readability)
  combined_cumhaz <- patchwork::wrap_plots(all_cumhaz_plots, ncol = 2)
  
  # Save combined PDF
  combined_filename <- file.path(fig_dir, "fig_cumhaz_ALL_scenarios_1_5.pdf")
  ggplot2::ggsave(
    filename = combined_filename,
    plot = combined_cumhaz,
    width = 20,
    height = 6 * ceiling(length(all_cumhaz_plots) / 2),
    units = "in",
    limitsize = FALSE
  )
  
  message("  ✓ Combined cumulative baseline hazard plot saved: fig_cumhaz_ALL_scenarios_1_5.pdf")
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
    
    # 1,         "alpha",          0.012,
    # 2,         "alpha",          0.012,
    # 3,         "alpha",          0.012,
    # 4,         "alpha",          0.012,
    # 5,         "alpha",          0.012,
    1,         "alpha",          0.000, # try no association 20260113
    2,         "alpha",          0.000, # try no association 20260113
    3,         "alpha",          0.000, # try no association 20260113
    4,         "alpha",          0.000, # try no association 20260113
    5,         "alpha",          0.000, # try no association 20260113
    
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
    
    # 1,         "alpha",          0.012,
    # 2,         "alpha",          0.012,
    # 3,         "alpha",          0.012,
    # 4,         "alpha",          0.012,
    # 5,         "alpha",          0.012,
    1,         "alpha",          0.000, # try no association 20260113
    2,         "alpha",          0.000, # try no association 20260113
    3,         "alpha",          0.000, # try no association 20260113
    4,         "alpha",          0.000, # try no association 20260113
    5,         "alpha",          0.000, # try no association 20260113
    
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
    
    Scenario  = factor(scenario),
    Cens      = factor(format_cens_table(cens)),
    k         = factor(k_bases)
  )

# Sanity checks
counts <- perf_sim %>%
  count(Parameter, Scenario, Cens, k) %>%
  arrange(Parameter, Cens, Scenario, k)
print(counts)

if (nrow(perf_sim) == 0) {
  stop("perf_sim is empty. Check all_filtered_data$config format, key_vars names, and that scenarios 1–5 exist.")
}

# Optional diagnostic: check beta_0 (intercept) centering
beta0_check <- perf_sim %>%
  dplyr::filter(as.character(Parameter) == "beta[0]") %>%
  group_by(Cens, Scenario, k) %>%
  summarise(
    mean_of_means   = mean(mean),
    sd_of_means     = sd(mean),
    median_abs_err  = median(abs_err),
    .groups = "drop"
  ) %>%
  arrange(Cens, Scenario, k)
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
  facet_grid(Cens ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Bias (posterior mean − truth)",
    title = "Bias distributions across simulations (scenarios 1–5), by censoring"
  )

# Absolute error distribution (NOT RMSE)
p_abserr <- ggplot(perf_sim, aes(x = Scenario, y = abs_err)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_grid(Cens ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Absolute error |posterior mean − truth|",
    title = "Absolute error distributions across simulations (scenarios 1–5), by censoring"
  )

# True RMSE is a group summary (one value per Scenario × Cens × Parameter × k)
rmse_sum <- perf_sim %>%
  group_by(Parameter, Scenario, Cens, k) %>%
  summarise(
    rmse = sqrt(mean(sq_error)),
    .groups = "drop"
  )

p_rmse <- ggplot(rmse_sum, aes(x = Scenario, y = rmse, group = k)) +
  geom_point(position = position_dodge(width = 0.35)) +
  facet_grid(Cens ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "RMSE across simulations",
    title = "RMSE (summary) by scenario (scenarios 1–5), by censoring"
  )

# Coverage summary + CI (binomial normal approx)
cov_sum <- perf_sim %>%
  group_by(Parameter, Scenario, Cens, k) %>%
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
  facet_grid(Cens ~ Parameter, scales = "free_y", labeller = facet_labeller) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Coverage (mean ± 1.96×SE)",
    title = "Coverage across simulations (scenarios 1–5), by censoring"
  )

# ------------------------------------------------------------
# SAVE as PNGs
# ------------------------------------------------------------
out_dir <- here("Bunya/ENZAMET/figures/performance")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "bias_distributions_scen1_5.png"),    p_bias,   width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "abs_error_distributions_scen1_5.png"), p_abserr, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "rmse_summary_scen1_5.png"),         p_rmse,   width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "coverage_scen1_5.png"),             p_cov,    width = 14, height = 8, dpi = 300)

message("Saved plots to: ", out_dir)

# ============================================================
# OPTIONAL (RECOMMENDED) — avoid mixing distributions
# If you can rebuild all_filtered_data, keep config_full as well:
#
# all_filtered_data <- purrr::map2(
#   results_by_config, names(results_by_config),
#   ~ dplyr::mutate(
#       .x$filtered_data,
#       config_full = .y,
#       config      = shorten_config(.y)
#     )
# ) %>% purrr::list_rbind()
#
# Then parse dist from config_full using parse_config(), and filter:
#   filter(dist == "loglogistic")
# ============================================================
