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
# Options: "BP", "GB", "GB_Quantile", "GP"
# Example: models_to_plot <- c("BP", "GP")  # for BP and GP
#          models_to_plot <- c("BP", "GB")  # for BP and GB
# models_to_plot <- c("BP")  # <-- CHANGE THIS to select models
models_to_plot <- c("GB")
# =============================================================================

##########


# File operations and setup
##########

bunya <- FALSE

if (bunya==TRUE) {
  out_dir <- here::here("JoMoNoPH-share", "Bunya", "ENZAMET", "output")
} else {
  out_dir <- here::here("Bunya", "ENZAMET", "output")
}

# out_dir <- "C:/Users/uqamar43/OneDrive - The University of Queensland/02 shared HERA ULTRA/04 Stat Methods/08 JoMoNoPH/Biometrical-Journal/Sim-Results"


##########


# Helper functions
##########

# --- Helper: make short, readable labels for facets/plots
shorten_config <- function(cfg) {
  # Updated pattern to handle BP, GB, GB_Quantile, GP and scenario
  parts <- stringr::str_match(
    cfg,
    "^([A-Za-z_]+)_k(\\d+)_n(\\d+)_cens(\\d+)_([A-Za-z]+)_scen(\\d+)_?(FU\\d+)?$"
  )
  if (is.na(parts[1, 1])) return(cfg)
  
  bf   <- parts[2]
  k    <- parts[3]
  n    <- parts[4]
  cens <- parts[5]
  dist_letter <- toupper(substr(parts[6], 1, 1))  # "L" for loglogistic
  scen <- parts[7]
  
  paste0(bf, " k=", k, ", c=", cens, "%, S", scen)
}


# Function to parse a config string 
parse_config <- function(config_str) {
  # Updated pattern to handle BP, GB, GB_Quantile, GP and scenario
  parts <- stringr::str_match(config_str, "^([A-Za-z_]+)_k(\\d+)_n(\\d+)_cens(\\d+)_([A-Za-z]+)_scen(\\d+)_?(FU\\d+)?$")
  
  if (is.na(parts[1, 1])) return(NULL)
  
  list(
    bf = parts[2],
    k_bases = as.numeric(parts[3]),
    n = as.numeric(parts[4]),
    cens = as.numeric(parts[5]),
    dist = parts[6],
    scenario = as.numeric(parts[7]),
    fu = parts[8]
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
  # sigma_0 = 15, sigma_1 = 0.20, alpha_AFT = 0.015
  # Scenario-specific:
  # Scenario 1: beta_2 = 0.00, log_AF = 0.00
  # Scenario 2: beta_2 = 0.04, log_AF = 0.90
  # Scenario 3: beta_2 = -0.04, log_AF = -0.90
  # ---- TRUE VALUES (scenarios 1–5 only), from supplied true_params ----
  stopifnot(scenario %in% 1:5)
  
  true_values <- c(
    "beta_long[1]"  = 73,     # beta_0
    "beta_long[2]"  = -0.04,  # beta_1
    "beta_long[3]"  = c(0.00, 0.04, -0.04, 0.04, -0.04)[scenario],  # beta_2 (scen 1–5)
    
    "sd_1_long[1]"  = 15,     # sigma_b0
    "sd_1_long[2]"  = 0.2,    # sigma_b1
    "sigma_long"    = 12,     # sigma_epsilon
    
    "alpha"         = 0.012,  # alpha (scen 1–5)
    "beta_surv[1]"  = c(0.00, 0.90, 0.90, -0.90, -0.90)[scenario]   # gamma_1 (scen 1–5)
  )
  
  
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
    summary       = combined_summary # aggregated table (Stan only)
  )
}

##########


# Configuration parsing and execution
##########

# ===== Run across all configs =====

# Get all unique config stubs from files (strip last 9 chars before .rds)
filenames <- list.files(out_dir, pattern = "\\.rds$")
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

# Build a long table of summaries with the config label attached
all_summary_tables <- purrr::map2_dfr(
  results_by_config, names(results_by_config),
  ~ dplyr::mutate(.x$summary, config = .y),
  .id = NULL
)

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
reference_lines <- tibble::tibble(
  variable = c("beta_long[1]","beta_long[2]","beta_long[3]",
               "sd_1_long[1]","sd_1_long[2]","sigma_long",
               "beta_surv[1]","alpha"),
  reference = c(73, -0.04, NA, 15, 0.20, 12, NA, 0.015)
)

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
      variable == "beta_long[3]" & scenario == 4 ~ 0.08,
      variable == "beta_long[3]" & scenario == 5 ~ -0.08,
      variable == "beta_long[3]" & scenario == 6 ~ 0.12,
      variable == "beta_long[3]" & scenario == 7 ~ -0.12,
      variable == "beta_surv[1]" & scenario == 1 ~ 0.00,
      variable == "beta_surv[1]" & scenario == 2 ~ 0.90,
      variable == "beta_surv[1]" & scenario == 3 ~ -0.90,
      variable == "beta_surv[1]" & scenario == 4 ~ 1.80,
      variable == "beta_surv[1]" & scenario == 5 ~ -1.80,
      variable == "beta_surv[1]" & scenario == 6 ~ 2.70,
      variable == "beta_surv[1]" & scenario == 7 ~ -2.70,
      TRUE ~ reference
    ),
    bias      = mean - reference,
    rel_bias  = if_else(reference != 0, bias / reference, NA_real_),
    rmse      = sqrt(sd^2 + bias^2)   # RMSE^2 = Var + Bias^2
  )

results_tab2 <- results_tab

nice_tab_wide <- results_tab2 %>%
  filter(
    variable %in% key_vars,
    source == "Stan posterior"
  ) %>%
  mutate(
    Parameter = variable,
    Basis     = bf,
    k         = as.integer(k_bases),
    n         = as.integer(n),
    Cens      = as.integer(cens),
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
  filter(
    stringr::str_detect(variable, "^gamma\\["),
    source == "Stan posterior"
  ) %>%
  mutate(
    Parameter = variable,
    Basis     = bf,
    k         = as.integer(k_bases),
    n         = as.integer(n),
    Cens      = as.integer(cens),
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
  filter(
    stringr::str_detect(variable, "^(gp_length_scale|gp_marginal_sd)$"),
    source == "Stan posterior"
  ) %>%
  mutate(
    Parameter = variable,
    Basis     = bf,
    k         = as.integer(k_bases),
    n         = as.integer(n),
    Cens      = as.integer(cens),
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

# Desired facet order for configs (ENZAMET specific: n=1100, loglogistic, scenarios 1-7, censoring 0%/50%)
# Format: "BF k=K, c=C%, SK" where BF=BP/GB, K=5/7/9, C=0/50, SK=scenario
# Total: 7 scenarios × 2 censoring levels × 2 basis types × 3 k values = 84 configurations
desired_order <- c(
  # 0% Censoring, Scenarios 1-7
  "BP k=5, c=0%, S1",  "GB k=5, c=0%, S1",
  "BP k=7, c=0%, S1",  "GB k=7, c=0%, S1",
  "BP k=9, c=0%, S1",  "GB k=9, c=0%, S1",
  
  "BP k=5, c=0%, S2",  "GB k=5, c=0%, S2",
  "BP k=7, c=0%, S2",  "GB k=7, c=0%, S2",
  "BP k=9, c=0%, S2",  "GB k=9, c=0%, S2",
  
  "BP k=5, c=0%, S3",  "GB k=5, c=0%, S3",
  "BP k=7, c=0%, S3",  "GB k=7, c=0%, S3",
  "BP k=9, c=0%, S3",  "GB k=9, c=0%, S3",
  
  "BP k=5, c=0%, S4",  "GB k=5, c=0%, S4",
  "BP k=7, c=0%, S4",  "GB k=7, c=0%, S4",
  "BP k=9, c=0%, S4",  "GB k=9, c=0%, S4",
  
  "BP k=5, c=0%, S5",  "GB k=5, c=0%, S5",
  "BP k=7, c=0%, S5",  "GB k=7, c=0%, S5",
  "BP k=9, c=0%, S5",  "GB k=9, c=0%, S5",
  
  "BP k=5, c=0%, S6",  "GB k=5, c=0%, S6",
  "BP k=7, c=0%, S6",  "GB k=7, c=0%, S6",
  "BP k=9, c=0%, S6",  "GB k=9, c=0%, S6",
  
  "BP k=5, c=0%, S7",  "GB k=5, c=0%, S7",
  "BP k=7, c=0%, S7",  "GB k=7, c=0%, S7",
  "BP k=9, c=0%, S7",  "GB k=9, c=0%, S7",
  
  # 50% Censoring, Scenarios 1-7
  "BP k=5, c=50%, S1",  "GB k=5, c=50%, S1",
  "BP k=7, c=50%, S1",  "GB k=7, c=50%, S1",
  "BP k=9, c=50%, S1",  "GB k=9, c=50%, S1",
  
  "BP k=5, c=50%, S2",  "GB k=5, c=50%, S2",
  "BP k=7, c=50%, S2",  "GB k=7, c=50%, S2",
  "BP k=9, c=50%, S2",  "GB k=9, c=50%, S2",
  
  "BP k=5, c=50%, S3",  "GB k=5, c=50%, S3",
  "BP k=7, c=50%, S3",  "GB k=7, c=50%, S3",
  "BP k=9, c=50%, S3",  "GB k=9, c=50%, S3",
  
  "BP k=5, c=50%, S4",  "GB k=5, c=50%, S4",
  "BP k=7, c=50%, S4",  "GB k=7, c=50%, S4",
  "BP k=9, c=50%, S4",  "GB k=9, c=50%, S4",
  
  "BP k=5, c=50%, S5",  "GB k=5, c=50%, S5",
  "BP k=7, c=50%, S5",  "GB k=7, c=50%, S5",
  "BP k=9, c=50%, S5",  "GB k=9, c=50%, S5",
  
  "BP k=5, c=50%, S6",  "GB k=5, c=50%, S6",
  "BP k=7, c=50%, S6",  "GB k=7, c=50%, S6",
  "BP k=9, c=50%, S6",  "GB k=9, c=50%, S6",
  
  "BP k=5, c=50%, S7",  "GB k=5, c=50%, S7",
  "BP k=7, c=50%, S7",  "GB k=7, c=50%, S7",
  "BP k=9, c=50%, S7",  "GB k=9, c=50%, S7"
)

all_filtered_data$config <- factor(all_filtered_data$config, levels = desired_order)

# Reference lines for key parameters
# Note: beta_long[3] and beta_surv[1] will be added dynamically per scenario in plots
reference_lines_plot <- tibble::tibble(
  variable = c("beta_long[1]","beta_long[2]",
               "sd_1_long[1]","sd_1_long[2]","sigma_long",
               "alpha"),
  reference = c(73, -0.04, 15, 0.20, 12, 0.015)
)

# Scenario-specific reference values for faceting
scenario_refs <- tibble::tibble(
  variable = rep(c("beta_long[3]", "beta_surv[1]"), each = 7),
  scenario = rep(1:7, 2),
  reference = c(0.00, 0.04, -0.04, 0.08, -0.08, 0.12, -0.12,  # beta_long[3] for scenarios 1-7
                0.00, 0.90, -0.90, 1.80, -1.80, 2.70, -2.70)   # beta_surv[1] for scenarios 1-7
)

# Shared palette - for BP/GB combinations
# Generate 84 distinct colors (7 scenarios × 2 censoring × 2 basis × 3 k values)
custom_palette <- rep(c(
  "#A6CEE3", "#1F78B4",  # BP k=5, GB k=5
  "#CAB2D6", "#6A3D9A",  # BP k=7, GB k=7
  "#FCCDE5", "#BC80BD"   # BP k=9, GB k=9
), times = 14)  # Repeat for 7 scenarios × 2 censoring levels

##########


# Figure generation
##########

### Stan posterior means

# Desired facet order
var_order <- reference_lines$variable


### Loglogistic distribution figures (ENZAMET only uses loglogistic)

message("Generating loglogistic figures for ENZAMET scenarios...")
message("Total rows in all_filtered_data: ", nrow(all_filtered_data))
message("Unique configs: ", paste(unique(all_filtered_data$config), collapse = ", "))

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
    title = "Posterior Means - ENZAMET (Loglogistic Hazard, 7 Scenarios, 0%/50% Censoring)"
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

message("✓ Saved ENZAMET figure to: ", fig_dir)
message("  - fig_ENZAMET_loglogistic.pdf")

} else {
  message("No data available for figure generation")
}

##########

# Baseline hazard plots
##########
# 
# Plot baseline hazards for the first 5 simulations using the basis weights
# (gamma parameters) to visualize how well the basis approximation captures
# the underlying hazard function.
#

# Function to create Gaussian basis functions
create_gaussian_basis <- function(t, n_bases = 5, t_max = 72) {
  # Create evenly spaced centers across time range
  centers <- seq(0, t_max, length.out = n_bases)
  bandwidth <- (t_max / n_bases) * 0.5  # Overlap between bases
  
  # Create basis matrix
  basis_matrix <- matrix(0, nrow = length(t), ncol = n_bases)
  for (j in 1:n_bases) {
    basis_matrix[, j] <- exp(-((t - centers[j])^2) / (2 * bandwidth^2))
  }
  
  return(basis_matrix)
}

# Function to create Bernstein polynomial basis (b_(k,m))
create_bernstein_basis <- function(t, n_bases = 5, t_max = 72) {
  # Normalize time to [0, 1]
  t_norm <- t / t_max
  t_norm <- pmin(pmax(t_norm, 0), 1)  # Clamp to [0, 1]
  
  # Bernstein polynomial basis of degree n_bases - 1
  m <- n_bases - 1
  basis_matrix <- matrix(0, nrow = length(t), ncol = n_bases)
  
  for (k in 0:m) {
    basis_matrix[, k + 1] <- choose(m, k) * t_norm^k * (1 - t_norm)^(m - k)
  }
  
  return(basis_matrix)
}

# Function to extract gamma weights from simulation
extract_gamma_weights <- function(sim_data, n_bases = 5) {
  gamma_vars <- paste0("gamma[", 1:n_bases, "]")
  gamma_df <- sim_data %>%
    dplyr::filter(variable %in% gamma_vars) %>%
    dplyr::arrange(variable)
  
  if (nrow(gamma_df) < n_bases) {
    return(NULL)
  }
  
  return(gamma_df$mean)
}

# Select configurations to plot - dynamically discover from available files
# Get all unique configs from output directory
all_rds_files <- list.files(out_dir, pattern = "\\.rds$", full.names = FALSE)
all_configs <- substr(tools::file_path_sans_ext(all_rds_files),
                      1, nchar(tools::file_path_sans_ext(all_rds_files)) - 9) %>%
  unique()

# Filter for selected models with k_bases to plot baseline hazards
# ENZAMET: n=1100, cens=0 or 50, loglogistic, FU72, scenarios 1-5 ONLY
# Use the models_to_plot vector defined at the top of the script
models_pattern <- paste0("^(", paste(models_to_plot, collapse = "|"), ")_k[579]_n1100_cens(0|50)_loglogistic_scen[1-5]_FU72$")
configs_to_plot <- all_configs[grepl(models_pattern, all_configs)]

if (length(configs_to_plot) == 0) {
  message("No matching configs found for baseline hazard plots. Available configs:")
  message(paste(head(all_configs, 10), collapse = "\n"))
  message("\nLooking for models: ", paste(models_to_plot, collapse = ", "))
}

# Constructing fixed basis knots for Gaussian Bases
build_fixed_basis_knots <- function(m) {
  mu_fix    <- numeric(m)
  sigma_fix <- rep(2.0 / (3.0 * (m - 1.0)), m)
  
  mu_fix[1] <- 0
  mu_fix[m] <- 1
  
  if (m > 2) {
    mu_fix[2:(m-1)] <- seq(1, m - 2) / (m - 1)
  }
  
  list(mu = mu_fix, sigma = sigma_fix)
}

# Create the gaussian bases for Stan-like computation
create_gaussian_basis_stan <- function(u, mu, sigma) {
  stopifnot(length(mu) == length(sigma))
  m <- length(mu)
  
  basis_mat <- sapply(seq_len(m), function(k) {
    dnorm(u, mean = mu[k], sd = sigma[k])  # Normal pdf, like exp(normal_lpdf)
  })
  # n_grid x m
  basis_mat
}

# Create Bernstein polynomial basis functions
create_bernstein_basis <- function(u, m) {
  # Bernstein polynomial basis: B_k(u) = dbeta(u | k, m-k+1)
  # where k = 1, 2, ..., m
  basis_mat <- sapply(seq_len(m), function(k) {
    dbeta(u, shape1 = k, shape2 = m - k + 1)
  })
  # n_grid x m
  basis_mat
}

reconstruct_h0_u <- function(gamma_weights, m, n_grid = 200, basis_type = "GB") {
  u_grid <- seq(0, 1, length.out = n_grid)
  
  if (basis_type == "GB") {
    knots <- build_fixed_basis_knots(m)
    basis <- create_gaussian_basis_stan(u_grid, knots$mu, knots$sigma)
  } else if (basis_type == "BP") {
    basis <- create_bernstein_basis(u_grid, m)
  } else {
    stop("basis_type must be 'GB' or 'BP'")
  }
  
  h0_u  <- as.vector(basis %*% gamma_weights)  # n_grid x 1
  
  tibble::tibble(u = u_grid, h0 = h0_u)
}




message("\n=== BASELINE HAZARD ESTIMATION ===")
message("Generating baseline hazard plots for first 5 simulations...")
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
  
  # Read first batch and extract first 5 simulations
  sims <- readRDS(rds_files[1])
  n_sims_to_plot <- min(5, length(sims))
  
  if (n_sims_to_plot == 0) {
    message("No simulations in file: ", rds_files[1])
    next
  }
  
  # Extract model type from summaries
  model_type <- sims[[1]]$joint_summary$model[1]
  message("  Model type from data: ", model_type)
  message("  Plotting ", n_sims_to_plot, " simulations")
  
  # Time grid for plotting (ENZAMET uses FU=72 months)
  t_grid <- seq(0, 72, length.out = 200)
  
  # Create basis matrix based on model type
  if (model_type == "GB" || model_type == "GB_Quantile") {
    knots       <- build_fixed_basis_knots(k_bases)
    u_grid      <- t_grid / max(t_grid)
    basis_matrix <- create_gaussian_basis_stan(u_grid, knots$mu, knots$sigma)
    basis_name  <- if (model_type == "GB") "Gaussian Basis" else "GB Quantile"
  } else if (model_type == "BP") {
    u_grid      <- t_grid / max(t_grid)
    basis_matrix <- create_bernstein_basis(u_grid, m = k_bases)
    basis_name   <- "Bernstein Polynomial"
  } else if (model_type == "GP") {
    u_grid      <- seq(0, 1, length.out = k_bases)
    basis_matrix <- NULL
    basis_name   <- "Gaussian Process"
  } else {
    message("Unknown model type: ", model_type)
    next
  }
  
  
  # Initialize data frame for plotting
  hazard_data <- data.frame()
  true_hazard_data <- data.frame()
  
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
  
  # Extract gamma weights (for BP/GB) or f_gp values (for GP) and compute baseline hazards
  for (sim_idx in 1:n_sims_to_plot) {
    if (model_type == "GP") {
      # Extract f_gp (log hazard at knots) for GP model
      f_gp_vars <- paste0("f_gp[", 1:k_bases, "]")
      f_gp_df <- sims[[sim_idx]]$joint_summary %>%
        dplyr::filter(variable %in% f_gp_vars) %>%
        dplyr::arrange(variable)
      
      if (nrow(f_gp_df) < k_bases) next
      
      f_gp_values <- f_gp_df$mean
      
      # Linear interpolation of log hazard, then exponentiate
      u_grid_fine <- seq(0, 1, length.out = length(t_grid))
      log_h0_interp <- approx(u_grid, f_gp_values, xout = u_grid_fine, rule = 2)$y
      h0_t <- exp(log_h0_interp)
      
    } else {
      # BP/GB models use gamma weights
      gamma_weights <- extract_gamma_weights(sims[[sim_idx]]$joint_summary,
                                             n_bases = k_bases)
      if (is.null(gamma_weights)) next
      
      h0_t <- basis_matrix %*% gamma_weights
    }
    
    sim_data <- data.frame(
      time   = t_grid,
      hazard = as.vector(h0_t),
      sim_id = factor(sim_idx),
      type   = "Estimated"
    )
    
    hazard_data <- rbind(hazard_data, sim_data)
  }
  
  
  # Generate true baseline hazard for each simulation
  # ENZAMET uses loglogistic distribution
  for (sim_idx in 1:n_sims_to_plot) {
    # For Loglogistic: h0(t) = (shape/scale) * (t/scale)^(shape-1) / (1 + (t/scale)^shape)
    true_h0 <- (loglogistic_shape / loglogistic_scale) * 
               (t_grid / loglogistic_scale)^(loglogistic_shape - 1) / 
               (1 + (t_grid / loglogistic_scale)^loglogistic_shape)
    
    # Scale true hazard by dividing by its maximum
    max_h0 <- max(true_h0, na.rm = TRUE)
    if (!is.na(max_h0) && max_h0 > 0) {
      true_h0 <- true_h0 / max_h0
    }
    
    true_sim_data <- data.frame(
      time = t_grid,
      hazard = true_h0,
      sim_id = factor(sim_idx),
      type = "True"
    )
    
    true_hazard_data <- rbind(true_hazard_data, true_sim_data)
  }
  
  # Combine estimated and true hazards
  all_hazard_data <- rbind(hazard_data, true_hazard_data)
  
  # Get covariate effects from true parameter values for AFT model
  beta_surv_1 <- 0.25
  X_ref <- 1  # Reference covariate value
  time_acceleration <- exp(-beta_surv_1 * X_ref)  # AFT time scaling
  
  # Compute cumulative hazards via trapezoidal integration
  cumhaz_data <- data.frame()
  for (sim_id_val in unique(all_hazard_data$sim_id)) {
    for (type_val in c("Estimated", "True")) {
      subset_df <- all_hazard_data %>%
        dplyr::filter(sim_id == sim_id_val, type == type_val) %>%
        dplyr::arrange(time)
      
      if (nrow(subset_df) > 1) {
        # For AFT model, transform time scale: t_transformed = t * exp(-beta * X)
        time_transformed <- subset_df$time * time_acceleration
        
        # Trapezoidal rule on transformed time scale
        H_cumsum <- cumsum(c(0, diff(time_transformed) * 
                                (subset_df$hazard[-1] + subset_df$hazard[-nrow(subset_df)]) / 2))
        
        cumhaz_df <- data.frame(
          time = subset_df$time,  # Original time scale for plotting
          cumhaz = H_cumsum,
          sim_id = sim_id_val,
          type = type_val
        )
        cumhaz_data <- rbind(cumhaz_data, cumhaz_df)
      }
    }
  }
  
  # Create cumulative hazard plot
  if (nrow(cumhaz_data) > 0) {
    
    p_cumhaz <- ggplot2::ggplot() +
      # Plot true cumulative hazards (black, thicker)
      ggplot2::geom_line(
        data = cumhaz_data %>% dplyr::filter(type == "True"),
        ggplot2::aes(x = time, y = cumhaz, group = sim_id),
        color = "black",
        linewidth = 1.2,
        linetype = "solid"
      ) +
      # Plot estimated cumulative hazards (colored by sim_id)
      ggplot2::geom_line(
        data = cumhaz_data %>% dplyr::filter(type == "Estimated"),
        ggplot2::aes(x = time, y = cumhaz, group = sim_id, color = sim_id),
        linewidth = 0.8,
        linetype = "dashed"
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
          "4" = "#984EA3", "5" = "#FF7F00"
        ),
        name = "Estimated (Sim ID)"
      ) +
      ggplot2::labs(
        title = paste0("Cumulative Baseline Hazard: ENZAMET Scenario ", scenario, " (", basis_name, ", k=", k_bases, ")"),
        subtitle = paste0("n=", config_parsed$n, ", ", config_parsed$cens, "% censoring, ", 
                          "loglogistic hazard, FU=72 months. ",
                          "Black solid = True, Colored dashed = Estimated"),
        x = "Time (months)",
        y = "Cumulative Baseline Hazard H_0(t)"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold", size = 11)
      )
    
    print(p_cumhaz)
    
    # Save cumulative hazard plot
    fig_name_cumhaz <- paste0("fig_cumhaz_", gsub("_FU72", "", config))
    
    ggplot2::ggsave(filename = file.path(fig_dir, paste0(fig_name_cumhaz, ".pdf")),
                    plot = p_cumhaz, width = 10, height = 6, units = "in")
    
    # Store plot for compilation
    all_cumhaz_plots[[config]] <- p_cumhaz
    
    message("  ✓ Cumulative baseline hazard plot created and saved: ", fig_name_cumhaz)
    
  } else {
    message("  ✗ No cumulative baseline hazard data to plot for: ", config)
  }
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

# ------------------------------------------------------------
# Parse metadata from SHORT config label: "BP k=5, c=0%, S1"
# ------------------------------------------------------------
parse_short_config <- function(lbl) {
  m <- stringr::str_match(lbl, "^([A-Za-z_]+)\\s+k=(\\d+),\\s+c=(\\d+)%?,\\s+S(\\d+)$")
  if (is.na(m[1,1])) {
    return(tibble(bf = NA_character_, k_bases = NA_real_, cens = NA_real_, scenario = NA_real_))
  }
  tibble(
    bf       = m[1,2],
    k_bases  = as.numeric(m[1,3]),
    cens     = as.numeric(m[1,4]),
    scenario = as.numeric(m[1,5])
  )
}

config_info <- tibble(config = unique(all_filtered_data$config)) %>%
  mutate(parsed = map(config, parse_short_config)) %>%
  tidyr::unnest(parsed)

# ------------------------------------------------------------
# Build per-simulation performance table (creates distributions)
# ------------------------------------------------------------
perf_sim <- all_filtered_data %>%
  filter(variable %in% key_vars) %>%
  left_join(config_info, by = "config") %>%
  filter(!is.na(scenario), scenario %in% 1:5) %>%
  left_join(truth_15, by = c("scenario", "variable")) %>%
  filter(!is.na(truth)) %>%
  mutate(
    err      = mean - truth,
    bias     = err,
    abs_err  = abs(err),   # distribution across sims (absolute error)
    sq_error = err^2,
    covered  = as.integer(q2.5 <= truth & truth <= q97.5),
    
    # Use manuscript-friendly labels for facets
    Parameter = factor(variable, levels = key_vars, labels = unname(param_labels[key_vars])),
    
    Scenario  = factor(scenario),
    Cens      = factor(cens),
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
  filter(as.character(Parameter) == "beta[0]") %>%
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
    se = sqrt(coverage * (1 - coverage) / n),
    lo = pmax(0, coverage - 1.96 * se),
    hi = pmin(1, coverage + 1.96 * se),
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
