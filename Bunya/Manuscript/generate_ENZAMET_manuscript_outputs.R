#!/usr/bin/env Rscript

## Manuscript-facing tables and figures for the ENZAMET simulations.
##
## This script reads summary-level simulation outputs only. It intentionally
## does not read per-simulation or raw posterior objects, so it should remain
## fast and easy to rerun while the manuscript tables evolve.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1] %||% file.path(getwd(), "Bunya", "Manuscript"))
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "Bunya", "ENZAMET"))) {
  repo_root <- normalizePath(getwd(), mustWork = FALSE)
}

enzamet_dir <- file.path(repo_root, "Bunya", "ENZAMET")
manuscript_dir <- Sys.getenv(
  "ENZ_MANUSCRIPT_DIR",
  file.path(repo_root, "Bunya", "Manuscript")
)
table_dir <- file.path(manuscript_dir, "tables")
figure_dir <- file.path(manuscript_dir, "figures")
csv_dir <- file.path(manuscript_dir, "derived_csv")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  proposed_40 = file.path(
    enzamet_dir,
    "summary_tables_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo"
  ),
  jmbayes2_40 = file.path(
    enzamet_dir,
    "summary_tables_20260602_ENZAMET-JMbayes2-DM-all40-1000rep-randSDsummary"
  ),
  alpha_maxF125 = file.path(
    enzamet_dir,
    "summary_tables_20260528_ENZAMET-alphaSens-LL-maxF125-altDM-gamma5-1000rep-zeroprior-initinfo"
  ),
  alpha_maxF150 = file.path(
    enzamet_dir,
    "summary_tables_20260528_ENZAMET-alphaSens-LL-maxF150-altDM-gamma5-1000rep-zeroprior-initinfo"
  ),
  ## Baseline alpha-prior setting in the main 40-config run.
  ## User note, 2026-06-04: treat this setting as maxF200 in the
  ## alpha-prior sensitivity tables and plots.
  alpha_maxF200 = file.path(
    enzamet_dir,
    "summary_tables_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo"
  ),
  alpha_maxF300 = file.path(
    enzamet_dir,
    "summary_tables_20260528_ENZAMET-alphaSens-LL-maxF300-altDM-gamma5-1000rep-zeroprior-initinfo"
  )
)

required_files <- c(
  file.path(paths$proposed_40, "jm_parameter_summary_long_loglogistic_weibull.csv"),
  file.path(paths$jmbayes2_40, "jmbayes2_comparison_summary_loglogistic_weibull.csv"),
  file.path(paths$alpha_maxF125, "jm_parameter_summary_long_loglogistic.csv"),
  file.path(paths$alpha_maxF150, "jm_parameter_summary_long_loglogistic.csv"),
  file.path(paths$alpha_maxF200, "jm_parameter_summary_long_loglogistic_weibull.csv"),
  file.path(paths$alpha_maxF300, "jm_parameter_summary_long_loglogistic.csv")
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop(
    "Missing required summary files:\n",
    paste(missing_required, collapse = "\n"),
    call. = FALSE
  )
}

read_summary_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

format_number <- function(x, digits = 3) {
  out <- rep("NA", length(x))
  ok <- !is.na(x) & is.finite(x)
  if (any(ok)) {
    formatted <- formatC(x[ok], format = "f", digits = digits)
    ## Keep signs for tiny negative values after ordinary fixed-decimal
    ## rounding, e.g. -4e-05 -> -0.0000 with digits = 4.
    negative_zero <- x[ok] < 0 &
      round(abs(x[ok]), digits = digits) == 0 &
      !startsWith(formatted, "-")
    formatted[negative_zero] <- paste0("-", formatted[negative_zero])
    out[ok] <- formatted
  }
  out
}

format_table <- function(df, digits = 3, column_digits = list()) {
  out <- df
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      dg <- column_digits[[nm]] %||% digits
      if (grepl("^n($|_)|N cfg|N configs|N bias|Scenario|FU|Sample size|Basis|Follow-up|maxF", nm)) {
        out[[nm]] <- format_number(out[[nm]], digits = 0)
      } else {
        out[[nm]] <- format_number(out[[nm]], digits = dg)
      }
    } else {
      out[[nm]] <- ifelse(is.na(out[[nm]]), "NA", as.character(out[[nm]]))
    }
  }
  out
}

write_latex_table <- function(df, path, caption, label, digits = 3,
                              column_digits = list(), align = NULL,
                              table_env = TRUE, font_size = "\\small") {
  display <- format_table(df, digits = digits, column_digits = column_digits)
  display[] <- lapply(display, latex_escape)
  display_names <- names(display)
  display_names[display_names == "Scenario"] <- "Scen."

  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", max(ncol(display) - 1, 0)), collapse = ""))
  }

  lines <- character()
  if (table_env) {
    lines <- c(
      lines,
      "\\begin{table}[!htbp]",
      "\\centering",
      font_size,
      paste0("\\caption{", latex_escape(caption), "}"),
      paste0("\\label{", label, "}")
    )
  }
  lines <- c(
    lines,
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(latex_escape(display_names), collapse = " & "),
    "\\\\",
    "\\midrule"
  )
  body <- apply(display, 1, function(row) paste(row, collapse = " & "))
  lines <- c(lines, paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}")
  if (table_env) {
    lines <- c(lines, "\\end{table}")
  }
  writeLines(lines, path)
  message("Wrote LaTeX table: ", path)
}

parameter_math_label <- function(variable) {
  dplyr::recode(
    variable,
    "beta_surv[1]" = "$\\gamma$",
    "alpha" = "$\\alpha$",
    "beta_long[1]" = "$\\beta_0$",
    "beta_long[2]" = "$\\beta_1$",
    "beta_long[3]" = "$\\beta_2$",
    "sd_1_long[1]" = "$\\sigma_{b0}$",
    "sd_1_long[2]" = "$\\sigma_{b1}$",
    "sigma_long" = "$\\sigma_{\\epsilon}$",
    .default = variable
  )
}

write_alpha_prior_estimate_table <- function(df, path, caption, label,
                                             prior_labels, cm_labels,
                                             prior_display_labels = prior_labels,
                                             cm_display_labels = cm_labels,
                                             super_header = NULL,
                                             na_text = "NA",
                                             note = paste0(
                                               "\\footnotesize Notes: $\\text{CM}_1$ = administrative censoring; ",
                                               "$\\text{CM}_2$ = controlled 50\\% censoring. ",
                                               "Entries are means of posterior means across simulations."
                                             )) {
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{llr", paste(rep("r", length(prior_labels) * length(cm_labels)), collapse = ""), "}"),
    "\\toprule"
  )

  n_estimate_cols <- length(prior_labels) * length(cm_labels)
  if (!is.null(super_header)) {
    header_super <- c(
      "Scen.", "Par.", "True",
      paste0("\\multicolumn{", n_estimate_cols, "}{c}{", super_header, "}")
    )
    lines <- c(
      lines,
      paste(header_super, collapse = " & "), "\\\\",
      paste0("\\cmidrule(lr){4-", 3 + n_estimate_cols, "}")
    )
  }

  header_1 <- c(
    if (is.null(super_header)) c("Scen.", "Par.", "True") else c("", "", ""),
    vapply(
      prior_display_labels,
      function(pr) paste0("\\multicolumn{", length(cm_labels), "}{c}{", pr, "}"),
      character(1)
    )
  )
  lines <- c(lines, paste(header_1, collapse = " & "), "\\\\")

  cm_start <- 4
  cmid <- unlist(lapply(seq_along(prior_labels), function(i) {
    start <- cm_start + (i - 1) * length(cm_labels)
    end <- start + length(cm_labels) - 1
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }))
  lines <- c(lines, paste(cmid, collapse = " "))

  header_2 <- c("", "", "", rep(cm_display_labels, times = length(prior_labels)))
  lines <- c(lines, paste(header_2, collapse = " & "), "\\\\", "\\midrule")

  for (scenario_id in unique(df$Scenario)) {
    scenario_rows <- df[df$Scenario == scenario_id, , drop = FALSE]
    for (i in seq_len(nrow(scenario_rows))) {
      row <- scenario_rows[i, ]
      scenario_cell <- ifelse(i == 1, as.character(scenario_id), "")
      estimate_cells <- unlist(lapply(prior_labels, function(pr) {
        vapply(cm_labels, function(cm) {
          col_nm <- paste0(pr, "_", cm)
          value <- if (col_nm %in% names(row)) row[[col_nm]] else NA_real_
          if (is.na(value) || !is.finite(value)) na_text else format_number(value, digits = 3)
        }, character(1))
      }))
      line <- paste(
        c(
          latex_escape(scenario_cell),
          row$Par,
          format_number(row$True, digits = 3),
          estimate_cells
        ),
        collapse = " & "
      )
      lines <- c(lines, paste0(line, " \\\\"))
    }
    if (scenario_id != tail(unique(df$Scenario), 1)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.5ex}",
    "\\begin{minipage}{0.95\\linewidth}",
    note,
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, path)
  message("Wrote LaTeX table: ", path)
}

write_distribution_model_estimate_table <- function(df, path, caption, label,
                                                    distribution_labels,
                                                    cm_labels,
                                                    model_labels,
                                                    distribution_display_labels = distribution_labels,
                                                    cm_display_labels = cm_labels,
                                                    model_display_labels = model_labels,
                                                    na_text = "--",
                                                    note = paste0(
                                                      "\\footnotesize Notes: $\\text{CM}_1$ = administrative censoring; ",
                                                      "$\\text{CM}_2$ = controlled 50\\% censoring. ",
                                                      "Entries are means of point estimates across simulations."
                                                    )) {
  n_model <- length(model_labels)
  n_cm <- length(cm_labels)
  n_distribution <- length(distribution_labels)
  n_estimate_cols <- n_distribution * n_cm * n_model

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\scriptsize",
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", label, "}"),
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{llr", paste(rep("r", n_estimate_cols), collapse = ""), "}"),
    "\\toprule"
  )

  header_1 <- c(
    "Scen.", "Par.", "True",
    vapply(
      distribution_display_labels,
      function(distribution_label) paste0("\\multicolumn{", n_cm * n_model, "}{c}{", distribution_label, "}"),
      character(1)
    )
  )
  lines <- c(lines, paste(header_1, collapse = " & "), "\\\\")

  distribution_start <- 4
  distribution_cmid <- unlist(lapply(seq_along(distribution_labels), function(i) {
    start <- distribution_start + (i - 1) * n_cm * n_model
    end <- start + n_cm * n_model - 1
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }))
  lines <- c(lines, paste(distribution_cmid, collapse = " "))

  header_2 <- c(
    "", "", "",
    unlist(lapply(seq_along(distribution_labels), function(i) {
      vapply(
        cm_display_labels,
        function(cm_label) paste0("\\multicolumn{", n_model, "}{c}{", cm_label, "}"),
        character(1)
      )
    }))
  )
  lines <- c(lines, paste(header_2, collapse = " & "), "\\\\")

  cm_start <- 4
  cm_cmid <- unlist(lapply(seq_len(n_distribution * n_cm), function(i) {
    start <- cm_start + (i - 1) * n_model
    end <- start + n_model - 1
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }))
  lines <- c(lines, paste(cm_cmid, collapse = " "))

  header_3 <- c("", "", "", rep(model_display_labels, times = n_distribution * n_cm))
  lines <- c(lines, paste(header_3, collapse = " & "), "\\\\", "\\midrule")

  for (scenario_id in unique(df$Scenario)) {
    scenario_rows <- df[df$Scenario == scenario_id, , drop = FALSE]
    for (i in seq_len(nrow(scenario_rows))) {
      row <- scenario_rows[i, ]
      scenario_cell <- ifelse(i == 1, as.character(scenario_id), "")
      estimate_cells <- unlist(lapply(distribution_labels, function(distribution_label) {
        unlist(lapply(cm_labels, function(cm) {
          vapply(model_labels, function(model_label) {
            col_nm <- paste(distribution_label, cm, model_label, sep = "_")
            value <- if (col_nm %in% names(row)) row[[col_nm]] else NA_real_
            if (is.na(value) || !is.finite(value)) na_text else format_number(value, digits = 3)
          }, character(1))
        }))
      }))
      line <- paste(
        c(
          latex_escape(scenario_cell),
          row$Par,
          format_number(row$True, digits = 3),
          estimate_cells
        ),
        collapse = " & "
      )
      lines <- c(lines, paste0(line, " \\\\"))
    }
    if (scenario_id != tail(unique(df$Scenario), 1)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.5ex}",
    "\\begin{minipage}{0.95\\linewidth}",
    note,
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, path)
  message("Wrote LaTeX table: ", path)
}

write_alpha_prior_detail_table <- function(df, path, caption, label,
                                           prior_labels, metric_labels,
                                           censoring_label,
                                           note = NULL,
                                           prior_display_labels = prior_labels) {
  n_metric <- length(metric_labels)
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\scriptsize",
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", label, "}"),
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{ll", paste(rep("r", length(prior_labels) * n_metric), collapse = ""), "}"),
    "\\toprule"
  )

  header_1 <- c(
    "Scen.", "Par.",
    vapply(
      prior_display_labels,
      function(pr) paste0("\\multicolumn{", n_metric, "}{c}{", pr, "}"),
      character(1)
    )
  )
  lines <- c(lines, paste(header_1, collapse = " & "), "\\\\")

  cm_start <- 3
  cmid <- unlist(lapply(seq_along(prior_labels), function(i) {
    start <- cm_start + (i - 1) * n_metric
    end <- start + n_metric - 1
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }))
  lines <- c(lines, paste(cmid, collapse = " "))

  header_2 <- c("", "", rep(names(metric_labels), times = length(prior_labels)))
  lines <- c(lines, paste(latex_escape(header_2), collapse = " & "), "\\\\", "\\midrule")

  metric_digits <- c(
    Bias = 4,
    ESD = 4,
    PSD = 4,
    RMSE = 4,
    CP = 3
  )

  for (scenario_id in unique(df$Scenario)) {
    scenario_rows <- df[df$Scenario == scenario_id, , drop = FALSE]
    for (i in seq_len(nrow(scenario_rows))) {
      row <- scenario_rows[i, ]
      scenario_cell <- ifelse(i == 1, as.character(scenario_id), "")
      stat_cells <- unlist(lapply(prior_labels, function(pr) {
        vapply(names(metric_labels), function(metric_name) {
          col_nm <- paste0(pr, "_", metric_labels[[metric_name]])
          format_number(row[[col_nm]], digits = metric_digits[[metric_name]])
        }, character(1))
      }))
      line <- paste(
        c(
          latex_escape(scenario_cell),
          row$Par,
          stat_cells
        ),
        collapse = " & "
      )
      lines <- c(lines, paste0(line, " \\\\"))
    }
    if (scenario_id != tail(unique(df$Scenario), 1)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.5ex}",
    "\\begin{minipage}{0.95\\linewidth}",
    note %||% paste0("\\footnotesize Notes: ", censoring_label, ". ESD = empirical SD; PSD = posterior SD; CP = coverage probability. Values summarize the 1000 simulation replicates within each log-logistic data-generation scenario."),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, path)
  message("Wrote LaTeX table: ", path)
}

write_csv <- function(df, path) {
  write.csv(df, path, row.names = FALSE)
  message("Wrote CSV: ", path)
}

clean_config <- function(config) {
  config <- sub("^JMbayes2_All_", "", config)
  config <- sub("^BP_", "", config)
  config
}

parse_dist_code <- function(code) {
  family <- ifelse(grepl("^LL", code), "Log-logistic",
                   ifelse(grepl("^W", code), "Weibull", NA_character_))
  shape_code <- sub("_.*$", "", code)
  shape_code <- sub("^[A-Z]+", "", shape_code)
  shape <- suppressWarnings(as.numeric(shape_code) / 100)
  scale <- suppressWarnings(as.numeric(sub("^.*_", "", code)) / 10)
  dist_short <- ifelse(
    family == "Log-logistic",
    sprintf("LL%.2f", shape),
    sprintf("W%.2f", shape)
  )
  dist_label <- sprintf("%s(%.2f, %.1f)", family, shape, scale)
  data.frame(
    family = family,
    baseline_shape = shape,
    baseline_scale = scale,
    dist_short = dist_short,
    dist_label = dist_label,
    stringsAsFactors = FALSE
  )
}

parse_config <- function(configs) {
  configs <- unique(configs)
  clean <- clean_config(configs)
  pattern <- "^k([0-9]+)_n([0-9]+)_cens([^_]+)_([A-Z]+[0-9]+_[0-9]+)_scen([0-9]+)_FU([0-9]+)$"
  matches <- regexec(pattern, clean)
  parts <- regmatches(clean, matches)
  bad <- lengths(parts) == 0
  if (any(bad)) {
    stop("Could not parse config(s): ", paste(configs[bad], collapse = ", "), call. = FALSE)
  }

  parsed <- lapply(seq_along(configs), function(i) {
    p <- parts[[i]]
    dist_info <- parse_dist_code(p[5])
    cens_raw <- p[4]
    scenario <- as.integer(p[6])
    data.frame(
      config = configs[i],
      config_clean = clean[i],
      k_bases = as.integer(p[2]),
      n = as.integer(p[3]),
      censoring = ifelse(cens_raw == "admin", "Admin", paste0(cens_raw, "%")),
      censoring_order = ifelse(cens_raw == "admin", 2, 1),
      dist_code = p[5],
      family = dist_info$family,
      baseline_shape = dist_info$baseline_shape,
      baseline_scale = dist_info$baseline_scale,
      dist_short = dist_info$dist_short,
      dist_label = dist_info$dist_label,
      scenario = scenario,
      fu = as.integer(p[7]),
      alpha_aft_truth = 0.012,
      beta_surv_aft_truth = c(0, 0.9, 0.9, -0.9, -0.9)[scenario],
      beta_long3_truth = c(0, 0.04, -0.04, 0.04, -0.04)[scenario],
      stringsAsFactors = FALSE
    )
  })

  bind_rows(parsed) %>%
    mutate(
      config_id = sprintf("%s | %s | S%d", dist_short, censoring, scenario),
      design_order = censoring_order * 1000 + match(dist_code, c("LL120_230", "W090_380", "W130_380", "W210_380")) * 10 + scenario
    ) %>%
    arrange(design_order)
}

parse_short_config <- function(configs) {
  configs <- as.character(configs)
  parsed <- lapply(configs, function(lbl) {
    pattern <- "^([A-Za-z0-9_() -]+)\\s+k=([0-9]+),\\s+c=([A-Za-z0-9]+)%?(?:,\\s+d=([^,]+))?,\\s+S([0-9]+)$"
    m <- regexec(pattern, lbl)
    p <- regmatches(lbl, m)[[1]]
    if (length(p) == 0) {
      return(data.frame(
        config_short = lbl,
        model_short = NA_character_,
        k_bases_short = NA_integer_,
        cens_short = NA_character_,
        dist_short_from_label = NA_character_,
        scenario_short = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      config_short = lbl,
      model_short = p[2],
      k_bases_short = as.integer(p[3]),
      cens_short = p[4],
      dist_short_from_label = p[5],
      scenario_short = as.integer(p[6]),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(parsed)
}

parameter_label <- function(variable) {
  dplyr::recode(
    variable,
    "alpha" = "alpha",
    "alpha_tilde" = "alpha_tilde",
    "beta_surv[1]" = "beta_surv[1]",
    "beta_long[1]" = "beta_long[1]",
    "beta_long[2]" = "beta_long[2]",
    "beta_long[3]" = "beta_long[3]",
    "sd_1_long[1]" = "sd_1_long[1]",
    "sd_1_long[2]" = "sd_1_long[2]",
    "sigma_long" = "sigma_long",
    .default = variable
  )
}

parameter_plot_label <- function(variable) {
  dplyr::recode(
    variable,
    "alpha" = "alpha",
    "beta_surv[1]" = "gamma",
    "beta_long[1]" = "beta[0]",
    "beta_long[2]" = "beta[1]",
    "beta_long[3]" = "beta[2]",
    "sd_1_long[1]" = "sigma[b0]",
    "sd_1_long[2]" = "sigma[b1]",
    "sigma_long" = "sigma[epsilon]",
    .default = variable
  )
}

summarise_config_rows <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      `N configs` = n(),
      `N bias-defined` = sum(!is.na(bias)),
      `Mean bias` = safe_mean(bias),
      `Median bias` = safe_median(bias),
      `Mean abs bias` = safe_mean(abs(bias)),
      `Mean RMSE` = safe_mean(rmse),
      `Mean coverage` = safe_mean(coverage),
      `Mean CI width` = safe_mean(mean_ci_width),
      `Mean post SD` = safe_mean(mean_post_sd),
      .groups = "drop"
    )
}

proposed_40 <- read_summary_csv(
  file.path(paths$proposed_40, "jm_parameter_summary_long_loglogistic_weibull.csv")
) %>%
  mutate(Method = "sAFT(BP)-JM")

jmbayes2_40 <- read_summary_csv(
  file.path(paths$jmbayes2_40, "jmbayes2_comparison_summary_loglogistic_weibull.csv")
) %>%
  mutate(Method = source)

config_meta <- parse_config(unique(c(proposed_40$config, jmbayes2_40$config)))

all_summary <- bind_rows(
  proposed_40 %>%
    select(any_of(names(jmbayes2_40)), everything()) %>%
    mutate(Method = "sAFT(BP)-JM"),
  jmbayes2_40
) %>%
  left_join(config_meta, by = "config") %>%
  mutate(
    Method = factor(
      Method,
      levels = c("sAFT(BP)-JM", "Cox-JM", "Weibull-JM", "loglogistic-JM", "LMM")
    ),
    Parameter = parameter_label(variable),
    PlotParameter = parameter_plot_label(variable)
  )

method_palette <- c(
  "sAFT(BP)-JM" = "#1B4D89",
  "Cox-JM" = "#C7472B",
  "Weibull-JM" = "#3D7D3A",
  "loglogistic-JM" = "#7B4AB8",
  "LMM" = "#555555"
)

## Table 1: simulation design.
simulation_design <- config_meta %>%
  filter(grepl("^BP_", config) | !duplicated(config_clean)) %>%
  distinct(config_clean, .keep_all = TRUE) %>%
  arrange(design_order) %>%
  transmute(
    Config = config_id,
    `Event-time distribution` = dist_label,
    `Sample size` = n,
    `Basis functions` = k_bases,
    Censoring = censoring,
    `Follow-up` = fu,
    Scenario = scenario,
    `True beta_long[3]` = beta_long3_truth,
    `True alpha (AFT)` = alpha_aft_truth,
    `True beta_surv[1] (AFT)` = beta_surv_aft_truth
  )

write_csv(simulation_design, file.path(csv_dir, "table1_simulation_design.csv"))
write_latex_table(
  simulation_design,
  file.path(table_dir, "table1_simulation_design.tex"),
  caption = "ENZAMET simulation configurations.",
  label = "tab:enzamet-simulation-design",
  digits = 3,
  column_digits = list(`True alpha (AFT)` = 3, `True beta_surv[1] (AFT)` = 1, `True beta_long[3]` = 2),
  align = "llrrlrrrrr"
)

## Table 2: proposed sAFT(BP)-JM key-parameter performance.
longitudinal_vars <- c(
  "beta_long[1]", "beta_long[2]", "beta_long[3]",
  "sd_1_long[1]", "sd_1_long[2]", "sigma_long"
)
key_vars <- c(longitudinal_vars, "alpha", "beta_surv[1]")
cm_labels <- c("CM1", "CM2")
cm_display_labels <- c("$\\text{CM}_1$", "$\\text{CM}_2$")

proposed_key_performance <- proposed_40 %>%
  filter(variable %in% key_vars) %>%
  mutate(Parameter = parameter_label(variable)) %>%
  summarise_config_rows(c("Parameter")) %>%
  arrange(match(Parameter, parameter_label(key_vars)))

write_csv(proposed_key_performance, file.path(csv_dir, "table2_saft_bp_key_performance.csv"))
write_latex_table(
  proposed_key_performance,
  file.path(table_dir, "table2_saft_bp_key_performance.tex"),
  caption = "Average key-parameter performance of the proposed sAFT(BP)-JM across the 40 ENZAMET configurations.",
  label = "tab:saft-bp-key-performance",
  digits = 4,
  column_digits = list(`Mean coverage` = 3, `Mean CI width` = 3, `Mean post SD` = 4),
  align = "lrrrrrrrrr"
)

saft_bp_distribution_labels <- c("LL1.20", "W0.90", "W1.30", "W2.10")
saft_bp_distribution_note <- paste0(
  "\\footnotesize Notes: $\\text{CM}_1$ = administrative censoring; ",
  "$\\text{CM}_2$ = controlled 50\\% censoring. ",
  "LL1.20 = log-logistic(1.20, 23.0); W0.90, W1.30, and W2.10 = Weibull(shape, 38.0). ",
  "Entries are means of posterior means across simulations."
)
saft_bp_distribution_detail_note <- function(cm_label) {
  paste0(
    "\\footnotesize Notes: ", cm_label, ". ",
    "LL1.20 = log-logistic(1.20, 23.0); W0.90, W1.30, and W2.10 = Weibull(shape, 38.0). ",
    "ESD = empirical SD; PSD = posterior SD; CP = coverage probability. ",
    "Values summarize the 1000 simulation replicates within each scenario and data-generation distribution."
  )
}

saft_bp_distribution_estimate_cols <- paste0(
  rep(saft_bp_distribution_labels, each = length(cm_labels)),
  "_",
  rep(cm_labels, times = length(saft_bp_distribution_labels))
)

proposed_distribution_source <- proposed_40 %>%
  left_join(config_meta, by = "config") %>%
  filter(variable %in% key_vars) %>%
  mutate(
    Scenario = scenario,
    Par = parameter_math_label(variable),
    Parameter = factor(variable, levels = key_vars),
    Distribution = factor(dist_short, levels = saft_bp_distribution_labels),
    CM = ifelse(censoring == "Admin", "CM1", "CM2")
  )

saft_bp_distribution_estimate_rows <- proposed_distribution_source %>%
  group_by(Scenario, Parameter, Par, truth, Distribution, CM) %>%
  summarise(Estimate = safe_mean(mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = c(Distribution, CM),
    values_from = Estimate,
    names_glue = "{Distribution}_{CM}"
  ) %>%
  arrange(Scenario, Parameter) %>%
  select(Scenario, Par, True = truth, dplyr::all_of(saft_bp_distribution_estimate_cols))

write_csv(
  saft_bp_distribution_estimate_rows,
  file.path(csv_dir, "table2_1_saft_bp_distribution_estimates.csv")
)
write_alpha_prior_estimate_table(
  saft_bp_distribution_estimate_rows,
  file.path(table_dir, "table2_1_saft_bp_distribution_estimates.tex"),
  caption = "sAFT(BP)-JM estimates by data-generation distribution.",
  label = "tab:saft-bp-distribution-estimates",
  prior_labels = saft_bp_distribution_labels,
  cm_labels = cm_labels,
  cm_display_labels = cm_display_labels,
  note = saft_bp_distribution_note
)

saft_lmm_distribution_model_labels <- c("sAFT", "LMM")
saft_lmm_distribution_model_display_labels <- c("sAFT", "LMM")
saft_lmm_distribution_vars <- c(
  "beta_long[1]", "beta_long[2]", "beta_long[3]",
  "alpha", "beta_surv[1]"
)
saft_lmm_distribution_estimate_cols <- unlist(lapply(
  saft_bp_distribution_labels,
  function(distribution_label) {
    unlist(lapply(cm_labels, function(cm_label) {
      paste(distribution_label, cm_label, saft_lmm_distribution_model_labels, sep = "_")
    }))
  }
))
saft_lmm_distribution_note <- paste0(
  "\\footnotesize Notes: $\\text{CM}_1$ = administrative censoring; ",
  "$\\text{CM}_2$ = controlled 50\\% censoring. ",
  "LL1.20 = log-logistic(1.20, 23.0); W0.90, W1.30, and W2.10 = Weibull(shape, 38.0). ",
  "sAFT denotes the proposed sAFT(BP)-JM. ",
  "Entries are means of point estimates across simulations; ",
  "dashes indicate parameters not fitted by the standalone LMM."
)

saft_lmm_distribution_source <- all_summary %>%
  filter(
    variable %in% saft_lmm_distribution_vars,
    as.character(Method) %in% c("sAFT(BP)-JM", "LMM"),
    Method == "sAFT(BP)-JM" | variable %in% longitudinal_vars
  ) %>%
  mutate(
    Model = ifelse(Method == "sAFT(BP)-JM", "sAFT", "LMM"),
    Model = factor(Model, levels = saft_lmm_distribution_model_labels),
    Scenario = scenario,
    Par = parameter_math_label(variable),
    Parameter = factor(variable, levels = saft_lmm_distribution_vars),
    Distribution = factor(dist_short, levels = saft_bp_distribution_labels),
    CM = ifelse(censoring == "Admin", "CM1", "CM2")
  )

saft_lmm_distribution_estimate_rows <- saft_lmm_distribution_source %>%
  group_by(Scenario, Parameter, Par, truth, Distribution, CM, Model) %>%
  summarise(Estimate = safe_mean(mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = c(Distribution, CM, Model),
    values_from = Estimate,
    names_glue = "{Distribution}_{CM}_{Model}"
  ) %>%
  arrange(Scenario, Parameter) %>%
  select(Scenario, Par, True = truth, dplyr::any_of(saft_lmm_distribution_estimate_cols))

missing_saft_lmm_distribution_cols <- setdiff(
  saft_lmm_distribution_estimate_cols,
  names(saft_lmm_distribution_estimate_rows)
)
if (length(missing_saft_lmm_distribution_cols) > 0) {
  saft_lmm_distribution_estimate_rows[missing_saft_lmm_distribution_cols] <- NA_real_
}
saft_lmm_distribution_estimate_rows <- saft_lmm_distribution_estimate_rows %>%
  select(Scenario, Par, True, dplyr::all_of(saft_lmm_distribution_estimate_cols))

write_csv(
  saft_lmm_distribution_estimate_rows,
  file.path(csv_dir, "table2_1_expanded_saft_lmm_distribution_estimates.csv")
)
write_distribution_model_estimate_table(
  saft_lmm_distribution_estimate_rows,
  file.path(table_dir, "table2_1_expanded_saft_lmm_distribution_estimates.tex"),
  caption = "sAFT(BP)-JM and LMM estimates by data-generation distribution.",
  label = "tab:saft-lmm-distribution-estimates",
  distribution_labels = saft_bp_distribution_labels,
  cm_labels = cm_labels,
  model_labels = saft_lmm_distribution_model_labels,
  cm_display_labels = cm_display_labels,
  model_display_labels = saft_lmm_distribution_model_display_labels,
  note = saft_lmm_distribution_note
)

saft_bp_distribution_detail_metrics <- c(
  Bias = "bias",
  ESD = "empirical_sd",
  PSD = "mean_post_sd",
  RMSE = "rmse",
  CP = "coverage"
)
saft_bp_distribution_detail_cols <- unlist(lapply(
  saft_bp_distribution_labels,
  function(distribution_label) paste0(distribution_label, "_", unname(saft_bp_distribution_detail_metrics))
))

make_saft_bp_distribution_detail_rows <- function(censoring_value) {
  proposed_distribution_source %>%
    filter(censoring == censoring_value) %>%
    group_by(Scenario, Parameter, Par, truth, Distribution) %>%
    summarise(
      mean = safe_mean(mean),
      bias = safe_mean(bias),
      empirical_sd = safe_mean(empirical_sd),
      mean_post_sd = safe_mean(mean_post_sd),
      rmse = safe_mean(rmse),
      coverage = safe_mean(coverage),
      mean_ci_width = safe_mean(mean_ci_width),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Distribution,
      values_from = c(mean, bias, empirical_sd, mean_post_sd, rmse, coverage, mean_ci_width),
      names_glue = "{Distribution}_{.value}"
    ) %>%
    arrange(Scenario, Parameter) %>%
    select(Scenario, Par, dplyr::all_of(saft_bp_distribution_detail_cols))
}

saft_bp_distribution_detail_admin <- make_saft_bp_distribution_detail_rows("Admin")
saft_bp_distribution_detail_cens50 <- make_saft_bp_distribution_detail_rows("50%")

write_csv(
  saft_bp_distribution_detail_admin,
  file.path(csv_dir, "table2_2_saft_bp_distribution_details_CM1.csv")
)
write_alpha_prior_detail_table(
  saft_bp_distribution_detail_admin,
  file.path(table_dir, "table2_2_saft_bp_distribution_details_CM1.tex"),
  caption = "Detailed sAFT(BP)-JM performance by data-generation distribution under administrative censoring.",
  label = "tab:saft-bp-distribution-details-cm1",
  prior_labels = saft_bp_distribution_labels,
  metric_labels = saft_bp_distribution_detail_metrics,
  censoring_label = "$\\text{CM}_1$ = administrative censoring",
  note = saft_bp_distribution_detail_note("$\\text{CM}_1$ = administrative censoring")
)

write_csv(
  saft_bp_distribution_detail_cens50,
  file.path(csv_dir, "table2_3_saft_bp_distribution_details_CM2.csv")
)
write_alpha_prior_detail_table(
  saft_bp_distribution_detail_cens50,
  file.path(table_dir, "table2_3_saft_bp_distribution_details_CM2.tex"),
  caption = "Detailed sAFT(BP)-JM performance by data-generation distribution under controlled 50 percent censoring.",
  label = "tab:saft-bp-distribution-details-cm2",
  prior_labels = saft_bp_distribution_labels,
  metric_labels = saft_bp_distribution_detail_metrics,
  censoring_label = "$\\text{CM}_2$ = controlled 50 percent censoring",
  note = saft_bp_distribution_detail_note("$\\text{CM}_2$ = controlled 50 percent censoring")
)

## Table 3: longitudinal comparison across all fitted models.
longitudinal_comparison <- all_summary %>%
  filter(variable %in% longitudinal_vars) %>%
  mutate(Parameter = parameter_label(variable)) %>%
  summarise_config_rows(c("Method", "Parameter")) %>%
  arrange(Method, match(Parameter, parameter_label(longitudinal_vars)))

write_csv(longitudinal_comparison, file.path(csv_dir, "table3_longitudinal_model_comparison.csv"))
write_latex_table(
  longitudinal_comparison,
  file.path(table_dir, "table3_longitudinal_model_comparison.tex"),
  caption = "Longitudinal-parameter performance compared across the proposed model, JMbayes2 joint models, and LMM.",
  label = "tab:longitudinal-model-comparison",
  digits = 4,
  column_digits = list(`Mean coverage` = 3, `Mean CI width` = 3, `Mean post SD` = 4),
  align = "llrrrrrrrrr"
)

## Screenshot-style compact table: proposed BP joint model and LMM under
## log-logistic data generation. Survival-side parameters are not fitted by
## the standalone LMM, so those cells are intentionally left as dashes.
bp_lmm_loglogistic_vars <- c(
  "beta_long[1]", "beta_long[2]", "beta_long[3]",
  "alpha", "beta_surv[1]"
)
bp_lmm_loglogistic_labels <- c("BP", "LMM")
bp_lmm_loglogistic_cols <- paste0(
  rep(bp_lmm_loglogistic_labels, each = length(cm_labels)),
  "_",
  rep(cm_labels, times = length(bp_lmm_loglogistic_labels))
)
bp_lmm_parameter_math_label <- function(variable) {
  out <- parameter_math_label(variable)
  out[variable == "beta_surv[1]"] <- "$\\gamma$"
  out
}

bp_lmm_loglogistic_source <- all_summary %>%
  filter(
    family == "Log-logistic",
    variable %in% bp_lmm_loglogistic_vars,
    Method %in% c("sAFT(BP)-JM", "LMM"),
    Method == "sAFT(BP)-JM" | variable %in% longitudinal_vars
  ) %>%
  mutate(
    Model = ifelse(Method == "sAFT(BP)-JM", "BP", "LMM"),
    Model = factor(Model, levels = bp_lmm_loglogistic_labels),
    Scenario = scenario,
    Par = bp_lmm_parameter_math_label(variable),
    Parameter = factor(variable, levels = bp_lmm_loglogistic_vars),
    CM = ifelse(censoring == "Admin", "CM1", "CM2")
  )

bp_lmm_loglogistic_estimate_rows <- bp_lmm_loglogistic_source %>%
  group_by(Scenario, Parameter, Par, truth, Model, CM) %>%
  summarise(Estimate = safe_mean(mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = c(Model, CM),
    values_from = Estimate,
    names_glue = "{Model}_{CM}"
  ) %>%
  arrange(Scenario, Parameter) %>%
  select(Scenario, Par, True = truth, dplyr::all_of(bp_lmm_loglogistic_cols))

write_csv(
  bp_lmm_loglogistic_estimate_rows,
  file.path(csv_dir, "table3_0_loglogistic_bp_lmm_estimates.csv")
)
write_alpha_prior_estimate_table(
  bp_lmm_loglogistic_estimate_rows,
  file.path(table_dir, "table3_0_loglogistic_bp_lmm_estimates.tex"),
  caption = "Simulation results under log-logistic data generation for the proposed BP joint model and LMM.",
  label = "tab:loglogistic-bp-lmm-estimates",
  prior_labels = bp_lmm_loglogistic_labels,
  cm_labels = cm_labels,
  cm_display_labels = cm_display_labels,
  super_header = "Log-logistic baseline",
  na_text = "--",
  note = paste0(
    "\\footnotesize Notes: $\\text{CM}_1$ = administrative censoring; ",
    "$\\text{CM}_2$ = controlled 50\\% censoring. ",
    "BP denotes the proposed sAFT(BP)-JM. ",
    "Entries are means of point estimates across the 1000 simulations in each scenario; ",
    "dashes indicate survival-side parameters not fitted by the standalone LMM."
  )
)

longitudinal_ll_model_labels <- c("LMM", "sAFT(BP)-JM", "Cox-JM")
longitudinal_ll_note <- paste0(
  "\\footnotesize Notes: $\\text{CM}_1$ = administrative censoring; ",
  "$\\text{CM}_2$ = controlled 50\\% censoring. ",
  "Only log-logistic data-generation configurations are included. ",
  "Entries are means of point estimates across simulations."
)
longitudinal_ll_detail_note <- function(cm_label) {
  paste0(
    "\\footnotesize Notes: ", cm_label, ". ",
    "Only log-logistic data-generation configurations are included. ",
    "ESD = empirical SD; PSD = posterior SD; CP = coverage probability. ",
    "Values summarize the 1000 simulation replicates within each scenario and fitted model. ",
    "For LMM variance parameters, unavailable interval-based Bayesian summaries remain NA."
  )
}

longitudinal_ll_estimate_cols <- paste0(
  rep(longitudinal_ll_model_labels, each = length(cm_labels)),
  "_",
  rep(cm_labels, times = length(longitudinal_ll_model_labels))
)

longitudinal_ll_source <- all_summary %>%
  filter(
    family == "Log-logistic",
    variable %in% longitudinal_vars,
    as.character(Method) %in% longitudinal_ll_model_labels
  ) %>%
  mutate(
    Model = factor(as.character(Method), levels = longitudinal_ll_model_labels),
    Scenario = scenario,
    Par = parameter_math_label(variable),
    Parameter = factor(variable, levels = longitudinal_vars),
    CM = ifelse(censoring == "Admin", "CM1", "CM2")
  )

longitudinal_ll_estimate_rows <- longitudinal_ll_source %>%
  group_by(Scenario, Parameter, Par, truth, Model, CM) %>%
  summarise(Estimate = safe_mean(mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = c(Model, CM),
    values_from = Estimate,
    names_glue = "{Model}_{CM}"
  ) %>%
  arrange(Scenario, Parameter) %>%
  select(Scenario, Par, True = truth, dplyr::all_of(longitudinal_ll_estimate_cols))

write_csv(
  longitudinal_ll_estimate_rows,
  file.path(csv_dir, "table3_1_longitudinal_LL_model_estimates.csv")
)
write_alpha_prior_estimate_table(
  longitudinal_ll_estimate_rows,
  file.path(table_dir, "table3_1_longitudinal_LL_model_estimates.tex"),
  caption = "Longitudinal estimates under log-logistic data generation.",
  label = "tab:longitudinal-ll-model-estimates",
  prior_labels = longitudinal_ll_model_labels,
  cm_labels = cm_labels,
  cm_display_labels = cm_display_labels,
  note = longitudinal_ll_note
)

longitudinal_ll_detail_metrics <- c(
  Bias = "bias",
  ESD = "empirical_sd",
  PSD = "mean_post_sd",
  RMSE = "rmse",
  CP = "coverage"
)
longitudinal_ll_detail_cols <- unlist(lapply(
  longitudinal_ll_model_labels,
  function(model_label) paste0(model_label, "_", unname(longitudinal_ll_detail_metrics))
))

make_longitudinal_ll_detail_rows <- function(censoring_value) {
  longitudinal_ll_source %>%
    filter(censoring == censoring_value) %>%
    group_by(Scenario, Parameter, Par, truth, Model) %>%
    summarise(
      mean = safe_mean(mean),
      bias = safe_mean(bias),
      empirical_sd = safe_mean(empirical_sd),
      mean_post_sd = safe_mean(mean_post_sd),
      rmse = safe_mean(rmse),
      coverage = safe_mean(coverage),
      mean_ci_width = safe_mean(mean_ci_width),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Model,
      values_from = c(mean, bias, empirical_sd, mean_post_sd, rmse, coverage, mean_ci_width),
      names_glue = "{Model}_{.value}"
    ) %>%
    arrange(Scenario, Parameter) %>%
    select(Scenario, Par, dplyr::all_of(longitudinal_ll_detail_cols))
}

longitudinal_ll_detail_admin <- make_longitudinal_ll_detail_rows("Admin")
longitudinal_ll_detail_cens50 <- make_longitudinal_ll_detail_rows("50%")

write_csv(
  longitudinal_ll_detail_admin,
  file.path(csv_dir, "table3_2_longitudinal_LL_model_details_CM1.csv")
)
write_alpha_prior_detail_table(
  longitudinal_ll_detail_admin,
  file.path(table_dir, "table3_2_longitudinal_LL_model_details_CM1.tex"),
  caption = "Detailed longitudinal performance under log-logistic data generation and administrative censoring.",
  label = "tab:longitudinal-ll-model-details-cm1",
  prior_labels = longitudinal_ll_model_labels,
  metric_labels = longitudinal_ll_detail_metrics,
  censoring_label = "$\\text{CM}_1$ = administrative censoring",
  note = longitudinal_ll_detail_note("$\\text{CM}_1$ = administrative censoring")
)

write_csv(
  longitudinal_ll_detail_cens50,
  file.path(csv_dir, "table3_3_longitudinal_LL_model_details_CM2.csv")
)
write_alpha_prior_detail_table(
  longitudinal_ll_detail_cens50,
  file.path(table_dir, "table3_3_longitudinal_LL_model_details_CM2.tex"),
  caption = "Detailed longitudinal performance under log-logistic data generation and controlled 50 percent censoring.",
  label = "tab:longitudinal-ll-model-details-cm2",
  prior_labels = longitudinal_ll_model_labels,
  metric_labels = longitudinal_ll_detail_metrics,
  censoring_label = "$\\text{CM}_2$ = controlled 50 percent censoring",
  note = longitudinal_ll_detail_note("$\\text{CM}_2$ = controlled 50 percent censoring")
)

## Table 4: alpha bias comparison.
alpha_summary_source <- all_summary %>%
  filter(variable == "alpha") %>%
  mutate(
    `Config group` = dist_label,
    `Truth scale` = ifelse(
      Method == "sAFT(BP)-JM",
      "AFT",
      ifelse(family == "Weibull", "converted PH", "undefined for LL DGP")
    )
  )

alpha_bias_overall <- alpha_summary_source %>%
  mutate(
    `Config group` = "Overall",
    `Truth scale` = ifelse(
      Method == "sAFT(BP)-JM",
      "AFT",
      "converted PH where defined"
    )
  ) %>%
  summarise_config_rows(c("Method", "Config group", "Truth scale"))

alpha_bias_by_distribution <- alpha_summary_source %>%
  summarise_config_rows(c("Method", "Config group", "Truth scale"))

alpha_bias_comparison <- bind_rows(alpha_bias_overall, alpha_bias_by_distribution) %>%
  arrange(Method, factor(`Config group`, levels = c("Overall", sort(unique(alpha_summary_source$`Config group`)))))

write_csv(alpha_bias_comparison, file.path(csv_dir, "table4_alpha_bias_model_comparison.csv"))
write_latex_table(
  alpha_bias_comparison,
  file.path(table_dir, "table4_alpha_bias_model_comparison.tex"),
  caption = "Alpha-bias comparison across fitted joint models. JMbayes2 alpha truth is the converted PH-scale truth for Weibull-generated data and is undefined for log-logistic-generated data.",
  label = "tab:alpha-bias-model-comparison",
  digits = 4,
  column_digits = list(`Mean coverage` = 3, `Mean CI width` = 4, `Mean post SD` = 4),
  align = "lllrrrrrrrrr"
)

## Table 5: alpha-prior sensitivity, LL-only configs.
read_alpha_sensitivity <- function(folder, maxF, label, file_name) {
  df <- read_summary_csv(file.path(folder, file_name)) %>%
    left_join(parse_config(unique(.$config)), by = "config") %>%
    filter(family == "Log-logistic") %>%
    mutate(
      `Alpha prior maxF` = maxF,
      `Prior setting` = label,
      Method = "sAFT(BP)-JM",
      Parameter = parameter_label(variable)
    )
  df
}

alpha_sensitivity <- bind_rows(
  read_alpha_sensitivity(paths$alpha_maxF125, 125, "maxF125", "jm_parameter_summary_long_loglogistic.csv"),
  read_alpha_sensitivity(paths$alpha_maxF150, 150, "maxF150", "jm_parameter_summary_long_loglogistic.csv"),
  read_alpha_sensitivity(paths$alpha_maxF200, 200, "maxF200 baseline", "jm_parameter_summary_long_loglogistic_weibull.csv"),
  read_alpha_sensitivity(paths$alpha_maxF300, 300, "maxF300", "jm_parameter_summary_long_loglogistic.csv")
)

alpha_prior_vars <- c(
  longitudinal_vars,
  "alpha", "beta_surv[1]"
)
alpha_prior_labels <- c("maxF125", "maxF150", "maxF200", "maxF300")
alpha_prior_sd_values <- log(c(1.25, 1.5, 2, 3)) / 1.96
names(alpha_prior_sd_values) <- alpha_prior_labels
alpha_prior_sd_labels <- format_number(alpha_prior_sd_values, digits = 3)
alpha_prior_display_labels <- paste0(
  "$\\sigma_{\\alpha}=", alpha_prior_sd_labels, "$"
)
alpha_prior_note <- paste0(
  "\\footnotesize Notes: $\\sigma_{\\alpha}$ denotes the prior SD for $\\alpha$, ",
  "with $\\sigma_{\\alpha}=\\log(F_{\\max})/1.96$ for ",
  "$F_{\\max}=1.25,1.5,2,3$. $\\text{CM}_1$ = administrative censoring; ",
  "$\\text{CM}_2$ = controlled 50\\% censoring. Entries are means of posterior means across simulations."
)
alpha_prior_detail_note <- function(cm_label) {
  paste0(
    "\\footnotesize Notes: ", cm_label, ". ",
    "$\\sigma_{\\alpha}$ denotes the prior SD for $\\alpha$, ",
    "with $\\sigma_{\\alpha}=\\log(F_{\\max})/1.96$ for ",
    "$F_{\\max}=1.25,1.5,2,3$. ",
    "ESD = empirical SD; PSD = posterior SD; CP = coverage probability. ",
    "Values summarize the 1000 simulation replicates within each log-logistic data-generation scenario."
  )
}
alpha_prior_censor_labels <- cm_labels
alpha_prior_censor_display_labels <- cm_display_labels

alpha_prior_sensitivity_table <- alpha_sensitivity %>%
  filter(variable %in% alpha_prior_vars) %>%
  summarise_config_rows(c("Alpha prior maxF", "Prior setting", "Parameter")) %>%
  mutate(
    `Alpha prior SD` = alpha_prior_sd_values[paste0("maxF", `Alpha prior maxF`)]
  ) %>%
  arrange(`Alpha prior maxF`, match(Parameter, alpha_prior_vars)) %>%
  select(`Alpha prior SD`, Parameter, everything(), -`Alpha prior maxF`, -`Prior setting`)

write_csv(alpha_prior_sensitivity_table, file.path(csv_dir, "table5_alpha_prior_sensitivity_LL.csv"))
write_latex_table(
  alpha_prior_sensitivity_table,
  file.path(table_dir, "table5_alpha_prior_sensitivity_LL.tex"),
  caption = "Alpha-prior sensitivity for survival-association and longitudinal parameters in the 10 log-logistic ENZAMET configurations. The prior SD is sigma_alpha = log(F_max) / 1.96.",
  label = "tab:alpha-prior-sensitivity",
  digits = 4,
  column_digits = list(`Alpha prior SD` = 3, `Mean coverage` = 3, `Mean CI width` = 4, `Mean post SD` = 4),
  align = "rlrrrrrrrrr"
)

alpha_prior_estimate_cols <- paste0(
  rep(alpha_prior_labels, each = length(alpha_prior_censor_labels)),
  "_",
  rep(alpha_prior_censor_labels, times = length(alpha_prior_labels))
)

alpha_prior_estimate_rows <- alpha_sensitivity %>%
  filter(variable %in% alpha_prior_vars) %>%
  mutate(
    Scenario = scenario,
    Par = parameter_math_label(variable),
    Parameter = factor(variable, levels = alpha_prior_vars),
    Prior = paste0("maxF", `Alpha prior maxF`),
    CM = ifelse(censoring == "Admin", "CM1", "CM2")
  ) %>%
  group_by(Scenario, Parameter, Par, truth, Prior, CM) %>%
  summarise(Estimate = safe_mean(mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = c(Prior, CM),
    values_from = Estimate,
    names_glue = "{Prior}_{CM}"
  ) %>%
  arrange(Scenario, Parameter) %>%
  select(Scenario, Par, True = truth, dplyr::all_of(alpha_prior_estimate_cols))

write_csv(alpha_prior_estimate_rows, file.path(csv_dir, "table5_1_alpha_prior_sensitivity_LL_estimates.csv"))
write_alpha_prior_estimate_table(
  alpha_prior_estimate_rows,
  file.path(table_dir, "table5_1_alpha_prior_sensitivity_LL_estimates.tex"),
  caption = "Alpha-prior sensitivity estimates in the log-logistic data-generation settings.",
  label = "tab:alpha-prior-sensitivity-estimates",
  prior_labels = alpha_prior_labels,
  prior_display_labels = alpha_prior_display_labels,
  cm_labels = alpha_prior_censor_labels,
  cm_display_labels = alpha_prior_censor_display_labels,
  note = alpha_prior_note
)

alpha_prior_detail_metrics <- c(
  Bias = "bias",
  ESD = "empirical_sd",
  PSD = "mean_post_sd",
  RMSE = "rmse",
  CP = "coverage"
)
alpha_prior_detail_cols <- unlist(lapply(
  alpha_prior_labels,
  function(prior_label) paste0(prior_label, "_", unname(alpha_prior_detail_metrics))
))

make_alpha_prior_detail_rows <- function(censoring_value) {
  alpha_sensitivity %>%
    filter(variable %in% alpha_prior_vars, censoring == censoring_value) %>%
    mutate(
      Scenario = scenario,
      Par = parameter_math_label(variable),
      Parameter = factor(variable, levels = alpha_prior_vars),
      Prior = paste0("maxF", `Alpha prior maxF`)
    ) %>%
    group_by(Scenario, Parameter, Par, truth, Prior) %>%
    summarise(
      mean = safe_mean(mean),
      bias = safe_mean(bias),
      empirical_sd = safe_mean(empirical_sd),
      mean_post_sd = safe_mean(mean_post_sd),
      rmse = safe_mean(rmse),
      coverage = safe_mean(coverage),
      mean_ci_width = safe_mean(mean_ci_width),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Prior,
      values_from = c(mean, bias, empirical_sd, mean_post_sd, rmse, coverage, mean_ci_width),
      names_glue = "{Prior}_{.value}"
    ) %>%
    arrange(Scenario, Parameter) %>%
    select(Scenario, Par, dplyr::all_of(alpha_prior_detail_cols))
}

alpha_prior_detail_admin <- make_alpha_prior_detail_rows("Admin")
alpha_prior_detail_cens50 <- make_alpha_prior_detail_rows("50%")

write_csv(alpha_prior_detail_admin, file.path(csv_dir, "table5_2_alpha_prior_sensitivity_LL_details_CM1.csv"))
write_alpha_prior_detail_table(
  alpha_prior_detail_admin,
  file.path(table_dir, "table5_2_alpha_prior_sensitivity_LL_details_CM1.tex"),
  caption = "Detailed alpha-prior sensitivity statistics for log-logistic generated data under administrative censoring.",
  label = "tab:alpha-prior-sensitivity-details-cm1",
  prior_labels = alpha_prior_labels,
  prior_display_labels = alpha_prior_display_labels,
  metric_labels = alpha_prior_detail_metrics,
  censoring_label = "$\\text{CM}_1$ = administrative censoring",
  note = alpha_prior_detail_note("$\\text{CM}_1$ = administrative censoring")
)

write_csv(alpha_prior_detail_cens50, file.path(csv_dir, "table5_3_alpha_prior_sensitivity_LL_details_CM2.csv"))
write_alpha_prior_detail_table(
  alpha_prior_detail_cens50,
  file.path(table_dir, "table5_3_alpha_prior_sensitivity_LL_details_CM2.tex"),
  caption = "Detailed alpha-prior sensitivity statistics for log-logistic generated data under controlled 50 percent censoring.",
  label = "tab:alpha-prior-sensitivity-details-cm2",
  prior_labels = alpha_prior_labels,
  prior_display_labels = alpha_prior_display_labels,
  metric_labels = alpha_prior_detail_metrics,
  censoring_label = "$\\text{CM}_2$ = controlled 50 percent censoring",
  note = alpha_prior_detail_note("$\\text{CM}_2$ = controlled 50 percent censoring")
)

## Figures.
save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

## Replacement for the older five separate scenario figures in the draft:
## one matrix of comparative per-simulation boxplots under log-logistic data
## generation. This intentionally reads compact per-simulation summaries.
proposed_per_sim_path <- file.path(
  paths$proposed_40,
  "jm_posterior_summaries_per_sim_loglogistic_weibull.rds"
)
jmbayes2_per_sim_path <- file.path(
  paths$jmbayes2_40,
  "jmbayes2_comparison_per_sim_loglogistic_weibull.rds"
)
if (file.exists(proposed_per_sim_path) && file.exists(jmbayes2_per_sim_path)) {
  proposed_per_sim <- readRDS(proposed_per_sim_path)
  jmbayes2_per_sim <- readRDS(jmbayes2_per_sim_path)
  boxplot_vars <- c(
    "beta_long[1]", "beta_long[2]", "beta_long[3]",
    "alpha", "beta_surv[1]"
  )
  boxplot_plot_labels <- parameter_plot_label(boxplot_vars)
  boxplot_plot_labels[boxplot_vars == "beta_surv[1]"] <- "gamma"
  boxplot_parameter_levels <- unname(boxplot_plot_labels)

  boxplot_truth <- expand.grid(
    Scenario = 1:5,
    variable = boxplot_vars,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      truth = case_when(
        variable == "beta_long[1]" ~ 73,
        variable == "beta_long[2]" ~ -0.04,
        variable == "beta_long[3]" ~ c(0, 0.04, -0.04, 0.04, -0.04)[Scenario],
        variable == "alpha" ~ 0.012,
        variable == "beta_surv[1]" ~ c(0, 0.9, 0.9, -0.9, -0.9)[Scenario],
        TRUE ~ NA_real_
      ),
      PlotParameter = parameter_plot_label(variable),
      PlotParameter = ifelse(variable == "beta_surv[1]", "gamma", PlotParameter),
      PlotParameter = factor(PlotParameter, levels = boxplot_parameter_levels),
      ScenarioLabel = factor(paste0("S", Scenario), levels = paste0("S", 1:5))
    )

  proposed_boxplot <- proposed_per_sim %>%
    filter(variable %in% boxplot_vars) %>%
    left_join(parse_short_config(unique(.$config)), by = c("config" = "config_short")) %>%
    filter(
      dist_short_from_label == "LL1.20",
      scenario_short %in% 1:5,
      cens_short %in% c("admin", "50")
    ) %>%
    transmute(
      variable,
      estimate = mean,
      Scenario = scenario_short,
      Censoring = ifelse(cens_short == "admin", "CM1", "CM2"),
      Method = "sAFT-JM",
      Estimator = paste("sAFT-JM", Censoring)
    )

  lmm_boxplot <- jmbayes2_per_sim %>%
    filter(source == "LMM", variable %in% boxplot_vars) %>%
    left_join(parse_short_config(unique(.$cfg_short)), by = c("cfg_short" = "config_short")) %>%
    filter(
      dist_short_from_label == "LL1.20",
      scenario_short %in% 1:5,
      cens_short %in% c("admin", "50")
    ) %>%
    transmute(
      variable,
      estimate = mean,
      Scenario = scenario_short,
      Censoring = ifelse(cens_short == "admin", "CM1", "CM2"),
      Method = "LMM",
      Estimator = paste("LMM", Censoring)
    )

  scenario_boxplot_data <- bind_rows(proposed_boxplot, lmm_boxplot) %>%
    mutate(
      Method = factor(Method, levels = c("sAFT-JM", "LMM")),
      Estimator = factor(
        Estimator,
        levels = c("sAFT-JM CM1", "sAFT-JM CM2", "LMM CM1", "LMM CM2")
      ),
      PlotParameter = parameter_plot_label(variable),
      PlotParameter = ifelse(variable == "beta_surv[1]", "gamma", PlotParameter),
      PlotParameter = factor(PlotParameter, levels = boxplot_parameter_levels),
      ScenarioLabel = factor(paste0("S", Scenario), levels = paste0("S", 1:5))
    )

  boxplot_fill_values <- c(
    "sAFT-JM CM1" = "#A6CEE3",
    "sAFT-JM CM2" = "#1F78B4",
    "LMM CM1" = "#FDB462",
    "LMM CM2" = "#D95F02"
  )
  boxplot_linetype_values <- c("sAFT-JM" = "solid", "LMM" = "solid")

  fig_scenario_coefficient_boxplot_matrix <- ggplot(
    scenario_boxplot_data,
    aes(x = Estimator, y = estimate, fill = Estimator, linetype = Method)
  ) +
    geom_boxplot(
      outlier.size = 0.15,
      linewidth = 0.28,
      width = 0.7,
      alpha = 0.85
    ) +
    geom_hline(
      data = boxplot_truth,
      aes(yintercept = truth),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.25
    ) +
    facet_grid(
      PlotParameter ~ ScenarioLabel,
      scales = "free_y",
      labeller = labeller(PlotParameter = label_parsed)
    ) +
    scale_fill_manual(
      name = "Configuration",
      values = boxplot_fill_values,
      labels = names(boxplot_fill_values),
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = "Model",
      values = boxplot_linetype_values,
      labels = names(boxplot_linetype_values),
      drop = FALSE,
      guide = "none"
    ) +
    guides(
      fill = guide_legend(order = 1)
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 11),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      x = NULL,
      y = "Estimate",
      title = "Simulation estimates under log-logistic data generation",
      subtitle = "Rows are coefficients and columns are scenarios; dashed horizontal lines indicate true values."
    )

  save_figure(
    fig_scenario_coefficient_boxplot_matrix,
    "fig_scenario_coefficient_boxplot_matrix_loglogistic",
    13,
    10
  )

  boxplot_truth_all40 <- expand.grid(
    Distribution = saft_bp_distribution_labels,
    Scenario = 1:5,
    variable = boxplot_vars,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      truth = case_when(
        variable == "beta_long[1]" ~ 73,
        variable == "beta_long[2]" ~ -0.04,
        variable == "beta_long[3]" ~ c(0, 0.04, -0.04, 0.04, -0.04)[Scenario],
        variable == "alpha" ~ 0.012,
        variable == "beta_surv[1]" ~ c(0, 0.9, 0.9, -0.9, -0.9)[Scenario],
        TRUE ~ NA_real_
      ),
      Distribution = factor(Distribution, levels = saft_bp_distribution_labels),
      PlotParameter = parameter_plot_label(variable),
      PlotParameter = ifelse(variable == "beta_surv[1]", "gamma", PlotParameter),
      PlotParameter = factor(PlotParameter, levels = boxplot_parameter_levels),
      ScenarioLabel = factor(paste0("S", Scenario), levels = paste0("S", 1:5))
    )

  proposed_boxplot_all40 <- proposed_per_sim %>%
    filter(variable %in% boxplot_vars) %>%
    left_join(parse_short_config(unique(.$config)), by = c("config" = "config_short")) %>%
    filter(
      dist_short_from_label %in% saft_bp_distribution_labels,
      scenario_short %in% 1:5,
      cens_short %in% c("admin", "50")
    ) %>%
    transmute(
      variable,
      estimate = mean,
      Distribution = factor(dist_short_from_label, levels = saft_bp_distribution_labels),
      Scenario = scenario_short,
      Censoring = ifelse(cens_short == "admin", "CM1", "CM2"),
      Method = "sAFT-JM",
      Estimator = paste("sAFT-JM", Censoring)
    )

  lmm_boxplot_all40 <- jmbayes2_per_sim %>%
    filter(source == "LMM", variable %in% boxplot_vars) %>%
    left_join(parse_short_config(unique(.$cfg_short)), by = c("cfg_short" = "config_short")) %>%
    filter(
      dist_short_from_label %in% saft_bp_distribution_labels,
      scenario_short %in% 1:5,
      cens_short %in% c("admin", "50")
    ) %>%
    transmute(
      variable,
      estimate = mean,
      Distribution = factor(dist_short_from_label, levels = saft_bp_distribution_labels),
      Scenario = scenario_short,
      Censoring = ifelse(cens_short == "admin", "CM1", "CM2"),
      Method = "LMM",
      Estimator = paste("LMM", Censoring)
    )

  scenario_boxplot_data_all40 <- bind_rows(proposed_boxplot_all40, lmm_boxplot_all40) %>%
    mutate(
      Method = factor(Method, levels = c("sAFT-JM", "LMM")),
      Estimator = factor(
        Estimator,
        levels = c("sAFT-JM CM1", "sAFT-JM CM2", "LMM CM1", "LMM CM2")
      ),
      PlotParameter = parameter_plot_label(variable),
      PlotParameter = ifelse(variable == "beta_surv[1]", "gamma", PlotParameter),
      PlotParameter = factor(PlotParameter, levels = boxplot_parameter_levels),
      ScenarioLabel = factor(paste0("S", Scenario), levels = paste0("S", 1:5))
    )

  boxplot_truth_all40_panels <- boxplot_truth_all40 %>%
    distinct(Scenario, variable, truth, PlotParameter, ScenarioLabel)

  fig_scenario_coefficient_boxplot_matrix_all40 <- ggplot(
    scenario_boxplot_data_all40,
    aes(
      x = Distribution,
      y = estimate,
      fill = Estimator,
      linetype = Method,
      group = interaction(Distribution, Estimator)
    )
  ) +
    geom_boxplot(
      outlier.size = 0.08,
      linewidth = 0.22,
      width = 0.72,
      alpha = 0.86,
      position = position_dodge2(width = 0.88, preserve = "single", padding = 0.08)
    ) +
    geom_hline(
      data = boxplot_truth_all40_panels,
      aes(yintercept = truth),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.22
    ) +
    facet_grid(
      PlotParameter ~ ScenarioLabel,
      scales = "free_y",
      labeller = labeller(PlotParameter = label_parsed)
    ) +
    scale_fill_manual(
      name = "Model/censoring",
      values = boxplot_fill_values,
      labels = names(boxplot_fill_values),
      drop = FALSE
    ) +
    scale_linetype_manual(
      name = "Model",
      values = boxplot_linetype_values,
      labels = names(boxplot_linetype_values),
      drop = FALSE,
      guide = "none"
    ) +
    guides(
      fill = guide_legend(order = 1)
    ) +
    theme_bw(base_size = 7.2) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.2),
      axis.ticks.x = element_line(linewidth = 0.2),
      strip.text.x = element_text(size = 7.2),
      strip.text.y = element_text(size = 6.8),
      panel.grid.minor = element_blank(),
      panel.spacing.x = grid::unit(0.45, "lines"),
      panel.spacing.y = grid::unit(0.35, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 7.2),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8.2)
    ) +
    labs(
      x = "Data-generation distribution",
      y = "Estimate",
      title = "Simulation estimates across data-generation distributions",
      subtitle = "Each panel compares four data-generation distributions; dashed horizontal lines indicate true values."
    )

  save_figure(
    fig_scenario_coefficient_boxplot_matrix_all40,
    "fig_scenario_coefficient_boxplot_matrix_all40",
    14,
    10.5
  )

  slide_treatment_vars <- c("beta_long[3]", "beta_surv[1]")
  slide_treatment_parameter_levels <- parameter_plot_label(slide_treatment_vars)
  slide_treatment_parameter_levels[slide_treatment_vars == "beta_surv[1]"] <- "gamma"

  slide_treatment_truth <- boxplot_truth_all40 %>%
    filter(variable %in% slide_treatment_vars) %>%
    distinct(Scenario, variable, truth, PlotParameter, ScenarioLabel)

  slide_treatment_data <- proposed_boxplot_all40 %>%
    filter(variable %in% slide_treatment_vars) %>%
    mutate(
      Censoring = factor(Censoring, levels = cm_labels),
      PlotParameter = parameter_plot_label(variable),
      PlotParameter = ifelse(variable == "beta_surv[1]", "gamma", PlotParameter),
      PlotParameter = factor(PlotParameter, levels = slide_treatment_parameter_levels),
      ScenarioLabel = factor(paste0("S", Scenario), levels = paste0("S", 1:5))
    )

  fig_slide_saft_treatment_effect_recovery_all40 <- ggplot(
    slide_treatment_data,
    aes(
      x = Distribution,
      y = estimate,
      fill = Censoring,
      group = interaction(Distribution, Censoring)
    )
  ) +
    geom_boxplot(
      outlier.size = 0.08,
      linewidth = 0.32,
      width = 0.68,
      alpha = 0.9,
      position = position_dodge2(width = 0.72, preserve = "single", padding = 0.08)
    ) +
    geom_hline(
      data = slide_treatment_truth,
      aes(yintercept = truth),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      linewidth = 0.28
    ) +
    facet_grid(
      PlotParameter ~ ScenarioLabel,
      scales = "free_y",
      labeller = labeller(PlotParameter = label_parsed)
    ) +
    scale_fill_manual(
      name = "Censoring",
      values = c("CM1" = "#A6CEE3", "CM2" = "#1F78B4"),
      labels = cm_labels,
      drop = FALSE
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8.2),
      axis.text.y = element_text(size = 8.6),
      axis.title = element_text(size = 10),
      axis.ticks.x = element_line(linewidth = 0.25),
      strip.text.x = element_text(size = 9.6),
      strip.text.y = element_text(size = 9.6),
      panel.grid.minor = element_blank(),
      panel.spacing.x = grid::unit(0.5, "lines"),
      panel.spacing.y = grid::unit(0.45, "lines"),
      legend.position = "bottom",
      legend.title = element_text(size = 9.5),
      legend.text = element_text(size = 9.5),
      plot.margin = margin(t = 3, r = 8, b = 3, l = 3)
    ) +
    labs(
      x = "Data-generation distribution",
      y = "Estimate"
    )

  save_figure(
    fig_slide_saft_treatment_effect_recovery_all40,
    "fig_slide_saft_treatment_effect_recovery_all40",
    12.8,
    5.4
  )
} else {
  warning(
    "Skipping scenario-coefficient boxplot matrix because per-simulation files are missing:\n",
    proposed_per_sim_path, "\n", jmbayes2_per_sim_path,
    call. = FALSE
  )
}

plot_config_order <- simulation_design$Config

alpha_plot_data <- all_summary %>%
  filter(variable == "alpha") %>%
  mutate(
    config_id = factor(config_id, levels = plot_config_order),
    Method = factor(Method, levels = names(method_palette))
  )

fig_alpha_bias_by_config <- alpha_plot_data %>%
  filter(!is.na(bias)) %>%
  ggplot(aes(x = config_id, y = bias, color = Method, group = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  scale_color_manual(values = method_palette, drop = FALSE) +
  facet_grid(. ~ censoring, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7)
  ) +
  labs(
    x = "Configuration",
    y = expression(alpha~"bias"),
    color = "Model",
    title = expression(alpha~"bias by ENZAMET configuration"),
    subtitle = expression("JMbayes2 "*alpha*" bias is shown only where an exact converted PH-scale truth is available.")
  )

ggsave(
  file.path(figure_dir, "fig1_alpha_bias_by_config.pdf"),
  fig_alpha_bias_by_config,
  width = 14,
  height = 6,
  units = "in"
)
ggsave(
  file.path(figure_dir, "fig1_alpha_bias_by_config.png"),
  fig_alpha_bias_by_config,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300
)

alpha_heatmap_grid <- alpha_plot_data %>%
  select(config_id, Method, bias) %>%
  complete(config_id = factor(plot_config_order, levels = plot_config_order),
           Method = factor(names(method_palette), levels = names(method_palette))) %>%
  mutate(
    config_id = factor(config_id, levels = rev(plot_config_order)),
    Method = factor(Method, levels = names(method_palette))
  )

fig_alpha_bias_heatmap <- ggplot(alpha_heatmap_grid, aes(x = Method, y = config_id, fill = bias)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "white",
    high = "#C53030",
    midpoint = 0,
    na.value = "grey85"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    x = "Model",
    y = "Configuration",
    fill = expression(alpha~"bias"),
    title = expression(alpha*"-bias heatmap across models and configurations"),
    subtitle = expression("Grey cells indicate undefined "*alpha*" truth for that fitted-model/data-generation scale.")
  )

ggsave(
  file.path(figure_dir, "fig2_alpha_bias_heatmap.pdf"),
  fig_alpha_bias_heatmap,
  width = 8.5,
  height = 10,
  units = "in"
)
ggsave(
  file.path(figure_dir, "fig2_alpha_bias_heatmap.png"),
  fig_alpha_bias_heatmap,
  width = 8.5,
  height = 10,
  units = "in",
  dpi = 300
)

longitudinal_plot_data <- all_summary %>%
  filter(variable %in% longitudinal_vars, !is.na(bias)) %>%
  mutate(
    Method = factor(Method, levels = names(method_palette)),
    PlotParameter = factor(PlotParameter, levels = parameter_plot_label(longitudinal_vars))
  )

fig_longitudinal_bias <- ggplot(longitudinal_plot_data, aes(x = Method, y = bias, fill = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  geom_boxplot(outlier.size = 0.7, width = 0.72) +
  scale_fill_manual(values = method_palette, drop = FALSE) +
  facet_wrap(~ PlotParameter, scales = "free_y", ncol = 3, labeller = label_parsed) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Bias",
    title = "Longitudinal-parameter bias across fitted models",
    subtitle = "Each box summarizes the 40 configuration-level simulation summaries."
  )

ggsave(
  file.path(figure_dir, "fig3_longitudinal_bias_comparison.pdf"),
  fig_longitudinal_bias,
  width = 11,
  height = 7,
  units = "in"
)
ggsave(
  file.path(figure_dir, "fig3_longitudinal_bias_comparison.png"),
  fig_longitudinal_bias,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300
)

coverage_plot_data <- all_summary %>%
  filter(variable %in% key_vars, !is.na(coverage)) %>%
  mutate(
    Method = factor(Method, levels = names(method_palette)),
    PlotParameter = factor(PlotParameter, levels = parameter_plot_label(key_vars))
  )

fig_coverage <- ggplot(coverage_plot_data, aes(x = Method, y = coverage, color = Method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey35") +
  geom_point(position = position_jitter(width = 0.08, height = 0), alpha = 0.65, size = 1.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, color = "black") +
  scale_color_manual(values = method_palette, drop = FALSE) +
  facet_wrap(~ PlotParameter, ncol = 4, labeller = label_parsed) +
  coord_cartesian(ylim = c(0.75, 1.0)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Coverage",
    title = "Interval coverage for key parameters",
    subtitle = "Points are configuration-level summaries; black bars are averages across configurations."
  )

ggsave(
  file.path(figure_dir, "fig4_key_parameter_coverage.pdf"),
  fig_coverage,
  width = 12,
  height = 7,
  units = "in"
)
ggsave(
  file.path(figure_dir, "fig4_key_parameter_coverage.png"),
  fig_coverage,
  width = 12,
  height = 7,
  units = "in",
  dpi = 300
)

alpha_sensitivity_plot <- alpha_prior_sensitivity_table %>%
  filter(Parameter == "alpha") %>%
  select(
    `Alpha prior SD`,
    `Mean bias`,
    `Mean RMSE`,
    `Mean coverage`,
    `Mean CI width`,
    `Mean post SD`
  ) %>%
  pivot_longer(
    cols = c(`Mean bias`, `Mean RMSE`, `Mean coverage`, `Mean CI width`, `Mean post SD`),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Mean bias", "Mean RMSE", "Mean coverage", "Mean CI width", "Mean post SD")
    )
  )

fig_alpha_prior_sensitivity <- ggplot(
  alpha_sensitivity_plot,
  aes(x = `Alpha prior SD`, y = Value)
) +
  geom_line(color = "#1B4D89", linewidth = 0.8) +
  geom_point(color = "#1B4D89", size = 2) +
  geom_vline(xintercept = alpha_prior_sd_values["maxF200"], linetype = "dotted", color = "grey35") +
  scale_x_continuous(breaks = alpha_prior_sd_values, labels = alpha_prior_sd_labels) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 11) +
  labs(
    x = expression("Prior SD for "*alpha),
    y = NULL,
    title = expression(alpha*"-prior sensitivity in log-logistic generated configurations"),
    subtitle = "The dotted line marks the baseline prior SD from the main 40-configuration simulation."
  )

ggsave(
  file.path(figure_dir, "fig5_alpha_prior_sensitivity_LL.pdf"),
  fig_alpha_prior_sensitivity,
  width = 10,
  height = 6.5,
  units = "in"
)
ggsave(
  file.path(figure_dir, "fig5_alpha_prior_sensitivity_LL.png"),
  fig_alpha_prior_sensitivity,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300
)

## Manuscript-table companion figures.
cm_plot_labels <- c(
  "CM1" = "CM1: administrative",
  "CM2" = "CM2: controlled 50%"
)

## Companion to Tables 2.1-2.3: sAFT(BP)-JM across DGP distributions.
saft_bp_distribution_bias_plot_data <- proposed_distribution_source %>%
  mutate(
    Scenario = factor(Scenario),
    PlotParameter = factor(parameter_plot_label(variable), levels = parameter_plot_label(key_vars)),
    Distribution = factor(as.character(Distribution), levels = saft_bp_distribution_labels),
    CM = factor(CM, levels = cm_labels, labels = cm_plot_labels[cm_labels])
  )

fig6_saft_bp_distribution_bias_heatmap <- ggplot(
  saft_bp_distribution_bias_plot_data,
  aes(x = Scenario, y = PlotParameter, fill = bias)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "white",
    high = "#C53030",
    midpoint = 0,
    na.value = "grey85"
  ) +
  facet_grid(CM ~ Distribution) +
  scale_y_discrete(labels = scales::parse_format()) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8)
  ) +
  labs(
    x = "Scenario",
    y = "Parameter",
    fill = "Bias",
    title = "sAFT(BP)-JM bias across data-generation distributions",
    subtitle = "Cells summarize the 1000-replicate mean bias for each parameter and scenario."
  )

save_figure(fig6_saft_bp_distribution_bias_heatmap, "fig6_saft_bp_distribution_bias_heatmap", 12, 6.8)

saft_bp_distribution_relative_bias_plot_data <- saft_bp_distribution_bias_plot_data %>%
  mutate(
    RelativeBias = ifelse(is.na(truth) | truth == 0, NA_real_, bias / abs(truth))
  )

cox_jm_distribution_relative_bias_plot_data <- all_summary %>%
  filter(as.character(Method) == "Cox-JM", variable %in% key_vars) %>%
  mutate(
    Scenario = factor(scenario),
    PlotParameter = factor(parameter_plot_label(variable), levels = parameter_plot_label(key_vars)),
    Distribution = factor(dist_short, levels = saft_bp_distribution_labels),
    CM = ifelse(censoring == "Admin", "CM1", "CM2"),
    CM = factor(CM, levels = cm_labels, labels = cm_plot_labels[cm_labels]),
    RelativeBias = ifelse(is.na(truth) | truth == 0, NA_real_, bias / abs(truth))
  )

relative_bias_limit <- safe_max(abs(c(
  saft_bp_distribution_relative_bias_plot_data$RelativeBias,
  cox_jm_distribution_relative_bias_plot_data$RelativeBias
)))

fig6b_saft_bp_distribution_relative_bias_heatmap <- ggplot(
  saft_bp_distribution_relative_bias_plot_data,
  aes(x = Scenario, y = PlotParameter, fill = RelativeBias)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "white",
    high = "#C53030",
    midpoint = 0,
    limits = c(-relative_bias_limit, relative_bias_limit),
    na.value = "grey85"
  ) +
  facet_grid(CM ~ Distribution) +
  scale_y_discrete(labels = scales::parse_format()) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8)
  ) +
  labs(
    x = "Scenario",
    y = "Parameter",
    fill = "Relative bias",
    title = "sAFT(BP)-JM relative bias across data-generation distributions",
    subtitle = "Relative bias is bias divided by abs(truth); grey cells indicate zero or unavailable truth."
  )

save_figure(
  fig6b_saft_bp_distribution_relative_bias_heatmap,
  "fig6b_saft_bp_distribution_relative_bias_heatmap",
  12,
  6.8
)

fig6c_cox_jm_distribution_relative_bias_heatmap <- ggplot(
  cox_jm_distribution_relative_bias_plot_data,
  aes(x = Scenario, y = PlotParameter, fill = RelativeBias)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "white",
    high = "#C53030",
    midpoint = 0,
    limits = c(-relative_bias_limit, relative_bias_limit),
    na.value = "grey85"
  ) +
  facet_grid(CM ~ Distribution) +
  scale_y_discrete(labels = scales::parse_format()) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8)
  ) +
  labs(
    x = "Scenario",
    y = "Parameter",
    fill = "Relative bias",
    title = "Cox-JM relative bias across data-generation distributions",
    subtitle = "Relative bias is bias divided by abs(truth); grey cells indicate zero or unavailable Cox-scale truth."
  )

save_figure(
  fig6c_cox_jm_distribution_relative_bias_heatmap,
  "fig6c_cox_jm_distribution_relative_bias_heatmap",
  12,
  6.8
)

fig7_saft_bp_distribution_esd_psd <- ggplot(
  saft_bp_distribution_bias_plot_data,
  aes(x = empirical_sd, y = mean_post_sd, color = Distribution, shape = CM)
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey35") +
  geom_point(size = 2, alpha = 0.85) +
  facet_wrap(~ PlotParameter, scales = "free", ncol = 4, labeller = label_parsed) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") +
  labs(
    x = "Empirical SD",
    y = "Posterior SD",
    color = "Distribution",
    shape = "Censoring",
    title = "sAFT(BP)-JM uncertainty calibration across distributions",
    subtitle = "Points near the dashed line indicate agreement between empirical and posterior SD."
  )

save_figure(fig7_saft_bp_distribution_esd_psd, "fig7_saft_bp_distribution_esd_psd", 11, 7.2)

## Companion to Tables 3.1-3.3: longitudinal comparison under LL DGP.
longitudinal_ll_model_plot_data <- longitudinal_ll_source %>%
  mutate(
    Model = factor(as.character(Model), levels = longitudinal_ll_model_labels),
    CM = factor(CM, levels = cm_labels, labels = cm_plot_labels[cm_labels]),
    PlotParameter = factor(parameter_plot_label(variable), levels = parameter_plot_label(longitudinal_vars))
  ) %>%
  group_by(Model, CM, PlotParameter) %>%
  summarise(
    Bias = safe_mean(bias),
    RMSE = safe_mean(rmse),
    Coverage = safe_mean(coverage),
    .groups = "drop"
  )

fig8_longitudinal_LL_bias_rmse <- longitudinal_ll_model_plot_data %>%
  pivot_longer(cols = c(Bias, RMSE), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Bias", "RMSE"))) %>%
  ggplot(aes(x = Model, y = Value, color = CM, group = CM)) +
  geom_hline(data = data.frame(Metric = "Bias"), aes(yintercept = 0),
             inherit.aes = FALSE, linetype = "dashed", color = "grey35") +
  geom_point(position = position_dodge(width = 0.45), size = 2) +
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.6) +
  facet_grid(Metric ~ PlotParameter, scales = "free_y", labeller = labeller(PlotParameter = label_parsed)) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Censoring",
    title = "Longitudinal bias and RMSE under log-logistic data generation",
    subtitle = "Values are averaged across the five log-logistic scenarios."
  )

save_figure(fig8_longitudinal_LL_bias_rmse, "fig8_longitudinal_LL_bias_rmse", 13, 6.8)

fig_slide_longitudinal_LL_beta_bias_rmse <- longitudinal_ll_model_plot_data %>%
  filter(PlotParameter %in% parameter_plot_label(c("beta_long[1]", "beta_long[2]", "beta_long[3]"))) %>%
  pivot_longer(cols = c(Bias, RMSE), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = c("Bias", "RMSE")),
    PlotParameter = factor(
      as.character(PlotParameter),
      levels = parameter_plot_label(c("beta_long[1]", "beta_long[2]", "beta_long[3]"))
    ),
    Model = factor(
      as.character(Model),
      levels = longitudinal_ll_model_labels,
      labels = c("LMM", "sAFT-JM", "Cox-JM")
    ),
    CM = factor(as.character(CM), levels = cm_plot_labels[cm_labels], labels = cm_labels)
  ) %>%
  ggplot(aes(x = Model, y = Value, color = CM, group = CM)) +
  geom_hline(
    data = data.frame(Metric = "Bias"),
    aes(yintercept = 0),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.35
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.8) +
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.75) +
  facet_grid(Metric ~ PlotParameter, scales = "free_y", labeller = labeller(PlotParameter = label_parsed)) +
  scale_color_manual(
    name = "Censoring",
    values = c("CM1" = "#A6CEE3", "CM2" = "#1F78B4"),
    drop = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(1.0, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 10.5),
    strip.text.x = element_text(size = 10.5),
    strip.text.y = element_text(size = 10.5),
    plot.margin = margin(t = 3, r = 8, b = 3, l = 3)
  ) +
  labs(
    x = NULL,
    y = NULL
  )

save_figure(
  fig_slide_longitudinal_LL_beta_bias_rmse,
  "fig_slide_longitudinal_LL_beta_bias_rmse",
  11.6,
  5.2
)

fig9_longitudinal_LL_coverage <- ggplot(
  longitudinal_ll_model_plot_data %>% filter(!is.na(Coverage)),
  aes(x = PlotParameter, y = Coverage, color = Model, shape = CM)
) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey35") +
  geom_point(position = position_dodge(width = 0.45), size = 2.2) +
  scale_x_discrete(labels = scales::parse_format()) +
  coord_cartesian(ylim = c(0.90, 0.98)) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    x = "Parameter",
    y = "Coverage",
    color = "Model",
    shape = "Censoring",
    title = "Longitudinal interval coverage under log-logistic data generation",
    subtitle = "Coverage for LMM variance parameters is unavailable and omitted."
  )

save_figure(fig9_longitudinal_LL_coverage, "fig9_longitudinal_LL_coverage", 9.5, 5.5)

## Companion to Tables 5.1-5.3: alpha-prior sensitivity.
alpha_prior_gamma_alpha_plot_data <- alpha_sensitivity %>%
  filter(variable %in% c("beta_surv[1]", "alpha")) %>%
  mutate(
    ScenarioLabel = factor(paste0("S", scenario), levels = paste0("S", 1:5)),
    Parameter = factor(parameter_plot_label(variable), levels = parameter_plot_label(c("alpha", "beta_surv[1]"))),
    `Alpha prior SD` = alpha_prior_sd_values[paste0("maxF", `Alpha prior maxF`)],
    CM = ifelse(censoring == "Admin", "CM1", "CM2"),
    CM = factor(CM, levels = cm_labels, labels = cm_plot_labels[cm_labels])
  ) %>%
  group_by(`Alpha prior SD`, ScenarioLabel, Parameter, CM, truth) %>%
  summarise(Estimate = safe_mean(mean), .groups = "drop")

fig10_alpha_prior_gamma_alpha_estimates <- ggplot(
  alpha_prior_gamma_alpha_plot_data,
  aes(x = `Alpha prior SD`, y = Estimate, color = CM, group = CM)
) +
  geom_hline(aes(yintercept = truth), linetype = "dashed", color = "black") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = alpha_prior_sd_values,
    labels = alpha_prior_sd_labels,
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  facet_grid(Parameter ~ ScenarioLabel, scales = "free_y", labeller = labeller(Parameter = label_parsed)) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.spacing.x = grid::unit(0.8, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
  ) +
  labs(
    x = expression("Prior SD for "*alpha),
    y = "Mean estimate",
    color = "Censoring",
    title = expression(alpha~"and"~gamma~"estimates across"~alpha*"-prior settings"),
    subtitle = "Dashed lines show the true value in each scenario."
  )

save_figure(fig10_alpha_prior_gamma_alpha_estimates, "fig10_alpha_prior_gamma_alpha_estimates", 13, 6.5)

alpha_prior_bias_coverage_plot_data <- alpha_sensitivity %>%
  filter(variable %in% alpha_prior_vars) %>%
  mutate(
    PlotParameter = factor(parameter_plot_label(variable), levels = parameter_plot_label(alpha_prior_vars)),
    `Alpha prior SD` = alpha_prior_sd_values[paste0("maxF", `Alpha prior maxF`)],
    `Alpha prior SD label` = factor(
      format_number(`Alpha prior SD`, digits = 3),
      levels = alpha_prior_sd_labels
    ),
    CM = ifelse(censoring == "Admin", "CM1", "CM2"),
    CM = factor(CM, levels = cm_labels, labels = cm_plot_labels[cm_labels])
  ) %>%
  group_by(`Alpha prior SD label`, PlotParameter, CM) %>%
  summarise(
    Bias = safe_mean(bias),
    `Coverage - 0.95` = safe_mean(coverage) - 0.95,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Bias, `Coverage - 0.95`), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Bias", "Coverage - 0.95")))

fig11_alpha_prior_bias_coverage_heatmap <- ggplot(
  alpha_prior_bias_coverage_plot_data,
  aes(x = `Alpha prior SD label`, y = PlotParameter, fill = Value)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "white",
    high = "#C53030",
    midpoint = 0,
    na.value = "grey85"
  ) +
  facet_grid(Metric ~ CM, scales = "free") +
  scale_y_discrete(labels = scales::parse_format()) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 9)) +
  labs(
    x = expression("Prior SD for "*alpha),
    y = "Parameter",
    fill = "Value",
    title = expression(alpha*"-prior sensitivity for bias and coverage"),
    subtitle = "Coverage values are shown as coverage minus 0.95."
  )

save_figure(fig11_alpha_prior_bias_coverage_heatmap, "fig11_alpha_prior_bias_coverage_heatmap", 10.5, 6.8)

message("\nDone.")
message("Tables:  ", table_dir)
message("Figures: ", figure_dir)
message("CSVs:    ", csv_dir)
