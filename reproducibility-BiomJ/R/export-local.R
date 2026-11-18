# Editor-facing reduced-run records, separate from all manuscript/full-run data.
audit_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(file), fileext = ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, file)) {
    file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
  }
}

audit_record_one <- function(root, row, out, result) {
  meta_names <- c("collection", "full_config_id", "replicate_id", "unique_seed", "baseline_label", "scenario", "cens_label", "n_patients", "max_FU", "alpha_max_f")
  meta <- row[, meta_names, drop = FALSE]
  annotate <- function(x) {
    x <- as.data.frame(x)
    if (!nrow(x)) return(x)
    x <- x[, !names(x) %in% names(meta), drop = FALSE]
    cbind(meta[rep(1, nrow(x)), , drop = FALSE], x, row.names = NULL)
  }
  if (row$collection == "comparator") {
    x <- as.data.frame(result$comparison_tidy)
    x$method <- x$model
    # Use neutral terminology: LMM SD entries are frequentist standard errors,
    # whereas JM entries summarize posterior uncertainty.
    x$uncertainty_sd <- x$sd
    x$interval_type <- ifelse(x$method == "LMM", "confidence", "credible")
    runtime <- as.data.frame(result$runtimes)
    status <- as.data.frame(result$model_status)
    status$method <- status$model
    audit_write_csv(annotate(status), file.path(out, "local-model-status.csv"))
  } else {
    x <- read.csv(file.path(out, "posterior-summary-original-scale.csv"), check.names = FALSE)
    x$method <- "sAFT-JM"
    x$uncertainty_sd <- x$sd
    x$interval_type <- "credible"
    runtime <- as.data.frame(result$runtimes)
    d <- readRDS(file.path(out, "fit-audit.rds"))$diagnostics
    chain <- data.frame(chain_id = seq_along(d$num_divergent),
      num_divergent = d$num_divergent, num_max_treedepth = d$num_max_treedepth,
      ebfmi = d$ebfmi)
    audit_write_csv(annotate(chain), file.path(out, "local-chain-diagnostics.csv"))
    ic <- cbind(as.data.frame(result$loo_joint), as.data.frame(result$waic_joint))
    if (nrow(ic)) audit_write_csv(annotate(ic), file.path(out, "local-information-criteria.csv"))
  }
  keep <- intersect(c("variable", "method", "mean", "median", "uncertainty_sd", "q2.5", "q97.5", "rhat", "ess_bulk", "ess_tail", "ess", "mcse_mean", "mcse_q2.5", "mcse_q97.5", "interval_type"), names(x))
  audit_write_csv(annotate(x[, keep, drop = FALSE]), file.path(out, "local-parameter-records.csv"))
  audit_write_csv(annotate(runtime), file.path(out, "local-runtime-records.csv"))
  comparison <- file.path(out, "comparison-with-production.csv")
  if (file.exists(comparison)) audit_write_csv(annotate(read.csv(comparison, check.names = FALSE)), file.path(out, "local-reference-comparison.csv"))
}

audit_export_local <- function(root, quiet = TRUE) {
  audit_require(c("dplyr", "tidyr"))
  base <- file.path(root, "output", "local-simulations")
  plan_file <- file.path(base, "run-plan.csv")
  if (!file.exists(plan_file)) {
    message("No local fitting plan saved yet. Use --mode=plan or start a selected run.")
    return(invisible(NULL))
  }
  plan <- read.csv(plan_file, stringsAsFactors = FALSE)
  tables <- list(parameters = list(), runtimes = list(), diagnostics = list(), model_status = list(), comparisons = list(), information_criteria = list())
  for (j in seq_len(nrow(plan))) {
    out <- file.path(base, "fits", plan$output_key[j])
    marker <- file.path(out, "completed.rds")
    failure <- file.path(out, "failure.rds")
    plan$status[j] <- if (file.exists(marker)) "completed" else if (file.exists(failure)) "failed" else if (file.exists(file.path(out, "configuration.rds"))) "started_or_interrupted" else "not_started"
    plan$error[j] <- if (plan$status[j] == "failed") readRDS(failure)$message else NA_character_
    if (plan$status[j] != "completed") next
    filenames <- c(parameters = "local-parameter-records.csv", runtimes = "local-runtime-records.csv", diagnostics = "local-chain-diagnostics.csv", model_status = "local-model-status.csv", comparisons = "local-reference-comparison.csv", information_criteria = "local-information-criteria.csv")
    for (name in names(filenames)) {
      f <- file.path(out, filenames[[name]])
      if (file.exists(f)) {
        x <- read.csv(f, check.names = FALSE)
        for (col in intersect(c("collection", "cens_label", "baseline_label", "method", "variable", "model", "status", "error_message", "retry_error_message", "interval_type"), names(x))) x[[col]] <- as.character(x[[col]])
        tables[[name]][[length(tables[[name]]) + 1L]] <- x
      }
    }
  }
  out <- file.path(base, "editor-records")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  merged <- lapply(tables, dplyr::bind_rows)
  export_names <- c(parameters = "replicate_parameter_records.csv", runtimes = "replicate_runtime_records.csv", diagnostics = "stan_chain_diagnostics.csv", model_status = "comparator_model_status.csv", comparisons = "comparison_with_historical_production.csv", information_criteria = "saft_information_criteria.csv")
  for (name in names(merged)) audit_write_csv(merged[[name]], file.path(out, export_names[[name]]))
  audit_write_csv(plan, file.path(out, "progress.csv"))
  counts <- as.data.frame(table(plan$collection, plan$status), responseName = "dataset_runs")
  names(counts)[1:2] <- c("collection", "status")
  audit_write_csv(counts, file.path(out, "progress_counts.csv"))
  p <- merged$parameters
  summary <- data.frame()
  if (nrow(p)) {
    summary <- p |>
      dplyr::group_by(collection, full_config_id, baseline_label, scenario, cens_label, method, variable) |>
      dplyr::summarise(n_available = sum(is.finite(mean)),
        mean_estimate = if (any(is.finite(mean))) base::mean(mean[is.finite(mean)]) else NA_real_,
        median_estimate = if (any(is.finite(mean))) stats::median(mean[is.finite(mean)]) else NA_real_,
        sd_across_estimates = stats::sd(mean[is.finite(mean)]), .groups = "drop")
  }
  audit_write_csv(summary, file.path(out, "reduced_run_summary.csv"))
  saveRDS(c(list(plan = plan, reduced_summary = summary), merged), file.path(out, "local-audit-records.rds"))
  if (!quiet) print(counts[counts$dataset_runs > 0, ], row.names = FALSE)
  invisible(counts)
}
