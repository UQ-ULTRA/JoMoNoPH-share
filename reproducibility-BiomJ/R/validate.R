# Re-runnable, non-sampling package checks. No access to the original repo needed.
audit_validate <- function(root) {
  rows <- list()
  add <- function(item, status, detail) {
    rows[[length(rows) + 1]] <<- data.frame(item = item, status = status, detail = detail)
  }
  scripts <- list.files(root, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  invisible(lapply(scripts, parse))
  add("R syntax", "PASS", paste(length(scripts), "R scripts parsed"))
  stan <- file.path(root, "production", "Stan", "Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan")
  provenance <- read.csv(file.path(root, "documentation", "source-provenance.csv"))
  stan_row <- provenance$role == "unchanged Stan model"
  stopifnot(sum(stan_row) == 1L,
    unname(tools::md5sum(stan)) == provenance$before_md5[stan_row])
  add("Stan source consistency", "PASS", "One shared Stan source; original model checksum unchanged")
  grid <- audit_grid(root)
  stopifnot(nrow(grid) == 40, !anyDuplicated(grid$full_config_id))
  add("Configuration grid", "PASS", "40 unique original configuration IDs")
  cols <- audit_collections(root)
  expected <- c(10000, 2500, 2500, 2500, 4000)
  counts <- vapply(cols$input, function(p) length(list.files(p, pattern = "\\.rds$")), integer(1))
  add("Merged archives", if (all(counts == expected)) "PASS" else "FAIL", paste(counts, collapse = ", "))
  plan <- audit_simulation_plan(root, "all", replicates = 1:5)
  planned <- as.integer(table(factor(plan$collection, levels = cols$collection)))
  stopifnot(nrow(plan) == 550L, identical(planned, c(200L, 50L, 50L, 50L, 200L)), !anyDuplicated(plan$output_key))
  add("Available selection space", "PASS", "All configuration/replicate mappings verified; running all 550 candidate combinations is not required")
  selected_file <- file.path(root, "output", "local-simulations", "run-plan.csv")
  selected <- if (file.exists(selected_file)) read.csv(selected_file) else plan[FALSE, ]
  for (id in unique(selected$collection)) {
    keys <- selected$output_key[selected$collection == id]
    n <- sum(file.exists(file.path(root, "output", "local-simulations", "fits", keys, "completed.rds")))
    add(paste("Recorded local checks", id), if (n == length(keys)) "EXECUTED" else "PENDING", paste(n, "of", length(keys), "recorded dataset checkpoints completed; inspect diagnostics separately"))
  }
  # Inspect current local fits only; validation never changes saved fit records.
  for (j in seq_len(nrow(selected))) {
    p <- selected[j, , drop = FALSE]
    folder <- file.path(root, "output", "local-simulations", "fits", p$output_key)
    fit_file <- file.path(folder, "fit-audit.rds")
    if (!file.exists(file.path(folder, "completed.rds")) || !file.exists(fit_file)) next
    x <- readRDS(fit_file)
    input <- readRDS(file.path(folder, "inputs-and-initial-values.rds"))
    draws <- readRDS(file.path(folder, "coefficient-draws-centered-scale.rds"))
    recalculated <- audit_parameter_summary(draws, input$longitudinal_response_center)
    same <- isTRUE(all.equal(as.data.frame(x$posterior_summary),
      as.data.frame(recalculated), check.attributes = FALSE, tolerance = 1e-10))
    add(paste("Saved summary", p$output_key), if (same) "PASS" else "FAIL",
      "Compared saved coefficient summaries with saved posterior draws")
    sm <- as.data.frame(recalculated)
    div <- sum(x$diagnostics$num_divergent)
    depth <- sum(x$diagnostics$num_max_treedepth)
    rh <- max(sm$rhat, na.rm = TRUE)
    add(paste("Sampling diagnostics", p$output_key),
      if (div > 0 || depth > 0 || rh >= 1.01) "REVIEW" else "CHECKED",
      sprintf("%d divergences; %d treedepth hits; max coefficient R-hat %.4f; %.1f sampling min",
        div, depth, rh, x$sampling_seconds / 60))
  }
  refdir <- file.path(root, "archived-production-results", "auxiliary", "manuscript-outputs")
  for (type in c("derived_csv", "tables")) {
    files <- list.files(file.path(refdir, type), pattern = if (type == "derived_csv") "\\.csv$" else "\\.tex$")
    exists <- file.exists(file.path(root, "output", "manuscript", type, files))
    equal <- unname(tools::md5sum(file.path(refdir, type, files))) ==
      unname(tools::md5sum(file.path(root, "output", "manuscript", type, files)))
    status <- if (all(exists) && all(equal)) "PASS" else "FAIL"
    add(paste("Manuscript", type), status, paste(sum(equal, na.rm = TRUE), "of", length(files), "reference files byte-identical"))
  }
  for (i in seq_len(nrow(cols))) {
    f <- if (cols$collection[i] == "comparator") "jmbayes2_comparison_summary_loglogistic_weibull.csv" else if (cols$collection[i] == "main") "jm_parameter_summary_long_loglogistic_weibull.csv" else "jm_parameter_summary_long_loglogistic.csv"
    fresh <- file.path(cols$regenerated_summary[i], f)
    if (file.exists(fresh)) {
      a <- read.csv(file.path(cols$saved_summary[i], f), check.names = FALSE)
      b <- read.csv(fresh, check.names = FALSE)
      # Sorting can differ across dplyr versions. Compare after stable key sort.
      keys <- intersect(c("config", "model", "method", "variable"), names(a))
      if (length(keys)) {
        a <- a[do.call(order, a[keys]), , drop = FALSE]
        b <- b[do.call(order, b[keys]), , drop = FALSE]
        rownames(a) <- rownames(b) <- NULL
      }
      same <- isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))
      add(paste("Reaggregate", cols$collection[i]), if (same) "PASS" else "REVIEW", if (same) "Full regenerated parameter table agrees within 1e-10" else paste(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE), collapse = "; "))
    } else add(paste("Reaggregate", cols$collection[i]), "PENDING", "Full saved-batch aggregation not yet completed locally")
  }
  mock_fit <- file.path(root, "output", "case-study", "fit-audit.rds")
  if (file.exists(mock_fit)) {
    folder <- dirname(mock_fit)
    x <- readRDS(mock_fit)
    input <- readRDS(file.path(folder, "inputs-and-initial-values.rds"))
    draws <- readRDS(file.path(folder, "coefficient-draws-centered-scale.rds"))
    recalculated <- audit_parameter_summary(draws, input$longitudinal_response_center)
    same <- isTRUE(all.equal(as.data.frame(x$posterior_summary),
      as.data.frame(recalculated), check.attributes = FALSE, tolerance = 1e-10))
    add("Mock saved summary", if (same) "PASS" else "FAIL",
      "Compared saved coefficient summaries with saved posterior draws")
  }
  add("Mock case-study fit", if (file.exists(mock_fit)) "EXECUTED" else "PENDING", "Inspect output/case-study diagnostics and original supplied mock reference table")
  out <- file.path(root, "documentation", "validation")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  report <- do.call(rbind, rows)
  write.csv(report, file.path(out, "checks.csv"), row.names = FALSE)
  lines <- c("# Validation Report", "", paste("Recorded:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
    "These checks distinguish completed execution from acceptable sampling diagnostics.",
    "EXECUTED is not an automatic declaration of convergence or statistical equivalence.", "",
    paste0("- **", report$status, "** ", report$item, ": ", report$detail), "",
    "## Still Requires Author Sign-Off", "",
    "- Confirm final manuscript/supplement numbering and coverage of every reported result.",
    "- Have a coauthor test the final package on another machine.",
    "- Confirm data/code redistribution permissions and the large-archive delivery URL.",
    "- Review all new-fit diagnostics before claiming successful numerical reproduction.")
  writeLines(lines, file.path(root, "documentation", "validation-report.md"))
  audit_environment(out)
  print(report, row.names = FALSE)
  if (any(report$status == "FAIL")) stop("One or more package validation checks failed.")
}
