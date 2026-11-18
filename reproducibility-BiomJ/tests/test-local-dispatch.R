#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "audit.R"))
source(file.path(root, "R", "simulation.R"))
source(file.path(root, "R", "export-local.R"))
# Regression: changing a build-only flag must not invalidate prepared inputs,
# but must invalidate an actual Stan fit. Relevant identity differences remain
# strict in both modes, and the recorded identity must not be modified.
saved_identity <- list(seed = 1000001, source_md5 = "original", no_pch = "true")
current_identity <- saved_identity
current_identity$no_pch <- "false"
stopifnot(!length(audit_checkpoint_differences(saved_identity, current_identity, preparation_only = TRUE)),
  !length(audit_checkpoint_differences(saved_identity, current_identity, comparator = TRUE)),
  identical(audit_checkpoint_differences(saved_identity, current_identity), "no_pch"),
  identical(saved_identity$no_pch, "true"))
current_identity$source_md5 <- "changed"
stopifnot(identical(audit_checkpoint_differences(saved_identity, current_identity, preparation_only = TRUE), "source_md5"))
current_identity$seed <- 1000002
stopifnot(setequal(audit_checkpoint_differences(saved_identity, current_identity, preparation_only = TRUE), c("seed", "source_md5")))
# Verified source cleanup must not accept a different model or configuration grid.
saved_identity <- list(seed = 1000001, source_md5 = c("old-driver", "old-generator", "stan"))
current_identity <- list(seed = 1000001, source_md5 = c("driver", "generator", "stan", "grid"))
aliases <- data.frame(before_md5 = c("old-driver", "old-generator"),
  after_md5 = c("driver", "generator"), extra_md5 = c("grid", ""))
stopifnot(!length(audit_checkpoint_differences(saved_identity, current_identity, source_aliases = aliases)),
  identical(saved_identity$source_md5, c("old-driver", "old-generator", "stan")))
for (field in 3:4) {
  changed <- current_identity
  changed$source_md5[field] <- "unverified"
  stopifnot(identical(audit_checkpoint_differences(saved_identity, changed, source_aliases = aliases), "source_md5"))
}
for (mode in c("prepare", "spot-check", "sensitivity", "comparator")) {
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(root, "main.R")), paste0("--mode=", mode)), stdout = TRUE, stderr = TRUE))
  stopifnot(!is.null(attr(output, "status")), attr(output, "status") != 0,
    any(grepl("Select settings explicitly", output, fixed = TRUE)))
}
plan <- audit_simulation_plan(root, "all")
stopifnot(nrow(plan) == 550L, !anyDuplicated(plan$output_key))
stopifnot(all(startsWith(plan$reference_file, "archived-production-results/")),
  all(file.exists(file.path(root, plan$reference_file))))
stopifnot(all(table(plan$collection)[c("main", "alpha125", "alpha150", "alpha300", "comparator")] == c(200, 50, 50, 50, 200)))
stopifnot(all(table(plan$collection, plan$full_config_id)["main", ] == 5))
stopifnot(all(plan$unique_seed == 1e6 * plan$full_config_id + plan$replicate_id))
stopifnot(all(plan$full_config_id[grepl("^alpha", plan$collection)] %in% 1:10))
stopifnot(nrow(audit_simulation_plan(root, "all", replicates = 1:10)) == 1100L)
bad <- try(audit_simulation_plan(root, "alpha125", config_ids = 11), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
for (id in unique(plan$collection)) {
  p <- plan[plan$collection == id & plan$full_config_id == 1 & plan$replicate_id == 5, , drop = FALSE]
  ref <- audit_plan_reference(root, p)
  stopifnot(ref$unique_seed == 1000005)
  stopifnot(p$reference_batch == if (id == "comparator") 1L else 2L)
  stopifnot(p$reference_position == if (id == "comparator") 5L else 1L)
}
# Synthetic export fixtures test mixed numeric/admin censor labels, duplicate
# metadata names and absent uncertainty columns. They are not model results.
fixture <- tempfile("biomj-local-export-test-")
dir.create(fixture, recursive = TRUE, showWarnings = FALSE)
p <- audit_simulation_plan(root, "main", config_ids = c(1, 2), replicates = 1)
base <- file.path(fixture, "output", "local-simulations")
audit_write_csv(p, file.path(base, "run-plan.csv"))
for (j in 1:2) {
  out <- file.path(base, "fits", p$output_key[j])
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  audit_write_csv(data.frame(variable = "alpha", mean = 0.01, sd = 0.003,
    rhat = 1.001, mcse_mean = 0.0001), file.path(out, "posterior-summary-original-scale.csv"))
  saveRDS(list(diagnostics = list(num_divergent = rep(0L, 4), num_max_treedepth = rep(0L, 4), ebfmi = rep(.8, 4))), file.path(out, "fit-audit.rds"))
  fake <- list(runtimes = data.frame(unique_seed = p$unique_seed[j], total_time = 60))
  audit_record_one(root, p[j, , drop = FALSE], out, fake)
  saveRDS(list(test_fixture_only = TRUE), file.path(out, "completed.rds"))
}
audit_export_local(fixture)
x <- read.csv(file.path(base, "editor-records", "replicate_parameter_records.csv"))
stopifnot(nrow(x) == 2L, !anyDuplicated(names(x)), setequal(x$cens_label, c("50", "admin")))
s <- read.csv(file.path(base, "editor-records", "progress.csv"))
stopifnot(all(s$status == "completed"))
unlink(fixture, recursive = TRUE)
cat("PASS: collection counts, seeds, reference slots, selection and CSV exports. No MCMC run.\n")
