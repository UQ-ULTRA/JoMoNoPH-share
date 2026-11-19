#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)

# Disposable fixtures exercise incomplete installations without MCMC or real results.
local({
  parent <- tempfile("biomj-overall-check-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  fixture <- function(name) {
    path <- file.path(parent, name)
    dir.create(file.path(path, "case-study", "source"), recursive = TRUE)
    dir.create(file.path(path, "production", "Stan"), recursive = TRUE)
    saveRDS(list(fixture_only = TRUE), file.path(path, "case-study/source/mock-data.rds"))
    file.create(file.path(path, "case-study/source/fit-case-study-model.R"))
    file.create(file.path(path,
      "production/Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"))
    e <- new.env(parent = globalenv())
    sys.source(file.path(root, "R/audit.R"), envir = e)
    sys.source(file.path(root, "R/session.R"), envir = e)
    sys.source(file.path(root, "R/entry.R"), envir = e)
    e$requireNamespace <- function(...) FALSE
    e$audit_grid <- function(root) data.frame(full_config_id = 1:40)
    e$audit_environment <- function(out) {
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      writeLines("Fixture environment only", file.path(out, "environment.txt"))
    }
    list(root = path, env = e)
  }
  run <- function(f, ...) {
    notes <- character()
    result <- NULL
    invisible(capture.output(suppressMessages(
      result <- withCallingHandlers(f$env$audit_run(f$root, task = "check", ...),
        warning = function(w) {
          notes <<- c(notes, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    )))
    list(report = result, warnings = notes)
  }
  status <- function(x, item) x$report$status[match(item, x$report$item)]

  f <- fixture("missing-archive-and-packages")
  x <- run(f, cmdstan = "/missing/requested/cmdstan")
  stopifnot(length(x$warnings) == 1L,
    grepl("Overall check completed with unresolved items", x$warnings),
    status(x, "Archive batches main") == "WARNING",
    status(x, "R package cmdstanr") == "WARNING",
    status(x, "CmdStan installation") == "WARNING",
    status(x, "Basic C++ toolchain") == "WARNING",
    status(x, "Mock data") == "PASS",
    status(x, "Configuration grid") == "PASS",
    status(x, "Environment record") == "PASS")
  saved <- read.csv(file.path(f$root, "output/check/check-summary.csv"))
  stopifnot(identical(saved, x$report),
    file.exists(file.path(f$root, "output/check/configurations.csv")),
    file.exists(file.path(f$root, "output/check/packages.csv")),
    file.exists(file.path(f$root, "output/check/archive-inventory.csv")))

  f <- fixture("missing-inputs-grid-and-environment")
  unlink(file.path(f$root, "case-study/source/mock-data.rds"))
  f$env$audit_grid <- function(root) stop("Fixture missing grid")
  f$env$audit_environment <- function(out) stop("Fixture environment unavailable")
  x <- run(f)
  stopifnot(length(x$warnings) == 1L,
    status(x, "Mock data") == "WARNING",
    status(x, "Configuration grid") == "WARNING",
    status(x, "Environment record") == "WARNING",
    status(x, "Mock fitting script") == "PASS",
    status(x, "Stan source") == "PASS",
    status(x, "Package inventory export") == "PASS")

  f <- fixture("environment-warning")
  f$env$audit_environment <- function(out) warning("Fixture metadata warning")
  x <- run(f)
  stopifnot(length(x$warnings) == 1L,
    status(x, "Environment record") == "WARNING",
    grepl("Fixture metadata warning", x$warnings))

  f <- fixture("namespace-error")
  attempted <- character()
  f$env$requireNamespace <- function(package, ...) {
    attempted <<- c(attempted, package)
    if (package == "dplyr") stop("Fixture namespace error")
    FALSE
  }
  x <- run(f)
  stopifnot(length(x$warnings) == 1L, "coda" %in% attempted,
    status(x, "R package dplyr") == "WARNING",
    status(x, "Environment record") == "PASS")

  f <- fixture("blocked-output")
  dir.create(file.path(f$root, "output"))
  file.create(file.path(f$root, "output/check"))
  x <- run(f)
  stopifnot(length(x$warnings) == 1L,
    status(x, "Output directory") == "WARNING",
    status(x, "Configuration grid") == "PASS",
    status(x, "Check summary export") == "WARNING",
    grepl("summary could not be saved", x$warnings))

  # Fake batch listings test inventory counts without creating 21,500 result files.
  f <- fixture("complete-inventory")
  cols <- f$env$audit_collections(f$root)
  counts <- c(10000L, 2500L, 2500L, 2500L, 4000L)
  for (p in c(cols$input, cols$saved_summary)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  f$env$list.files <- function(path, ...) {
    i <- match(path, cols$input)
    if (!is.na(i)) rep("fixture.rds", counts[i]) else base::list.files(path, ...)
  }
  x <- run(f)
  stopifnot(all(x$report$status[grepl("^Archive batches|^Saved summaries", x$report$item)]
    == "PASS"), length(x$warnings) == 1L)

  # Exercise the no-warning branch when the local software/toolchain are available.
  f$env$requireNamespace <- base::requireNamespace
  x <- run(f)
  if (all(x$report$status == "PASS")) {
    stopifnot(!length(x$warnings))
    before <- cmdstanr::cmdstan_path()
    y <- run(f, cmdstan = file.path(parent, "nonexistent-cmdstan"))
    stopifnot(length(y$warnings) == 1L,
      status(y, "CmdStan installation") == "WARNING",
      status(y, "Basic C++ toolchain") == "PASS",
      status(y, "Environment record") == "PASS",
      identical(cmdstanr::cmdstan_path(), before))
    cat("PASS: complete preflight returns without warnings.\n")
  } else {
    stopifnot(length(x$warnings) == 1L)
    cat("CHECKED: local software gaps reported without stopping the inventory.\n")
  }
})
cat("PASS: overall check continues after missing components and reports one final warning. No MCMC.\n")
