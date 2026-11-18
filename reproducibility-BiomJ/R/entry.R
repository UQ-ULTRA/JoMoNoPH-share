# Shared task dispatch for the RStudio script and optional command-line entry.
# Statistical calculations remain in the existing audit/production functions.
audit_run <- function(root, task = "plan", collection = "main", config_ids = NULL,
                      replicates = 1:5, parallel_chains = 4L, cmdstan = "",
                      no_pch = FALSE, resume = TRUE, summaries = "saved") {
  restore_session <- audit_session()
  on.exit(restore_session(), add = TRUE)
  aliases <- c(fit = "spot-check", records = "export-local")
  if (!is.character(task) || length(task) != 1L || is.na(task)) stop("Choose one task.")
  if (task %in% names(aliases)) task <- unname(aliases[task])
  tasks <- c("check", "results", "aggregate", "plan", "prepare", "spot-check",
    "sensitivity", "comparator", "export-local", "local-status",
    "case-study", "manifest", "validate")
  if (!task %in% tasks) stop("Unknown task: ", task)
  if (!is.character(collection) || length(collection) != 1L || is.na(collection)) stop("Choose one collection (or alpha/all).")
  if (!collection %in% c("main", "alpha125", "alpha150", "alpha300", "alpha", "comparator", "all")) stop("Unknown collection: ", collection)
  ids <- function(x, name, upper) {
    if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < 1 | x > upper) || anyDuplicated(x)) {
      stop(name, " must contain distinct integer IDs from 1 to ", upper, ".")
    }
    as.integer(x)
  }
  explicit_all <- identical(config_ids, "all")
  if (explicit_all) config_ids <- NULL
  if (!is.null(config_ids)) config_ids <- ids(config_ids, "config_ids", 40L)
  replicates <- ids(replicates, "replicates", 1000L)
  if (task %in% c("prepare", "spot-check", "sensitivity", "comparator") && is.null(config_ids) && !explicit_all) {
    stop("Select config_ids explicitly before preparing or fitting. Use task = 'plan' to list all candidates.")
  }
  parallel_chains <- ids(parallel_chains, "parallel_chains", 4L)
  if (length(parallel_chains) != 1L) stop("parallel_chains must be one integer.")
  for (value in list(no_pch, resume)) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("no_pch and resume must be TRUE or FALSE.")
  }
  if (!is.character(cmdstan) || length(cmdstan) != 1L || is.na(cmdstan)) stop("cmdstan must be a path or an empty string.")
  if (!is.character(summaries) || length(summaries) != 1L || is.na(summaries) ||
    !summaries %in% c("saved", "regenerated")) stop("summaries must be saved or regenerated.")

  # The session guard restores build and prior options after every task.
  Sys.setenv(BIOMJ_DISABLE_PCH = if (no_pch) "true" else "false")
  if (nzchar(cmdstan)) cmdstanr::set_cmdstan_path(path.expand(cmdstan))

  if (task == "check") {
    audit_check(root)
  } else if (task == "results") {
    audit_results(root, summaries)
  } else if (task == "aggregate") {
    audit_aggregate(root, collection)
  } else if (task %in% c("plan", "prepare", "spot-check", "sensitivity", "comparator")) {
    audit_simulations(root, collection, config_ids, replicates, task,
      parallel_chains, resume = resume)
  } else if (task %in% c("export-local", "local-status")) {
    audit_export_local(root, quiet = FALSE)
  } else if (task == "case-study") {
    audit_case_study(root, parallel_chains)
  } else if (task == "manifest") {
    audit_manifest(root)
  } else if (task == "validate") {
    validation <- new.env(parent = environment())
    sys.source(file.path(root, "R", "validate.R"), envir = validation)
    validation$audit_validate(root)
  }
}
