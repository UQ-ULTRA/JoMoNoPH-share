# Local adapters around the numerically verified production definitions.
# Statistical calculations are retained; cluster dispatch is not included.
audit_collections <- function(root) {
  suffix <- c(
    main = "20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo",
    alpha125 = "20260528_ENZAMET-alphaSens-LL-maxF125-altDM-gamma5-1000rep-zeroprior-initinfo",
    alpha150 = "20260528_ENZAMET-alphaSens-LL-maxF150-altDM-gamma5-1000rep-zeroprior-initinfo",
    alpha300 = "20260528_ENZAMET-alphaSens-LL-maxF300-altDM-gamma5-1000rep-zeroprior-initinfo",
    comparator = "20260602_ENZAMET-JMbayes2-DM-all40-1000rep-randSDsummary"
  )
  data.frame(collection = names(suffix), suffix = unname(suffix),
    input = file.path(root, "archived-production-results", paste0("output_", suffix)),
    saved_summary = file.path(root, "archived-production-results", "auxiliary", substr(suffix, 1, 8), paste0("summary_tables_", suffix)),
    regenerated_summary = file.path(root, "output", "aggregated", paste0("summary_tables_", suffix)),
    stringsAsFactors = FALSE)
}

audit_parameters <- c("beta_long", "beta_surv", "alpha", "alpha_tilde", "gamma", "sigma_long", "sd_1_long", "L_1_long")

audit_parameter_summary <- function(draws, center) {
  sm <- posterior::summarise_draws(draws, posterior::default_summary_measures(),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail,
    quantiles = ~ posterior::quantile2(., probs = c(.025, .975)),
    mcse_mean = posterior::mcse_mean,
    mcse_q025 = ~ posterior::mcse_quantile(., probs = .025),
    mcse_q975 = ~ posterior::mcse_quantile(., probs = .975))
  shift <- sm$variable == "beta_long[1]"
  for (nm in intersect(c("mean", "median", "q5", "q95", "q2.5", "q97.5"), names(sm))) sm[[nm]][shift] <- sm[[nm]][shift] + center
  sm
}

audit_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these R packages first (see README): ", paste(missing, collapse = ", "))
}

audit_command <- function(command, args = character()) {
  tryCatch(suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE)), error = function(e) conditionMessage(e))
}

audit_environment <- function(out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  versions <- as.data.frame(installed.packages()[, c("Package", "Version", "Built")], stringsAsFactors = FALSE)
  write.csv(versions, file.path(out, "installed-packages.csv"), row.names = FALSE)
  stan_path <- if (requireNamespace("cmdstanr", quietly = TRUE)) tryCatch(cmdstanr::cmdstan_path(), error = function(e) NA_character_) else NA_character_
  details <- list(recorded_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    process_id = Sys.getpid(), command = commandArgs(FALSE),
    system = Sys.info(), platform = R.version, RNGkind = RNGkind(),
    os_version = if (Sys.info()[["sysname"]] == "Darwin") audit_command("sw_vers") else Sys.info(),
    locale = Sys.getlocale(), capabilities = capabilities(),
    cpu = if (Sys.info()[["sysname"]] == "Darwin") audit_command("sysctl", c("-n", "machdep.cpu.brand_string", "hw.model", "hw.memsize"))
    else if (Sys.info()[["sysname"]] == "Windows") audit_command("powershell", c("-NoProfile", "-Command", shQuote("Get-CimInstance Win32_Processor | Select-Object Name; Get-CimInstance Win32_ComputerSystem | Select-Object Model,TotalPhysicalMemory")))
    else audit_command("sh", c("-c", shQuote("uname -a; lscpu; free -b"))),
    # Motherboard identifiers may be unavailable in a VM or cluster allocation.
    motherboard = if (Sys.info()[["sysname"]] == "Darwin") audit_command("sysctl", c("-n", "hw.model")) else "Not automatically available; record manually if required.",
    cmdstan_path = stan_path,
    cmdstan_version = if (!is.na(stan_path)) as.character(cmdstanr::cmdstan_version()) else NA_character_,
    compiler_candidates = lapply(c("c++", "g++", "clang++", "make"), function(x) audit_command(x, "--version")),
    R_CPP_config = lapply(c("CXX", "CXX14", "CXX17", "CXXFLAGS", "CXX14FLAGS", "CXX17FLAGS", "LDFLAGS"), function(x) audit_command(file.path(R.home("bin"), "R"), c("CMD", "config", x))),
    compiler_environment = Sys.getenv(c("CXX", "CC", "CXXFLAGS", "CPPFLAGS", "LDFLAGS", "MAKE", "STAN_NUM_THREADS")))
  details$loaded_package_provenance <- lapply(loadedNamespaces(), function(p) {
    d <- utils::packageDescription(p)
    list(package = p, loaded_version = as.character(getNamespaceVersion(p)),
      repository = d$Repository, remote_repo = d$RemoteRepo, remote_ref = d$RemoteRef,
      remote_sha = d$RemoteSha, github_sha = d$GithubSHA1)
  })
  if (!is.na(stan_path)) {
    deps <- file.path(stan_path, "stan", "lib", "stan_math", "lib")
    details$included_libraries <- list.files(deps, pattern = "boost|eigen|tbb|sundials", ignore.case = TRUE)
    for (relative in c("make/local", "makefile", "stan/lib/stan_math/make/compiler_flags")) {
      f <- file.path(stan_path, relative)
      if (file.exists(f)) file.copy(f, file.path(out, gsub("/", "_", relative)), overwrite = TRUE)
    }
  }
  saveRDS(details, file.path(out, "environment.rds"))
  capture.output(str(details, max.level = 3), file = file.path(out, "environment.txt"))
}

# Load local fitting definitions into a private environment. No cluster dispatch.
# Inputs: package root and comparator flag. Returns: functions/settings environment.
audit_load_production <- function(root, comparator = FALSE) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(survival)
    library(MASS)
    library(nlme)
    library(cmdstanr)
  })
  e <- new.env(parent = globalenv())
  e$repo_root <- file.path(root, "production")
  e$stan_model_path <- file.path(root, "production", "Stan",
    "Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan")
  folder <- file.path(root, "production", "Bunya", "ENZAMET")
  sys.source(file.path(folder, "Data-Gen-ENZAMET.R"), envir = e)
  sys.source(file.path(folder, if (comparator) "ENZAMET-JMbayes2-DM.R" else
    "ENZAMET-DM.R"), envir = e)
  e
}

# Return the fixed 40-row configuration data frame; no fitting or file writes.
audit_grid <- function(root) {
  e <- new.env(parent = baseenv())
  sys.source(file.path(root, "production", "Bunya", "ENZAMET",
    "configuration-grid.R"), envir = e)
  stopifnot(nrow(e$config_grid) == 40L)
  e$config_grid
}

audit_model <- function(root, out, parallel_chains, prepare = FALSE) {
  function(stan_file, ...) {
    if (!file.exists(stan_file)) stop("Missing Stan source: ", stan_file)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    if (prepare) {
      dummy <- file.path(out, "preparation-only-not-a-Stan-binary")
      file.create(dummy)
      model <- NULL
    } else {
      start <- Sys.time()
      build_args <- list(...)
      if (Sys.getenv("BIOMJ_DISABLE_PCH", "false") == "true") {
        build_args$cpp_options$PRECOMPILED_HEADERS <- FALSE
      }
      build_log <- capture.output(model <- do.call(cmdstanr::cmdstan_model,
        c(list(stan_file = stan_file, quiet = FALSE), build_args)))
      writeLines(build_log, file.path(out, "compile.log"))
      saveRDS(list(seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
        stan_source = basename(stan_file), source_md5 = unname(tools::md5sum(stan_file)),
        executable_md5 = unname(tools::md5sum(model$exe_file())),
        cpp_options = build_args$cpp_options), file.path(out, "build.rds"))
      links <- if (Sys.info()[["sysname"]] == "Darwin") audit_command("otool", c("-L", shQuote(model$exe_file())))
      else if (Sys.info()[["sysname"]] == "Linux") audit_command("ldd", shQuote(model$exe_file())) else "See executable metadata and local compiler toolchain."
      writeLines(links, file.path(out, "linked-libraries.txt"))
    }
    list(exe_file = function() if (prepare) dummy else model$exe_file(),
      sample = function(...) {
        call <- list(...)
        caller <- parent.frame()
        center <- get0("longitudinal_response_center", envir = caller, inherits = TRUE, ifnotfound = 0)
        # Changing simultaneous chains changes resource use, not the four-chain
        # statistical settings. Retain all four exact initial-value lists.
        call$parallel_chains <- parallel_chains
        call$chain_ids <- seq_len(call$chains)
        call$output_dir <- file.path(out, "stan-csv")
        call$output_basename <- "audit"
        dir.create(call$output_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(list(data = call$data, init = call$init,
          longitudinal_response_center = center, call = call,
          RNGkind = RNGkind(), R_rng_state_at_sampling = .Random.seed), file.path(out, "inputs-and-initial-values.rds"))
        cmdstanr::write_stan_json(call$data, file.path(out, "stan-data.json"))
        for (chain in seq_along(call$init)) {
          cmdstanr::write_stan_json(call$init[[chain]], file.path(out, sprintf("initial-values-chain%d.json", chain)))
        }
        if (prepare) stop(structure(list(message = "Prepared inputs without sampling.", call = NULL), class = c("audit_prepared", "error", "condition")))
        start <- Sys.time()
        fit <- do.call(model$sample, call)
        elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
        draws <- fit$draws(variables = audit_parameters)
        sm <- audit_parameter_summary(draws, center)
        write.csv(sm, file.path(out, "posterior-summary-original-scale.csv"), row.names = FALSE)
        saveRDS(list(sampling_seconds = elapsed, diagnostics = fit$diagnostic_summary(),
          metadata = fit$metadata(), posterior_summary = sm), file.path(out, "fit-audit.rds"))
        # Selected coefficient draws are small; full raw CSVs remain available
        # for diagnostics without saving an R6 object with temporary-file paths.
        saveRDS(draws, file.path(out, "coefficient-draws-centered-scale.rds"))
        writeLines(readLines(fit$output_files()[1], n = 120, warn = FALSE), file.path(out, "stan-csv-header.txt"))
        fit
      })
  }
}

audit_compare_reference <- function(root, out, reference_result) {
  if (!file.exists(file.path(out, "inputs-and-initial-values.rds"))) return(invisible(NULL))
  ref <- reference_result
  inputs <- readRDS(file.path(out, "inputs-and-initial-values.rds"))
  center_diff <- inputs$longitudinal_response_center - ref$longitudinal_response_center
  beta_init <- inputs$data$beta_surv_init_mean
  check <- data.frame(item = c("response_center_difference", "exponential_beta_init_difference"),
    difference = c(center_diff, as.numeric(beta_init - ref$prior_info$beta_surv_init_mean)))
  write.csv(check, file.path(out, "data-regeneration-check.csv"), row.names = FALSE)
  if (any(!is.finite(check$difference)) || any(abs(check$difference) > 1e-7)) stop("Generated data/staging fits do not match the selected historical replicate; inspect ", out)
  fields <- intersect(names(inputs$init[[1]]), names(ref$init_info$values[[1]]))
  fields <- fields[vapply(fields, function(n) is.numeric(inputs$init[[1]][[n]]) && is.numeric(ref$init_info$values[[1]][[n]]), logical(1))]
  init_diff <- do.call(rbind, lapply(seq_along(inputs$init), function(chain) {
    data.frame(chain = chain, parameter = fields, max_absolute_difference = vapply(fields,
      function(n) max(abs(as.numeric(inputs$init[[chain]][[n]]) - as.numeric(ref$init_info$values[[chain]][[n]]))), numeric(1)))
  }))
  write.csv(init_diff, file.path(out, "initial-values-regeneration-check.csv"), row.names = FALSE)
  if (any(!is.finite(init_diff$max_absolute_difference)) || any(init_diff$max_absolute_difference > 1e-7)) stop("Historical compact initial values do not match: ", out)
  smfile <- file.path(out, "posterior-summary-original-scale.csv")
  if (file.exists(smfile)) {
    fresh <- read.csv(smfile, check.names = FALSE)
    old <- as.data.frame(ref$joint_summary)
    x <- merge(fresh, old[, intersect(c("variable", "mean", "sd", "ess_bulk", "rhat"), names(old))], by = "variable", suffixes = c("_fresh", "_production"))
    x$mean_difference <- x$mean_fresh - x$mean_production
    # Historical summaries lack mcse_mean. Do not invent exact MCSEs or mark a
    # fit as passing solely because rounded point estimates happen to match.
    write.csv(x, file.path(out, "comparison-with-production.csv"), row.names = FALSE)
  }
}

audit_case_study <- function(root, parallel_chains = 4) {
  audit_require(c("tidyverse", "survival", "nlme", "cmdstanr", "posterior"))
  out <- file.path(root, "output", "case-study")
  source_dir <- file.path(root, "case-study", "source")
  audit_environment(file.path(out, "environment"))
  # This supplied example used R's default RNG, unlike the furrr simulation.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  e <- new.env(parent = globalenv())
  e$cmdstan_model <- function(stan_file, ...) audit_model(root, out, parallel_chains)(file.path(root, "production", "Stan", basename(stan_file)), ...)
  e$readRDS <- function(file, ...) base::readRDS(file.path(source_dir, basename(file)), ...)
  e$write.csv <- function(x, file, ...) utils::write.csv(x, file.path(out, basename(file)), ...)
  e$writeLines <- function(text, con, ...) base::writeLines(text, file.path(out, basename(con)), ...)
  e$saveRDS <- function(object, file, ...) {
    if (basename(file) == "mock_fit_objects.rds") {
      # Raw Stan CSVs and selected draws are saved by audit_model(). Do not
      # serialize an R6 fit pointing to temporary files that disappear later.
      object$fit <- NULL
    }
    base::saveRDS(object, file.path(out, basename(file)), ...)
  }
  start <- Sys.time()
  for (expr in parse(file.path(source_dir, "fit-case-study-model.R"))) {
    eval(expr, e)
  }
  write.csv(data.frame(elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs"))), file.path(out, "elapsed.csv"), row.names = FALSE)
}

# Return a private script environment with explicit dplyr bindings, rather than
# changing the user's global conflicted-package preferences. No file writes.
audit_analysis_environment <- function() {
  e <- new.env(parent = globalenv())
  for (name in c("select", "filter", "mutate", "arrange", "summarise")) {
    e[[name]] <- getExportedValue("dplyr", name)
  }
  e
}

audit_results <- function(root, summaries = "saved") {
  if (!summaries %in% c("saved", "regenerated")) stop("--summaries must be saved or regenerated.")
  Sys.setenv(BIOMJ_ROOT = root, BIOMJ_SUMMARIES = summaries,
    ENZ_MANUSCRIPT_DIR = file.path(root, "output", "manuscript"))
  audit_environment(file.path(root, "output", "manuscript", "environment"))
  sys.source(file.path(root, "analysis", "generate_manuscript_outputs.R"), envir = audit_analysis_environment())
}

audit_aggregate <- function(root, collection) {
  Sys.setenv(BIOMJ_ROOT = root)
  collections <- audit_collections(root)
  if (collection == "all") collection <- collections$collection
  if (any(!collection %in% collections$collection)) stop("Unknown --collection.")
  for (id in collection) {
    x <- collections[collections$collection == id, ]
    Sys.setenv(SIM_OUTPUT_DIR = x$input, SIM_SUMMARY_DIR = x$regenerated_summary,
      SIM_FIGURE_DIR = file.path(root, "output", "aggregated", paste0("figures_", x$suffix)),
      SIM_BATCH_SIZE = if (id == "comparator") "10" else "4", SIMS_PER_CONFIG = "1000")
    f <- if (id == "comparator") "SumSimRes-ENZAMET-JMbayes2-DM.R" else "SumSimRes-ENZAMET-DM.R"
    fig <- Sys.getenv("SIM_FIGURE_DIR")
    dir.create(fig, recursive = TRUE, showWarnings = FALSE)
    # Production scripts print diagnostic plots as well as calling ggsave().
    # An explicit device prevents R from creating Rplots.pdf in the caller's cwd.
    grDevices::pdf(file.path(fig, "console-plots.pdf"))
    tryCatch(sys.source(file.path(root, "analysis", f), envir = audit_analysis_environment()),
      finally = grDevices::dev.off())
  }
}

audit_check <- function(root) {
  x <- audit_collections(root)
  x$expected_batches <- c(10000L, 2500L, 2500L, 2500L, 4000L)
  x$found_batches <- vapply(x$input, function(p) length(list.files(p, "\\.rds$", recursive = FALSE)), integer(1))
  x$summary_available <- dir.exists(x$saved_summary)
  print(x[, c("collection", "expected_batches", "found_batches", "summary_available")], row.names = FALSE)
  if (any(x$expected_batches != x$found_batches) || !all(x$summary_available)) stop("Result archive incomplete.")
  if (!file.exists(file.path(root, "case-study", "source", "mock-data.rds"))) stop("Missing mock data.")
  # Include tidyverse components used by the scripts in the loaded session record.
  p <- c("tidyverse", "dplyr", "ggplot2", "purrr", "readr", "stringr", "tibble",
         "tidyr", "forcats", "lubridate", "survival", "nlme", "MASS", "cmdstanr",
         "posterior", "loo", "patchwork", "JMbayes2", "coda")
  versions <- data.frame(package = p, installed = vapply(p, requireNamespace, logical(1), quietly = TRUE))
  print(versions, row.names = FALSE)
  audit_environment(file.path(root, "output", "check"))
  write.csv(audit_grid(root), file.path(root, "output", "check", "configurations.csv"), row.names = FALSE)
  message("Archive count check passed. This does not certify every fit's diagnostics.")
}

audit_manifest <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[!file.info(files)$isdir & !grepl("/output/|/files_manifest.csv$|\\.DS_Store$", files)]
  message("Hashing ", length(files), " files; the result archive can take several minutes.")
  rel <- substring(files, nchar(root) + 2)
  tab <- data.frame(path = rel, bytes = file.info(files)$size, md5 = unname(tools::md5sum(files)))
  write.csv(tab, file.path(root, "files_manifest.csv"), row.names = FALSE)
}
