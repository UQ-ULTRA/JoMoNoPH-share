# Local dispatch only: the data-generation and fitting functions remain the
# numerically verified local production functions loaded by audit_load_production().
audit_parse_ids <- function(value, name, upper, allow_all = FALSE) {
  if (allow_all && identical(value, "all")) return(NULL)
  pieces <- strsplit(value, ",", fixed = TRUE)[[1]]
  if (!length(pieces) || any(!grepl("^[1-9][0-9]*$", pieces))) stop("--", name, " must be comma-separated positive integer IDs.")
  ids <- suppressWarnings(as.integer(pieces))
  if (anyNA(ids) || any(ids > upper) || anyDuplicated(ids)) stop("Invalid or repeated --", name, " IDs; maximum is ", upper, ".")
  ids
}

audit_simulation_plan <- function(root, collection = "main", config_ids = NULL,
                                  replicates = 1:5) {
  collections <- audit_collections(root)
  if (identical(collection, "all")) collection <- collections$collection
  if (identical(collection, "alpha")) collection <- c("alpha125", "alpha150", "alpha300")
  if (!length(collection) || any(!collection %in% collections$collection)) stop("Unknown collection. Use main, alpha125, alpha150, alpha300, alpha, comparator or all.")
  grid <- audit_grid(root)
  pieces <- lapply(collection, function(id) {
    g <- grid
    if (grepl("^alpha", id)) g <- g[g$baseline_label == "LL120_230", , drop = FALSE]
    if (!is.null(config_ids)) {
      if (any(!config_ids %in% g$full_config_id)) stop("Some --config-ids are not available in collection ", id, ". Alpha collections use full-grid IDs 1 to 10.")
      g <- g[match(config_ids, g$full_config_id), , drop = FALSE]
    }
    p <- g[rep(seq_len(nrow(g)), each = length(replicates)), , drop = FALSE]
    p$collection <- id
    p$replicate_id <- rep(replicates, times = nrow(g))
    p$unique_seed <- 1e6 * p$full_config_id + p$replicate_id
    p$alpha_max_f <- unname(c(main = 2, alpha125 = 1.25, alpha150 = 1.5, alpha300 = 3, comparator = NA_real_)[id])
    batch_size <- if (id == "comparator") 10L else 4L
    p$reference_batch <- ceiling(p$replicate_id / batch_size)
    p$reference_position <- (p$replicate_id - 1L) %% batch_size + 1L
    prefix <- if (id == "comparator") "JMbayes2_All" else "BP"
    filename <- sprintf("%s_k%d_n%d_cens%s_%s_scen%d_FU%d_batch%03d.rds", prefix,
      p$k_bases, p$n_patients, p$cens_label, p$baseline_label, p$scenario, p$max_FU, p$reference_batch)
    archive <- basename(collections$input[match(id, collections$collection)])
    p$reference_file <- file.path("archived-production-results", archive, filename)
    p$output_key <- file.path(id, sprintf("config%03d", p$full_config_id), sprintf("rep%03d", p$replicate_id))
    p
  })
  plan <- do.call(rbind, pieces)
  rownames(plan) <- NULL
  stopifnot(!anyDuplicated(plan$output_key))
  plan
}

audit_plan_reference <- function(root, row) {
  f <- file.path(root, row$reference_file)
  if (!file.exists(f)) stop("Missing production reference: ", f)
  x <- readRDS(f)[[row$reference_position]]
  if (is.null(x) || !identical(as.numeric(x$unique_seed), as.numeric(row$unique_seed))) stop("Missing/wrong production replicate: ", f)
  x
}

audit_comparator_check <- function(out, reference, fresh = NULL) {
  if (is.null(fresh)) fresh <- readRDS(file.path(out, "comparator-prepared.rds"))
  d <- data.frame(item = "response_center_difference", difference = fresh$longitudinal_response_center - reference$longitudinal_response_center)
  write.csv(d, file.path(out, "data-regeneration-check.csv"), row.names = FALSE)
  if (!is.finite(d$difference) || abs(d$difference) > 1e-7) stop("Comparator dataset/staging check differs from production: ", out)
  if (!is.null(fresh$comparison_tidy)) {
    a <- as.data.frame(fresh$comparison_tidy)
    b <- as.data.frame(reference$comparison_tidy)
    keys <- c("model", "variable")
    a <- a[a$variable %in% c("beta_long[1]", "beta_long[2]", "beta_long[3]", "sigma_long", "sd_1_long[1]", "sd_1_long[2]", "alpha", "beta_surv[1]"), ]
    x <- merge(a[, c(keys, "mean")], b[, c(keys, "mean")], by = keys, all = TRUE, suffixes = c("_fresh", "_production"))
    x$mean_difference <- x$mean_fresh - x$mean_production
    write.csv(x, file.path(out, "comparison-with-production.csv"), row.names = FALSE)
  }
}

# Build flags cannot change preparation-only data or initial values because
# preparation never compiles/runs Stan. Keep checking them for actual Stan fits.
# Compare copies only: the saved checkpoint retains its original provenance.
audit_checkpoint_differences <- function(saved, current, preparation_only = FALSE,
                                         comparator = FALSE, source_aliases = NULL) {
  if (!is.list(saved) || !is.list(current)) return("identity")
  # Accept only explicitly verified source-only cleanup mappings. Do not modify
  # the historical checkpoint, its software versions, priors, seed or data IDs.
  if (!is.null(source_aliases) && !identical(saved$source_md5, current$source_md5)) {
    previous <- saved$source_md5
    idx <- match(previous, source_aliases$before_md5)
    mapped <- previous
    mapped[!is.na(idx)] <- source_aliases$after_md5[idx[!is.na(idx)]]
    if ("extra_md5" %in% names(source_aliases)) {
      extra <- source_aliases$extra_md5[idx[!is.na(idx)]]
      extra <- unique(extra[!is.na(extra) & nzchar(extra)])
      mapped <- c(mapped, extra)
    }
    if (identical(mapped, current$source_md5)) saved$source_md5 <- mapped
  }
  fields <- union(names(saved), names(current))
  if (preparation_only || comparator) fields <- setdiff(fields, "no_pch")
  fields[!vapply(fields, function(field) identical(saved[[field]], current[[field]]), logical(1))]
}

audit_simulations <- function(root, collection = "main", config_ids = NULL,
                              replicates = 1:5, mode = "spot-check",
                              parallel_chains = 4, resume = TRUE,
                              continue_on_error = FALSE) {
  plan <- audit_simulation_plan(root, collection, config_ids, replicates)
  plan_dir <- file.path(root, "output", "simulation-plans")
  dir.create(plan_dir, recursive = TRUE, showWarnings = FALSE)
  audit_write_csv(audit_grid(root), file.path(plan_dir, "configuration_candidates.csv"))
  write.csv(plan, file.path(plan_dir, paste0(collection, ".csv")), row.names = FALSE)
  print(as.data.frame(table(plan$collection), responseName = "dataset_runs"), row.names = FALSE)
  message("Plan: ", nrow(plan), " dataset runs; comparator runs each fit three JMs plus an LMM.")
  if (mode == "plan") return(invisible(plan))
  if (mode != "prepare") {
    active_file <- file.path(root, "output", "local-simulations", "run-plan.csv")
    if (file.exists(active_file)) {
      active <- read.csv(active_file, stringsAsFactors = FALSE)
      active <- active[!active$output_key %in% plan$output_key, , drop = FALSE]
      plan_union <- rbind(active[, names(plan), drop = FALSE], plan)
    } else plan_union <- plan
    audit_write_csv(plan_union, active_file)
  }
  audit_require(c("tidyverse", "survival", "nlme", "MASS", "cmdstanr", "posterior", "loo"))
  if ("comparator" %in% plan$collection) audit_require("JMbayes2")
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  Sys.setenv(SIM_STAN_MODEL_FILE = file.path(root, "production", "Stan", "Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"),
    SIM_ALPHA_DELTA_Y = "1.96", SIM_JMB2_N_CHAINS = "4", SIM_JMB2_CORES = as.character(parallel_chains),
    SIM_JMB2_N_ITER = "3500", SIM_JMB2_N_BURNIN = "500", SIM_JMB2_N_THIN = "1", SIM_JMB2_MODELS = "all",
    SIM_JMB2_SAVE_MODE = "summary", SIM_JMB2_RETRY_CORES = "0")
  for (id in unique(plan$collection)) {
    run_root <- file.path(root, "output", "local-simulations", id)
    audit_environment(file.path(run_root, "environment"))
    for (j in which(plan$collection == id)) {
      p <- plan[j, , drop = FALSE]
      out <- file.path(root, "output", "local-simulations", if (mode == "prepare") "prepared" else "fits", p$output_key)
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      error <- tryCatch(
        {
          marker <- file.path(out, "completed.rds")
          fingerprint <- unname(tools::md5sum(c(file.path(root, "production", "Bunya", "ENZAMET", if (id == "comparator") "ENZAMET-JMbayes2-DM.R" else "ENZAMET-DM.R"),
            file.path(root, "production", "Bunya", "ENZAMET", "Data-Gen-ENZAMET.R"),
            file.path(root, "production", "Stan", "Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"),
            file.path(root, "production", "Bunya", "ENZAMET", "configuration-grid.R"))))
          identity <- list(collection = id, config = p$full_config_id, replicate = p$replicate_id,
            seed = p$unique_seed, alpha_max_f = p$alpha_max_f, source_md5 = fingerprint,
            R_version = as.character(getRversion()),
            package_versions = vapply(c("nlme", "survival", "MASS", if (id == "comparator") "JMbayes2" else c("cmdstanr", "posterior")),
              function(x) as.character(utils::packageVersion(x)), character(1)),
            cmdstan_version = if (id == "comparator") NA_character_ else as.character(cmdstanr::cmdstan_version()),
            no_pch = Sys.getenv("BIOMJ_DISABLE_PCH", "false"))
          if (resume && file.exists(marker)) {
            saved <- readRDS(marker)
            expected <- if (mode == "prepare") if (id == "comparator") "comparator-prepared.rds" else "inputs-and-initial-values.rds" else "production-format-result.rds"
            aliases_file <- file.path(root, "documentation", "verified-source-cleanup.csv")
            aliases <- if (file.exists(aliases_file)) read.csv(aliases_file) else NULL
            differences <- audit_checkpoint_differences(saved$identity, identity,
              preparation_only = mode == "prepare", comparator = id == "comparator",
              source_aliases = aliases)
            missing <- !file.exists(file.path(out, expected))
            if (length(differences) || missing) {
              detail <- c(if (length(differences)) paste0("different fields: ", paste(differences, collapse = ", ")),
                if (missing) paste0("missing file: ", expected))
              stop("Incomplete or incompatible checkpoint at ", out, " (", paste(detail, collapse = "; "),
                "). Inspect it first; to regenerate explicitly, set resume <- FALSE in run_checks.R (terminal: --resume=false).")
            }
            message("Skipping completed ", p$output_key)
            next
          }
          reference <- audit_plan_reference(root, p)
          write.csv(p, file.path(out, "reference-provenance.csv"), row.names = FALSE)
          if (id != "comparator") Sys.setenv(SIM_ALPHA_MAX_F = as.character(p$alpha_max_f))
          e <- audit_load_production(root, comparator = id == "comparator")
          e$parallel_chains <- parallel_chains
          if (id != "comparator") e$cmdstan_model <- audit_model(root, out, parallel_chains, prepare = mode == "prepare")
          if (id == "comparator") {
            original_fit <- e$fit_one_jmb2_model
            # Record the generated observations and LMM staging inputs for every
            # comparator run. Preparation stops before the first JM fit.
            e$fit_one_jmb2_model <- function(model, fit_long, surv_data, unique_seed, sim_id, longitudinal_response_center = 0) {
              saveRDS(list(longitudinal = nlme::getData(fit_long), survival = surv_data,
                unique_seed = unique_seed, longitudinal_response_center = longitudinal_response_center,
                fixed_effects = nlme::fixef(fit_long), fixed_vcov = stats::vcov(fit_long),
                random_effects = nlme::ranef(fit_long), residual_sd = fit_long$sigma,
                random_effect_covariance = as.matrix(nlme::getVarCov(fit_long)),
                jm_seed_by_model = unique_seed + setNames(1:3, c("Cox", "Weibull", "loglogistic")),
                jm_control = list(n_iter = 3500L, n_burnin = 500L, n_thin = 1L,
                  n_chains = 4L, cores = parallel_chains, save_random_effects = FALSE,
                  save_logLik_contributions = FALSE)), file.path(out, "comparator-prepared.rds"))
              if (mode == "prepare") stop(structure(list(message = "Prepared comparator inputs without sampling.", call = NULL), class = c("audit_prepared", "error", "condition")))
              original_fit(model, fit_long, surv_data, unique_seed, sim_id, longitudinal_response_center)
            }
          }
          call <- list(i = p$replicate_id, n_patients = p$n_patients, lambda_c = p$lambda_c,
            aft_mode = p$aft_mode, max_FU = p$max_FU, config_hash = p$full_config_id,
            k_bases = p$k_bases, scenario = p$scenario, baseline_label = p$baseline_label,
            cfg_weibull_shape = p$weibull_shape, cfg_weibull_scale = p$weibull_scale,
            cfg_loglogistic_shape = p$loglogistic_shape, cfg_loglogistic_scale = p$loglogistic_scale)
          saveRDS(call, file.path(out, "configuration.rds"))
          message("Running ", mode, ": ", p$output_key, "; original seed ", p$unique_seed)
          start <- Sys.time()
          # furrr resets its RNG stream for each production replicate. Do the same
          # even if a preceding comparator package changed the parent RNG kind.
          RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
          warning_file <- file.path(out, "warnings.log")
          base::cat("\nAttempt started: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", file = warning_file, append = TRUE, sep = "")
          result <- tryCatch(withCallingHandlers(do.call(e$run_simulation, call), warning = function(w) {
            text <- conditionMessage(w)
            base::cat(text, "\n", file = warning_file, append = TRUE)
            message("WARNING [", p$output_key, "]: ", text)
            # Evidence is retained both per replicate and immediately in the run
            # log, rather than queued until the multi-day top-level call returns.
            invokeRestart("muffleWarning")
          }), error = function(err) {
            cause <- err
            while (!is.null(cause)) {
              if (inherits(cause, "audit_prepared")) return(NULL)
              cause <- cause$parent
            }
            stop(err)
          })
          if (mode != "prepare") {
            if (is.null(result)) stop("Fit returned NULL: ", out)
            saveRDS(result, file.path(out, "production-format-result.rds"))
          }
          if (id == "comparator") audit_comparator_check(out, reference, result)
          else audit_compare_reference(root, out, reference_result = reference)
          if (mode != "prepare") audit_record_one(root, p, out, result)
          write.csv(data.frame(elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs"))), file.path(out, "elapsed.csv"), row.names = FALSE)
          # This records completed execution, not acceptable diagnostics or success
          # of every comparator model. Its original model_status is retained.
          saveRDS(list(identity = identity, completed_at = Sys.time(), preparation_only = mode == "prepare"), marker)
          NULL
        },
        error = function(err) err)
      if (inherits(error, "error")) {
        saveRDS(list(message = conditionMessage(error), at = Sys.time(), plan = p), file.path(out, "failure.rds"))
        message("FAILED ", p$output_key, ": ", conditionMessage(error))
        if (!continue_on_error) {
          if (mode != "prepare") audit_export_local(root)
          stop(error)
        }
      }
      if (mode != "prepare" && exists("audit_export_local", mode = "function")) audit_export_local(root)
      invisible(gc(verbose = FALSE))
    }
  }
  if (mode != "prepare" && exists("audit_export_local", mode = "function")) audit_export_local(root)
  invisible(plan)
}
