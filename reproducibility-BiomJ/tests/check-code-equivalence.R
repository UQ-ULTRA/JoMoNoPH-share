# Optional author-environment regression check, not a cross-platform MCMC test.
# Compare against the recorded pre-cleanup hashes without running MCMC.
stage <- "after"
root <- if (file.exists("R/audit.R")) normalizePath(".") else normalizePath("reproducibility-BiomJ")
source(file.path(root, "R", "audit.R"))
source(file.path(root, "R", "simulation.R"))
out <- file.path(root, "output", "compliance-checks")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
hash <- function(x) {
  f <- tempfile(tmpdir = out)
  on.exit(unlink(f))
  saveRDS(x, f, version = 3)
  unname(tools::md5sum(f))
}
rows <- list()
record <- function(kind, config, rep, value) {
  rows[[length(rows) + 1L]] <<- data.frame(kind = kind, config = config,
    replicate = rep, md5 = hash(value))
}
grid <- audit_grid(root)
for (j in seq_len(nrow(grid))) for (rep in c(1L, 56L)) {
  p <- grid[j, ]
  e <- audit_load_production(root)
  list2env(e$simulation_parameters(p$scenario), envir = e)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  d <- e$simulate_joint_dataset(e$D, e$beta_0, e$beta_1, e$beta_2, e$sigma_e,
    e$log_HR, 99, e$log_AF, e$alpha_AFT, p$weibull_shape, p$weibull_scale,
    p$loglogistic_shape, p$loglogistic_scale, e$visit, 1e6 * j + rep,
    p$n_patients, p$max_FU, p$lambda_c, p$aft_mode, "value")
  record("generated-data", j, rep, list(data = d, RNG = .Random.seed))
  message(stage, " data config ", j, " rep ", rep)
}

# Capture exactly what reaches Stan, including all four random-effect starts.
dummy <- tempfile(tmpdir = out)
file.create(dummy)
on_stop <- function() stop(structure(list(message = "Captured inputs", call = NULL),
  class = c("equivalence_capture", "error", "condition")))
captured <- NULL
capture_model <- function(...) list(exe_file = function() dummy, sample = function(...) {
  x <- list(...)
  captured <<- list(data = x$data, init = x$init, RNG = .Random.seed,
    center = get0("longitudinal_response_center", parent.frame(), inherits = TRUE))
  on_stop()
})
cases <- rbind(data.frame(config = 1:40, rep = 1L, factor = 2),
  data.frame(config = 1L, rep = 56L, factor = c(1.25, 1.5, 3)))
for (k in seq_len(nrow(cases))) {
  s <- cases[k, ]
  p <- grid[s$config, ]
  Sys.setenv(SIM_ALPHA_MAX_F = s$factor, SIM_ALPHA_DELTA_Y = "1.96",
    SIM_STAN_MODEL_FILE = file.path(root, "production/Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"))
  e <- audit_load_production(root)
  e$cmdstan_model <- capture_model
  captured <- NULL
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  tryCatch(e$run_simulation(s$rep, p$n_patients, p$lambda_c, p$aft_mode, p$max_FU,
    s$config, k_bases = p$k_bases, scenario = p$scenario, baseline_label = p$baseline_label,
    cfg_weibull_shape = p$weibull_shape, cfg_weibull_scale = p$weibull_scale,
    cfg_loglogistic_shape = p$loglogistic_shape, cfg_loglogistic_scale = p$loglogistic_scale),
  error = function(err) if (is.null(captured)) stop(err))
  stopifnot(!is.null(captured))
  record(paste0("stan-inputs-alpha", s$factor), s$config, s$rep, captured)
  message(stage, " Stan inputs config ", s$config, " alpha ", s$factor)
}
unlink(dummy)
result <- do.call(rbind, rows)
write.csv(result, file.path(out, paste0(stage, "-hashes.csv")), row.names = FALSE)
if (stage == "after") {
  before <- read.csv(file.path(out, "before-hashes.csv"))
  stopifnot(identical(before, result))
  message("PASS: generated data, Stan inputs/initial values and RNG states are identical.")
}
