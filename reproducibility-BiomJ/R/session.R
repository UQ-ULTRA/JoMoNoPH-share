# Capture task settings; return a function that restores them on exit.
# No files are written. Restore RNG only after all inputs/results are saved.
audit_session <- function() {
  old_env <- Sys.getenv()
  old_options <- options()
  old_rng <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  old_plan <- if ("future" %in% loadedNamespaces()) future::plan() else NULL
  old_cmdstan <- if (requireNamespace("cmdstanr", quietly = TRUE)) {
    tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  } else NULL
  function() {
    current <- Sys.getenv()
    keys <- grep("^(SIM_|BIOMJ_|ENZ_MANUSCRIPT_DIR$)", union(names(current), names(old_env)), value = TRUE)
    Sys.unsetenv(setdiff(keys, names(old_env)))
    keep <- intersect(keys, names(old_env))
    if (length(keep)) do.call(Sys.setenv, as.list(old_env[keep]))
    if (!is.null(old_plan)) future::plan(old_plan)
    new_options <- setdiff(names(options()), names(old_options))
    if (length(new_options)) options(setNames(rep(list(NULL), length(new_options)), new_options))
    options(old_options)
    if (!is.null(old_cmdstan)) suppressMessages(cmdstanr::set_cmdstan_path(old_cmdstan))
    do.call(RNGkind, as.list(old_rng))
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}
