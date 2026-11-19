#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)

# Stubs exercise routing and validation without compiling or sampling any model.
e <- new.env(parent = globalenv())
sys.source(file.path(root, "R", "session.R"), envir = e)
sys.source(file.path(root, "R", "entry.R"), envir = e)
e$audit_simulations <- function(root, collection, config_ids, replicates, mode,
                                parallel_chains, resume) {
  list(root = root, collection = collection, config_ids = config_ids,
    replicates = replicates, mode = mode, parallel_chains = parallel_chains,
    resume = resume, no_pch = Sys.getenv("BIOMJ_DISABLE_PCH"))
}
e$audit_check <- function(root, cmdstan = "") "check"
e$audit_results <- function(root, summaries) summaries
e$audit_aggregate <- function(root, collection) collection
e$audit_export_local <- function(root, quiet) list(task = "records", quiet = quiet)
e$audit_case_study <- function(root, parallel_chains) parallel_chains
bad <- function(...) inherits(try(e$audit_run(root, ...), silent = TRUE), "try-error")
for (collection in c("main", "alpha125", "alpha150", "alpha300", "comparator")) {
  x <- e$audit_run(root, "fit", collection, c(1, 2), c(7, 42), 2L)
  stopifnot(x$mode == "spot-check", x$collection == collection,
    identical(x$config_ids, c(1L, 2L)), identical(x$replicates, c(7L, 42L)),
    x$parallel_chains == 2L, x$resume, x$no_pch == "false")
}
stopifnot(e$audit_run(root, "check") == "check",
  e$audit_run(root, "results", summaries = "regenerated") == "regenerated",
  e$audit_run(root, "aggregate", collection = "all") == "all",
  e$audit_run(root, "records")$task == "records",
  e$audit_run(root, "case-study", parallel_chains = 2L) == 2L,
  is.null(e$audit_run(root, "fit", config_ids = "all")$config_ids),
  is.null(e$audit_run(root, "plan")$config_ids),
  bad(task = "fit"), bad(task = "prepare"),
  bad(config_ids = numeric()), bad(config_ids = c(1, 1)),
  bad(config_ids = NA), bad(config_ids = 1.5), bad(config_ids = 41),
  bad(replicates = 1001), bad(replicates = c(1, 1)),
  bad(parallel_chains = 5), bad(parallel_chains = c(1, 2)),
  bad(resume = NA), bad(no_pch = "true"), bad(summaries = NA),
  bad(collection = "typo"), bad(task = "typo"))

local({
  old <- Sys.getenv("BIOMJ_DISABLE_PCH", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("BIOMJ_DISABLE_PCH") else
    Sys.setenv(BIOMJ_DISABLE_PCH = old))
  Sys.unsetenv("BIOMJ_DISABLE_PCH")
  stopifnot(e$audit_run(root, no_pch = TRUE)$no_pch == "true",
    is.na(Sys.getenv("BIOMJ_DISABLE_PCH", unset = NA_character_)))
  Sys.setenv(BIOMJ_DISABLE_PCH = "previous-setting")
  e$audit_run(root, no_pch = FALSE)
  stopifnot(Sys.getenv("BIOMJ_DISABLE_PCH") == "previous-setting")
})

# Errors must also restore temporary settings and the user's RNG stream.
local({
  set.seed(901)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  env_before <- Sys.getenv("SIM_BIOMJ_TEST", unset = NA_character_)
  option_before <- getOption("biomj_test_option")
  previous <- e$audit_simulations
  e$audit_simulations <- function(...) {
    Sys.setenv(SIM_BIOMJ_TEST = "temporary")
    options(biomj_test_option = "temporary")
    RNGkind("L'Ecuyer-CMRG")
    runif(3)
    stop("Deliberate restoration test")
  }
  stopifnot(inherits(try(e$audit_run(root, "plan"), silent = TRUE), "try-error"),
    identical(.Random.seed, seed_before), identical(RNGkind(), kind_before),
    identical(Sys.getenv("SIM_BIOMJ_TEST", unset = NA_character_), env_before),
    identical(getOption("biomj_test_option"), option_before))
  e$audit_simulations <- previous
})

# Execute the entry script in both supported project directories, forcing its
# task to plan in memory. Users may have edited it to fit/prepare; tests must
# never launch those tasks or overwrite the user's choice in the source file.
local({
  cwd <- getwd()
  on.exit(setwd(cwd))
  expressions <- parse(file.path(root, "run_checks.R"))
  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("task"))) expressions[[i]] <- quote(task <- "plan")
  }
  for (wd in c(root, dirname(root))) {
    setwd(wd)
    session <- new.env(parent = globalenv())
    invisible(capture.output(eval(expressions, envir = session)))
    expected <- session$audit_simulation_plan(root, session$collection, session$config_ids, session$replicates)
    stopifnot(identical(session$audit_result, expected), nrow(session$configurations) == 40L,
      identical(getwd(), wd), session$task == "plan")
  }
})
# The case-study table must also work after a simulation has loaded MASS.
# Evaluate only its formatting expressions with fixtures, not fitting or I/O.
local({
  tab_env <- new.env(parent = asNamespace("dplyr"))
  tab_env$tibble <- tibble::tibble
  tab_env$replace_na <- tidyr::replace_na
  expressions <- parse(file.path(root, "case-study", "source", "fit-case-study-model.R"))
  wanted <- c("format_est_ci", "uncentre", "labels", "table1_update")
  assignments <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    is.symbol(x[[2]]) && as.character(x[[2]]) %in% wanted, as.list(expressions))
  for (x in assignments[-length(assignments)]) eval(x, tab_env)
  tab_env$longitudinal_response_center <- 73
  tab_env$saft <- tibble::tibble(variable = names(tab_env$labels), mean = 1,
    lower = 0, upper = 2)
  tab_env$lmm <- tab_env$saft[1:6, ]
  table_expr <- assignments[[length(assignments)]]
  tab_env$select <- dplyr::select
  expected <- eval(table_expr, tab_env)
  tab_env$select <- MASS::select
  actual <- eval(table_expr, tab_env)
  stopifnot(identical(actual, expected), nrow(actual) == 13L,
    identical(names(actual), c("parameter", "variable", "sAFT_JM", "LMM")))
})
cat("PASS: RStudio routing, settings restoration and case-study table with masked select(). No MCMC run.\n")
