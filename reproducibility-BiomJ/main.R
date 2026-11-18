#!/usr/bin/env Rscript
# One entry point; paths are relative to this file, not to the working directory.
args <- commandArgs(trailingOnly = TRUE)
allowed <- c("mode", "collection", "config-ids", "replicates", "n-reps", "parallel-chains",
  "cmdstan", "no-pch", "summaries", "resume")
if (any(!grepl("^--[^=]+=", args)) || any(!sub("^--([^=]+)=.*$", "\\1", args) %in% allowed)) {
  stop("Unknown argument. Use named --option=value arguments documented in README.md.")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(dirname(script), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "audit.R"), local = TRUE)
source(file.path(root, "R", "simulation.R"), local = TRUE)
source(file.path(root, "R", "export-local.R"), local = TRUE)
source(file.path(root, "R", "entry.R"), local = TRUE)
source(file.path(root, "R", "session.R"), local = TRUE)
opt <- function(name, default = "") {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x) > 1) stop("Repeated option: ", name)
  if (length(x)) sub(paste0("^--", name, "="), "", x) else default
}
mode <- opt("mode", "check")
if (mode %in% c("prepare", "spot-check", "sensitivity", "comparator") && !any(grepl("^--config-ids=", args))) {
  stop("Select settings explicitly with --config-ids=1 (or other candidate IDs). Use --mode=plan to inspect settings without fitting.")
}
collection <- opt("collection", if (mode == "sensitivity") "alpha" else if (mode == "comparator") "comparator" else "main")
config_ids <- audit_parse_ids(opt("config-ids", "all"), "config-ids", 40L, allow_all = TRUE)
if (any(grepl("^--replicates=", args)) && any(grepl("^--n-reps=", args))) stop("Use --n-reps or --replicates, not both.")
n_reps <- audit_parse_ids(opt("n-reps", "5"), "n-reps", 1000L)
if (length(n_reps) != 1L) stop("--n-reps must be one integer.")
replicates <- audit_parse_ids(opt("replicates", paste(seq_len(n_reps), collapse = ",")), "replicates", 1000L)
resume_value <- opt("resume", "true")
if (!resume_value %in% c("true", "false")) stop("--resume must be true or false.")
chain_value <- opt("parallel-chains", "4")
if (!grepl("^[1-4]$", chain_value)) stop("--parallel-chains must be an integer from 1 to 4.")
parallel_chains <- as.integer(chain_value)
if (is.na(parallel_chains) || !parallel_chains %in% 1:4) stop("--parallel-chains must be 1, 2, 3 or 4.")
cmdstan <- opt("cmdstan", Sys.getenv("CMDSTAN", ""))
pch_value <- opt("no-pch", "false")
if (!pch_value %in% c("true", "false")) stop("--no-pch must be true or false.")
# The RStudio and terminal routes call the same dispatcher and fitting engine.
# The CLI parser already requires an explicit selection for fitting; explicit
# --config-ids=all is retained as an explicit all-settings selection.
if (mode %in% c("prepare", "spot-check", "sensitivity", "comparator") && is.null(config_ids)) {
  config_ids <- "all"
}
audit_run(root, mode, collection, config_ids, replicates, parallel_chains,
  cmdstan, pch_value == "true", resume_value == "true", opt("summaries", "saved"))
