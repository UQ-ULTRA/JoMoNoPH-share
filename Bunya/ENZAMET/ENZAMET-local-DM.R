#!/usr/bin/env Rscript

# Local launcher for ENZAMET-DM.R.
#
# This script emulates the SLURM environment used by ENZAMET-DM.R so one batch
# can be run locally on macOS while producing the same batch-level .rds structure
# expected by Bunya/SumSimRes-ENZAMET-DM.R.
#
# Examples:
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 1 --cpus 4
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 126 --cpus 8 --output-dir Bunya/ENZAMET/output_local
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 1 --batch-size 1 --sims-per-config 1
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 1 --batch-size 1 --sims-per-config 1 --seed-offset 4
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 1 --batch-size 1 --sims-per-config 1 --fixed-seed 5
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 1 --force-recompile
#   Rscript Bunya/ENZAMET/ENZAMET-local-DM.R --task-id 1 --stan-model-file Stan/Bernstein-Polynomials-JM-Hist-Penalty.stan --force-recompile

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1]]
}

has_flag <- function(flag) {
  flag %in% args
}

as_int_arg <- function(flag, default) {
  value <- arg_value(flag, default)
  value <- suppressWarnings(as.integer(value))
  if (is.na(value)) {
    stop(sprintf("Argument %s must be an integer.", flag), call. = FALSE)
  }
  value
}

file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_args) > 0) {
  this_file <- sub("^--file=", "", file_args[[1]])
} else {
  this_file <- "Bunya/ENZAMET/ENZAMET-local-DM.R"
}
this_file <- normalizePath(this_file, mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(this_file), "..", ".."), mustWork = TRUE)

task_id <- as_int_arg("--task-id", Sys.getenv("LOCAL_TASK_ID", unset = "1"))
task_offset <- as_int_arg("--offset", Sys.getenv("LOCAL_TASK_OFFSET", unset = "0"))
cpus <- as_int_arg("--cpus", Sys.getenv("LOCAL_CPUS", unset = "4"))

default_output_dir <- file.path(repo_root, "Bunya", "ENZAMET", "output_local")
output_dir <- arg_value("--output-dir", Sys.getenv("SIM_OUTPUT_DIR", unset = default_output_dir))
output_dir <- normalizePath(output_dir, mustWork = FALSE)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

default_cmdstan_path <- Sys.getenv("CMDSTAN", unset = "")
if (!nzchar(default_cmdstan_path)) {
  cmdstan_candidates <- list.dirs(path.expand("~/.cmdstan"), recursive = FALSE, full.names = TRUE)
  cmdstan_candidates <- cmdstan_candidates[grepl("cmdstan-[0-9]", basename(cmdstan_candidates))]
  default_cmdstan_path <- if (length(cmdstan_candidates) > 0) {
    sort(cmdstan_candidates, decreasing = TRUE)[[1]]
  } else {
    "~/.cmdstan/cmdstan-2.37.0"
  }
}
cmdstan_path <- arg_value("--cmdstan-path", default_cmdstan_path)
cmdstan_path <- path.expand(cmdstan_path)

Sys.setenv(
  SLURM_ARRAY_TASK_ID = as.character(task_id),
  SLURM_ARRAY_TASK_OFFSET = as.character(task_offset),
  SLURM_CPUS_PER_TASK = as.character(cpus),
  SIM_OUTPUT_DIR = output_dir,
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

pipeline_file <- file.path(repo_root, "Bunya", "ENZAMET", "ENZAMET-DM.R")
pipeline <- readLines(pipeline_file, warn = FALSE)

# Make the Bunya-oriented script portable to the current checkout.
pipeline <- gsub("~/JoMoNoPH-share", repo_root, pipeline, fixed = TRUE)
pipeline <- gsub("~/.cmdstan/cmdstan-2.37.0", cmdstan_path, pipeline, fixed = TRUE)

# Local timing instrumentation. The Bunya script is kept unchanged; these edits
# are applied only to the in-memory local copy of the pipeline.
rm_line <- grep("rm(list = ls())", pipeline, fixed = TRUE)
if (length(rm_line) == 1) {
  pipeline <- append(
    pipeline,
    c(
      ".local_batch_start <- Sys.time()",
      "message(\"Local batch started: \", format(.local_batch_start, \"%Y-%m-%d %H:%M:%S\"))"
    ),
    after = rm_line
  )
}

rep_elapsed_line <- grep('cat(sprintf("Trial %d elapsed: %.1f min\\n",', pipeline, fixed = TRUE)
if (length(rep_elapsed_line) == 1) {
  pipeline[rep_elapsed_line] <- "    .local_rep_elapsed_secs <- as.numeric(difftime(Sys.time(), t0, units = \"secs\"))"
  pipeline[rep_elapsed_line + 1] <- "    message(sprintf(\"Replicate/trial %d finished; elapsed: %.1f min (%.1f sec)\","
  pipeline <- append(
    pipeline,
    "                    .x, .local_rep_elapsed_secs / 60, .local_rep_elapsed_secs))",
    after = rep_elapsed_line + 1
  )
}

save_line <- grep("saveRDS(results, file = out_file)", pipeline, fixed = TRUE)
if (length(save_line) == 1) {
  pipeline <- append(
    pipeline,
    c(
      ".local_batch_elapsed_secs <- as.numeric(difftime(Sys.time(), .local_batch_start, units = \"secs\"))",
      "message(sprintf(\"Local batch total elapsed: %.1f min (%.1f sec)\",",
      "                .local_batch_elapsed_secs / 60, .local_batch_elapsed_secs))"
    ),
    after = save_line
  )
}

# Optional local-only overrides. Leave unset to match ENZAMET-DM.R exactly.
batch_size <- arg_value("--batch-size", NULL)
if (!is.null(batch_size)) {
  batch_size <- as_int_arg("--batch-size", NA_integer_)
  pipeline <- sub("batch_size <- 8", sprintf("batch_size <- %d", batch_size), pipeline, fixed = TRUE)
}

sims_per_config <- arg_value("--sims-per-config", NULL)
if (!is.null(sims_per_config)) {
  sims_per_config <- as_int_arg("--sims-per-config", NA_integer_)
  pipeline <- sub("sims_per_config <- 1000", sprintf("sims_per_config <- %d", sims_per_config), pipeline, fixed = TRUE)
}

safe_elapsed <- arg_value("--safe-elapsed", NULL)
if (!is.null(safe_elapsed)) {
  safe_elapsed <- as_int_arg("--safe-elapsed", NA_integer_)
  pipeline <- sub("safe_elapsed <- 14400", sprintf("safe_elapsed <- %d", safe_elapsed), pipeline, fixed = TRUE)
}

seed_offset <- arg_value("--seed-offset", NULL)
if (!is.null(seed_offset)) {
  seed_offset <- as_int_arg("--seed-offset", NA_integer_)
  pipeline <- sub(
    "trial_ids <- ((batch_id - 1) * batch_size + 1):(batch_id * batch_size)",
    sprintf("trial_ids <- (((batch_id - 1) * batch_size + 1):(batch_id * batch_size)) + %d", seed_offset),
    pipeline,
    fixed = TRUE
  )
}

fixed_seed <- arg_value("--fixed-seed", NULL)
if (!is.null(fixed_seed)) {
  fixed_seed <- as_int_arg("--fixed-seed", NA_integer_)
  pipeline <- sub(
    "unique_seed <- 1e6 * config_hash + i",
    sprintf("unique_seed <- %d", fixed_seed),
    pipeline,
    fixed = TRUE
  )
}

if (has_flag("--force-recompile")) {
  pipeline <- gsub("force_recompile = FALSE", "force_recompile = TRUE", pipeline, fixed = TRUE)
}

stan_model_file <- arg_value("--stan-model-file", NULL)
if (!is.null(stan_model_file)) {
  stan_model_file <- path.expand(stan_model_file)
  if (!grepl("^/", stan_model_file)) {
    stan_model_file <- file.path(repo_root, stan_model_file)
  }
  stan_model_file <- normalizePath(stan_model_file, mustWork = TRUE)
  pipeline <- sub(
    'stan_model_file_joint <- ".*Stan/Bernstein-Polynomials-JM-Hist\\.stan"',
    sprintf('stan_model_file_joint <- "%s"', stan_model_file),
    pipeline
  )
}

message("Local ENZAMET-DM launcher")
message("  repo_root: ", repo_root)
message("  task_id: ", task_id)
message("  task_offset: ", task_offset)
message("  cpus: ", cpus)
message("  output_dir: ", output_dir)
message("  cmdstan_path: ", cmdstan_path)
if (!is.null(seed_offset)) {
  message("  seed_offset: ", seed_offset)
}
if (!is.null(fixed_seed)) {
  message("  fixed_seed: ", fixed_seed)
}
if (!is.null(stan_model_file)) {
  message("  stan_model_file: ", stan_model_file)
}
message("  force_recompile: ", has_flag("--force-recompile"))

eval(parse(text = pipeline), envir = globalenv())
