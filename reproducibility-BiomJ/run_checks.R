# START HERE: open reproducibility-BiomJ.Rproj, then open this script.
# Run sections 1-4 in order, or click Source to run the selected task.
# The default task checks files/software; it does NOT compile or fit a model.
# No packages are installed automatically. See README.md for installation.

# 1. Choose what to do ------------------------------------------------------
# Tasks: "check", "plan", "prepare", "fit", "records", "case-study",
#        "results", "aggregate", "validate".
# See README.md for each task's inputs, runtime and output locations.
task <- "check"

# Collections: "main", "alpha125", "alpha150", "alpha300", "comparator".
# "alpha" combines the three sensitivity collections; "all" combines all five.
collection <- "main"

# 2. Choose simulation settings --------------------------------------------
# Original IDs: main/comparator 1:40; alpha sensitivity 1:10 (log-logistic).
# Use c(1, 2) for multiple configurations. Configuration 1 is an example,
# not a prescribed selection. Inspect the candidate table in section 4.
# NULL lists all candidates in a plan only; "all" explicitly selects every
# available configuration (not recommended for a small audit fit).
config_ids <- 1L

# Five datasets PER selected configuration. Increase to 1:10 if desired,
# or use c(7, 42) to check specific original replicates. Maximum ID: 1000.
replicates <- 1:5

# All four chains are retained. This controls simultaneous chains only.
# Datasets are processed sequentially; no SLURM or terminal is needed.
parallel_chains <- 4L
resume <- TRUE        # Skip compatible completed local checks.
cmdstan <- ""         # Use configured CmdStan; or supply its installation path.
# This Mac's cached CmdStan header predates an Apple Clang update. Bypass it
# for compilation; the model and sampling settings are unchanged. Other
# machines may also use TRUE, or use FALSE with compatible cached headers.
no_pch <- TRUE
summaries <- "saved"  # For task="results"; "regenerated" after aggregation.

# 3. Load the existing implementation --------------------------------------
# With the supplied R project open, leave package_dir unchanged. The parent
# JoMoNoPH-share project also works. Else set an explicit package folder here.
package_dir <- "."
if (!file.exists(file.path(package_dir, "R", "audit.R")) &&
    file.exists(file.path(package_dir, "reproducibility-BiomJ", "R", "audit.R"))) {
  package_dir <- file.path(package_dir, "reproducibility-BiomJ")
}
if (!file.exists(file.path(package_dir, "R", "audit.R"))) {
  stop("Open reproducibility-BiomJ.Rproj, or set package_dir to the extracted package folder.")
}
root <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
for (helper in c("audit.R", "simulation.R", "export-local.R", "session.R", "entry.R")) {
  source(file.path(root, "R", helper), local = TRUE)
}

# 4. Inspect the available settings, then run the selected task -------------
# View(configurations) in RStudio to see the full grid. Changing collection,
# config_ids or replicates selects datasets; it does not alter the model.
configurations <- audit_grid(root)
print(configurations[, c("full_config_id", "baseline_label", "scenario",
                          "cens_label", "max_FU")], row.names = FALSE)
message("Selected task: ", task, "; collection: ", collection)
audit_result <- audit_run(root, task, collection, config_ids, replicates,
                          parallel_chains, cmdstan, no_pch, resume, summaries)

# After inspecting task="plan", change task to "prepare" (no MCMC) or "fit",
# then rerun sections 1-4. Five sAFT-JM fits may take several hours.
# Use task="records" after fitting to combine estimates, diagnostics, runtimes
# and comparisons with the matching saved production replicates into CSV/RDS.
# All generated files stay under this package's output/; supplied references
# and the original production files are not overwritten.
