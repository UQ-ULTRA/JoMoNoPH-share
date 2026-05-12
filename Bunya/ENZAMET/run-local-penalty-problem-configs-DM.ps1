$ErrorActionPreference = "Stop"

# Windows PowerShell version for replicating the previous penalty-model local run.
# Run from the repository root.
#
# This intentionally does NOT pass adapt_delta or max_treedepth, so CmdStan uses
# its defaults. This is for reproducing the earlier penalty-model diagnostic run.
#
# Config/task mapping with batch_size = 1 and sims_per_config = 1:
#   task 1  = LL120_230, scenario 1, 50% censoring
#   task 16 = W090_380,  scenario 3, administrative censoring
#   task 20 = W090_380,  scenario 5, administrative censoring

$CPUS = if ($env:CPUS) { $env:CPUS } else { "8" }
$FIXED_SEED = if ($env:FIXED_SEED) { $env:FIXED_SEED } else { "5" }

$REPO_ROOT = (Get-Location).Path
$OUTPUT_DIR = if ($env:OUTPUT_DIR) {
    $env:OUTPUT_DIR
} else {
    Join-Path $REPO_ROOT "Bunya/ENZAMET/output_local_penalty_fixed_seed$FIXED_SEED"
}

$STAN_FILE = Join-Path $REPO_ROOT "Stan/Bernstein-Polynomials-JM-Hist-Penalty.stan"
$LOCAL_RUNNER = Join-Path $REPO_ROOT "Bunya/ENZAMET/ENZAMET-local-DM.R"

New-Item -ItemType Directory -Force -Path $OUTPUT_DIR | Out-Null

Write-Host "Output directory: $OUTPUT_DIR"
Write-Host "Stan model: $STAN_FILE"
Write-Host "Fixed seed: $FIXED_SEED"
Write-Host "CPUs per run: $CPUS"
Write-Host "Sampling controls: CmdStan defaults"

$commonArgs = @(
    "--batch-size", "1",
    "--sims-per-config", "1",
    "--fixed-seed", $FIXED_SEED,
    "--cpus", $CPUS,
    "--safe-elapsed", "21600",
    "--output-dir", $OUTPUT_DIR,
    "--stan-model-file", $STAN_FILE
)

Write-Host "Running task 1 with forced compilation..."
Rscript $LOCAL_RUNNER --task-id 1 @commonArgs --force-recompile

Write-Host "Running task 16..."
Rscript $LOCAL_RUNNER --task-id 16 @commonArgs

Write-Host "Running task 20..."
Rscript $LOCAL_RUNNER --task-id 20 @commonArgs

Write-Host "Finished penalty-model reruns."
