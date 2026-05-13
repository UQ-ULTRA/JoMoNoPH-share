#!/usr/bin/env bash
set -euo pipefail

# Local rerun of the three one-replicate JM configs that showed severe
# convergence problems under the standard non-penalised Bernstein model.
#
# This variant isolates the effect of mean-centering the longitudinal response:
#   - Stan model: Bernstein-Polynomials-JM-Hist.stan
#   - NUTS controls: production pipeline settings (max_treedepth = 12; no adapt_delta override)
#   - Response centering: subtract the observed mean of Y_obs within each simulated dataset
#
# Config/task mapping with batch_size = 1 and sims_per_config = 1:
#   task 1  = LL120_230, scenario 1, 50% censoring
#   task 16 = W090_380,  scenario 3, administrative censoring
#   task 20 = W090_380,  scenario 5, administrative censoring

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CPUS="${CPUS:-8}"
FIXED_SEED="${FIXED_SEED:-5}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/Bunya/ENZAMET/output_local_standard_centeredY_default_fixed_seed${FIXED_SEED}}"
STAN_FILE="${REPO_ROOT}/Stan/Bernstein-Polynomials-JM-Hist.stan"
LOCAL_RUNNER="${REPO_ROOT}/Bunya/ENZAMET/ENZAMET-local-DM.R"

mkdir -p "${OUTPUT_DIR}"

common_args=(
  --batch-size 1
  --sims-per-config 1
  --fixed-seed "${FIXED_SEED}"
  --cpus "${CPUS}"
  --safe-elapsed 21600
  --center-longitudinal-response
  --output-dir "${OUTPUT_DIR}"
  --stan-model-file "${STAN_FILE}"
)

if [[ " ${common_args[*]} " != *" --center-longitudinal-response "* ]]; then
  echo "Internal error: centered-Y run is missing --center-longitudinal-response." >&2
  exit 1
fi

echo "Output directory: ${OUTPUT_DIR}"
echo "Stan model: ${STAN_FILE}"
echo "Fixed seed: ${FIXED_SEED}"
echo "CPUs per run: ${CPUS}"
echo "Longitudinal response centering: observed mean of Y_obs within each dataset"
echo "Sampling controls: production pipeline settings, no adapt_delta override, max_treedepth=12"

echo "Running task 1 with forced compilation..."
Rscript "${LOCAL_RUNNER}" --task-id 1 "${common_args[@]}" --force-recompile

echo "Running task 16..."
Rscript "${LOCAL_RUNNER}" --task-id 16 "${common_args[@]}"

echo "Running task 20..."
Rscript "${LOCAL_RUNNER}" --task-id 20 "${common_args[@]}"

echo "Finished standard centered-Y reruns."
