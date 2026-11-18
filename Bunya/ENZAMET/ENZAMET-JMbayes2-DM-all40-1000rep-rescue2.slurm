#!/bin/bash -l
# Second rescue for remaining JMbayes2 Weibull-JM failures.
#
# The first rescue fixed the missing batches and most model-level failures.
# These remaining targets all failed quickly in the Weibull-JM fit with:
#   4 nodes produced errors; first error: argument is of length zero
#
# This rescue uses the updated R driver with SIM_JMB2_RETRY_CORES=1, so a
# JMbayes2 fit that fails with cores=4 is retried once with cores=1.
# Output is written to a separate -rescue2 folder and does not overwrite the
# original production or first rescue output folders.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-7
#SBATCH --mem=12GB
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-JMB2-R2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Original global sim_ids under:
#   SIM_BASELINE_FILTER=, SIMS_PER_CONFIG=1000, SIM_BATCH_SIZE=10.
GLOBAL_SIM_IDS=(
  3620
  3633
  3638
  3649
  3804
  3844
  3895
)

RESCUE_LABELS=(
  "cfg37_W210_scen4_cens50_batch020_weibull_retry"
  "cfg37_W210_scen4_cens50_batch033_weibull_retry"
  "cfg37_W210_scen4_cens50_batch038_weibull_retry"
  "cfg37_W210_scen4_cens50_batch049_weibull_retry"
  "cfg39_W210_scen5_cens50_batch004_weibull_retry"
  "cfg39_W210_scen5_cens50_batch044_weibull_retry"
  "cfg39_W210_scen5_cens50_batch095_weibull_retry"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
GLOBAL_SIM_ID="${GLOBAL_SIM_IDS[$IDX]}"
RESCUE_LABEL="${RESCUE_LABELS[$IDX]}"

export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))
export ENZ_JMB2_SLICE_LABEL=RESCUE2

export JOB_OUTPUT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
JOB_OUTPUT_DATE_FILE="$HOME/JoMoNoPH-share/Bunya/ENZAMET/.out_date_${JOB_OUTPUT_ID}"
JOB_OUTPUT_DATE_LOCK="${JOB_OUTPUT_DATE_FILE}.lock"
if [ ! -s "$JOB_OUTPUT_DATE_FILE" ]; then
  if mkdir "$JOB_OUTPUT_DATE_LOCK" 2>/dev/null; then
    date +%Y%m%d > "$JOB_OUTPUT_DATE_FILE"
    rmdir "$JOB_OUTPUT_DATE_LOCK"
  else
    while [ ! -s "$JOB_OUTPUT_DATE_FILE" ]; do sleep 1; done
  fi
fi

export JOB_OUTPUT_DATE="$(cat "$JOB_OUTPUT_DATE_FILE")"
export JOB_OUTPUT_DIR="$HOME/out_${JOB_OUTPUT_DATE}_${JOB_OUTPUT_ID}"
export JOB_OUTPUT_TASK_SUFFIX="${SLURM_ARRAY_TASK_ID:+_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$JOB_OUTPUT_DIR"
exec > "$JOB_OUTPUT_DIR/${SLURM_JOB_NAME}_${JOB_OUTPUT_ID}${JOB_OUTPUT_TASK_SUFFIX}.out" 2>&1

module purge
module load r/4.4.0-gfbf-2023a

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

export SIM_OUTPUT_DATE="$JOB_OUTPUT_DATE"
export SIM_OUTPUT_JOB_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"

# Match the production grid exactly.
export SIM_BASELINE_FILTER=
export SIM_BATCH_SIZE=10
export SIMS_PER_CONFIG=1000

# More generous than production; the failures were fast, but the retry may run
# four JMbayes2 chains sequentially for a failed model.
export SIM_SAFE_ELAPSED=1800

export SIM_JMB2_MODELS=all
export SIM_JMB2_SAVE_MODE=summary
export SIM_JMB2_N_ITER=3500
export SIM_JMB2_N_BURNIN=500
export SIM_JMB2_N_THIN=1
export SIM_JMB2_N_CHAINS=4
export SIM_JMB2_CORES=4
export SIM_JMB2_RETRY_CORES=1

export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260601_ENZAMET-JMbayes2-DM-all40-1000rep-rescue2"
mkdir -p "$SIM_OUTPUT_DIR"

echo "JMbayes2 all40 production rescue2"
echo "Target global sim_id: $GLOBAL_SIM_ID"
echo "Rescue label: $RESCUE_LABEL"
echo "Dynamic SLURM_ARRAY_TASK_OFFSET: $SLURM_ARRAY_TASK_OFFSET"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Array Task Offset: $SLURM_ARRAY_TASK_OFFSET"
echo "JMbayes2 slice label: $ENZ_JMB2_SLICE_LABEL"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: JMbayes2 ENZAMET comparison; rescue2; all 40 configs; 1000 reps/config"
echo "SIM_BASELINE_FILTER: $SIM_BASELINE_FILTER"
echo "SIM_BATCH_SIZE: $SIM_BATCH_SIZE"
echo "SIMS_PER_CONFIG: $SIMS_PER_CONFIG"
echo "SIM_SAFE_ELAPSED: $SIM_SAFE_ELAPSED"
echo "SIM_JMB2_MODELS: $SIM_JMB2_MODELS"
echo "SIM_JMB2_SAVE_MODE: $SIM_JMB2_SAVE_MODE"
echo "SIM_JMB2_N_ITER: $SIM_JMB2_N_ITER"
echo "SIM_JMB2_N_BURNIN: $SIM_JMB2_N_BURNIN"
echo "SIM_JMB2_N_CHAINS: $SIM_JMB2_N_CHAINS"
echo "SIM_JMB2_CORES: $SIM_JMB2_CORES"
echo "SIM_JMB2_RETRY_CORES: $SIM_JMB2_RETRY_CORES"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript "$HOME/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R"
