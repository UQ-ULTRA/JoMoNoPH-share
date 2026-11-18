#!/bin/bash
# Self-contained launcher body for the JMbayes2 all-40-config 1000-replicate run.
#
# Design:
#   - all 40 ENZAMET configs from Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R
#   - 1000 replicates/config, 10 replicates/batch = 100 batches/config
#   - 40 configs x 100 batches/config = 4000 global array tasks
#   - 5 Slurm wrappers split the 4000 tasks as 800 local array tasks each
#   - each task requests 8 CPUs so two replicates can run in parallel, with
#     4 JMbayes2 cores per replicate

if [ "${SLURM_ARRAY_TASK_OFFSET+x}" != "x" ]; then
  echo "SLURM_ARRAY_TASK_OFFSET must be set by the production .slurm wrapper."
  exit 1
fi

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

# Empty baseline filter means use the full 40-config grid from the R driver.
export SIM_BASELINE_FILTER=
export SIM_BATCH_SIZE=10
export SIMS_PER_CONFIG=1000

# Timeout is per replicate. With two parallel replicates/wave, 10 replicates
# run in five waves; 480 seconds gives up to 40 minutes of fitting time under
# the 45-minute Slurm limit.
export SIM_SAFE_ELAPSED=480

export SIM_JMB2_MODELS=all
export SIM_JMB2_SAVE_MODE=summary

# Keep JMbayes2's default iteration/burn-in/thinning settings for posterior
# summaries so the comparator is not disadvantaged by a shorter MCMC schedule.
export SIM_JMB2_N_ITER=3500
export SIM_JMB2_N_BURNIN=500
export SIM_JMB2_N_THIN=1

# JMbayes2 defaults to 3 chains; use 4 chains here to align diagnostics with
# the Stan sAFT-JM runs while otherwise retaining JMbayes2's MCMC defaults.
# The fitted object is still discarded in summary mode; only compact posterior
# summary rows, including derived sd_1_long rows, are saved.
export SIM_JMB2_N_CHAINS=4
export SIM_JMB2_CORES=4

export SIM_OUTPUT_DIR="${SIM_OUTPUT_DIR_OVERRIDE:-$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260602_ENZAMET-JMbayes2-DM-all40-1000rep-randSDsummary}"
mkdir -p "$SIM_OUTPUT_DIR"

if [ -n "${SIM_JMB2_FULL_FIT_DIR_OVERRIDE:-}" ]; then
  export SIM_JMB2_FULL_FIT_DIR="$SIM_JMB2_FULL_FIT_DIR_OVERRIDE"
fi

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Array Task Offset: $SLURM_ARRAY_TASK_OFFSET"
if [ -n "${ENZ_JMB2_SLICE_LABEL:-}" ]; then
  echo "JMbayes2 slice label: $ENZ_JMB2_SLICE_LABEL"
fi
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: JMbayes2 ENZAMET comparison; all 40 configs; 1000 reps/config"
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
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript "$HOME/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R"
