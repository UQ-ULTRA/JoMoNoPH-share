#!/bin/bash
# Shared launcher body for the 40-config x 1000-replicate alt-DM manuscript run
# with zero-centered beta priors and compact saved initial-value metadata.
#
# Each account/part .slurm file sets SLURM_ARRAY_TASK_OFFSET before sourcing
# this file. With SIM_BATCH_SIZE=4 and --array=1-500, each submission covers
# 500 array tasks x 4 replicates/task = 2000 fits.

if [ -z "${SLURM_ARRAY_TASK_OFFSET:-}" ]; then
  echo "SLURM_ARRAY_TASK_OFFSET must be set by the account/part .slurm file."
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
export SIM_BATCH_SIZE=4
export SIMS_PER_CONFIG=1000
# Timeout is per replicate. With 8 cores, 4 replicates run in two waves of two;
# 2h55m preserves a saveRDS buffer before the 6h walltime. Rescue jobs can set
# SIM_SAFE_ELAPSED_OVERRIDE before sourcing this file if needed.
if [ -n "${SIM_SAFE_ELAPSED_OVERRIDE:-}" ]; then
  export SIM_SAFE_ELAPSED="$SIM_SAFE_ELAPSED_OVERRIDE"
else
  export SIM_SAFE_ELAPSED=9600
fi
export SIM_STAN_MODEL_FILE="$HOME/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"
export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo"
mkdir -p "$SIM_OUTPUT_DIR"

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Array Task Offset: $SLURM_ARRAY_TASK_OFFSET"
echo "Account label: ${ENZ_ACCOUNT_LABEL:-NA}"
echo "Part ID: ${ENZ_PART_ID:-NA}"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: optimized BP alt-DM sAFT-JM manuscript run; gamma prior SD 5; zero beta priors; saved compact init_info; max_treedepth 10; 40 configs x 1000 reps/config"
echo "SIM_BATCH_SIZE: $SIM_BATCH_SIZE"
echo "SIMS_PER_CONFIG: $SIMS_PER_CONFIG"
echo "SIM_SAFE_ELAPSED: $SIM_SAFE_ELAPSED"
echo "SIM_STAN_MODEL_FILE: $SIM_STAN_MODEL_FILE"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript "$HOME/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM.R"
