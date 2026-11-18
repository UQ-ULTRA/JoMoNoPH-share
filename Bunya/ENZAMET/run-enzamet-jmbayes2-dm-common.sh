#!/bin/bash
# Shared launcher body for ENZAMET JMbayes2 comparison runs.
#
# The R driver contains the full 40-config data-generation grid. This launcher
# defaults to the 10 log-logistic configs by setting SIM_BASELINE_FILTER=LL120_230.
# Override variables in the calling .slurm file or with exported *_OVERRIDE values.

if [ -z "${SLURM_ARRAY_TASK_OFFSET:-}" ]; then
  export SLURM_ARRAY_TASK_OFFSET=0
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

if [ "${SIM_BASELINE_FILTER_OVERRIDE+x}" = "x" ]; then
  export SIM_BASELINE_FILTER="$SIM_BASELINE_FILTER_OVERRIDE"
else
  export SIM_BASELINE_FILTER=LL120_230
fi
export SIM_BATCH_SIZE="${SIM_BATCH_SIZE_OVERRIDE:-1}"
export SIMS_PER_CONFIG="${SIMS_PER_CONFIG_OVERRIDE:-2}"
export SIM_SAFE_ELAPSED="${SIM_SAFE_ELAPSED_OVERRIDE:-21600}"

export SIM_JMB2_MODELS="${SIM_JMB2_MODELS_OVERRIDE:-all}"
export SIM_JMB2_SAVE_MODE="${SIM_JMB2_SAVE_MODE_OVERRIDE:-summary}"
export SIM_JMB2_N_ITER="${SIM_JMB2_N_ITER_OVERRIDE:-3500}"
export SIM_JMB2_N_BURNIN="${SIM_JMB2_N_BURNIN_OVERRIDE:-500}"
export SIM_JMB2_N_THIN="${SIM_JMB2_N_THIN_OVERRIDE:-1}"
export SIM_JMB2_N_CHAINS="${SIM_JMB2_N_CHAINS_OVERRIDE:-3}"
export SIM_JMB2_CORES="${SIM_JMB2_CORES_OVERRIDE:-3}"

if [ -n "${SIM_OUTPUT_DIR_OVERRIDE:-}" ]; then
  export SIM_OUTPUT_DIR="$SIM_OUTPUT_DIR_OVERRIDE"
else
  if [ -z "$SIM_BASELINE_FILTER" ]; then
    SIM_OUTPUT_LABEL="${SIM_OUTPUT_LABEL_OVERRIDE:-ENZAMET-JMbayes2-DM-all40-pilot}"
  else
    SIM_OUTPUT_LABEL="${SIM_OUTPUT_LABEL_OVERRIDE:-ENZAMET-JMbayes2-DM-LL-pilot}"
  fi
  export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_${JOB_OUTPUT_DATE}_${SIM_OUTPUT_LABEL}"
fi
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
echo "Model: JMbayes2 ENZAMET comparison"
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
