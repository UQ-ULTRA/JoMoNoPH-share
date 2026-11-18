#!/bin/bash -l
# Rescue run for JMbayes2 all-40-config 1000-replicate production batches.
#
# Targets:
#   - 4 missing .rds batches
#   - 25 saved batches with at least one model_status != "ok"
#
# The local array index maps to an explicit original global sim_id. The R driver
# computes global sim_id as SLURM_ARRAY_TASK_ID + SLURM_ARRAY_TASK_OFFSET, so the
# dynamic offset below preserves the original config/batch mapping exactly.
#
# Rescue output is written to a separate -rescue folder and does not overwrite
# the original production output folder.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-29
#SBATCH --mem=12GB
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-JMB2-R
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Original global sim_ids under:
#   SIM_BASELINE_FILTER=, SIMS_PER_CONFIG=1000, SIM_BATCH_SIZE=10.
GLOBAL_SIM_IDS=(
  24
  87
  167
  290
  378
  656
  820
  1000
  1271
  1494
  1510
  1807
  1852
  2043
  2048
  2520
  2644
  3091
  3185
  3227
  3262
  3266
  3620
  3633
  3638
  3649
  3804
  3844
  3895
)

RESCUE_LABELS=(
  "missing_cfg01_batch024_serverSocket"
  "modelFail_cfg01_batch087"
  "modelFail_cfg02_batch067"
  "modelFail_cfg03_batch090"
  "missing_cfg04_batch078_serverSocket"
  "modelFail_cfg07_batch056"
  "modelFail_cfg09_batch020"
  "modelFail_cfg10_batch100"
  "modelFail_cfg13_batch071"
  "modelFail_cfg15_batch094"
  "modelFail_cfg16_batch010"
  "missing_cfg19_batch007_serverSocket"
  "missing_cfg19_batch052_serverSocket"
  "modelFail_cfg21_batch043_two_reps"
  "modelFail_cfg21_batch048"
  "modelFail_cfg26_batch020"
  "modelFail_cfg27_batch044"
  "modelFail_cfg31_batch091"
  "modelFail_cfg32_batch085_timeout"
  "modelFail_cfg33_batch027"
  "modelFail_cfg33_batch062"
  "modelFail_cfg33_batch066"
  "modelFail_cfg37_batch020"
  "modelFail_cfg37_batch033"
  "modelFail_cfg37_batch038"
  "modelFail_cfg37_batch049"
  "modelFail_cfg39_batch004"
  "modelFail_cfg39_batch044"
  "modelFail_cfg39_batch095"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
GLOBAL_SIM_ID="${GLOBAL_SIM_IDS[$IDX]}"
RESCUE_LABEL="${RESCUE_LABELS[$IDX]}"

export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))
export ENZ_JMB2_SLICE_LABEL=RESCUE

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

# Rescue gives more time than production: 15 minutes per replicate. With two
# parallel replicates/wave and 10 replicates/batch, this allows about 75 minutes
# of model-fitting time under the 2-hour Slurm walltime.
export SIM_SAFE_ELAPSED=900

export SIM_JMB2_MODELS=all
export SIM_JMB2_SAVE_MODE=summary
export SIM_JMB2_N_ITER=3500
export SIM_JMB2_N_BURNIN=500
export SIM_JMB2_N_THIN=1
export SIM_JMB2_N_CHAINS=4
export SIM_JMB2_CORES=4

export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260601_ENZAMET-JMbayes2-DM-all40-1000rep-rescue"
mkdir -p "$SIM_OUTPUT_DIR"

echo "JMbayes2 all40 production rescue"
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
echo "Model: JMbayes2 ENZAMET comparison; rescue; all 40 configs; 1000 reps/config"
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
