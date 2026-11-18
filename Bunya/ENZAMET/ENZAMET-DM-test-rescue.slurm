#!/bin/bash -l
# Rescue the 40-config/10-replicate ENZAMET test batches from job 24604129
# where one or more trials hit the 300-minute R internal elapsed-time limit.
# ENZAMET-DM.R skips complete .rds files and overwrites incomplete/NULL ones.
# Array IDs are the original global task IDs because this test run uses no offset.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=20
#SBATCH --array=7,10,12,20,24,25,26,27,31,39,41,48,51,52,59,64,68,71,72,79
#SBATCH --mem=16GB
#SBATCH --time=09:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
##SBATCH --reservation=ding_diag
#SBATCH --job-name=ENZAMET-TSTR
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

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
export SLURM_ARRAY_TASK_OFFSET=0
export SIM_OUTPUT_DATE="$JOB_OUTPUT_DATE"
export SIM_OUTPUT_JOB_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
export SIM_BATCH_SIZE=5
export SIMS_PER_CONFIG=10
export SIM_SAFE_ELAPSED=32400
export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_test_10rep_W210_centered_td12"
mkdir -p "$SIM_OUTPUT_DIR"

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Array Task Offset: $SLURM_ARRAY_TASK_OFFSET"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: BP sAFT-JM test rescue run for job 24604129 internal-timeout batches"
echo "SIM_BATCH_SIZE: $SIM_BATCH_SIZE"
echo "SIMS_PER_CONFIG: $SIMS_PER_CONFIG"
echo "SIM_SAFE_ELAPSED: $SIM_SAFE_ELAPSED"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM.R
