#!/bin/bash -l
# Rerun the latest part-4 walltime failures from job 24438528 plus batches
# with the old 170-minute R internal elapsed-time-limit signature.
# ENZAMET-DM.R skips complete .rds files and overwrites incomplete/NULL ones.
# Local array IDs use offset 2813, mapping to global tasks 2814-3750.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=16
#SBATCH --array=4,9,12,14,16,24,25,26,27,31,36,39,41,42,44,45,47,50,52,53,56,58,59,60,61,62,63,64,65,66,68,69,70,72,73,74,75,77,78,81,82,89,91,92,93,95,97,99,100,101,103,105,107,108,109,112,114,116,122,124,126,140,141,142,145,146,149,150,151,157,159,161,162,165,166,167,172,173,174,175,179,182,184,188,191,197,199,200,205,207,211,212,222,225,229,231,233,234,235,237,238,239,243,246,251,253,255,259,261,263,265,266,268,273,276,277,279,280,283,284,286,287,290,291,292,293,296,299,301,303,305,306,309,314,316,317,322,330,335,337,338,340,345,348,350,352,355,357,361,363,364,368,369,373,376,377,382,383,384,385,386,387,393,401,406,407,410,413,414,415,417,419,420,423,424,427,431,432,433,434,435,437,614,635,817,828,861,867,869,880,882,883
#SBATCH --mem=10GB
#SBATCH --time=09:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
##SBATCH --reservation=ding_diag
#SBATCH --job-name=ENZAMET-P4R
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
export SLURM_ARRAY_TASK_OFFSET=2813
export SIM_OUTPUT_DATE="$JOB_OUTPUT_DATE"
export SIM_OUTPUT_JOB_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output"
mkdir -p "$SIM_OUTPUT_DIR"

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Array Task Offset: $SLURM_ARRAY_TASK_OFFSET"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: BP sAFT-JM rescue run; part 4/4 missing walltime batches"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM.R
