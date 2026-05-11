#!/bin/bash -l
# Min P2 rescue run.
# Reruns the latest P2 walltime failures from job 24444774 plus batches with
# the old 170-minute R internal elapsed-time-limit signature.
# ENZAMET-DM.R skips complete .rds files and overwrites incomplete/NULL ones.
# Local array IDs use offset 938, mapping to global tasks 939-1876.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=16
#SBATCH --array=5,7,12,17,28,32,41,42,45,47,131,196,202,221,222,223,224,232,240,243,262,271,275,281,290,296,298,299,306,332,360,406,433,443,445,447,448,449,453,454,455,456,458,460,461,465,470,478,485,494,495,496,498,499,501,503,504,505,506,507,513,522,524,529,544,552,553,554,556,562,586,589,590,595,601,604,605,607,609,610,613,620,621,640,659,672,673,675,683,684,689,691,692,695,697,699,700,702,707,713,715,717,725,726,734,735,736,738,740,741,743,746,753,757,758,759,761,762,765,771,774,775,776,777,778,779,790,792,793,794,798,799,803,811,818,822,825,826,834,838,841,853,855,869,875,885,887,904,906,907,915,924,926,935,937,938
#SBATCH --mem=10GB
#SBATCH --time=09:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
##SBATCH --reservation=ding_diag
#SBATCH --job-name=ENZAMET-P2R
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
export SLURM_ARRAY_TASK_OFFSET=938
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
echo "Model: BP sAFT-JM rescue run; Min part 2/4"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM.R
