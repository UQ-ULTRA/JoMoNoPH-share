#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --job-name=CompStan
#SBATCH --account=a_uqccr
#SBATCH --time=00:10:00
#SBATCH --partition=general
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export JOB_OUTPUT_DATE="$(date +%Y%m%d)"
export JOB_OUTPUT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
export JOB_OUTPUT_DIR="$HOME/out_${JOB_OUTPUT_DATE}_${JOB_OUTPUT_ID}"
export JOB_OUTPUT_TASK_SUFFIX="${SLURM_ARRAY_TASK_ID:+_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$JOB_OUTPUT_DIR"
exec > "$JOB_OUTPUT_DIR/${SLURM_JOB_NAME}_${JOB_OUTPUT_ID}${JOB_OUTPUT_TASK_SUFFIX}.out" 2>&1

# --- Load required modules ---
module purge
module load r/4.4.0-gfbf-2023a # load the same version in the ENZAMET.slurm

# --- Run R script ---

Rscript compile-stan-models.R
