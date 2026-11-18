#!/bin/bash -l
# High-time rescue for the one incomplete production batch with two timed-out
# replicates:
#   cfg35 batch171, BP_k5_n1100_cens50_W210_380_scen3_FU120,
#   reps 681 and 683.
#
# This reruns the original 4-replicate batch and overwrites the incomplete .rds
# in the production output folder, preserving the original batch numbering.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-1
#SBATCH --mem=8GB
#SBATCH --time=15:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-ZP-TO-W
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

GLOBAL_SIM_ID=8671
RESCUE_LABEL="cfg35_batch171_cens50_W210_380_scen3_reps681_683"

JOB_OUTPUT_DATE="$(date +%Y%m%d)"
JOB_OUTPUT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
JOB_OUTPUT_DIR="$HOME/out_${JOB_OUTPUT_DATE}_${JOB_OUTPUT_ID}_timeout_rescue"
mkdir -p "$JOB_OUTPUT_DIR"
exec > "$JOB_OUTPUT_DIR/${SLURM_JOB_NAME}_${JOB_OUTPUT_ID}_${SLURM_ARRAY_TASK_ID}.out" 2>&1

module purge
module load r/4.4.0-gfbf-2023a

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ENZAMET-DM.R computes global sim_id as SLURM_ARRAY_TASK_ID +
# SLURM_ARRAY_TASK_OFFSET. This maps array task 1 to original global sim_id 8671.
export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))

export SIM_BATCH_SIZE=4
export SIMS_PER_CONFIG=1000
export SIM_SAFE_ELAPSED=24600
export SIM_STAN_MODEL_FILE="$HOME/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"
export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo"

echo "Timeout rescue: high-time W210 scen3 batch171"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Target global sim_id: $GLOBAL_SIM_ID"
echo "Rescue label: $RESCUE_LABEL"
echo "Dynamic SLURM_ARRAY_TASK_OFFSET: $SLURM_ARRAY_TASK_OFFSET"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "SIM_BATCH_SIZE: $SIM_BATCH_SIZE"
echo "SIMS_PER_CONFIG: $SIMS_PER_CONFIG"
echo "SIM_SAFE_ELAPSED: $SIM_SAFE_ELAPSED"
echo "SIM_STAN_MODEL_FILE: $SIM_STAN_MODEL_FILE"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript "$HOME/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-DM.R"
