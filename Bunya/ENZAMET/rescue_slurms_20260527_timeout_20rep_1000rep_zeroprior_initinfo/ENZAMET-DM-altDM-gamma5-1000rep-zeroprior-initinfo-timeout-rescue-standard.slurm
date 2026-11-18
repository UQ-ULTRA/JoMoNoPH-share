#!/bin/bash -l
# Rescue for 18 incomplete full-batch files from the 40-config x 1000-replicate
# zero-prior/init-info production run. Each array task reruns one original
# 4-replicate batch and overwrites the incomplete .rds in the production output
# folder. This keeps SumSimRes-ENZAMET-DM.R compatible because rescued files keep
# the original batch numbering and list length.
#
# The batch with two timed-out replicates is handled separately by the
# high-time rescue slurm in this folder.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-18
#SBATCH --mem=8GB
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-ZP-TO-R
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Original global sim_ids under SIM_BATCH_SIZE=4. These correspond to the
# incomplete batch files containing one timed-out replicate each.
GLOBAL_SIM_IDS=(
  255
  261
  326
  330
  406
  735
  736
  2830
  3803
  3814
  5705
  5829
  6029
  8668
  8694
  8698
  8837
  8849
)

RESCUE_LABELS=(
  "cfg02_batch005_censadmin_LL120_230_scen1_rep17"
  "cfg02_batch011_censadmin_LL120_230_scen1_rep43"
  "cfg02_batch076_censadmin_LL120_230_scen1_rep303"
  "cfg02_batch080_censadmin_LL120_230_scen1_rep317"
  "cfg02_batch156_censadmin_LL120_230_scen1_rep621"
  "cfg03_batch235_cens50_LL120_230_scen2_rep939"
  "cfg03_batch236_cens50_LL120_230_scen2_rep943"
  "cfg12_batch080_censadmin_W090_380_scen1_rep317"
  "cfg16_batch053_censadmin_W090_380_scen3_rep210"
  "cfg16_batch064_censadmin_W090_380_scen3_rep253"
  "cfg23_batch205_cens50_W130_380_scen2_rep819"
  "cfg24_batch079_censadmin_W130_380_scen2_rep313"
  "cfg25_batch029_cens50_W130_380_scen3_rep116"
  "cfg35_batch168_cens50_W210_380_scen3_rep671"
  "cfg35_batch194_cens50_W210_380_scen3_rep773"
  "cfg35_batch198_cens50_W210_380_scen3_rep791"
  "cfg36_batch087_censadmin_W210_380_scen3_rep345"
  "cfg36_batch099_censadmin_W210_380_scen3_rep393"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
GLOBAL_SIM_ID="${GLOBAL_SIM_IDS[$IDX]}"
RESCUE_LABEL="${RESCUE_LABELS[$IDX]}"

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
# SLURM_ARRAY_TASK_OFFSET. Keep the real array index for logging and set the
# offset dynamically to target the original full-batch global sim_id.
export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))

export SIM_BATCH_SIZE=4
export SIMS_PER_CONFIG=1000
export SIM_SAFE_ELAPSED=24600
export SIM_STAN_MODEL_FILE="$HOME/JoMoNoPH-share/Stan/Bernstein-Polynomials-JM-Hist-alt-DM-Gamma5.stan"
export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260525_ENZAMET-altDM-gamma5-1000rep-zeroprior-initinfo"

echo "Timeout rescue: standard single-timeout batches"
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
