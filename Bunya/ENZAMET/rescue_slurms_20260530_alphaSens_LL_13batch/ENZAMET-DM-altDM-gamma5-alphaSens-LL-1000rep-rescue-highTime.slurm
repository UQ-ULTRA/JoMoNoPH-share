#!/bin/bash -l
# High-time rescue for the slowest remaining LL120_230 alpha-prior sensitivity
# batches. These had incomplete chains around 900-1400/2000 under the original
# time cap, so they get 5h50m per replicate and a 12h30m walltime.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4
#SBATCH --mem=8GB
#SBATCH --time=12:30:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-AS-RH
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Original global sim_ids under:
#   SIM_BASELINE_FILTER=LL120_230, SIMS_PER_CONFIG=1000, SIM_BATCH_SIZE=4.
GLOBAL_SIM_IDS=(
  1262
  1315
  1047
  1482
)

ALPHA_LABELS=(
  maxF125
  maxF125
  maxF150
  maxF300
)

ALPHA_MAX_F=(
  1.25
  1.25
  1.5
  3
)

RESCUE_LABELS=(
  "maxF125_cfg06_scen3_admin_batch012_slow"
  "maxF125_cfg06_scen3_admin_batch065_slow_two_nulls"
  "maxF150_cfg05_scen3_cens50_batch047_slow_two_nulls"
  "maxF300_cfg06_scen3_admin_batch232_slow_two_nulls"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
GLOBAL_SIM_ID="${GLOBAL_SIM_IDS[$IDX]}"
RESCUE_LABEL="${RESCUE_LABELS[$IDX]}"

export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))
export ENZ_ACCOUNT_LABEL=RESCUE
export ENZ_ALPHA_LABEL="${ALPHA_LABELS[$IDX]}"
export SIM_ALPHA_MAX_F="${ALPHA_MAX_F[$IDX]}"

# High-time rescue: 5h50m per replicate. With 8 cores, 4 replicates run in two
# waves of two, so the 12h30m walltime leaves roughly 50 minutes for overhead/save.
export SIM_SAFE_ELAPSED_OVERRIDE=21000

echo "Alpha-sensitivity high-time rescue batch"
echo "Target global sim_id: $GLOBAL_SIM_ID"
echo "Rescue label: $RESCUE_LABEL"
echo "Dynamic SLURM_ARRAY_TASK_OFFSET: $SLURM_ARRAY_TASK_OFFSET"
echo "ENZ_ALPHA_LABEL: $ENZ_ALPHA_LABEL"
echo "SIM_ALPHA_MAX_F: $SIM_ALPHA_MAX_F"
echo "SIM_SAFE_ELAPSED_OVERRIDE: $SIM_SAFE_ELAPSED_OVERRIDE"

source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-alphaSens-LL-1000rep-common.sh"
