#!/bin/bash -l
# Standard-time rescue for remaining LL120_230 alpha-prior sensitivity batches.
# These targets either failed before Stan started, were close to completion under
# the original time cap, or had all chains finish but still saved a NULL entry.
# Each array task reruns one original 4-replicate batch and overwrites the
# missing/incomplete .rds file in the original alpha-specific output folder.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-9
#SBATCH --mem=8GB
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-AS-RS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Original global sim_ids under:
#   SIM_BASELINE_FILTER=LL120_230, SIMS_PER_CONFIG=1000, SIM_BATCH_SIZE=4.
GLOBAL_SIM_IDS=(
  522
  1261
  1264
  1398
  742
  861
  902
  1110
  1123
)

ALPHA_LABELS=(
  maxF125
  maxF125
  maxF125
  maxF125
  maxF300
  maxF300
  maxF300
  maxF300
  maxF300
)

ALPHA_MAX_F=(
  1.25
  1.25
  1.25
  1.25
  3
  3
  3
  3
  3
)

RESCUE_LABELS=(
  "maxF125_cfg03_scen2_cens50_batch022_missing_serverSocket"
  "maxF125_cfg06_scen3_admin_batch011_near_finish"
  "maxF125_cfg06_scen3_admin_batch014_near_finish"
  "maxF125_cfg06_scen3_admin_batch148_near_finish"
  "maxF300_cfg03_scen2_cens50_batch242_near_finish"
  "maxF300_cfg04_scen2_admin_batch111_near_finish"
  "maxF300_cfg04_scen2_admin_batch152_near_finish"
  "maxF300_cfg05_scen3_cens50_batch110_finished_but_null"
  "maxF300_cfg05_scen3_cens50_batch123_finished_but_null"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
GLOBAL_SIM_ID="${GLOBAL_SIM_IDS[$IDX]}"
RESCUE_LABEL="${RESCUE_LABELS[$IDX]}"

export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))
export ENZ_ACCOUNT_LABEL=RESCUE
export ENZ_ALPHA_LABEL="${ALPHA_LABELS[$IDX]}"
export SIM_ALPHA_MAX_F="${ALPHA_MAX_F[$IDX]}"

# Standard rescue: 4h50m per replicate. With 8 cores, 4 replicates run in two
# waves of two, so the 10h walltime leaves roughly 20 minutes for overhead/save.
export SIM_SAFE_ELAPSED_OVERRIDE=17400

echo "Alpha-sensitivity standard rescue batch"
echo "Target global sim_id: $GLOBAL_SIM_ID"
echo "Rescue label: $RESCUE_LABEL"
echo "Dynamic SLURM_ARRAY_TASK_OFFSET: $SLURM_ARRAY_TASK_OFFSET"
echo "ENZ_ALPHA_LABEL: $ENZ_ALPHA_LABEL"
echo "SIM_ALPHA_MAX_F: $SIM_ALPHA_MAX_F"
echo "SIM_SAFE_ELAPSED_OVERRIDE: $SIM_SAFE_ELAPSED_OVERRIDE"

source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-alphaSens-LL-1000rep-common.sh"
