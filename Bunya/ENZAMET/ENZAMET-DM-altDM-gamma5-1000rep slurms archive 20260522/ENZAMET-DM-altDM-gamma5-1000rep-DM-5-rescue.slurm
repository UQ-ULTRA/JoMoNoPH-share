#!/bin/bash -l
# DM part 5 rescue: re-run global array tasks 9509, 9608, and 9737.
# These correspond to config 39 batches 9, 108, and 237:
#   BP_k5_n1100_cens50_W210_380_scen5_FU120_batch009.rds
#   BP_k5_n1100_cens50_W210_380_scen5_FU120_batch108.rds
#   BP_k5_n1100_cens50_W210_380_scen5_FU120_batch237.rds
# Keep all array indices below Bunya's 0-999 array-index limit.
# The shared launcher writes to the common output folder and will overwrite only
# if the existing batch output is incomplete.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=9,108,237
#SBATCH --mem=8GB
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-G5-DM5-R
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=9500
export ENZ_ACCOUNT_LABEL=DM
export ENZ_PART_ID=5-rescue
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-common.sh"
