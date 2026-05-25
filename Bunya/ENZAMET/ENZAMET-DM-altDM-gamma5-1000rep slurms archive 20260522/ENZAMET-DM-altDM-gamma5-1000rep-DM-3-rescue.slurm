#!/bin/bash -l
# DM part 3 rescue: re-run global array task 5940.
# This corresponds to config 24, batch 190:
#   BP_k5_n1100_censadmin_W130_380_scen2_FU120_batch190.rds
# Keep the array index below Bunya's 0-999 array-index limit.
# The shared launcher writes to the common output folder and will overwrite only
# if the existing batch output is incomplete. This batch previously hit the
# 2h50m R internal limit, so this rescue uses a longer 6h50m per-replicate
# limit and a 10h walltime.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=440
#SBATCH --mem=8GB
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-G5-DM3-R
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=5500
export ENZ_ACCOUNT_LABEL=DM
export ENZ_PART_ID=3-rescue
export SIM_SAFE_ELAPSED_OVERRIDE=24600
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-common.sh"
