#!/bin/bash -l
# Rescue for DM part 2, array task 49 only.
# This reruns global array task 3549:
#   Config ID 15, BP_k5_n1100_cens50_W090_380_scen3_FU120, batch 49
# The original task failed before R started because srun could not contact the
# SLURM controller, so this uses the same output directory and normal settings.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=49
#SBATCH --mem=8GB
#SBATCH --time=05:30:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-ZP-DM2R
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=3500
export ENZ_ACCOUNT_LABEL=DM
export ENZ_PART_ID=2-rescue-batch049
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-zeroprior-initinfo-common.sh"
