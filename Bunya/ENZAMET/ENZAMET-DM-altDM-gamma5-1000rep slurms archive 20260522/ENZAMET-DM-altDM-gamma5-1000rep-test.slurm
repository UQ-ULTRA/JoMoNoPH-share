#!/bin/bash -l
# Smoke test for the 40-config x 1000-replicate alt-DM manuscript run.
# This runs only two array tasks from the AM-1 slice:
#   global task 1 -> config 1, batch 1, reps 1-4
#   global task 2 -> config 1, batch 2, reps 5-8
# It uses the same shared launcher and resource settings as the production files.

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-2
#SBATCH --mem=8GB
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-G5-TEST
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=0
export ENZ_ACCOUNT_LABEL=TEST
export ENZ_PART_ID=test
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-common.sh"
