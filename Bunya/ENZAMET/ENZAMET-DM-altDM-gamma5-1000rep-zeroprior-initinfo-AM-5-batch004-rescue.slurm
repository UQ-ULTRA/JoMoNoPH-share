#!/bin/bash -l
# Rescue for AM part 5, array task 254 only.
# This reruns global array task 8254:
#   Config ID 34, BP_k5_n1100_censadmin_W210_380_scen2_FU120, batch 4
# The original task failed before sampling because future could not open a
# server socket, so this uses the same output directory and normal settings.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=254
#SBATCH --mem=8GB
#SBATCH --time=05:30:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-ZP-AM5R
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=8000
export ENZ_ACCOUNT_LABEL=AM
export ENZ_PART_ID=5-rescue-batch004
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-zeroprior-initinfo-common.sh"
