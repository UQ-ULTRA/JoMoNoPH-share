#!/bin/bash -l
# DM part 5: global array tasks 9501-10000, configs 39-40.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-500
#SBATCH --mem=8GB
#SBATCH --time=05:30:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-ZP-DM5
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=9500
export ENZ_ACCOUNT_LABEL=DM
export ENZ_PART_ID=5
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-zeroprior-initinfo-common.sh"
