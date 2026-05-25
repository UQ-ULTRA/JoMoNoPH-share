#!/bin/bash -l
# MZ part 4: global array tasks 6501-7000, configs 27-28.
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
#SBATCH --job-name=ENZ-ZP-MZ4
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=6500
export ENZ_ACCOUNT_LABEL=MZ
export ENZ_PART_ID=4
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-zeroprior-initinfo-common.sh"
