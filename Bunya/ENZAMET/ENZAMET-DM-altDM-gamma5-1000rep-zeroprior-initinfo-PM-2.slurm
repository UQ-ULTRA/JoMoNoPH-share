#!/bin/bash -l
# PM part 2: global array tasks 3001-3500, configs 13-14.
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
#SBATCH --job-name=ENZ-ZP-PM2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=3000
export ENZ_ACCOUNT_LABEL=PM
export ENZ_PART_ID=2
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-1000rep-zeroprior-initinfo-common.sh"
