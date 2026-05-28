#!/bin/bash -l
# Alpha sensitivity max_f=1.5, MZ slice: global array tasks 626-1250.
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=8
#SBATCH --array=1-625
#SBATCH --mem=8GB
#SBATCH --time=05:30:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-A150-MZ
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=625
export ENZ_ACCOUNT_LABEL=MZ
export ENZ_ALPHA_LABEL=maxF150
export SIM_ALPHA_MAX_F=1.5
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-alphaSens-LL-1000rep-common.sh"
