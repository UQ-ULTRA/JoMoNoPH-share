#!/bin/bash -l
# Alpha sensitivity max_f=1.25, PM slice: global array tasks 1251-1875.
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
#SBATCH --job-name=ENZ-A125-PM
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=1250
export ENZ_ACCOUNT_LABEL=PM
export ENZ_ALPHA_LABEL=maxF125
export SIM_ALPHA_MAX_F=1.25
source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-altDM-gamma5-alphaSens-LL-1000rep-common.sh"
