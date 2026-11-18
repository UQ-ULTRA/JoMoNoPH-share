#!/bin/bash -l
# JMbayes2 all-40-config 1000-replicate production slice 01:
# global array tasks 1-800.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-800
#SBATCH --mem=12GB
#SBATCH --time=00:45:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-JMB2-01
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=0
export ENZ_JMB2_SLICE_LABEL=01

source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-jmbayes2-dm-all40-1000rep-common.sh"
