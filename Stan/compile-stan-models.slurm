#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G
#SBATCH --job-name=CompStan
#SBATCH --account=a_uqccr
#SBATCH --time=00:10:00
#SBATCH --partition=general

# --- Load required modules ---
module purge
module load r

# --- Run R script ---

Rscript compile-stan-models.R
