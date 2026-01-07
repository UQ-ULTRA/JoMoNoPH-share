#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
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
