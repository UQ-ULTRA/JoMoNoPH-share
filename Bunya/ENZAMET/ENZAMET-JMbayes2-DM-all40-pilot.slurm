#!/bin/bash -l
# Pilot run for JMbayes2 comparison timing/storage across all 40 ENZAMET configs.
#
# Defaults here:
#   - all 40 configs: 4 baselines x scenarios 1:5 x censoring 50/admin
#   - 10 replicates per config
#   - 5 replicates per batch, run sequentially within each task with 4 JMbayes2 cores
#   - all three JMbayes2 models in each replicate: Cox, Weibull, loglogistic
#
# Array size:
#   40 configs x 2 batches/config = 80 tasks.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-80
#SBATCH --mem=8GB
#SBATCH --time=00:45:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-JMB2-40
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=0

# Empty baseline filter means use the full 40-config grid from the R driver.
export SIM_BASELINE_FILTER_OVERRIDE=
export SIM_OUTPUT_LABEL_OVERRIDE=ENZAMET-JMbayes2-DM-all40-pilot-batch5
export SIM_BATCH_SIZE_OVERRIDE=5
export SIMS_PER_CONFIG_OVERRIDE=10
# 8 minutes per replicate gives each 5-replicate batch up to 40 minutes of
# model-fitting time, leaving roughly 5 minutes for overhead under the 45-minute
# Slurm wall time.
export SIM_SAFE_ELAPSED_OVERRIDE=480
export SIM_JMB2_MODELS_OVERRIDE=all
export SIM_JMB2_SAVE_MODE_OVERRIDE=summary
export SIM_JMB2_N_ITER_OVERRIDE=3500
export SIM_JMB2_N_BURNIN_OVERRIDE=500

# JMbayes2 defaults to 3 chains; use 4 chains here to align diagnostics with
# the Stan sAFT-JM runs while keeping JMbayes2's default iteration/burn-in
# lengths, which may be better tuned for its MCMC scheme.
export SIM_JMB2_N_CHAINS_OVERRIDE=4
export SIM_JMB2_CORES_OVERRIDE=4

source "$HOME/JoMoNoPH-share/Bunya/ENZAMET/run-enzamet-jmbayes2-dm-common.sh"
