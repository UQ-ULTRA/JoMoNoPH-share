#!/bin/bash -l
# Pilot run for JMbayes2 comparison timing/storage.
#
# Defaults from run-enzamet-jmbayes2-dm-common.sh:
#   - 10 log-logistic configs: LL120_230 x scenarios 1:5 x censoring 50/admin
#   - 2 replicates per config
#   - 1 replicate per batch
#   - all three JMbayes2 models in each replicate: Cox, Weibull, loglogistic
#
# Array size:
#   10 configs x 2 batches/config = 20 tasks.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-20
#SBATCH --mem=24GB
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-JMB2-LL
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export SLURM_ARRAY_TASK_OFFSET=0

# Keep these here so the pilot is self-documenting. Override before sbatch if
# you want to split models or increase replicates; remember to update --array.
export SIM_BASELINE_FILTER_OVERRIDE=LL120_230
export SIM_BATCH_SIZE_OVERRIDE=1
export SIMS_PER_CONFIG_OVERRIDE=2
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
