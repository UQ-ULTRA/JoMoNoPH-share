#!/bin/bash -l
# Rescue run for the JMbayes2 all-40-config 1000-replicate randSDsummary run.
#
# Targets from the 2026-06-02/03 production rerun:
#   - 1 missing .rds batch caused by serverSocket failure before fitting started
#   - 22 saved batches with at least one model_status != "ok"
#
# The local array index maps to an explicit original global sim_id. The R driver
# computes global sim_id as SLURM_ARRAY_TASK_ID + SLURM_ARRAY_TASK_OFFSET, so the
# dynamic offset below preserves the original config/batch mapping exactly.
#
# Rescue output is written to a separate -rescue folder and does not overwrite
# the production output folder.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-23
#SBATCH --mem=12GB
#SBATCH --time=03:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZ-JMB2-RS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Original global sim_ids under:
#   SIM_BASELINE_FILTER=, SIMS_PER_CONFIG=1000, SIM_BATCH_SIZE=10.
GLOBAL_SIM_IDS=(
  534
  879
  893
  1300
  1502
  1517
  1606
  1835
  2301
  2302
  2775
  2929
  3140
  3283
  3440
  3620
  3633
  3638
  3649
  3723
  3804
  3844
  3895
)

RESCUE_LABELS=(
  "cfg06_LL120_scen3_admin_batch034_model_fail"
  "cfg09_LL120_scen5_cens50_batch079_missing_socket"
  "cfg09_LL120_scen5_cens50_batch093_model_fail"
  "cfg13_W090_scen2_cens50_batch100_model_fail"
  "cfg16_W090_scen3_admin_batch002_model_fail"
  "cfg16_W090_scen3_admin_batch017_model_fail"
  "cfg17_W090_scen4_cens50_batch006_model_fail"
  "cfg19_W090_scen5_cens50_batch035_model_fail"
  "cfg24_W130_scen2_admin_batch001_model_fail"
  "cfg24_W130_scen2_admin_batch002_model_fail"
  "cfg28_W130_scen4_admin_batch075_model_fail"
  "cfg30_W130_scen5_admin_batch029_model_fail"
  "cfg32_W210_scen1_admin_batch040_timeout"
  "cfg33_W210_scen2_cens50_batch083_model_fail"
  "cfg35_W210_scen3_cens50_batch040_model_fail"
  "cfg37_W210_scen4_cens50_batch020_model_fail"
  "cfg37_W210_scen4_cens50_batch033_model_fail"
  "cfg37_W210_scen4_cens50_batch038_model_fail"
  "cfg37_W210_scen4_cens50_batch049_model_fail"
  "cfg38_W210_scen4_admin_batch023_model_fail"
  "cfg39_W210_scen5_cens50_batch004_model_fail"
  "cfg39_W210_scen5_cens50_batch044_model_fail"
  "cfg39_W210_scen5_cens50_batch095_model_fail"
)

IDX=$((SLURM_ARRAY_TASK_ID - 1))
GLOBAL_SIM_ID="${GLOBAL_SIM_IDS[$IDX]}"
RESCUE_LABEL="${RESCUE_LABELS[$IDX]}"

export SLURM_ARRAY_TASK_OFFSET=$((GLOBAL_SIM_ID - SLURM_ARRAY_TASK_ID))
export ENZ_JMB2_SLICE_LABEL=RANDSD-RESCUE

export JOB_OUTPUT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
JOB_OUTPUT_DATE_FILE="$HOME/JoMoNoPH-share/Bunya/ENZAMET/.out_date_${JOB_OUTPUT_ID}"
JOB_OUTPUT_DATE_LOCK="${JOB_OUTPUT_DATE_FILE}.lock"
if [ ! -s "$JOB_OUTPUT_DATE_FILE" ]; then
  if mkdir "$JOB_OUTPUT_DATE_LOCK" 2>/dev/null; then
    date +%Y%m%d > "$JOB_OUTPUT_DATE_FILE"
    rmdir "$JOB_OUTPUT_DATE_LOCK"
  else
    while [ ! -s "$JOB_OUTPUT_DATE_FILE" ]; do sleep 1; done
  fi
fi

export JOB_OUTPUT_DATE="$(cat "$JOB_OUTPUT_DATE_FILE")"
export JOB_OUTPUT_DIR="$HOME/out_${JOB_OUTPUT_DATE}_${JOB_OUTPUT_ID}"
export JOB_OUTPUT_TASK_SUFFIX="${SLURM_ARRAY_TASK_ID:+_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$JOB_OUTPUT_DIR"
exec > "$JOB_OUTPUT_DIR/${SLURM_JOB_NAME}_${JOB_OUTPUT_ID}${JOB_OUTPUT_TASK_SUFFIX}.out" 2>&1

module purge
module load r/4.4.0-gfbf-2023a

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

export SIM_OUTPUT_DATE="$JOB_OUTPUT_DATE"
export SIM_OUTPUT_JOB_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"

# Match the production grid exactly.
export SIM_BASELINE_FILTER=
export SIM_BATCH_SIZE=10
export SIMS_PER_CONFIG=1000

# Rescue is intentionally more generous than production. The production run used
# 480 seconds per replicate under a 45-minute Slurm limit. Here 1800 seconds
# allows slow or retried replicates to finish, and the 3-hour walltime covers
# 10 sequential replicates with margin.
export SIM_SAFE_ELAPSED=1800

export SIM_JMB2_MODELS=all
export SIM_JMB2_SAVE_MODE=summary

# Keep the same JMbayes2 MCMC settings as the randSDsummary production run.
export SIM_JMB2_N_ITER=3500
export SIM_JMB2_N_BURNIN=500
export SIM_JMB2_N_THIN=1
export SIM_JMB2_N_CHAINS=4

# Use one replicate worker per Slurm task to reduce socket pressure. Each
# JMbayes2 fit still starts with 4 cores/chains; if that fails, retry once with
# cores=1 to rescue socket-related or parallel-backend failures.
export SIM_JMB2_CORES=4
export SIM_JMB2_RETRY_CORES=1

export SIM_OUTPUT_DIR="$HOME/JoMoNoPH-share/Bunya/ENZAMET/output_20260602_ENZAMET-JMbayes2-DM-all40-1000rep-randSDsummary-rescue"
mkdir -p "$SIM_OUTPUT_DIR"

echo "JMbayes2 all40 randSDsummary rescue"
echo "Target global sim_id: $GLOBAL_SIM_ID"
echo "Rescue label: $RESCUE_LABEL"
echo "Dynamic SLURM_ARRAY_TASK_OFFSET: $SLURM_ARRAY_TASK_OFFSET"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Array Task Offset: $SLURM_ARRAY_TASK_OFFSET"
echo "JMbayes2 slice label: $ENZ_JMB2_SLICE_LABEL"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: JMbayes2 ENZAMET comparison; randSDsummary rescue; all 40 configs; 1000 reps/config"
echo "SIM_BASELINE_FILTER: $SIM_BASELINE_FILTER"
echo "SIM_BATCH_SIZE: $SIM_BATCH_SIZE"
echo "SIMS_PER_CONFIG: $SIMS_PER_CONFIG"
echo "SIM_SAFE_ELAPSED: $SIM_SAFE_ELAPSED"
echo "SIM_JMB2_MODELS: $SIM_JMB2_MODELS"
echo "SIM_JMB2_SAVE_MODE: $SIM_JMB2_SAVE_MODE"
echo "SIM_JMB2_N_ITER: $SIM_JMB2_N_ITER"
echo "SIM_JMB2_N_BURNIN: $SIM_JMB2_N_BURNIN"
echo "SIM_JMB2_N_CHAINS: $SIM_JMB2_N_CHAINS"
echo "SIM_JMB2_CORES: $SIM_JMB2_CORES"
echo "SIM_JMB2_RETRY_CORES: $SIM_JMB2_RETRY_CORES"
echo "Output directory: $SIM_OUTPUT_DIR"
echo "SLURM .out directory: $JOB_OUTPUT_DIR"

srun --cpu-bind=cores Rscript "$HOME/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET-JMbayes2-DM.R"
