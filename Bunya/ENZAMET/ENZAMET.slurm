#!/bin/bash -l 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=20
## change to 40 for larger simulation
#SBATCH --array=1-100
#SBATCH --mem=20GB
## change to 40GB for larger simulation
#SBATCH --time=2:00:00
## change to 8 for larger simulation
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=ENZAMET

module purge
module load r/4.4.0-gfbf-2023a 

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Host: $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Model: BP/GB + LMM comparison"

srun --cpu-bind=cores Rscript ~/JoMoNoPH-share/Bunya/ENZAMET/ENZAMET.R