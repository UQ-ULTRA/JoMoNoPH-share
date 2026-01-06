#!/bin/bash -l 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=bun[083-115,126-143]
#SBATCH --exclude=bun091
#SBATCH --cpus-per-task=40
#SBATCH --array=1-100
#SBATCH --mem=40GB
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH --account=a_uqccr
#SBATCH --job-name=I2a-ENZAMET

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

srun --cpu-bind=cores Rscript /home/uqpmahe1/JoMoNoPH/Bunya/ENZAMET/I2a-ENZAMET.R