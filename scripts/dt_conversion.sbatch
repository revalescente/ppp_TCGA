#!/bin/bash
#SBATCH --job-name=dt_conversion
#SBATCH --ntasks=1                      # 1 task (the R script itself)
#SBATCH --cpus-per-task=16              # Number of CPU cores for mclapply
#SBATCH --mem=50G                       # Total memory (adjust based on data size)
#SBATCH --time=02:00:00                 # Max runtime (HH:MM:SS)
#SBATCH --output=/projects/shared/TCGA_data/logs/slurm-%j.out
#SBATCH --error=/projects/shared/TCGA_data/logs/slurm-%j.err

# ---------------------------------------------------------
# Environment Construction
# ---------------------------------------------------------

# Clear any loaded modules to ensure a clean environment
module purge
module load micromamba/rpy-sc

# ---------------------------------------------------------
# Execution
# ---------------------------------------------------------

echo "Job started on $(hostname) at $(date)"
echo "Allocated $SLURM_CPUS_PER_TASK cores."

# Execute the R script
Rscript /home/vreffo/ppp_TCGA/src/pipeline/dt_creation.R

echo "Job completed at $(date)"