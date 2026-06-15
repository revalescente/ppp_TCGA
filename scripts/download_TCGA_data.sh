#!/bin/bash
#SBATCH --job-name=download_TCGA
#SBATCH --partition=cpu
#SBATCH --time=06-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --output=ppp_TCGA/outputs/slurm-%j.out
#SBATCH --error=ppp_TCGA/logs/slurm-%j.err

# specify the paths of the out and err to have them created in a sensible position

set -euo pipefail


echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start time: $(date)"

# Load the required environment
module purge
module load micromamba/rpy-ml

Rscript ppp_TCGA/src/TCGA_data_download.R
