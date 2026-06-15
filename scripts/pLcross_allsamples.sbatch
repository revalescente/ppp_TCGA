#!/bin/bash
#SBATCH --job-name=pLcross_allsamples
#SBATCH --array=1-116%10
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=01-00:00:00
#SBATCH --output=/projects/shared/TCGA_data/logs/slurm-%j.out
#SBATCH --error=/projects/shared/TCGA_data/logs/slurm-%j.err

module purge
module load micromamba/rpy-sc

INPUT_LIST=/projects/shared/TCGA_data/missing_tcga.txt

INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")

# Extract the base name of the input to construct the expected output filename
BASENAME=$(basename "$INPUT_FILE")
NAME="${BASENAME%.*}"

OUTPUT_DIR=/projects/shared/TCGA_data/results/Lfun_KDE

# Define the exact file you expect the R script to create.
# CHANGE '.rds' to the actual file extension your R script saves!
EXPECTED_OUTPUT="${OUTPUT_DIR}/${NAME}.rds"

mkdir -p "$OUTPUT_DIR"

# Check if the output file already exists
if [ -f "$EXPECTED_OUTPUT" ]; then
    echo "Skipping task ${SLURM_ARRAY_TASK_ID}: ${EXPECTED_OUTPUT} already exists."
    exit 0
fi

echo "Processing task ${SLURM_ARRAY_TASK_ID}: $INPUT_FILE"

micromamba run Rscript ppp_TCGA/src/pipeline/Lcross_estimate.R --sample_dir "$INPUT_FILE" \
--save_dir "$OUTPUT_DIR" --AppMag 40 --nthreads 2
