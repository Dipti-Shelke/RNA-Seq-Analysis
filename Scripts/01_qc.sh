#!/bin/bash
#SBATCH --job-name=fastqc_analysis
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/data/projects/dipti/16s_projects/FTP_11052026B/logs/fastqc_%j.out
#SBATCH --error=/data/projects/dipti/16s_projects/FTP_11052026B/logs/fastqc_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=diptigulumbe@chromosomelabs.com

# =============================================================================
# FASTQC ANALYSIS — scratch-based workflow
# Reads input from NAS, processes on local scratch, copies results back to NAS
# =============================================================================

# =============================================================================
# INPUT / OUTPUT DIRECTORIES
# =============================================================================

INPUT_DIR="/data/projects/dipti/16s_projects/FTP_11052026B"
RESULT_DIR="/data/projects/dipti/16s_projects/FTP_11052026B/fastqc_results"
LOG_DIR="/data/projects/dipti/16s_projects/FTP_11052026B/logs"

THREADS=$SLURM_CPUS_PER_TASK

# =============================================================================
# SETUP SCRATCH
# =============================================================================

SCRATCH="/scratch/${USER}_${SLURM_JOB_ID}"

mkdir -p "$SCRATCH/input"
mkdir -p "$SCRATCH/output"
mkdir -p "$RESULT_DIR"
mkdir -p "$LOG_DIR"

echo "================================================="
echo "                FASTQC ANALYSIS                  "
echo "================================================="
echo "Job ID        : $SLURM_JOB_ID"
echo "Node          : $SLURMD_NODENAME"
echo "Start Time    : $(date)"
echo "Input Dir     : $INPUT_DIR"
echo "Result Dir    : $RESULT_DIR"
echo "Scratch Dir   : $SCRATCH"
echo "Threads       : $THREADS"
echo "================================================="

# =============================================================================
# ACTIVATE CONDA ENVIRONMENT
# =============================================================================

echo "[$(date)] Activating conda environment..."

source /local_conda/miniconda3/etc/profile.d/conda.sh
conda activate trimQC

echo "[$(date)] FastQC version:"
fastqc --version

# =============================================================================
# COPY ALL FASTQ FILES TO SCRATCH AUTOMATICALLY
# =============================================================================

echo "[$(date)] Copying FASTQ files to scratch..."

FASTQ_COUNT=$(ls "$INPUT_DIR"/*.fastq.gz 2>/dev/null | wc -l)

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $INPUT_DIR"
    exit 1
fi

cp "$INPUT_DIR"/*.fastq.gz "$SCRATCH/input/"

echo "[$(date)] Total FASTQ files copied: $FASTQ_COUNT"

# =============================================================================
# RUN FASTQC
# =============================================================================

echo "[$(date)] Running FastQC..."

fastqc \
    --threads "$THREADS" \
    --outdir "$SCRATCH/output" \
    --format fastq \
    "$SCRATCH/input/"*.fastq.gz

EXIT_CODE=$?

# =============================================================================
# CHECK FASTQC STATUS
# =============================================================================

if [ $EXIT_CODE -ne 0 ]; then
    echo "================================================="
    echo "ERROR: FastQC failed with exit code $EXIT_CODE"
    echo "================================================="

    rm -rf "$SCRATCH"
    exit $EXIT_CODE
fi

echo "[$(date)] FastQC completed successfully."

# =============================================================================
# COPY RESULTS BACK
# =============================================================================

echo "[$(date)] Copying results back to NAS..."

cp -r "$SCRATCH/output/"* "$RESULT_DIR/"

echo "[$(date)] Results copied successfully."

# =============================================================================
# CLEANUP
# =============================================================================

echo "[$(date)] Cleaning scratch directory..."

rm -rf "$SCRATCH"

# =============================================================================
# JOB COMPLETE
# =============================================================================

echo "================================================="
echo "FASTQC JOB COMPLETED SUCCESSFULLY"
echo "Completion Time : $(date)"
echo "Results Location: $RESULT_DIR"
echo "================================================="
