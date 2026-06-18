#!/bin/bash
#SBATCH --job-name=sortmerna_rRNAFree
#SBATCH --partition=all
#SBATCH --nodelist=bionode1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/data/projects/dipti/Test/logs/slurm_sortmerna_%j.out
#SBATCH --error=/data/projects/dipti/Test/logs/slurm_sortmerna_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=diptigulumbe@chromosomelabs.com

set -euo pipefail

PROJECT_DIR="/data/projects/dipti/Test"

INPUT_DIR="${PROJECT_DIR}/02_Trim_reads"
DB_DIR="/data/projects/dipti/Test/db"

OUTPUT_DIR="${PROJECT_DIR}/03_Sortmerna"
RRNAFREE_DIR="${OUTPUT_DIR}/rRNAFree"
SAMPLE_LOG_DIR="${OUTPUT_DIR}/sample_logs"
LOG_DIR="${PROJECT_DIR}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-25}"

# =============================================================================
# SELF-SUBMIT MODE: create output folders before SLURM starts
# =============================================================================

if [ -z "${SLURM_JOB_ID:-}" ]; then
    mkdir -p "$RRNAFREE_DIR"
    mkdir -p "$SAMPLE_LOG_DIR"
    mkdir -p "$LOG_DIR"

    echo "Output folders created:"
    echo "$RRNAFREE_DIR"
    echo "$LOG_DIR"

    echo "Submitting SortMeRNA job to SLURM..."
    sbatch "$0"
    exit 0
fi

mkdir -p "$RRNAFREE_DIR"
mkdir -p "$SAMPLE_LOG_DIR"
mkdir -p "$LOG_DIR"

# =============================================================================
# SCRATCH SETUP
# =============================================================================

SCRATCH="/scratch/${USER}_${SLURM_JOB_ID}_sortmerna"

SCRATCH_INPUT="${SCRATCH}/input"
SCRATCH_RRNAFREE="${SCRATCH}/rRNAFree"
SCRATCH_LOGS="${SCRATCH}/sortmerna_logs"
SCRATCH_WORKDIR="${SCRATCH}/workdir"

mkdir -p "$SCRATCH_INPUT"
mkdir -p "$SCRATCH_RRNAFREE"
mkdir -p "$SCRATCH_LOGS"
mkdir -p "$SCRATCH_WORKDIR"

# =============================================================================
# DATABASE REFERENCES
# =============================================================================

REFS=(
    "${DB_DIR}/silva-bac-16s-id90.fasta"
    "${DB_DIR}/silva-bac-23s-id98.fasta"
    "${DB_DIR}/silva-arc-16s-id95.fasta"
    "${DB_DIR}/silva-arc-23s-id98.fasta"
    "${DB_DIR}/silva-euk-18s-id95.fasta"
    "${DB_DIR}/silva-euk-28s-id98.fasta"
    "${DB_DIR}/rfam-5s-database-id98.fasta"
    "${DB_DIR}/rfam-5.8s-database-id98.fasta"
)

echo "================================================="
echo "              SORTMERNA rRNA REMOVAL             "
echo "================================================="
echo "Job ID          : ${SLURM_JOB_ID}"
echo "Node            : ${SLURMD_NODENAME:-NA}"
echo "Start Time      : $(date)"
echo "Project Dir     : $PROJECT_DIR"
echo "Input Dir       : $INPUT_DIR"
echo "Database Dir    : $DB_DIR"
echo "Output Dir      : $OUTPUT_DIR"
echo "rRNAFree Output : $RRNAFREE_DIR"
echo "SLURM Log Output: $LOG_DIR"
echo "Sample Log Dir  : $SAMPLE_LOG_DIR"
echo "Scratch Dir     : $SCRATCH"
echo "Threads         : $THREADS"
echo "================================================="

# =============================================================================
# CHECK INPUT FILES
# =============================================================================

FASTQ_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $INPUT_DIR"
    rm -rf "$SCRATCH"
    exit 1
fi

echo "[$(date)] FASTQ files found in input folder: $FASTQ_COUNT"

# =============================================================================
# CHECK DATABASE FILES
# =============================================================================

echo "[$(date)] Checking SortMeRNA database files..."

for ref in "${REFS[@]}"; do
    if [ ! -f "$ref" ]; then
        echo "ERROR: Missing database file: $ref"
        rm -rf "$SCRATCH"
        exit 1
    fi
done

echo "[$(date)] All database files found."

# =============================================================================
# ACTIVATE CONDA ENVIRONMENT
# =============================================================================

echo "[$(date)] Activating conda environment: sortmerna"

source /local_conda/miniconda3/etc/profile.d/conda.sh
conda activate sortmerna

echo "[$(date)] SortMeRNA path:"
which sortmerna

echo "[$(date)] SortMeRNA version:"
sortmerna --version

# =============================================================================
# COPY INPUT READS TO SCRATCH
# =============================================================================

echo "[$(date)] Copying input FASTQ files to scratch..."
cp "$INPUT_DIR"/*.fastq.gz "$SCRATCH_INPUT/"
echo "[$(date)] FASTQ files copied to scratch."

cd "$SCRATCH_INPUT"

# =============================================================================
# AUTOMATIC PAIRED-END DETECTION
# =============================================================================

shopt -s nullglob

R1_FILES=( *R1*.fastq.gz )

if [ "${#R1_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No R1 FASTQ files found in $INPUT_DIR"
    rm -rf "$SCRATCH"
    exit 1
fi

echo "[$(date)] Paired-end samples detected: ${#R1_FILES[@]}"

# =============================================================================
# RUN SORTMERNA
# =============================================================================

echo "[$(date)] Starting SortMeRNA..."

for R1 in "${R1_FILES[@]}"
do
    R2="${R1/R1/R2}"

    if [ ! -f "$R2" ]; then
        echo "ERROR: Matching R2 file not found for R1 file: $R1"
        echo "Expected R2 file: $R2"
        cp -r "$SCRATCH_LOGS/"* "$SAMPLE_LOG_DIR/" 2>/dev/null || true
        exit 1
    fi

    sample="$R1"
    sample="${sample%.fastq.gz}"
    sample="${sample/_R1_001/}"
    sample="${sample/_R1_trimmed/}"
    sample="${sample/_R1/}"

    SAMPLE_WORKDIR="${SCRATCH_WORKDIR}/${sample}"
    mkdir -p "$SAMPLE_WORKDIR"

    echo "================================================="
    echo "Processing sample : $sample"
    echo "R1                : $R1"
    echo "R2                : $R2"
    echo "Sample log        : ${SCRATCH_LOGS}/sortmerna_${sample}.log"
    echo "Start             : $(date)"
    echo "================================================="

    REF_ARGS=()
    for ref in "${REFS[@]}"; do
        REF_ARGS+=(--ref "$ref")
    done

    set +e

    sortmerna \
        "${REF_ARGS[@]}" \
        --reads "$R1" \
        --reads "$R2" \
        --paired_in \
        --fastx \
        --out2 \
        --other "${SCRATCH_RRNAFREE}/${sample}_rRNAFree" \
        --threads "$THREADS" \
        --workdir "$SAMPLE_WORKDIR" \
        > "${SCRATCH_LOGS}/sortmerna_${sample}.log" 2>&1

    EXIT_CODE=$?

    set -e

    if [ "$EXIT_CODE" -ne 0 ]; then
        echo "================================================="
        echo "ERROR: SortMeRNA failed for sample: $sample"
        echo "Exit code: $EXIT_CODE"
        echo "================================================="

        echo "Copying sample logs to:"
        echo "$SAMPLE_LOG_DIR"

        cp -r "$SCRATCH_LOGS/"* "$SAMPLE_LOG_DIR/" 2>/dev/null || true

        echo "Last 80 lines of SortMeRNA sample log:"
        echo "-------------------------------------------------"
        tail -n 80 "${SCRATCH_LOGS}/sortmerna_${sample}.log" || true
        echo "-------------------------------------------------"

        echo "Scratch kept for debugging:"
        echo "$SCRATCH"

        exit "$EXIT_CODE"
    fi

    echo "Completed sample : $sample"
    echo "End              : $(date)"
done

echo "[$(date)] SortMeRNA completed for all paired-end samples."

# =============================================================================
# MOVE ONLY REQUIRED OUTPUT FROM SCRATCH
# =============================================================================

echo "[$(date)] Moving only rRNAFree reads and SortMeRNA logs..."

cp -r "$SCRATCH_RRNAFREE/"* "$RRNAFREE_DIR/" 2>/dev/null || true
cp -r "$SCRATCH_LOGS/"* "$SAMPLE_LOG_DIR/" 2>/dev/null || true

echo "[$(date)] Required output files copied successfully."

# =============================================================================
# FINAL OUTPUT SUMMARY
# =============================================================================

echo "================================================="
echo "SORTMERNA OUTPUT SUMMARY"
echo "================================================="
echo "rRNAFree reads saved in:"
echo "$RRNAFREE_DIR"
echo ""
echo "Total rRNAFree output files:"
find "$RRNAFREE_DIR" -type f | wc -l

echo ""
echo "Sample-wise SortMeRNA logs saved in:"
echo "$SAMPLE_LOG_DIR"
echo ""
echo "Total SortMeRNA log files:"
find "$SAMPLE_LOG_DIR" -type f | wc -l
echo "================================================="

# =============================================================================
# CLEAN SCRATCH ONLY AFTER SUCCESS
# =============================================================================

echo "[$(date)] Removing scratch directory after successful run..."
rm -rf "$SCRATCH"

echo "================================================="
echo "SORTMERNA JOB COMPLETED SUCCESSFULLY"
echo "Completion Time  : $(date)"
echo "rRNAFree Output  : $RRNAFREE_DIR"
echo "SLURM Log Output : $LOG_DIR"
echo "Sample Log Output: $SAMPLE_LOG_DIR"
echo "================================================="
