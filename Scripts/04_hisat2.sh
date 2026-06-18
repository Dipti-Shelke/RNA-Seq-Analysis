#!/bin/bash
#SBATCH --job-name=hisat2_mapping
#SBATCH --partition=all
#SBATCH --nodelist=bionode1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/data/projects/dipti/Test/logs/slurm_hisat2_%j.out
#SBATCH --error=/data/projects/dipti/Test/logs/slurm_hisat2_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=diptigulumbe@chromosomelabs.com

set -euo pipefail

# =============================================================================
# PROJECT PATHS
# =============================================================================

PROJECT_DIR="/data/projects/dipti/Test"

# SortMeRNA rRNA-free reads
INPUT_DIR="${PROJECT_DIR}/eukaryotic/03_Sortmerna"

# Genome FASTA for building index
GENOME_FASTA="${PROJECT_DIR}/genome/GRCm39_genomic.fna"

# Index will be built here: genome/mm_index/mm39_index
INDEX_DIR="${PROJECT_DIR}/genome/mm_index"
HISAT2_INDEX="${INDEX_DIR}/mm39_index"

# HISAT2 output
OUTPUT_DIR="${PROJECT_DIR}/04_hisat2"
BAM_DIR="${OUTPUT_DIR}"
SAMPLE_LOG_DIR="${OUTPUT_DIR}/sample_logs"
LOG_DIR="${PROJECT_DIR}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-30}"

# =============================================================================
# SELF-SUBMIT MODE
# =============================================================================

if [ -z "${SLURM_JOB_ID:-}" ]; then
    mkdir -p "$INDEX_DIR"
    mkdir -p "$BAM_DIR"
    mkdir -p "$SAMPLE_LOG_DIR"
    mkdir -p "$LOG_DIR"

    echo "Folders created:"
    echo "  $INDEX_DIR"
    echo "  $BAM_DIR"
    echo "  $SAMPLE_LOG_DIR"
    echo "  $LOG_DIR"
    echo ""
    echo "Submitting HISAT2 index + mapping job to SLURM..."
    sbatch "$0"
    exit 0
fi

mkdir -p "$INDEX_DIR"
mkdir -p "$BAM_DIR"
mkdir -p "$SAMPLE_LOG_DIR"
mkdir -p "$LOG_DIR"

# =============================================================================
# SCRATCH SETUP
# =============================================================================

SCRATCH="/scratch/${USER}_${SLURM_JOB_ID}_hisat2"
SCRATCH_INPUT="${SCRATCH}/input"
SCRATCH_BAM="${SCRATCH}/bam"
SCRATCH_LOGS="${SCRATCH}/sample_logs"
SCRATCH_TMP="${SCRATCH}/tmp"

mkdir -p "$SCRATCH_INPUT" "$SCRATCH_BAM" "$SCRATCH_LOGS" "$SCRATCH_TMP"

cleanup_on_error() {
    echo "ERROR: Job failed. Copying available output before exit..."
    cp -r "$SCRATCH_BAM"/  "$BAM_DIR/"        2>/dev/null || true
    cp -r "$SCRATCH_LOGS"/ "$SAMPLE_LOG_DIR/" 2>/dev/null || true
    echo "Scratch kept for debugging: $SCRATCH"
}
trap cleanup_on_error ERR

echo "================================================="
echo "         HISAT2 INDEX BUILD + MAPPING            "
echo "================================================="
echo "Job ID           : ${SLURM_JOB_ID}"
echo "Node             : ${SLURMD_NODENAME:-NA}"
echo "Start Time       : $(date)"
echo "Project Dir      : $PROJECT_DIR"
echo "Input Dir        : $INPUT_DIR"
echo "Genome FASTA     : $GENOME_FASTA"
echo "Index Dir        : $INDEX_DIR"
echo "HISAT2 Index     : $HISAT2_INDEX"
echo "Output Dir       : $OUTPUT_DIR"
echo "Scratch Dir      : $SCRATCH"
echo "Threads          : $THREADS"
echo "================================================="

# =============================================================================
# ACTIVATE CONDA
# =============================================================================

echo "[$(date)] Activating conda environment: mapper"
source /local_conda/miniconda3/etc/profile.d/conda.sh
conda activate mapper

echo "[$(date)] hisat2     : $(which hisat2)  | $(hisat2 --version | head -n 1)"
echo "[$(date)] samtools   : $(which samtools) | $(samtools --version | head -n 1)"

# =============================================================================
# PHASE 1 — BUILD HISAT2 INDEX (skip if already built)
# =============================================================================

echo ""
echo "================================================="
echo "PHASE 1 : HISAT2 INDEX BUILD"
echo "================================================="

if ls "${HISAT2_INDEX}".*.ht2 >/dev/null 2>&1 || \
   ls "${HISAT2_INDEX}".*.ht2l >/dev/null 2>&1; then

    echo "[$(date)] Index already exists — skipping build."
    echo "  Index prefix : $HISAT2_INDEX"
    ls -lh "${INDEX_DIR}/"
else
    echo "[$(date)] No index found. Building HISAT2 index..."

    # Check genome FASTA exists
    if [ ! -f "$GENOME_FASTA" ]; then
        echo "ERROR: Genome FASTA not found: $GENOME_FASTA"
        rm -rf "$SCRATCH"
        exit 1
    fi

    echo "[$(date)] Genome FASTA : $GENOME_FASTA ($(du -sh "$GENOME_FASTA" | cut -f1))"
    echo "[$(date)] Index prefix : $HISAT2_INDEX"
    echo "[$(date)] Threads      : $THREADS"

    INDEX_LOG="${LOG_DIR}/hisat2_build_${SLURM_JOB_ID}.log"

    hisat2-build \
        -p "$THREADS" \
        "$GENOME_FASTA" \
        "$HISAT2_INDEX" \
        > "$INDEX_LOG" 2>&1

    if [ $? -ne 0 ]; then
        echo "ERROR: hisat2-build failed. See log: $INDEX_LOG"
        tail -n 40 "$INDEX_LOG" || true
        rm -rf "$SCRATCH"
        exit 1
    fi

    echo "[$(date)] Index build complete."
    echo "Index files created:"
    ls -lh "${INDEX_DIR}/"
fi

# Verify index is usable
if ! ls "${HISAT2_INDEX}".*.ht2 >/dev/null 2>&1 && \
   ! ls "${HISAT2_INDEX}".*.ht2l >/dev/null 2>&1; then
    echo "ERROR: Index still not found after build step."
    echo "  Expected: ${HISAT2_INDEX}.1.ht2 or ${HISAT2_INDEX}.1.ht2l"
    rm -rf "$SCRATCH"
    exit 1
fi

# =============================================================================
# PHASE 2 — CHECK INPUT READS
# =============================================================================

echo ""
echo "================================================="
echo "PHASE 2 : HISAT2 MAPPING"
echo "================================================="

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    rm -rf "$SCRATCH"
    exit 1
fi

# Count all fastq/fq files (parentheses required around -o)
FASTQ_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type f \
    \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | wc -l)

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $INPUT_DIR"
    echo "Contents:"
    ls -lh "$INPUT_DIR" || true
    rm -rf "$SCRATCH"
    exit 1
fi

echo "[$(date)] FASTQ files found: $FASTQ_COUNT"
echo "Files:"
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) \
    | sort | while read -r f; do echo "  $(basename "$f")"; done

# =============================================================================
# COPY INPUT READS TO SCRATCH
# =============================================================================

echo "[$(date)] Copying reads to scratch..."
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) \
    -exec cp {} "$SCRATCH_INPUT/" \;
echo "[$(date)] Copy done."

cd "$SCRATCH_INPUT"

# =============================================================================
# PAIRED-END DETECTION
# Supports all common naming conventions including SortMeRNA output
#
#   Pattern                          R1 suffix               R2 suffix
#   -------                          ---------               ---------
#   SortMeRNA non_rRNA  :  _1.non_rRNA.fastq.gz   _2.non_rRNA.fastq.gz
#   SortMeRNA rRNAFree  :  _rRNAFree_fwd.fq.gz    _rRNAFree_rev.fq.gz
#   Standard _R1/_R2    :  _R1.fastq.gz            _R2.fastq.gz
#   Standard _R1/_R2    :  _R1.fq.gz               _R2.fq.gz
#   Illumina _001       :  _R1_001.fastq.gz         _R2_001.fastq.gz
#   Illumina _001       :  _R1_001.fq.gz            _R2_001.fq.gz
#   Simple 1/2          :  _1.fastq.gz              _2.fastq.gz
#   Simple 1/2          :  _1.fq.gz                 _2.fq.gz
# =============================================================================

shopt -s nullglob

R1_FILES=(
    *_1.non_rRNA.fastq.gz
    *_1.non_rRNA.fq.gz
    *_rRNAFree_fwd.fq.gz
    *_rRNAFree_fwd.fastq.gz
    *_R1_001.fastq.gz
    *_R1_001.fq.gz
    *_R1.fastq.gz
    *_R1.fq.gz
    *_1.fastq.gz
    *_1.fq.gz
)

if [ "${#R1_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No R1 files detected in scratch input."
    echo "Files present:"
    ls -lh "$SCRATCH_INPUT"
    echo ""
    echo "Supported R1 naming patterns:"
    echo "  *_1.non_rRNA.fastq.gz  (SortMeRNA non_rRNA output)"
    echo "  *_rRNAFree_fwd.fq.gz   (SortMeRNA rRNAFree output)"
    echo "  *_R1.fastq.gz / *_R1.fq.gz"
    echo "  *_R1_001.fastq.gz / *_R1_001.fq.gz"
    echo "  *_1.fastq.gz / *_1.fq.gz"
    rm -rf "$SCRATCH"
    exit 1
fi

echo "[$(date)] Paired-end samples detected: ${#R1_FILES[@]}"
for f in "${R1_FILES[@]}"; do echo "  R1: $f"; done

# =============================================================================
# RUN HISAT2 → SAMTOOLS VIEW → SAMTOOLS SORT
# =============================================================================

echo ""
echo "[$(date)] Starting HISAT2 mapping..."

for R1 in "${R1_FILES[@]}"
do
    # ------------------------------------------------------------------
    # Match R1 pattern → derive R2 filename and clean sample name
    # ------------------------------------------------------------------
    if   [[ "$R1" == *_1.non_rRNA.fastq.gz ]]; then
        R2="${R1/_1.non_rRNA.fastq.gz/_2.non_rRNA.fastq.gz}"
        sample="${R1%_1.non_rRNA.fastq.gz}"

    elif [[ "$R1" == *_1.non_rRNA.fq.gz ]]; then
        R2="${R1/_1.non_rRNA.fq.gz/_2.non_rRNA.fq.gz}"
        sample="${R1%_1.non_rRNA.fq.gz}"

    elif [[ "$R1" == *_rRNAFree_fwd.fq.gz ]]; then
        R2="${R1/_rRNAFree_fwd.fq.gz/_rRNAFree_rev.fq.gz}"
        sample="${R1%_rRNAFree_fwd.fq.gz}"

    elif [[ "$R1" == *_rRNAFree_fwd.fastq.gz ]]; then
        R2="${R1/_rRNAFree_fwd.fastq.gz/_rRNAFree_rev.fastq.gz}"
        sample="${R1%_rRNAFree_fwd.fastq.gz}"

    elif [[ "$R1" == *_R1_001.fastq.gz ]]; then
        R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
        sample="${R1%_R1_001.fastq.gz}"

    elif [[ "$R1" == *_R1_001.fq.gz ]]; then
        R2="${R1/_R1_001.fq.gz/_R2_001.fq.gz}"
        sample="${R1%_R1_001.fq.gz}"

    elif [[ "$R1" == *_R1.fastq.gz ]]; then
        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        sample="${R1%_R1.fastq.gz}"

    elif [[ "$R1" == *_R1.fq.gz ]]; then
        R2="${R1/_R1.fq.gz/_R2.fq.gz}"
        sample="${R1%_R1.fq.gz}"

    elif [[ "$R1" == *_1.fastq.gz ]]; then
        R2="${R1/_1.fastq.gz/_2.fastq.gz}"
        sample="${R1%_1.fastq.gz}"

    elif [[ "$R1" == *_1.fq.gz ]]; then
        R2="${R1/_1.fq.gz/_2.fq.gz}"
        sample="${R1%_1.fq.gz}"

    else
        echo "ERROR: Cannot parse R1 filename: $R1"
        exit 1
    fi

    if [ ! -f "$R2" ]; then
        echo "ERROR: R2 file not found for R1: $R1"
        echo "  Expected R2: $R2"
        echo "  Files in scratch:"
        ls "$SCRATCH_INPUT"
        exit 1
    fi

    SAMPLE_LOG="${SCRATCH_LOGS}/hisat2_${sample}.log"
    SORT_LOG="${SCRATCH_LOGS}/samtools_${sample}.log"
    BAM_OUT="${SCRATCH_BAM}/${sample}.bam"

    echo "================================================="
    echo "Sample   : $sample"
    echo "R1       : $R1"
    echo "R2       : $R2"
    echo "BAM out  : $BAM_OUT"
    echo "Start    : $(date)"
    echo "================================================="

    set +e

    hisat2 \
        -x "$HISAT2_INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS" \
        2> "$SAMPLE_LOG" \
    | samtools view -@ "$THREADS" -bS - 2>> "$SORT_LOG" \
    | samtools sort -@ "$THREADS" \
                    -T "${SCRATCH_TMP}/${sample}_sorttmp" \
                    -o "$BAM_OUT" - 2>> "$SORT_LOG"

    # Save ALL exit codes at once before PIPESTATUS is overwritten
    PIPE_STATUS=( "${PIPESTATUS[@]}" )
    HISAT2_EC="${PIPE_STATUS[0]}"
    VIEW_EC="${PIPE_STATUS[1]:-0}"
    SORT_EC="${PIPE_STATUS[2]:-0}"

    set -e

    if [ "$HISAT2_EC" -ne 0 ] || [ "$VIEW_EC" -ne 0 ] || [ "$SORT_EC" -ne 0 ]; then
        echo "ERROR: Pipeline failed for sample: $sample"
        echo "  hisat2 exit        : $HISAT2_EC"
        echo "  samtools view exit : $VIEW_EC"
        echo "  samtools sort exit : $SORT_EC"
        echo "--- HISAT2 log (last 60 lines) ---"
        tail -n 60 "$SAMPLE_LOG" || true
        echo "--- Samtools log (last 60 lines) ---"
        tail -n 60 "$SORT_LOG"   || true
        exit 1
    fi

    if [ ! -s "$BAM_OUT" ]; then
        echo "ERROR: BAM file empty or missing: $BAM_OUT"
        exit 1
    fi

    echo "[$(date)] Indexing BAM: $sample"
    samtools index -@ "$THREADS" "$BAM_OUT" >> "$SORT_LOG" 2>&1

    # Print alignment summary from HISAT2 log
    echo "--- Alignment summary for $sample ---"
    cat "$SAMPLE_LOG" || true
    echo "--- End summary ---"

    echo "Done: $sample  ($(date))"
done

echo ""
echo "[$(date)] All samples mapped successfully."

# =============================================================================
# COPY OUTPUT FROM SCRATCH TO FINAL DESTINATION
# =============================================================================

echo "[$(date)] Copying BAM files and logs to final output..."
cp -r "$SCRATCH_BAM"/  "$BAM_DIR/"        2>/dev/null || true
cp -r "$SCRATCH_LOGS"/ "$SAMPLE_LOG_DIR/" 2>/dev/null || true
echo "[$(date)] Copy complete."

# =============================================================================
# FINAL SUMMARY
# =============================================================================

echo ""
echo "================================================="
echo "FINAL OUTPUT SUMMARY"
echo "================================================="
echo "Index location   : $INDEX_DIR"
ls -lh "${INDEX_DIR}/" | grep -v "^total" || true
echo ""
echo "BAM directory    : $BAM_DIR"
echo "BAM files        : $(find "$BAM_DIR" -maxdepth 1 -name "*.bam"     | wc -l)"
echo "BAI files        : $(find "$BAM_DIR" -maxdepth 1 -name "*.bam.bai" | wc -l)"
echo ""
echo "Sample logs      : $SAMPLE_LOG_DIR"
echo "Log files        : $(find "$SAMPLE_LOG_DIR" -type f | wc -l)"
echo "================================================="

# =============================================================================
# CLEANUP SCRATCH
# =============================================================================

trap - ERR

echo "[$(date)] Removing scratch..."
rm -rf "$SCRATCH"

echo ""
echo "================================================="
echo "JOB COMPLETED SUCCESSFULLY"
echo "Completion Time  : $(date)"
echo "BAM Output       : $BAM_DIR"
echo "================================================="
