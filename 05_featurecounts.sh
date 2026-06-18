#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --partition=all
#SBATCH --nodelist=bionode1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=/data/projects/dipti/Test/logs/slurm_featurecounts_%j.out
#SBATCH --error=/data/projects/dipti/Test/logs/slurm_featurecounts_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=diptigulumbe@chromosomelabs.com

set -euo pipefail

# =============================================================================
# PROJECT PATHS
# =============================================================================

PROJECT_DIR="/data/projects/dipti/Test"

# Sorted BAM files from HISAT2
INPUT_DIR="${PROJECT_DIR}/04_hisat2/bam"

# GTF annotation
GTF="${PROJECT_DIR}/genome/GRCm39_genomic.gtf"

# featureCounts output
OUTPUT_DIR="${PROJECT_DIR}/05_featurecounts"
SAMPLE_LOG_DIR="${OUTPUT_DIR}/sample_logs"
LOG_DIR="${PROJECT_DIR}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-25}"

# =============================================================================
# featureCounts PARAMETERS — EUKARYOTIC
# -p               paired-end mode
# --countReadPairs count fragments (pairs), not individual reads
# -B               both mates must map
# -C               exclude chimeric fragments
# -s 2             reverse-stranded library (dUTP / TruSeq stranded)
# -t exon          feature type to count
# -g gene_id       group by gene_id attribute
# -F GTF           annotation format
# =============================================================================

# =============================================================================
# SELF-SUBMIT MODE
# =============================================================================

if [ -z "${SLURM_JOB_ID:-}" ]; then
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$SAMPLE_LOG_DIR"
    mkdir -p "$LOG_DIR"

    echo "Folders created:"
    echo "  $OUTPUT_DIR"
    echo "  $SAMPLE_LOG_DIR"
    echo "  $LOG_DIR"
    echo ""
    echo "Submitting featureCounts job to SLURM..."
    sbatch "$0"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$SAMPLE_LOG_DIR"
mkdir -p "$LOG_DIR"

# =============================================================================
# SCRATCH SETUP
# =============================================================================

SCRATCH="/scratch/${USER}_${SLURM_JOB_ID}_featurecounts"
SCRATCH_BAM="${SCRATCH}/bam"
SCRATCH_OUT="${SCRATCH}/counts"
SCRATCH_LOGS="${SCRATCH}/sample_logs"

mkdir -p "$SCRATCH_BAM" "$SCRATCH_OUT" "$SCRATCH_LOGS"

cleanup_on_error() {
    echo "ERROR: featureCounts failed."
    echo "Copying available output before exit..."
    cp -r "$SCRATCH_OUT"/  "$OUTPUT_DIR/"     2>/dev/null || true
    cp -r "$SCRATCH_LOGS"/ "$SAMPLE_LOG_DIR/" 2>/dev/null || true
    echo "Scratch kept for debugging: $SCRATCH"
}
trap cleanup_on_error ERR

echo "================================================="
echo "         FEATURECOUNTS — EUKARYOTIC              "
echo "================================================="
echo "Job ID           : ${SLURM_JOB_ID}"
echo "Node             : ${SLURMD_NODENAME:-NA}"
echo "Start Time       : $(date)"
echo "Project Dir      : $PROJECT_DIR"
echo "Input BAM Dir    : $INPUT_DIR"
echo "GTF File         : $GTF"
echo "Output Dir       : $OUTPUT_DIR"
echo "Scratch Dir      : $SCRATCH"
echo "Threads          : $THREADS"
echo "-------------------------------------------------"
echo "featureCounts parameters:"
echo "  -p               paired-end"
echo "  --countReadPairs count fragments not reads"
echo "  -B               both mates must map"
echo "  -C               exclude chimeric fragments"
echo "  -s 2             reverse-stranded"
echo "  -t exon          feature type"
echo "  -g gene_id       group by gene_id"
echo "  -F GTF           annotation format"
echo "================================================="

# =============================================================================
# ACTIVATE CONDA
# =============================================================================

echo "[$(date)] Activating conda environment: mapper"
source /local_conda/miniconda3/etc/profile.d/conda.sh
conda activate mapper

echo "[$(date)] featureCounts : $(which featureCounts)"
echo "[$(date)] Version       : $(featureCounts -v 2>&1 | head -n 1)"

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

# Input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: BAM input directory not found: $INPUT_DIR"
    rm -rf "$SCRATCH"
    exit 1
fi

# BAM files
BAM_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.bam" | wc -l)
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "ERROR: No BAM files found in $INPUT_DIR"
    echo "Contents:"
    ls -lh "$INPUT_DIR" || true
    rm -rf "$SCRATCH"
    exit 1
fi

echo "[$(date)] BAM files found: $BAM_COUNT"
find "$INPUT_DIR" -maxdepth 1 -type f -name "*.bam" | sort | \
    while read -r f; do echo "  $(basename "$f")  ($(du -sh "$f" | cut -f1))"; done

# GTF file
if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF file not found: $GTF"
    echo "Expected: $GTF"
    echo "Contents of genome dir:"
    ls -lh "${PROJECT_DIR}/genome/" || true
    rm -rf "$SCRATCH"
    exit 1
fi
echo "[$(date)] GTF found: $GTF  ($(du -sh "$GTF" | cut -f1))"

# Check / create BAM indexes
echo "[$(date)] Checking BAM index files (.bai)..."
while IFS= read -r bam; do
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ]; then
        echo "  No index for $(basename "$bam") — indexing now..."
        samtools index -@ "$THREADS" "$bam"
    fi
done < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.bam" | sort)
echo "[$(date)] BAM index check done."

# =============================================================================
# COPY BAM + BAI TO SCRATCH
# =============================================================================

echo "[$(date)] Copying BAM files to scratch..."
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.bam" -o -name "*.bam.bai" \) \
    -exec cp {} "$SCRATCH_BAM/" \;
echo "[$(date)] Copy done."

cd "$SCRATCH_BAM"

mapfile -t BAM_FILES < <(find . -maxdepth 1 -type f -name "*.bam" | sort)

if [ "${#BAM_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No BAM files found in scratch: $SCRATCH_BAM"
    rm -rf "$SCRATCH"
    exit 1
fi

echo "[$(date)] Samples to quantify: ${#BAM_FILES[@]}"
for f in "${BAM_FILES[@]}"; do echo "  $(basename "$f")"; done

# =============================================================================
# RUN featureCounts
# All samples in one call → one matrix, most efficient
# =============================================================================

RAW_OUT="${SCRATCH_OUT}/featureCounts.txt"
FC_LOG="${SCRATCH_LOGS}/featurecounts_run.log"

echo ""
echo "[$(date)] Running featureCounts..."

set +e

featureCounts \
    -T  "$THREADS" \
    -p  \
    --countReadPairs \
    -B  \
    -C  \
    -s  2 \
    -F  GTF \
    -t  exon \
    -g  gene_id \
    -a  "$GTF" \
    -o  "$RAW_OUT" \
    "${BAM_FILES[@]}" \
    > "$FC_LOG" 2>&1

FC_EC=$?
set -e

if [ "$FC_EC" -ne 0 ]; then
    echo "ERROR: featureCounts failed (exit code: $FC_EC)"
    echo "--- Last 80 lines of log ---"
    tail -n 80 "$FC_LOG" || true
    exit 1
fi

if [ ! -s "$RAW_OUT" ]; then
    echo "ERROR: Count matrix empty or not created: $RAW_OUT"
    exit 1
fi

echo "[$(date)] featureCounts finished."

# Print assignment summary
echo ""
echo "--- featureCounts assignment summary ---"
grep -A 30 "Summary" "$FC_LOG" || tail -n 30 "$FC_LOG" || true
echo "--- End ---"

# =============================================================================
# CONVERT COUNT MATRIX → CSV
# featureCounts outputs tab-separated with a comment header line (#)
# Steps:
#   1. Strip the leading comment line (starts with #)
#   2. Clean column headers: remove path prefix and .bam suffix
#   3. Convert tabs → commas
# =============================================================================

echo ""
echo "[$(date)] Converting count matrix to CSV..."

COUNTS_CSV="${SCRATCH_OUT}/featureCounts.csv"

awk 'NR==1 && /^#/ { next }          # skip comment line
     NR==1 {                          # header row — clean column names
         printf "%s", $1
         for (i=2; i<=6; i++) printf ",%s", $i
         for (i=7; i<=NF; i++) {
             col = $i
             sub(".*/", "", col)      # strip path
             sub("\\.bam$", "", col)  # strip .bam
             printf ",%s", col
         }
         print ""
         next
     }
     {                                # data rows — tab to comma
         printf "%s", $1
         for (i=2; i<=NF; i++) printf ",%s", $i
         print ""
     }' "$RAW_OUT" > "$COUNTS_CSV"

echo "[$(date)] featureCounts.csv written."

# =============================================================================
# CONVERT SUMMARY FILE → CSV
# .summary is tab-separated: Status | sample1.bam | sample2.bam ...
# Clean paths and convert to CSV
# =============================================================================

SUMMARY_TSV="${RAW_OUT}.summary"
SUMMARY_CSV="${SCRATCH_OUT}/featureCounts_summary.csv"

if [ -f "$SUMMARY_TSV" ]; then
    echo "[$(date)] Converting summary to CSV..."

    awk 'NR==1 {
             printf "%s", $1
             for (i=2; i<=NF; i++) {
                 col = $i
                 sub(".*/", "", col)
                 sub("\\.bam$", "", col)
                 printf ",%s", col
             }
             print ""
             next
         }
         {
             printf "%s", $1
             for (i=2; i<=NF; i++) printf ",%s", $i
             print ""
         }' "$SUMMARY_TSV" > "$SUMMARY_CSV"

    echo "[$(date)] featureCounts_summary.csv written."
else
    echo "WARNING: Summary file not found: $SUMMARY_TSV"
fi

# Preview outputs
echo ""
echo "--- Preview: featureCounts.csv (first 3 data rows) ---"
head -n 4 "$COUNTS_CSV" || true
echo ""
echo "--- Preview: featureCounts_summary.csv ---"
cat "$SUMMARY_CSV" || true
echo "--- End preview ---"

# =============================================================================
# COPY OUTPUT TO FINAL DESTINATION
# =============================================================================

echo ""
echo "[$(date)] Copying outputs to $OUTPUT_DIR ..."
cp -r "$SCRATCH_OUT"/  "$OUTPUT_DIR/"     2>/dev/null || true
cp -r "$SCRATCH_LOGS"/ "$SAMPLE_LOG_DIR/" 2>/dev/null || true
echo "[$(date)] Copy complete."

# =============================================================================
# FINAL SUMMARY
# =============================================================================

echo ""
echo "================================================="
echo "FEATURECOUNTS OUTPUT SUMMARY"
echo "================================================="
echo "Output directory      : $OUTPUT_DIR"
echo ""
echo "Files saved:"
find "$OUTPUT_DIR" -maxdepth 1 -type f | sort | \
    while read -r f; do
        echo "  $(basename "$f")  ($(du -sh "$f" | cut -f1))"
    done
echo ""
echo "Log files:"
find "$SAMPLE_LOG_DIR" -type f | sort | \
    while read -r f; do
        echo "  $(basename "$f")  ($(du -sh "$f" | cut -f1))"
    done
echo "================================================="

# =============================================================================
# CLEANUP
# =============================================================================

trap - ERR

echo "[$(date)] Removing scratch..."
rm -rf "$SCRATCH"

echo ""
echo "================================================="
echo "FEATURECOUNTS JOB COMPLETED SUCCESSFULLY"
echo "Completion Time  : $(date)"
echo "Counts CSV       : ${OUTPUT_DIR}/featureCounts.csv"
echo "Summary CSV      : ${OUTPUT_DIR}/featureCounts_summary.csv"
echo "================================================="
