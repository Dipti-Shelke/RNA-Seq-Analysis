#!/bin/bash

set -euo pipefail

# ============================
# fastp paired-end trimming script
# ============================

OUTDIR="02_Trim_reads"
THREADS=20

mkdir -p "$OUTDIR"

LOGFILE="${OUTDIR}/fastp_run.log"

echo "fastp trimming started on: $(date)" > "$LOGFILE"
echo "Input directory: $(pwd)" >> "$LOGFILE"
echo "Output directory: $OUTDIR" >> "$LOGFILE"
echo "Threads: $THREADS" >> "$LOGFILE"
echo "----------------------------------------" >> "$LOGFILE"

if ! command -v fastp &> /dev/null; then
    echo "ERROR: fastp is not installed or not found in PATH" | tee -a "$LOGFILE"
    exit 1
fi

for R1 in *_R1_001.fastq.gz; do

    [[ -e "$R1" ]] || {
        echo "ERROR: No *_R1_001.fastq.gz files found." | tee -a "$LOGFILE"
        exit 1
    }

    echo "Processing $R1 ..."
    echo "Processing $R1" >> "$LOGFILE"

    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: Matching R2 file not found for $R1. Skipping." | tee -a "$LOGFILE"
        continue
    fi

    SAMPLE=$(echo "$R1" | sed -E 's/_S[0-9]+//g; s/_L[0-9]+//g; s/_R1_001\.fastq\.gz$//')

    OUT_R1="${OUTDIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    OUT_R2="${OUTDIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    HTML="${OUTDIR}/${SAMPLE}_fastp.html"
    JSON="${OUTDIR}/${SAMPLE}_fastp.json"

    echo "Sample name: $SAMPLE" >> "$LOGFILE"
    echo "R1 file: $R1" >> "$LOGFILE"
    echo "R2 file: $R2" >> "$LOGFILE"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT_R1" \
        -O "$OUT_R2" \
        --detect_adapter_for_pe \
        --cut_mean_quality 30 \
        --length_required 50 \
        -q 30 \
        --thread "$THREADS" \
        -h "$HTML" \
        -j "$JSON" \
        2>> "$LOGFILE"

    echo "Finished: $SAMPLE" >> "$LOGFILE"
    echo "----------------------------------------" >> "$LOGFILE"

done

echo "fastp trimming completed on: $(date)" >> "$LOGFILE"
echo "All files processed successfully."
