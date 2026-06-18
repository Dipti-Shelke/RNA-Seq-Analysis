#!/bin/bash
# =============================================================================
# RNA-SEQ MASTER PIPELINE — SELF-SUBMITTING + SINGLE SLURM JOB
#
# HOW IT WORKS
# ────────────────────────────────────────────────────────────────────────────
#  bash rnaseq_pipeline.sh        ← you run this once from the terminal
#       │
#       ├─ NOT inside SLURM yet  → creates ALL output/log folders first,
#       │                           then calls  sbatch rnaseq_pipeline.sh
#       │                           and exits.
#       │
#       └─ INSIDE SLURM job     → runs all 5 steps sequentially.
#
# This two-phase design avoids the "log dir doesn't exist" SLURM crash.
# =============================================================================

# ── SLURM directives (only used when sbatch reads the file) ─────────────────
#SBATCH --job-name=rnaseq_pipeline
#SBATCH --partition=all
#SBATCH --nodelist=bionode1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --output=/data/projects/dipti/Test/final_test/rnaseq_output/logs/rnaseq_pipeline_%j.out
#SBATCH --error=/data/projects/dipti/Test/final_test/rnaseq_output/logs/rnaseq_pipeline_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=diptigulumbe@chromosomelabs.com

set -euo pipefail

# =============================================================================
# PATHS
# =============================================================================

WORK_DIR="/data/projects/dipti/Test/final_test"
RAW_FASTQ_DIR="${WORK_DIR}/00_Raw_reads"

GENOME_FASTA="/data/projects/dipti/Test/genome/GRCm39_genomic.fna"
GTF="/data/projects/dipti/Test/genome/GRCm39_genomic.gtf"
INDEX_DIR="/data/projects/dipti/Test/genome/mm_index"
HISAT2_INDEX="${INDEX_DIR}/mm39_index"
DB_DIR="/data/projects/dipti/Test/db"

OUTPUT_BASE="${WORK_DIR}/rnaseq_output"
LOG_DIR="${OUTPUT_BASE}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-30}"

# =============================================================================
# COLOUR HELPERS  (auto-disabled when not a terminal / inside SLURM log)
# =============================================================================

if [ -t 1 ]; then          # stdout is a real terminal
    C_RESET='\033[0m'
    C_BOLD='\033[1m'
    C_CYAN='\033[1;36m'
    C_GREEN='\033[1;32m'
    C_YELLOW='\033[1;33m'
    C_RED='\033[1;31m'
    C_BLUE='\033[1;34m'
else                        # inside SLURM log file — no escape codes
    C_RESET=''; C_BOLD=''; C_CYAN=''; C_GREEN=''
    C_YELLOW=''; C_RED=''; C_BLUE=''
fi

step_banner () {
    local num="$1" total="$2" name="$3"
    echo ""
    echo -e "${C_CYAN}${C_BOLD}╔══════════════════════════════════════════════════════════╗${C_RESET}"
    echo -e "${C_CYAN}${C_BOLD}  ▶  STEP ${num}/${total} — ${name}${C_RESET}"
    echo -e "${C_CYAN}${C_BOLD}     Started : $(date '+%Y-%m-%d %H:%M:%S')${C_RESET}"
    echo -e "${C_CYAN}${C_BOLD}╚══════════════════════════════════════════════════════════╝${C_RESET}"
}

step_done () {
    local num="$1" name="$2"
    echo -e "${C_GREEN}${C_BOLD}  ✔  STEP ${num} COMPLETE — ${name} — $(date '+%H:%M:%S')${C_RESET}"
    echo ""
}

info  () { echo -e "${C_BLUE}[INFO]${C_RESET}  $*"; }
warn  () { echo -e "${C_YELLOW}[WARN]${C_RESET}  $*"; }
error () { echo -e "${C_RED}[ERROR]${C_RESET} $*" >&2; }

# =============================================================================
# PHASE 1 — PRE-SLURM  (script is launched with plain  bash  by the user)
#
# $SLURM_JOB_ID is only set once SLURM actually starts running the job.
# If it's empty we are still on the login node → create folders, then sbatch.
# =============================================================================

if [ -z "${SLURM_JOB_ID:-}" ]; then

    echo -e "${C_BOLD}"
    echo "  ┌─────────────────────────────────────────────────────┐"
    echo "  │        RNA-SEQ PIPELINE  —  LAUNCH SEQUENCE         │"
    echo "  └─────────────────────────────────────────────────────┘"
    echo -e "${C_RESET}"

    # ── 1. Create every folder SLURM will need for its log file ──────────────
    info "Creating output directory tree …"

    mkdir -p "${LOG_DIR}"
    mkdir -p "${OUTPUT_BASE}/01_qc"
    mkdir -p "${OUTPUT_BASE}/02_fastp"
    mkdir -p "${OUTPUT_BASE}/03_sortmerna/rRNAFree"
    mkdir -p "${OUTPUT_BASE}/03_sortmerna/sample_logs"
    mkdir -p "${OUTPUT_BASE}/04_hisat2/sample_logs"
    mkdir -p "${OUTPUT_BASE}/05_featurecounts/sample_logs"

    echo -e "${C_GREEN}  ✔  Folders created:${C_RESET}"
    echo "       ${OUTPUT_BASE}/"
    echo "       ├── logs/"
    echo "       ├── 01_qc/"
    echo "       ├── 02_fastp/"
    echo "       ├── 03_sortmerna/  (rRNAFree/ + sample_logs/)"
    echo "       ├── 04_hisat2/     (sample_logs/)"
    echo "       └── 05_featurecounts/ (sample_logs/)"
    echo ""

    # ── 2. Quick sanity checks before wasting queue time ─────────────────────
    info "Checking raw FASTQ files …"
    FASTQ_COUNT=$(ls "${RAW_FASTQ_DIR}"/*.fastq.gz 2>/dev/null | wc -l)
    if [ "$FASTQ_COUNT" -eq 0 ]; then
        error "No *.fastq.gz files found in: ${RAW_FASTQ_DIR}"
        error "Check the path and try again."
        exit 1
    fi
    echo -e "${C_GREEN}  ✔  ${FASTQ_COUNT} FASTQ file(s) found in ${RAW_FASTQ_DIR}${C_RESET}"
    echo ""

    # ── 3. Submit self to SLURM ───────────────────────────────────────────────
    info "Submitting pipeline to SLURM …"
    JOB_OUT=$(sbatch "$0")          # $0 = this script
    JOB_ID=$(echo "${JOB_OUT}" | awk '{print $NF}')

    echo ""
    echo -e "${C_GREEN}${C_BOLD}  ✔  Job submitted successfully!${C_RESET}"
    echo ""
    echo "  Job ID   : ${JOB_ID}"
    echo "  Log file : ${LOG_DIR}/rnaseq_pipeline_${JOB_ID}.out"
    echo "  Err file : ${LOG_DIR}/rnaseq_pipeline_${JOB_ID}.err"
    echo ""
    echo "  Monitor with:"
    echo "    squeue -j ${JOB_ID}"
    echo "    tail -f ${LOG_DIR}/rnaseq_pipeline_${JOB_ID}.out"
    echo ""
    exit 0
fi

# =============================================================================
# PHASE 2 — INSIDE SLURM  (all steps run here)
# =============================================================================

# Ensure folders exist even if the script is re-submitted directly via sbatch
mkdir -p "${LOG_DIR}" \
         "${OUTPUT_BASE}/01_qc" \
         "${OUTPUT_BASE}/02_fastp" \
         "${OUTPUT_BASE}/03_sortmerna/rRNAFree" \
         "${OUTPUT_BASE}/03_sortmerna/sample_logs" \
         "${OUTPUT_BASE}/04_hisat2/sample_logs" \
         "${OUTPUT_BASE}/05_featurecounts/sample_logs"

# ── Conda ─────────────────────────────────────────────────────────────────────
source /local_conda/miniconda3/etc/profile.d/conda.sh

# ── Pipeline header ───────────────────────────────────────────────────────────
echo "=========================================================="
echo "         RNA-SEQ PIPELINE STARTED                        "
echo "=========================================================="
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Node       : ${SLURMD_NODENAME}"
echo "Start Time : $(date)"
echo "Work Dir   : ${WORK_DIR}"
echo "Output     : ${OUTPUT_BASE}"
echo "Threads    : ${THREADS}"
echo "=========================================================="

PIPELINE_START=$(date +%s)

# =============================================================================
# STEP 1 — FastQC
# =============================================================================

step_banner 1 5 "FastQC Quality Control"
info "→ Running FastQC on all raw reads"

conda activate trimQC
info "FastQC version: $(fastqc --version)"

FASTQ_COUNT=$(ls "${RAW_FASTQ_DIR}"/*.fastq.gz 2>/dev/null | wc -l)
if [ "$FASTQ_COUNT" -eq 0 ]; then
    error "No FASTQ files found in ${RAW_FASTQ_DIR}"
    exit 1
fi
info "Files to process: ${FASTQ_COUNT}"

SCRATCH_QC="/scratch/${USER}_${SLURM_JOB_ID}_qc"
mkdir -p "${SCRATCH_QC}/input" "${SCRATCH_QC}/output"

info "Copying FASTQ files to scratch …"
cp "${RAW_FASTQ_DIR}"/*.fastq.gz "${SCRATCH_QC}/input/"

info "Launching FastQC …"
fastqc \
    --threads "${THREADS}" \
    --outdir  "${SCRATCH_QC}/output" \
    --format  fastq \
    "${SCRATCH_QC}/input/"*.fastq.gz

cp -r "${SCRATCH_QC}/output/"* "${OUTPUT_BASE}/01_qc/"
rm -rf "${SCRATCH_QC}"

info "Results saved → ${OUTPUT_BASE}/01_qc/"
step_done 1 "FastQC"

# =============================================================================
# STEP 2 — fastp Adapter Trimming
# =============================================================================

step_banner 2 5 "fastp Adapter Trimming"
info "→ Trimming adapters and low-quality bases"

conda activate trimQC
info "fastp version: $(fastp --version 2>&1 | head -n1)"

cd "${RAW_FASTQ_DIR}"

SAMPLE_NUM=0
for R1 in *_R1_001.fastq.gz; do

    [[ -e "$R1" ]] || { error "No *_R1_001.fastq.gz files found."; exit 1; }

    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    if [ ! -f "$R2" ]; then
        warn "R2 not found for ${R1} — skipping."
        continue
    fi

    SAMPLE=$(echo "$R1" | sed -E 's/_S[0-9]+//g; s/_L[0-9]+//g; s/_R1_001\.fastq\.gz$//')
    SAMPLE_NUM=$(( SAMPLE_NUM + 1 ))
    info "  [Sample ${SAMPLE_NUM}] Trimming: ${SAMPLE}"

    fastp \
        -i  "$R1" \
        -I  "$R2" \
        -o  "${OUTPUT_BASE}/02_fastp/${SAMPLE}_R1_trimmed.fastq.gz" \
        -O  "${OUTPUT_BASE}/02_fastp/${SAMPLE}_R2_trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --cut_mean_quality 30 \
        --length_required 50 \
        -q  30 \
        --thread "${THREADS}" \
        -h  "${OUTPUT_BASE}/02_fastp/${SAMPLE}_fastp.html" \
        -j  "${OUTPUT_BASE}/02_fastp/${SAMPLE}_fastp.json" \
        2>> "${OUTPUT_BASE}/02_fastp/fastp_run.log"

    info "  [Sample ${SAMPLE_NUM}] Done: ${SAMPLE}"
done

cd "${WORK_DIR}"
info "All ${SAMPLE_NUM} sample(s) trimmed → ${OUTPUT_BASE}/02_fastp/"
step_done 2 "fastp"

# =============================================================================
# STEP 3 — SortMeRNA rRNA Removal
# =============================================================================

step_banner 3 5 "SortMeRNA rRNA Removal"
info "→ Removing ribosomal RNA reads"

conda activate sortmerna
info "SortMeRNA version: $(sortmerna --version 2>&1 | head -n1)"

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

for ref in "${REFS[@]}"; do
    [ -f "$ref" ] || { error "Missing DB file: ${ref}"; exit 1; }
done

SCRATCH_SMR="/scratch/${USER}_${SLURM_JOB_ID}_sortmerna"
mkdir -p "${SCRATCH_SMR}/input" "${SCRATCH_SMR}/rRNAFree" \
         "${SCRATCH_SMR}/logs"  "${SCRATCH_SMR}/workdir"

cp "${OUTPUT_BASE}/02_fastp/"*_trimmed.fastq.gz "${SCRATCH_SMR}/input/"

shopt -s nullglob
R1_FILES=( "${SCRATCH_SMR}/input/"*R1*.fastq.gz )
[ "${#R1_FILES[@]}" -gt 0 ] || { error "No R1 files in SortMeRNA input."; exit 1; }

REF_ARGS=()
for ref in "${REFS[@]}"; do REF_ARGS+=(--ref "$ref"); done

SAMPLE_NUM=0
for R1 in "${R1_FILES[@]}"; do

    R2="${R1/R1/R2}"
    [ -f "$R2" ] || { error "R2 not found for ${R1}"; exit 1; }

    SAMPLE=$(basename "$R1")
    SAMPLE="${SAMPLE%_R1_trimmed.fastq.gz}"
    SAMPLE_NUM=$(( SAMPLE_NUM + 1 ))
    info "  [Sample ${SAMPLE_NUM}] rRNA removal: ${SAMPLE}"

    SAMPLE_WORKDIR="${SCRATCH_SMR}/workdir/${SAMPLE}"
    mkdir -p "${SAMPLE_WORKDIR}"

    sortmerna \
        "${REF_ARGS[@]}" \
        --reads "$R1" \
        --reads "$R2" \
        --paired_in \
        --fastx \
        --out2 \
        --other "${SCRATCH_SMR}/rRNAFree/${SAMPLE}_rRNAFree" \
        --threads "${THREADS}" \
        --workdir "${SAMPLE_WORKDIR}" \
        > "${SCRATCH_SMR}/logs/sortmerna_${SAMPLE}.log" 2>&1

    info "  [Sample ${SAMPLE_NUM}] Done: ${SAMPLE}"
done

cp -r "${SCRATCH_SMR}/rRNAFree/"* "${OUTPUT_BASE}/03_sortmerna/rRNAFree/"
cp -r "${SCRATCH_SMR}/logs/"*     "${OUTPUT_BASE}/03_sortmerna/sample_logs/"
rm -rf "${SCRATCH_SMR}"

info "All ${SAMPLE_NUM} sample(s) processed → ${OUTPUT_BASE}/03_sortmerna/rRNAFree/"
step_done 3 "SortMeRNA"

# =============================================================================
# STEP 4 — HISAT2 Alignment
# =============================================================================

step_banner 4 5 "HISAT2 Alignment"
info "→ Aligning rRNA-free reads to the mouse genome (GRCm39)"

conda activate mapper
info "HISAT2  : $(hisat2 --version | head -n1)"
info "samtools: $(samtools --version | head -n1)"

# Check / build index
if ls "${HISAT2_INDEX}".*.ht2 >/dev/null 2>&1 || \
   ls "${HISAT2_INDEX}".*.ht2l >/dev/null 2>&1; then
    info "Existing HISAT2 index found — skipping build."
else
    info "No index found — building from ${GENOME_FASTA} …"
    [ -f "${GENOME_FASTA}" ] || { error "Genome FASTA not found: ${GENOME_FASTA}"; exit 1; }
    hisat2-build -p "${THREADS}" "${GENOME_FASTA}" "${HISAT2_INDEX}" \
        > "${LOG_DIR}/hisat2_build_${SLURM_JOB_ID}.log" 2>&1
    info "Index build complete."
fi

SCRATCH_H2="/scratch/${USER}_${SLURM_JOB_ID}_hisat2"
mkdir -p "${SCRATCH_H2}/input" "${SCRATCH_H2}/bam" \
         "${SCRATCH_H2}/logs"  "${SCRATCH_H2}/tmp"

# Copy ALL fastq/fq files from SortMeRNA output (handles .fq.gz and .fastq.gz)
shopt -s nullglob
SMR_FILES=( "${OUTPUT_BASE}/03_sortmerna/rRNAFree/"*.fastq.gz
            "${OUTPUT_BASE}/03_sortmerna/rRNAFree/"*.fq.gz )
[ "${#SMR_FILES[@]}" -gt 0 ] || { error "No reads found in 03_sortmerna/rRNAFree/ — check SortMeRNA step."; exit 1; }
info "Copying ${#SMR_FILES[@]} file(s) from SortMeRNA output to scratch …"
cp "${SMR_FILES[@]}" "${SCRATCH_H2}/input/"

R1_FILES=( "${SCRATCH_H2}/input/"*_rRNAFree_fwd.fastq.gz
           "${SCRATCH_H2}/input/"*_rRNAFree_fwd.fq.gz
           "${SCRATCH_H2}/input/"*_1.non_rRNA.fastq.gz
           "${SCRATCH_H2}/input/"*_1.non_rRNA.fq.gz
           "${SCRATCH_H2}/input/"*_R1_trimmed.fastq.gz
           "${SCRATCH_H2}/input/"*_R1_trimmed.fq.gz )

[ "${#R1_FILES[@]}" -gt 0 ] || { error "No R1 files found in HISAT2 input."; exit 1; }

SAMPLE_NUM=0
for R1 in "${R1_FILES[@]}"; do

    if   [[ "$R1" == *_rRNAFree_fwd.fastq.gz ]]; then
        R2="${R1/_rRNAFree_fwd.fastq.gz/_rRNAFree_rev.fastq.gz}"
        SAMPLE=$(basename "${R1%_rRNAFree_fwd.fastq.gz}")
    elif [[ "$R1" == *_rRNAFree_fwd.fq.gz ]]; then
        R2="${R1/_rRNAFree_fwd.fq.gz/_rRNAFree_rev.fq.gz}"
        SAMPLE=$(basename "${R1%_rRNAFree_fwd.fq.gz}")
    elif [[ "$R1" == *_1.non_rRNA.fastq.gz ]]; then
        R2="${R1/_1.non_rRNA.fastq.gz/_2.non_rRNA.fastq.gz}"
        SAMPLE=$(basename "${R1%_1.non_rRNA.fastq.gz}")
    elif [[ "$R1" == *_R1_trimmed.fastq.gz ]]; then
        R2="${R1/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}"
        SAMPLE=$(basename "${R1%_R1_trimmed.fastq.gz}")
    fi

    [ -f "$R2" ] || { error "R2 not found for ${R1}"; exit 1; }

    SAMPLE_NUM=$(( SAMPLE_NUM + 1 ))
    BAM_OUT="${SCRATCH_H2}/bam/${SAMPLE}.bam"
    info "  [Sample ${SAMPLE_NUM}] Aligning: ${SAMPLE}"

    set +e
    hisat2 \
        -x "${HISAT2_INDEX}" \
        -1 "$R1" \
        -2 "$R2" \
        -p "${THREADS}" \
        2> "${SCRATCH_H2}/logs/hisat2_${SAMPLE}.log" \
    | samtools view -@ "${THREADS}" -bS - \
    | samtools sort -@ "${THREADS}" \
                    -T "${SCRATCH_H2}/tmp/${SAMPLE}_sort" \
                    -o "$BAM_OUT" -

    PIPE_STATUS=( "${PIPESTATUS[@]}" )
    set -e

    if [ "${PIPE_STATUS[0]}" -ne 0 ] || \
       [ "${PIPE_STATUS[1]:-0}" -ne 0 ] || \
       [ "${PIPE_STATUS[2]:-0}" -ne 0 ]; then
        error "HISAT2 pipeline failed for ${SAMPLE}"
        echo "  hisat2: ${PIPE_STATUS[0]}  view: ${PIPE_STATUS[1]:-0}  sort: ${PIPE_STATUS[2]:-0}"
        tail -n 30 "${SCRATCH_H2}/logs/hisat2_${SAMPLE}.log" || true
        exit 1
    fi

    samtools index -@ "${THREADS}" "$BAM_OUT"

    # Print alignment summary to terminal / log
    echo "--- Alignment summary: ${SAMPLE} ---"
    cat "${SCRATCH_H2}/logs/hisat2_${SAMPLE}.log"
    info "  [Sample ${SAMPLE_NUM}] Done: ${SAMPLE}"
done

cp -r "${SCRATCH_H2}/bam/"*  "${OUTPUT_BASE}/04_hisat2/"
cp -r "${SCRATCH_H2}/logs/"* "${OUTPUT_BASE}/04_hisat2/sample_logs/"
rm -rf "${SCRATCH_H2}"

info "All ${SAMPLE_NUM} sample(s) aligned → ${OUTPUT_BASE}/04_hisat2/"
step_done 4 "HISAT2"

# =============================================================================
# STEP 5 — featureCounts Quantification
# =============================================================================

step_banner 5 5 "featureCounts Gene Quantification"
info "→ Counting reads per gene from BAM files"

conda activate mapper
info "featureCounts: $(featureCounts -v 2>&1 | head -n1)"

[ -f "${GTF}" ] || { error "GTF not found: ${GTF}"; exit 1; }

BAM_COUNT=$(find "${OUTPUT_BASE}/04_hisat2" -maxdepth 1 -name "*.bam" | wc -l)
[ "$BAM_COUNT" -gt 0 ] || { error "No BAM files found in ${OUTPUT_BASE}/04_hisat2"; exit 1; }
info "BAM files to quantify: ${BAM_COUNT}"

SCRATCH_FC="/scratch/${USER}_${SLURM_JOB_ID}_featurecounts"
mkdir -p "${SCRATCH_FC}/bam" "${SCRATCH_FC}/counts" "${SCRATCH_FC}/logs"

find "${OUTPUT_BASE}/04_hisat2" -maxdepth 1 \
    \( -name "*.bam" -o -name "*.bam.bai" \) \
    -exec cp {} "${SCRATCH_FC}/bam/" \;

cd "${SCRATCH_FC}/bam"
mapfile -t BAM_FILES < <(find . -maxdepth 1 -name "*.bam" | sort)

RAW_OUT="${SCRATCH_FC}/counts/featureCounts.txt"
FC_LOG="${SCRATCH_FC}/logs/featurecounts_run.log"

info "Running featureCounts on ${#BAM_FILES[@]} BAM file(s) …"

set +e
featureCounts \
    -T  "${THREADS}" \
    -p  \
    --countReadPairs \
    -B  \
    -C  \
    -s  2 \
    -F  GTF \
    -t  exon \
    -g  gene_id \
    -a  "${GTF}" \
    -o  "${RAW_OUT}" \
    "${BAM_FILES[@]}" \
    > "$FC_LOG" 2>&1
FC_EC=$?
set -e

[ "$FC_EC" -eq 0 ] || { error "featureCounts failed."; tail -n 40 "$FC_LOG"; exit 1; }
[ -s "${RAW_OUT}" ]  || { error "Count matrix is empty."; exit 1; }

echo "--- featureCounts summary ---"
grep -A 30 "Summary" "$FC_LOG" || tail -n 30 "$FC_LOG"

# ── Convert raw counts to CSV ─────────────────────────────────────────────────
info "Converting count matrix to CSV …"
COUNTS_CSV="${SCRATCH_FC}/counts/featureCounts.csv"
awk 'NR==1 && /^#/ { next }
     NR==1 {
         printf "%s", $1
         for (i=2; i<=6; i++) printf ",%s", $i
         for (i=7; i<=NF; i++) {
             col = $i; sub(".*/","",col); sub("\\.bam$","",col)
             printf ",%s", col
         }
         print ""; next
     }
     { printf "%s", $1; for (i=2; i<=NF; i++) printf ",%s", $i; print "" }
    ' "${RAW_OUT}" > "${COUNTS_CSV}"

# ── Convert summary to CSV ────────────────────────────────────────────────────
SUMMARY_TSV="${RAW_OUT}.summary"
SUMMARY_CSV="${SCRATCH_FC}/counts/featureCounts_summary.csv"
if [ -f "${SUMMARY_TSV}" ]; then
    awk 'NR==1 {
             printf "%s", $1
             for (i=2; i<=NF; i++) {
                 col=$i; sub(".*/","",col); sub("\\.bam$","",col)
                 printf ",%s", col
             }
             print ""; next
         }
         { printf "%s", $1; for (i=2; i<=NF; i++) printf ",%s", $i; print "" }
        ' "${SUMMARY_TSV}" > "${SUMMARY_CSV}"
fi

cp -r "${SCRATCH_FC}/counts/"* "${OUTPUT_BASE}/05_featurecounts/"
cp -r "${SCRATCH_FC}/logs/"*   "${OUTPUT_BASE}/05_featurecounts/sample_logs/"
rm -rf "${SCRATCH_FC}"

info "Results saved → ${OUTPUT_BASE}/05_featurecounts/"
step_done 5 "featureCounts"

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
HOURS=$(( ELAPSED / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo -e "${C_GREEN}${C_BOLD}"
echo "  ╔══════════════════════════════════════════════════════════╗"
echo "  ║      RNA-SEQ PIPELINE COMPLETED SUCCESSFULLY  ✔         ║"
echo "  ╚══════════════════════════════════════════════════════════╝"
echo -e "${C_RESET}"
echo "Job ID          : ${SLURM_JOB_ID}"
echo "End Time        : $(date)"
echo "Total Runtime   : ${HOURS}h ${MINS}m ${SECS}s"
echo ""
echo "Output summary:"
echo "  01_qc           : $(find "${OUTPUT_BASE}/01_qc"           -type f | wc -l) files"
echo "  02_fastp        : $(find "${OUTPUT_BASE}/02_fastp"         -type f | wc -l) files"
echo "  03_sortmerna    : $(find "${OUTPUT_BASE}/03_sortmerna"     -type f | wc -l) files"
echo "  04_hisat2       : $(find "${OUTPUT_BASE}/04_hisat2"        -maxdepth 1 -name '*.bam' | wc -l) BAM files"
echo "  05_featurecounts: featureCounts.csv + featureCounts_summary.csv"
echo ""
echo "All results in  : ${OUTPUT_BASE}/"
echo "=========================================================="
