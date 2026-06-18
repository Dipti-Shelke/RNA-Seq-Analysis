# RNA-Seq-Analysis
Robust end-to-end RNA-Seq pipeline from raw reads to publication-ready figures, covering QC, trimming, filtering, rRNA removal, alignment, sorting, and featureCounts — delivering clean gene count matrices, differential expression results, and easy-to-interpret biological visualizations built from 5 years of multi-organism transcriptomics expertise


<h1 align="center">🧬 RNA-Seq Analysis Pipeline</h1>

<p align="center">
  <img src="https://img.shields.io/badge/Tools-STAR%20|%20DESeq2%20|%20edgeR-blue?style=flat-square"/>
  <img src="https://img.shields.io/badge/Language-R%20|%20Python%20|%20Bash-green?style=flat-square"/>
  <img src="https://img.shields.io/badge/HPC-SLURM-orange?style=flat-square"/>
  <img src="https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square"/>
</p>

<p align="center">
  <b>End-to-end RNA-Seq pipeline built from 40+ real-world projects across human, mouse, plant & custom organisms.</b><br/>
  Covers QC → Alignment → Quantification → Differential Expression → Pathway Enrichment.
</p>

---

## 🔬 What This Pipeline Does

- Performs quality control and adapter trimming on raw FASTQ files
- Aligns reads to reference genomes using STAR (splice-aware)
- Quantifies gene-level expression using featureCounts
- Identifies differentially expressed genes with DESeq2 and edgeR
- Runs GO and KEGG pathway enrichment using clusterProfiler
- Produces publication-ready plots: volcano plots, heatmaps, PCA

---

## 🗺️ Pipeline Workflow

```
Raw FASTQ
  ↓
FastQC + MultiQC        → Quality Reports
  ↓
Trimmomatic             → Adapter & Quality Trimming
  ↓
STAR                    → Genome Alignment (.bam)
  ↓
featureCounts           → Gene-level Read Counts
  ↓
DESeq2 / edgeR          → Differential Expression
  ↓
clusterProfiler         → GO & KEGG Pathway Enrichment
  ↓
ggplot2 / pheatmap      → Volcano Plots | Heatmaps | PCA
```

---

## 🛠️ Tools Used

| Step | Tool |
|------|------|
| Quality Control | FastQC, MultiQC |
| Trimming | Trimmomatic |
| Alignment | STAR v2.7.10 |
| Quantification | featureCounts (Subread) |
| DE Analysis | DESeq2, edgeR |
| Pathway Enrichment | clusterProfiler, fgsea |
| Visualization | ggplot2, pheatmap, EnhancedVolcano |
| Cluster | SLURM (HPC job arrays) |

---

## 📁 Folder Structure

```
RNA-Seq-Analysis/
├── scripts/
│   ├── 01_QC_fastqc.sh
│   ├── 02_trimming_trimmomatic.sh
│   ├── 03_alignment_STAR.sh
│   ├── 04_featureCounts.sh
│   ├── 05_DESeq2_analysis.R
│   └── 06_pathway_enrichment.R
├── envs/
│   └── rnaseq_env.yml
├── config/
│   └── params.yaml
└── results/
    ├── volcano_plot.png
    ├── heatmap.png
    └── PCA.png
```

---

## ⚡ Quick Start

```bash
# 1. Clone
git clone https://github.com/Dipti-Shelke/RNA-Seq-Analysis.git
cd RNA-Seq-Analysis

# 2. Setup environment
conda env create -f envs/rnaseq_env.yml
conda activate rnaseq

# 3. Run steps
bash scripts/01_QC_fastqc.sh /path/to/fastq/
bash scripts/02_trimming_trimmomatic.sh
bash scripts/03_alignment_STAR.sh /path/to/genome/ /path/to/gtf/
bash scripts/04_featureCounts.sh
Rscript scripts/05_DESeq2_analysis.R
```

---

## 📊 Sample Outputs

- ✅ DE gene table with log2FC, p-value, and adjusted p-value
- ✅ Volcano plot highlighting significant genes
- ✅ Heatmap of top DEGs across conditions
- ✅ PCA plot for sample QC and clustering
- ✅ GO/KEGG enrichment bar and dot plots

---

## 💡 Highlights

- 40+ RNA-Seq projects across human, mouse, plant, and microbial organisms
- Experience with paired-end, stranded, and single-end libraries
- Multi-factorial and time-series experimental designs
- HPC/SLURM execution with job arrays for large cohorts

---

## 👩‍🔬 Author

**Dipti Shelke** | Bioinformatics Scientist

[![GitHub](https://img.shields.io/badge/GitHub-Dipti--Shelke-black?logo=github&style=flat-square)](https://github.com/Dipti-Shelke)
