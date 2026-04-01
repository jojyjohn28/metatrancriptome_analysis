#!/bin/bash
#SBATCH --job-name=fastqc_metatx
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 2 — FastQC: Raw read quality assessment
# Metatranscriptomics Series · Campbell Lab · 2026
# Run BEFORE trimming (raw) and AFTER preprocessing (clean)
# Usage: sbatch 01_fastqc.slurm
# =============================================================

module load fastqc

# ── directories ──────────────────────────────────────────────
RAW_DIR="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/raw_data"
OUT_RAW="fastqc_raw"
OUT_CLEAN="fastqc_clean"

mkdir -p "$OUT_RAW" "$OUT_CLEAN" logs

# ── sample list ───────────────────────────────────────────────
# List only R1 files; R2 is derived automatically
SAMPLES=($(ls ${RAW_DIR}/*_R1.paired.fq | xargs -n1 basename | sed 's/_R1.paired.fq//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "================================================"
echo "Running FastQC on: $SAMPLE"
echo "Start: $(date)"
echo "================================================"

# ── run on raw reads ──────────────────────────────────────────
fastqc \
  "${RAW_DIR}/${SAMPLE}_R1.paired.fq" \
  "${RAW_DIR}/${SAMPLE}_R2.paired.fq" \
  -o "$OUT_RAW" \
  -t "${SLURM_CPUS_PER_TASK}"

echo "Raw FastQC done. Start: $(date)"

# ── run on clean mRNA reads (after fastp + SortMeRNA) ─────────
# Uncomment once preprocessing is complete
# fastqc \
#   "mRNA_reads/${SAMPLE}_mRNA_fwd.fq" \
#   "mRNA_reads/${SAMPLE}_mRNA_rev.fq" \
#   -o "$OUT_CLEAN" \
#   -t "${SLURM_CPUS_PER_TASK}"

echo "================================================"
echo "FastQC complete for: $SAMPLE"
echo "End: $(date)"
echo "================================================"
