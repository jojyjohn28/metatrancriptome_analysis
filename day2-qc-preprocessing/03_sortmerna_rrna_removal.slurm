#!/bin/bash
#SBATCH --job-name=sortmerna_metatx
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --output=logs/sortmerna_%A_%a.out
#SBATCH --error=logs/sortmerna_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 2 — SortMeRNA: ribosomal RNA removal
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Input:  trimmed paired FASTQ (from fastp)
# Output: mRNA-enriched reads (--other) + rRNA reads (--aligned)
#
# Databases required (set SILVA_DB path below):
#   silva-bac-16s-id90.fasta
#   silva-arc-16s-id95.fasta
#   silva-euk-18s-id95.fasta
#   rfam-5s-database-id98.fasta
#
# Download databases:
#   wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
#   tar -xzf database.tar.gz
#
# Usage: sbatch 03_sortmerna_rrna_removal.slurm
# =============================================================

module load sortmerna

# ── paths — edit these ────────────────────────────────────────
TRIM_DIR="trimmed_reads"
RRNA_DIR="rRNA_reads"
MRNA_DIR="mRNA_reads"
TMP_DIR="sortmerna_tmp"

# Path to SortMeRNA database directory
SILVA_DB="/project/bcampb7/camplab/databases/sortmerna_db"

mkdir -p "$RRNA_DIR" "$MRNA_DIR" "$TMP_DIR" logs

# ── sample list ───────────────────────────────────────────────
SAMPLES=($(ls ${TRIM_DIR}/*_R1.trimmed.fq | xargs -n1 basename | sed 's/_R1.trimmed.fq//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# SortMeRNA requires its own working directory per sample
WORKDIR="${TMP_DIR}/${SAMPLE}"
mkdir -p "$WORKDIR"

echo "================================================"
echo "Running SortMeRNA on: $SAMPLE"
echo "Start: $(date)"
echo "================================================"

sortmerna \
  --ref "${SILVA_DB}/silva-bac-16s-id90.fasta" \
  --ref "${SILVA_DB}/silva-arc-16s-id95.fasta" \
  --ref "${SILVA_DB}/silva-euk-18s-id95.fasta" \
  --ref "${SILVA_DB}/rfam-5s-database-id98.fasta" \
  --reads "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" \
  --reads "${TRIM_DIR}/${SAMPLE}_R2.trimmed.fq" \
  --aligned  "${RRNA_DIR}/${SAMPLE}_rRNA" \
  --other    "${MRNA_DIR}/${SAMPLE}_mRNA" \
  --fastx \
  --paired-in \
  --out2 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --workdir "$WORKDIR"

# ── read count summary ────────────────────────────────────────
TRIM_COUNT=$(cat "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" | wc -l)
MRNA_COUNT=$(cat "${MRNA_DIR}/${SAMPLE}_mRNA_fwd.fq"   | wc -l)
RRNA_COUNT=$(cat "${RRNA_DIR}/${SAMPLE}_rRNA_fwd.fq"   | wc -l)

TRIM_PAIRS=$((TRIM_COUNT / 4))
MRNA_PAIRS=$((MRNA_COUNT / 4))
RRNA_PAIRS=$((RRNA_COUNT / 4))

echo "------------------------------------------------"
echo "Sample:           $SAMPLE"
echo "Trimmed pairs:    $TRIM_PAIRS"
echo "rRNA pairs:       $RRNA_PAIRS"
echo "mRNA pairs:       $MRNA_PAIRS"
echo "mRNA retention:   $(echo "scale=1; $MRNA_PAIRS * 100 / $TRIM_PAIRS" | bc)%"
echo "------------------------------------------------"

# ── clean up SortMeRNA tmp (large KV store) ───────────────────
rm -rf "$WORKDIR"

echo "================================================"
echo "SortMeRNA complete for: $SAMPLE"
echo "End: $(date)"
echo "================================================"
