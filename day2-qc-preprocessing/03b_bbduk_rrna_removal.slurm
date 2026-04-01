#!/bin/bash
#SBATCH --job-name=bbduk_metatx
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=3:00:00
#SBATCH --output=logs/bbduk_%A_%a.out
#SBATCH --error=logs/bbduk_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 2 — bbduk (BBTools): rRNA removal (alternative to SortMeRNA)
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Use this when:
#   - Processing many samples quickly (bbduk is faster)
#   - Your system is well-characterised (standard rRNA databases)
#
# Use SortMeRNA (03_sortmerna_rrna_removal.slurm) when:
#   - Working with novel/understudied environments
#   - You need maximum sensitivity for divergent rRNA
#
# Input:  trimmed paired FASTQ (from fastp)
# Output: mRNA-enriched R1/R2 reads
#
# Usage: sbatch 03b_bbduk_rrna_removal.slurm
# =============================================================

module load bbtools

# ── paths — edit these ────────────────────────────────────────
TRIM_DIR="trimmed_reads"
MRNA_DIR="mRNA_reads_bbduk"
RRNA_DIR="rRNA_reads_bbduk"

# ribokmers.fa.gz ships with BBTools installation
# Typical path on Palmetto:
RKMER="/software/external/bbtools/38.xx/resources/ribokmers.fa.gz"

mkdir -p "$MRNA_DIR" "$RRNA_DIR" logs

# ── sample list ───────────────────────────────────────────────
SAMPLES=($(ls ${TRIM_DIR}/*_R1.trimmed.fq | xargs -n1 basename | sed 's/_R1.trimmed.fq//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "================================================"
echo "Running bbduk rRNA removal on: $SAMPLE"
echo "Start: $(date)"
echo "================================================"

bbduk.sh \
  in1="${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" \
  in2="${TRIM_DIR}/${SAMPLE}_R2.trimmed.fq" \
  out1="${MRNA_DIR}/${SAMPLE}_R1.mRNA.fq" \
  out2="${MRNA_DIR}/${SAMPLE}_R2.mRNA.fq" \
  outm1="${RRNA_DIR}/${SAMPLE}_R1.rRNA.fq" \
  outm2="${RRNA_DIR}/${SAMPLE}_R2.rRNA.fq" \
  ref="$RKMER" \
  k=31 \
  ordered=t \
  threads="${SLURM_CPUS_PER_TASK}" \
  -Xmx28g

# ── read count summary ────────────────────────────────────────
TRIM_COUNT=$(cat "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" | wc -l)
MRNA_COUNT=$(cat "${MRNA_DIR}/${SAMPLE}_R1.mRNA.fq"    | wc -l)

TRIM_PAIRS=$((TRIM_COUNT / 4))
MRNA_PAIRS=$((MRNA_COUNT / 4))

echo "------------------------------------------------"
echo "Sample:           $SAMPLE"
echo "Trimmed pairs:    $TRIM_PAIRS"
echo "mRNA pairs:       $MRNA_PAIRS"
echo "mRNA retention:   $(echo "scale=1; $MRNA_PAIRS * 100 / $TRIM_PAIRS" | bc)%"
echo "------------------------------------------------"

echo "================================================"
echo "bbduk complete for: $SAMPLE"
echo "End: $(date)"
echo "================================================"
