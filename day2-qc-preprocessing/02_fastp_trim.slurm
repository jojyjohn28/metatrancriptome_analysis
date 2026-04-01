#!/bin/bash
#SBATCH --job-name=fastp_metatx
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=3:00:00
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 2 — fastp: Adapter trimming & quality filtering
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Input:  raw paired-end FASTQ files
# Output: trimmed R1/R2 + orphaned reads + HTML/JSON report
#
# Usage: sbatch 02_fastp_trim.slurm
# =============================================================

module load fastp

# ── directories ──────────────────────────────────────────────
RAW_DIR="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/raw_data"
TRIM_DIR="trimmed_reads"
REPORT_DIR="fastp_reports"

mkdir -p "$TRIM_DIR" "$REPORT_DIR" logs

# ── sample list ───────────────────────────────────────────────
SAMPLES=($(ls ${RAW_DIR}/*_R1.paired.fq | xargs -n1 basename | sed 's/_R1.paired.fq//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "================================================"
echo "Running fastp on: $SAMPLE"
echo "Start: $(date)"
echo "================================================"

fastp \
  -i "${RAW_DIR}/${SAMPLE}_R1.paired.fq" \
  -I "${RAW_DIR}/${SAMPLE}_R2.paired.fq" \
  -o "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" \
  -O "${TRIM_DIR}/${SAMPLE}_R2.trimmed.fq" \
  --unpaired1 "${TRIM_DIR}/${SAMPLE}_unpaired.fq" \
  --unpaired2 "${TRIM_DIR}/${SAMPLE}_unpaired.fq" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --length_required 50 \
  --thread "${SLURM_CPUS_PER_TASK}" \
  --html "${REPORT_DIR}/${SAMPLE}_fastp.html" \
  --json "${REPORT_DIR}/${SAMPLE}_fastp.json"

# ── sanity check: confirm output exists and is non-empty ──────
if [[ -s "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" ]]; then
  RAW_COUNT=$(cat "${RAW_DIR}/${SAMPLE}_R1.paired.fq" | wc -l)
  TRIM_COUNT=$(cat "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fq" | wc -l)
  echo "Raw reads (R1):     $((RAW_COUNT / 4))"
  echo "Trimmed reads (R1): $((TRIM_COUNT / 4))"
else
  echo "ERROR: trimmed output not found or empty for $SAMPLE"
  exit 1
fi

echo "================================================"
echo "fastp complete for: $SAMPLE"
echo "End: $(date)"
echo "================================================"
