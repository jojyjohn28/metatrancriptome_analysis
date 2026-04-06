#!/bin/bash
#SBATCH --job-name=bowtie2_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --output=logs/bowtie2_%A_%a.out
#SBATCH --error=logs/bowtie2_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 4 — Bowtie2 mapping (one sample per array task)
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Input:  mRNA_reads/${SAMPLE}_mRNA_fwd.fq
#         mRNA_reads/${SAMPLE}_mRNA_rev.fq
# Output: sam/${SAMPLE}.sam
#         logs/${SAMPLE}_bowtie2.log  (alignment statistics)
#
# Prerequisites:
#   ✓ Day 2: clean mRNA reads in mRNA_reads/ (R1 == R2 count)
#   ✓ Day 3: Bowtie2 index built (01_build_mag_index.slurm)
#
# Usage: sbatch 01_bowtie2_mapping.slurm
# =============================================================

module load bowtie2/2.5.4

# ── paths — edit these ────────────────────────────────────────
MRNA_DIR="mRNA_reads"
INDEX="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/MAGS/bowtie_index/master_MAGs_index"
SAM_DIR="sam"

mkdir -p "$SAM_DIR" logs

# ── sample list ───────────────────────────────────────────────
SAMPLES=($(ls ${MRNA_DIR}/*_mRNA_fwd.fq | xargs -n1 basename | sed 's/_mRNA_fwd.fq//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

R1="${MRNA_DIR}/${SAMPLE}_mRNA_fwd.fq"
R2="${MRNA_DIR}/${SAMPLE}_mRNA_rev.fq"
SAM="${SAM_DIR}/${SAMPLE}.sam"
LOG="logs/${SAMPLE}_bowtie2.log"

echo "=============================================="
echo " Bowtie2 mapping · $SAMPLE"
echo " Start: $(date)"
echo "=============================================="

# ── pre-flight: check reads are synced ───────────────────────
R1_N=$(cat "$R1" | wc -l | awk '{print int($1/4)}')
R2_N=$(cat "$R2" | wc -l | awk '{print int($1/4)}')

if [ "$R1_N" -ne "$R2_N" ]; then
  echo "ERROR: R1/R2 mismatch — R1=$R1_N R2=$R2_N"
  echo "Run Day 2 repair scripts before mapping."
  exit 1
fi
echo "✓ Input reads synchronized: $R1_N pairs"

# ── run Bowtie2 ───────────────────────────────────────────────
bowtie2 \
  -x "$INDEX" \
  -1 "$R1" \
  -2 "$R2" \
  -p "${SLURM_CPUS_PER_TASK}" \
  --no-unal \
  --very-sensitive \
  -S "$SAM" \
  2> "$LOG"

# Print alignment summary
echo ""
echo "Alignment summary:"
cat "$LOG"

ALIGN_RATE=$(grep "overall alignment rate" "$LOG" | awk '{print $1}')
echo ""
echo "=============================================="
echo " $SAMPLE → $ALIGN_RATE overall alignment rate"
echo " SAM: $SAM"
echo " Log: $LOG"
echo " End: $(date)"
echo "=============================================="
