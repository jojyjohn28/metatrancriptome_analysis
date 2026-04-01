#!/bin/bash
#SBATCH --job-name=repair_pairs
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=4:00:00
#SBATCH --output=logs/repair_%A_%a.out
#SBATCH --error=logs/repair_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 2 — repair.sh (BBTools): Fix mismatched paired-end reads
# Metatranscriptomics Series · Campbell Lab · 2026
#
# When to run this script:
#   After SortMeRNA (or any step), check if R1 and R2 read
#   counts are equal. If they differ, run this script BEFORE
#   proceeding to mapping. Mismatched pairs will cause:
#     - Bowtie2: "reads file 2 contains different number of reads"
#     - featureCounts: silent zeros in count table
#     - SAMtools flagstat: 0 properly paired
#
# How it works:
#   repair.sh loads all read IDs from both files into memory,
#   matches them by ID, and writes:
#     - out1/out2: synchronized matched pairs (use for mapping)
#     - outs:      singletons — reads whose mate was lost
#                  (save separately; can map as single-end)
#
# Memory: repair.sh loads IDs into RAM.
#   64 GB handles ~50M read pairs comfortably.
#   For larger files increase --mem and -Xmx proportionally.
#   Rule of thumb: -Xmx = --mem minus 4 GB
#
# Usage: sbatch 04_repair_pairs.slurm
# =============================================================

module load bbtools

# ── paths — edit these ────────────────────────────────────────
# Input: your mRNA reads after SortMeRNA
MRNA_DIR="mRNA_reads"

# Output directories
REPAIRED_DIR="repaired_reads"
SINGLETON_DIR="singleton_reads"

mkdir -p "$REPAIRED_DIR" "$SINGLETON_DIR" logs

# ── sample list ───────────────────────────────────────────────
SAMPLES=($(ls ${MRNA_DIR}/*_mRNA_fwd.fq | xargs -n1 basename | sed 's/_mRNA_fwd.fq//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

R1="${MRNA_DIR}/${SAMPLE}_mRNA_fwd.fq"
R2="${MRNA_DIR}/${SAMPLE}_mRNA_rev.fq"

echo "=============================================="
echo " repair.sh · Sample: $SAMPLE"
echo " Start: $(date)"
echo "=============================================="

# ── pre-repair check: show the mismatch ──────────────────────
R1_PRE=$(cat "$R1" | wc -l | awk '{print int($1/4)}')
R2_PRE=$(cat "$R2" | wc -l | awk '{print int($1/4)}')

echo ""
echo "Before repair:"
echo "  R1 reads: $R1_PRE"
echo "  R2 reads: $R2_PRE"

if [ "$R1_PRE" -eq "$R2_PRE" ]; then
  echo ""
  echo "  ✓ Files already synchronized — repair not needed."
  echo "  Copying to repaired directory as-is."
  cp "$R1" "${REPAIRED_DIR}/${SAMPLE}_R1.repaired.fq"
  cp "$R2" "${REPAIRED_DIR}/${SAMPLE}_R2.repaired.fq"
  touch "${SINGLETON_DIR}/${SAMPLE}_singletons.fq"
  echo "=============================================="
  echo " No repair needed for: $SAMPLE"
  echo " End: $(date)"
  echo "=============================================="
  exit 0
fi

DIFF=$(( R1_PRE - R2_PRE ))
echo "  ✗ MISMATCH detected: difference of $DIFF reads"
echo "  → Running repair.sh ..."
echo ""

# ── run repair.sh ─────────────────────────────────────────────
repair.sh \
  in1="$R1" \
  in2="$R2" \
  out1="${REPAIRED_DIR}/${SAMPLE}_R1.repaired.fq" \
  out2="${REPAIRED_DIR}/${SAMPLE}_R2.repaired.fq" \
  outs="${SINGLETON_DIR}/${SAMPLE}_singletons.fq" \
  repair=t \
  threads="${SLURM_CPUS_PER_TASK}" \
  -Xmx60g

# ── post-repair check: confirm synchronization ────────────────
R1_POST=$(cat "${REPAIRED_DIR}/${SAMPLE}_R1.repaired.fq" | wc -l | awk '{print int($1/4)}')
R2_POST=$(cat "${REPAIRED_DIR}/${SAMPLE}_R2.repaired.fq" | wc -l | awk '{print int($1/4)}')
SING_N=$(cat  "${SINGLETON_DIR}/${SAMPLE}_singletons.fq"  | wc -l | awk '{print int($1/4)}')

echo ""
echo "After repair:"
echo "  R1 repaired: $R1_POST"
echo "  R2 repaired: $R2_POST"
echo "  Singletons:  $SING_N  (saved to ${SINGLETON_DIR}/)"

if [ "$R1_POST" -eq "$R2_POST" ]; then
  echo ""
  echo "  ✓ Repair successful — R1 and R2 synchronized."
else
  echo ""
  echo "  ✗ ERROR: still mismatched after repair. Check input files."
  exit 1
fi

# ── summary ───────────────────────────────────────────────────
LOST=$(( R1_PRE - R1_POST ))

echo ""
echo "----------------------------------------------"
echo " Sample:          $SAMPLE"
echo " Original pairs:  $R1_PRE"
echo " Repaired pairs:  $R1_POST"
echo " Singletons:      $SING_N"
echo " Lost (neither):  $LOST"
echo "----------------------------------------------"
echo " Repaired reads → ${REPAIRED_DIR}/"
echo " Singletons     → ${SINGLETON_DIR}/"
echo " → Use repaired reads for Bowtie2 mapping (Day 4)"
echo " → Singletons can be mapped as single-end if needed"
echo "=============================================="
echo " repair.sh complete: $(date)"
echo "=============================================="
