#!/bin/bash
#SBATCH --job-name=check_pairs
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=0:30:00
#SBATCH --output=logs/check_pairs_%j.out
#SBATCH --error=logs/check_pairs_%j.err

# =============================================================
# Day 2 — Paired-end read count validation (run before mapping)
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Run this after EVERY preprocessing step that writes FASTQ:
#   - After fastp
#   - After SortMeRNA / bbduk
#   - After repair.sh
#   - Before Bowtie2 (Day 4)
#
# Output: a clear PASS/FAIL table for all samples.
# Any FAIL must be resolved before proceeding to mapping.
#
# Usage: sbatch 00_check_pair_counts.slurm
#   or interactively: bash 00_check_pair_counts.slurm
# =============================================================

# ── set directory to check ────────────────────────────────────
# Change this to the directory containing your current mRNA reads
CHECK_DIR="mRNA_reads"   # or "repaired_reads", "synced_reads", etc.

# ── file suffix pattern ───────────────────────────────────────
# Adjust to match your filename convention
R1_SUFFIX="_mRNA_fwd.fq"
R2_SUFFIX="_mRNA_rev.fq"

# ── report file ───────────────────────────────────────────────
REPORT="logs/pair_count_report_$(date +%Y%m%d_%H%M%S).txt"
mkdir -p logs

echo "=============================================="  | tee "$REPORT"
echo " Paired-end read count validation"              | tee -a "$REPORT"
echo " Directory: $CHECK_DIR"                         | tee -a "$REPORT"
echo " Run: $(date)"                                  | tee -a "$REPORT"
echo "=============================================="  | tee -a "$REPORT"
printf " %-30s %10s %10s %8s\n" "Sample" "R1 reads" "R2 reads" "Status" | tee -a "$REPORT"
echo "----------------------------------------------"  | tee -a "$REPORT"

PASS=0
FAIL=0
FAIL_LIST=()

for R1 in ${CHECK_DIR}/*${R1_SUFFIX}; do
  [ -f "$R1" ] || continue

  SAMPLE=$(basename "$R1" "$R1_SUFFIX")
  R2="${CHECK_DIR}/${SAMPLE}${R2_SUFFIX}"

  # Check R2 exists
  if [ ! -f "$R2" ]; then
    printf " %-30s %10s %10s %8s\n" "$SAMPLE" "N/A" "MISSING" "FAIL" | tee -a "$REPORT"
    FAIL=$((FAIL + 1))
    FAIL_LIST+=("$SAMPLE (R2 missing)")
    continue
  fi

  R1_N=$(cat "$R1" | wc -l | awk '{print int($1/4)}')
  R2_N=$(cat "$R2" | wc -l | awk '{print int($1/4)}')

  if [ "$R1_N" -eq "$R2_N" ]; then
    printf " %-30s %10s %10s %8s\n" "$SAMPLE" "$R1_N" "$R2_N" "✓ PASS" | tee -a "$REPORT"
    PASS=$((PASS + 1))
  else
    DIFF=$(( R1_N - R2_N ))
    if [ "$DIFF" -lt 0 ]; then DIFF=$(( -DIFF )); fi
    PCT=$(echo "scale=2; $DIFF * 100 / $R1_N" | bc)
    printf " %-30s %10s %10s %8s\n" "$SAMPLE" "$R1_N" "$R2_N" "✗ FAIL" | tee -a "$REPORT"
    echo "   → difference: $DIFF reads (${PCT}%)" | tee -a "$REPORT"
    FAIL=$((FAIL + 1))
    FAIL_LIST+=("$SAMPLE (diff=${DIFF}, ${PCT}%)")
  fi
done

echo ""                                               | tee -a "$REPORT"
echo "=============================================="  | tee -a "$REPORT"
echo " PASS: $PASS samples"                           | tee -a "$REPORT"
echo " FAIL: $FAIL samples"                           | tee -a "$REPORT"

if [ "${#FAIL_LIST[@]}" -gt 0 ]; then
  echo ""                                             | tee -a "$REPORT"
  echo " Failed samples:"                             | tee -a "$REPORT"
  for F in "${FAIL_LIST[@]}"; do
    echo "   ✗ $F"                                    | tee -a "$REPORT"
  done
  echo ""                                             | tee -a "$REPORT"
  echo " Action required:"                            | tee -a "$REPORT"
  echo "   Difference <1%  → 04b_fastp_paired_sync.slurm"   | tee -a "$REPORT"
  echo "   Difference ≥1%  → 04_repair_bbtools.slurm"       | tee -a "$REPORT"
  echo "   DO NOT proceed to mapping (Day 4) with FAIL samples" | tee -a "$REPORT"
  echo "=============================================="  | tee -a "$REPORT"
  exit 1
else
  echo ""                                             | tee -a "$REPORT"
  echo " ✓ All samples synchronized — safe to proceed"| tee -a "$REPORT"
  echo " → Next: Day 4 mapping (bowtie_mapping.slurm)"| tee -a "$REPORT"
  echo "=============================================="  | tee -a "$REPORT"
  exit 0
fi
