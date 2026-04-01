#!/bin/bash
#SBATCH --job-name=fastp_repair
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --output=logs/fastp_repair_%A_%a.out
#SBATCH --error=logs/fastp_repair_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 2 — fastp: Paired-end repair / re-sync (alternative to repair.sh)
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Two use cases for this script:
#
# USE CASE 1 — PREVENTION (run during initial trimming):
#   Use --unpaired1/--unpaired2 so fastp NEVER produces
#   mismatched output in the first place. Orphaned mates are
#   written to separate files, not mixed into paired output.
#   → Best practice: always run fastp this way.
#
# USE CASE 2 — RE-SYNC AFTER MISMATCH:
#   If you already have mismatched R1/R2 files (e.g. after
#   SortMeRNA or a failed pipeline step), fastp can re-sync
#   them using --reads_to_process and strict paired mode.
#   However, for severe mismatches repair.sh is more reliable.
#   Use fastp re-sync for minor mismatches (<1% difference).
#
# Choose MODE below: "prevent" or "resync"
# =============================================================

module load fastp

MODE="prevent"   # set to "prevent" or "resync"

# ── paths — edit these ────────────────────────────────────────
RAW_DIR="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/raw_data"
MRNA_DIR="mRNA_reads"          # used in resync mode
SYNCED_DIR="synced_reads"
UNPAIRED_DIR="unpaired_reads"
REPORT_DIR="qc_reports/fastp"

mkdir -p "$SYNCED_DIR" "$UNPAIRED_DIR" "$REPORT_DIR" logs

# ── sample list ───────────────────────────────────────────────
if [ "$MODE" = "prevent" ]; then
  SAMPLES=($(ls ${RAW_DIR}/*_R1.paired.fq | xargs -n1 basename | sed 's/_R1.paired.fq//'))
else
  SAMPLES=($(ls ${MRNA_DIR}/*_mRNA_fwd.fq | xargs -n1 basename | sed 's/_mRNA_fwd.fq//'))
fi

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=============================================="
echo " fastp paired-end sync · Mode: $MODE"
echo " Sample: $SAMPLE"
echo " Start: $(date)"
echo "=============================================="

# ──────────────────────────────────────────────────────────────
# MODE 1: PREVENTION — strict paired output during trimming
# ──────────────────────────────────────────────────────────────
if [ "$MODE" = "prevent" ]; then

  R1_IN="${RAW_DIR}/${SAMPLE}_R1.paired.fq"
  R2_IN="${RAW_DIR}/${SAMPLE}_R2.paired.fq"

  echo ""
  echo "Running fastp in strict paired mode (prevention)..."

  fastp \
    -i "$R1_IN" \
    -I "$R2_IN" \
    -o "${SYNCED_DIR}/${SAMPLE}_R1.trimmed.fq" \
    -O "${SYNCED_DIR}/${SAMPLE}_R2.trimmed.fq" \
    --unpaired1 "${UNPAIRED_DIR}/${SAMPLE}_R1.unpaired.fq" \
    --unpaired2 "${UNPAIRED_DIR}/${SAMPLE}_R2.unpaired.fq" \
    --failed_out "${UNPAIRED_DIR}/${SAMPLE}_failed.fq" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 50 \
    --correction \
    --thread "${SLURM_CPUS_PER_TASK}" \
    --html "${REPORT_DIR}/${SAMPLE}_fastp.html" \
    --json "${REPORT_DIR}/${SAMPLE}_fastp.json"

  # Note on --correction:
  # fastp can attempt to correct mismatched bases in overlapping
  # paired-end reads. Optional but useful for short-insert libraries.

fi

# ──────────────────────────────────────────────────────────────
# MODE 2: RESYNC — re-process mismatched files through fastp
# ──────────────────────────────────────────────────────────────
if [ "$MODE" = "resync" ]; then

  R1_IN="${MRNA_DIR}/${SAMPLE}_mRNA_fwd.fq"
  R2_IN="${MRNA_DIR}/${SAMPLE}_mRNA_rev.fq"

  R1_PRE=$(cat "$R1_IN" | wc -l | awk '{print int($1/4)}')
  R2_PRE=$(cat "$R2_IN" | wc -l | awk '{print int($1/4)}')

  echo ""
  echo "Pre-resync counts:"
  echo "  R1: $R1_PRE"
  echo "  R2: $R2_PRE"

  if [ "$R1_PRE" -eq "$R2_PRE" ]; then
    echo "  ✓ Already synchronized — no action needed."
    cp "$R1_IN" "${SYNCED_DIR}/${SAMPLE}_R1.synced.fq"
    cp "$R2_IN" "${SYNCED_DIR}/${SAMPLE}_R2.synced.fq"
    exit 0
  fi

  DIFF=$(( R1_PRE - R2_PRE ))
  echo "  ✗ Mismatch: $DIFF reads difference"

  # For large mismatches (>1%): use repair.sh instead (04_repair_bbtools.slurm)
  PCT_DIFF=$(echo "scale=2; ($DIFF < 0 ? -$DIFF : $DIFF) * 100 / $R1_PRE" | bc)
  echo "  Difference: ${PCT_DIFF}%"

  if (( $(echo "$PCT_DIFF > 1.0" | bc -l) )); then
    echo ""
    echo "  WARNING: mismatch >1% — recommend using repair.sh instead."
    echo "  (04_repair_bbtools.slurm handles large mismatches more reliably)"
    echo "  Continuing with fastp resync anyway — check output carefully."
  fi

  echo ""
  echo "Running fastp resync..."

  fastp \
    -i "$R1_IN" \
    -I "$R2_IN" \
    -o "${SYNCED_DIR}/${SAMPLE}_R1.synced.fq" \
    -O "${SYNCED_DIR}/${SAMPLE}_R2.synced.fq" \
    --unpaired1 "${UNPAIRED_DIR}/${SAMPLE}_R1.unpaired.fq" \
    --unpaired2 "${UNPAIRED_DIR}/${SAMPLE}_R2.unpaired.fq" \
    --detect_adapter_for_pe \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --thread "${SLURM_CPUS_PER_TASK}" \
    --html "${REPORT_DIR}/${SAMPLE}_fastp_resync.html" \
    --json "${REPORT_DIR}/${SAMPLE}_fastp_resync.json"

  # Note: disable trimming/filtering flags used here because
  # reads are already trimmed. We only want the pairing fix,
  # not additional quality changes.

fi

# ── post-run check: verify synchronization ────────────────────
if [ "$MODE" = "prevent" ]; then
  OUT_R1="${SYNCED_DIR}/${SAMPLE}_R1.trimmed.fq"
  OUT_R2="${SYNCED_DIR}/${SAMPLE}_R2.trimmed.fq"
else
  OUT_R1="${SYNCED_DIR}/${SAMPLE}_R1.synced.fq"
  OUT_R2="${SYNCED_DIR}/${SAMPLE}_R2.synced.fq"
fi

R1_POST=$(cat "$OUT_R1" | wc -l | awk '{print int($1/4)}')
R2_POST=$(cat "$OUT_R2" | wc -l | awk '{print int($1/4)}')
UNP_N=$(cat "${UNPAIRED_DIR}/${SAMPLE}_R1.unpaired.fq" | wc -l | awk '{print int($1/4)}')

echo ""
if [ "$R1_POST" -eq "$R2_POST" ]; then
  echo "  ✓ Output synchronized: $R1_POST pairs"
else
  echo "  ✗ ERROR: output still mismatched (R1=$R1_POST R2=$R2_POST)"
  echo "  → Run 04_repair_bbtools.slurm for rescue"
  exit 1
fi

echo ""
echo "----------------------------------------------"
echo " Sample:          $SAMPLE"
echo " Paired output:   $R1_POST reads each"
echo " Unpaired R1:     $UNP_N"
echo " Mode:            $MODE"
echo "----------------------------------------------"
echo " Paired reads  → ${SYNCED_DIR}/"
echo " Unpaired      → ${UNPAIRED_DIR}/"
echo " Reports       → ${REPORT_DIR}/"
echo " → Proceed to SortMeRNA (03_sortmerna.slurm)"
echo "   then check counts again before mapping"
echo "=============================================="
echo " fastp sync complete: $(date)"
echo "=============================================="
