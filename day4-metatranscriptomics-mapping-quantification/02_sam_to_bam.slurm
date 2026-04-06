#!/bin/bash
#SBATCH --job-name=sam_to_bam
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=3:00:00
#SBATCH --output=logs/sam_to_bam_%A_%a.out
#SBATCH --error=logs/sam_to_bam_%A_%a.err
#SBATCH --array=0-2          # adjust to number_of_samples - 1

# =============================================================
# Day 4 — SAM → sorted BAM → flagstat
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Steps:
#   1. SAM → BAM (samtools view)
#   2. Sort BAM by coordinate (samtools sort)
#   3. Index sorted BAM (samtools index)
#   4. Run flagstat for alignment diagnostics
#   5. Remove unsorted BAM and SAM to save disk space
#
# Input:  sam/${SAMPLE}.sam
# Output: bam/${SAMPLE}_sorted.bam
#         bam/${SAMPLE}_sorted.bam.bai
#         flagstat/${SAMPLE}_flagstat.txt
#
# Usage: sbatch 02_sam_to_bam.slurm
# =============================================================

module load samtools

# ── paths ─────────────────────────────────────────────────────
SAM_DIR="sam"
BAM_DIR="bam"
FLAGSTAT_DIR="flagstat"

mkdir -p "$BAM_DIR" "$FLAGSTAT_DIR" logs

# ── sample list ───────────────────────────────────────────────
SAMPLES=($(ls ${SAM_DIR}/*.sam | xargs -n1 basename | sed 's/.sam//'))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

SAM="${SAM_DIR}/${SAMPLE}.sam"
BAM="${BAM_DIR}/${SAMPLE}.bam"
SORTED="${BAM_DIR}/${SAMPLE}_sorted.bam"
FLAGSTAT="${FLAGSTAT_DIR}/${SAMPLE}_flagstat.txt"

echo "=============================================="
echo " SAM → sorted BAM · $SAMPLE"
echo " Start: $(date)"
echo "=============================================="

# ── check SAM exists ─────────────────────────────────────────
if [ ! -f "$SAM" ]; then
  echo "ERROR: SAM not found: $SAM"
  echo "Run 01_bowtie2_mapping.slurm first."
  exit 1
fi

# ── Step 1: SAM → BAM ────────────────────────────────────────
echo "[1/4] Converting SAM → BAM..."
samtools view -bS -@ "${SLURM_CPUS_PER_TASK}" "$SAM" > "$BAM"
echo "  ✓ BAM created"

# ── Step 2: Sort ──────────────────────────────────────────────
echo "[2/4] Sorting by coordinate..."
samtools sort \
  -@ "${SLURM_CPUS_PER_TASK}" \
  -m 3G \
  -o "$SORTED" \
  "$BAM"
echo "  ✓ Sorted BAM: $SORTED"

# ── Step 3: Index ─────────────────────────────────────────────
echo "[3/4] Indexing sorted BAM..."
samtools index "$SORTED"
echo "  ✓ Index: ${SORTED}.bai"

# ── Step 4: flagstat ──────────────────────────────────────────
echo "[4/4] Running flagstat..."
samtools flagstat "$SORTED" | tee "$FLAGSTAT"

# ── cleanup ───────────────────────────────────────────────────
rm -f "$BAM" "$SAM"
echo ""
echo "  ✓ Cleaned SAM and unsorted BAM"

# ── summary ───────────────────────────────────────────────────
MAPPED=$(grep "primary mapped" "$FLAGSTAT" | awk '{print $1}')
SINGLETONS=$(grep "singletons" "$FLAGSTAT" | awk '{print $1}')
PAIRED_OK=$(grep "properly paired" "$FLAGSTAT" | awk '{print $1}')

echo ""
echo "=============================================="
echo " $SAMPLE flagstat summary"
echo "----------------------------------------------"
echo " Primary mapped:   $MAPPED"
echo " Properly paired:  $PAIRED_OK  (0% expected — see Day 4 blog)"
echo " Singletons:       $SINGLETONS"
echo " Output: $SORTED"
echo " End: $(date)"
echo "=============================================="
