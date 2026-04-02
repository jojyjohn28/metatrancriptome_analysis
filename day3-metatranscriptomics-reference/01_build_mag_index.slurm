#!/bin/bash
#SBATCH --job-name=build_mag_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --output=logs/build_mag_index_%j.out
#SBATCH --error=logs/build_mag_index_%j.err

# =============================================================
# Day 3 — Build Bowtie2 index from MAGs (master_mags.fa)
# Metatranscriptomics Series · Campbell Lab · 2026
#
# What this script does:
#   1. Checks that all MAG FASTA files are present
#   2. Concatenates MAGs into a single master reference
#   3. Reports contig count per MAG and total
#   4. Builds Bowtie2 index (master_MAGs_index)
#   5. Verifies all 6 index files were created
#
# Run this ONCE before mapping (Day 4).
# Reuse the same index for all samples.
#
# Memory: for 4 MAGs (~few MB) — 8 GB is plenty.
#         For large references (>4 GB) use --large-index flag
#         and increase --mem to 200+ GB.
#
# Usage: sbatch 01_build_mag_index.slurm
# =============================================================

module load bowtie2/2.5.4

# ── paths — edit these ────────────────────────────────────────
MAG_DIR="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/MAGS"
INDEX_DIR="${MAG_DIR}/bowtie_index"
MASTER_FA="${MAG_DIR}/master_mags.fa"
INDEX_PREFIX="${INDEX_DIR}/master_MAGs_index"

mkdir -p "$INDEX_DIR" logs

echo "=============================================="
echo " Day 3 — Build Bowtie2 index from MAGs"
echo " Start: $(date)"
echo "=============================================="

# ── Step 1: verify MAG files exist ───────────────────────────
echo ""
echo "[1/4] Checking MAG files in: $MAG_DIR"

MAG_FILES=($(ls ${MAG_DIR}/*.fasta 2>/dev/null || ls ${MAG_DIR}/*.fa 2>/dev/null))

if [ ${#MAG_FILES[@]} -eq 0 ]; then
  echo "ERROR: No .fasta or .fa files found in $MAG_DIR"
  echo "       Expected: MAG1.fasta MAG2.fasta MAG3.fasta MAG4.fasta"
  exit 1
fi

echo "Found ${#MAG_FILES[@]} MAG file(s):"
TOTAL_CONTIGS=0
for F in "${MAG_FILES[@]}"; do
  N=$(grep -c ">" "$F")
  TOTAL_CONTIGS=$((TOTAL_CONTIGS + N))
  echo "  $(basename $F): $N contigs"
done
echo "Total contigs: $TOTAL_CONTIGS"

# ── Step 2: concatenate MAGs ──────────────────────────────────
echo ""
echo "[2/4] Concatenating MAGs → master_mags.fa"

cat "${MAG_FILES[@]}" > "$MASTER_FA"

# Verify count matches
MASTER_N=$(grep -c ">" "$MASTER_FA")
echo "Contigs in master_mags.fa: $MASTER_N"

if [ "$MASTER_N" -ne "$TOTAL_CONTIGS" ]; then
  echo "ERROR: contig count mismatch after concatenation"
  echo "       Expected $TOTAL_CONTIGS, found $MASTER_N"
  exit 1
fi
echo "✓ Concatenation verified"

# ── Step 3: build Bowtie2 index ───────────────────────────────
echo ""
echo "[3/4] Building Bowtie2 index"
echo "Index prefix: $INDEX_PREFIX"
echo ""

bowtie2-build \
  --threads "${SLURM_CPUS_PER_TASK}" \
  "$MASTER_FA" \
  "$INDEX_PREFIX"

# Note: add --large-index for references >4 GB:
# bowtie2-build --large-index --threads 16 reference.fa index_prefix

# ── Step 4: verify index files ────────────────────────────────
echo ""
echo "[4/4] Verifying index files"

EXPECTED_FILES=(
  "${INDEX_PREFIX}.1.bt2"
  "${INDEX_PREFIX}.2.bt2"
  "${INDEX_PREFIX}.3.bt2"
  "${INDEX_PREFIX}.4.bt2"
  "${INDEX_PREFIX}.rev.1.bt2"
  "${INDEX_PREFIX}.rev.2.bt2"
)

ALL_OK=true
for F in "${EXPECTED_FILES[@]}"; do
  if [ -f "$F" ]; then
    SIZE=$(du -sh "$F" | cut -f1)
    echo "  ✓ $(basename $F)  ($SIZE)"
  else
    echo "  ✗ MISSING: $(basename $F)"
    ALL_OK=false
  fi
done

if [ "$ALL_OK" = true ]; then
  echo ""
  echo "=============================================="
  echo " Index build complete ✓"
  echo " Total contigs indexed: $MASTER_N"
  echo " Index location: $INDEX_DIR/"
  echo " Index prefix:   $INDEX_PREFIX"
  echo ""
  echo " → Ready for Day 4 mapping"
  echo "   Use prefix: $INDEX_PREFIX"
  echo " End: $(date)"
  echo "=============================================="
else
  echo ""
  echo "ERROR: Some index files are missing. Re-run this script."
  exit 1
fi
