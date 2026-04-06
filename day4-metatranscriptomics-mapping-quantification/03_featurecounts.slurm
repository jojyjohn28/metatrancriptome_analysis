#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --output=logs/featurecounts_%j.out
#SBATCH --error=logs/featurecounts_%j.err

# =============================================================
# Day 4 — featureCounts gene expression quantification
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Counts reads overlapping annotated CDS features in the GFF.
# Runs on ALL sorted BAM files at once — produces a single
# count matrix (genes × samples) ready for DESeq2 / TPM.
#
# Input:  bam/*_sorted.bam  (all samples)
#         MAGS/master_MAGs.fixed.gff
# Output: counts/all_samples_counts.txt
#         counts/all_samples_counts.txt.summary
#
# GFF requirement:
#   The GFF must have:
#     - Column 3 (feature type): CDS  (we use -t CDS)
#     - Column 9 (attributes):   ID=  (we use -g ID)
#   If your GFF uses 'gene' instead of 'CDS', change -t to 'gene'
#   If your GFF uses 'Name' instead of 'ID', change -g to 'Name'
#   Check your GFF: head -20 master_MAGs.fixed.gff
#
# Usage: sbatch 03_featurecounts.slurm
# =============================================================

module load subread

# ── paths — edit these ────────────────────────────────────────
BAM_DIR="bam"
GFF="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/MAGS/master_MAGs.fixed.gff"
OUT_DIR="counts"

mkdir -p "$OUT_DIR" logs

echo "=============================================="
echo " featureCounts quantification"
echo " Start: $(date)"
echo "=============================================="

# ── collect all sorted BAM files ─────────────────────────────
BAM_FILES=($(ls ${BAM_DIR}/*_sorted.bam))

if [ ${#BAM_FILES[@]} -eq 0 ]; then
  echo "ERROR: No sorted BAM files found in $BAM_DIR"
  echo "Run 02_sam_to_bam.slurm for all samples first."
  exit 1
fi

echo ""
echo "BAM files to quantify:"
for B in "${BAM_FILES[@]}"; do
  echo "  $(basename $B)"
done
echo ""

# ── check GFF ────────────────────────────────────────────────
if [ ! -f "$GFF" ]; then
  echo "ERROR: GFF not found: $GFF"
  echo "Provide DRAM or Prokka annotation GFF for your MAGs."
  exit 1
fi

GFF_FEATURES=$(grep -v "^#" "$GFF" | awk '{print $3}' | sort | uniq -c | sort -rn)
echo "GFF feature types (for reference):"
echo "$GFF_FEATURES"
echo ""

# ── run featureCounts ─────────────────────────────────────────
# -T 8    : threads
# -p      : paired-end mode (counts fragments, not reads)
# -B      : both mates must map (remove if you want singletons counted)
# -C      : exclude chimeric pairs (mates on different contigs)
# -t CDS  : count CDS features only
# -g ID   : use ID attribute as gene name
# -a      : annotation GFF
# -o      : output count table
# remaining args: all BAM files

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -p \
  -B \
  -C \
  -t CDS \
  -g ID \
  -a "$GFF" \
  -o "${OUT_DIR}/all_samples_counts.txt" \
  "${BAM_FILES[@]}"

# ── check output ──────────────────────────────────────────────
if [ ! -f "${OUT_DIR}/all_samples_counts.txt" ]; then
  echo "ERROR: featureCounts output not produced. Check logs."
  exit 1
fi

# ── print summary ─────────────────────────────────────────────
echo ""
echo "featureCounts summary:"
cat "${OUT_DIR}/all_samples_counts.txt.summary"

# Count genes with at least 1 read across any sample
TOTAL_GENES=$(tail -n +3 "${OUT_DIR}/all_samples_counts.txt" | wc -l)
EXPRESSED=$(tail -n +3 "${OUT_DIR}/all_samples_counts.txt" | \
  awk 'BEGIN{n=0}{for(i=7;i<=NF;i++){if($i>0){n++;break}}}END{print n}')
ZERO=$(( TOTAL_GENES - EXPRESSED ))

echo ""
echo "=============================================="
echo " Quantification complete"
echo "----------------------------------------------"
echo " Total annotated genes:  $TOTAL_GENES"
echo " Genes with ≥1 read:     $EXPRESSED"
echo " Genes with 0 reads:     $ZERO"
echo ""

if (( $(echo "$EXPRESSED * 100 / $TOTAL_GENES < 10" | bc -l) )); then
  echo " ℹ  Low fraction of expressed genes — expected"
  echo "    with a small MAG reference and low alignment rate."
  echo "    Non-zero genes ARE your signal. See Day 5."
fi

echo ""
echo " Count table: ${OUT_DIR}/all_samples_counts.txt"
echo " Summary:     ${OUT_DIR}/all_samples_counts.txt.summary"
echo " → Load into R for Day 5: DESeq2 and TPM normalization"
echo " End: $(date)"
echo "=============================================="
