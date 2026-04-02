#!/bin/bash
#SBATCH --job-name=build_catalog_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=6:00:00
#SBATCH --output=logs/build_catalog_index_%j.out
#SBATCH --error=logs/build_catalog_index_%j.err

# =============================================================
# Day 3 — Build Bowtie2 index from a gene catalog
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Gene catalog mapping: an alternative to genome-level mapping.
# Instead of full genomes, you map against a non-redundant
# catalog of predicted protein-coding genes.
#
# Advantages vs genome mapping:
#   + Much broader functional coverage
#   + Captures genes from unbinned/unassembled organisms
#   + Directly links reads to functional annotations (KEGG, COG)
#
# Disadvantages:
#   - Taxonomic resolution requires extra steps
#   - You know the function; you need provenance to know the organism
#
# Commonly used gene catalogs:
#   OM-RGC v2 (marine):    https://www.ebi.ac.uk/biostudies/studies/S-BSST297
#   IGC v2 (human gut):    https://db.cngb.org/microbiome/igc/
#   Soil gene catalog:     (various publications)
#   Custom (from your own metagenomes): see Prodigal + CD-HIT workflow below
#
# Usage: sbatch 03_build_catalog_index.slurm
# =============================================================

module load bowtie2/2.5.4

# ── paths — edit for your catalog ────────────────────────────
CATALOG_FA="/project/bcampb7/camplab/databases/OM-RGC_v2/OM-RGC_v2_genes.fa"
INDEX_DIR="/project/bcampb7/camplab/databases/OM-RGC_v2/bowtie_index"
INDEX_PREFIX="${INDEX_DIR}/omrgc_v2_index"

mkdir -p "$INDEX_DIR" logs

echo "=============================================="
echo " Day 3 — Build gene catalog Bowtie2 index"
echo " Catalog: $CATALOG_FA"
echo " Start: $(date)"
echo "=============================================="

# ── check catalog file ────────────────────────────────────────
if [ ! -f "$CATALOG_FA" ]; then
  echo "ERROR: Gene catalog not found: $CATALOG_FA"
  echo ""
  echo "  To build a custom gene catalog from your metagenomes:"
  echo "  ── Step 1: predict genes with Prodigal ─────────────"
  echo "    prodigal -i assembly.fa -a proteins.faa -d genes.fna -f gff -p meta"
  echo ""
  echo "  ── Step 2: cluster at 95% identity with CD-HIT ─────"
  echo "    cd-hit-est -i all_samples_genes.fna \\"
  echo "      -o catalog_95.fa \\"
  echo "      -c 0.95 -aS 0.9 -n 8 -d 0 -T 16 -M 60000"
  echo ""
  echo "  ── Step 3: annotate with KEGG / COG ────────────────"
  echo "    # Use DRAM, eggNOG-mapper, or HMMER"
  echo ""
  exit 1
fi

CATALOG_SIZE=$(du -sh "$CATALOG_FA" | cut -f1)
GENE_COUNT=$(grep -c ">" "$CATALOG_FA")

echo "Gene catalog: $CATALOG_FA"
echo "File size:    $CATALOG_SIZE"
echo "Gene count:   $GENE_COUNT"
echo ""

# ── choose index format based on file size ────────────────────
# Files < 4 GB: standard .bt2 format
# Files >= 4 GB: --large-index (.bt2l format)
SIZE_BYTES=$(stat -c%s "$CATALOG_FA")
FOUR_GB=4294967296

if [ "$SIZE_BYTES" -gt "$FOUR_GB" ]; then
  echo "Reference >4 GB — using --large-index"
  LARGE_FLAG="--large-index"
  EXT="bt2l"
else
  echo "Reference <4 GB — using standard index"
  LARGE_FLAG=""
  EXT="bt2"
fi

# ── build index ───────────────────────────────────────────────
echo ""
echo "Building Bowtie2 index..."

bowtie2-build \
  $LARGE_FLAG \
  --threads "${SLURM_CPUS_PER_TASK}" \
  "$CATALOG_FA" \
  "$INDEX_PREFIX"

# ── verify ────────────────────────────────────────────────────
echo ""
echo "Verifying index files:"

ALL_OK=true
for SUFFIX in "1.${EXT}" "2.${EXT}" "3.${EXT}" "4.${EXT}" "rev.1.${EXT}" "rev.2.${EXT}"; do
  F="${INDEX_PREFIX}.${SUFFIX}"
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
  echo " Gene catalog index build complete ✓"
  echo " Genes indexed: $GENE_COUNT"
  echo " Index prefix:  $INDEX_PREFIX"
  echo ""
  echo " IMPORTANT notes for gene catalog mapping (Day 4):"
  echo "   Add --no-unal to bowtie2 — suppresses unmapped"
  echo "   reads in SAM output (critical for large catalogs)"
  echo ""
  echo "   featureCounts annotation:"
  echo "   Use a BED or SAF file derived from catalog coords"
  echo "   rather than a GFF — gene catalogs lack GFF structure"
  echo " End: $(date)"
  echo "=============================================="
else
  echo "ERROR: Index incomplete. Check logs and retry."
  exit 1
fi
