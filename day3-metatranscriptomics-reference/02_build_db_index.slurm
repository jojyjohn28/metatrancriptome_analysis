#!/bin/bash
#SBATCH --job-name=build_db_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=220GB
#SBATCH --time=8:00:00
#SBATCH --output=logs/build_db_index_%j.out
#SBATCH --error=logs/build_db_index_%j.err

# =============================================================
# Day 3 — Build Bowtie2 index from a reference genome database
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Use this script when mapping against a reference database
# instead of (or in addition to) your own MAGs. Examples:
#   - GTDB-Tk representative genomes
#   - NCBI RefSeq bacterial genomes
#   - Environment-specific MAG catalogs (TARA Oceans, etc.)
#
# KEY DIFFERENCES from MAG index script:
#   - Requires --large-index flag (reference >4 GB)
#   - Much higher memory requirement (200+ GB)
#   - Much longer runtime (2–8 hours)
#   - Index files will be .bt2l (large format), not .bt2
#
# DOWNLOAD OPTIONS (uncomment the one you need):
#
#   GTDB-Tk representative genomes:
#     https://data.gtdb.ecogenomics.org/releases/latest/
#     → genomic_files_reps/gtdb_genomes_reps.tar.gz
#
#   NCBI RefSeq bacteria:
#     datasets download genome taxon bacteria \
#       --reference-only --include genome
#
#   TARA Oceans MAG catalog (marine):
#     https://www.ebi.ac.uk/biostudies/studies/S-BSST297
#
# Usage: sbatch 02_build_db_index.slurm
# =============================================================

module load bowtie2/2.5.4

# ── paths — edit these ────────────────────────────────────────
# Set DB_FA to your concatenated reference database FASTA
DB_FA="/project/bcampb7/camplab/databases/gtdb_representatives.fa"
INDEX_DIR="/project/bcampb7/camplab/databases/gtdb_bowtie_index"
INDEX_PREFIX="${INDEX_DIR}/gtdb_index"

mkdir -p "$INDEX_DIR" logs

echo "=============================================="
echo " Day 3 — Build reference database Bowtie2 index"
echo " Reference: $DB_FA"
echo " Start: $(date)"
echo "=============================================="

# ── check reference file ──────────────────────────────────────
if [ ! -f "$DB_FA" ]; then
  echo "ERROR: Reference file not found: $DB_FA"
  echo "       Download and concatenate your reference genomes first."
  echo ""
  echo "  Example — GTDB representative genomes:"
  echo "    cd /path/to/databases"
  echo "    wget https://data.gtdb.ecogenomics.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz"
  echo "    tar -xzf gtdb_genomes_reps.tar.gz"
  echo "    cat gtdb_genomes_reps/**/*.fna > gtdb_representatives.fa"
  exit 1
fi

DB_SIZE=$(du -sh "$DB_FA" | cut -f1)
DB_CONTIGS=$(grep -c ">" "$DB_FA")
echo "Reference file:   $DB_FA"
echo "File size:        $DB_SIZE"
echo "Contig count:     $DB_CONTIGS"
echo ""

# ── build large-format Bowtie2 index ─────────────────────────
echo "Building Bowtie2 index (large format — this takes time)..."
echo ""

# --large-index is REQUIRED for references larger than 4 GB
# It produces .bt2l files instead of .bt2
bowtie2-build \
  --large-index \
  --threads "${SLURM_CPUS_PER_TASK}" \
  "$DB_FA" \
  "$INDEX_PREFIX"

# ── verify index files ────────────────────────────────────────
echo ""
echo "Verifying index files:"

EXPECTED=(
  "${INDEX_PREFIX}.1.bt2l"
  "${INDEX_PREFIX}.2.bt2l"
  "${INDEX_PREFIX}.3.bt2l"
  "${INDEX_PREFIX}.4.bt2l"
  "${INDEX_PREFIX}.rev.1.bt2l"
  "${INDEX_PREFIX}.rev.2.bt2l"
)

ALL_OK=true
for F in "${EXPECTED[@]}"; do
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
  echo " Database index build complete ✓"
  echo " Contigs indexed: $DB_CONTIGS"
  echo " Index location:  $INDEX_DIR/"
  echo " Index prefix:    $INDEX_PREFIX"
  echo ""
  echo " NOTE: Use this prefix in Day 4 bowtie2 mapping"
  echo "       -x $INDEX_PREFIX"
  echo ""
  echo " Expected alignment rate vs MAG index:"
  echo "   MAG index (4 MAGs):    ~2%    (few organisms)"
  echo "   GTDB index (all reps): ~15-40% (broader coverage)"
  echo "   Actual rate depends on your community composition."
  echo " End: $(date)"
  echo "=============================================="
else
  echo "ERROR: Index build incomplete. Check logs and retry."
  exit 1
fi
