#!/bin/bash
#SBATCH --job-name=inspect_reference
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=0:30:00
#SBATCH --output=logs/inspect_reference_%j.out
#SBATCH --error=logs/inspect_reference_%j.err

# =============================================================
# Day 3 — Reference inspection & sanity checks
# Metatranscriptomics Series · Campbell Lab · 2026
#
# Run this BEFORE building your index to understand what is
# in your reference. Produces a summary report covering:
#   - Number of contigs / sequences
#   - Total reference length (bp)
#   - Contig length distribution (min, max, N50, mean)
#   - GC content (proxy for taxonomic breadth)
#   - Duplicate sequence ID check (will break Bowtie2)
#
# Works on any FASTA: MAGs, genome database, gene catalog.
#
# Usage: sbatch 00_inspect_reference.slurm
#   or:  bash 00_inspect_reference.slurm (interactive)
# =============================================================

# ── set your reference FASTA here ────────────────────────────
REF_FA="/project/bcampb7/camplab/MICRO_8130_2026/metatranscriptome/MAGS/master_mags.fa"

# Label for the report
REF_LABEL="master_mags (4 MAGs)"

REPORT="logs/reference_inspection_$(date +%Y%m%d_%H%M%S).txt"
mkdir -p logs

echo "=============================================="  | tee "$REPORT"
echo " Reference inspection report"                   | tee -a "$REPORT"
echo " File:  $REF_FA"                                | tee -a "$REPORT"
echo " Label: $REF_LABEL"                             | tee -a "$REPORT"
echo " Run:   $(date)"                                | tee -a "$REPORT"
echo "=============================================="  | tee -a "$REPORT"

# ── basic checks ──────────────────────────────────────────────
if [ ! -f "$REF_FA" ]; then
  echo "ERROR: File not found: $REF_FA" | tee -a "$REPORT"
  exit 1
fi

FILE_SIZE=$(du -sh "$REF_FA" | cut -f1)
echo "File size:     $FILE_SIZE"                      | tee -a "$REPORT"

# ── contig count ──────────────────────────────────────────────
N_CONTIGS=$(grep -c ">" "$REF_FA")
echo "Contig count:  $N_CONTIGS"                      | tee -a "$REPORT"

# ── total length + GC + length stats via awk ─────────────────
echo ""                                               | tee -a "$REPORT"
echo "Computing sequence statistics..."               | tee -a "$REPORT"

awk '
/^>/ {
  if (seq != "") {
    len = length(seq)
    total_len += len
    n_seqs++
    lens[n_seqs] = len
    # GC count
    gc = gsub(/[GCgc]/, "", seq)
    total_gc += gc
  }
  seq = ""
  next
}
{ seq = seq $0 }
END {
  if (seq != "") {
    len = length(seq)
    total_len += len
    n_seqs++
    lens[n_seqs] = len
    gc = gsub(/[GCgc]/, "", seq)
    total_gc += gc
  }

  # Sort lengths
  n = asort(lens)

  min_len = lens[1]
  max_len = lens[n]
  mean_len = total_len / n_seqs

  # N50
  half = total_len / 2
  cumsum = 0
  n50 = 0
  for (i = n; i >= 1; i--) {
    cumsum += lens[i]
    if (cumsum >= half) { n50 = lens[i]; break }
  }

  gc_pct = (total_gc / total_len) * 100

  printf "\n"
  printf "Total length:  %d bp (%.1f Mbp)\n", total_len, total_len/1e6
  printf "Min length:    %d bp\n", min_len
  printf "Max length:    %d bp\n", max_len
  printf "Mean length:   %.0f bp\n", mean_len
  printf "N50:           %d bp\n", n50
  printf "GC content:    %.1f%%\n", gc_pct
}
' "$REF_FA" | tee -a "$REPORT"

# ── duplicate ID check ────────────────────────────────────────
echo ""                                               | tee -a "$REPORT"
echo "Checking for duplicate sequence IDs..."         | tee -a "$REPORT"

DUP_COUNT=$(grep "^>" "$REF_FA" | sort | uniq -d | wc -l)

if [ "$DUP_COUNT" -eq 0 ]; then
  echo "  ✓ No duplicate IDs found — safe to index"   | tee -a "$REPORT"
else
  echo "  ✗ FOUND $DUP_COUNT duplicate sequence IDs!" | tee -a "$REPORT"
  echo "    Duplicate IDs will cause Bowtie2 to fail." | tee -a "$REPORT"
  echo "    Fix with: awk script to rename duplicates" | tee -a "$REPORT"
  echo "    or: seqkit rename -n master_mags.fa > master_mags_renamed.fa" | tee -a "$REPORT"
fi

# ── index status check ────────────────────────────────────────
echo ""                                               | tee -a "$REPORT"
echo "Checking for existing Bowtie2 index..."         | tee -a "$REPORT"

INDEX_DIR="$(dirname $REF_FA)/bowtie_index"
INDEX_BASE="$(basename $REF_FA .fa)"

BT2_COUNT=$(ls ${INDEX_DIR}/${INDEX_BASE}*.bt2 2>/dev/null | wc -l)
BT2L_COUNT=$(ls ${INDEX_DIR}/${INDEX_BASE}*.bt2l 2>/dev/null | wc -l)

if [ "$BT2_COUNT" -eq 6 ]; then
  echo "  ✓ Standard index (.bt2) found — $BT2_COUNT files" | tee -a "$REPORT"
  echo "    Location: ${INDEX_DIR}/"                         | tee -a "$REPORT"
elif [ "$BT2L_COUNT" -eq 6 ]; then
  echo "  ✓ Large index (.bt2l) found — $BT2L_COUNT files"  | tee -a "$REPORT"
  echo "    Location: ${INDEX_DIR}/"                         | tee -a "$REPORT"
else
  echo "  → No index found. Run 01_build_mag_index.slurm"   | tee -a "$REPORT"
fi

# ── summary ───────────────────────────────────────────────────
echo ""                                               | tee -a "$REPORT"
echo "=============================================="  | tee -a "$REPORT"
echo " Inspection complete"                           | tee -a "$REPORT"
echo " Full report saved to: $REPORT"                 | tee -a "$REPORT"
echo ""                                               | tee -a "$REPORT"

if [ "$N_CONTIGS" -lt 100 ]; then
  echo " ⚠ Very few contigs ($N_CONTIGS). Expect very low"  | tee -a "$REPORT"
  echo "   alignment rates (<1%). Consider adding more"      | tee -a "$REPORT"
  echo "   MAGs or using a reference database."              | tee -a "$REPORT"
elif [ "$N_CONTIGS" -lt 5000 ]; then
  echo " ℹ Modest reference size ($N_CONTIGS contigs)."     | tee -a "$REPORT"
  echo "   Alignment rate will reflect reference coverage."  | tee -a "$REPORT"
  echo "   Low alignment is expected and informative."       | tee -a "$REPORT"
else
  echo " ✓ Reference is large enough for broad mapping."    | tee -a "$REPORT"
fi

echo "=============================================="  | tee -a "$REPORT"
echo " End: $(date)"                                  | tee -a "$REPORT"
echo "=============================================="  | tee -a "$REPORT"
