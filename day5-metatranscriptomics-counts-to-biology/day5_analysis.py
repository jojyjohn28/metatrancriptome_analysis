#!/usr/bin/env python3
"""
Day 5 — From Counts to Biology (Python version)
Metatranscriptomics Series · Campbell Lab · 2026

Covers:
  1. TPM normalization from featureCounts output
  2. DNA:RNA ratio (if paired metagenome TPM available)
  3. Top expressed gene heatmap
  4. Per-sample expression summary statistics
  5. Functional gene family expression (with annotation table)

Input:
  counts/all_samples_counts.txt     - featureCounts output (Day 4)
  counts/DNA_TPM_matrix.csv         - from paired metagenome (optional)
  metadata/sample_metadata.csv      - SampleID, Site, Season, Salinity
  MAGS/master_MAGs_annotations.tsv  - DRAM annotation (optional)

Usage:
  python day5_analysis.py
  (or: Rscript day5_analysis.R for the full DESeq2 workflow)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for HPC
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# ── Output directories ────────────────────────────────────────
OUT_DIR  = Path("results")
PLOT_DIR = Path("plots")
OUT_DIR.mkdir(exist_ok=True)
PLOT_DIR.mkdir(exist_ok=True)

# ── File paths ────────────────────────────────────────────────
COUNTS_FILE   = "counts/all_samples_counts.txt"
DNA_TPM_FILE  = "counts/DNA_TPM_matrix.csv"
METADATA_FILE = "metadata/sample_metadata.csv"
ANNOT_FILE    = "MAGS/master_MAGs_annotations.tsv"

print("=" * 55)
print(" Day 5 — Counts to Biology (Python)")
print("=" * 55)

# ──────────────────────────────────────────────────────────────
# 1. Load featureCounts output
# ──────────────────────────────────────────────────────────────
print("\n[1/5] Loading count matrix...")

try:
    counts_raw = pd.read_csv(COUNTS_FILE, sep="\t", skiprows=1, index_col=0)
except FileNotFoundError:
    sys.exit(f"ERROR: Count file not found: {COUNTS_FILE}\nRun Day 4 featureCounts first.")

gene_lengths = counts_raw["Length"]

# featureCounts: Chr, Start, End, Strand, Length are cols 1-5 (after Geneid)
# Sample count columns start at index 5 (0-based)
sample_cols  = [c for c in counts_raw.columns if c not in
                ["Chr","Start","End","Strand","Length"]]
count_matrix = counts_raw[sample_cols].copy()

# Clean up sample names (featureCounts appends full path)
count_matrix.columns = [Path(c).stem.replace("_sorted","") for c in count_matrix.columns]

print(f"  Genes:   {count_matrix.shape[0]}")
print(f"  Samples: {count_matrix.shape[1]}")
print(f"  Samples: {list(count_matrix.columns)}")

# ──────────────────────────────────────────────────────────────
# 2. TPM normalization
# ──────────────────────────────────────────────────────────────
print("\n[2/5] Computing TPM...")

def compute_tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """
    counts  : genes × samples DataFrame (raw integers)
    lengths : gene lengths in bp (Series, same index as counts)
    returns : TPM DataFrame (same shape as counts)
    """
    # Step 1: RPK — correct for gene length
    rpk = counts.div(lengths / 1000, axis=0)
    # Step 2: scale to per-million
    tpm = rpk.div(rpk.sum(axis=0) / 1e6, axis=1)
    return tpm

TPM = compute_tpm(count_matrix, gene_lengths)

print("  Column sums (should all = 1,000,000):")
for col, val in TPM.sum().items():
    print(f"    {col}: {val:,.0f}")

TPM.to_csv(OUT_DIR / "TPM_matrix.csv")
print(f"  Saved: {OUT_DIR / 'TPM_matrix.csv'}")

# ──────────────────────────────────────────────────────────────
# 3. DNA:RNA ratio
# ──────────────────────────────────────────────────────────────
print("\n[3/5] DNA:RNA ratio...")

if Path(DNA_TPM_FILE).exists():
    DNA_TPM = pd.read_csv(DNA_TPM_FILE, index_col=0)

    shared_genes   = TPM.index.intersection(DNA_TPM.index)
    shared_samples = TPM.columns.intersection(DNA_TPM.columns)

    print(f"  Shared genes:   {len(shared_genes)}")
    print(f"  Shared samples: {len(shared_samples)}")

    RNA_sub = TPM.loc[shared_genes, shared_samples]
    DNA_sub = DNA_TPM.loc[shared_genes, shared_samples]

    # log2 ratio with pseudocount to avoid log(0)
    ratio = np.log2((RNA_sub + 1) / (DNA_sub + 1))
    ratio.to_csv(OUT_DIR / "DNA_RNA_ratio.csv")
    print(f"  Saved: {OUT_DIR / 'DNA_RNA_ratio.csv'}")

    # Interpret ratio values
    for col in ratio.columns:
        pos = (ratio[col] > 1).sum()
        neg = (ratio[col] < -1).sum()
        print(f"    {col}: {pos} genes actively upregulated (ratio>1), "
              f"{neg} genes silenced (ratio<-1)")

    # Boxplot of ratio distribution per sample
    fig, ax = plt.subplots(figsize=(max(6, len(shared_samples) * 0.9), 4))
    ratio.boxplot(ax=ax, grid=False)
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, label="ratio = 0")
    ax.set_title("DNA:RNA ratio distribution per sample")
    ax.set_ylabel("log2(RNA_TPM+1 / DNA_TPM+1)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "DNA_RNA_ratio_boxplot.pdf")
    plt.close()
    print(f"  Plot saved: {PLOT_DIR / 'DNA_RNA_ratio_boxplot.pdf'}")

else:
    print(f"  DNA TPM file not found: {DNA_TPM_FILE}")
    print("  To enable DNA:RNA ratio:")
    print("  1. Map paired metagenome reads to master_mags.fa")
    print("  2. Run featureCounts on metagenome BAMs")
    print("  3. Compute TPM the same way as metatranscriptome")
    print("  4. Save as counts/DNA_TPM_matrix.csv")

# ──────────────────────────────────────────────────────────────
# 4. Expression summary + top-50 heatmap
# ──────────────────────────────────────────────────────────────
print("\n[4/5] Expression summary and heatmap...")

total_genes = count_matrix.shape[0]
expressed   = (count_matrix.sum(axis=1) > 0).sum()
mean_tpm    = TPM.mean(axis=1)
top10       = mean_tpm.nlargest(10)

print(f"  Total genes:              {total_genes}")
print(f"  Genes with ≥1 count:      {expressed}")
print(f"  Zero-count genes:         {total_genes - expressed}")
print(f"  Median mean TPM:          {mean_tpm.median():.2f}")
print(f"\n  Top 10 expressed genes (mean TPM):")
for gene, val in top10.items():
    print(f"    {gene[:50]:<50} {val:.1f}")

# Save expression summary
summary_df = pd.DataFrame({
    "mean_TPM":   TPM.mean(axis=1),
    "max_TPM":    TPM.max(axis=1),
    "n_samples_expressed": (TPM > 0).sum(axis=1)
}).sort_values("mean_TPM", ascending=False)
summary_df.to_csv(OUT_DIR / "gene_expression_summary.csv")
print(f"\n  Summary saved: {OUT_DIR / 'gene_expression_summary.csv'}")

# Top-50 heatmap
top50_genes = mean_tpm.nlargest(50).index
top50_tpm   = np.log2(TPM.loc[top50_genes] + 1)

fig, ax = plt.subplots(figsize=(max(8, len(TPM.columns) * 1.0), 12))
sns.heatmap(
    top50_tpm,
    ax=ax,
    cmap="YlOrRd",
    xticklabels=True,
    yticklabels=False,
    linewidths=0,
    cbar_kws={"label": "log2(TPM + 1)", "shrink": 0.6}
)
ax.set_title("Top 50 expressed genes — log2(TPM+1)", fontsize=13, pad=12)
ax.set_xlabel("Sample")
plt.tight_layout()
plt.savefig(PLOT_DIR / "top50_TPM_heatmap.pdf", bbox_inches="tight")
plt.close()
print(f"  Heatmap saved: {PLOT_DIR / 'top50_TPM_heatmap.pdf'}")

# ──────────────────────────────────────────────────────────────
# 5. Functional gene family expression
# ──────────────────────────────────────────────────────────────
print("\n[5/5] Functional gene family expression...")

if Path(ANNOT_FILE).exists():
    annots = pd.read_csv(ANNOT_FILE, sep="\t", low_memory=False)
    annots.columns = annots.columns.str.strip()

    tpm_annot = TPM.copy()
    tpm_annot.index.name = "Gene"
    tpm_annot = tpm_annot.reset_index()

    if "Gene" in annots.columns and "gene_description" in annots.columns:
        tpm_annot = tpm_annot.merge(
            annots[["Gene","gene_description"]].drop_duplicates("Gene"),
            on="Gene", how="left"
        )
    else:
        print("  Warning: annotation file missing 'Gene' or 'gene_description' columns")
        print(f"  Columns found: {list(annots.columns[:8])}")

    sample_cols_tpm = [c for c in tpm_annot.columns
                       if c not in ["Gene","gene_description"]]

    # Target gene families — edit this dict to add your genes of interest
    gene_families = {
        "dsrA (sulfate reduction)":     ["dsrA","dissimilatory sulfite reductase"],
        "alkB (alkane degradation)":    ["alkB","alkane monooxygenase","alkane 1-monooxygenase"],
        "mcrA (methanogenesis)":        ["mcrA","methyl-coenzyme M reductase"],
        "amoA (nitrification)":         ["amoA","ammonia monooxygenase"],
        "hzo (anammox)":                ["hzo","hydrazine oxidoreductase"],
    }

    results_families = {}
    for family_label, keywords in gene_families.items():
        pattern = "|".join(keywords)
        if "gene_description" in tpm_annot.columns:
            subset = tpm_annot[tpm_annot["gene_description"].str.contains(
                pattern, case=False, na=False)]
        else:
            subset = tpm_annot[tpm_annot["Gene"].str.contains(
                pattern, case=False, na=False)]

        if len(subset) > 0:
            family_sum = subset[sample_cols_tpm].sum(axis=0)
            results_families[family_label] = family_sum
            print(f"  {family_label}: {len(subset)} genes found")
        else:
            print(f"  {family_label}: not found in annotations")

    if results_families:
        func_df = pd.DataFrame(results_families).T
        func_df.to_csv(OUT_DIR / "functional_gene_expression.csv")
        print(f"  Functional expression saved: {OUT_DIR / 'functional_gene_expression.csv'}")

        # Heatmap of functional families across samples
        if func_df.shape[0] > 1:
            fig, ax = plt.subplots(figsize=(max(7, len(sample_cols_tpm) * 0.9), 4))
            sns.heatmap(
                np.log2(func_df + 1),
                ax=ax, cmap="Greens", annot=True, fmt=".1f",
                linewidths=0.5, cbar_kws={"label": "log2(sum TPM + 1)"}
            )
            ax.set_title("Community-level functional gene expression")
            ax.set_xlabel("Sample")
            plt.tight_layout()
            plt.savefig(PLOT_DIR / "functional_gene_heatmap.pdf", bbox_inches="tight")
            plt.close()
            print(f"  Functional heatmap saved: {PLOT_DIR / 'functional_gene_heatmap.pdf'}")

else:
    print(f"  Annotation file not found: {ANNOT_FILE}")
    print("  Provide DRAM annotation TSV to enable functional analysis.")

# ── Final summary ─────────────────────────────────────────────
print("\n" + "=" * 55)
print(" Day 5 analysis complete")
print("-" * 55)
print(f" TPM matrix:         {OUT_DIR / 'TPM_matrix.csv'}")
print(f" Expression summary: {OUT_DIR / 'gene_expression_summary.csv'}")
print(f" Plots:              {PLOT_DIR}/")
print(f"\n For DESeq2 differential expression:")
print(f"   Rscript day5_analysis.R")
print("=" * 55)
