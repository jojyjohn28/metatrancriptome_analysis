# Day 5 — From Counts to Biology 🌊

### Normalization · DESeq2 · DNA:RNA Ratios · Ecological Interpretation

This is the final day of the series. The count matrix from Day 4 is your raw material — this day transforms it into biological insight using TPM normalization, differential expression analysis, and ecological interpretation frameworks.

📘 **Blog post:** https://jojyjohn28.github.io/blog/day5-metatranscriptomics-counts-to-biology/

---

## 🎯 What You Will Do

- Normalize raw counts to TPM for within-sample gene expression comparisons
- Run DESeq2 differential expression across three environmental factors (Site, Season, Salinity)
- Generate PCA, sample distance heatmaps, and volcano plots
- Calculate DNA:RNA ratios using paired metagenome data (optional but recommended)
- Quantify community-level functional gene expression (dsrA, alkB, and others)

---

## 📁 Folder Structure

```
Day5_Counts_to_Biology/
├── scripts/
│   ├── day5_analysis.R          # DESeq2 · PCA · volcano plots · DNA:RNA ratio
│   ├── day5_analysis.py         # TPM · heatmaps · functional gene expression
│   └── install_packages.R       # Install all required R packages (run once)
└── README.md
```

---

## 🧰 Tools Used

| Tool                     | Purpose                                          |
| ------------------------ | ------------------------------------------------ |
| DESeq2 (R)               | Differential expression · PCA · variance heatmap |
| ggplot2 + ggrepel (R)    | Volcano plots with gene labels                   |
| pheatmap (R)             | Sample distance and DNA:RNA ratio heatmaps       |
| Python (pandas, seaborn) | TPM normalization · expression heatmaps          |

---

## 🚀 Quick Start

**Step 1 — Install R packages (first time only):**

```bash
Rscript scripts/install_packages.R
```

**Step 2 — Run with toy data (test that everything works):**

```r
# In R — uses built-in toy data, no real data needed
source("scripts/Deseq_with_toydata.R")
```

**Step 3 — Run with your real data:**

```bash
# On HPC (runs both R and Python):
sbatch scripts/day5_analysis.slurm

# Locally:
Rscript scripts/day5_analysis.R
python  scripts/day5_analysis.py
```

---

## 📥 Input Files Required

| File                               | Description             | Where it comes from |
| ---------------------------------- | ----------------------- | ------------------- |
| `counts/all_samples_counts.txt`    | Raw count matrix        | Day 4 featureCounts |
| `metadata/sample_metadata.csv`     | Sample information      | You create this     |
| `counts/DNA_TPM_matrix.csv`        | DNA TPM from metagenome | Day 4 (optional)    |
| `MAGS/master_MAGs_annotations.tsv` | Gene annotations        | DRAM (optional)     |

**Metadata format** (`sample_metadata.csv`):

```
SampleID,       Site, Season, Salinity
CP_Spr15G08,    PA,   Spring, Low
CP_Spr15L08,    PA,   Spring, High
CP_Spr31G08,    FL,   Spring, Low
```

---

## 📤 Output Files

| File                                     | Description                              |
| ---------------------------------------- | ---------------------------------------- |
| `results/TPM_matrix.csv`                 | TPM-normalized expression matrix         |
| `results/DESeq2_all_factors.csv`         | DESeq2 results for all three factors     |
| `results/DNA_RNA_ratio.csv`              | log2 DNA:RNA ratio per gene per sample   |
| `results/functional_gene_expression.csv` | Community-level pathway expression       |
| `plots/PCA_samples.pdf`                  | PCA of samples (VST-transformed)         |
| `plots/sample_distance_heatmap.pdf`      | Sample-to-sample distance matrix         |
| `plots/volcano_all_factors.pdf`          | Volcano plots (Site · Season · Salinity) |
| `plots/volcano_labeled.pdf`              | Volcano plots with top 5 gene labels     |
| `plots/DNA_RNA_ratio_heatmap.pdf`        | Heatmap of top variable ratio genes      |
| `plots/top50_TPM_heatmap.pdf`            | Top 50 expressed genes (log2 TPM)        |

---

## 📢 Notes

- DESeq2 requires **raw integer counts** — do not pass TPM or RPKM to DESeq2
- The DNA:RNA ratio requires a paired metagenome mapped to the same reference (Day 3 index); skip if unavailable
- With a small MAG reference (~2% alignment rate), expect few significant genes — non-zero counts are real signal
- For co-expression network analysis of these results, see the [WGCNA post](https://jojyjohn28.github.io/series/)
- For multi-omics integration, see the [multi-omics post](https://jojyjohn28.github.io/series/)

---

Last updated: April 3, 2026\_
