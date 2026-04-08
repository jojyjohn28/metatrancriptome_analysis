# Metatranscriptomics — From Raw RNA → Biology 🧬

### A 5-Day Metatranscriptomics Learning Series

This repository accompanies a **5-day blog series on metatranscriptomics**, guiding you from raw RNA sequencing reads to biological interpretation of microbial community activity. Each day focuses on a key step in the workflow, with detailed blog explanations, reproducible scripts, and toy datasets for hands-on practice.

📘 **Full blog series:** https://jojyjohn28.github.io/series/metatranscriptome/

---

## 🚀 What You Will Learn

By the end of this series, you will be able to:

- Understand what metatranscriptomics measures and how it differs from metagenomics
- Perform RNA-seq quality control, adapter trimming, and rRNA removal
- Diagnose and repair mismatched paired-end read files before mapping
- Choose the right reference strategy — MAGs, reference genomes, or gene catalogs
- Build a Bowtie2 index and map clean mRNA reads to a reference
- Quantify gene expression using featureCounts and interpret low alignment rates
- Normalize counts to TPM, run DESeq2 differential expression, and compute DNA:RNA ratios
- Connect expression data to ecological questions about active microbial communities

---

## 📁 Repository Structure

```
metatranscriptome_analysis/
│
├── day1-toy-data/
│   ├── toy_data/                        # Toy FASTQ files for pipeline testing
│   ├── bowtie_index/                    # Pre-built Bowtie2 index for toy MAGs
│   └── README.md
│
├── day2-qc-preprocessing/
│   ├── 00_check_pair_counts.slurm       # Validate R1/R2 read count sync
│   ├── 01_fastqc.slurm                  # FastQC on raw reads
│   ├── 02_fastp_trim.slurm              # Adapter trimming and quality filtering
│   ├── 03_sortmerna_rrna_removal.slurm  # rRNA removal (SortMeRNA)
│   ├── 03b_bbduk_rrna_removal.slurm     # rRNA removal alternative (bbduk)
│   ├── 04_repair_bbtools.slurm          # Fix mismatched paired-end files (repair.sh)
│   ├── 04b_fastp_paired_sync.slurm      # Fix/prevent pairing issues (fastp)
│   └── README.md
│
├── day3-metatranscriptomics-reference/
│   ├── 00_inspect_reference.slurm       # Reference QC: contig stats + duplicate ID check
│   ├── 01_build_mag_index.slurm         # Concatenate MAGs + build Bowtie2 index
│   ├── 02_build_db_index.slurm          # Build index from reference genome database
│   ├── 03_build_catalog_index.slurm     # Build index from gene catalog (OM-RGC etc.)
│   └── README.md
│
├── day4-metatranscriptomics-mapping-quantification/
│   ├── 01_bowtie2_mapping.slurm         # Map mRNA reads with Bowtie2 (array job)
│   ├── 02_sam_to_bam.slurm              # SAM → sorted BAM + flagstat (array job)
│   ├── 03_featurecounts.slurm           # Count reads per gene across all samples
│   └── README.md
│
├── day5-metatranscriptomics-counts-to-biology/
│   ├── Deseq_with_toydata.R             # DESeq2 workflow with toy count data
│   ├── TPM_normalization.R              # TPM calculation + DESeq2 + volcano plots
│   ├── day5_analysis.py                 # Python: TPM · DNA:RNA ratio · heatmaps
│   ├── Differential_abundance_cazy_all.pdf  # Example output figure
│   └── README.md
│
└── README.md                            # This file
```

---

## 📅 Series Overview

| Day       | Topic                    | Key Tools                                      | Scripts            |
| --------- | ------------------------ | ---------------------------------------------- | ------------------ |
| **Day 1** | Introduction + toy data  | —                                              | toy_data/          |
| **Day 2** | QC & preprocessing       | FastQC · fastp · SortMeRNA · bbduk · repair.sh | 8 SLURM scripts    |
| **Day 3** | Reference strategy       | Bowtie2-build                                  | 4 SLURM scripts    |
| **Day 4** | Mapping & quantification | Bowtie2 · SAMtools · featureCounts             | 3 SLURM scripts    |
| **Day 5** | Counts to biology        | DESeq2 · R · Python                            | 3 analysis scripts |

---

## 🧰 Tools Used

| Tool                             | Purpose                                          |
| -------------------------------- | ------------------------------------------------ |
| FastQC                           | Raw read quality assessment                      |
| fastp                            | Adapter trimming and quality filtering           |
| SortMeRNA                        | rRNA removal (multi-database)                    |
| bbduk / repair.sh (BBTools)      | rRNA removal (alternative) + paired-end repair   |
| Bowtie2                          | Read alignment to reference                      |
| SAMtools                         | SAM/BAM conversion, sorting, and diagnostics     |
| featureCounts (Subread)          | Gene-level read quantification                   |
| DESeq2 (R)                       | Differential expression analysis                 |
| ggplot2 / ggrepel / pheatmap (R) | Visualization                                    |
| Python (pandas, seaborn)         | TPM normalization, heatmaps, functional analysis |

---

## 🖥️ HPC Requirements

All SLURM scripts are written for the **Palmetto cluster (Clemson University)** but are straightforward to adapt for any SLURM-based HPC system. Adjust the following in each script as needed:

- `module load` — replace with your cluster's module names
- File paths under `/project/bcampb7/` — replace with your project directory
- `--mem` and `--cpus-per-task` — adjust to your allocation

---

## ⚡ Quick Start with Toy Data

If you want to run the full pipeline end-to-end with the included toy data before using your own:

```bash
# Clone the repository
git clone https://github.com/jojyjohn28/metatranscriptomics.git
cd metatranscriptomics

# Day 2: QC check on toy reads
bash day2-qc-preprocessing/00_check_pair_counts.slurm

# Day 3: The Bowtie2 index is pre-built in day1-toy-data/bowtie_index/
# No need to rebuild for the toy data

# Day 4: Map toy reads
sbatch day4-metatranscriptomics-mapping-quantification/01_bowtie2_mapping.slurm

# Day 5: Run DESeq2 with toy counts
Rscript day5-metatranscriptomics-counts-to-biology/Deseq_with_toydata.R
```

---

## 📊 Example Output

The file `day5-metatranscriptomics-counts-to-biology/Differential_abundance_cazy_all.pdf` shows an example volcano plot output from the DESeq2 analysis — faceted by environmental factor (Site, Season, Salinity) with top significant genes labeled.

---

## 📢 Notes

- Toy datasets are provided for learning purposes only — replace all input paths when running on your own data
- The pre-built Bowtie2 index in `day1-toy-data/bowtie_index/` is for the toy MAGs only; build your own index using Day 3 scripts for real data
- All SLURM array jobs use `--array=0-2` by default (3 samples); adjust the upper index to match your sample count
- Scripts are numbered to indicate the recommended run order within each day

---

## 🤝 Acknowledgment

I thank my PI, **Dr. Barbara Campbell**, for her continued support and guidance while learning and implementing metatranscriptomics analysis in our daily research.

---

## ⭐ Found this helpful?

Consider starring the repository or sharing it with others learning metatranscriptomics. Feedback, issues, and pull requests are welcome.

---

Last updated: April 8 2026
This series is completed.
