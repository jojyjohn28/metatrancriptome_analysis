# Metatranscriptomics — From Raw RNA → Biology 🧬

### A 5-Day Metatranscriptomics Learning Series

This repository accompanies a **5-day blog series on metatranscriptomics**, guiding you from raw RNA sequencing reads to biological interpretation of microbial community activity. Each day focuses on a key step in the workflow, with detailed blog explanations, reproducible scripts, and toy datasets for quick testing.

📘 **Full blog series:** https://jojyjohn28.github.io/series/metatranscriptome/

---

## 🚀 What You Will Learn

By the end of this series, you will be able to:

- Perform RNA-seq quality control and rRNA removal
- Choose the right reference strategy (MAGs vs. reference genomes vs. gene catalogs)
- Map metatranscriptomic reads with Bowtie2 and process alignments with SAMtools
- Quantify gene expression using featureCounts
- Normalize and interpret expression data in an ecological context (TPM, DESeq2, DNA:RNA ratios)

---

## 📁 Repository Structure

```
.
├── Day1_Introduction/
│   ├── toy_data/
│   └── README.md
├── Day2_QC_Preprocessing/
│   ├── scripts/
│   └── README.md
├── Day3_Reference_Strategy/
│   ├── scripts/
│   └── README.md
├── Day4_Mapping_Quantification/
│   ├── scripts/
│   └── README.md
├── Day5_Counts_to_Biology/
│   ├── scripts/
│   └── README.md
└── README.md
```

> Toy datasets are provided for learning purposes. Replace input paths with your own data when adapting to your project.

---

## 🧰 Tools Used

| Tool                    | Purpose                              |
| ----------------------- | ------------------------------------ |
| FastQC / fastp          | Quality control and adapter trimming |
| SortMeRNA / bbduk       | rRNA removal                         |
| Bowtie2                 | Read mapping                         |
| SAMtools                | SAM/BAM processing                   |
| featureCounts (Subread) | Gene-level quantification            |
| DESeq2 (R)              | Differential expression analysis     |

---

## 📢 Notes

- All paths are simplified for clarity; adapt to your HPC directory structure as needed
- SLURM batch scripts are included for HPC users (Palmetto cluster examples)
- This is an actively maintained repository — scripts and data will be updated as the series progresses

---

## 🤝 Acknowledgment

I thank my PI, **Dr. Barbara Campbell**, for her continued support and guidance while learning and implementing metatranscriptomics analysis in our daily research.

---

## ⭐ Found this helpful?

Consider starring the repository or sharing it with others learning metatranscriptomics. Feedback and contributions are welcome.

---

_Last updated: March 30, 2026_
