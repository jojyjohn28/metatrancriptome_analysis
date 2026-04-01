# 🧬 Day 2 — QC & Preprocessing for Metatranscriptomics

This folder contains all scripts used in **Day 2 of the Metatranscriptomics Series**, focusing on **RNA quality control and preprocessing**.

📘 Blog post: **[Day 2](https://jojyjohn28.github.io/blog/day2-metatranscriptomics-qc-preprocessing/)**

---

## 🎯 Goal of Day 2

Prepare raw RNA-seq reads for downstream analysis by:

- Assessing read quality
- Removing adapters and low-quality bases
- Eliminating ribosomal RNA (rRNA)
- Ensuring paired-end consistency

> **Garbage RNA in → garbage expression out**

---

## 📂 Workflow Overview

Raw Reads
↓
FastQC (diagnose quality)
↓
fastp (trim adapters & low-quality bases)
↓
SortMeRNA / bbduk (remove rRNA)
↓
Repair paired reads (if needed)
↓
Clean mRNA reads → ready for mapping

---

## 📜 Scripts in This Folder

### 🔹 00_check_pair_counts.slurm

Check whether R1 and R2 files have **equal read counts**.

- Detects pairing issues early
- Prevents mapping and counting errors

---

### 🔹 01_fastqc.slurm

Run **FastQC** on raw reads.

- Evaluates:
  - Per-base quality
  - Adapter contamination
  - Duplication levels

📌 Always run this before trimming

---

### 🔹 02_fastp_trim.slurm

Trim adapters and low-quality bases using **fastp**.

- Automatic adapter detection
- Quality filtering (Phred-based)
- Outputs:
  - Paired reads
  - Unpaired reads

---

### 🔹 03_sortmerna_rrna_removal.slurm

Remove rRNA reads using **SortMeRNA**.

- Uses SILVA and Rfam databases
- Outputs:
  - rRNA reads (for QC)
  - mRNA-enriched reads (used downstream)

📌 Recommended for accuracy

---

### 🔹 03b_bbduk_rrna_removal.slurm

Alternative rRNA removal using **bbduk (BBTools)**.

- Faster than SortMeRNA
- Suitable for large datasets
- Slightly less sensitive

📌 Choose this for speed or large-scale processing

---

### 🔹 04b_fastp_paired_sync.slurm

Ensure paired-end consistency during preprocessing.

- Prevents mismatches caused by trimming
- Routes orphan reads separately

📌 Preventive approach (recommended)

---

### 🔹 04_repair_bbtools.slurm

Repair mismatched paired-end reads using **repair.sh (BBTools)**.

- Synchronizes R1 and R2 files
- Outputs:
  - Repaired paired reads
  - Singleton reads

📌 Use only if mismatch is detected

---

## ⚠️ Critical Checks

Run pairing checks at two stages:

1. After **fastp**
2. After **rRNA removal**

Example:

```bash
# Count reads (lines / 4)
echo "R1:" $(($(wc -l < R1.fq)/4))
echo "R2:" $(($(wc -l < R2.fq)/4))
```

If counts differ:
→ run 04_repair_bbtools.slurm

📊 Expected Results

After preprocessing:

✅ Adapter sequences removed
✅ Improved read quality
✅ Significant reduction in dataset size (rRNA removal)
✅ Paired-end consistency maintained

Typical retention:
| Step | Retention |
| -------------- | --------- |
| Raw → Trimmed | ~95% |
| Trimmed → mRNA | ~5–40% |

💡 Notes
Scripts are designed for SLURM HPC environments
Modify file paths for your system
Toy dataset available in Day 1 folder

🚀 Next Step

➡️ Proceed to Day 3 — Reference Strategy

Your reference database defines your biological interpretation.
