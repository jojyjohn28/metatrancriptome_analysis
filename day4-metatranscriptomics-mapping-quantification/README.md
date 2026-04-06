# 🧬 Day 4 — Mapping & Quantification

This folder contains SLURM scripts used in **Day 4 of the Metatranscriptomics Series**, where RNA reads are mapped to the reference and converted into gene expression counts.

📘 Blog post: **[Day 4](https://jojyjohn28.github.io/blog/day4-metatranscriptomics-mapping-quantification/)** :contentReference[oaicite:0]{index=0}

---

## 🎯 Goal of Day 4

Convert clean RNA reads into a **gene count table** by:

- Mapping reads to the reference (Bowtie2)
- Converting alignment files (SAM → BAM)
- Sorting and indexing alignments (SAMtools)
- Quantifying gene expression (featureCounts)

---

## 🔄 Workflow Overview

Clean mRNA Reads
↓
Bowtie2 Mapping
↓
SAM → BAM Conversion
↓
Sorting & Indexing
↓
featureCounts
↓
Gene Count Table

---

## 📜 Scripts in This Folder

### 🔹 `01_bowtie2_mapping.slurm`

Map paired-end RNA reads to the reference using **Bowtie2**.

- Uses indexed reference from Day 3
- Outputs alignment file (`.sam`)
- Includes alignment statistics (log file)

📌 Expect low alignment (~2%) when using small MAG references — this is normal

---

### 🔹 `02_sam_to_bam.slurm`

Convert and process alignment files using **SAMtools**.

- Converts `.sam` → `.bam`
- Sorts alignments
- Indexes BAM file

📌 Required for downstream quantification

---

### 🔹 `03_featurecounts.slurm`

Quantify gene expression using **featureCounts**.

- Uses GFF annotation
- Counts reads per gene
- Outputs count table for downstream analysis

📌 Output is used in Day 5 for normalization and interpretation

---

## ⚠️ Important Notes

- Low alignment rates are expected with limited references
- Zero counts do not indicate failure — they reflect reference coverage
- Ensure:
  - Paired reads are synchronized (Day 2)
  - Reference index exists (Day 3)

---

## 🚀 Next Step

➡️ Proceed to **Day 5 — From Counts to Biology**

This is where raw counts become biological insight.

---
