# 🧬 Day 3 — Reference Strategy: Your Database Defines Your Biology

This folder contains SLURM scripts used in **Day 3 of the Metatranscriptomics Series**, focusing on **building reference databases for mapping RNA reads**.

📘 Blog post: **[Day 3](https://jojyjohn28.github.io/blog/day3-metatranscriptomics-reference-strategy/)** :contentReference[oaicite:0]{index=0}

---

## 🎯 Goal of Day 3

Prepare a reference database for read mapping by:

- Understanding different reference strategies
- Building MAG-based references
- Preparing database-level references
- Creating gene catalog indices

> **You can only detect what exists in your reference**

---

## 🔄 Workflow Overview

Clean mRNA Reads
↓
Choose Reference Strategy
↓
MAGs / Genomes / Gene Catalog
↓
Build Bowtie2 Index
↓
Reference Ready → Mapping (Day 4)

---

## 🧠 Core Concept

Mapping is not just a technical step — it is a **biological filter**.

- Unmapped reads ≠ no expression
- Unmapped reads = **not represented in your reference**

Your results are only as complete as your reference.

---

## 📜 Scripts in This Folder

### 🔹 00_inspect_reference.slurm

Inspect input reference files before building indices.

- Counts contigs
- Checks file integrity
- Helps estimate reference coverage

📌 Always run before index building

---

### 🔹 01_build_mag_index.slurm

Build a **Bowtie2 index from MAGs**.

- Concatenates MAG FASTA files
- Creates `master_mags.fa`
- Builds Bowtie2 index

📌 High specificity, low coverage  
📌 Used in this series

---

### 🔹 02_build_db_index.slurm

Build index from **reference genome databases** (e.g., GTDB, RefSeq).

- Handles large FASTA files
- Uses `--large-index` when needed

📌 Medium coverage, medium specificity  
📌 Useful for broader taxonomic profiling

---

### 🔹 03_build_catalog_index.slurm

Build index from **gene catalogs**.

- Uses non-redundant gene collections
- Optimized for functional profiling

📌 High coverage, function-focused  
📌 Limited taxonomic resolution

---

## ⚖️ Reference Strategy Comparison

| Strategy          | Coverage | Specificity | Best Use             |
| ----------------- | -------- | ----------- | -------------------- |
| MAGs              | Low      | High        | Target organisms     |
| Reference genomes | Medium   | Medium      | Broad surveys        |
| Gene catalogs     | High     | Low–Medium  | Functional profiling |

---

## 📊 Expected Outcomes

After this step:

- ✅ A valid Bowtie2 index (6 `.bt2` files)
- ✅ Reference ready for mapping
- ✅ Understanding of what your dataset can and cannot detect

---

## 💡 Notes

- MAG-based reference (used here) will result in **low alignment (~2%)**
- This is expected due to limited community coverage
- This is a **biological limitation**, not a pipeline failure

---

## 📜 Scripts in This Folder

### 🔹 `00_inspect_reference.slurm`

Inspect reference FASTA files before building any index.

This script is used to:

- verify that input reference files exist
- count the number of contigs/sequences
- check file size and general structure
- confirm that the reference is ready for indexing

📌 Use this first to sanity-check your reference before running Bowtie2 index building.

---

### 🔹 `01_build_mag_index.slurm`

Build a Bowtie2 index from a small set of MAGs.

This script is used to:

- combine selected MAG FASTA files into a single master reference
- generate a `master_mags.fa` file
- build the Bowtie2 index for MAG-based mapping

This is the reference strategy used in this tutorial series because it gives:

- high specificity to your system
- organism-resolved expression
- a simple and fast index for teaching purposes

📌 Best when you want to map reads to genomes reconstructed directly from your own metagenomes.

---

### 🔹 `02_build_db_index.slurm`

Build a Bowtie2 index from a larger genome reference database.

This script is used to:

- prepare and index a genome collection such as GTDB, RefSeq, or another curated reference set
- handle larger input FASTA files
- optionally use Bowtie2 large-index settings when required

This strategy gives:

- broader taxonomic coverage than a small MAG set
- improved alignment rates in many cases
- lower specificity than sample-derived MAGs

📌 Best when you want broader community coverage beyond your own reconstructed genomes.

---

### 🔹 `03_build_catalog_index.slurm`

Build a Bowtie2 index from a gene catalog.

This script is used to:

- prepare a non-redundant gene catalog FASTA
- build an index for function-focused mapping
- support expression analysis at the gene level rather than genome level

This strategy gives:

- the broadest functional coverage
- strong detection of expressed genes across diverse communities
- weaker direct taxonomic resolution unless the catalog is already annotated with source genomes

📌 Best when your main question is about **which functions are active**, rather than exactly **which organism is expressing them**.

---

## 🚀 Next Step

➡️ Proceed to **Day 4 — Mapping & Quantification**

This is where reads are aligned and expression is quantified.

---
