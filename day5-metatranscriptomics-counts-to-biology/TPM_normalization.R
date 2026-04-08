# ── TPM calculation from featureCounts output ─────────────────

library(dplyr)
library(tidyr)
library(ggplot2)

# Load the featureCounts table
# Skip the first comment line (starts with #)
counts_raw <- read.table("counts/all_samples_counts.txt",
                         header = TRUE,
                         skip = 1,
                         sep = "\t",
                         check.names = FALSE)

# Extract gene lengths and count matrix
gene_ids    <- counts_raw$Geneid
gene_length <- counts_raw$Length

# Sample columns start at column 7 in featureCounts output
count_matrix <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(count_matrix) <- gene_ids

# ── Calculate TPM ──────────────────────────────────────────────
# Step 1: RPK = counts / (length in kb)
RPK <- sweep(count_matrix, 1, gene_length / 1000, FUN = "/")

# Step 2: TPM = RPK / (sum of RPK per sample) * 1e6
TPM <- sweep(RPK, 2, colSums(RPK) / 1e6, FUN = "/")

# Verify: each column should sum to ~1,000,000
colSums(TPM)

# Save TPM table
write.csv(as.data.frame(TPM), "counts/TPM_matrix.csv")