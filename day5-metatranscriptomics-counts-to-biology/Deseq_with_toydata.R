# Toy DESeq2 example (CAZy-like gene table)
# Includes: PCA + sample-distance heatmap + volcano + faceted volcano (top 5 labels)

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(pheatmap)

set.seed(1)

########################
# 1) Create TOY count matrix (genes x samples)
########################
ngenes <- 80
nsamp  <- 12

samples <- paste0("S", 1:nsamp)
genes <- paste0("GH", sample(1:50, ngenes, replace=TRUE), "_toy(", sample(LETTERS, ngenes, TRUE), ")")

countData <- matrix(rnbinom(ngenes*nsamp, mu=50, size=1), nrow=ngenes)
rownames(countData) <- genes
colnames(countData) <- samples

########################
# 2) Create TOY metadata (samples x factors)
########################
colData <- data.frame(
  SampleID  = samples,
  Bay       = rep(c("Chesapeake","Delaware"), each = 6),
  Season    = rep(c("Spring","Summer"), times = 6),
  Salinity  = rep(c("Low","High"), times = 6),
  stringsAsFactors = FALSE
)
rownames(colData) <- colData$SampleID
colData$Bay <- factor(colData$Bay)
colData$Season <- factor(colData$Season)
colData$Salinity <- factor(colData$Salinity)

# reference levels
colData$Bay <- relevel(colData$Bay, "Chesapeake")
colData$Season <- relevel(colData$Season, "Spring")
colData$Salinity <- relevel(colData$Salinity, "Low")

########################
# 3) Add a small "true signal" so plots look realistic
########################
# 10 genes respond to Salinity (High up)
sal_genes <- 1:10
countData[sal_genes, colData$Salinity == "High"] <- countData[sal_genes, colData$Salinity == "High"] + 80

# 8 genes respond to Bay (Delaware up)
bay_genes <- 11:18
countData[bay_genes, colData$Bay == "Delaware"] <- countData[bay_genes, colData$Bay == "Delaware"] + 60

# 8 genes respond to Season (Summer up)
sea_genes <- 19:26
countData[sea_genes, colData$Season == "Summer"] <- countData[sea_genes, colData$Season == "Summer"] + 60

########################
# 4) QC: VST -> PCA + sample-distance heatmap
########################
dds_qc <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds_qc, blind = TRUE)

# PCA
pca_df <- plotPCA(vsd, intgroup = c("Season","Bay","Salinity"), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
ggplot(pca_df, aes(PC1, PC2, color = Season, shape = Bay)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA (toy data, VST)")

# Sample distance heatmap
sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
pheatmap(sampleDistMatrix,
         annotation_col = colData[, c("Season","Bay","Salinity"), drop=FALSE],
         main = "Sample distances (toy data, VST)")

########################
# 5) DESeq2 for each factor (one at a time)
########################
run_one <- function(factor_name) {
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData,
                                design = as.formula(paste0("~ ", factor_name)))
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  res$Gene <- rownames(res)
  res$Factor <- factor_name
  res
}

res_all <- bind_rows(
  run_one("Salinity"),
  run_one("Season"),
  run_one("Bay")
) %>% filter(!is.na(padj))

########################
# 6) Volcano plot by factor (all genes)
########################
ggplot(res_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  facet_wrap(~ Factor, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Differential abundance by factor (toy data)",
       x = "Log2 fold change", y = "-Log10 adjusted p-value") +
  theme(legend.position = "bottom")

########################
# 7) Faceted volcano with top 5 significant labels per factor
########################
res_all_clean <- res_all %>%
  mutate(
    ShortName = str_extract(Gene, "(?<=\\().+?(?=\\))"),
    ShortName = ifelse(is.na(ShortName), Gene, ShortName)
  )

top_genes <- res_all_clean %>%
  filter(padj < 0.05) %>%
  group_by(Factor) %>%
  slice_min(order_by = padj, n = 5) %>%
  ungroup()

ggplot(res_all_clean, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = ShortName),
                  size = 3, max.overlaps = 50, segment.color = "grey50") +
  facet_wrap(~ Factor, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw() +
  labs(title = "Faceted volcano (toy data) — top 5 labels per factor",
       x = "Log2 fold change", y = "-Log10 adjusted p-value",
       color = "padj < 0.05") +
  theme(legend.position = "bottom")
#########################################
#end