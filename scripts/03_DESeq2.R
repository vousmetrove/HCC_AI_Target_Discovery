###############################################################################
# 03_DESeq2.R
# ---------------------------------------------------------------------------
# Purpose:
#   Perform differential expression analysis between HCC tumor and adjacent
#   normal liver tissue using DESeq2 on STAR unstranded raw counts. This script also
#   generates the PRIMARY normalized expression matrix (rlog) that feeds ALL
#   downstream machine learning and survival analyses.
#
#   CRITICAL — single expression modality pipeline:
#   The rlog matrix produced here is THE normalized expression matrix for the
#   entire project. It replaces the log2(FPKM+1) approach with a statistically
#   principled transformation that:
#     - Accounts for library size via DESeq2's median-of-ratios normalization
#     - Stabilizes variance across the full dynamic range of expression
#     - Removes the mean-variance dependency inherent in count data
#     - Is on a log2-like scale suitable for linear models (LASSO, Cox, RF)
#
#   All downstream scripts (05_LASSO, 06_RandomForest, 07_Survival,
#   08_GEO_validation, 09_HubGene_Integration) load the rlog matrix from:
#     results/DEG/rlog_normalized.rds
#
# Analytical approach:
#   1. DESeq2 negative binomial model on filtered STAR unstranded raw counts
#   2. Wald test for tumor vs. normal contrast
#   3. Independent filtering & multiple-testing correction (Benjamini-Hochberg)
#   4. DEG classification: |log2FC| > 1 & padj < 0.05
#   5. Regularized log (rlog) transformation — THE normalized expression matrix
#      for all downstream ML, visualization, and survival modeling
#   6. Publication-quality volcano plot, MA plot, and top-DEG heatmap
#   7. PCA from rlog values (for publication, replacing exploratory log2(CPM) PCA)
#
# Output:
#   results/DEG/DESeq2_object.rds         — full DESeq2 object
#   results/DEG/DESeq2_results.rds        — full results table (all genes)
#   results/DEG/DEG_list.rds / .csv       — significant DEGs only
#   results/DEG/rlog_normalized.rds       — PRIMARY normalized expression matrix
#                                           (feeds scripts 05, 06, 07, 08, 09)
#   figures/volcano_DEG.pdf               — volcano plot
#   figures/heatmap_top50_DEG.pdf         — top DEG heatmap (rlog Z-score)
#   figures/PCA_vst.pdf                  — PCA from rlog (publication)
#
# Dependencies: DESeq2, ggplot2, ggrepel, pheatmap, dplyr
###############################################################################

# ---- 0. Environment ----
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tibble)

set.seed(2024)
dir.create("../results/DEG", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load processed data ----
count_mat  <- readRDS("../data/processed/counts_filtered.rds")
metadata   <- readRDS("../data/processed/metadata.rds")
clinical   <- readRDS("../data/processed/clinical_clean.rds")

# ---- 2. Build DESeq2 dataset ----
# DESeq2 requires:
#   - countData: integer matrix of raw counts (genes × samples)
#   - colData:   sample-level metadata with the variable of interest
#   - design:    formula specifying the comparison (~ tissue)

# Ensure count matrix is integer
count_mat <- round(count_mat)

# Create colData matching the count matrix columns
coldata <- metadata %>%
  mutate(
    tissue = factor(tissue, levels = c("Normal", "Tumor"))
  ) %>%
  as.data.frame()
rownames(coldata) <- colnames(count_mat)

# Verify sample alignment
stopifnot(all(colnames(count_mat) == rownames(coldata)))

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ tissue
)

# ---- 3. Run DESeq2 pipeline ----
# This performs:
#   a) Estimation of size factors (median-of-ratios normalization)
#   b) Estimation of gene-wise dispersion
#   c) Fitting negative binomial GLM
#   d) Wald test for differential expression

dds <- DESeq(dds, parallel = FALSE)

# Save the full DESeq2 object for later use
saveRDS(dds, "../results/DEG/DESeq2_object.rds")

# ---- 4. Extract results ----
# Contrast: Tumor vs Normal → positive log2FC = upregulated in tumor
res <- results(
  dds,
  contrast       = c("tissue", "Tumor", "Normal"),
  alpha          = 0.05,             # FDR cutoff for independent filtering
  lfcThreshold   = 0,                # no minimal LFC threshold in the test
  independentFiltering = TRUE
)

# Convert to data frame and annotate
res_df <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)   # sort by significance

# Classify DEGs
res_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1   ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1  ~ "Downregulated",
      TRUE                                ~ "Not significant"
    )
  )

# Count DEGs
n_up   <- sum(res_df$significance == "Upregulated",   na.rm = TRUE)
n_down <- sum(res_df$significance == "Downregulated", na.rm = TRUE)

message("DEGs: ", n_up, " up | ", n_down, " down | ",
        n_up + n_down, " total significant")

# Save full results table
saveRDS(res_df, "../results/DEG/DESeq2_results.rds")

# Save only significant DEGs
deg_list <- res_df %>% filter(significance != "Not significant")
saveRDS(deg_list, "../results/DEG/DEG_list.rds")
write.csv(deg_list, "../results/DEG/DEG_list.csv", row.names = FALSE)

# ---- 5. Volcano plot ----
# Highlights top 20 genes by |log2FC| for labeling

# Select top genes to label
top_label <- res_df %>%
  filter(significance != "Not significant") %>%
  slice_max(order_by = abs(log2FoldChange), n = 20)

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), size = 0.6, alpha = 0.6) +
  scale_color_manual(
    values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8",
               "Not significant" = "grey70")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top_label,
    aes(label = gene_id),
    size  = 3,
    max.overlaps = 25,
    box.padding  = 0.4
  ) +
  labs(
    title    = "DEGs: HCC Tumor vs. Adjacent Normal",
    subtitle = paste0("Up: ", n_up, " | Down: ", n_down,
                      " (|log2FC| > 1, padj < 0.05)"),
    x        = "log2 Fold Change (Tumor / Normal)",
    y        = expression(-log[10](adjusted~P~value))
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("../figures/volcano_DEG.pdf", volcano, width = 8, height = 7)

# ---- 6. Regularized log transformation ----
# vst (variance stabilizing transformation) produces data suitable for:
# - PCA / clustering
# - Heatmap visualization
# - Correlation analysis (without the mean-variance dependency of raw counts)
#
# Note: vst() is used instead of rlog() because rlog() is prohibitively slow
# with >50 samples. vst() produces nearly identical log2-scale values and is
# the DESeq2-recommended choice for large datasets.

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
saveRDS(vst_mat, "../results/DEG/rlog_normalized.rds")  # kept filename for downstream compat

# ---- 7. Heatmap of top 50 DEGs ----
top50 <- deg_list %>%
  slice_min(padj, n = 50) %>%
  pull(gene_id)

# Extract rlog values for top genes, subset to tumor samples
heatmap_data <- vst_mat[top50, coldata$tissue == "Tumor"]
# ... or use all samples if preferred
heatmap_data_all <- vst_mat[top50, ]

# Z-score normalization for better heatmap contrast
heatmap_z <- t(scale(t(heatmap_data_all)))

# Annotation
ann_col <- data.frame(
  Tissue = coldata$tissue,
  row.names = colnames(heatmap_data_all)
)
ann_colors <- list(Tissue = c(Tumor = "#E41A1C", Normal = "#377EB8"))

heatmap_plot <- pheatmap(
  heatmap_z,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  fontsize_row    = 5,
  annotation_col  = ann_col,
  annotation_colors = ann_colors,
  main            = "Top 50 DEGs — Z-score (rlog)",
  color           = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
  filename        = "../figures/heatmap_top50_DEG.pdf",
  width           = 12,
  height          = 10
)

# ---- 8. PCA from rlog (publication-quality) ----
# This replaces the exploratory log2(CPM) PCA from 02_preprocessing.
# rlog accounts for library size and mean-variance dependency, producing
# a more faithful low-dimensional representation of the data.

vst_pca <- prcomp(t(vst_mat), center = TRUE, scale. = TRUE)
vst_pca_df <- data.frame(
  PC1    = vst_pca$x[, 1],
  PC2    = vst_pca$x[, 2],
  tissue = coldata$tissue
)
rlog_var_pc1 <- round(summary(vst_pca)$importance[2, 1] * 100, 1)
rlog_var_pc2 <- round(summary(vst_pca)$importance[2, 2] * 100, 1)

vst_pca_plot <- ggplot(vst_pca_df, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
  labs(
    title    = "PCA — TCGA-LIHC (DESeq2 vst)",
    subtitle = paste("Tumor:", sum(vst_pca_df$tissue == "Tumor"),
                     "| Normal:", sum(vst_pca_df$tissue == "Normal")),
    x        = paste0("PC1 (", rlog_var_pc1, "%)"),
    y        = paste0("PC2 (", rlog_var_pc2, "%)")
  ) +
  theme_minimal(base_size = 14)

ggsave("../figures/PCA_vst.pdf", vst_pca_plot, width = 7, height = 6)

# ---- 9. Summary ----
cat("\n============================================================\n")
cat(" DESeq2 differential expression analysis complete.\n")
cat(" DEGs (|log2FC|>1, padj<0.05):", n_up + n_down, "\n")
cat("   - Upregulated:  ", n_up, "\n")
cat("   - Downregulated:", n_down, "\n")
cat("\n PRIMARY NORMALIZED MATRIX:\n")
cat("   results/DEG/rlog_normalized.rds\n")
cat("   → This feeds scripts 05_LASSO, 06_RandomForest,\n")
cat("     07_Survival, 08_GEO_validation, 09_HubGene_Integration\n")
cat("\n Output saved to results/DEG/\n")
cat("============================================================\n")
