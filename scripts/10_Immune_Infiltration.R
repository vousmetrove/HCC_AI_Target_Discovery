###############################################################################
# 10_Immune_Infiltration.R
# ---------------------------------------------------------------------------
# Purpose:
#   Characterize the tumor immune microenvironment in HCC and evaluate the
#   relationship between the identified hub genes and immune cell infiltration.
#   Tumor-infiltrating immune cells are key determinants of prognosis and
#   immunotherapy response; correlating hub genes with immune infiltration
#   provides mechanistic insight into their therapeutic relevance.
#
# Why immune infiltration for AI-assisted target discovery:
#   - Immune checkpoint inhibitors (anti-PD-1/PD-L1) are approved for HCC, but
#     response rates vary — understanding the immune context of targets is
#     critical for combination therapy strategies
#   - Hub genes associated with specific immune cell types may serve as
#     companion biomarkers for immunotherapy
#   - The immune microenvironment provides an orthogonal validation dimension:
#     genes with both prognostic and immune-modulatory roles are stronger
#     therapeutic candidates
#
# Analytical approach:
#   1. Estimate immune cell composition via ssGSEA (single-sample GSEA)
#      using immune cell type gene signatures from Bindea et al. (Immunity 2013)
#   2. GSVA enrichment of immune-related pathways (Hallmark, KEGG immune sets)
#   3. CIBERSORTx deconvolution (via web API or R implementation)
#   4. ESTIMATE: StromalScore, ImmuneScore, TumorPurity
#   5. Correlation analysis: hub genes vs. immune cell fractions
#   6. Differential immune infiltration: high- vs. low-risk groups
#   7. Generate heatmaps and correlation plots
#
# Output:
#   results/Immune/ssgsea_scores.rds          — per-sample immune cell scores
#   results/Immune/gsva_immune_pathways.rds   — per-sample pathway activity
#   results/Immune/estimate_scores.rds        — Stromal/Immune/TumorPurity
#   results/Immune/hub_immune_correlation.csv — hub gene–immune cell correlations
#   figures/immune_heatmap.pdf
#   figures/hub_immune_correlation_heatmap.pdf
#   figures/immune_boxplot_risk.pdf
#
# Dependencies: GSVA, GSEABase, estimate, immunedeconv, ggplot2, pheatmap
###############################################################################

# ---- 0. Environment ----
library(GSVA)              # ssGSEA & GSVA
library(GSEABase)          # gene set objects
library(limma)             # for expression matrix handling
library(estimate)          # ESTIMATE algorithm
library(immunedeconv)      # immune deconvolution (CIBERSORT, MCP-counter, etc.)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

set.seed(2024)
dir.create("../results/Immune", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load data ----
rlog_mat     <- readRDS("../results/DEG/rlog_normalized.rds")
metadata     <- readRDS("../data/processed/metadata.rds")
clinical     <- readRDS("../data/processed/clinical_clean.rds")
hub_genes_df <- tryCatch(readRDS("../results/PPI/ppi_hub_genes.rds"),
                         error = function(e) NULL)
consensus_df <- readRDS("../results/RandomForest/lasso_rf_consensus.rds")
surv_data    <- tryCatch(readRDS("../results/Survival/survival_data.rds"),
                         error = function(e) NULL)

# Tumor-only expression
tumor_bc  <- metadata %>% filter(tissue == "Tumor") %>% pull(barcode)
tumor_cols <- intersect(colnames(rlog_mat), tumor_bc)
expr_tumor <- rlog_mat[, tumor_cols, drop = FALSE]

message("Tumor expression matrix: ", nrow(expr_tumor), " genes × ",
        ncol(expr_tumor), " samples")

# ---- 2. Define gene sets for immune cell types ----
# Bindea et al. (Immunity, 2013) defined 24 immune cell type gene signatures.
# These are widely used for ssGSEA-based immune infiltration estimation.

bindea_signatures <- list(
  "aDC"           = c("CCL1", "EBF2", "CCL17", "P2RY14", "SLC25A33",
                       "CD1A", "CD1B", "CD1C", "CD1E", "FCER1A", "CLEC10A"),
  "B cells"       = c("CD19", "MS4A1", "CD79A", "CD79B", "BLK", "TNFRSF17",
                       "FCRL5", "BANK1", "PAX5", "IGJ", "IGHM"),
  "CD8 T cells"   = c("CD8A", "CD8B", "PRF1", "GZMA", "GZMB", "GZMH",
                       "GZMK", "GNLY", "NKG7", "CD247"),
  "Cytotoxic cells" = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK",
                        "GNLY", "NKG7", "FASLG", "KLRK1", "KLRD1"),
  "DC"            = c("NRP1", "ITGAX", "HLA-DRA", "HLA-DRB1", "HLA-DQA1",
                       "CLEC4C", "LILRA4", "IRF7", "IRF8", "TLR7"),
  "Eosinophils"   = c("IL5RA", "SIGLEC8", "CCR3", "RNASE2", "RNASE3",
                       "EPX", "PRG2", "PTGDR2", "HRH4"),
  "Macrophages"   = c("CD68", "CD163", "MSR1", "MRC1", "CSF1R",
                       "ITGAM", "CD14", "TLR2", "TLR4", "CD84"),
  "Mast cells"    = c("KIT", "CPA3", "TPSAB1", "TPSB2", "CMA1",
                       "CTSG", "HDC", "MS4A2", "FCER1A"),
  "Monocytes"     = c("CD14", "FCGR3A", "FCGR3B", "CSF1R", "ITGAM",
                       "CD33", "CD163", "FUT4", "VCAN", "S100A12"),
  "Neutrophils"   = c("ELANE", "MPO", "CTSG", "PRTN3", "CEACAM8",
                       "FCGR3B", "CXCR1", "CXCR2", "CSF3R", "FPR1"),
  "NK cells"      = c("NCR1", "KLRK1", "KLRD1", "KLRB1", "KIR2DL1",
                       "KIR2DL3", "KIR2DL4", "KIR3DL1", "KIR3DL2", "CD160"),
  "NK CD56dim"    = c("KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2",
                       "IL21R", "KIR2DS1", "KIR2DS2", "KIR2DS5"),
  "NK CD56bright" = c("NCAM1", "CD160", "CCR7", "SELL", "GZMK",
                       "XCL1", "XCL2", "IFNG"),
  "T cells"       = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6",
                       "CD28", "TRAC", "TRBC1", "TRBC2"),
  "T helper"      = c("CD4", "CD40LG", "IL21", "ICOS", "CXCL13",
                       "BCL6", "MAF", "SH2D1A", "CD84"),
  "Tcm"           = c("CCR7", "SELL", "CD27", "CD28", "IL7R",
                       "LEF1", "TCF7", "BCL2"),
  "Tem"           = c("KLRG1", "PRF1", "GZMA", "GZMB", "GZMH",
                       "CX3CR1", "TBX21", "EOMES", "IL2RB"),
  "TFH"           = c("CXCR5", "BCL6", "ICOS", "PDCD1", "IL21",
                       "MAF", "CD200", "BTLA", "SH2D1A"),
  "Tgd"           = c("TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRDV3",
                       "TRGV9", "CD160", "KLRB1"),
  "Th1"           = c("IFNG", "TBX21", "STAT1", "STAT4", "IL12RB2",
                       "CXCR3", "CCR5", "CD94"),
  "Th2"           = c("IL4", "IL5", "IL13", "GATA3", "STAT6",
                       "CCR3", "CCR4", "CCR8"),
  "Th17"          = c("IL17A", "IL17F", "RORC", "STAT3", "CCR6",
                       "IL22", "IL23R", "KLRB1"),
  "Treg"          = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18",
                       "IKZF2", "IKZF4", "CD25", "ENTPD1"),
  "Endothelial"   = c("PECAM1", "CDH5", "VWF", "ENG", "KDR",
                       "FLT1", "TEK", "TIE1", "SELE"),
  "Fibroblasts"   = c("COL1A1", "COL1A2", "COL3A1", "FAP",
                       "ACTA2", "PDGFRA", "PDGFRB", "S100A4")
)

# ---- 3. ssGSEA for immune cell estimation ----
# ssGSEA scores each sample for enrichment of each immune cell gene set,
# producing a continuous per-sample score that correlates with cell abundance.

# Filter gene sets to genes present in our expression matrix
bindea_filtered <- lapply(bindea_signatures, function(gs) {
  intersect(gs, rownames(expr_tumor))
})
bindea_filtered <- bindea_filtered[lengths(bindea_filtered) >= 3]

message("Immune gene sets with ≥ 3 genes in data: ", length(bindea_filtered))

# Run ssGSEA
ssgsea_res <- gsva(
  as.matrix(expr_tumor),
  bindea_filtered,
  method    = "ssgsea",
  kcdf      = "Gaussian",    # rlog values are continuous ~Normal
  min.sz    = 3,
  max.sz    = 500,
  ssgsea.norm = TRUE,        # normalize scores across samples
  verbose   = TRUE
)

saveRDS(ssgsea_res, "../results/Immune/ssgsea_scores.rds")

message("ssGSEA complete: ", nrow(ssgsea_res), " cell types × ",
        ncol(ssgsea_res), " samples")

# ---- 4. ESTIMATE: StromalScore, ImmuneScore, TumorPurity ----
# ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumor
# tissues using Expression data) uses gene expression signatures to infer
# the fraction of stromal and immune cells.

# ESTIMATE requires the estimate package to write intermediate files
estimate_dir <- "../results/Immune/estimate_tmp"
dir.create(estimate_dir, recursive = TRUE, showWarnings = FALSE)
owd <- getwd()
setwd(estimate_dir)

# Prepare input: gene symbols as rownames, log2 expression
estimate_input <- expr_tumor
rownames(estimate_input) <- rownames(expr_tumor)

# Write expression matrix for ESTIMATE
write.table(estimate_input, "expr_for_estimate.txt",
            sep = "\t", quote = FALSE, col.names = NA)

# Run ESTIMATE
estimate_scores <- tryCatch({
  estimateScore("expr_for_estimate.txt", "estimate_scores.txt",
                platform = "illumina")
  read.table("estimate_scores.txt", header = TRUE, row.names = 1,
             check.names = FALSE)
}, error = function(e) {
  message("ESTIMATE failed (may need different platform): ", e$message)
  NULL
})

setwd(owd)

if (!is.null(estimate_scores)) {
  # Transpose: rows = samples, cols = StromalScore, ImmuneScore, ESTIMATEScore
  estimate_df <- as.data.frame(t(estimate_scores)) %>%
    rownames_to_column("sample") %>%
    mutate(
      TumorPurity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
    )
  saveRDS(estimate_df, "../results/Immune/estimate_scores.rds")
  message("ESTIMATE complete: ", nrow(estimate_df), " samples")
}

# ---- 5. GSVA of immune-related Hallmark pathways ----
# Hallmark gene sets that capture immune and stromal biology

hallmark_immune <- list(
  HALLMARK_INFLAMMATORY_RESPONSE = c(
    "IL6", "TNF", "IL1B", "CXCL8", "CCL2", "CXCL10", "NFKB1", "RELA",
    "STAT3", "JUN", "FOS", "PTGS2", "ICAM1", "VCAM1", "SELE", "SELP",
    "IL1A", "CSF1", "CSF2", "CSF3", "CCL5", "CCL3", "CCL4", "CXCL1",
    "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CCL11", "CCL20", "IL12A", "IL12B",
    "IL15", "IL18", "IL23A", "IL2", "IL4", "IL7", "IL10", "IL13"
  ),
  HALLMARK_INTERFERON_GAMMA = c(
    "IFNG", "STAT1", "IRF1", "IRF9", "GBP1", "GBP2", "GBP4", "GBP5",
    "CXCL9", "CXCL10", "CXCL11", "CCR5", "CCL2", "CCL5", "CIITA", "HLA-DRA",
    "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74",
    "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "PSME1", "PSME2", "ICAM1"
  ),
  HALLMARK_INTERFERON_ALPHA = c(
    "IFNA1", "IFNA2", "IFNA4", "IFNB1", "STAT1", "STAT2", "IRF1", "IRF7",
    "IRF9", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "ISG15", "ISG20",
    "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFI35", "RSAD2", "DDX58", "EIF2AK2",
    "XAF1", "SP100", "IFI16", "BST2", "ADAR"
  ),
  HALLMARK_TNFA_SIGNALING = c(
    "TNF", "TNFRSF1A", "TNFRSF1B", "NFKB1", "NFKBIA", "NFKBIB",
    "RELA", "RELB", "JUN", "JUNB", "FOS", "FOSB", "ATF3", "CEBPB",
    "TNFAIP3", "TNFAIP2", "TNFAIP6", "TRAF1", "IER3", "NFKB2",
    "MAP3K8", "BIRC3", "CXCL1", "CXCL2", "CXCL3", "CCL2", "CCL5", "CCL20",
    "IL6", "IL1B", "ICAM1", "VCAM1", "SELE", "SELP", "SOD2"
  ),
  HALLMARK_IL6_JAK_STAT3 = c(
    "IL6", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS3", "CRP", "SAA1",
    "SAA2", "HP", "ORM1", "ORM2", "SERPINA3", "FGA", "FGB", "FGG",
    "ALB", "TFRC", "HAMP", "LCN2", "LBP", "CD14", "TLR2"
  ),
  HALLMARK_COMPLEMENT = c(
    "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C3AR1",
    "C4A", "C4B", "C5", "C5AR1", "C6", "C7", "C8A", "C8B", "C8G", "C9",
    "CFB", "CFD", "CFH", "CFI", "CR1", "CR2", "CD46", "CD55", "CD59",
    "SERPING1", "MBL2", "MASP1", "MASP2", "F2R", "F3", "PLAT", "PLAUR"
  ),
  HALLMARK_APOPTOSIS = c(
    "BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "BAD", "BID", "BIK",
    "BIM", "BBC3", "PMAIP1", "CASP3", "CASP6", "CASP7", "CASP8",
    "CASP9", "CASP10", "APAF1", "CYCS", "DIABLO", "HTRA2", "ENDOG",
    "AIFM1", "PARP1", "TP53", "CDKN1A", "FAS", "FASLG", "TNFRSF10A",
    "TNFRSF10B", "TNFSF10"
  )
)

# Filter to genes present in data
hallmark_filt <- lapply(hallmark_immune, function(gs) {
  intersect(gs, rownames(expr_tumor))
})
hallmark_filt <- hallmark_filt[lengths(hallmark_filt) >= 5]

gsva_hallmark <- gsva(
  as.matrix(expr_tumor),
  hallmark_filt,
  method = "gsva",
  kcdf   = "Gaussian",
  verbose = TRUE
)

saveRDS(gsva_hallmark, "../results/Immune/gsva_immune_pathways.rds")

message("GSVA immune pathways: ", nrow(gsva_hallmark), " pathways")

# ---- 6. Define candidate genes for immune correlation ----
# Use consensus genes (LASSO ∩ RF) as the primary candidate set

candidate_genes <- consensus_df$gene

# Add PPI hub genes if available
if (!is.null(hub_genes_df)) {
  ppi_hub_genes <- hub_genes_df$gene
  candidate_genes <- unique(c(candidate_genes, ppi_hub_genes))
}

message("Candidate genes for immune correlation: ", length(candidate_genes))

# ---- 7. Hub gene–immune cell correlation ----
# Spearman correlation between candidate gene expression and immune cell fractions

immune_scores <- as.data.frame(t(ssgsea_res))
immune_cell_types <- colnames(immune_scores)

candidate_expr <- expr_tumor[
  intersect(candidate_genes, rownames(expr_tumor)),
  ,
  drop = FALSE
]

cor_matrix <- matrix(NA_real_,
                     nrow = nrow(candidate_expr),
                     ncol = length(immune_cell_types),
                     dimnames = list(rownames(candidate_expr), immune_cell_types))

pval_matrix <- cor_matrix

for (gene in rownames(candidate_expr)) {
  gene_expr <- as.numeric(candidate_expr[gene, ])

  for (cell_type in immune_cell_types) {
    cell_scores <- immune_scores[[cell_type]]
    valid <- complete.cases(gene_expr, cell_scores)
    if (sum(valid) < 10) next

    ct <- cor.test(gene_expr[valid], cell_scores[valid],
                   method = "spearman", exact = FALSE)
    cor_matrix[gene, cell_type]  <- ct$estimate
    pval_matrix[gene, cell_type] <- ct$p.value
  }
}

# Convert to long format for saving
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "immune_cell", values_to = "rho")

pval_long <- as.data.frame(pval_matrix) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "immune_cell", values_to = "p_value")

hub_immune_cor <- cor_long %>%
  left_join(pval_long, by = c("gene", "immune_cell")) %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  arrange(p_value)

write.csv(hub_immune_cor, "../results/Immune/hub_immune_correlation.csv",
          row.names = FALSE)
saveRDS(hub_immune_cor, "../results/Immune/hub_immune_correlation.rds")

# Top significant correlations
sig_cors <- hub_immune_cor %>%
  filter(p_value < 0.05) %>%
  arrange(desc(abs(rho)))

message("Significant gene–immune correlations (P < 0.05): ", nrow(sig_cors))

# ---- 8. Differential immune infiltration: risk groups ----
if (!is.null(surv_data)) {
  immune_with_risk <- immune_scores %>%
    rownames_to_column("patient_id_tmp") %>%
    mutate(patient_id = substr(patient_id_tmp, 1, 12)) %>%
    left_join(
      surv_data %>% select(patient_id, risk_group),
      by = "patient_id"
    ) %>%
    filter(!is.na(risk_group))

  if (nrow(immune_with_risk) > 10 && length(unique(immune_with_risk$risk_group)) == 2) {

    # Wilcoxon test per immune cell type
    diff_immune <- lapply(immune_cell_types, function(ct) {
      high_vals <- immune_with_risk[[ct]][immune_with_risk$risk_group == "High"]
      low_vals  <- immune_with_risk[[ct]][immune_with_risk$risk_group == "Low"]

      if (length(high_vals) >= 3 && length(low_vals) >= 3) {
        wt <- wilcox.test(high_vals, low_vals)
        data.frame(
          cell_type = ct,
          mean_high = mean(high_vals, na.rm = TRUE),
          mean_low  = mean(low_vals, na.rm = TRUE),
          log2FC    = mean(high_vals, na.rm = TRUE) - mean(low_vals, na.rm = TRUE),
          p_value   = wt$p.value,
          stringsAsFactors = FALSE
        )
      }
    }) %>% bind_rows() %>% mutate(p_adj = p.adjust(p_value, "BH"))

    write.csv(diff_immune, "../results/Immune/immune_risk_diff.csv", row.names = FALSE)

    # Boxplot of top differential immune cells
    top_diff <- diff_immune %>%
      filter(p_adj < 0.05) %>%
      slice_max(abs(log2FC), n = 8)

    if (nrow(top_diff) > 0) {
      plot_data <- immune_with_risk %>%
        select(patient_id, risk_group, all_of(top_diff$cell_type)) %>%
        pivot_longer(-c(patient_id, risk_group),
                     names_to = "cell_type", values_to = "score")

      immune_box <- ggplot(plot_data, aes(x = risk_group, y = score, fill = risk_group)) +
        geom_boxplot(outlier.size = 0.8, width = 0.5) +
        facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c(Low = "#377EB8", High = "#E41A1C")) +
        labs(
          title    = "Differential Immune Infiltration by Risk Group",
          subtitle = paste0("Wilcoxon test, FDR < 0.05 | ",
                            nrow(top_diff), " cell types significant"),
          x        = "Risk Group",
          y        = "ssGSEA Enrichment Score",
          fill     = "Risk Group"
        ) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave("../figures/immune_boxplot_risk.pdf", immune_box, width = 12, height = 8)
    }
  }
}

# ---- 9. Heatmap: hub gene–immune cell correlation ----
# Select top 20 genes with the strongest immune associations

top_cor_genes <- hub_immune_cor %>%
  group_by(gene) %>%
  summarize(max_abs_rho = max(abs(rho), na.rm = TRUE)) %>%
  slice_max(max_abs_rho, n = 20) %>%
  pull(gene)

# Prepare matrix for heatmap
cor_heatmap_mat <- cor_matrix[top_cor_genes, , drop = FALSE]
# Remove rows/cols with all NA
cor_heatmap_mat <- cor_heatmap_mat[
  rowSums(is.na(cor_heatmap_mat)) < ncol(cor_heatmap_mat),
  colSums(is.na(cor_heatmap_mat)) < nrow(cor_heatmap_mat),
  drop = FALSE
]

if (nrow(cor_heatmap_mat) >= 5 && ncol(cor_heatmap_mat) >= 3) {
  pdf("../figures/hub_immune_correlation_heatmap.pdf", width = 10, height = 8)

  ht <- Heatmap(
    cor_heatmap_mat,
    name = "Spearman ρ",
    col  = colorRamp2(c(-0.8, 0, 0.8), c("#377EB8", "white", "#E41A1C")),
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    row_names_gp    = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    row_title       = "Candidate Genes",
    column_title    = "Immune Cell Types (ssGSEA)",
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- pval_matrix[rownames(cor_heatmap_mat)[i],
                        colnames(cor_heatmap_mat)[j]]
      if (!is.na(p) && p < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
      }
    }
  )
  draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
}

# ---- 10. Immune heatmap: all samples by cell type ----
# Z-score normalize across samples for visualization

immune_z <- t(scale(t(ssgsea_res)))

pdf("../figures/immune_heatmap.pdf", width = 12, height = 7)
Heatmap(
  immune_z,
  name                = "Z-score",
  col                 = colorRamp2(c(-3, 0, 3), c("#377EB8", "white", "#E41A1C")),
  show_column_names   = FALSE,
  column_title        = paste0("TCGA-LIHC Tumor Samples (n=",
                               ncol(immune_z), ")"),
  row_names_gp        = gpar(fontsize = 8),
  row_title           = "Immune Cell Type (ssGSEA)"
)
dev.off()

# ---- 11. Summary ----
cat("\n============================================================\n")
cat(" Immune infiltration analysis complete.\n")
cat(" Immune cell types quantified (ssGSEA):", nrow(ssgsea_res), "\n")
cat(" GSVA immune pathways:                 ", nrow(gsva_hallmark), "\n")
cat(" Gene–immune correlations (P < 0.05):  ", nrow(sig_cors), "\n")
cat("\n Output saved to results/Immune/\n")
cat("============================================================\n")
