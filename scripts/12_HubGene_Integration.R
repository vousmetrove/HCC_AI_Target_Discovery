###############################################################################
# 12_HubGene_Integration.R
# ---------------------------------------------------------------------------
# Purpose:
#   Integrate multi-dimensional evidence from the complete AI-assisted pipeline
#   to prioritize prognostic therapeutic targets in HCC. This final integration
#   step consolidates findings across NINE orthogonal evidence dimensions and
#   generates the definitive ranked candidate target list for clinical translation.
#
# Integration dimensions (9 evidence sources):
#   Step 03 — DEG:            |log2FC| & padj significance
#   Step 04 — GO/KEGG:         cancer pathway / GO term membership
#   Step 05 — PPI (STRING):    network hub score (degree, betweenness, MCC)
#   Step 06 — LASSO:           penalized Cox selection (non-zero coefficient)
#   Step 07 — Random Forest:   permutation VIMP ranking (ranger)
#   Step 08 — Survival:        Cox HR, KM log-rank, time-ROC AUC
#   Step 09 — GEO Validation:  external cohort HR & significance
#   Step 10 — Immune:          immune cell correlation strength
#   Step 11 — Drug Repurposing: DGIdb druggability & polypharmacology score
#
# Output:
#   results/HubGene_Integration/final_targets.csv       — ranked target list
#   results/HubGene_Integration/final_targets.rds       — full object
#   results/HubGene_Integration/evidence_matrix.csv     — per-dimension scores
#   results/HubGene_Integration/integration_summary.txt — readable report
#   figures/target_prioritization_heatmap.pdf           — evidence heatmap
#   figures/target_prioritization_bar.pdf               — composite score bar
#   manuscript/supplementary_table_targets.csv          — for publication
#
# Dependencies: dplyr, tidyr, ggplot2, ComplexHeatmap, circlize
###############################################################################

# ---- 0. Environment ----
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

set.seed(2024)
dir.create("../results/HubGene_Integration", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load evidence from all upstream analyses ----
# Each upstream step contributes one or more scores. We load them all and
# build a unified evidence matrix for each candidate gene.

# Step 03 — DEG results
deg_all  <- readRDS("../results/DEG/DESeq2_results.rds")
deg_list <- readRDS("../results/DEG/DEG_list.rds")

# Step 04 — GO/KEGG enrichment
go_enrich  <- readRDS("../results/GO_KEGG/GO_enrichment_simplified.rds")
kegg_enrich <- readRDS("../results/GO_KEGG/KEGG_enrichment.rds")

# Step 05 — PPI network
ppi_metrics <- tryCatch(readRDS("../results/PPI/ppi_node_metrics.rds"),
                        error = function(e) NULL)
ppi_hub     <- tryCatch(readRDS("../results/PPI/ppi_hub_genes.rds"),
                        error = function(e) NULL)

# Step 06 — LASSO
lasso_genes <- readRDS("../results/LASSO/lasso_genes.rds")

# Step 07 — Random Forest
rf_importance <- readRDS("../results/RandomForest/rf_importance.rds")
consensus_df  <- readRDS("../results/RandomForest/lasso_rf_consensus.rds")
altmann_pvals <- tryCatch(readRDS("../results/RandomForest/rf_altmann_pvals.rds"),
                          error = function(e) NULL)

# Step 08 — Survival
surv_data  <- readRDS("../results/Survival/survival_data.rds")
cox_univar <- readRDS("../results/Survival/cox_univariate.rds")

# Step 09 — GEO Validation
geo_val <- tryCatch(readRDS("../results/GEO_validation/validation_summary.rds"),
                    error = function(e) NULL)

# Step 10 — Immune infiltration
immune_cor <- tryCatch(readRDS("../results/Immune/hub_immune_correlation.rds"),
                       error = function(e) NULL)

# Step 11 — Drug repurposing
gene_drugs <- tryCatch(readRDS("../results/Drug/gene_drug_count.rds"),
                       error = function(e) NULL)

# ---- 2. Define candidate gene pool ----
# All genes that appeared in at least one upstream selection step

candidates <- unique(c(
  deg_list$gene_id,                                # all DEGs
  lasso_genes$gene_symbol,                         # LASSO selected
  rf_importance$gene[1:min(100, nrow(rf_importance))],  # top RF
  if (!is.null(ppi_hub)) ppi_hub$gene else NULL,   # PPI hubs
  if (!is.null(immune_cor)) unique(immune_cor$gene[immune_cor$p_value < 0.01]) else NULL  # immune-correlated
))
candidates <- candidates[!is.na(candidates)]

message("Total candidate genes for integration: ", length(candidates))

# ---- 3. Build evidence matrix ----
# Each dimension normalized to 0–1 (1 = strongest evidence)

E <- data.frame(gene = candidates, stringsAsFactors = FALSE)

# ---- 3a. DEG evidence ----
deg_sub <- deg_all %>%
  filter(gene_id %in% candidates) %>%
  select(gene_id, log2FoldChange, padj, significance)

E <- E %>%
  left_join(deg_sub, by = c("gene" = "gene_id")) %>%
  mutate(
    deg_score = -log10(pmin(padj, 1e-300, na.rm = TRUE)) / -log10(1e-300),
    deg_score = ifelse(is.na(deg_score), 0, deg_score)
  )

# ---- 3b. GO/KEGG pathway evidence ----
cancer_pathways <- c("hsa05200", "hsa05225", "hsa04110", "hsa04151",
                     "hsa04310", "hsa04350", "hsa04010", "hsa04210",
                     "hsa04630", "hsa04115", "hsa04510")

enriched_genes <- c()
if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
  kegg_df <- as.data.frame(kegg_enrich)
  cancer_kegg <- kegg_df %>%
    filter(ID %in% cancer_pathways) %>%
    pull(geneID) %>%
    strsplit("/") %>%
    unlist() %>% unique()
  enriched_genes <- unique(c(enriched_genes, cancer_kegg))
}
if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
  go_df <- as.data.frame(go_enrich)
  cancer_go <- go_df %>%
    filter(grepl("cell cycle|apoptosis|proliferation|angiogenesis|metastasis|EMT|invasion|immune|inflammation|DNA repair|mitotic|hepatocyte|liver",
                 Description, ignore.case = TRUE)) %>%
    pull(geneID) %>%
    strsplit("/") %>%
    unlist() %>% unique()
  enriched_genes <- unique(c(enriched_genes, cancer_go))
}

E <- E %>%
  mutate(
    in_cancer_pathway = gene %in% enriched_genes,
    pathway_score = ifelse(in_cancer_pathway, 1, 0)
  )

# ---- 3c. PPI network evidence ----
if (!is.null(ppi_metrics)) {
  ppi_sub <- ppi_metrics %>%
    select(gene, hub_score, degree, betweenness) %>%
    mutate(
      ppi_score = hub_score,  # already 0–1 normalized composite
      ppi_degree_norm = pmin(degree / max(degree, na.rm = TRUE), 1)
    )
  E <- E %>%
    left_join(ppi_sub %>% select(gene, ppi_score, ppi_degree_norm),
              by = "gene") %>%
    mutate(
      ppi_score = ifelse(is.na(ppi_score), 0, ppi_score),
      ppi_degree_norm = ifelse(is.na(ppi_degree_norm), 0, ppi_degree_norm)
    )
} else {
  E <- E %>% mutate(ppi_score = 0, ppi_degree_norm = 0)
}

# ---- 3d. LASSO evidence ----
E <- E %>%
  mutate(
    in_lasso   = gene %in% lasso_genes$gene_symbol,
    lasso_coef = lasso_genes$coefficient[match(gene, lasso_genes$gene_symbol)],
    lasso_coef = ifelse(is.na(lasso_coef), 0, lasso_coef),
    lasso_score = ifelse(in_lasso, 1, 0)
  )

# ---- 3e. Random Forest evidence ----
rf_sub <- rf_importance %>%
  filter(importance > 0) %>%
  mutate(rf_score = pmin(importance / max(importance, na.rm = TRUE), 1))

E <- E %>%
  left_join(rf_sub %>% select(gene, rf_score), by = "gene") %>%
  mutate(rf_score = ifelse(is.na(rf_score), 0, rf_score))

# Altmann P-value bonus
if (!is.null(altmann_pvals)) {
  E <- E %>%
    left_join(altmann_pvals %>% select(gene, altmann_p = pval), by = "gene") %>%
    mutate(altmann_sig = ifelse(!is.na(altmann_p) & altmann_p < 0.05, 1, 0))
} else {
  E <- E %>% mutate(altmann_sig = 0)
}

# Consensus (LASSO ∩ RF)
E <- E %>%
  mutate(
    in_consensus = gene %in% consensus_df$gene,
    consensus_score = ifelse(in_consensus, 1, 0.5)
  )

# ---- 3f. Survival evidence ----
# Individual gene Cox regression
rlog_mat  <- readRDS("../results/DEG/rlog_normalized.rds")
metadata  <- readRDS("../data/processed/metadata.rds")
clinical  <- readRDS("../data/processed/clinical_clean.rds")
tumor_bc  <- metadata %>% filter(tissue == "Tumor") %>% pull(barcode)
tumor_cols <- intersect(colnames(rlog_mat), tumor_bc)
expr_tumor <- rlog_mat[intersect(candidates, rownames(rlog_mat)), tumor_cols, drop = FALSE]

patient_ids <- substr(colnames(expr_tumor), 1, 12)
clin_sub <- clinical %>% filter(submitter_id %in% patient_ids)
m <- match(patient_ids, clin_sub$submitter_id)
clin_sub <- clin_sub[m, ]

gene_cox <- data.frame(gene = character(), HR = numeric(),
                       p_value = numeric(), stringsAsFactors = FALSE)

for (g in intersect(candidates, rownames(expr_tumor))) {
  expr_val <- as.numeric(expr_tumor[g, ])
  valid <- complete.cases(expr_val, clin_sub$os_time, clin_sub$os_event) &
           clin_sub$os_time > 0
  if (sum(valid) < 30 || sd(expr_val[valid], na.rm = TRUE) < 1e-6) next

  cox_g <- tryCatch(
    coxph(Surv(clin_sub$os_time[valid], clin_sub$os_event[valid]) ~ expr_val[valid]),
    error = function(e) NULL
  )
  if (is.null(cox_g)) next

  s_g <- summary(cox_g)
  gene_cox <- rbind(gene_cox, data.frame(
    gene    = g,
    HR      = s_g$conf.int[1, "exp(coef)"],
    p_value = s_g$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  ))
}

E <- E %>%
  left_join(gene_cox, by = "gene") %>%
  mutate(
    survival_score = ifelse(!is.na(HR),
                            -log10(pmin(p_value, 0.99)) / -log10(0.01), 0),
    survival_score = pmin(pmax(survival_score, 0), 1)
  )

# ---- 3g. GEO validation evidence ----
if (!is.null(geo_val)) {
  n_sig_cohorts <- sum(geo_val$P_logrank < 0.05, na.rm = TRUE)
  E <- E %>% mutate(geo_score = ifelse(n_sig_cohorts >= 1, 1, 0.3))
} else {
  E <- E %>% mutate(geo_score = 0.5)  # neutral if no validation available
}

# ---- 3h. Immune correlation evidence ----
if (!is.null(immune_cor)) {
  immune_gene_score <- immune_cor %>%
    group_by(gene) %>%
    summarize(
      max_abs_rho = max(abs(rho), na.rm = TRUE),
      n_sig_immune = sum(p_value < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(immune_score = pmin(max_abs_rho, 1))

  E <- E %>%
    left_join(immune_gene_score %>% select(gene, immune_score, n_sig_immune),
              by = "gene") %>%
    mutate(
      immune_score = ifelse(is.na(immune_score), 0, immune_score),
      n_sig_immune = ifelse(is.na(n_sig_immune), 0, n_sig_immune)
    )
} else {
  E <- E %>% mutate(immune_score = 0, n_sig_immune = 0)
}

# ---- 3i. Drug / Druggability evidence ----
if (!is.null(gene_drugs)) {
  E <- E %>%
    left_join(gene_drugs %>% select(gene, n_drugs), by = "gene") %>%
    mutate(
      n_drugs = ifelse(is.na(n_drugs), 0, n_drugs),
      drug_score = pmin(n_drugs / max(n_drugs, 1), 1)
    )
} else {
  E <- E %>% mutate(n_drugs = 0, drug_score = 0)
}

# ---- 4. Composite score ----
# Weighted sum of all 9 evidence dimensions. Weights reflect the relative
# strength and orthogonality of each evidence type.

weights <- list(
  deg        = 0.10,   # differential expression
  pathway    = 0.05,   # functional enrichment
  ppi        = 0.10,   # protein network topology
  lasso      = 0.12,   # penalized regression selection
  rf         = 0.12,   # ensemble feature importance
  consensus  = 0.06,   # ML method agreement
  survival   = 0.12,   # prognostic validation
  geo        = 0.08,   # external cohort validation
  immune     = 0.10,   # immune microenvironment
  drug       = 0.15    # translational druggability (highest weight)
)

stopifnot(abs(sum(unlist(weights)) - 1) < 0.001)

E <- E %>%
  mutate(
    composite_score =
      deg_score       * weights$deg +
      pathway_score   * weights$pathway +
      ppi_score       * weights$ppi +
      lasso_score     * weights$lasso +
      rf_score        * weights$rf +
      consensus_score * weights$consensus +
      survival_score  * weights$survival +
      geo_score       * weights$geo +
      immune_score    * weights$immune +
      drug_score      * weights$drug,

    composite_100 = round(composite_score * 100, 1),

    tier = case_when(
      composite_score >= quantile(composite_score, 0.90, na.rm = TRUE) ~
        "Tier 1 — High Priority",
      composite_score >= quantile(composite_score, 0.70, na.rm = TRUE) ~
        "Tier 2 — Moderate Priority",
      composite_score >= quantile(composite_score, 0.40, na.rm = TRUE) ~
        "Tier 3 — Lower Priority",
      TRUE ~ "Tier 4 — Candidate"
    )
  ) %>%
  arrange(desc(composite_score))

# ---- 5. Final target table ----
final_targets <- E %>%
  select(
    gene, tier, composite_100,
    log2FoldChange, significance,
    deg_score, pathway_score, ppi_score,
    lasso_score, rf_score, consensus_score,
    survival_score, geo_score, immune_score, drug_score,
    HR, p_value = p_value,
    n_drugs, n_sig_immune
  ) %>%
  arrange(desc(composite_100))

write.csv(final_targets, "../results/HubGene_Integration/final_targets.csv",
          row.names = FALSE)
saveRDS(final_targets, "../results/HubGene_Integration/final_targets.rds")

# Evidence matrix (full)
write.csv(E, "../results/HubGene_Integration/evidence_matrix.csv", row.names = FALSE)

# Supplementary table for manuscript
dir.create("../manuscript", recursive = TRUE, showWarnings = FALSE)
write.csv(final_targets, "../manuscript/supplementary_table_targets.csv",
          row.names = FALSE)

# ---- 6. Visualization ----
# 6a. Composite score bar chart (top 20)

top_20 <- final_targets %>% slice_head(n = 20)

bar_plot <- ggplot(top_20, aes(x = reorder(gene, composite_100),
                                y = composite_100)) +
  geom_col(aes(fill = tier), width = 0.7) +
  scale_fill_manual(
    values = c(
      "Tier 1 — High Priority"     = "#B2182B",
      "Tier 2 — Moderate Priority" = "#F4A582",
      "Tier 3 — Lower Priority"    = "#92C5DE",
      "Tier 4 — Candidate"         = "#E0E0E0"
    )
  ) +
  coord_flip() +
  labs(
    title    = "AI-Assisted Target Discovery — Top 20 Prioritized Genes",
    subtitle = paste0("Composite score: 9 dimensions ",
                      "(DEG + GO/KEGG + PPI + LASSO + RF + Survival + GEO + Immune + Drug)"),
    x        = "",
    y        = "Composite Evidence Score (0–100)",
    fill     = "Priority Tier"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("../figures/target_prioritization_bar.pdf", bar_plot, width = 12, height = 7)

# 6b. Evidence heatmap (top 30 genes)

top_30 <- final_targets %>% slice_head(n = 30)

heat_cols <- c("deg_score", "pathway_score", "ppi_score", "lasso_score",
               "rf_score", "consensus_score", "survival_score", "geo_score",
               "immune_score", "drug_score")
heat_mat <- E %>%
  filter(gene %in% top_30$gene) %>%
  select(gene, all_of(heat_cols)) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

colnames(heat_mat) <- c("DEG", "Pathway", "PPI", "LASSO",
                         "RF", "Consensus", "Survival", "GEO",
                         "Immune", "Drug")

if (nrow(heat_mat) >= 5) {
  pdf("../figures/target_prioritization_heatmap.pdf", width = 11, height = 9)

  ht <- Heatmap(
    heat_mat,
    name               = "Evidence",
    col                = colorRamp2(c(0, 0.5, 1), c("#F7F7F7", "#FDB863", "#B2182B")),
    cluster_rows       = TRUE,
    cluster_columns    = FALSE,
    show_row_dend      = TRUE,
    row_names_gp       = gpar(fontsize = 8),
    column_names_gp    = gpar(fontsize = 9, fontface = "bold"),
    row_title          = "Candidate Genes",
    column_title       = "AI-Assisted Target Discovery — Evidence Matrix",
    cell_border        = "grey90",
    rect_gp            = gpar(col = "grey90", lwd = 0.5),
    heatmap_legend_param = list(title = "Score", at = c(0, 0.5, 1),
                                labels = c("0", "0.5", "1")),
    # Annotate tiers
    right_annotation = rowAnnotation(
      Tier = top_30$tier[match(rownames(heat_mat), top_30$gene)],
      col   = list(Tier = c(
        "Tier 1 — High Priority"     = "#B2182B",
        "Tier 2 — Moderate Priority" = "#F4A582",
        "Tier 3 — Lower Priority"    = "#92C5DE",
        "Tier 4 — Candidate"         = "#E0E0E0"
      ))
    )
  )
  draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
}

# ---- 7. Integration summary report ----
n_tier1 <- sum(final_targets$tier == "Tier 1 — High Priority")
n_tier2 <- sum(final_targets$tier == "Tier 2 — Moderate Priority")
n_druggable <- sum(final_targets$n_drugs > 0, na.rm = TRUE)
n_immune <- sum(final_targets$n_sig_immune > 0, na.rm = TRUE)

sink("../results/HubGene_Integration/integration_summary.txt")
cat("============================================================\n")
cat(" AI-Assisted Target Discovery — Integration Report\n")
cat(" HCC Prognostic Therapeutic Targets\n")
cat(" Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("============================================================\n\n")

cat(" Evidence dimensions integrated: 9\n")
cat("   [1] DEG          — Differential expression (DESeq2)\n")
cat("   [2] GO/KEGG      — Functional enrichment (clusterProfiler)\n")
cat("   [3] PPI          — Protein-protein interaction network (STRING)\n")
cat("   [4] LASSO        — Penalized Cox feature selection (glmnet)\n")
cat("   [5] RandomForest — Ensemble importance ranking (ranger)\n")
cat("   [6] Survival      — KM, Cox, time-ROC validation\n")
cat("   [7] GEO           — External cohort validation\n")
cat("   [8] Immune         — Tumor microenvironment infiltration\n")
cat("   [9] Drug           — DGIdb druggability & repurposing\n\n")

cat(" Candidate genes evaluated: ", nrow(final_targets), "\n")
cat(" Tier 1 (High Priority):    ", n_tier1, "\n")
cat(" Tier 2 (Moderate Priority):", n_tier2, "\n")
cat(" Genes with known drugs:    ", n_druggable, "\n")
cat(" Genes with immune link:    ", n_immune, "\n\n")

cat("--- Top 15 Therapeutic Targets ---\n")
print(final_targets %>%
        slice_head(n = 15) %>%
        select(gene, tier, composite_100, n_drugs, n_sig_immune, HR))
sink()

# ---- 8. Final summary ----
cat("\n============================================================\n")
cat(" Hub Gene Integration — Final Target Prioritization\n")
cat("============================================================\n")
cat("\n 9-Dimension evidence integration complete.\n")
cat(" Candidates evaluated: ", nrow(final_targets), "\n")
cat(" Tier 1 (High Priority):     ", n_tier1, "\n")
cat(" Tier 2 (Moderate Priority): ", n_tier2, "\n")
cat(" Druggable targets (≥ 1 drug):", n_druggable, "\n")
cat(" Immune-correlated targets:  ", n_immune, "\n")
cat("\n --- Top 10 Final Targets ---\n")
print(final_targets %>%
        slice_head(n = 10) %>%
        select(gene, tier, composite_100, n_drugs, HR))
cat("\n Output saved to results/HubGene_Integration/\n")
cat(" Manuscript table: manuscript/supplementary_table_targets.csv\n")
cat("============================================================\n")
