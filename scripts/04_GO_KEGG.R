###############################################################################
# 04_GO_KEGG.R
# ---------------------------------------------------------------------------
# Purpose:
#   Perform functional enrichment analysis on the differentially expressed
#   genes (DEGs) identified in Step 03. This reveals the biological processes,
#   molecular functions, cellular components (GO terms), and signaling pathways
#   (KEGG) that are dysregulated in HCC, providing mechanistic context for
#   the candidate therapeutic targets.
#
# Enrichment methods:
#   1. Over-Representation Analysis (ORA) — clusterProfiler::enrichGO / enrichKEGG
#      Tests whether DEGs are enriched in specific GO/KEGG categories relative
#      to a background gene set (all expressed genes).
#   2. Gene Set Enrichment Analysis (GSEA) — clusterProfiler::gseGO / gseKEGG
#      Rank-based approach using the full gene list ordered by log2FC; detects
#      coordinated up- or down-regulation of pathway members.
#
# Visualizations:
#   - Dot plot of top enriched terms
#   - Enrichment network (cnetplot)
#   - Heatmap of gene-to-term membership
#   - Ridge plot for GSEA
#
# Output:
#   results/GO_KEGG/GO_enrichment.rds        — ORA results for GO
#   results/GO_KEGG/KEGG_enrichment.rds      — ORA results for KEGG
#   results/GO_KEGG/GSEA_GO.rds              — GSEA results for GO
#   results/GO_KEGG/GSEA_KEGG.rds            — GSEA results for KEGG
#   figures/dotplot_GO.pdf / dotplot_KEGG.pdf
#   figures/cnetplot_GO.pdf / cnetplot_KEGG.pdf
#
# Dependencies: clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, ggplot2
###############################################################################

# ---- 0. Environment ----
library(clusterProfiler)
library(org.Hs.eg.db)      # Human gene annotation database
library(DOSE)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

set.seed(2024)
dir.create("../results/GO_KEGG", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load DEG results ----
deg_list <- readRDS("../results/DEG/DEG_list.rds")     # significant DEGs
res_df   <- readRDS("../results/DEG/DESeq2_results.rds")  # full results w/ log2FC

# ---- 2. Gene ID conversion ----
# Convert Ensembl gene IDs to Entrez IDs required by clusterProfiler.
# Many TCGA expression matrices use Ensembl IDs; clusterProfiler's enrich
# functions need Entrez IDs for database lookup.

# Extract the gene ID column (may need to strip version numbers: ENSGxxx.1 → ENSGxxx)
gene_ids <- gsub("\\..*", "", res_df$gene_id)

# Map to Entrez and SYMBOL using org.Hs.eg.db
gene_map <- bitr(
  gene_ids,
  fromType = "ENSEMBL",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
)

# Merge with DEG results
res_annotated <- res_df %>%
  mutate(ENSEMBL = gsub("\\..*", "", gene_id)) %>%
  inner_join(gene_map, by = "ENSEMBL") %>%
  distinct(ENTREZID, .keep_all = TRUE)   # remove duplicate gene mappings

# ---- 3. Background gene set ----
# The background should be all expressed genes (those retained after filtering
# in Step 02). Using all genes in the DESeq2 results satisfies this.

background_entrez <- unique(
  bitr(gsub("\\..*", "", res_df$gene_id),
       fromType = "ENSEMBL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db
  )$ENTREZID
)

message("Background genes (Entrez): ", length(background_entrez))

# ---- 4. DEG lists by direction ----
up_genes   <- res_annotated %>% filter(significance == "Upregulated")   %>% pull(ENTREZID)
down_genes <- res_annotated %>% filter(significance == "Downregulated") %>% pull(ENTREZID)
all_deg    <- c(up_genes, down_genes)

message("DEGs mapped to Entrez — Up: ", length(up_genes),
        " | Down: ", length(down_genes))

# ---- 5. GO Over-Representation Analysis ----
# Tests all three GO sub-ontologies: BP (Biological Process),
# CC (Cellular Component), MF (Molecular Function).

go_enrich <- enrichGO(
  gene          = all_deg,
  universe      = background_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",            # BP + CC + MF
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE,            # convert Entrez → SYMBOL in output
  minGSSize     = 10,              # minimum gene set size
  maxGSSize     = 500              # maximum gene set size
)

# Simplify: remove redundant GO terms using semantic similarity
go_simplified <- simplify(go_enrich, cutoff = 0.7, by = "p.adjust")

message("GO terms enriched: ", nrow(go_enrich), " total, ",
        nrow(go_simplified), " after simplification")

saveRDS(go_enrich,     "../results/GO_KEGG/GO_enrichment.rds")
saveRDS(go_simplified, "../results/GO_KEGG/GO_enrichment_simplified.rds")

# ---- 6. KEGG Pathway Over-Representation Analysis ----

kegg_enrich <- enrichKEGG(
  gene          = all_deg,
  universe      = background_entrez,
  organism      = "hsa",           # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 10,
  maxGSSize     = 500
)

message("KEGG pathways enriched: ", nrow(kegg_enrich))
saveRDS(kegg_enrich, "../results/GO_KEGG/KEGG_enrichment.rds")

# ---- 7. GSEA (Gene Set Enrichment Analysis) ----
# Uses the full gene list ranked by log2FC, providing more sensitivity than ORA
# because it considers fold-change magnitude and direction.

# Prepare ranked gene list
ranked_genes <- res_annotated %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange)
names(ranked_genes) <- res_annotated %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(ENTREZID)

# GSEA — GO
gsea_go <- gseGO(
  geneList      = ranked_genes,
  ont           = "ALL",
  OrgDb         = org.Hs.eg.db,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  eps           = 0,              # no p-value boundary
  seed          = 2024
)
saveRDS(gsea_go, "../results/GO_KEGG/GSEA_GO.rds")

# GSEA — KEGG
gsea_kegg <- gseKEGG(
  geneList      = ranked_genes,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  eps           = 0,
  seed          = 2024
)
saveRDS(gsea_kegg, "../results/GO_KEGG/GSEA_KEGG.rds")

# ---- 8. Visualizations ----
# 8a. Dot plot — top 15 GO terms per sub-ontology, top 20 KEGG pathways

if (nrow(go_simplified) > 0) {
  # GO dot plot — split by ontology
  go_simplified@result <- go_simplified@result %>%
    mutate(ONTOLOGY = recode(ONTOLOGY, BP = "Biological Process",
                            CC = "Cellular Component", MF = "Molecular Function"))

  dot_go <- dotplot(go_simplified, showCategory = 15, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    labs(title = "GO Enrichment Analysis", x = "Gene Ratio") +
    theme_minimal(base_size = 11)
  ggsave("../figures/dotplot_GO.pdf", dot_go, width = 9, height = 10)
}

if (nrow(kegg_enrich) > 0) {
  dot_kegg <- dotplot(kegg_enrich, showCategory = 20) +
    labs(title = "KEGG Pathway Enrichment", x = "Gene Ratio") +
    theme_minimal(base_size = 11)
  ggsave("../figures/dotplot_KEGG.pdf", dot_kegg, width = 9, height = 7)
}

# 8b. Cnet plot — gene-to-term network for top enriched terms

if (nrow(go_simplified) > 0) {
  cnet_go <- cnetplot(
    go_simplified,
    showCategory   = 5,
    foldChange     = ranked_genes,
    circular       = TRUE,
    colorEdge      = TRUE,
    cex_label_gene = 0.6,
    cex_label_category = 0.8
  ) +
    labs(title = "GO Gene-Concept Network")
  ggsave("../figures/cnetplot_GO.pdf", cnet_go, width = 10, height = 9)
}

# 8c. GSEA ridge plot
if (nrow(as.data.frame(gsea_kegg)) > 0) {
  ridge_kegg <- ridgeplot(gsea_kegg, showCategory = 15) +
    labs(title = "GSEA — KEGG Pathways",
         x = "log2 Fold Change (Tumor vs. Normal)") +
    theme_minimal()
  ggsave("../figures/ridgeplot_GSEA_KEGG.pdf", ridge_kegg, width = 10, height = 8)
}

# ---- 9. Write enrichment tables ----
# Save as CSV for manuscript tables / supplementary
if (nrow(go_simplified) > 0) {
  go_df <- as.data.frame(go_simplified) %>%
    arrange(p.adjust)
  write.csv(go_df, "../results/GO_KEGG/GO_enrichment_table.csv", row.names = FALSE)
}

if (nrow(kegg_enrich) > 0) {
  kegg_df <- as.data.frame(kegg_enrich) %>%
    arrange(p.adjust)
  write.csv(kegg_df, "../results/GO_KEGG/KEGG_enrichment_table.csv", row.names = FALSE)
}

# ---- 10. Summary ----
cat("\n============================================================\n")
cat(" GO & KEGG enrichment analysis complete.\n")
cat(" GO  terms enriched:", nrow(go_simplified), "\n")
cat(" KEGG pathways:     ", nrow(kegg_enrich), "\n")
cat(" Output saved to results/GO_KEGG/\n")
cat("============================================================\n")
