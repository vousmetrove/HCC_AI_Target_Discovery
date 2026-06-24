###############################################################################
# 07_RandomForest.R
# ---------------------------------------------------------------------------
# Purpose:
#   Apply Random Forest (ranger) to rank DEGs by their prognostic importance
#   and complement the LASSO selection with an ensemble-based, non-linear
#   feature importance assessment. RF captures complex gene-gene interactions
#   and nonlinear expression-survival relationships that LASSO (being a linear
#   model) may miss.
#
# Why Random Forest for target discovery:
#   - Handles p ≫ n without overfitting (ensemble of decision trees)
#   - Built-in variable importance (permutation-based, impurity-based)
#   - Captures non-linear associations and high-order interactions
#   - Robust to outliers and non-normal expression distributions
#   - Provides a complementary gene ranking to LASSO (consensus = robust targets)
#
# Why ranger (not randomForestSRC):
#   - ranger is a fast C++ implementation suitable for high-dimensional data
#   - Native survival forest support via Surv(time, event) ~ . formula
#   - Permutation importance (altmann = TRUE corrects for bias toward high-cardinality features)
#   - Designed for feature selection workflows, not just survival modeling
#
# Expression matrix: DESeq2 rlog (same as LASSO — single normalized matrix)
#   Loaded from results/DEG/rlog_normalized.rds
#
# Input from previous steps:
#   - results/DEG/rlog_normalized.rds   (Step 03 — DESeq2 normalized matrix)
#   - results/DEG/DEG_list.rds          (Step 03 — significant DEGs)
#   - results/LASSO/lasso_genes.rds     (Step 06 — LASSO-selected genes)
#   - results/PPI/hub_genes.rds         (Step 05 — PPI hub genes, optional)
#
# Workflow:
#   1. Load DEG rlog expression + survival data (tumor samples only)
#   2. Hyperparameter tuning: grid search mtry, min.node.size, num.trees
#   3. Train survival Random Forest with ranger
#   4. Extract permutation-based variable importance (Altmann P-values)
#   5. Rank genes by importance; select top prognostic candidates
#   6. Evaluate model via OOB C-index
#   7. Intersection analysis with LASSO + PPI hub genes
#
# Output:
#   results/RandomForest/rf_model.rds          — fitted ranger object
#   results/RandomForest/rf_importance.csv     — full ranked gene importance
#   results/RandomForest/rf_top30.csv          — top 30 genes
#   results/RandomForest/lasso_rf_consensus.rds / .csv
#   figures/rf_vimp.pdf                        — top 30 VIMP bar plot
#
# Dependencies: ranger, survival, dplyr, ggplot2
###############################################################################

# ---- 0. Environment ----
library(ranger)
library(survival)
library(dplyr)
library(tibble)
library(ggplot2)

set.seed(2024)
dir.create("../results/RandomForest", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load and prepare data ----
# Using DESeq2 rlog as the single normalized expression matrix (no FPKM)
rlog_mat  <- readRDS("../results/DEG/rlog_normalized.rds")
clinical  <- readRDS("../data/processed/clinical_clean.rds")
deg_list  <- readRDS("../results/DEG/DEG_list.rds")
metadata  <- readRDS("../data/processed/metadata.rds")

# Tumor-only: survival data only available for tumor patients
tumor_barcodes <- metadata %>% filter(tissue == "Tumor") %>% pull(barcode)
tumor_cols <- intersect(colnames(rlog_mat), tumor_barcodes)
deg_genes  <- intersect(deg_list$gene_id, rownames(rlog_mat))
expr_tumor <- rlog_mat[deg_genes, tumor_cols]

# Transpose → patients × genes
expr_mat <- t(expr_tumor)

# Merge with survival
patient_ids <- substr(rownames(expr_mat), 1, 12)
clin <- clinical %>%
  filter(submitter_id %in% patient_ids) %>%
  select(submitter_id, os_time, os_event)

match_idx <- match(patient_ids, clin$submitter_id)
has_clin  <- !is.na(match_idx)
expr_mat  <- expr_mat[has_clin, ]
clin      <- clin[match_idx[has_clin], ]

message("Analysis set: ", nrow(expr_mat), " patients × ", ncol(expr_mat), " DEGs")

# ---- 2. Clean gene names for ranger formula compatibility ----
# ranger uses formula interface; special characters in gene names cause errors.

original_names <- colnames(expr_mat)
clean_names    <- make.names(original_names, unique = TRUE)
colnames(expr_mat) <- clean_names

name_map <- data.frame(
  original = original_names,
  clean    = clean_names,
  stringsAsFactors = FALSE
)

# ---- 3. Build analysis data frame ----
rf_data <- data.frame(
  os_time  = clin$os_time,
  os_event = clin$os_event,
  expr_mat,
  check.names = FALSE
)

# ---- 4. Hyperparameter tuning via grid search ----
# ranger key parameters:
#   mtry:          number of variables randomly sampled at each split
#                  (default: sqrt(p) for classification, p/3 for survival)
#   min.node.size: minimum number of observations in terminal node
#                  (larger = more regularization, less overfitting)
#   num.trees:     number of trees (more = more stable importance, slower)
#
# We evaluate combinations by OOB C-index (Harrell's concordance).

p <- ncol(expr_mat)
mtry_values  <- unique(round(c(sqrt(p), p / 3, p / 5, p / 10)))
mtry_values  <- mtry_values[mtry_values >= 2 & mtry_values <= p]
node_values  <- c(3, 5, 10, 20)
n_trees      <- 1000   # fixed for importance stability

tune_grid <- expand.grid(
  mtry          = mtry_values,
  min.node.size = node_values,
  cindex        = NA_real_,
  stringsAsFactors = FALSE
)

message("Tuning grid: ", nrow(tune_grid), " combinations")

for (i in seq_len(nrow(tune_grid))) {
  fit_tune <- ranger(
    Surv(os_time, os_event) ~ .,
    data          = rf_data,
    num.trees     = n_trees,
    mtry          = tune_grid$mtry[i],
    min.node.size = tune_grid$min.node.size[i],
    importance    = "none",           # skip importance for speed during tuning
    seed          = 2024,
    verbose       = FALSE
  )
  tune_grid$cindex[i] <- fit_tune$prediction.error
  message("  mtry=", tune_grid$mtry[i],
          " nodesize=", tune_grid$min.node.size[i],
          " → OOB C-index=", round(1 - fit_tune$prediction.error, 4))
}

# Select best parameters (lowest prediction error = highest C-index)
best_row <- which.min(tune_grid$cindex)
best_mtry  <- tune_grid$mtry[best_row]
best_nodes <- tune_grid$min.node.size[best_row]

message("Best: mtry=", best_mtry, " nodesize=", best_nodes,
        " C-index=", round(1 - tune_grid$cindex[best_row], 4))

# ---- 5. Final model with permutation importance ----
# importance = "permutation":   drop in predictive performance when variable is permuted
# importance = "altmann":       permutation importance with Altmann's P-value
#                               (corrects for bias; more robust but slower)
#
# We use "permutation" for ranking; for publication we add Altmann P-values
# on the top 100 genes as a robustness check.

rf_model <- ranger(
  Surv(os_time, os_event) ~ .,
  data           = rf_data,
  num.trees      = n_trees,
  mtry           = best_mtry,
  min.node.size  = best_nodes,
  importance     = "permutation",
  seed           = 2024,
  verbose        = TRUE
)

saveRDS(rf_model, "../results/RandomForest/rf_model.rds")

# OOB C-index: ranger reports 1 - C-index as prediction.error for survival
oob_cindex <- 1 - rf_model$prediction.error
message("Final OOB C-index: ", round(oob_cindex, 4))

# ---- 6. Variable Importance extraction ----
# Positive VIMP → gene contributes to survival prediction
# Negative VIMP → permuting the gene improves prediction (noise)

vimp_raw <- rf_model$variable.importance
vimp_sorted <- sort(vimp_raw, decreasing = TRUE)

vimp_df <- data.frame(
  gene_clean = names(vimp_sorted),
  importance = as.numeric(vimp_sorted)
) %>%
  left_join(name_map, by = c("gene_clean" = "clean")) %>%
  mutate(gene = ifelse(is.na(original), gene_clean, original)) %>%
  select(gene, importance) %>%
  arrange(desc(importance))

# Filter positive importance
vimp_df_pos <- vimp_df %>% filter(importance > 0)

write.csv(vimp_df, "../results/RandomForest/rf_importance.csv", row.names = FALSE)
saveRDS(vimp_df, "../results/RandomForest/rf_importance.rds")

message("Genes with positive importance: ", nrow(vimp_df_pos),
        " / ", nrow(vimp_df))

# ---- 7. Altmann P-values for top genes (robustness) ----
# Altmann's method uses permutations of the outcome to compute a P-value
# for each variable's importance, correcting for multiple testing implicitly.
# Run on the top 100 genes (full set would be computationally expensive).

top100_genes <- vimp_df_pos %>% slice_head(n = 100) %>% pull(gene_clean)

if (length(top100_genes) >= 10) {
  # Subset data to top 100 genes + survival columns
  rf_data_top100 <- rf_data[, c("os_time", "os_event", top100_genes)]

  rf_altmann <- ranger(
    Surv(os_time, os_event) ~ .,
    data           = rf_data_top100,
    num.trees      = 500,
    mtry           = max(2, floor(sqrt(length(top100_genes)))),
    min.node.size  = best_nodes,
    importance     = "altmann",
    seed           = 2024,
    verbose        = FALSE
  )

  alt_pvals <- rf_altmann$variable.importance.pval
  alt_imp   <- rf_altmann$variable.importance

  altmann_df <- data.frame(
    gene_clean = names(alt_imp),
    importance = as.numeric(alt_imp),
    pval       = as.numeric(alt_pvals),
    stringsAsFactors = FALSE
  ) %>%
    left_join(name_map, by = c("gene_clean" = "clean")) %>%
    mutate(gene = ifelse(is.na(original), gene_clean, original)) %>%
    select(gene, importance, pval) %>%
    arrange(pval)

  write.csv(altmann_df, "../results/RandomForest/rf_altmann_pvals.csv",
            row.names = FALSE)
  saveRDS(altmann_df, "../results/RandomForest/rf_altmann_pvals.rds")

  message("Altmann P < 0.05 genes: ", sum(altmann_df$pval < 0.05, na.rm = TRUE))
}

# ---- 8. Top gene selection ----
top_n  <- min(30, nrow(vimp_df_pos))
top_rf <- vimp_df_pos %>% slice_head(n = top_n)
write.csv(top_rf, "../results/RandomForest/rf_top30.csv", row.names = FALSE)

# ---- 9. Visualization ----
# Top 30 variable importance bar plot

vimp_plot <- ggplot(top_rf, aes(x = reorder(gene, importance), y = importance)) +
  geom_col(aes(fill = importance), width = 0.65) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C") +
  coord_flip() +
  labs(
    title    = "Random Forest (ranger) — Top 30 Prognostic Genes",
    subtitle = paste0("OOB C-index = ", round(oob_cindex, 4),
                      " | ntree = ", n_trees,
                      " | mtry = ", best_mtry),
    x        = "Gene",
    y        = "Variable Importance (permutation VIMP)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("../figures/rf_vimp.pdf", vimp_plot, width = 8, height = 7)

# ---- 10. Intersection: LASSO ∩ RF ∩ PPI ----
# Three-way consensus: genes selected by multiple orthogonal methods are the
# highest-confidence therapeutic targets.

lasso_genes <- readRDS("../results/LASSO/lasso_genes.rds")
lasso_set   <- lasso_genes$gene_symbol

# PPI hub genes (if available)
ppi_hub_set <- tryCatch({
  ppi_hub <- readRDS("../results/PPI/ppi_hub_genes.rds")
  ppi_hub$gene
}, error = function(e) NULL)

# Top RF genes (comparable number to LASSO)
rf_top_set <- vimp_df_pos %>%
  slice_head(n = max(30, nrow(lasso_genes))) %>%
  pull(gene)

# LASSO ∩ RF
consensus_lr <- intersect(lasso_set, rf_top_set)

# LASSO ∩ RF ∩ PPI (three-way)
if (!is.null(ppi_hub_set)) {
  consensus_three <- intersect(consensus_lr, ppi_hub_set)
  message("LASSO ∩ RF ∩ PPI: ", length(consensus_three), " genes")
} else {
  consensus_three <- character(0)
}

# Build consensus table
consensus_df <- data.frame(
  gene          = consensus_lr,
  lasso_coef    = lasso_genes$coefficient[match(consensus_lr, lasso_genes$gene_symbol)],
  rf_importance = vimp_df$importance[match(consensus_lr, vimp_df$gene)],
  in_ppi        = consensus_lr %in% (ppi_hub_set %||% character(0)),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(rf_importance))

write.csv(consensus_df, "../results/RandomForest/lasso_rf_consensus.csv", row.names = FALSE)
saveRDS(consensus_df, "../results/RandomForest/lasso_rf_consensus.rds")

# Three-way consensus
if (length(consensus_three) > 0) {
  three_way_df <- consensus_df %>% filter(in_ppi)
  write.csv(three_way_df, "../results/RandomForest/three_way_consensus.csv", row.names = FALSE)
  saveRDS(three_way_df, "../results/RandomForest/three_way_consensus.rds")
}

message("LASSO ∩ RF consensus: ", nrow(consensus_df), " genes")
print(consensus_df %>% select(gene, lasso_coef, rf_importance, in_ppi))

# ---- 11. Summary ----
cat("\n============================================================\n")
cat(" Random Forest (ranger) feature selection complete.\n")
cat(" Total DEGs evaluated:  ", ncol(expr_mat), "\n")
cat(" OOB C-index:           ", round(oob_cindex, 4), "\n")
cat(" Best mtry / nodesize:  ", best_mtry, "/", best_nodes, "\n")
cat(" Positive VIMP genes:   ", nrow(vimp_df_pos), "\n")
cat(" LASSO ∩ RF consensus:  ", nrow(consensus_df), "genes\n")
if (!is.null(ppi_hub_set)) {
  cat(" LASSO ∩ RF ∩ PPI:     ", length(consensus_three), "genes\n")
}
cat(" Output saved to results/RandomForest/\n")
cat("============================================================\n")
