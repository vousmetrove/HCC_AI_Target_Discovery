###############################################################################
# 09_GEO_validation.R
# ---------------------------------------------------------------------------
# Purpose:
#   Validate the TCGA-derived prognostic gene signature in independent external
#   cohorts from the Gene Expression Omnibus (GEO). External validation is
#   essential to confirm that the identified targets are generalizable and not
#   artifacts of the TCGA dataset. This step strengthens the translational
#   relevance of the findings.
#
# Validation cohorts (priority order for HCC):
#   1. GSE14520 — 221 HCC tumors + paired non-tumor (Roessler et al.)
#   2. GSE76427 — 115 HCC tumors
#   3. ICGC-LIRI-JP — 240 HCC tumors (via ICGC portal or GEO mirror)
#   4. GSE10143  — 80 HCC tumors
#
# Validation approach:
#   1. Download & preprocess GEO expression matrices (microarray or RNA-seq)
#   2. Map probe IDs to gene symbols; collapse multi-probe genes by max mean
#   3. Compute risk score using TCGA-derived coefficients (frozen model)
#   4. Stratify patients into high- vs. low-risk (median split)
#   5. Kaplan-Meier survival analysis
#   6. Time-dependent ROC (validation AUC)
#   7. Meta-analysis across cohorts (forest plot of HRs)
#
# Key concept: "Frozen model" — the coefficients and gene set are fixed from
# the TCGA training; no refitting on external data. This tests generalizability.
#
# Output:
#   results/GEO_validation/GEO_metadata.rds          — cohort-level summary
#   results/GEO_validation/GSE14520_validation.rds   — per-cohort results
#   results/GEO_validation/validation_km.pdf         — combined KM curves
#   results/GEO_validation/validation_meta_forest.pdf — meta-analysis forest
#   results/GEO_validation/validation_roc.pdf        — external AUC comparison
#
# Dependencies: GEOquery, limma, survival, survminer, timeROC, metafor, dplyr
###############################################################################

# ---- 0. Environment ----
library(GEOquery)
library(limma)            # for microarray normalization & probe collapsing
library(survival)
library(survminer)
library(timeROC)
library(metafor)          # meta-analysis (random-effects)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)

set.seed(2024)
dir.create("../results/GEO_validation", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load TCGA signature (frozen model) ----
consensus_df <- readRDS("../results/RandomForest/lasso_rf_consensus.rds")
lasso_genes  <- readRDS("../results/LASSO/lasso_genes.rds")

# Define the frozen signature: genes and their LASSO coefficients
sig_genes <- if (nrow(consensus_df) > 0) consensus_df$gene else lasso_genes$gene_symbol
sig_coefs <- lasso_genes %>%
  filter(gene_symbol %in% sig_genes) %>%
  select(gene_symbol, coefficient)

# If any consensus genes don't have LASSO coefs, assign sign by log2FC
deg_all <- readRDS("../results/DEG/DESeq2_results.rds")
for (g in sig_genes) {
  if (!(g %in% sig_coefs$gene_symbol)) {
    deg_row <- deg_all %>% filter(gene_id == g)
    sign_val <- if (nrow(deg_row) > 0) sign(deg_row$log2FoldChange[1]) else 0
    sig_coefs <- rbind(sig_coefs, data.frame(gene_symbol = g, coefficient = sign_val))
  }
}

message("Frozen signature: ", nrow(sig_coefs), " genes")

# ---- 2. GEO dataset definitions ----
# List of validation datasets with accession IDs and metadata

geo_datasets <- list(
  GSE14520 = list(
    accession = "GSE14520",
    platform  = "GPL3921",    # Affymetrix HT Human Genome U133A Array
    n_samples = 221,          # tumor samples
    reference = "Roessler et al., 2010"
  ),
  GSE76427 = list(
    accession = "GSE76427",
    platform  = "GPL10558",   # Illumina HumanHT-12 V4.0
    n_samples = 115,
    reference = "Grinchuk et al., 2018"
  )
)

# ---- 3. Validation function ----
# Modular function to process one GEO dataset end-to-end

validate_cohort <- function(dataset_info, sig_coefs) {

  accession <- dataset_info$accession
  message("\n>>> Processing ", accession, " ...")

  # 3a. Download GEO Series
  gse <- tryCatch(
    getGEO(accession, GSEMatrix = TRUE, getGPL = TRUE),
    error = function(e) { message("  Failed to download ", accession); return(NULL) }
  )
  if (is.null(gse)) return(NULL)

  # Extract expression matrix
  if (length(gse) > 1) {
    # Multiple platforms — pick the mRNA/miRNA one with most features
    n_features <- sapply(gse, function(x) nrow(exprs(x)))
    expr_set    <- gse[[which.max(n_features)]]
  } else {
    expr_set <- gse[[1]]
  }

  expr_mat <- exprs(expr_set)

  # 3b. Extract phenotype (survival) data
  pheno <- pData(expr_set)

  # Identify survival columns (variable names vary across GEO series)
  os_time_col <- grep("overall survival|os\\.time|survival time|days",
                       colnames(pheno), ignore.case = TRUE, value = TRUE)
  os_event_col <- grep("status|event|death|outcome|vital",
                        colnames(pheno), ignore.case = TRUE, value = TRUE)

  # Filter out non-tumor samples if annotation exists
  tissue_col <- grep("tissue|source|histolog", colnames(pheno),
                     ignore.case = TRUE, value = TRUE)

  message("  Samples: ", ncol(expr_mat))

  # 3c. Map probes to gene symbols
  # Use feature data from the ExpressionSet, or use annotation package
  feat_data <- fData(expr_set)

  # Try common gene symbol columns
  symbol_col <- grep("gene.*symbol|symbol|gene name",
                      colnames(feat_data), ignore.case = TRUE, value = TRUE)
  if (length(symbol_col) == 0) {
    # Fallback: use the row names (may already be symbols for some platforms)
    message("  No gene symbol column found — using rownames")
    gene_symbols <- rownames(expr_mat)
  } else {
    gene_symbols <- feat_data[[symbol_col[1]]]
  }

  # Handle multi-probe genes: collapse by taking the probe with max mean expression
  gene_list <- split(seq_len(nrow(expr_mat)), gene_symbols)
  gene_list <- gene_list[names(gene_list) != "" & !is.na(names(gene_list))]

  # Collapse: for each gene, take probe(s) with maximum mean expression
  collapsed_expr <- t(sapply(gene_list, function(idx) {
    if (length(idx) == 1) {
      expr_mat[idx, ]
    } else {
      # Use limma::avereps or take max-mean probe
      sub_mat    <- expr_mat[idx, , drop = FALSE]
      row_means  <- rowMeans(sub_mat, na.rm = TRUE)
      sub_mat[which.max(row_means), ]
    }
  }))

  message("  Genes after collapse: ", nrow(collapsed_expr))

  # 3d. Match signature genes to this dataset's gene space
  common_genes <- intersect(sig_coefs$gene_symbol, rownames(collapsed_expr))
  missing_genes <- setdiff(sig_coefs$gene_symbol, rownames(collapsed_expr))

  message("  Signature genes found: ", length(common_genes),
          " / ", nrow(sig_coefs), " (missing: ", length(missing_genes), ")")

  if (length(common_genes) < 2) {
    message("  Insufficient gene overlap — skipping cohort")
    return(NULL)
  }

  # Subset and compute risk score
  expr_sig_sub <- collapsed_expr[common_genes, , drop = FALSE]
  coefs_sub    <- sig_coefs$coefficient[match(common_genes, sig_coefs$gene_symbol)]

  # Compute risk score (note: may need to handle scale differences between
  # RNA-seq and microarray — we z-score within cohort before scoring)
  expr_z <- t(scale(t(expr_sig_sub)))
  risk_score <- as.vector(t(expr_z) %*% coefs_sub)

  # 3e. Extract survival — manual column parsing is dataset-specific
  # For robustness, try multiple column naming conventions

  # Survival time: try known column names
  os_time <- NA
  time_candidates <- c(
    "os.time", "overall_survival", "survival_time",
    "days_to_death", "follow_up_time", "overall survival (months):ch1",
    "survival time (months):ch1", "overall survival:ch1"
  )
  for (col in time_candidates) {
    if (col %in% colnames(pheno)) {
      os_time <- as.numeric(as.character(pheno[[col]]))
      break
    }
  }
  # Fallback: use the first column matching "survival" or "time"
  if (all(is.na(os_time))) {
    time_hits <- grep("survival|os\\.|follow|month|day",
                      colnames(pheno), ignore.case = TRUE, value = TRUE)
    if (length(time_hits) > 0) {
      os_time <- as.numeric(as.character(pheno[[time_hits[1]]]))
    }
  }

  # Convert months to days if needed
  if (!all(is.na(os_time)) && max(os_time, na.rm = TRUE) < 200) {
    os_time <- os_time * 30.44   # months → days
    message("  Survival times appear to be in months — converted to days")
  }

  # Event indicator
  os_event <- NA
  event_candidates <- c(
    "os.event", "vital_status", "death", "event",
    "status", "outcome", "overall survival event:ch1",
    "death due to tumor:ch1", "recurrence:ch1"
  )
  for (col in event_candidates) {
    if (col %in% colnames(pheno)) {
      val <- as.character(pheno[[col]])
      os_event <- ifelse(grepl("dead|event|yes|recurrence|1|positive",
                                val, ignore.case = TRUE), 1, 0)
      break
    }
  }
  # Fallback
  if (all(is.na(os_event))) {
    event_hits <- grep("event|status|death|alive|dead|vital",
                       colnames(pheno), ignore.case = TRUE, value = TRUE)
    if (length(event_hits) > 0) {
      val <- as.character(pheno[[event_hits[1]]])
      os_event <- ifelse(grepl("dead|1|event|yes|positive",
                                val, ignore.case = TRUE), 1, 0)
    }
  }

  # 3f. Filter valid survival records
  valid <- complete.cases(os_time, os_event) & os_time > 0 & os_event %in% c(0, 1)
  if (sum(valid) < 20) {
    message("  Insufficient valid survival records — skipping")
    return(NULL)
  }

  surv_df <- data.frame(
    sample     = colnames(expr_sig_sub)[valid],
    os_time    = os_time[valid],
    os_event   = os_event[valid],
    risk_score = risk_score[valid],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      risk_group = if_else(risk_score > median(risk_score),
                          "High", "Low"),
      risk_group = factor(risk_group, levels = c("Low", "High")),
      cohort     = accession
    )

  # 3g. Kaplan-Meier
  km_fit <- survfit(Surv(os_time, os_event) ~ risk_group, data = surv_df)
  logrank_p <- surv_pvalue(km_fit)$pval

  # 3h. Cox regression
  cox_fit <- coxph(Surv(os_time, os_event) ~ risk_score, data = surv_df)
  cox_summary <- summary(cox_fit)

  # 3i. Time-dependent ROC (1-year, 3-year)
  time_pts <- c(365, 1095)
  names(time_pts) <- c("1-year", "3-year")
  roc_val <- tryCatch({
    timeROC(
      T     = surv_df$os_time,
      delta = surv_df$os_event,
      marker = surv_df$risk_score,
      cause  = 1,
      times  = time_pts,
      ROC    = TRUE
    )
  }, error = function(e) NULL)

  # 3j. Compile results
  result <- list(
    accession   = accession,
    n_patients  = nrow(surv_df),
    n_genes     = length(common_genes),
    missing_pct = length(missing_genes) / nrow(sig_coefs),
    logrank_p   = logrank_p,
    cox_HR      = cox_summary$conf.int[1, "exp(coef)"],
    cox_HR_lower = cox_summary$conf.int[1, "lower .95"],
    cox_HR_upper = cox_summary$conf.int[1, "upper .95"],
    cox_p       = cox_summary$coefficients[1, "Pr(>|z|)"],
    auc_1yr     = if (!is.null(roc_val)) roc_val$AUC[1] else NA,
    auc_3yr     = if (!is.null(roc_val)) roc_val$AUC[2] else NA,
    surv_df     = surv_df,
    km_fit      = km_fit
  )

  message("  N = ", nrow(surv_df),
          " | Log-rank P = ", format.pval(logrank_p, digits = 2),
          " | HR = ", round(result$cox_HR, 2))
  return(result)
}

# ---- 4. Run validation across cohorts ----
validation_results <- lapply(geo_datasets, validate_cohort, sig_coefs = sig_coefs)
validation_results <- validation_results[!sapply(validation_results, is.null)]

message("\nSuccessfully validated in ", length(validation_results),
        " GEO cohorts")

saveRDS(validation_results, "../results/GEO_validation/validation_results.rds")

# ---- 5. Meta-analysis forest plot ----
if (length(validation_results) >= 2) {
  # Combine HRs across cohorts with random-effects meta-analysis
  meta_data <- data.frame(
    study    = sapply(validation_results, `[[`, "accession"),
    n        = sapply(validation_results, `[[`, "n_patients"),
    HR       = sapply(validation_results, `[[`, "cox_HR"),
    CI_lower = sapply(validation_results, `[[`, "cox_HR_lower"),
    CI_upper = sapply(validation_results, `[[`, "cox_HR_upper"),
    stringsAsFactors = FALSE
  )

  # Add TCGA as the discovery cohort
  surv_tcga <- readRDS("../results/Survival/survival_data.rds")
  cox_tcga <- coxph(Surv(os_time, os_event) ~ risk_score, data = surv_tcga)
  cox_tcga_s <- summary(cox_tcga)

  meta_data <- rbind(
    data.frame(
      study    = "TCGA-LIHC (discovery)",
      n        = nrow(surv_tcga),
      HR       = cox_tcga_s$conf.int[1, "exp(coef)"],
      CI_lower = cox_tcga_s$conf.int[1, "lower .95"],
      CI_upper = cox_tcga_s$conf.int[1, "upper .95"],
      stringsAsFactors = FALSE
    ),
    meta_data
  )

  # Random-effects meta-analysis
  meta_res <- rma(
    yi   = log(meta_data$HR),
    sei  = (log(meta_data$CI_upper) - log(meta_data$CI_lower)) / (2 * 1.96),
    slab = meta_data$study,
    method = "REML"
  )

  saveRDS(meta_res, "../results/GEO_validation/meta_analysis.rds")

  # Forest plot
  pdf("../figures/validation_meta_forest.pdf", width = 10, height = 5)
  forest(
    meta_res,
    xlab         = "Hazard Ratio (log scale)",
    header       = "Cohort",
    mlab         = "Random-effects summary",
    refline      = 1,
    atransf      = exp,
    at           = log(c(0.5, 1, 2, 5, 10)),
    main         = "Meta-Analysis — Prognostic Gene Signature Across Cohorts"
  )
  dev.off()

  message("Meta-analysis HR: ", round(exp(meta_res$b), 2),
          " [95% CI: ", round(exp(meta_res$ci.lb), 2),
          "–", round(exp(meta_res$ci.ub), 2),
          "] P = ", format.pval(meta_res$pval, digits = 3))
}

# ---- 6. Combined validation KM plot ----
if (length(validation_results) >= 1) {

  # Combine survival data across validation cohorts
  all_surv <- bind_rows(lapply(validation_results, `[[`, "surv_df"))

  km_combined <- survfit(Surv(os_time, os_event) ~ risk_group, data = all_surv)

  km_val_plot <- ggsurvplot(
    km_combined,
    data             = all_surv,
    pval             = TRUE,
    conf.int         = TRUE,
    palette          = c("#377EB8", "#E41A1C"),
    xlab             = "Overall Survival (days)",
    title            = "External Validation — Combined GEO Cohorts",
    legend.title     = "Risk Group",
    legend.labs      = c("Low Risk", "High Risk"),
    ggtheme          = theme_minimal(base_size = 12)
  )

  pdf("../figures/validation_km_combined.pdf", width = 8, height = 7)
  print(km_val_plot, newpage = FALSE)
  dev.off()
}

# ---- 7. Validation AUC summary ----
auc_summary <- data.frame(
  Cohort     = sapply(validation_results, `[[`, "accession"),
  N          = sapply(validation_results, `[[`, "n_patients"),
  AUC_1year  = sapply(validation_results, `[[`, "auc_1yr"),
  AUC_3year  = sapply(validation_results, `[[`, "auc_3yr"),
  HR         = sapply(validation_results, `[[`, "cox_HR"),
  P_logrank  = sapply(validation_results, `[[`, "logrank_p"),
  stringsAsFactors = FALSE
)
write.csv(auc_summary, "../results/GEO_validation/validation_summary.csv",
          row.names = FALSE)
saveRDS(auc_summary, "../results/GEO_validation/validation_summary.rds")

print(auc_summary)

# ---- 8. Summary ----
cat("\n============================================================\n")
cat(" GEO external validation complete.\n")
cat(" Cohorts validated:", length(validation_results), "\n")
for (v in validation_results) {
  cat("   ", v$accession, ": N=", v$n_patients,
      "  HR=", round(v$cox_HR, 2),
      "  P=", format.pval(v$logrank_p, digits = 2), "\n")
}
cat(" Output saved to results/GEO_validation/\n")
cat("============================================================\n")
