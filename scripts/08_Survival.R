###############################################################################
# 08_Survival.R
# ---------------------------------------------------------------------------
# Purpose:
#   Validate the prognostic power of the consensus gene signature
#   (LASSO Step 06 ∩ RF Step 07 ∩ PPI Step 05) through comprehensive
#   survival analysis. Establishes the clinical relevance of identified
#   targets using multiple complementary approaches.
#
# Analyses performed:
#   1. Risk score construction for the final gene signature
#   2. Kaplan-Meier survival curves (high- vs. low-risk groups)
#   3. Log-rank test for survival difference significance
#   4. Univariate Cox proportional hazards regression per gene & risk score
#   5. Multivariate Cox regression adjusting for clinical covariates
#      (stage, grade, age, sex) — tests independence of the signature
#   6. Time-dependent ROC curves (1-, 3-, 5-year) via timeROC
#   7. Nomogram integrating gene signature with clinical variables
#   8. Calibration plots to assess nomogram accuracy
#
# Key clinical questions answered:
#   - Does the gene signature stratify patients into distinct prognostic groups?
#   - Is the signature an independent predictor beyond standard clinical factors?
#   - How well does the model predict survival at specific time horizons?
#   - Can we build a clinically actionable nomogram?
#
# Output:
#   results/Survival/km_curves.pdf / km_curves.rds
#   results/Survival/cox_univariate.rds / cox_multivariate.rds
#   results/Survival/time_roc.rds / time_roc.pdf
#   results/Survival/nomogram.rds / nomogram.pdf
#   results/Survival/calibration.pdf
#   results/Survival/forest_plot.pdf
#
# Dependencies: survival, survminer, timeROC, rms, ggplot2, dplyr
###############################################################################

# ---- 0. Environment ----
library(survival)
library(survminer)
library(timeROC)
library(rms)           # nomogram & calibration (loads after survival to avoid masking)
library(ggplot2)
library(dplyr)
library(tibble)

set.seed(2024)
dir.create("../results/Survival", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load data and gene signature ----
# Using DESeq2 rlog as the single normalized expression matrix (no FPKM)
rlog_mat     <- readRDS("../results/DEG/rlog_normalized.rds")
clinical     <- readRDS("../data/processed/clinical_clean.rds")
metadata     <- readRDS("../data/processed/metadata.rds")
consensus_df <- readRDS("../results/RandomForest/lasso_rf_consensus.rds")
lasso_genes  <- readRDS("../results/LASSO/lasso_genes.rds")

# ---- 2. Prepare tumor expression with survival ----
# Using DESeq2 rlog matrix (single normalized expression modality)
tumor_barcodes <- metadata %>% filter(tissue == "Tumor") %>% pull(barcode)
tumor_cols     <- intersect(colnames(rlog_mat), tumor_barcodes)

# Expression of consensus genes (or LASSO genes if consensus is empty)
sig_genes <- if (nrow(consensus_df) > 0) consensus_df$gene else lasso_genes$gene_symbol
sig_genes <- intersect(sig_genes, rownames(rlog_mat))

message("Signature genes for survival analysis: ", length(sig_genes))
print(sig_genes)

# Subset rlog expression
expr_sig  <- t(rlog_mat[sig_genes, tumor_cols, drop = FALSE])
patient_ids <- substr(rownames(expr_sig), 1, 12)

# Merge with clinical
clin <- clinical %>%
  filter(submitter_id %in% patient_ids) %>%
  select(
    submitter_id, os_time, os_event,
    age_years, gender,
    stage_simple, ajcc_pathologic_stage
  )

match_idx <- match(patient_ids, clin$submitter_id)
has_clin  <- !is.na(match_idx)
expr_sig  <- expr_sig[has_clin, , drop = FALSE]
clin      <- clin[match_idx[has_clin], ]

# ---- 3. Compute risk score ----
# Risk score = Σ (expression_i × LASSO coefficient_i)
# Using LASSO coefficients as weights; if a gene was only in RF, we use
# the sign of the correlation with survival.

coef_lookup <- lasso_genes %>%
  filter(gene_symbol %in% sig_genes) %>%
  select(gene_symbol, coefficient)

# For genes in consensus but not LASSO, assign weight by direction
all_weights <- data.frame(
  gene   = sig_genes,
  weight = rep(NA_real_, length(sig_genes)),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(all_weights))) {
  g <- all_weights$gene[i]
  if (g %in% coef_lookup$gene_symbol) {
    all_weights$weight[i] <- coef_lookup$coefficient[coef_lookup$gene_symbol == g]
  } else {
    # Estimate direction: correlate with DEG log2FC (positive → risk, negative → protective)
    deg <- readRDS("../results/DEG/DESeq2_results.rds")
    deg_gene <- deg %>% filter(gene_id == g)
    if (nrow(deg_gene) > 0) {
      all_weights$weight[i] <- sign(deg_gene$log2FoldChange)
    } else {
      all_weights$weight[i] <- 0
    }
  }
}

# Calculate risk scores
risk_score <- as.vector(expr_sig %*% all_weights$weight)

# ---- 4. Survival data assembly ----
surv_df <- data.frame(
  patient_id = rownames(expr_sig),
  os_time    = clin$os_time,
  os_event   = clin$os_event,
  risk_score = risk_score,
  age        = clin$age_years,
  gender     = clin$gender,
  stage      = clin$stage_simple,
  stringsAsFactors = FALSE
) %>%
  # Stratify into high- vs. low-risk by median
  mutate(
    risk_group = if_else(risk_score > median(risk_score, na.rm = TRUE),
                         "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low", "High"))
  )

saveRDS(surv_df, "../results/Survival/survival_data.rds")

# ---- 5. Kaplan-Meier survival analysis ----
# Overall survival stratified by risk group

km_fit <- survfit(Surv(os_time, os_event) ~ risk_group, data = surv_df)

# Log-rank test
logrank_p <- surv_pvalue(km_fit)$pval
message("Log-rank P = ", format.pval(logrank_p, digits = 3))

# KM plot
km_plot <- ggsurvplot(
  km_fit,
  data             = surv_df,
  pval             = TRUE,
  pval.size        = 4.5,
  conf.int         = TRUE,
  conf.int.alpha   = 0.15,
  risk.table       = TRUE,
  risk.table.height = 0.25,
  palette          = c("#377EB8", "#E41A1C"),
  xlab             = "Overall Survival (days)",
  ylab             = "Survival Probability",
  title            = "KM Curve — Gene Signature Risk Stratification",
  legend.title     = "Risk Group",
  legend.labs      = c("Low Risk", "High Risk"),
  ggtheme          = theme_minimal(base_size = 12)
)

pdf("../figures/km_curves.pdf", width = 8, height = 7)
print(km_plot, newpage = FALSE)
dev.off()

saveRDS(km_fit, "../results/Survival/km_fit.rds")

# KM for individual stages (subgroup analysis)
if (length(unique(surv_df$stage[!is.na(surv_df$stage)])) >= 2) {
  km_by_stage <- ggsurvplot(
    survfit(Surv(os_time, os_event) ~ stage, data = surv_df),
    data      = surv_df,
    pval      = TRUE,
    xlab      = "Overall Survival (days)",
    title     = "KM by AJCC Stage (All Patients)",
    palette   = "jco",
    ggtheme   = theme_minimal(base_size = 12)
  )
  pdf("../figures/km_by_stage.pdf", width = 8, height = 7)
  print(km_by_stage, newpage = FALSE)
  dev.off()
}

# ---- 6. Cox proportional hazards regression ----
# 6a. Univariate Cox — for each gene, the risk score, and clinical variables

univariate_results <- data.frame(
  variable  = character(),
  HR        = numeric(),
  CI_lower  = numeric(),
  CI_upper  = numeric(),
  p_value   = numeric(),
  stringsAsFactors = FALSE
)

# Risk score (continuous)
cox_rs <- coxph(Surv(os_time, os_event) ~ risk_score, data = surv_df)
s_rs   <- summary(cox_rs)
univariate_results <- rbind(univariate_results, data.frame(
  variable = "Risk Score (continuous)",
  HR       = s_rs$conf.int[1, "exp(coef)"],
  CI_lower = s_rs$conf.int[1, "lower .95"],
  CI_upper = s_rs$conf.int[1, "upper .95"],
  p_value  = s_rs$coefficients[1, "Pr(>|z|)"]
))

# Risk group (binary)
cox_rg <- coxph(Surv(os_time, os_event) ~ risk_group, data = surv_df)
s_rg   <- summary(cox_rg)
univariate_results <- rbind(univariate_results, data.frame(
  variable = "Risk Group (High vs. Low)",
  HR       = s_rg$conf.int[1, "exp(coef)"],
  CI_lower = s_rg$conf.int[1, "lower .95"],
  CI_upper = s_rg$conf.int[1, "upper .95"],
  p_value  = s_rg$coefficients[1, "Pr(>|z|)"]
))

# Age (continuous)
cox_age <- coxph(Surv(os_time, os_event) ~ age, data = surv_df)
s_age   <- summary(cox_age)
univariate_results <- rbind(univariate_results, data.frame(
  variable = "Age",
  HR       = s_age$conf.int[1, "exp(coef)"],
  CI_lower = s_age$conf.int[1, "lower .95"],
  CI_upper = s_age$conf.int[1, "upper .95"],
  p_value  = s_age$coefficients[1, "Pr(>|z|)"]
))

# Stage (if available)
stage_clean <- surv_df %>% filter(!is.na(stage), stage != "")
if (nrow(stage_clean) > 10) {
  stage_clean$stage_factor <- factor(stage_clean$stage, ordered = FALSE)
  cox_stage <- coxph(Surv(os_time, os_event) ~ stage_factor, data = stage_clean)
  s_stage   <- summary(cox_stage)
  # Use likelihood ratio test for multi-level factor
  univariate_results <- rbind(univariate_results, data.frame(
    variable = "AJCC Stage",
    HR       = NA,
    CI_lower = NA,
    CI_upper = NA,
    p_value  = s_stage$logtest["pvalue"]
  ))
}

saveRDS(univariate_results, "../results/Survival/cox_univariate.rds")
write.csv(univariate_results, "../results/Survival/cox_univariate.csv", row.names = FALSE)

# 6b. Multivariate Cox — risk score adjusted for clinical covariates
# Tests whether the gene signature is an INDEPENDENT prognostic factor.

cox_multi <- coxph(
  Surv(os_time, os_event) ~ risk_score + age + stage,
  data = surv_df %>% filter(!is.na(stage))
)
s_multi <- summary(cox_multi)

# Extract coefficients
multi_coefs <- data.frame(
  variable  = names(coef(cox_multi)),
  HR        = exp(coef(cox_multi)),
  CI_lower  = exp(confint(cox_multi)[, 1]),
  CI_upper  = exp(confint(cox_multi)[, 2]),
  p_value   = s_multi$coefficients[, "Pr(>|z|)"]
)

saveRDS(cox_multi,  "../results/Survival/cox_multivariate.rds")
saveRDS(multi_coefs, "../results/Survival/cox_multivariate_coefs.rds")
write.csv(multi_coefs, "../results/Survival/cox_multivariate.csv", row.names = FALSE)

message("Multivariate Cox — Risk Score HR: ",
        round(multi_coefs$HR[multi_coefs$variable == "risk_score"], 2),
        " P = ", format.pval(multi_coefs$p_value[multi_coefs$variable == "risk_score"]))

# ---- 7. Forest plot ----
# Combined univariate + multivariate results

forest_df <- bind_rows(
  univariate_results %>% mutate(model = "Univariate"),
  multi_coefs %>% rename(CI_lower = CI_lower, CI_upper = CI_upper,
                          HR = HR, p_value = p_value) %>%
    mutate(model = "Multivariate")
) %>%
  filter(variable != "(Intercept)")

forest_plot <- ggplot(forest_df, aes(x = HR, y = variable, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                  height = 0.2, position = position_dodge(width = 0.5)) +
  scale_x_log10() +
  scale_color_manual(values = c("Univariate" = "#377EB8", "Multivariate" = "#E41A1C")) +
  labs(
    title    = "Forest Plot — Prognostic Factors in HCC",
    subtitle = "Hazard Ratio (HR) with 95% CI",
    x        = "Hazard Ratio (log scale)",
    y        = "",
    color    = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("../figures/forest_plot.pdf", forest_plot, width = 9, height = 5)

# ---- 8. Time-dependent ROC analysis ----
# Evaluates the predictive accuracy of the risk score at specific follow-up
# times (1, 3, and 5 years). Unlike standard ROC, timeROC accounts for censoring.

time_points <- c(365, 1095, 1825)   # 1, 3, 5 years
names(time_points) <- c("1-year", "3-year", "5-year")

roc_res <- timeROC(
  T         = surv_df$os_time,
  delta     = surv_df$os_event,
  marker    = surv_df$risk_score,
  cause     = 1,
  weighting = "marginal",          # accounts for censoring distribution
  times     = time_points,
  ROC       = TRUE,
  iid       = TRUE                 # for CI estimation
)

saveRDS(roc_res, "../results/Survival/time_roc.rds")

# AUC values
for (i in seq_along(time_points)) {
  message(names(time_points)[i], " AUC: ",
          round(roc_res$AUC[i], 4),
          " [95% CI: ", round(confint(roc_res)$CI_AUC[i, 1], 4), "–",
          round(confint(roc_res)$CI_AUC[i, 2], 4), "]")
}

# ROC plot
pdf("../figures/time_roc.pdf", width = 7, height = 6)
plot(roc_res, time = time_points[1], col = "#377EB8", lwd = 2,
     title = "Time-Dependent ROC — Gene Signature Risk Score")
plot(roc_res, time = time_points[2], col = "#4DAF4A", lwd = 2, add = TRUE)
plot(roc_res, time = time_points[3], col = "#E41A1C", lwd = 2, add = TRUE)
legend("bottomright",
       legend = paste0(names(time_points), " AUC = ", round(roc_res$AUC, 3)),
       col    = c("#377EB8", "#4DAF4A", "#E41A1C"),
       lwd    = 2, bty = "n")
dev.off()

# ---- 9. Nomogram construction ----
# Integrates the gene signature risk score with clinical variables into a
# single graphical tool for individualized survival prediction.

if (length(unique(surv_df$stage[!is.na(surv_df$stage)])) >= 2) {

  # Use rms package for nomogram
  ddist <- datadist(surv_df %>% filter(!is.na(stage)));
  options(datadist = "ddist")

  cph_fit <- cph(
    Surv(os_time, os_event) ~ risk_score + age + stage,
    data     = surv_df %>% filter(!is.na(stage)),
    surv     = TRUE,
    x        = TRUE,
    y        = TRUE,
    time.inc = 365 * 3       # 3-year survival prediction
  )

  nom <- nomogram(
    cph_fit,
    fun  = function(x) 1 - x,   # convert to survival probability
    lp   = FALSE,
    funlabel      = "3-Year Survival Probability",
    conf.int      = c(0.1, 0.5, 0.9),
    fun.at        = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  )

  pdf("../figures/nomogram.pdf", width = 9, height = 7)
  plot(nom, xfrac = 0.25, lmgp = 0.3)
  dev.off()

  saveRDS(nom,    "../results/Survival/nomogram.rds")
  saveRDS(cph_fit, "../results/Survival/cph_model.rds")

  # ---- 10. Calibration plot ----
  # Compares nomogram-predicted survival with observed survival.

  cal <- calibrate(
    cph_fit,
    cmethod = "KM",
    method  = "boot",
    u       = 365 * 3,    # 3-year calibration
    m       = 50,          # subgroups of ~50 patients
    B       = 200          # bootstrap resamples
  )

  pdf("../figures/calibration.pdf", width = 6, height = 6)
  plot(cal, xlab = "Nomogram-Predicted 3-Year Survival",
       ylab = "Observed 3-Year Survival",
       main = "Calibration Plot")
  abline(0, 1, col = "grey60", lty = 2)
  dev.off()

  saveRDS(cal, "../results/Survival/calibration.rds")
}

# ---- 11. Summary ----
cat("\n============================================================\n")
cat(" Survival analysis complete.\n")
cat(" Signature genes:     ", length(sig_genes), "\n")
cat(" Log-rank P:          ", format.pval(logrank_p, digits = 3), "\n")
cat(" C-index (risk score):", round(s_rs$concordance[1], 4), "\n")
cat(" 1-year AUC:          ", round(roc_res$AUC[1], 4), "\n")
cat(" 3-year AUC:          ", round(roc_res$AUC[2], 4), "\n")
cat(" 5-year AUC:          ", round(roc_res$AUC[3], 4), "\n")
cat(" Output saved to results/Survival/\n")
cat("============================================================\n")
