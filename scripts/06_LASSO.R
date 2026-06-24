###############################################################################
# 06_LASSO.R
# ---------------------------------------------------------------------------
# Purpose:
#   Apply LASSO (Least Absolute Shrinkage and Selection Operator) penalized
#   Cox regression to the DEG expression matrix to select a compact, highly
#   prognostic gene signature for HCC. LASSO performs simultaneous feature
#   selection and regularization, shrinking the coefficients of non-informative
#   genes to exactly zero, yielding a sparse model with maximal predictive power.
#
# Why LASSO for biomarker discovery:
#   - High-dimensional setting (p ≫ n): thousands of DEGs but ~370 patients
#   - Automatic feature selection: irrelevant genes are eliminated
#   - Prevents overfitting via L1-penalty (λ tuning by cross-validation)
#   - Interpretable: each retained gene has a non-zero coefficient
#
# Expression matrix: DESeq2 rlog (variance-stabilized, log2-like scale)
#   Loaded from results/DEG/rlog_normalized.rds — the single normalized
#   expression matrix for the entire pipeline (no FPKM).
#
# Input from previous steps:
#   - results/DEG/rlog_normalized.rds  (Step 03 — DESeq2 normalized matrix)
#   - results/DEG/DEG_list.rds         (Step 03 — significant DEGs)
#   - results/PPI/hub_genes.rds        (Step 05 — PPI hub genes, optional)
#
# Workflow:
#   1. Load DESeq2 rlog expression matrix (primary normalized data)
#   2. Merge rlog expression with survival data (tumor samples only)
#   3. Split into training (70%) and internal test (30%) sets
#   4. 10-fold cross-validated LASSO-Cox via glmnet::cv.glmnet
#   5. Select optimal λ (lambda.min or lambda.1se)
#   6. Extract non-zero coefficient genes
#   7. Compute patient risk scores
#
# Output:
#   results/LASSO/lasso_model.rds         — fitted cv.glmnet object
#   results/LASSO/lasso_genes.rds         — selected gene list with coefficients
#   results/LASSO/lasso_coefficient_plot.pdf
#   results/LASSO/lasso_cv_plot.pdf
#
# Dependencies: glmnet, survival, dplyr, ggplot2, caret
###############################################################################

# ---- 0. Environment ----
library(glmnet)
library(survival)
library(dplyr)
library(tibble)
library(ggplot2)
library(caret)

set.seed(2024)
dir.create("../results/LASSO", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load data ----
# rlog_normalized.rds is the primary normalized expression matrix from DESeq2.
# It is variance-stabilized and on a log2-like scale — ideal for linear models
# (LASSO, Cox regression) and machine learning (Random Forest).
rlog_mat  <- readRDS("../results/DEG/rlog_normalized.rds")
clinical  <- readRDS("../data/processed/clinical_clean.rds")
deg_list  <- readRDS("../results/DEG/DEG_list.rds")
metadata  <- readRDS("../data/processed/metadata.rds")

# ---- 2. Prepare expression matrix (DEGs × tumor samples) ----
# Subset to DEGs and tumor samples only (LASSO is for prognostic modeling;
# normal tissue doesn't have survival outcomes)

tumor_barcodes <- metadata %>% filter(tissue == "Tumor") %>% pull(barcode)

# Match expression columns to tumor barcodes
tumor_cols <- intersect(colnames(rlog_mat), tumor_barcodes)

# Subset rlog expression to DEGs in tumor samples
deg_genes <- deg_list$gene_id
common_genes <- intersect(deg_genes, rownames(rlog_mat))

expr_tumor <- rlog_mat[common_genes, tumor_cols]   # DEGs × tumor samples

# Transpose → patients × genes (orientation needed for glmnet)
expr_mat <- t(expr_tumor)   # rows = patients, cols = genes

message("Expression matrix (tumor only, DESeq2 rlog): ", nrow(expr_mat),
        " patients × ", ncol(expr_mat), " DEGs")

# ---- 3. Merge with survival data ----
# Align expression patients with clinical records

patient_ids <- substr(rownames(expr_mat), 1, 12)
clin <- clinical %>%
  filter(submitter_id %in% patient_ids) %>%
  select(submitter_id, os_time, os_event)

# Match order
match_idx <- match(patient_ids, clin$submitter_id)

# Remove patients without clinical match
has_clin <- !is.na(match_idx)
expr_mat <- expr_mat[has_clin, ]
clin     <- clin[match_idx[has_clin], ]

stopifnot(nrow(expr_mat) == nrow(clin))

message("Final analysis set: ", nrow(expr_mat), " patients with survival data")

# ---- 4. Train-test split ----
# 70% training for model building, 30% for internal validation

train_idx <- createDataPartition(
  clin$os_event,
  p     = 0.7,
  list  = FALSE
)[, 1]

x_train <- expr_mat[train_idx, ]
x_test  <- expr_mat[-train_idx, ]
y_train <- Surv(clin$os_time[train_idx], clin$os_event[train_idx])
y_test  <- Surv(clin$os_time[-train_idx], clin$os_event[-train_idx])

message("Training set: ", nrow(x_train), " | Test set: ", nrow(x_test))

# ---- 5. LASSO-Cox regression ----
# cv.glmnet performs k-fold cross-validation over a grid of λ values to
# find the penalty that minimizes partial likelihood deviance.
#
# family = "cox": Cox proportional hazards model
# alpha   = 1:    LASSO penalty (L1); alpha = 0 would be ridge (L2)
# nfolds  = 10:   10-fold CV

cv_fit <- cv.glmnet(
  x        = x_train,
  y        = y_train,
  family   = "cox",
  alpha    = 1,
  nfolds   = 10,
  type.measure = "deviance",
  maxit    = 10000        # increase iterations for convergence
)

saveRDS(cv_fit, "../results/LASSO/lasso_model.rds")

# ---- 6. Identify optimal λ ----
# Two common choices:
#   lambda.min  — λ that minimizes CV deviance (more genes, potentially better fit)
#   lambda.1se  — largest λ within 1 SE of the minimum (more regularization, fewer genes)
# We default to lambda.1se for a parsimonious, more reproducible signature.

lambda_opt <- cv_fit$lambda.1se
message("Optimal λ (1se): ", round(lambda_opt, 4))
message("Optimal λ (min): ", round(cv_fit$lambda.min, 4))

# ---- 7. Extract selected genes ----
# Coefficients at optimal λ — non-zero entries correspond to selected genes.

coefs <- coef(cv_fit, s = "lambda.1se")
coef_mat <- as.matrix(coefs)
selected_idx <- which(coef_mat[, 1] != 0)

lasso_genes <- data.frame(
  gene_symbol = names(coef_mat[selected_idx, 1]),
  coefficient = coef_mat[selected_idx, 1]
) %>%
  arrange(desc(abs(coefficient)))

message("Genes selected by LASSO: ", nrow(lasso_genes))
print(lasso_genes)

saveRDS(lasso_genes, "../results/LASSO/lasso_genes.rds")
write.csv(lasso_genes, "../results/LASSO/lasso_genes.csv", row.names = FALSE)

# ---- 8. Risk score computation ----
# Risk score = Σ (expression_i × β_i) for i in selected genes
# Higher risk score → worse prognosis

# Training set risk scores
risk_train <- predict(cv_fit, newx = x_train, s = "lambda.1se", type = "link")
risk_train <- as.vector(risk_train)

# Test set risk scores
risk_test <- predict(cv_fit, newx = x_test, s = "lambda.1se", type = "link")
risk_test <- as.vector(risk_test)

# Save risk scores
risk_df <- data.frame(
  patient_id = c(rownames(x_train), rownames(x_test)),
  risk_score = c(risk_train, risk_test),
  dataset    = c(rep("train", length(risk_train)), rep("test", length(risk_test)))
)
saveRDS(risk_df, "../results/LASSO/risk_scores.rds")

# ---- 9. Visualization ----
# 9a. Cross-validation curve

pdf("../figures/lasso_cv_plot.pdf", width = 8, height = 6)
plot(cv_fit)
title("LASSO Cox — 10-fold Cross-Validation", line = 2.5)
dev.off()

# 9b. Coefficient path plot
lasso_fit <- cv_fit$glmnet.fit

pdf("../figures/lasso_coefficient_path.pdf", width = 8, height = 6)
plot(lasso_fit, xvar = "lambda", label = FALSE)
abline(v = log(lambda_opt), lty = 2, col = "red")
title("LASSO Coefficient Path", line = 2.5)
mtext(paste0("Selected genes at λ.1se: ", nrow(lasso_genes)),
      side = 3, line = 0.5, cex = 0.9)
dev.off()

# 9c. Coefficient bar plot
coef_plot <- ggplot(lasso_genes, aes(x = reorder(gene_symbol, coefficient),
                                      y = coefficient)) +
  geom_col(aes(fill = coefficient > 0), width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("TRUE" = "Risk (HR > 1)", "FALSE" = "Protective (HR < 1)"),
                    guide = guide_legend(title = NULL)) +
  coord_flip() +
  labs(
    title    = "LASSO-Selected Prognostic Genes",
    subtitle = paste0("λ.1se = ", round(lambda_opt, 4),
                      " | ", nrow(lasso_genes), " genes"),
    x        = "Gene",
    y        = "LASSO Coefficient (log HR)"
  ) +
  theme_minimal(base_size = 12)

ggsave("../figures/lasso_coefficient_bar.pdf", coef_plot, width = 8, height = 5)

# ---- 10. Summary ----
cat("\n============================================================\n")
cat(" LASSO feature selection complete.\n")
cat(" Genes selected:", nrow(lasso_genes), "\n")
cat("   Training patients:", nrow(x_train), "\n")
cat("   Test patients:    ", nrow(x_test), "\n")
lasso_genes %>% mutate(cat = ifelse(coefficient > 0, "Risk+", "Risk−")) %>%
  count(cat) %>% { cat("   ", .$cat[1], ":", .$n[1], " | ", .$cat[2], ":", .$n[2], "\n") }
cat(" Output saved to results/LASSO/\n")
cat("============================================================\n")
