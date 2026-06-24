###############################################################################
# 01_download_TCGA.R
# ---------------------------------------------------------------------------
# Purpose:
#   Download Hepatocellular Carcinoma (HCC) transcriptomic and clinical data
#   from The Cancer Genome Atlas (TCGA) via the GDC Data Portal using the
#   TCGAbiolinks R/Bioconductor package.
#
# Data retrieved:
#   1. RNA-seq gene expression (HTSeq-Counts) — PRIMARY expression matrix
#      - Used for DESeq2 differential expression analysis (requires raw integer
#        counts for the negative binomial model)
#      - After DESeq2, the regularized log (rlog) transformation provides a
#        variance-stabilized, normalized expression matrix suitable for all
#        downstream analyses: machine learning (LASSO, Random Forest),
#        visualization (heatmaps, PCA), and survival modeling
#      - This keeps the entire pipeline on a single, coherent expression
#        modality — no FPKM nor cross-normalization artifacts
#   2. Clinical & follow-up data — survival time, vital status, tumor stage,
#      histologic grade, age, sex, and other prognostic covariates
#
# Project: TCGA-LIHC (Liver Hepatocellular Carcinoma)
# Samples: ~374 HCC tumor + ~50 adjacent normal tissues
#
# Output:
#   data/raw/TCGA_LIHC_counts.rds     — raw HTSeq count matrix (genes × samples)
#   data/raw/TCGA_LIHC_clinical.rds   — clinical annotations
#
# Why HTSeq-Counts (not FPKM) as primary:
#   - DESeq2 requires raw integer counts; FPKM breaks the count-variance
#     relationship the model relies on
#   - rlog/vst transformations from DESeq2 produce properly normalized values
#     that account for library size, composition bias, and mean-variance
#     dependency — superior to log2(FPKM+1) for downstream ML
#   - Using a single expression modality eliminates cross-normalization
#     discrepancies and simplifies the analytical provenance
#
# Dependencies: TCGAbiolinks, SummarizedExperiment, dplyr
###############################################################################

# ---- 0. Environment setup ----
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

set.seed(2024)
dir.create("../data/raw", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Query TCGA-LIHC RNA-seq (HTSeq-Counts) — PRIMARY ----
# GDC data category: Transcriptome Profiling
# Experimental strategy: RNA-Seq
# Workflow type:   HTSeq-Counts
#
# HTSeq-Counts produce integer read counts per gene. These are the required
# input for DESeq2's negative binomial GLM (counts ~ NB(μ, α)). DESeq2
# internally performs median-of-ratios normalization and produces rlog/vst
# transformations that serve as the normalized expression matrix for all
# downstream machine learning and visualization tasks.

query_counts <- GDCquery(
  project           = "TCGA-LIHC",
  data.category     = "Transcriptome Profiling",
  data.type         = "Gene Expression Quantification",
  workflow.type     = "STAR - Counts",
  sample.type       = c("Primary Tumor", "Solid Tissue Normal"),
)

# Download the data from GDC
GDCdownload(query_counts, method = "api", files.per.chunk = 10)

# Prepare into a SummarizedExperiment object
#   - assay: raw integer count matrix
#   - rowData: gene-level annotations (Ensembl ID, gene symbol, gene type)
#   - colData: sample-level metadata (barcode, sample type, patient ID)
counts_data <- GDCprepare(
  query_counts,
  summarizedExperiment = TRUE,
  save                 = TRUE,
  save.filename        = "../data/raw/TCGA_LIHC_counts.rds"
)

message("HTSeq-Count matrix dimensions: ", paste(dim(counts_data), collapse = " × "))

# ---- 2. Download clinical data ----
# Retrieves patient demographics, tumor staging (AJCC TNM), histologic grade,
# survival outcomes (OS, DSS, PFI), treatment history, and lab values.
#
# Two clinical tables are downloaded:
#   clinical:          curated clinical data (structured fields)
#   biospecimen:       supplemental biospecimen records (tissue source, etc.)

clinical_query <- GDCquery_clinic(
  project      = "TCGA-LIHC",
  type         = "clinical",
  save.csv     = FALSE
)

clinical_supplement <- GDCquery_clinic(
  project  = "TCGA-LIHC",
  type     = "biospecimen",
  save.csv = FALSE
)

# Save clinical data
saveRDS(clinical_query,      "../data/raw/TCGA_LIHC_clinical.rds")
saveRDS(clinical_supplement, "../data/raw/TCGA_LIHC_biospecimen.rds")

message("Clinical data: ", nrow(clinical_query), " patients")

# ---- 3. Summary ----
cat("\n============================================================\n")
cat(" TCGA-LIHC data download complete.\n")
cat("\n Primary expression matrix:\n")
cat("   data/raw/TCGA_LIHC_counts.rds    (HTSeq raw counts)\n")
cat("\n Clinical annotations:\n")
cat("   data/raw/TCGA_LIHC_clinical.rds  (clinical table)\n")
cat("\n Pipeline notes:\n")
cat("   - Raw counts → 02_preprocessing → DESeq2 (03)\n")
cat("   - DESeq2 rlog → normalized matrix for all downstream ML\n")
cat("============================================================\n")
