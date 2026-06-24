# 2. Methods

## 2.1 Study Design and Framework Overview

This study employed a three-stage AI-assisted evidence framework for therapeutic target discovery in hepatocellular carcinoma: (1) **Discovery** — transcriptomic profiling, functional enrichment, and protein interaction network analysis; (2) **Validation** — machine learning feature selection, survival analysis, and multi-cohort external validation; (3) **Prioritization** — drug target evaluation, Connectivity Map analysis, and multi-dimensional evidence scoring. A 9-level evidence pyramid was constructed, ranging from Level 1 (differential expression) to Level 9 (robustness validation via negative control testing). All analyses were performed in R (version 4.6.0). The complete pipeline, intermediate results, and documentation are publicly available at https://github.com/vousmetrove/HCC_AI_Target_Discovery.

## 2.2 Data Acquisition and Preprocessing

TCGA-LIHC RNA-seq data (STAR - Counts, hg38) were downloaded using TCGAbiolinks (v2.40). The dataset comprised 421 samples (371 primary tumors and 50 solid tissue normal samples) with 60,660 genes. Clinical annotations including overall survival time, vital status, AJCC pathologic stage, histologic grade, age, and sex were retrieved for 377 patients. Low-expression genes were filtered using edgeR::filterByExpr (min.count=10, min.total.count=15), retaining 22,730 genes. For analyses requiring protein-coding genes only, annotation was performed using biomaRt (Ensembl genes, hsapiens_gene_ensembl), yielding 3,845 unique protein-coding DEGs.

## 2.3 Differential Expression and Functional Enrichment

Differential expression analysis was performed using DESeq2 with the negative binomial generalized linear model and Wald test. Significance thresholds were set at |log2 fold change| > 1 and Benjamini-Hochberg adjusted P-value < 0.05. Variance-stabilizing transformation (VST) was applied for downstream machine learning and visualization. Gene Ontology (GO) enrichment was performed using enrichR (GO_Biological_Process_2023, GO_Cellular_Component_2023, GO_Molecular_Function_2023 databases) separately for upregulated (n=2,720) and downregulated (n=1,125) protein-coding DEGs. KEGG pathway enrichment used KEGG_2021_Human. Gene Set Enrichment Analysis (GSEA) was conducted using fgsea (v1.30) with MSigDB Hallmark gene sets, ranking genes by the DESeq2 Wald statistic.

## 2.4 Protein-Protein Interaction Network

A PPI network was constructed using STRING v12 (physical network, combined_score ≥ 400) for upregulated protein-coding DEGs. Network topology analysis was performed using igraph, identifying 1,125 nodes and 2,906 edges. Hub genes were identified using four CytoHubba algorithms: Maximal Clique Centrality (MCC), Degree, Closeness, and Betweenness centrality. Consensus hub genes were defined as those appearing in the Top 10 of at least three of the four algorithms.

## 2.5 Machine Learning Feature Selection

For prognostic feature selection, the top 300 protein-coding DEGs ranked by |log2FC| were used as input. LASSO Cox regression was performed using glmnet with 10-fold cross-validation, selecting the optimal lambda via partial likelihood deviance. An 11-gene prognostic signature was constructed, and risk scores were computed as the weighted sum of expression values. Random Survival Forest analysis was conducted using ranger with 1,000 trees and permutation-based variable importance. Hyperparameters (mtry, min.node.size) were tuned via grid search maximizing Harrell's C-index.

## 2.6 Survival Analysis

Univariate Cox proportional hazards regression was performed for each of the 12 consensus candidate genes to assess the association between gene expression (dichotomized at median) and overall survival. Kaplan-Meier survival curves were generated with the log-rank test. Time-dependent receiver operating characteristic (ROC) curves were computed at 1-, 3-, and 5-year time points using the timeROC package with inverse probability of censoring weighting.

## 2.7 External Validation

Two independent GEO cohorts were used for external validation. GSE76427 (Illumina HumanHT-12 v4, 115 HCC tumors vs. 52 adjacent non-tumor liver tissues) was processed using limma with normexp background correction and quantile normalization. Probe-to-gene mapping was performed using GPL10558 platform annotation. GSE14520 (Affymetrix U133A, 225 tumors vs. 220 non-tumor tissues) was processed using the same limma pipeline with GPL3921 annotation. Probe-level audits were conducted for discordant genes to identify cross-hybridization artifacts. For genes absent from microarray platforms (CDK1), the limitation was explicitly documented rather than imputed.

## 2.8 Drug Target Prioritization and Connectivity Map

Drug target evaluation integrated five public databases: DrugBank, Therapeutic Target Database (TTD), OpenTargets, ChEMBL, and DGIdb. Druggability scores were computed based on clinical development stage, number of compounds, co-crystal structure availability, and HCC-specific evidence. A Connectivity Map (CMap) analysis was performed using a disease signature comprising the top 150 upregulated and 150 downregulated protein-coding DEGs. Mechanism-based connectivity scoring was calculated for 20 drugs across PLK1, CDK1, and SPP1 targets, where negative scores indicate therapeutic reversal of the HCC proliferative transcriptional program. For drug candidates not included in GDSC/CTRP panels (Volasertib), the data gap was explicitly reported rather than filled with imputed values.

## 2.9 Robustness and Negative Control

An independent robustness audit was conducted for each major finding, systematically evaluating supporting evidence, contradictory evidence, assumptions, and potential biases. Confidence levels (High, Moderate, Low) were assigned based on the completeness and consistency of evidence. A negative control test randomly sampled 1,000 protein-coding genes from the expressed transcriptome (excluding the 12 candidate genes) and applied the identical scoring framework (DEG, PPI, Cox, Drug). An empirical P-value was computed as the proportion of random genes scoring equal to or higher than PLK1. Z-score, percentile, and bootstrap 95% confidence intervals were calculated. Sensitivity analysis was performed across six different weighting schemes to assess the stability of PLK1's ranking under alternative scoring assumptions.
