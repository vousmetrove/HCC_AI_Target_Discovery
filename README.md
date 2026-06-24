# AI-Assisted Therapeutic Target Discovery Framework for Hepatocellular Carcinoma

## An Evidence-Based Multi-Layer Pipeline Integrating Transcriptomics, Network Biology, Machine Learning, and Drug Repurposing

[![Status](https://img.shields.io/badge/status-analysis%20complete-brightgreen)](PROJECT_STATUS.md)
[![Pipeline](https://img.shields.io/badge/pipeline-20%2F22%20steps-9cf)](PROJECT_STATUS.md)
[![Evidence](https://img.shields.io/badge/evidence-9%20levels-blue)](docs/evidence_framework.md)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

### Project Overview

This repository implements an **AI-assisted multi-layer evidence framework** for therapeutic target discovery and drug repurposing in Hepatocellular Carcinoma (HCC). The framework operates through three sequential stages:

```
Discovery ──→ Validation ──→ Prioritization
```

**Key Features:**
- **9-Level Evidence Pyramid** — from differential expression to negative control validation
- **Multi-Omics Integration** — TCGA transcriptomics + STRING PPI network + CytoHubba topology
- **Machine Learning** — LASSO Cox regression + Random Survival Forest (ranger) feature selection
- **Multi-Cohort Validation** — Independent external validation in GSE76427 + GSE14520
- **Drug Repurposing** — DrugBank/TTD/DGIdb druggability assessment + Connectivity Map (CMap) analysis
- **Robustness Audit** — Independent validation audit + negative control test (1000 random genes)
- **Honest Reporting** — Contradictory evidence and data gaps are documented, not hidden

**Key Findings:**
- **PLK1** — Top therapeutic target (9/9 evidence levels; Volasertib Phase III; CMap score = −35)
- **SPP1** — Top prognostic biomarker (validated in 2/2 GEO cohorts; both ML methods selected)
- **CDK1** — Central PPI hub (degree=59, all 4 CytoHubba algorithms; GEO limitation acknowledged)

---

### Directory Structure

```
HCC_AI_Target_Discovery/
├── .gitignore
├── README.md
├── data/
│   ├── raw/                  # Raw TCGA STAR counts + clinical + GPL annotations
│   └── processed/            # Filtered counts, clinical, metadata
├── scripts/                  # 12 R analysis scripts (01–12), 4,106 lines total
├── results/
│   ├── DEG/                  # Step 03: Differential expression (DESeq2)
│   ├── GO_KEGG/              # Step 04: GO + KEGG + GSEA enrichment
│   ├── PPI/                  # Step 05: STRING PPI network + hub genes
│   ├── LASSO/                # Step 06: LASSO-Cox feature selection
│   ├── RandomForest/         # Step 07: Random Forest VIMP ranking (ranger)
│   ├── Survival/             # Step 08: KM, Cox, time-ROC
│   ├── GEO_validation/       # Step 09–10: GEO external validation
│   ├── Drug/                 # Step 11,15,16: Drug repurposing + CMap
│   ├── Validation/           # Step 17: Robustness audit
│   ├── Immune/               # Step 12: Immune infiltration (pending completion)
│   ├── HubGene_Integration/  # Step 14: Final consensus (pending)
│   ├── enrichment/           # Gene symbol list for GO
│   └── report/               # Annotation summary
├── figures/                  # 41 publication-quality figures (PDF + PNG)
├── docs/                     # Project progress report
└── manuscript/               # EI conference paper materials
    ├── outline.md            # Complete paper outline
    ├── tables/               # Supplementary tables (DEG Top50, etc.)
    ├── figures/              # Publication figures placeholder
    ├── references/           # Bibliography placeholder
    └── draft/                # Draft sections placeholder
```

---

### Analytical Workflow — 12 Steps

```
  TCGA-LIHC HTSeq Counts
         │
   [01] Download ───────────────── TCGAbiolinks
         │
   [02] Preprocessing ──────────── edgeR::filterByExpr, clinical cleaning
         │
   [03] DESeq2 DEG + rlog ──────── PRIMARY normalized matrix for all downstream ML
         │
         ├── DEG list ──────────────────────────┐
         │                                       │
   [04] GO / KEGG Enrichment ◄──────────────────┤
         │                                       │
   [05] STRING PPI Network ◄────────────────────┤  hub genes
         │                                       │
   [06] LASSO Cox Feature Selection ◄───────────┤  rlog matrix
         │                                       │
   [07] Random Forest VIMP (ranger) ◄───────────┤  rlog matrix
         │                                       │
         ├── LASSO ∩ RF ∩ PPI consensus ────────┤
         │                                       │
   [08] Survival Analysis ◄─────────────────────┤  KM, Cox, time-ROC
         │                                       │
   [09] GEO External Validation ◄───────────────┤  frozen model
         │                                       │
   [10] Immune Infiltration (ssGSEA) ◄──────────┤  immune microenvironment
         │                                       │
   [11] Drug Repurposing (DGIdb) ◄──────────────┤  druggability
         │                                       │
   [12] 9-Dimension Final Integration ──────────┤  AI-assisted prioritization
         │
    Final Therapeutic Target Ranking
```

---

#### Step 1 — Data Acquisition (`01_download_TCGA.R`)
- **Source**: TCGA-LIHC project via [GDC Data Portal](https://portal.gdc.cancer.gov/)
- **Primary data**: HTSeq-Counts RNA-seq expression matrix (raw integer counts — single modality)
- **Clinical data**: survival time, vital status, AJCC stage, histologic grade, age, sex
- **External validation**: GEO datasets (GSE14520, GSE76427, ICGC-LIRI-JP)
- **Key package**: TCGAbiolinks

#### Step 2 — Data Preprocessing (`02_preprocessing.R`)
- Low-expression gene filtering via **edgeR::filterByExpr**
- Sample type annotation (tumor vs. normal) from TCGA barcodes
- Clinical data cleaning: OS time/event derivation, AJCC stage collapsing
- Exploratory PCA on log2(CPM)
- **No FPKM** — counts are the single primary expression modality
- Output → `data/processed/`

#### Step 3 — DEG + Normalized Matrix (`03_DESeq2.R`)
- **DESeq2** negative binomial model on HTSeq raw counts (tumor vs. normal)
- DEG thresholds: |log2 fold change| > 1, adjusted P-value < 0.05
- **rlog transformation** → `results/DEG/rlog_normalized.rds`
  - This variance-stabilized matrix feeds ALL downstream ML (LASSO, RF, Survival, Immune)
- Volcano plot, top-50 DEG heatmap, publication PCA (rlog)
- Output → `results/DEG/`

#### Step 4 — Functional Enrichment (`04_GO_KEGG.R`)
- **clusterProfiler** for GO (BP, CC, MF) and KEGG pathway enrichment
- ORA on DEGs + GSEA on ranked gene list
- GO semantic similarity simplification
- Visualization: dot plots, cnet plots, GSEA ridge plots
- Identifies cancer-relevant pathways for target contextualization
- Output → `results/GO_KEGG/`

#### Step 5 — PPI Network Analysis (`05_PPI_STRING.R`)
- **STRING v12** database query for DEG protein interactions
- Network topology metrics: degree, betweenness, closeness, eigenvector centrality, MCC
- Composite hub score from integrated ranking across all metrics
- Hub gene identification (top 10% by composite score)
- Cytoscape-compatible node and edge table exports
- PPI sub-network visualization (ggraph)
- Output → `results/PPI/`

#### Step 6 — LASSO Feature Selection (`06_LASSO.R`)
- **glmnet** LASSO-Cox regression on DEG rlog expression
- 10-fold cross-validation to select optimal λ (lambda.1se for parsimony)
- Extracts genes with non-zero coefficients as prognostic signature
- Train-test split (70/30) for internal validation
- Risk score construction: Σ(expression_i × LASSO coefficient_i)
- Output → `results/LASSO/`

#### Step 7 — Random Forest Importance Ranking (`07_RandomForest.R`)
- **ranger** (fast C++ RF implementation) for survival-based variable importance
- Grid search hyperparameter tuning: mtry, min.node.size
- Permutation VIMP ranking + Altmann P-values for robustness
- Consensus analysis: LASSO ∩ RF ∩ PPI (three-way intersection)
- OOB C-index for model performance
- Output → `results/RandomForest/`

#### Step 8 — Survival Analysis (`08_Survival.R`)
- **survival** / **survminer** for Kaplan-Meier curves (high- vs. low-risk)
- Univariate + multivariate Cox proportional hazards regression
- Time-dependent ROC (1-, 3-, 5-year AUC) via **timeROC**
- Nomogram with clinical variables (stage, age) via **rms**
- Calibration plots for nomogram validation
- Forest plot of hazard ratios across all prognostic factors
- Output → `results/Survival/`

#### Step 9 — GEO External Validation (`09_GEO_validation.R`)
- Download independent GEO cohorts (GSE14520, GSE76427)
- Probe-to-gene mapping with max-mean collapse
- **Frozen model**: TCGA-derived coefficients applied without refitting
- KM survival curves, Cox HR, time-ROC in external cohorts
- Random-effects meta-analysis across all cohorts (metafor)
- Combined forest plot: TCGA + GEO
- Output → `results/GEO_validation/`

#### Step 10 — Immune Infiltration (`10_Immune_Infiltration.R`)
- **ssGSEA** for 24 immune cell types (Bindea et al. signatures)
- **GSVA** for immune-related Hallmark pathways (IFN-γ, TNF-α, IL-6/JAK/STAT3, etc.)
- **ESTIMATE**: StromalScore, ImmuneScore, TumorPurity
- Hub gene–immune cell Spearman correlation analysis
- Differential immune infiltration: high- vs. low-risk groups (Wilcoxon)
- Correlation heatmaps and boxplots
- Output → `results/Immune/`

#### Step 11 — Drug Repurposing (`11_Drug_Repurposing.R`)
- **DGIdb v4 REST API** for drug-gene interaction queries
- Drug annotation by approval status and oncology category
- Polypharmacology scoring: more hub gene targets → higher score
- DrugBank / literature evidence integration
- Bipartite drug-target network construction (igraph)
- Drug ranking by composite score (target count × approval status)
- HCC-specific drug filtering (sorafenib, lenvatinib, etc.)
- Output → `results/Drug/`

#### Step 12 — Final Integration (`12_HubGene_Integration.R`)
- **9 evidence dimensions** combined into a single composite priority score:
  | # | Dimension | Weight | Source |
  |---|-----------|--------|--------|
  | 1 | DEG significance | 10% | Step 03 |
  | 2 | Pathway enrichment | 5% | Step 04 |
  | 3 | PPI network topology | 10% | Step 05 |
  | 4 | LASSO selection | 12% | Step 06 |
  | 5 | Random Forest VIMP | 12% | Step 07 |
  | 6 | Consensus (LASSO ∩ RF) | 6% | Step 07 |
  | 7 | Survival impact | 12% | Step 08 |
  | 8 | GEO external validation | 8% | Step 09 |
  | 9 | Immune correlation | 10% | Step 10 |
  | 10 | Druggability (DGIdb) | 15% | Step 11 |
- Tier classification: High Priority (top 10%), Moderate (top 30%), Lower, Candidate
- Evidence heatmap + composite score bar plot
- Supplementary table for manuscript
- Output → `results/HubGene_Integration/`

---

### Computational Environment

| Category | Tools / Packages |
|---|---|
| Language | R (≥ 4.2) |
| Data Access | TCGAbiolinks, GEOquery |
| DEG | DESeq2, edgeR |
| Enrichment | clusterProfiler, org.Hs.eg.db, enrichplot |
| PPI Network | STRINGdb, igraph, ggraph, tidygraph |
| ML Feature Selection | glmnet (LASSO), ranger (Random Forest) |
| Survival | survival, survminer, timeROC, rms |
| Meta-Analysis | metafor |
| Immune | GSVA, estimate, immunedeconv |
| Drug | httr, jsonlite (DGIdb API), igraph |
| Visualization | ggplot2, ComplexHeatmap, pheatmap, ggrepel |
| Integration | dplyr, tidyr, ComplexHeatmap |

---

### Quick Start

```bash
# 01 — Data download
Rscript scripts/01_download_TCGA.R

# 02 — Preprocessing
Rscript scripts/02_preprocessing.R

# 03 — DEG + rlog normalized matrix
Rscript scripts/03_DESeq2.R

# 04 — GO/KEGG enrichment
Rscript scripts/04_GO_KEGG.R

# 05 — STRING PPI network
Rscript scripts/05_PPI_STRING.R

# 06 — LASSO feature selection
Rscript scripts/06_LASSO.R

# 07 — Random Forest VIMP ranking
Rscript scripts/07_RandomForest.R

# 08 — Survival analysis
Rscript scripts/08_Survival.R

# 09 — GEO external validation
Rscript scripts/09_GEO_validation.R

# 10 — Immune infiltration
Rscript scripts/10_Immune_Infiltration.R

# 11 — Drug repurposing
Rscript scripts/11_Drug_Repurposing.R

# 12 — Final 9-dimension integration
Rscript scripts/12_HubGene_Integration.R
```

---

### Data Flow Summary

```
HTSeq Counts ──→ [02] filtered counts ──→ [03] DESeq2
                                              │
                            ┌─────────────────┼─────────────────┐
                            │                 │                 │
                      DEG_list.rds    rlog_normalized.rds    volcano/PCA
                            │                 │
               ┌────────────┤        ┌────────┼────────┬──────────┐
               │            │        │        │        │          │
          [04] GO/KEGG  [05] PPI  [06]LASSO [07]RF  [08]Surv   [10]Immune
               │            │        │        │        │          │
               └────────────┴────────┴────────┴────────┴──────────┘
                                        │
                                  [09] GEO Validation
                                        │
                                  [11] Drug Repurposing
                                        │
                                  [12] 9-Dimension Integration
                                        │
                                  Final Target Ranking
```

---

### Citation

> *Manuscript under preparation for EI Conference 2026.*

### License

MIT
