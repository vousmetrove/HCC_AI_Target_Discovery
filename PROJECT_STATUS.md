# PROJECT_STATUS.md — Step-by-Step Completion Tracker

**Updated:** 2026-06-24 | **Overall Completion:** 18/23 steps (78%)

---

## Pipeline Execution Status

| Step | Task | Status | Date | Output Verified | Key File |
|------|------|--------|------|-----------------|----------|
| 01 | TCGA-LIHC Download | ✅ DONE | 06-23 | Yes | `data/raw/TCGA_LIHC_counts.rds` |
| 02 | Preprocessing | ✅ DONE | 06-23 | Yes | `data/processed/counts_filtered.rds` |
| 03 | DESeq2 DEG + VST + biomaRt | ✅ DONE | 06-23 | Yes | `results/DEG/DEG_list.rds` |
| 04 | GO + KEGG + GSEA Enrichment | ✅ DONE | 06-23 | Yes | `results/GO_KEGG/GO_enrichment_all.csv` |
| 05 | PPI Network + CytoHubba | ✅ DONE | 06-23 | Yes | `results/PPI/PPI_hub_genes_top10.csv` |
| 06 | LASSO Cox Feature Selection | ✅ DONE | 06-23 | Yes | `results/LASSO/lasso_genes.csv` |
| 07 | Random Forest (ranger) | ✅ DONE | 06-23 | Yes | `results/RandomForest/rf_importance.csv` |
| 08 | Survival Analysis (Cox + KM + timeROC) | ✅ DONE | 06-23 | Yes | `results/Survival/cox_univariate.csv` |
| 09 | GEO Validation (GSE76427) | ✅ DONE | 06-23 | Yes | `results/GEO_validation/GSE76427_validation_corrected.csv` |
| 10 | GEO Robustness (GAGE2A probe audit) | ✅ DONE | 06-24 | Yes | `results/GEO_validation/GAGE2A_probe_audit.csv` |
| 11 | Drug Target Prioritization | ✅ DONE | 06-24 | Yes | `results/Drug/drug_target_table.csv` |
| 12 | Immune Infiltration | ⚠️ PARTIAL | 06-24 | Limited | `results/Immune/` (marker genes filtered) |
| 13 | Molecular Docking Protocol | ✅ DONE | 06-24 | Protocol only | `results/Docking/docking_protocol.md` |
| 14 | Final Consensus Framework | ✅ DONE | 06-24 | Yes | `results/HubGene_Integration/` (7-dim scoring) |
| 15 | Drug Repurposing | ✅ DONE | 06-24 | Yes | `results/Drug/drug_repurposing_full.csv` |
| 16 | CMap Analysis | ✅ DONE | 06-24 | Yes | `results/Drug/cmap_connectivity_scores.csv` |
| 17 | Robustness Audit | ✅ DONE | 06-24 | Yes | `results/Validation/robustness_audit.md` |
| 18 | GSE14520 Validation | ✅ DONE | 06-24 | Yes | `results/GEO_validation/multi_cohort_validation.csv` |
| **19** | **ICGC-LIRI-JP Validation** | ⬜ PENDING | — | — | Requires portal download |
| **20** | **Immune: CIBERSORTx/xCell/EPIC** | ⬜ PENDING | — | — | Requires web submission |
| **21** | **Molecular Docking Execution** | ⬜ PENDING | — | — | Protocol ready; needs Vina |
| **22** | **Manuscript Writing** | ⬜ PENDING | — | — | `manuscript/outline.md` ready |
| **23** | **EI Conference Submission** | ⬜ PENDING | — | — | — |

---

## Evidence Summary

| Finding | Confidence | Evidence Dimensions |
|---------|------------|---------------------|
| PLK1 as therapeutic target | HIGH (0.85) | 7/7 positive |
| SPP1 as prognostic biomarker | HIGH (0.82) | Best GEO validation; 2/2 ML |
| CDK1 as mechanistic hub | MODERATE (0.68) | No GEO (probe absence) |
| BUB1B as emerging CIN target | LOW-MOD (0.49) | Preclinical drugs only |
| NR0B1 NOT targetable | CONFIRMED | Orphan nuclear receptor |
| GAGE2A deprioritized | CONFIRMED | GEO direction conflict |

## Bottlenecks

| # | Issue | Impact | Resolution |
|---|-------|--------|------------|
| 1 | ICGC portal access | GEO external validation incomplete | Download from dcc.icgc.org |
| 2 | Immune markers filtered | Step 12 incomplete | Use CIBERSORTx (web-based) |
| 3 | AutoDock Vina not installed | Step 13 pending | Install on workstation/HPC |
| 4 | No experimental validation | All findings computational | Frame as hypothesis-generating |

---

## Next Priority Task

**Step 19 — ICGC-LIRI-JP Validation** (P1)
- 240 HCC tumors with RNA-seq + survival data
- Available at https://dcc.icgc.org/projects/LIRI-JP
- Validate: PLK1, CDK1, SPP1
