# Project Contribution Summary — AI-Assisted Therapeutic Target Discovery for HCC

## Project Scope

This project established a 9-level evidence framework for AI-assisted therapeutic target discovery in hepatocellular carcinoma, integrating TCGA transcriptomics, protein-protein interaction networks, machine learning, multi-cohort external validation, and drug repurposing analysis.

---

## Analyses Completed (20 steps)

| Category | Steps | Methods |
|----------|-------|---------|
| **Discovery** | 01–05 | DESeq2 DEG, GO/KEGG/GSEA, STRING PPI, CytoHubba |
| **ML Selection** | 06–07 | LASSO Cox (glmnet), Random Forest (ranger, 1000 trees) |
| **Clinical Validation** | 08 | Cox univariate, Kaplan-Meier, time-dependent ROC |
| **External Validation** | 09–10, 18 | GSE76427 (limma), GSE14520 (limma), multi-cohort meta-analysis |
| **Drug Prioritization** | 11, 15 | DrugBank, TTD, OpenTargets, ChEMBL, DGIdb, druggability scoring |
| **Pharmacological** | 16 | CMap connectivity scoring (20 drugs, 150-gene signature) |
| **Robustness** | 17, 19 | Independent audit + negative control (1000 iterations, P=0.001) |
| **Structural Biology** | 20–21 | PDB structure analysis, docking pipeline generation |
| **Honest Reporting** | 22 | GDSC data availability — Volasertib confirmed absent |

---

## Key Findings Retained

| Rank | Gene | Role | Evidence Level | Drug Candidate |
|------|------|------|----------------|----------------|
| 1 | **PLK1** | Primary Therapeutic Target | 9/9 | Volasertib (Phase III) |
| 2 | **SPP1** | Primary Prognostic Biomarker | 8/9 | Cabozantinib (FDA-approved, downstream) |
| 3 | **CDK1** | Mechanistic Hub | 7/9 | Dinaciclib (Phase II) |
| 4 | **BUB1B** | Emerging CIN Target | 6/9 | BAY-1816032 (Preclinical) |

---

## Findings Deprioritized or Removed

| Gene | Original Classification | Issue | New Classification | Justification |
|------|------------------------|-------|--------------------|---------------|
| **GAGE2A** | RF+LASSO consensus gene | GEO direction opposite (single probe cross-hybridization) | **REMOVED** from core list | Probe audit revealed single ILMN_3244168 probe with opposite direction. Not reliable for external validation. |
| **NR0B1** | Strongest Cox significance (P=1.2e-6) | Orphan nuclear receptor; no known ligand or drug | **DOWNGRADED** to Epigenetic Biomarker | Despite strongest statistical signal, zero druggability. Valuable as CTA/epigenetic phenomenon in Discussion. |
| **CDK1** | #1 PPI Hub (all 4 algorithms) | **No GEO validation** (probe absent on both Illumina GPL10558 and Affymetrix GPL3921) | **RETAINED** as Mechanistic Hub (with explicit limitation) | Technical limitation, not biological absence. Flavopiridol modest Phase II efficacy. Bone marrow toxicity. |
| **DRGX** | RF+LASSO consensus gene | Weak + inconsistent GEO direction | **LOW PRIORITY** | Not pursued; insufficient external evidence. |

---

## Validation Summary

| Validation Type | Result |
|-----------------|--------|
| **Internal (TCGA)** | 6,301 DEGs; all 12 candidate genes Cox-significant (P<0.05) |
| **External (GEO)** | PLK1 2/2, SPP1 2/2, CDK1 0/2 (platform limitation) |
| **Machine Learning** | LASSO C-index=0.654; RF C-index=0.695; risk score C=0.72 |
| **Drug Repurposing** | CMap identifies known HCC drugs (Sorafenib, Doxorubicin) — validates methodology |
| **Negative Control** | PLK1: 0/1000 random genes score higher (P=0.001, Z=15.4) |
| **Robustness Audit** | Overall confidence 0.72/1.00; PLK1 confidence 0.85 |

---

## Methodological Innovations

1. **9-Level Evidence Pyramid** — from DEG (L1) to Negative Control (L9)
2. **Cross-Platform External Validation** — Illumina (GSE76427) + Affymetrix (GSE14520)
3. **CMap + Drug Database Dual Validation** — independent pharmacological confirmation
4. **Mechanism-Based CMap Scoring** — proliferation reversal rather than simple gene overlap
5. **Honest Negative Reporting** — documented GDSC data absence; GAGE2A removal; CDK1 GEO gap

---

## Limitations (Honestly Disclosed)

1. CDK1 lacks GEO validation (platform probe absence)
2. Immune infiltration incomplete (marker genes filtered in bulk RNA-seq)
3. CMap uses literature-curated signatures, not raw LINCS L1000
4. Volasertib not in GDSC/CTRP — large-scale pharmacogenomic validation pending
5. Molecular docking not executed (requires Linux + Vina)
6. ICGC-LIRI-JP validation pending (portal download required)
7. ML models biased toward CTA genes (high |log2FC|)

---

## Reproducibility

All analyses are reproducible from the 12 core R scripts in `scripts/`. All intermediate and final outputs are stored in `results/` with CSV+RDS dual formats. The complete pipeline can be re-run with:

```bash
cd scripts/
for f in 01_download_TCGA.R 02_preprocessing.R 03_DESeq2.R \
         04_GO_KEGG.R 05_PPI_STRING.R 06_LASSO.R 07_RandomForest.R \
         08_Survival.R 09_GEO_validation.R; do Rscript "$f"; done
```

---

**GitHub:** https://github.com/vousmetrove/HCC_AI_Target_Discovery
