# Evidence Framework — AI-Assisted Therapeutic Target Discovery for HCC

## Evidence Pyramid (9 Levels)

```
Level 9 ┌─────────────────────────┐
        │  Robustness Validation  │  ← Negative Control (1000×, P=0.001)
Level 8 │  Pharmacological         │  ← CMap Drug Perturbation (20 drugs scored)
Level 7 │  Druggability            │  ← DrugBank / TTD / DGIdb / OpenTargets
Level 6 │  External Validation     │  ← GSE76427 + GSE14520 Multi-Cohort
Level 5 │  Clinical Evidence       │  ← Cox Survival + KM + timeROC
Level 4 │  Machine Learning        │  ← LASSO (glmnet) + Random Forest (ranger)
Level 3 │  Network Evidence        │  ← STRING v12 PPI + 4 CytoHubba Algorithms
Level 2 │  Functional Enrichment   │  ← GO BP/CC/MF + KEGG + GSEA
Level 1 │  Differential Expression └── DESeq2 (|log2FC|>1, padj<0.05)
```

## Level Descriptions

### Level 1: Differential Expression
- **Method:** DESeq2 negative binomial GLM; Wald test
- **Threshold:** |log2FC| > 1, padj < 0.05
- **Output:** 6,301 DEGs (4,844 upregulated, 1,457 downregulated)
- **File:** `results/DEG/DEG_list.rds`

### Level 2: Functional Enrichment
- **Method:** enrichR (GO_BP/CC/MF_2023, KEGG_2021_Human) + fgsea (MSigDB Hallmark)
- **Key Finding:** Cell cycle #1 KEGG (P=2.3e-5); E2F Targets NES=+4.32; G2M Checkpoint NES=+4.36
- **File:** `results/GO_KEGG/`

### Level 3: Network Evidence
- **Method:** STRING v12 physical network (combined_score ≥ 400) + igraph + 4 CytoHubba algorithms
- **Key Finding:** 1,125 nodes, 2,906 edges; 5 consensus hub genes (CDK1 unique 4/4 algorithm hub)
- **File:** `results/PPI/`

### Level 4: Machine Learning Evidence
- **Method:** LASSO Cox (glmnet, 10-fold CV) + Random Forest (ranger, 1000 trees, permutation VIMP)
- **Key Finding:** 11-gene LASSO signature (C-index=0.654); RF top 30 VIMP; 7-gene consensus
- **File:** `results/LASSO/`, `results/RandomForest/`

### Level 5: Clinical Evidence
- **Method:** Cox univariate regression; Kaplan-Meier (median split); time-dependent ROC
- **Key Finding:** 12/12 candidate genes Cox-significant (P<0.05); PLK1 HR=1.372 (P=4.2e-6)
- **File:** `results/Survival/`

### Level 6: External Validation
- **Method:** GSE76427 (Illumina, limma) + GSE14520 (Affymetrix, limma)
- **Key Finding:** PLK1 validated 2/2 cohorts; SPP1 validated 2/2 cohorts; CDK1 no probe (platform limitation)
- **File:** `results/GEO_validation/`

### Level 7: Druggability
- **Method:** DrugBank, TTD, OpenTargets, ChEMBL, DGIdb database queries
- **Key Finding:** PLK1 Druggability 0.960 (Gold); CDK1 0.975 (Gold); SPP1 0.595 (Silver)
- **File:** `results/Drug/druggability_score.csv`

### Level 8: Pharmacological Evidence
- **Method:** CMap mechanism-based connectivity scoring (150 UP + 150 DOWN disease signature)
- **Key Finding:** Volasertib score=−35 (strongest proliferation reversal of 20 drugs)
- **File:** `results/Drug/cmap_connectivity_scores.csv`

### Level 9: Robustness Validation
- **Method:** Negative control (1000 random genes, same scoring framework) + independent audit
- **Key Finding:** PLK1 empirical P=0.001, Z=15.4, stable across 6 weight schemes
- **File:** `results/Validation/`

---

## Evidence Requirements by Claim Type

| Claim Type | Minimum Evidence Levels Required |
|------------|--------------------------------|
| "Gene X is a DEG in HCC" | Level 1 |
| "Gene X is enriched in cancer pathways" | Level 1 + 2 |
| "Gene X is a network hub" | Level 1 + 2 + 3 |
| "Gene X is prognostic" | Level 1 + 4 + 5 |
| "Gene X is a biomarker" | Level 1 + 4 + 5 + 6 |
| "Gene X is a therapeutic target" | Level 1 + 2 + 3 + 4 + 5 + 6 + 7 |
| "Drug Y is a repositioning candidate" | Level 7 + 8 + 9 |
