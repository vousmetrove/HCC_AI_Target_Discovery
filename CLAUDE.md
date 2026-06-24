# CLAUDE.md — AI-Assisted Therapeutic Target Discovery for HCC

## Project Identity

| Field | Value |
|-------|-------|
| **Project** | AI-Assisted Therapeutic Target Discovery Framework for HCC |
| **Full Title** | AI-Assisted Identification of Therapeutic Targets for Hepatocellular Carcinoma Through Integrated TCGA and GEO Transcriptomic Analysis |
| **Goal** | AI-assisted drug target identification + drug repurposing for HCC |
| **Target Venue** | EI-indexed conference paper |
| **Research Direction** | AI-assisted Drug Discovery |
| **Local Path** | `E:\analysis-bio\practice\HCC_AI_Target_Discovery` |
| **GitHub** | `git@github.com:vousmetrove/HCC_AI_Target_Discovery` |
| **Primary Language** | R (≥4.2), R 4.6.0 (Windows ucrt) |
| **Started** | 2026-06-23 |

---

## Scientific Design Philosophy

### Core Framework: Discovery → Validation → Prioritization

```
Stage 1: DISCOVERY           Stage 2: VALIDATION          Stage 3: PRIORITIZATION
─────────────────────────    ────────────────────────     ──────────────────────────
DEG Analysis (TCGA)          GEO External (GSE76427)     7-Dimension Evidence Scoring
GO/KEGG/GSEA Enrichment      GEO Multi-Cohort (GSE14520) Druggability Ranking
PPI Network + CytoHubba      Cox Survival Analysis       CMap Drug Repurposing
LASSO + Random Forest        Robustness Audit            Final Therapeutic Targets
                             Negative Control Test       Drug Repositioning Candidates
```

### Core Principles

1. **No single-analysis conclusions.** Every major finding must be supported by at least two independent evidence sources.
2. **Hub genes are NOT accepted solely by DEG ranking or PPI degree.** They must survive ML selection, survival validation, and external cohort confirmation.
3. **Every claim requires a complete evidence chain.** DEG → PPI → ML → Cox → GEO → Drug → CMap.
4. **Contradictory evidence is documented, not hidden.** GAGE2A removal, CDK1 GEO absence, NR0B1 undruggability are all explicitly recorded.
5. **Computational predictions are hypotheses, not clinical truths.** CMap suggests repositioning; it does not prove efficacy.
6. **Negative controls are mandatory.** Random gene testing confirms PLK1 is not an algorithmic artifact (Step 19, P=0.001).

### Evidence Levels (9-Tier Pyramid)

```
Level 9: Robustness Validation  ← Negative Control (1000×) + Independent Audit
Level 8: Pharmacological        ← CMap Drug Perturbation Signatures
Level 7: Druggability           ← DrugBank / TTD / DGIdb / OpenTargets
Level 6: External Validation    ← GSE76427 + GSE14520 Multi-Cohort
Level 5: Clinical Evidence      ← Cox Survival + KM + timeROC
Level 4: Machine Learning       ← LASSO (glmnet) + Random Forest (ranger)
Level 3: Network Evidence       ← STRING v12 PPI + 4 CytoHubba Algorithms
Level 2: Functional Enrichment  ← GO BP/CC/MF + KEGG + GSEA
Level 1: Differential Expression ← DESeq2 (|log2FC|>1, padj<0.05)
```

---

## Completed Steps (20/22)

| Step | Analysis | Status |
|------|----------|--------|
| 01-03 | TCGA Download + Preprocessing + DESeq2 DEG | ✅ |
| 04 | GO + KEGG + GSEA Enrichment | ✅ |
| 05 | STRING PPI Network + CytoHubba | ✅ |
| 06-07 | LASSO Cox + Random Forest (ranger) | ✅ |
| 08 | Cox Survival + KM + timeROC | ✅ |
| 09-10 | GEO Validation (GSE76427 + GSE14520) | ✅ |
| 11,15 | Drug Target Prioritization + Repurposing | ✅ |
| 16 | CMap Connectivity Analysis | ✅ |
| 17 | Robustness Audit | ✅ |
| 18 | Multi-Cohort Validation | ✅ |
| 19 | Negative Control Test (1000 iterations) | ✅ |
| 20-21 | Docking Reliability + Pipeline | ✅ |
| 22 | GDSC Drug Response (data not available) | ✅ |

---

## Current Top Targets

| Rank | Gene | Type | Evidence Level | Confidence | Drug (Phase) |
|------|------|------|----------------|------------|---------------|
| 1 | **PLK1** | Therapeutic Target | 9/9 | 0.85 | Volasertib (Phase III) |
| 2 | **SPP1** | Prognostic Biomarker | 8/9 | 0.82 | Cabozantinib (FDA, downstream) |
| 3 | **CDK1** | Mechanistic Hub | 7/9 | 0.68 | Dinaciclib (Phase II) |
| 4 | **BUB1B** | Emerging CIN Target | 6/9 | 0.49 | BAY-1816032 (Preclinical) |

## Deprioritized / Removed Genes

| Gene | Reason | Decision |
|------|--------|----------|
| GAGE2A | GEO direction opposite (single probe artifact) | Removed from core list |
| NR0B1 | Strongest Cox but orphan receptor, undruggable | Biomarker only |
| DRGX | Weak + inconsistent GEO direction | Low priority |

---

## Output File Locations
```
results/DEG/            # 6,301 DEGs; VST matrix; biomaRt annotations
results/GO_KEGG/        # GO + KEGG + GSEA enrichment
results/PPI/            # STRING v12 network; Cytoscape tables
results/LASSO/          # 11-gene LASSO signature
results/RandomForest/   # ranger VIMP ranking
results/Survival/       # Cox univariate + KM + timeROC
results/GEO_validation/ # GSE76427 + GSE14520 validation
results/Drug/           # Drug repurposing + CMap connectivity
results/Docking/        # PDB structures + docking pipeline
results/Validation/     # Robustness audit + negative control
results/Evidence/       # Gene evidence matrix
results/DrugResponse/   # GDSC data availability report
```

## Key References
- `docs/evidence_framework.md` — Evidence pyramid documentation
- `docs/project_summary.md` — Complete project contribution summary
- `docs/phd_application_summary.md` — PhD application statement
- `manuscript/storyline.md` — Discovery → Validation → Prioritization narrative
- `PROJECT_STATUS.md` — Step-by-step completion tracker
