# Independent Validation Audit — HCC_AI_Target_Discovery

**Date:** 2026-06-24 | **Auditor:** Systematic robustness assessment
**Scope:** All major findings from Steps 01–16

---

## A. PLK1 as Therapeutic Target

### Supporting Evidence (7/7 dimensions positive)

| # | Evidence | Source |
|---|----------|--------|
| 1 | DEG: log2FC=+3.49, padj=1.7e-66 | Step 03, `results/DEG/DEG_list.rds` |
| 2 | PPI Hub #3: Degree=47, 3/4 CytoHubba algorithms | Step 05, `results/PPI/ppi_node_metrics.rds` |
| 3 | TCGA Cox: HR=1.372 [1.20–1.57], P=4.2e-6 (highest HR among PPI genes) | Step 08, `results/Survival/cox_univariate.csv` |
| 4 | GEO validated: logFC=+0.158, P=2.6e-3, direction correct | Step 09, `results/GEO_validation/GSE76427_validation_corrected.csv` |
| 5 | Drug: Volasertib Phase III, BI2536 Phase II, Onvansertib Phase II | Step 11, `results/Drug/drug_target_table.csv` |
| 6 | CMap: Volasertib score=−35 (strongest proliferation reversal) | Step 16, `results/Drug/cmap_connectivity_scores.csv` |
| 7 | Druggability: 0.960 (Gold Tier) | Step 15, `results/Drug/druggability_score.csv` |

### Contradictory Evidence

| # | Issue | Severity |
|---|-------|----------|
| 1 | PLK1 is NOT HCC-specific — overexpressed in many cancers | Moderate — reduces novelty but not validity |
| 2 | GEO logFC is small (+0.158) — clinical significance of small effect? | Low — direction is correct; microarray dynamic range is limited |
| 3 | Not selected by LASSO or RF — ML methods favored CTA genes with higher |log2FC| | Moderate — PLK1 effect size is smaller than CTAs but more functionally relevant |
| 4 | Volasertib Phase III was in AML, not HCC — no HCC-specific Phase III data | Moderate — HCC Phase II exists but smaller |
| 5 | PLK1 inhibitors cause thrombocytopenia as on-target toxicity | Moderate — clinical limitation known from Phase I trials |

### Assumptions

1. PLK1 overexpression in HCC is functionally oncogenic (supported by literature but not directly tested here)
2. Volasertib's proliferation reversal in CMap translates to clinical efficacy (requires HCC-specific Phase III)
3. PLK1 expression measured by RNA-seq correlates with protein-level PLK1 activity (not verified by proteomics)

### Confidence: **HIGH** (0.85/1.00)

Justification: PLK1 has the most complete evidence package — positive in all 7 evaluation dimensions. Only CDK1 rivals it on druggability (0.975 vs 0.960), but PLK1 has stronger GEO validation and a more advanced drug (Phase III vs Phase II). The main limitation is lack of HCC-specific Phase III data for Volasertib.

### Potential Bias

- **Selection bias:** PLK1 was included as a PPI Hub gene candidate; other kinases with similar profiles may have been missed by not being in the PPI top 5.
- **Publication bias:** PLK1 has extensive literature; our confidence may be inflated by existing knowledge rather than novel discovery.
- **Platform bias:** GEO validation used GSE76427 only; second cohort (GSE14520) failed due to annotation issues.

---

## B. CDK1 as Therapeutic Target

### Supporting Evidence

| # | Evidence | Source |
|---|----------|--------|
| 1 | DEG: log2FC=+3.47, padj=1.4e-73 | Step 03 |
| 2 | PPI Hub #1: Degree=59, ALL 4 CytoHubba algorithms | Step 05 |
| 3 | TCGA Cox: HR=1.310 [1.15–1.50], P=8.1e-5 | Step 08 |
| 4 | Drug: Dinaciclib Phase II, Milciclib Phase II HCC, Flavopiridol Phase II HCC | Step 11 |
| 5 | CMap: Dinaciclib score=−24 (2nd strongest proliferation reversal) | Step 16 |
| 6 | Druggability: 0.975 (Gold Tier, highest score) | Step 15 |

### Contradictory Evidence

| # | Issue | Severity |
|---|-------|----------|
| 1 | **No GEO validation** — CDK1 has no probe on GPL10558 Illumina platform | **HIGH** — complete absence of external validation |
| 2 | CDK1 inhibitors (Flavopiridol) showed only MODEST efficacy in HCC Phase II | **HIGH** — clinical data contradicts strong computational prediction |
| 3 | On-target bone marrow toxicity — CDK1 is essential for hematopoiesis | Moderate — limits therapeutic window |
| 4 | Not selected by LASSO or RF — same CTA dominance issue as PLK1 | Low |
| 5 | Dinaciclib is a pan-CDK inhibitor, not CDK1-selective — therapeutic effect may not be CDK1-specific | Moderate |

### Assumptions

1. CDK1 network centrality (Degree=59) translates to functional importance (structural biology assumption)
2. Absence of GEO validation is purely technical (Illumina probe design) and not biological
3. Dinaciclib's CMap score reflects CDK1 inhibition rather than CDK2/5/9 off-target effects

### Confidence: **MODERATE** (0.68/1.00)

Justification: CDK1 has the strongest PPI evidence (unique 4/4 algorithm hub) and tied-highest druggability score (0.975). However, the COMPLETE ABSENCE of GEO validation and the MODEST clinical efficacy of Flavopiridol in HCC Phase II are significant red flags. The druggability score overestimates CDK1's clinical translatability because it penalizes GEO absence insufficiently.

### Potential Bias

- **Network bias:** CDK1's #1 PPI Hub status may be inflated by its role as a 'hub' in STRING networks (CDK1 interacts with many substrates). This is real biology but may overstate its position as a drug target relative to more selective targets.
- **Technology bias:** GEO absence is a TRUE limitation. We cannot rule out the possibility that CDK1 expression is less robust in external cohorts than our TCGA data suggests.
- **Overconfidence from druggability score:** The 0.975 score was driven by many compounds in clinical trials, but clinical FAILURES (Flavopiridol modest efficacy) were not penalized.

---

## C. SPP1 as Prognostic Biomarker

### Supporting Evidence

| # | Evidence | Source |
|---|----------|--------|
| 1 | Selected by BOTH LASSO and RF in top variable gene set | Step 06, 07 |
| 2 | TCGA Cox: HR=1.124 [1.07–1.18], P=9.8e-6 | Step 08 |
| 3 | GEO validated: logFC=+0.309, P=1.3e-11 (STRONGEST GEO signal) | Step 09 |
| 4 | 4/4 SPP1 probes show consistent positive direction in GSE76427 | Step 10, `results/GEO_validation/GAGE2A_probe_audit.csv` |
| 5 | Validated in 10+ independent HCC cohorts (PubMed literature) | Literature |
| 6 | Immune association: correlation with M2 macrophage infiltration (immunosuppressive) | Step 12 (limited data) |
| 7 | GEO2 (GSE14520): 225T/220N cohort available but annotation failed | Step 10 |

### Contradictory Evidence

| # | Issue | Severity |
|---|-------|----------|
| 1 | Drug development is PRECLINICAL only — not a drug target | Low for biomarker claim |
| 2 | Immune infiltration data is INCOMPLETE (marker genes filtered out) | Moderate — immune claim is weakly supported |
| 3 | Not in PPI top 5 hubs — low network centrality | Low — not required for biomarker status |
| 4 | Second GEO cohort (GSE14520) not successfully validated due to annotation failure | Low-Moderate |

### Assumptions

1. SPP1 RNA expression correlates with secreted osteopontin protein (not verified)
2. SPP1-high tumors have an immunosuppressive microenvironment (partially supported)
3. SPP1's prognostic value is independent of other clinical variables (multivariate Cox needed)

### Confidence: **HIGH** (0.82/1.00)

Justification: SPP1 has the strongest external validation evidence (GEO P=1.3e-11, consistent across 4 probes). Its role as a BIOMARKER (not drug target) is extremely well-supported. The main limitations are incomplete immune data and lack of multivariate Cox adjustment.

### Potential Bias

- **Literature bias:** SPP1 is well-studied in HCC — our findings confirm rather than discover. Novelty is low.
- **Incomplete validation:** GSE14520 failure is noted but not resolved.

---

## D. Volasertib as Repositioning Candidate

### Supporting Evidence
| # | Evidence |
|---|----------|
| 1 | CMap score: −35 (strongest of 20 drugs tested) |
| 2 | Phase III drug (most advanced among candidates) |
| 3 | Targets PLK1 (our #1 therapeutic target) |
| 4 | Reverses 15/25 HCC proliferation genes |
| 5 | HCC Phase II trial completed |

### Contradictory Evidence
| 1 | Phase III was in AML, NOT HCC |
| 2 | Thrombocytopenia as dose-limiting toxicity |
| 3 | CMap is computational — in vitro validation pending |

### Confidence: **MODERATE** (0.70/1.00)
### Potential Bias: CMap database is biased toward well-studied drugs; positive control overfitting possible.

---

## E. Dinaciclib as Repositioning Candidate

### Supporting Evidence
| 1 | CMap score: −24 (2nd strongest) |
| 2 | Phase II drug targeting CDK1 |
| 3 | HCC cell line validation |

### Contradictory Evidence
| 1 | NOT HCC-specific — pan-cancer Phase II |
| 2 | Pan-CDK (not CDK1-selective) — therapeutic effect may be multi-target |
| 3 | Similar pan-CDK inhibitor (Flavopiridol) showed modest HCC efficacy |
| 4 | CDK1 has NO GEO validation → weakness in the underlying target evidence |

### Confidence: **LOW-MODERATE** (0.55/1.00)
### Potential Bias: CMap may overestimate Dinaciclib's HCC effect by conflating CDK1/2/5/9 inhibition.

---

## F. CMap Conclusions — Overall Assessment

### Strengths
1. Volasertib's −35 score is internally consistent with PLK1's multi-dimensional evidence
2. Known HCC drugs (Sorafenib, Doxorubicin) are correctly identified — validates the methodology
3. CMap ranking correlates with drug development stage — not coincidental

### Weaknesses
1. CMap database is based on cell-line perturbation data — may not translate to human HCC
2. Limited to drugs in the CLUE LINCS database — novel compounds may be missed
3. Proliferation gene set used for scoring is heavily biased toward cell cycle genes
4. Connectivity scoring uses literature-curated drug-gene signatures, not raw LINCS data

### Confidence: **MODERATE** (0.65/1.00)

---

## Summary — All Findings Confidence Assessment

| # | Finding | Supporting | Contradictory | Confidence | Key Bias |
|---|---------|------------|---------------|------------|----------|
| A | PLK1 as therapeutic target | 7/7 dimensions | Non-specific to HCC; thrombocytopenia | **HIGH (0.85)** | Publication bias |
| B | CDK1 as therapeutic target | 6/7 dimensions | **No GEO validation**; modest Phase II | **MODERATE (0.68)** | Network centrality inflation |
| C | SPP1 as prognostic biomarker | Best GEO; both ML; 10+ cohorts | Incomplete immune data | **HIGH (0.82)** | Literature bias (low novelty) |
| D | Volasertib repositioning | CMap −35; Phase III; PLK1 target | AML not HCC; thrombocytopenia | **MODERATE (0.70)** | CMap overfitting |
| E | Dinaciclib repositioning | CMap −24; Phase II; CDK1 target | No GEO; pan-CDK; Flavopiridol weak | **LOW-MOD (0.55)** | Pan-CDK conflation |
| F | CMap overall conclusions | Identifies known HCC drugs correctly | Cell-line data; biased gene set | **MODERATE (0.65)** | Proliferation gene bias |

## Overall Project Confidence

| Dimension | Score |
|-----------|-------|
| DEG & pathway analysis | 0.90 — robust, well-validated |
| PPI network analysis | 0.85 — STRING v12 + 4 algorithms |
| ML feature selection | 0.70 — CTA dominance is a known limitation |
| Survival analysis | 0.80 — 12/12 genes Cox-significant |
| GEO external validation | 0.72 — 1/2 cohorts successful; CDK1 missing |
| Drug target prioritization | 0.75 — comprehensive but literature-based |
| CMap drug repurposing | 0.65 — methodology sound, limited by database |
| Immune infiltration | 0.40 — INCOMPLETE; marker genes filtered out |
| **PROJECT OVERALL** | **0.72** — Solid with acknowledged limitations |

## Critical Limitations Requiring Disclosure in Manuscript

1. **CDK1 has no external validation** (GEO). This must be explicitly stated.
2. **Immune infiltration analysis is incomplete** (filtered-out markers). Acknowledge as limitation.
3. **CMap uses literature-curated signatures**, not raw LINCS L1000 data. Disclose methodology.
4. **ML models were dominated by CTA genes** (high |log2FC|). Note this bias.
5. **GAGE2A was deprioritized** after GEO direction mismatch. Document the decision.
6. **NR0B1 has strongest Cox P but is not druggable.** Explain the discrepancy.

## Recommendations for Revision

1. **PLK1** — Retain as primary therapeutic target. Strongest overall evidence.
2. **SPP1** — Retain as primary prognostic biomarker. Emphasize biomarker role, not drug target.
3. **CDK1** — Down-rank to "mechanistic hub". Add GEO limitation prominently.
4. **Volasertib** — Present with caveat: Phase III in AML, HCC Phase II completed. Bridging needed.
5. **Dinaciclib** — Present as secondary candidate. Note pan-CDK limitation and Flavopiridol precedent.
6. **Immune analysis** — Either complete with CIBERSORTx or remove from abstract. Currently incomplete.

---
*Audit prepared for EI Conference 2026 submission. All limitations are honestly disclosed.*
*Data provenance: All evidence traceable to `results/` folder and specific CSV/RDS files.*
