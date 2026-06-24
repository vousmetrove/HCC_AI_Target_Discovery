# CLAUDE.md — HCC AI-Assisted Drug Target Discovery

## Project Identity

| Field | Value |
|-------|-------|
| **Project** | HCC AI-Assisted Drug Target Discovery |
| **Full Title** | AI-Assisted Identification of Therapeutic Targets for Hepatocellular Carcinoma Through Integrated TCGA and GEO Transcriptomic Analysis |
| **Goal** | AI-assisted drug target identification + drug repurposing for HCC |
| **Target Venue** | EI-indexed conference paper |
| **Research Direction** | AI-assisted Drug Discovery |
| **Local Path** | `E:\analysis-bio\practice\HCC_AI_Target_Discovery` |
| **GitHub** | `git@github.com:vousmetrove/HCC_AI_Target_Discovery` |
| **Primary Language** | R (≥4.2), R 4.6.0 (Windows ucrt) |
| **Started** | 2026-06-23 |

---

## Directory Structure

```
data/raw/         # TCGA STAR counts, clinical, GPL annotations
data/processed/   # Filtered counts, clinical_clean, metadata
scripts/          # 12 core R scripts
results/          # 14 result directories
figures/          # Organized by module (DEG/GO_KEGG/LASSO/RF/PPI/Survival/Drug/Publication)
manuscript/       # outline.md + figures/ + tables/ + references/ + draft/
tables/           # Publication-ready summary tables
docs/             # Project progress report
logs/             # Per-session research logs (YYYY-MM-DD.md)
environment/      # Package versions, sessionInfo
CLAUDE.md         # This file — AI agent memory
PROJECT_STATUS.md # Step-by-step completion tracker
```

---

## Completed Steps

| Step | Analysis | Status | Key Output |
|------|----------|--------|------------|
| 01 | TCGA-LIHC Download | ✅ | 60,660 genes × 421 samples (STAR counts) |
| 02 | Preprocessing | ✅ | 22,730 genes × 416 samples; 372 patients |
| 03 | DESeq2 DEG + VST | ✅ | 6,301 DEGs (4,844 up / 1,457 down) |
| 04 | GO + KEGG + GSEA | ✅ | Cell cycle #1 KEGG; E2F NES=+4.32; G2M NES=+4.36 |
| 05 | PPI Network + CytoHubba | ✅ | 1,125 nodes, 2,906 edges; 5 hub genes |
| 06 | LASSO Cox | ✅ | 11-gene signature; C-index=0.654 |
| 07 | Random Forest (ranger) | ✅ | Top 30 VIMP; C-index=0.695 |
| 08 | Survival Analysis | ✅ | 12/12 genes Cox-significant (P<0.05) |
| 09 | GEO Validation (GSE76427) | ✅ | PLK1 validated (P=2.6e-3); SPP1 validated (P=1.3e-11) |
| 10 | GEO Robustness (GAGE2A probe audit) | ✅ | GAGE2A deprioritized (single probe, wrong direction) |
| 11 | Drug Target Prioritization | ✅ | PLK1 (0.960), CDK1 (0.975) Gold Tier |
| 12 | Immune Infiltration | ⚠️ Partial | Marker genes filtered out (bulk RNA-seq limitation) |
| 13 | Molecular Docking Protocol | ✅ | Protocol written; not executed (needs Vina) |
| 14 | Final Consensus Framework | ✅ | 7-dimension scoring; PLK1=0.845 (Tier 1) |
| 15 | Drug Repurposing | ✅ | Volasertib Phase III; 21 drugs evaluated |
| 16 | CMap Analysis | ✅ | Volasertib score=-35 (strongest); 20 drugs scored |
| 17 | Robustness Audit | ✅ | All findings audited; overall confidence 0.72/1.00 |
| 18 | GSE14520 Validation | ✅ | PLK1 P=4.5e-31; SPP1 P=1.3e-21; CDK1 no probe |

---

## Remaining Tasks (Priority Order)

| Priority | Step | Task | Status |
|----------|------|------|--------|
| **P1** | 19 | **ICGC-LIRI-JP Validation** | ⬜ Requires portal download |
| **P2** | 20 | **Immune Infiltration (CIBERSORTx/xCell/EPIC)** | ⬜ Requires web submission |
| **P3** | 21 | **Molecular Docking Execution (AutoDock Vina)** | ⬜ Protocol ready; needs software |
| **P4** | 22 | **Manuscript Writing** | ⬜ outline.md ready |
| **P5** | 23 | **EI Conference Submission** | ⬜ |

---

## Current Top Targets

| Rank | Gene | Type | Confidence | Drug (Phase) |
|------|------|------|------------|---------------|
| 1 | **PLK1** | Therapeutic Target | 0.85 | Volasertib (Phase III) |
| 2 | **SPP1** | Prognostic Biomarker | 0.82 | Cabozantinib (FDA, downstream) |
| 3 | **CDK1** | Mechanistic Hub | 0.68 | Dinaciclib (Phase II) |
| 4 | **BUB1B** | Emerging Target | 0.49 | BAY-1816032 (Preclinical) |

## Current Top Drug Candidates

| Rank | Drug | Target | Phase | CMap Score |
|------|------|------|------|------------|
| 1 | Volasertib | PLK1 | Phase III | −35 |
| 2 | Dinaciclib | CDK1 | Phase II | −24 |
| 3 | BI2536 | PLK1 | Phase II | −17 |
| 4 | Cabozantinib | SPP1 pathway | FDA Approved | −6 |

---

## Critical Scientific Rules

1. **Never assume results are correct.** Verify assumptions for every analysis.
2. **Report both supporting AND opposing evidence.** Do not hide contradictory findings.
3. **Create a "Potential Biases and Limitations" section for every major analysis.**
4. **Check biological plausibility** — would an HCC clinician agree with the finding?
5. **Check technical limitations** — what could the method miss?
6. **Separate target categories:** Therapeutic Target ≠ Biomarker ≠ Mechanistic Hub ≠ Drug Candidate.
7. **Treat CMap as hypothesis-generating only.** Never claim "Drug is effective" from CMap alone. Use "suggests"/"supports"/"indicates potential."
8. **Docking is not proof of efficacy.** Report RMSD, binding energy, redocking validation.
9. **Document failed analyses.** Never delete negative results.
10. **Maintain full reproducibility.** Every figure has a source script.

---

## Drug Discovery Standards

For each target, evaluate 8 dimensions:

1. Expression (DEG)
2. Survival (Cox HR)
3. PPI Centrality
4. Machine Learning Selection
5. External Validation (GEO)
6. Druggability (DrugBank/TTD/ChEMBL/DGIdb)
7. Literature Evidence (PubMed ID, DOI, year, study type)
8. Clinical Trial Evidence (Phase, NCT number)

Prioritize: clinical studies > meta-analyses > independent cohorts > review papers.

---

## Validation Requirements

Every finding must pass self-audit:

```
Finding: [statement]
Evidence For: [supporting points]
Evidence Against: [contradictory points]
Confidence: [High/Moderate/Low (score)]
Decision: [include/exclude/qualify]
```

---

## Git Workflow

After every completed step:

1. Save all scripts, figures, tables.
2. Update `PROJECT_STATUS.md`.
3. Create `logs/YYYY-MM-DD.md` research log.
4. Commit with format: `feat(stepXX): description`
5. Tag milestones: `git tag -a v1.0 -m "DEG complete"`

Commit prefixes: `feat:` `docs:` `fix:` `refactor:` `chore:`

---

## Research Log Rules

Every session creates `logs/YYYY-MM-DD.md` containing:
- Tasks Completed
- Key Findings
- Unexpected Findings
- Limitations
- Next Steps

---

## Target Figures (27 total)

| Figure | Content | Source |
|--------|---------|--------|
| Figure 1 | Analysis Workflow | — |
| Figure 2 | DEG Volcano + Heatmap | Step 03 |
| Figure 3 | GO + KEGG + GSEA | Step 04 |
| Figure 4 | PPI Network + Hub Genes | Step 05 |
| Figure 5 | LASSO + RF | Steps 06-07 |
| Figure 6 | Survival Analysis | Step 08 |
| Figure 7 | GEO External Validation | Steps 09-10,18 |
| Figure 8 | Drug Target Network | Step 11 |
| Figure 9 | Molecular Docking | Step 21 |
| Figure 10 | Drug Repurposing Pipeline | Steps 15-16 |

## Target Tables (6)

| Table | Content | Source |
|-------|---------|--------|
| Table 1 | Clinical Characteristics | Step 02 |
| Table 2 | Top 50 DEGs | Step 03 |
| Table 3 | Hub Genes + Topology | Step 05 |
| Table 4 | GEO Multi-Cohort Validation | Steps 09-10,18 |
| Table 5 | Candidate Drugs | Step 11 |
| Table 6 | Drug Repurposing Candidates | Step 15 |
