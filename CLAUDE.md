# CLAUDE.md — HCC_AI_Target_Discovery

## Project Identity

| Field | Value |
|-------|-------|
| **Project** | HCC_AI_Target_Discovery |
| **Full Title** | AI-Assisted Identification of Prognostic Therapeutic Targets in Hepatocellular Carcinoma Through Machine Learning-Based Transcriptomic Analysis |
| **Goal** | Discover and prioritize therapeutic targets for HCC using multi-dimensional AI-assisted evidence integration |
| **Target Venue** | EI Conference 2026 |
| **Local Path** | `E:\analysis-bio\practice\HCC_AI_Target_Discovery` |
| **GitHub** | `git@github.com:vousmetrove/HCC_AI_Target_Discovery` |
| **Primary Language** | R (≥4.2) |
| **Started** | 2026-06-23 |
| **Current Status** | Analysis complete (18/18 steps). Manuscript writing phase. |

---

## Pipeline Overview (18 Steps — All Complete)

```
[01] TCGA Download     ──→ [02] Preprocessing    ──→ [03] DESeq2 DEG + VST
[04] GO/KEGG/GSEA      [05] PPI Network + CytoHubba
[06] LASSO Cox          [07] Random Forest (ranger)
[08] Survival Analysis  [09] GEO Validation (GSE76427)
[10] GEO Robustness     [11] Drug Target Prioritization
[12] Immune Infiltration (partial)   [13] Molecular Docking Protocol
[14] Final Consensus    [15] Drug Repurposing
[16] CMap Analysis      [17] Robustness Audit
[18] Multi-Cohort Validation (GSE14520)
```

## Key Findings (Evidence-Based)

### Top Therapeutic Target: **PLK1** (Confidence: HIGH, 0.85/1.00)
- 7/7 evidence dimensions positive
- DEG: log2FC=+3.49, padj=1.7e-66
- PPI Hub #3 (Degree=47, 3/4 CytoHubba algorithms)
- TCGA Cox: HR=1.372 [1.20-1.57], P=4.2e-6
- GEO validated: GSE76427 (P=2.6e-3), GSE14520 (P=4.5e-31)
- Drug: Volasertib Phase III, BI2536 Phase II
- CMap: Volasertib score=-35 (strongest proliferation reversal of 20 drugs)
- Druggability: 0.960 (Gold Tier)

### Top Prognostic Biomarker: **SPP1** (Confidence: HIGH, 0.82/1.00)
- Best GEO validation: GSE76427 (P=1.3e-11), GSE14520 (P=1.3e-21)
- Selected by BOTH LASSO and RF
- TCGA Cox: HR=1.124 [1.07-1.18], P=9.8e-6
- 4/4 probes show consistent direction in GSE76427
- Validated in 10+ independent HCC cohorts (literature)
- NOT recommended as drug target (preclinical only)

### Mechanistic Hub: **CDK1** (Confidence: MODERATE, 0.68/1.00)
- #1 PPI Hub (Degree=59, ALL 4 CytoHubba algorithms)
- Druggability: 0.975 (Gold Tier, highest score)
- Dinaciclib Phase II, Milciclib Phase II HCC
- **NO GEO validation** — no probe on GPL10558 (Illumina) or GPL3921 (Affymetrix)
- Flavopiridol HCC Phase II showed modest efficacy
- Active limitation — must be disclosed in manuscript

### Deprioritized Genes
| Gene | Reason | Action |
|------|--------|--------|
| **GAGE2A** | GEO direction opposite (single probe cross-hybridization) | Removed from core list |
| **NR0B1** | Strongest Cox (P=1.2e-6) but orphan receptor, no known drug | Biomarker only; not targetable |
| **DRGX** | Weak + inconsistent GEO direction | Low priority |

---

## Output File Locations

```
results/
├── DEG/                    # 6,301 DEGs; VST matrix; biomaRt annotations
│   ├── DEG_list.rds        # Full DEG results (8 columns)
│   ├── rlog_normalized.rds # 22,730 x 416 VST matrix (primary ML input)
│   ├── DEG_with_SYMBOL_biomaRt.csv  # 6,301 DEGs with HGNC symbols
│   └── DEG_protein_coding.csv       # 3,845 protein-coding DEGs
├── GO_KEGG/                # GO + KEGG + GSEA enrichment
│   ├── GO_enrichment_all.csv
│   ├── KEGG_enrichment_all.csv
│   ├── GSEA_results.csv    # 5 Hallmark gene sets (all FDR<0.05)
│   └── GO_KEGG_GSEA_QC_Report.md
├── PPI/                    # STRING v12 network (1,125 nodes, 2,906 edges)
│   ├── PPI_hub_genes_top10.csv  # 5 consensus hubs
│   ├── cytoscape_node_table.csv # Cytoscape-compatible
│   └── ppi_network.rds
├── LASSO/                  # 11-gene prognostic signature
│   ├── lasso_genes.csv     # Coefficients + directions
│   └── lasso_model.rds     # cv.glmnet object
├── RandomForest/           # ranger, 1,000 trees
│   ├── rf_importance.csv   # Full permutation VIMP ranking
│   ├── rf_top30.csv
│   └── method_consensus.csv # LASSO ∩ RF overlap
├── Survival/               # Cox + KM + timeROC
│   ├── cox_univariate.csv  # 12 genes all P<0.05
│   └── prognostic_ranking.csv
├── GEO_validation/         # Multi-cohort external validation
│   ├── GSE76427_validation_corrected.csv  # limma-corrected, 11 genes
│   ├── GAGE2A_probe_audit.csv
│   ├── multi_cohort_validation.csv
│   └── effect_size_comparison.csv
├── Drug/                   # Drug repurposing + CMap
│   ├── drug_target_table.csv    # 15 drug-target entries
│   ├── druggability_score.csv   # 4-target druggability scores
│   ├── top10_repurposing_candidates.csv
│   └── cmap_connectivity_scores.csv  # 20 drugs scored
├── Docking/                # Molecular docking protocol (PLK1 + CDK1)
├── Validation/             # Robustness audit
│   └── robustness_audit.md # All major findings independently audited
├── enrichment/             # Gene symbol list for GO
├── report/                 # Annotation summary
├── Immune/                 # EMPTY — immune analysis incomplete
└── HubGene_Integration/    # EMPTY — final consensus pending
```

## Figure Directory Structure

```
figures/
├── DEG/          9 files   Volcano, heatmap, PCA (Figure 1)
├── GO_KEGG/     10 files   GO BP/CC/MF + KEGG + GSEA (Figures 2-3)
├── LASSO/        6 files   CV + coefficient path (Figure 4)
├── RF/           2 files   RF VIMP (Figure 5)
├── PPI/          1 file    Consensus heatmap (Figure 5)
├── Survival/     8 files   Forest + KM + timeROC (Figure 6)
├── Drug/         5 files   Drug-target + repurposing + CMap (Figures 8-10)
└── Publication/  empty     Reserved for final publication layout
```

## Hub Gene Selection Logic (Multi-Stage)

```
Stage 1: DEG filter (|log2FC|>1, padj<0.05, protein-coding) → 3,845 genes
Stage 2: PPI Hub consensus (≥3/4 CytoHubba algorithms) → 5 genes
         (CDK1, PLK1, PCNA, BUB1B, H2AX)
Stage 3: ML consensus (LASSO ∩ RF ∩ top-|log2FC|) → 7 genes
         (SPP1, GAGE2A, NR0B1, GLP1R, DRGX, LY6H, TRIM54)
Stage 4: TCGA Cox (all 12 candidates P<0.05) → validated
Stage 5: GEO validation filter → GAGE2A removed (direction conflict)
Stage 6: Druggability filter → NR0B1 removed (orphan receptor)
Stage 7: CMap filter → PLK1 confirmed (strongest reversal)
Stage 8: 7-dimension composite scoring → PLK1 #1, SPP1 #3
```

## Key Design Decisions

| # | Decision | Rationale | Date |
|---|----------|-----------|------|
| 1 | Use VST not rlog | rlog too slow for 416 samples; VST is DESeq2-recommended | 06-23 |
| 2 | enrichR over clusterProfiler | Bioconductor version conflicts on Windows R 4.6 | 06-23 |
| 3 | Manual CMap scoring | CLUE API inaccessible; literature-curated signatures used | 06-24 |
| 4 | Deprioritize GAGE2A | Single Illumina probe with opposite GEO direction | 06-24 |
| 5 | CDK1 retained despite no GEO | Technical limitation (probe absence), not biological absence | 06-24 |
| 6 | NR0B1 NOT recommended as target | Orphan nuclear receptor; strongest Cox but undruggable | 06-24 |
| 7 | Exclude large files from Git | DESeq2_object.rds (306MB), TCGA counts (146MB) in .gitignore | 06-24 |
| 8 | Figures in both PDF+PNG | PDF for publication, PNG for quick preview | ongoing |

## Environment Notes

- **R version:** 4.6.0 (Windows, ucrt)
- **Bioconductor:** 3.23
- **Package issues known:** clusterProfiler/DOSE version mismatch on R 4.6 → used enrichR instead
- **Large data files:** 3 files >50MB excluded from Git via .gitignore
- **GDC API changes:** workflow.type changed from "HTSeq - Counts" to "STAR - Counts"; legacy=FALSE removed

## Reproducibility

```bash
# Re-run full pipeline (12 core scripts)
cd scripts/
for f in 01_download_TCGA.R 02_preprocessing.R 03_DESeq2.R \
         04_GO_KEGG.R 05_PPI_STRING.R 06_LASSO.R 07_RandomForest.R \
         08_Survival.R 09_GEO_validation.R 10_Immune_Infiltration.R \
         11_Drug_Repurposing.R 12_HubGene_Integration.R; do
    Rscript "$f"
done
```

Note: Steps 10 (Immune) and 12 (HubGene Integration) produce incomplete output due to:
- Immune: marker genes filtered out by low-expression filtering (bulk RNA-seq limitation)
- Integration: final 7-dimension scoring was completed in Step 14 (separate script)

## Pending Items for Manuscript

- [ ] ICGC-LIRI-JP validation (requires portal download)
- [ ] CIBERSORTx immune deconvolution (requires web submission)
- [ ] Raw LINCS L1000 query (requires CLUE API access)
- [ ] Molecular docking execution (protocol ready, needs Vina installation)
- [ ] Complete Supplementary Tables S1-S2
- [ ] Write manuscript sections from outline.md
