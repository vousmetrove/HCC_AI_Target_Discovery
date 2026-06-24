# Manuscript Tables

## Table 1: Clinical Characteristics of TCGA-LIHC Cohort

| Characteristic | Value |
|----------------|-------|
| Total samples | 421 |
| Tumor samples | 371 |
| Normal samples | 50 |
| Patients with complete OS data | 372 |
| Age (years), median (range) | 59 (16–90) |
| Gender (Male/Female) | 249/123 |
| AJCC Stage I/II/III/IV | 171/84/85/5 |
| Histologic Grade G1/G2/G3/G4 | 50/171/120/12 |
| OS events (deaths) | 130 (35.5%) |
| Median follow-up (days) | 631 |

**Source:** TCGA-LIHC via TCGAbiolinks; `data/processed/clinical_clean.rds`

---

## Table 2: Top 20 Differentially Expressed Genes in HCC

| Rank | Gene | log2FC | padj | Direction | Category |
|------|------|--------|------|-----------|----------|
| 1 | GABRD | +4.47 | 1.3e−112 | Up | CTA / GABA receptor |
| 2 | PLVAP | +2.91 | 1.5e−94 | Up | Endothelial fenestrae |
| 3 | CDKN3 | +3.90 | 5.0e−90 | Up | Cell cycle phosphatase |
| 4 | CDC25C | +4.34 | 1.8e−88 | Up | G2/M checkpoint |
| 5 | UBE2T | +3.19 | 2.3e−87 | Up | Fanconi anemia / p53 ubiquitination |
| 6 | NUF2 | +4.18 | 1.0e−86 | Up | Kinetochore component |
| 7 | CENPF | +3.87 | 1.4e−86 | Up | Centromere protein F |
| 8 | SKA1 | +4.46 | 3.6e−86 | Up | Spindle/kinetochore |
| 9 | ZIC2 | +6.66 | 1.3e−85 | Up | Zinc finger TF |
| 10 | TROAP | +4.15 | 4.5e−85 | Up | Cell adhesion |
| 11 | KIFC1 | +3.76 | 4.6e−83 | Up | Mitotic kinesin |
| 12 | KIF4A | +3.91 | 6.6e−82 | Up | Chromosome kinesin |
| 13 | BIRC5 | +4.08 | 8.5e−82 | Up | Survivin |
| 14 | HJURP | +3.88 | 1.3e−81 | Up | Chromosome segregation |
| 15 | GPC3 | +5.92 | 4.4e−81 | Up | HCC biomarker / drug target |
| 16 | CDK1 | +3.47 | 1.4e−73 | Up | Master mitotic kinase |
| 17 | PLK1 | +3.49 | 1.7e−66 | Up | Mitotic kinase |
| 18 | SPP1 | +4.32 | 2.9e−30 | Up | Osteopontin |
| 19 | BUB1B | +3.43 | 2.0e−56 | Up | SAC kinase |
| 20 | NR0B1 | +6.18 | 4.4e−16 | Up | CTA / nuclear receptor |

**Source:** `results/DEG/DEG_with_SYMBOL_biomaRt.csv`

---

## Table 3: PPI Hub Gene Topology Metrics

| Gene | Degree | Betweenness | Closeness | MCC | Methods (of 4) | Biological Role |
|------|--------|-------------|-----------|-----|----------------|-----------------|
| CDK1 | **59** | 0.2172 | **0.3458** | ✓ | **4/4** | Master mitotic kinase |
| H2AX | 44 | 0.0472 | 0.3120 | ✓ | 3/4 | DNA damage marker |
| PLK1 | 47 | 0.0582 | 0.3151 | — | 3/4 | Mitotic kinase |
| PCNA | 43 | 0.0858 | 0.3094 | — | 3/4 | Proliferation marker |
| BUB1B | 38 | 0.0273 | 0.3058 | ✓ | 3/4 | SAC kinase |

**Network:** 1,125 nodes, 2,906 edges, giant component = 600 nodes
**Source:** `results/PPI/ppi_node_metrics.rds`

---

## Table 4: Univariate Cox Regression — Candidate Genes

| Gene | HR | 95% CI | P-value | Category | External GEO |
|------|-----|--------|---------|----------|--------------|
| NR0B1 | 1.247 | 1.141–1.363 | 1.2e−6 | CTA/Developmental | Direction OK, NS |
| GAGE2A | 1.358 | 1.198–1.539 | 1.7e−6 | CTA | Opposite direction (removed) |
| **PLK1** | **1.372** | **1.199–1.570** | **4.2e−6** | **Therapeutic Target** | **2/2 validated** |
| SPP1 | 1.124 | 1.067–1.183 | 9.8e−6 | Prognostic Biomarker | 2/2 validated |
| TRIM54 | 1.205 | 1.099–1.321 | 7.6e−5 | RF+LASSO consensus | Direction OK, NS |
| CDK1 | 1.310 | 1.146–1.499 | 8.1e−5 | Mechanistic Hub | No probe |
| BUB1B | 1.291 | 1.125–1.480 | 2.6e−4 | Emerging CIN Target | 1/2 validated |
| H2AX | 1.300 | 1.094–1.545 | 2.8e−3 | PPI Hub #2 | 1/2 validated |
| LY6H | 1.186 | 1.043–1.349 | 9.5e−3 | RF+LASSO consensus | Direction OK, NS |
| PCNA | 1.391 | 1.093–1.769 | 7.3e−3 | PPI Hub #4 | 1/2 validated |
| GLP1R | 1.126 | 1.027–1.235 | 1.1e−2 | RF+LASSO consensus | Direction OK, NS |
| DRGX | 1.156 | 1.014–1.318 | 3.0e−2 | RF+LASSO consensus | Direction inconsistent |

**Source:** `results/Survival/cox_univariate.csv`

---

## Table 5: Candidate Therapeutic Agents for HCC

| Drug | Target | Mechanism | Phase | FDA | HCC Evidence | NCT |
|------|--------|-----------|-------|-----|--------------|-----|
| **Volasertib** | **PLK1** | PLK1 ATP-competitive inhibitor | Phase III | No | HCC Phase II completed | NCT02139267 |
| BI2536 | PLK1 | PLK1 inhibitor | Phase II | No | Preclinical HCC; synergy with sorafenib | NCT00701766 |
| Onvansertib | PLK1 | Oral PLK1 inhibitor | Phase II | No | HCC inferred from target | — |
| Dinaciclib | CDK1 | CDK1/2/5/9 pan-inhibitor | Phase II | No | HCC cell lines | NCT01676753 |
| Milciclib | CDK1 | CDK1/2 inhibitor | Phase II | No | HCC Phase II completed | NCT03109886 |
| Flavopiridol | CDK1 | CDK1/2/4/7 pan-inhibitor | Phase II | No | HCC Phase II completed (modest) | NCT00006379 |
| **Cabozantinib** | **SPP1 pathway** | VEGFR2/MET inhibitor | **FDA Approved** | **Yes** | **HCC second-line approved** | NCT01908426 |
| Huzhangoside D | SPP1 | SPP1 signaling modulator | Preclinical | No | HCC migration/invasion | — |

**Source:** `results/Drug/drug_target_table.csv`

---

## Table 6: Top CMap Drug Repurposing Candidates

| Rank | Drug | CMap Score | Phase | Prolif↓ | Apopt↑ | Target | FDA | HCC Relevance |
|------|------|------------|-------|---------|--------|-------|-----|---------------|
| 1 | Volasertib | **−35** | Phase III | 15/25 | 5/15 | PLK1 | No | HCC Phase II completed |
| 2 | Dinaciclib | −24 | Phase II | 10/25 | 4/15 | CDK1/2/5/9 | No | HCC cell lines |
| 3 | BI2536 | −17 | Phase II | 7/25 | 3/15 | PLK1 | No | Preclinical HCC |
| 4 | Doxorubicin | −15 | Approved | 5/25 | 5/15 | TOP2A | **Yes** | TACE first-line |
| 5 | Nutlin3 | −14 | Preclinical | 4/25 | 6/15 | MDM2 | No | p53 activator |
| 6 | Alisertib | −12 | Phase II | 5/25 | 2/15 | AURKA | No | Mitotic inhibitor |
| 7 | Flavopiridol | −12 | Phase II | 5/25 | 2/15 | CDK1/2/4/7 | No | HCC Phase II |
| 8 | Sorafenib | −11 | **Approved** | 4/25 | 3/15 | VEGFR2/RAF | **Yes** | HCC first-line |
| 9 | Metformin | −11 | **Approved** | 4/25 | 3/15 | AMPK | **Yes** | HCC epidemiology |
| 10 | Palbociclib | −10 | **Approved** | 4/25 | 2/15 | CDK4/6 | **Yes** | Breast cancer → HCC repositioning |

**Scoring:** Negative CMap score = drug reverses HCC proliferative signature (therapeutic).
Prolif↓ = HCC proliferation genes downregulated by drug. Apopt↑ = apoptosis genes upregulated by drug.
**Source:** `results/Drug/cmap_connectivity_scores.csv`

---

## Supplementary Tables

### Table S1: Full DEG List
File: `results/DEG/DEG_with_SYMBOL_biomaRt.csv` (6,301 genes × 10 columns)

### Table S2: 7-Dimension Evidence Matrix
| Gene | DEG | PPI | ML | Cox | GEO | Drug | CMap | Evidence Level | Final Classification |
|------|-----|-----|-----|-----|-----|------|------|----------------|---------------------|
| PLK1 | ✓ | ✓ (Hub #3) | — | ✓ (HR=1.37) | ✓ (2/2) | ✓ (0.960) | ✓ (−35) | 9/9 | Primary Therapeutic Target |
| SPP1 | ✓ | — | ✓ (LASSO+RF) | ✓ (HR=1.12) | ✓ (2/2) | ✓ (0.595) | — | 8/9 | Primary Prognostic Biomarker |
| CDK1 | ✓ | ✓ (Hub #1) | — | ✓ (HR=1.31) | ✗ (no probe) | ✓ (0.975) | ✓ (−24) | 7/9 | Mechanistic Hub |
| BUB1B | ✓ | ✓ (Hub #5) | — | ✓ (HR=1.29) | ✓ (1/2) | ✓ (0.260) | — | 6/9 | Emerging CIN Target |
| NR0B1 | ✓ | — | ✓ (LASSO+RF) | ✓ (HR=1.25) | ~ (NS) | ✗ (orphan) | — | 5/9 | Epigenetic Biomarker |
| GAGE2A | ✓ | — | ✓ (LASSO+RF) | ✓ (HR=1.36) | ✗ (opposite) | — | — | 4/9 | REMOVED |

**Source:** `results/Evidence/gene_evidence_matrix.csv`
