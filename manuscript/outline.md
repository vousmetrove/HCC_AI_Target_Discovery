# Manuscript Outline — EI Conference 2026

**Title (proposed):**
AI-Assisted Identification and Multi-Dimensional Prioritization of Therapeutic Targets in Hepatocellular Carcinoma Through Transcriptomic Network Analysis

**Authors:** Vousmetrove et al.
**Target Journal:** EI Conference Proceedings, 2026

---

## Abstract (target: 250 words)

- **Background:** HCC remains a leading cause of cancer mortality with limited systemic therapeutic options. Systematic multi-evidence target discovery is needed.
- **Methods:** TCGA-LIHC (371T/50N) DESeq2 → 6,301 DEGs. PPI network (STRING v12) + CytoHubba. LASSO + Random Forest. Connectivity Map drug repurposing. External validation in GSE76427 + GSE14520.
- **Results:** PLK1 emerged as top therapeutic target (7/7 evidence, CMap score −35, Volasertib Phase III). SPP1 validated as robust prognostic biomarker (GEO P=1.3e−11, HR=1.12). CDK1 identified as #1 PPI hub (degree=59) but lacked GEO probe availability.
- **Conclusion:** Multi-dimensional AI-assisted framework prioritizes PLK1 (therapeutic) and SPP1 (prognostic) for HCC.

**Keywords:** hepatocellular carcinoma; therapeutic target discovery; PLK1; SPP1; connectivity map; PPI network; external validation

---

## 1. Introduction

### 1.1 HCC Clinical Landscape
- Third leading cause of cancer mortality globally
- Current approved therapies: Sorafenib, Lenvatinib, Regorafenib, Cabozantinib, Atezolizumab+Bevacizumab
- Unmet need: drug resistance, limited response rates
- Rationale for computational target discovery

### 1.2 Transcriptomic Approaches to HCC
- TCGA-LIHC landmark studies (Wheeler et al. Cancer Cell 2017; TCGA Cell 2017)
- Hoshida molecular classification (S1/S2/S3)
- Gap: translation from DEG lists to actionable therapeutic targets

### 1.3 AI-Assisted Target Discovery Framework
- Multi-dimensional evidence integration
- 7-dimension framework: Expression + Network + ML + Survival + GEO + Drug + Immune

### 1.4 Study Objectives
1. Identify robust DEG signature of proliferative HCC
2. Construct PPI network for hub gene identification
3. Apply ML for prognostic feature selection
4. Validate in independent GEO cohorts
5. Prioritize therapeutic targets via drug repurposing + CMap
6. Establish multi-dimensional evidence scoring framework

---

## 2. Methods

### 2.1 Data Acquisition and Preprocessing (Steps 01–02)
- TCGA-LIHC: TCGAbiolinks → STAR counts (60,660 genes x 421 samples)
- Low-expression filtering: edgeR::filterByExpr (min.count=10, min.total.count=15)
- 22,730 genes retained; 371 tumors, 50 normal; 372 patients with complete OS data

### 2.2 Differential Expression Analysis (Step 03)
- DESeq2: negative binomial GLM; Wald test; |log2FC|>1, padj<0.05
- 6,301 DEGs (4,844 upregulated; 1,457 downregulated)
- VST normalization; biomaRt annotation: 5,212 unique HGNC symbols

### 2.3 GO/KEGG/GSEA Enrichment (Step 04)
- enrichR: GO_BP/CC/MF_2023, KEGG_2021_Human
- fgsea: MSigDB Hallmark gene sets
- Up/downregulated DEGs analyzed separately

### 2.4 PPI Network and Hub Gene Identification (Step 05)
- STRING v12: physical network, combined_score >= 0.4
- igraph: 1,125 nodes, 2,906 edges, 600-node giant component
- CytoHubba: MCC, Degree, Closeness, Betweenness consensus

### 2.5 Machine Learning Feature Selection (Steps 06–07)
- LASSO Cox: glmnet, 10-fold CV, 300 top-|log2FC| protein-coding DEGs
- Random Forest: ranger, 1,000 trees, permutation VIMP
- Cross-method consensus: RF intersect LASSO

### 2.6 Survival Analysis (Step 08)
- Cox univariate regression (HR, 95% CI); Kaplan-Meier (median split)
- Time-dependent ROC (1/3/5-year AUC); LASSO risk score analysis

### 2.7 GEO External Validation (Steps 09–10, 18)
- GSE76427: Illumina HumanHT-12 v4, 115T/52N, limma + neqc
- GSE14520: Affymetrix U133A, 225T/220N, limma
- Probe audit for GAGE2A direction discrepancy
- Multi-cohort validation forest plot

### 2.8 Drug Repurposing and CMap (Steps 11, 15–16)
- DrugBank, TTD, OpenTargets, ChEMBL, DGIdb
- CMap: 150 UP + 150 DOWN disease signature vs drug perturbation
- Mechanism-based connectivity scoring (proliferation reversal)

### 2.9 Robustness Audit (Step 17)
- Independent validation audit: supporting + contradictory evidence
- Confidence scoring: High / Moderate / Low for each finding

---

## 3. Results

### 3.1 DEG Landscape (Figure 1)
- 6,301 DEGs; top 50 all protein-coding
- GSEA: E2F Targets (NES=+4.32), G2M Checkpoint (NES=+4.36) — FDR<10^-50
- Proliferative HCC subtype dominance confirmed

### 3.2 Functional Enrichment (Figures 2–3)
- Cell cycle (#1 KEGG, P=2.3e-5); DNA replication (#2, P=1.8e-4)
- Complement cascade (#1 downregulated, P=2.7e-16)
- CTA-driven synaptic terms identified and contextualized

### 3.3 PPI Network (Table 3)
- 5 consensus hub genes: CDK1, H2AX, PLK1, PCNA, BUB1B
- CDK1: only gene in Top 10 of ALL 4 CytoHubba algorithms (degree=59)

### 3.4 ML Prognostic Signature (Figures 4–5)
- LASSO: 11-gene signature (C-index: Train=0.75, Test=0.65)
- RF: Top 30 permutation importance ranking
- RF + LASSO consensus: 7 genes

### 3.5 Survival Validation (Figure 6)
- All 12 candidate genes significant in Cox (P<0.05)
- Top prognostic: NR0B1 (HR=1.25, P=1.2e-6), PLK1 (HR=1.37, P=4.2e-6)
- LASSO risk score C-index = 0.72

### 3.6 GEO External Validation (Figure 11)
- PLK1: validated in 2/2 external cohorts (GSE76427 P=2.6e-3; GSE14520 P=4.5e-31)
- SPP1: validated in 2/2 external cohorts (GSE76427 P=1.3e-11; GSE14520 P=1.3e-21)
- CDK1: no probe on Illumina or Affymetrix platforms (technical limitation)
- Direction consistency: 4/4 (100%, excluding platform failures)

### 3.7 Drug Target Prioritization (Figure 8, Table 5)
- PLK1: Druggability 0.960 (Gold); Volasertib Phase III, BI2536 Phase II
- CDK1: Druggability 0.975 (Gold); Dinaciclib Phase II, Milciclib Phase II HCC
- SPP1: Druggability 0.595 (Silver); Cabozantinib FDA-approved (downstream)

### 3.8 CMap Drug Repurposing (Figures 9–10, Table 6)
- Volasertib: CMap score -35 (strongest); reverses 15/25 HCC proliferation genes
- Dinaciclib: CMap score -24; Doxorubicin: -15; Sorafenib: -11
- Known HCC drugs correctly identified — validates methodology

### 3.9 Final Consensus Framework
- PLK1: 0.845 (Tier 1) — primary therapeutic target
- CDK1: 0.735 (Tier 1) — core mechanistic hub (GEO limitation)
- SPP1: 0.675 (Tier 2) — primary prognostic biomarker
- BUB1B: 0.485 (Tier 2) — emerging CIN target

---

## 4. Discussion

### 4.1 PLK1 as a High-Confidence Therapeutic Target
- Most complete evidence package; 2/2 GEO cohorts validated
- CMap independently confirms PLK1 inhibition reverses HCC transcriptional program
- Volasertib: Phase III with completed HCC Phase II

### 4.2 SPP1 as a Robust Prognostic Biomarker
- Strongest GEO validation; consistent across 4 probes
- Bridges CTA/ECM/immune signaling
- Cabozantinib targets SPP1 downstream pathway

### 4.3 CDK1 — Central Hub with Validation Gap
- #1 PPI hub (all 4 algorithms)
- No GEO validation due to platform probe absence — technical, not biological
- Flavopiridol modest Phase II efficacy; bone marrow toxicity

### 4.4 CTA De-repression as HCC Epigenetic Phenomenon
- GABRD, GAGE2A, MAGEA family: extreme log2FC from promoter hypomethylation
- Synaptic/neural GO terms are CTA artifacts
- GAGE2A deprioritized after GEO direction discrepancy

### 4.5 CMap: Independent Pharmacological Validation
- Known HCC drugs correctly identified — validates methodology
- PLK1 inhibitors show strongest proliferation reversal

### 4.6 Limitations
1. CDK1 no GEO validation (explicitly stated)
2. Immune analysis incomplete (bulk RNA-seq marker limitation)
3. CMap uses curated signatures, not raw LINCS L1000
4. ML models biased toward CTA genes (high |log2FC|)
5. ICGC-LIRI-JP validation pending
6. Docking protocol defined but not executed

### 4.7 Future Directions
- HCC-specific Phase III for Volasertib with SPP1 stratification
- CIBERSORTx deconvolution for immune validation
- Raw LINCS L1000 query; molecular docking + MD simulation

---

## 5. Conclusion

1. A 7-dimension AI-assisted framework for therapeutic target prioritization in HCC is presented.
2. **PLK1** is the highest-confidence therapeutic target (all dimensions positive; Phase III drug).
3. **SPP1** is a robust prognostic biomarker with patient stratification potential.
4. **CDK1** is the central PPI hub but requires RNA-seq-based external validation.
5. CMap independently converges on PLK1 inhibition for reversing proliferative HCC.
6. This framework is generalizable to other cancers.

---

## Figures Index

| Figure | Content | Source Step |
|--------|---------|-------------|
| Figure 1 | Volcano + Heatmap + PCA | 03 |
| Figure 2 | GO BP/CC/MF + KEGG dot plots | 04 |
| Figure 3 | GSEA Hallmark bar plot | 04 |
| Figure 4 | LASSO CV + coefficient path | 06 |
| Figure 5 | RF VIMP + consensus heatmap | 07 |
| Figure 6 | Forest plot + KM + timeROC | 08 |
| Figure 7 | GEO KM curves | 09-10 |
| Figure 8 | Drug-Target Network | 11 |
| Figure 9 | Drug Repurposing Network | 15 |
| Figure 10 | CMap Drug Discovery Pipeline | 16 |
| Figure 11 | Multi-Cohort Validation Forest | 18 |

## Tables Index

| Table | Content | Source Step |
|-------|---------|-------------|
| Table 1 | DEG Summary Statistics | 03 |
| Table 2 | GO/KEGG/GSEA Enrichment | 04 |
| Table 3 | PPI Hub Gene Metrics | 05 |
| Table 4 | Cox Univariate + Prognostic Ranking | 08 |
| Table 5 | Candidate Therapeutic Agents | 11 |
| Table 6 | Top CMap Candidate Compounds | 16 |
| Table S1 | Full DEG List (6,301 genes) | Supplementary |
| Table S2 | 7-Dimension Evidence Matrix | 14 |

---

## Key References

1. Wheeler DA, et al. *Cell* 169(7):1327-1341, 2017.
2. Hoshida Y, et al. *Cancer Res* 69(18):7385-7392, 2009.
3. TCGA Network. *Cell* 169(7):1327-1341, 2017.
4. Zucman-Rossi J, et al. *Gastroenterology* 149(5):1226-1239, 2015.
5. Guichard C, et al. *Nat Genet* 44(6):694-698, 2012.
6. Cao K, et al. *Oncogene* 43:35-46, 2024.
7. Huang C, et al. *Cell Death Dis* 2026.
8. Liu LP, et al. *BBRC* 493(1):20-27, 2017.
9. Lamb J, et al. *Science* 313(5795):1929-1935, 2006.

---

## Draft Writing Status

| Section | Data Ready | Notes |
|---------|-----------|-------|
| Abstract | Yes | All findings complete |
| Introduction | Yes | Literature review needed |
| Methods | Yes | All steps documented |
| Results 3.1-3.3 | Yes | DEG + Enrichment + PPI |
| Results 3.4-3.6 | Yes | ML + Survival + GEO |
| Results 3.7-3.9 | Yes | Drug + CMap + Consensus |
| Discussion | Yes | Key targets + limitations |
| Conclusion | Yes | Final rankings ready |
| References | Partial | 9 key papers listed |
