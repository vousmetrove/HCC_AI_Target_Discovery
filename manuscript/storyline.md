# Manuscript Storyline — AI-Assisted Therapeutic Target Discovery for HCC

## Core Narrative: Discovery → Validation → Prioritization

---

## Stage 1: Discovery (What did we find?)

### 1.1 The HCC Transcriptomic Landscape
- TCGA-LIHC DESeq2 analysis: 6,301 DEGs (4,844 upregulated, 1,457 downregulated)
- Proliferative HCC subtype dominates: E2F Targets (NES=+4.32), G2M Checkpoint (NES=+4.36)
- CTA genes (GABRD, MAGEA family) show extreme upregulation due to genome-wide hypomethylation
- **Figure 1:** Volcano plot + heatmap + PCA

### 1.2 Functional Characterization
- Cell cycle (#1 KEGG, P=2.3e-5), DNA replication (#2, P=1.8e-4)
- Complement cascade downregulation (#1 down, P=2.7e-16) — loss of hepatocyte identity
- Synaptic/neural GO terms identified as CTA artifacts, not neuroendocrine differentiation
- **Figures 2–3:** GO + KEGG + GSEA

### 1.3 Network-Level Architecture
- PPI network: 1,125 nodes, 2,906 edges
- 5 hub genes identified via CytoHubba consensus (CDK1 unique 4/4 algorithm hub)
- SAC (Spindle Assembly Checkpoint) subnetwork: BUB1B → CDC20 → BUB1 → MAD2L1
- **Figure 4:** PPI network visualization

---

## Stage 2: Validation (Is it real?)

### 2.1 Machine Learning Confirmation
- LASSO Cox: 11-gene prognostic signature (C-index=0.654)
- Random Forest: Top 30 permutation VIMP ranking (C-index=0.695)
- RF ∩ LASSO consensus: 7 genes — independent ML agreement
- **Figure 5:** LASSO + RF results

### 2.2 Clinical Survival Validation
- All 12 candidate genes Cox-significant (P<0.05)
- PLK1: HR=1.372 [1.20–1.57], P=4.2e-6 — highest HR among PPI hub genes
- SPP1: HR=1.124 [1.07–1.18], P=9.8e-6 — selected by both ML methods
- **Figure 6:** Forest plot + KM curves + timeROC

### 2.3 External Cohort Validation
- GSE76427 (Illumina): PLK1 validated (P=2.6e-3), SPP1 validated (P=1.3e-11)
- GSE14520 (Affymetrix): PLK1 validated (P=4.5e-31), SPP1 validated (P=1.3e-21)
- CDK1: No probe on either platform — technical limitation, not biological
- GAGE2A: Single probe, direction opposite to TCGA — REMOVED from core list
- **Figures 7, 11:** GEO validation forest plot + KM curves

### 2.4 Robustness Confirmation
- Negative control: 0/1000 random genes score ≥ PLK1 (P=0.001, Z=15.4)
- PLK1 stable across all 6 weight schemes
- Independent audit confirms PLK1 confidence 0.85/1.00
- **Figure:** Negative control distribution plot

---

## Stage 3: Prioritization (Which target to pursue?)

### 3.1 Drug Target Evaluation
- PLK1: Druggability 0.960 (Gold); Volasertib Phase III, BI2536 Phase II, Onvansertib Phase II
- CDK1: Druggability 0.975 (Gold); Dinaciclib Phase II, Milciclib Phase II HCC
- SPP1: Druggability 0.595 (Silver); preclinical only — biomarker, not drug target
- **Figure 8:** Drug-target network

### 3.2 Pharmacological Validation (CMap)
- Volasertib: CMap score −35 (strongest of 20 drugs); reverses 15/25 proliferation genes
- Dinaciclib: CMap score −24; reverses 10/25 proliferation genes
- Known HCC drugs (Sorafenib, Doxorubicin) correctly identified — validates methodology
- **Figures 9–10:** Drug repurposing + CMap pipeline

### 3.3 Final 7-Dimension Evidence Integration
- PLK1: 0.845 (Tier 1) — Primary Therapeutic Target
- CDK1: 0.735 (Tier 1) — Mechanistic Hub (GEO limitation acknowledged)
- SPP1: 0.675 (Tier 2) — Primary Prognostic Biomarker
- BUB1B: 0.485 (Tier 2) — Emerging CIN Target

---

## Stage 4: Drug Repurposing (What can we reposition?)

### 4.1 Primary Repositioning: PLK1 → Volasertib
- Convergent evidence: DEG + PPI + Cox + GEO + Drug + CMap + Negative Control
- Volasertib Phase III (completed); HCC Phase II (completed)
- CMap independently confirms PLK1 inhibition reverses HCC program
- Limitation: Not in GDSC/CTRP (too recent for screening panels)

### 4.2 Secondary Repositioning: CDK1 → Dinaciclib
- #1 PPI Hub + Gold druggability
- Dinaciclib Phase II; Milciclib HCC Phase II completed (NCT03109886)
- Limitation: No GEO validation; Flavopiridol modest efficacy precedent

### 4.3 Biomarker-Guided Repositioning: SPP1 → Cabozantinib
- Cabozantinib is FDA-approved for HCC (second-line)
- Targets VEGFR2/MET — downstream of SPP1 signaling pathway
- SPP1 may serve as companion biomarker for Cabozantinib stratification

---

## The Paper's Central Argument

> *"Using a 9-level AI-assisted evidence framework integrating transcriptomics, network biology, machine learning, multi-cohort validation, and drug perturbation analysis, we identify **PLK1 as the highest-confidence therapeutic target** and **SPP1 as a robust prognostic biomarker** for proliferative hepatocellular carcinoma. The PLK1 inhibitor Volasertib (Phase III) emerges as a strong drug repositioning candidate, independently supported by a Connectivity Map score of −35. This framework demonstrates that systematic multi-dimensional evidence integration can prioritize therapeutic targets with higher confidence than single-dimension approaches."*
