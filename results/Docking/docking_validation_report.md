# Molecular Docking Reliability Assessment — HCC Therapeutic Targets

**Date:** 2026-06-24 | **Status:** Protocol & validation framework | **Execution:** Pending Vina installation

---

## 1. PDB Structure Verification (Corrected)

### PLK1 (Primary Target)

| PDB ID | Resolution | Co-crystal Ligand | Ligand Type | Fitness for Re-docking |
|--------|-----------|-------------------|-------------|----------------------|
| **2RKU** | **1.95 Å** | **BI-2536 (R78)** | Clinical PLK1 inhibitor (Phase II) | **BEST** — Highest resolution + clinically relevant ligand |
| 2OWB | 2.10 Å | PHA-680626 (626) | Experimental pyrrolo-pyrazole inhibitor | Good — alternative binding pose |
| 3FVH | 2.10 Å | Rigosertib (ON-01910) | PLK1/PI3K dual inhibitor (Phase III) | Alternative — different chemotype |

**Correction:** Step 13 protocol listed 2OWB as having BI-2536. **2OWB actually contains PHA-680626.** The correct structure for BI-2536 re-docking is **2RKU (1.95 Å)**.

### CDK1 (Secondary Target)

| PDB ID | Resolution | Co-crystal Ligand | Ligand Type | Fitness |
|--------|-----------|-------------------|-------------|---------|
| **4Y72** | 2.30 Å | LZ9 | ATP-competitive inhibitor | **Moderate** — resolution adequate but not optimal |
| **6GU2** | **1.65 Å** | RO-3306 | Selective CDK1 inhibitor | **Better** — higher resolution + selective inhibitor |

**Correction:** Step 13 protocol listed 4Y72 as having Dinaciclib. **4Y72 actually contains LZ9.** Dinaciclib co-crystal structures exist but were not identified in our initial search. **6GU2 (1.65 Å, RO-3306) is recommended as the primary CDK1 re-docking structure.**

---

## 2. Re-docking Validation Protocol

### 2.1 Redocking of Co-crystal Ligand (Self-docking)

```
Step 1: Extract co-crystal ligand from PDB
Step 2: Remove ligand from binding site
Step 3: Re-dock ligand into original binding site using AutoDock Vina
Step 4: Compare docked pose with experimental pose
Step 5: Calculate RMSD (root-mean-square deviation)
```

### 2.2 RMSD Calculation

```r
# Using Bio3D R package
library(bio3d)
exp_pdb  <- read.pdb("2RKU_ligand.pdb")     # experimental pose
dock_pdb <- read.pdb("2RKU_docked.pdb")     # docked pose
rmsd_val <- rmsd(exp_pdb$xyz, dock_pdb$xyz, fit=TRUE)
```

### 2.3 Acceptance Criteria

| RMSD | Interpretation | Action |
|------|---------------|--------|
| **< 1.0 Å** | Excellent — docking protocol faithfully reproduces experimental pose | Accept all results |
| **1.0–2.0 Å** | Good — protocol is reliable; minor conformational differences | Accept with note |
| **2.0–3.0 Å** | Borderline — protocol may have systematic errors | Flag results; report with caution |
| **> 3.0 Å** | Poor — protocol cannot reproduce known binding mode | **Reject docking results**; optimize protocol |

---

## 3. Ligand Library for Validation

### Tier 1: Co-crystal Ligands (Positive Controls)

| Ligand | Target | PDB | Expected RMSD | Expected Binding Energy (kcal/mol) |
|--------|--------|-----|---------------|-------------------------------------|
| BI-2536 (R78) | PLK1 | 2RKU | < 2.0 Å | −9.0 to −10.5 |
| RO-3306 | CDK1 | 6GU2 | < 2.0 Å | −8.0 to −9.5 |
| PHA-680626 | PLK1 | 2OWB | < 2.0 Å | −8.5 to −10.0 |

### Tier 2: Known Active Compounds (Expected ≥ Moderate Affinity)

| Ligand | Target | Known Activity | Expected Binding Energy |
|--------|--------|---------------|------------------------|
| Volasertib (BI-6727) | PLK1 | Phase III; IC50 = 0.87 nM | −9.5 to −11.0 |
| Onvansertib (PCM-075) | PLK1 | Phase II; oral PLK1 inhibitor | −9.0 to −10.5 |
| Dinaciclib | CDK1 | Phase II; CDK1 IC50 = 3 nM | −9.0 to −10.5 |
| Milciclib | CDK1 | Phase II HCC; CDK1/2 inhibitor | −8.5 to −10.0 |
| Rigosertib | PLK1 | Phase III; PLK1/PI3K dual | −8.5 to −10.0 |

### Tier 3: Negative Controls (Expected Weak/No Affinity)

| Ligand | Rationale | Expected Binding Energy |
|--------|-----------|------------------------|
| Glucose (GLC) | Small sugar; no kinase affinity | > −4.0 |
| Acetaminophen (TYP) | Common drug; not a kinase inhibitor | > −5.0 |
| Alanine (ALA) | Amino acid; no drug-like properties | > −3.0 |

---

## 4. Expected Docking Results Framework

### 4.1 PLK1 Single-Structure Validation

```
PDB: 2RKU (1.95 Å)
Receptor: PLK1 kinase domain (chain A)
Binding site center: Co-crystal BI-2536 centroid
Grid box: 22 × 22 × 22 Å, exhaustiveness = 32

┌─────────────────────┬──────────┬──────────┬──────────┐
│ Ligand              │ Vina Score│ RMSD     │ Status   │
├─────────────────────┼──────────┼──────────┼──────────┤
│ BI-2536 (re-dock)   │ _______  │ _______  │ ______   │
│ Volasertib          │ _______  │ N/A      │ ______   │
│ Onvansertib         │ _______  │ N/A      │ ______   │
│ Glucose (neg. ctrl) │ _______  │ N/A      │ ______   │
└─────────────────────┴──────────┴──────────┴──────────┘
```

### 4.2 CDK1 Single-Structure Validation

```
PDB: 6GU2 (1.65 Å)
Receptor: CDK1/CyclinB1/CKS2 (chain A)
Binding site center: Co-crystal RO-3306 centroid
Grid box: 22 × 22 × 22 Å, exhaustiveness = 32

┌─────────────────────┬──────────┬──────────┬──────────┐
│ Ligand              │ Vina Score│ RMSD     │ Status   │
├─────────────────────┼──────────┼──────────┼──────────┤
│ RO-3306 (re-dock)   │ _______  │ _______  │ ______   │
│ Dinaciclib          │ _______  │ N/A      │ ______   │
│ Milciclib           │ _______  │ N/A      │ ______   │
│ Acetaminophen (-)   │ _______  │ N/A      │ ______   │
└─────────────────────┴──────────┴──────────┴──────────┘
```

---

## 5. Autodock Vina Command Template

```bash
#!/bin/bash
# PLK1 + BI-2536 re-docking
# Source: PDB 2RKU, chain A

# 1. Prepare receptor
prepare_receptor -r 2RKU_chainA.pdb -o plk1_receptor.pdbqt \
    -A checkhydrogens -U nphs

# 2. Prepare BI-2536 ligand
obabel BI2536.sdf -O BI2536.pdbqt --gen3d

# 3. Run Vina
vina \
    --receptor plk1_receptor.pdbqt \
    --ligand BI2536.pdbqt \
    --center_x 15.2 --center_y 28.7 --center_z 12.4 \
    --size_x 22 --size_y 22 --size_z 22 \
    --exhaustiveness 32 \
    --num_modes 10 \
    --out results/2RKU_BI2536_redock.pdbqt \
    --log results/2RKU_BI2536_redock.log

# 4. Extract best pose (mode 1)
obabel results/2RKU_BI2536_redock.pdbqt -O results/best_pose.pdb -m 1

# 5. Calculate RMSD
echo "RMSD calculation: compare best_pose.pdb with experimental BI-2536"
```

---

## 6. Validation Decision Matrix

| Scenario | RMSD | Vina Score vs. Co-crystal | Action |
|----------|------|--------------------------|--------|
| A | < 1.0 | Better or equal | ✅ Protocol validated — publish all results |
| B | 1.0–2.0 | Better or equal | ✅ Protocol acceptable — report with RMSD note |
| C | 1.0–2.0 | Worse | ⚠️ Flag — re-optimize grid/protonation |
| D | 2.0–3.0 | Any | ⚠️ Borderline — re-optimize; report as preliminary |
| E | > 3.0 | Any | ❌ Protocol rejected — re-test with different parameters |

---

## 7. Publication-Ready Docking Figure Specifications

If validation passes (Scenario A or B):

1. **Figure 9a:** 3D binding pose — PLK1 surface + Volasertib (PyMOL)
2. **Figure 9b:** 2D interaction diagram — key H-bonds + hydrophobic contacts (LigPlot+)
3. **Figure 9c:** RMSD validation — overlay of docked BI-2536 vs. experimental (PyMOL)
4. **Figure 9d:** Binding energy comparison — bar chart of all tested ligands
5. **Supplementary Figure S1:** Full per-residue interaction table

---

## 8. Current Status

| Task | Status |
|------|--------|
| PDB structures identified | ✅ 2RKU (PLK1, 1.95Å), 6GU2 (CDK1, 1.65Å) |
| Protocol corrections applied | ✅ 2OWB→2RKU (BI-2536); 4Y72→6GU2 (RO-3306) |
| AutoDock Vina commands written | ✅ Ready for execution |
| Co-crystal ligand RMSD criteria | ✅ Defined (< 2.0 Å = pass) |
| Negative controls defined | ✅ Glucose, Acetaminophen, Alanine |
| Validation decision matrix | ✅ 5 scenarios A–E |
| **Actual docking execution** | ⬜ Requires Vina + MGLTools + OpenBabel installation |

---

## 9. Recommendations

1. **Install AutoDock Vina 1.2** on a Linux workstation or WSL2
2. **Execute re-docking validation first** (BI-2536 → 2RKU, RO-3306 → 6GU2)
3. **Only proceed to candidate docking AFTER RMSD < 2.0 Å is confirmed**
4. **Report RMSD in manuscript** — this is required for publication
5. If RMSD > 2.0 Å: optimize grid box, protonation states, and exhaustiveness
6. If RMSD remains > 2.0 Å after optimization: report as preliminary computational prediction only

---

**Note:** Until docking is executed, all drug binding predictions remain **computational hypotheses**. The CMap analysis (Step 16) provides independent pharmacological validation that does not depend on docking.
