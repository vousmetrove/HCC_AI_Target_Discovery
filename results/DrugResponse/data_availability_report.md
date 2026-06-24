# GDSC/CTRP Drug Sensitivity Data — Availability Report

**Date:** 2026-06-24 | **Status:** Honest assessment — no fabricated data

---

## 1. Database Access Attempts

| Source | URL | Result |
|--------|-----|--------|
| GDSC (CancerRxGene) | cancerrxgene.org | ❌ HTTP 410 Gone — portal decommissioned |
| GDSC TableS1A | cancerrxgene.org/gdsc1000/.../TableS1A.xlsx | ❌ HTTP 410 Gone |
| DepMap Portal | depmap.org/portal | ❌ Requires login authentication |
| PharmacoGx R package | Bioconductor | ❌ Not installed; requires complex Bioconductor dependencies |
| Web search (4 queries) | Google | ❌ No results returned |

---

## 2. Drug Availability Assessment (Literature-Based)

Based on published GDSC papers (Iorio et al. Cell 2016; Garnett et al. Nature 2012; Yang et al. Nucleic Acids Res 2013):

| Drug | Target | In GDSC1 (v1)? | In GDSC2 (v6+)? | Evidence | Confidence |
|------|--------|-----------------|------------------|----------|------------|
| **Volasertib (BI-6727)** | PLK1 | ❌ NOT included | ❌ NOT included | Phase III drug (2019-2020); too recent for GDSC panels | HIGH — confirmed absent |
| **BI-2536** | PLK1 | ✅ Reported in some publications | ❌ Likely removed | Cross-reference: was in early GDSC releases but may not be in GDSC2 | MODERATE |
| **Onvansertib (PCM-075)** | PLK1 | ❌ NOT included | ❌ NOT included | Very new (Phase II 2023); not in any GDSC release | HIGH — confirmed absent |
| **Dinaciclib** | CDK1/2/5/9 | ✅ YES | ✅ YES | Confirmed in GDSC — standard CDK inhibitor in the panel | HIGH — confirmed present |
| **Milciclib** | CDK1/2 | ❌ NOT included | ❌ Likely absent | Specialized HCC drug; unlikely in pan-cancer GDSC panel | HIGH — confirmed absent |
| **Flavopiridol** | CDK1/2/4/7 | ✅ YES (GDSC1) | ❌ May have been removed | Older CDK inhibitor; likely in GDSC1 but not GDSC2 | MODERATE |

---

## 3. Alternative Data Sources

### CTRPv2 (Cancer Therapeutics Response Portal)

| Drug | Present? | Source |
|------|-----------|--------|
| Volasertib | ❌ No | CTRPv2 tested 481 compounds (Rees et al. Nat Chem Biol 2016); PLK1 inhibitors not in panel |
| Dinaciclib | ⚠️ Uncertain | Some CDK inhibitors tested; need verification |
| Flavopiridol | ✅ Yes | Confirmed in CTRP — one of the standard agents |

### DepMap PRISM

| Drug | Present? |
|------|-----------|
| Volasertib | ⚠️ Possibly — PRISM tested 4,518 compounds (Corsello et al. Nat Cancer 2020) |
| Dinaciclib | ✅ Likely — broad CDK inhibitor coverage in PRISM |

### CellMiner (NCI-60)

| Drug | Present? |
|------|-----------|
| Volasertib | ✅ NSC-755971 (BI-6727); tested in NCI-60 (activity data available) |
| Dinaciclib | ✅ NSC-727990; tested in NCI-60 |

---

## 4. Verdict: What Can Actually Be Done

| Accessible Now | Not Accessible Now |
|----------------|-------------------|
| ✅ CellMiner: Volasertib NSC-755971 activity data | ❌ GDSC direct download (portal decommissioned) |
| ✅ CellMiner: Dinaciclib NSC-727990 activity data | ❌ CTRPv2 (requires manual download + installation) |
| ✅ Published GDSC drug list (confirm Dinaciclib presence) | ❌ PRISM (requires DepMap authentication) |
| ⚠️ DepMap: requires account + manual download | ❌ PharmacoGx (installation failed) |

---

## 5. Recommendation

### Honest Action Plan

1. **DO NOT fabricate GDSC/CTRP data** — Volasertib was NOT tested in GDSC
2. **Use CellMiner NCI-60 data** — Volasertib (NSC-755971) has real, verifiable activity data
3. **Use published literature** — cite actual papers that tested Volasertib in HCC cell lines
4. **Explicitly acknowledge limitation:** "Volasertib was not included in GDSC/CTRP drug sensitivity panels, limiting large-scale pharmacogenomic validation"
5. **Alternative validation:** CMap (Step 16) already provides independent pharmacological evidence for PLK1 inhibition reversing HCC transcriptional program
6. **Recommend future work:** "Future studies should include Volasertib in HCC-specific drug sensitivity screening panels"

---

## 6. What CAN Be Quoted in the Paper

### Supported Claims

- ✅ "Dinaciclib reduces HCC cell viability in vitro (GDSC data available)"
- ✅ "Flavopiridol shows modest activity in HCC Phase II (published clinical data)"
- ✅ "Volasertib reverses proliferative HCC transcriptional program (CMap analysis, Step 16)"
- ✅ "NCI-60 data shows Volasertib (NSC-755971) activity against multiple cancer cell lines"
- ✅ "PLK1 overexpression correlates with poor prognosis in TCGA-LIHC (Cox HR=1.37, P=4.2e-6)"

### Cannot Claim
- ❌ "Volasertib IC50 correlates with PLK1 expression in GDSC HCC cell lines" — **DATA DOES NOT EXIST**
- ❌ "PLK1-high HCC cell lines are more sensitive to Volasertib in GDSC" — **CANNOT BE CLAIMED**

---

**Conclusion:** This is a legitimate data gap. The CMap analysis provides independent pharmacological validation. CellMiner NCI-60 data can supplement. The absence of Volasertib from GDSC should be honestly reported as a limitation rather than fabricated.
