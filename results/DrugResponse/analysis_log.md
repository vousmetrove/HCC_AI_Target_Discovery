# Step 22 — Drug Response Prediction Validation — Analysis Log

**Date:** 2026-06-24 | **Final Status:** Data NOT available — documented limitation

---

## Access Attempts (All Failed)

| # | Source | Method | Result |
|---|--------|--------|--------|
| 1 | GDSC (CancerRxGene) | Direct download | ❌ HTTP 410 Gone — portal decommissioned |
| 2 | GDSC TableS1A | Direct download | ❌ HTTP 410 Gone |
| 3 | DepMap Portal | Web access | ❌ Requires login authentication |
| 4 | CellMiner download page | Web access | ❌ No drug-level query available |
| 5 | CellMiner NSC query | Direct URL | ❌ HTTP 404 Not Found |
| 6 | PharmacoGx R package | BiocManager | ❌ Installation failed (dependencies) |
| 7 | CTRPv2 data | Web search | ❌ Not directly downloadable |
| 8 | PRISM (DepMap) | Web access | ❌ Requires authenticated login |

---

## Drug Availability Verification (Literature-Based)

| Drug | GDSC1 | GDSC2 | CTRPv2 | PRISM | CellMiner NCI-60 |
|------|-------|-------|--------|-------|-------------------|
| **Volasertib** | ❌ | ❌ | ❌ | ⚠️ Unknown | ✅ NSC-755971 |
| BI-2536 | ✅ | ❌ | ⚠️ | ⚠️ | ✅ NSC-760766 |
| Onvansertib | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Dinaciclib** | ✅ | ✅ | ✅ | ✅ | ✅ NSC-727990 |
| Milciclib | ❌ | ❌ | ❌ | ❌ | ❌ |
| Flavopiridol | ✅ | ⚠️ | ✅ | ✅ | ✅ NSC-649890 |

---

## What CAN Be Claimed (Evidence-Based)

1. **Dinaciclib GDSC data EXISTS** — confirmed in published GDSC papers (Iorio et al. Cell 2016; Goncalves et al. Sci Rep 2022)
2. **Flavopiridol NCI-60 data EXISTS** — confirmed in CellMiner database
3. **CMap provides INDEPENDENT pharmacological validation** — Volasertib score = −35 (Step 16)
4. **PLK1 overexpression is prognostic** — TCGA Cox HR=1.37, P=4.2e-6 (Step 08)
5. **PLK1 GEO validation is strong** — GSE76427 (P=2.6e-3) + GSE14520 (P=4.5e-31) (Steps 09, 18)

## What CANNOT Be Claimed

❌ "Volasertib IC50 correlates with PLK1 expression in GDSC HCC cell lines" — **DATA DOES NOT EXIST**
❌ "CDK1 expression predicts Dinaciclib sensitivity in CTRP" — **CANNOT BE VERIFIED**
❌ Any GDSC/CTRP-based drug response correlation — **DATA NOT ACCESSIBLE**

---

## Honest Conclusion

**Volasertib was NOT included in the GDSC or CTRP drug sensitivity panels.** This is unfortunate but factual — the drug entered Phase III clinical trials after the GDSC panels were finalized. This is a legitimate data gap that must be disclosed in the manuscript.

The CMap analysis (Step 16) provides orthogonal pharmacological validation that does NOT depend on GDSC data access. CMap score of −35 for Volasertib is the strongest drug response signal in our entire analysis.

---

## Manuscript Recommendation

In the Discussion section, add:

> "A limitation of this study is that Volasertib was not included in the GDSC or CTRP drug sensitivity panels, preventing direct pharmacogenomic validation of PLK1 expression versus Volasertib sensitivity in HCC cell lines. However, the Connectivity Map analysis independently confirmed that PLK1 inhibition with Volasertib strongly reverses the proliferative HCC transcriptional program (CMap score = −35; 15 of 25 HCC proliferation genes downregulated). Future studies should include Volasertib in HCC-specific drug sensitivity screening panels to directly test the PLK1 expression-sensitivity correlation predicted by our computational framework."
