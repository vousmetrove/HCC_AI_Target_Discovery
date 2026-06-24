# DEG Annotation Summary — TCGA-LIHC HCC

**Generated:**  2026-06-23 22:02 

## Data Sources

- **DESeq2 DEG list:** `results/DEG/DEG_list.csv` (6,301 DEGs)
- **Annotation database:** biomaRt (Ensembl genes, hsapiens_gene_ensembl)
- **Fallback:** org.Hs.eg.db for initial SYMBOL mapping

## Statistics

| Metric | Count | Percentage |
|--------|-------|------------|
| Total DEGs | 6301 | 100% |
| With HGNC symbol | 5213 | 82.7% |
| Without HGNC symbol | 1088 | 17.3% |
| Protein-coding | 3871 | 61.4% |
| lncRNA | 1691 | 26.8% |

## Gene Biotype Distribution

| Biotype | Count |
|---------|-------|
| protein_coding | 3871 |
| lncRNA | 1691 |
| transcribed_unprocessed_pseudogene | 177 |
| processed_pseudogene | 159 |
| TEC | 85 |
| transcribed_processed_pseudogene | 80 |
| transcribed_unitary_pseudogene | 35 |
| IG_V_gene | 23 |
| unprocessed_pseudogene | 21 |
| misc_RNA | 10 |
| snoRNA | 9 |
| miRNA | 6 |
| snRNA | 6 |
| IG_C_gene | 5 |
| unitary_pseudogene | 3 |
| artifact | 2 |
| IG_V_pseudogene | 1 |

## Output Files

- **GO/KEGG gene list:** `results/enrichment/gene_symbol_for_GO.txt` (5212 unique symbols)
- **PPI gene list:** `results/PPI/PPI_gene_list.txt` (5212 unique symbols)
- **Protein-coding DEGs:** `results/DEG/DEG_protein_coding.csv` (3845 genes)
- **Full annotated DEGs:** `results/DEG/DEG_with_SYMBOL_biomaRt.csv` (6301 genes)

## Downstream Pipeline Input

- **Genes for GO/KEGG analysis (Step 04):** 5212
- **Genes for PPI analysis (Step 05):** 5212
- **Protein-coding genes for ML (Step 06/07):** 3845

---
**Pipeline:** HCC_AI_Target_Discovery | **Manuscript:** EI Conference 2026
