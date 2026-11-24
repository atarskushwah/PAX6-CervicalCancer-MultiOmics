# Supplementary Figure S1 — PAX6 annotation and STRING associations

## Scientific scope

Panels A–C summarise **single-gene Enrichr annotations** for PAX6 (GO Biological Process, GO Molecular Function and KEGG pathway membership). These are **not** enrichment results obtained from a differentially expressed gene set and do **not** establish activation of the displayed terms in cervical cancer.

Panel D refers to **predicted STRING functional associations**, not experimentally validated PAX6 protein binding or transcriptional control. The repository SVG includes a *schematic* first-shell partner diagram; it is not a pixel-for-pixel reproduction of the original STRING network. The manuscript's figure panel uses the archived STRING network image.

## Source of the corrected figure

The manuscript's Supplementary Figure S1 was updated on 24 November 2026 from the author-supplied `Figure4_Enrichment_STRING.pptx` and `Figure4_Enrichment_STRING.tiff`, combined with the network image already present in the author's supplementary-figure Rv2.0 materials. The user-provided TIFF/PPTX suffered from right-side clipping and crowded term labels; the corrected publication figure restructures the four panels without introducing any new gene association claims.

**Figure legend**

**Supplementary Figure S1.** PAX6 functional annotations and predicted protein functional associations. (A) Gene Ontology (GO) Biological Process terms, (B) GO Molecular Function terms and (C) KEGG pathway membership returned from an Enrichr query containing PAX6 alone. Horizontal distances show −log10 of the source-output nominal p-values, and the single queried gene overlaps each displayed term. These annotations are descriptive and do not demonstrate significantly activated pathways in cervical cancer. (D) Archived STRING network illustrating predicted functional associations of PAX6 with SOX1, SOX2, SIX3, NEUROG2, TBR1, MITF, MAF, INS, BANF1 and DCAF7. STRING edges denote functional association evidence, not direct physical interaction or causal regulation. The source graphic does not constitute an independent tumour-specific functional validation.

## Figure numbering correspondence

| Repository or original output | Revised manuscript designation |
|---|---|
| `Figure5_PAX6_Coexpression_Scatterplots` output from the archived R analysis | **Main Figure 4** |
| `Figure6_FiveGene_Correlation_Heatmap` output from the archived R analysis | **Supplementary Figure S2** |
| `Table4_Coexpression` output from the archived R analysis | **Main Table 3** |
| Enrichr/STRING annotation graphic | **Supplementary Figure S1** |

The historical filenames in `PAX6_Coexpression_Analysis.R` are preserved for reproducibility; use the manuscript labels when citing figures in the revised article.

## Public vs restricted data

The existing repository contains **TCGA-CESC public expression and manifest data only**. It must not include the confidential patient-level clinical survival, HPV, treatment-response or RT-qPCR datasets. Clinical survival estimates in the manuscript are analyses of a separate, restricted cohort (103 cases, 30 deaths); their individual-level inputs are not released here.

## Technical note

This repository provides the aligned TCGA matrix and downstream co-expression script. GO/KEGG annotations and STRING associations are portal-derived exports, not direct outputs of that co-expression script. Exact recreation of the four-panel source figure requires the original saved Enrichr/STRING export settings and records.

