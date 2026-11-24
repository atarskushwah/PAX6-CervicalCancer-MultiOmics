# PAX6 Cervical Cancer: Five-Gene Co-expression Analysis

This repository contains **one R script** that reproduces the five-gene co-expression analysis, Figures 5 and 6, Table 4, and the associated sensitivity analyses for the PAX6 cervical cancer study.

## Study scope

The analysis evaluates transcriptomic co-expression among five genes in the TCGA-CESC primary-tumour cohort:

- `PAX6`
- `SOX1`
- `SOX2`
- `NEUROG2`
- `SIX3`

The supplied cohort contains 304 primary-tumour samples from 304 patients. The accompanying GDC manifest identifies all 304 records as `Primary Tumor`, processed using the `STAR - Counts` workflow (workflow version `5d8c131bbff59fb0c969217fc1d44e6d1503cd1f`). The expression matrix and manifest contain the same 304 sample identifiers in the same row order.

## Repository contents

```text
PAX6-CervicalCancer-MultiOmics/
├── README.md
├── PAX6_Coexpression_Analysis.R
└── data/                                  # local input files; not committed
    ├── TCGA_304_FiveGene_log2TPM.csv
    └── TCGA_304_Aligned_Manifest.csv
```

The `Results/`, `Figures/` and `Tables/` directories are created automatically when the script runs. They do not need to be maintained as separate source-code modules.

## Input files

Two local CSV files are required.

### 1. Five-gene expression matrix

File: `data/TCGA_304_FiveGene_log2TPM.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | TCGA sample identifier |
| `PAX6` | PAX6 expression |
| `SOX1` | SOX1 expression |
| `SOX2` | SOX2 expression |
| `NEUROG2` | NEUROG2 expression |
| `SIX3` | SIX3 expression |

The supplied expression values are already on the **log2(TPM + 1)** scale. The script uses these values as supplied and does not apply a second log transformation. It checks that the matrix contains exactly 304 unique samples, all five gene columns are numeric and finite, and no values are missing or negative.

### 2. Aligned sample manifest

File: `data/TCGA_304_Aligned_Manifest.csv`

The script uses these fields:

- `sample_id`
- `patient_id`
- `sample_type`
- `workflow`
- `workflow_version`

The manifest is used to check sample identity, one patient per sample, primary-tumour sample type and workflow consistency. The script requires the expression matrix and manifest to contain exactly the same sample IDs. It aligns the manifest to the expression matrix by `sample_id` before validating the records.

**Data handling:** the two input CSVs are not included in the public repository by default. Keep them in the local `data/` folder. Do not commit clinical or patient-level data, credentials, or any file whose redistribution is restricted. The GDC manifest and expression matrix should only be redistributed if their source and applicable data-use terms permit it.

## Statistical methods

- Cohort: 304 TCGA-CESC primary-tumour samples.
- Expression scale: log2(TPM + 1), supplied as the input matrix.
- Genes: PAX6, SOX1, SOX2, NEUROG2 and SIX3.
- Primary method: two-sided Spearman rank correlation.
- Testing family: ten unique pairwise comparisons among the five genes.
- Multiple testing: Benjamini-Hochberg false-discovery-rate correction across the ten unique comparisons.
- Significance threshold: BH-adjusted `q < 0.05`.
- Confidence intervals: pointwise 95% percentile-bootstrap intervals from 5,000 non-parametric resamples per pair.
- Sensitivity analysis: Kendall's tau-b for all ten pairs, with BH correction across the ten Kendall tests.
- Influence analysis: leave-one-out Spearman correlations for PAX6-SOX2 and PAX6-SIX3, omitting each of the 304 samples once.

The mirrored cells in Figure 6 display the same pairwise association and are not counted as separate tests. The diagonal cells represent self-correlation (rho = 1.00).

The bootstrap uses a fixed seed (`20261124`, with a pair-specific offset) so the intervals are reproducible for this script version and R random-number generator implementation. Confidence intervals can differ slightly from intervals generated using a different seed or bootstrap implementation; the script does not overwrite calculated intervals with previously reported values.

## Principal findings

The supplied expression matrix reproduces the following Spearman coefficients and BH-adjusted q-values:

| Gene pair | Spearman rho | BH-adjusted q-value |
|---|---:|---:|
| PAX6-SOX1 | 0.011 | 0.845 |
| PAX6-SOX2 | 0.296 | 1.41 x 10^-6 |
| PAX6-NEUROG2 | 0.055 | 0.499 |
| PAX6-SIX3 | 0.274 | 6.41 x 10^-6 |
| SOX1-SOX2 | 0.038 | 0.561 |
| SOX1-NEUROG2 | 0.054 | 0.499 |
| SOX1-SIX3 | 0.081 | 0.316 |
| SOX2-NEUROG2 | 0.125 | 0.074 |
| SOX2-SIX3 | 0.160 | 0.0172 |
| NEUROG2-SIX3 | 0.046 | 0.530 |

Three associations are significant after BH correction: **PAX6-SOX2, PAX6-SIX3 and SOX2-SIX3**. SOX2-NEUROG2 has a nominal p-value below 0.05 but an adjusted q-value of 0.074 and is not significant after correction.

The results indicate modest positive co-expression associations. They do not establish direct transcriptional regulation, physical interaction or causality.

## Software requirements

- R (the preceding analysis was conducted using R 4.6.1; record the version used for each final run).
- R packages: `ggplot2` and `openxlsx`.

Install the packages in R:

```r
install.packages(c("ggplot2", "openxlsx"))
```

Base R functions are used for correlation tests, BH adjustment, bootstrap resampling and leave-one-out calculations.

## Run the analysis

From the repository root in R or RStudio:

```r
source("PAX6_Coexpression_Analysis.R")
```

Or from a terminal:

```bash
Rscript PAX6_Coexpression_Analysis.R
```

If the inputs are stored elsewhere, set their paths before running:

```r
Sys.setenv(
  PAX6_EXPRESSION_FILE = "/path/to/TCGA_304_FiveGene_log2TPM.csv",
  PAX6_MANIFEST_FILE = "/path/to/TCGA_304_Aligned_Manifest.csv"
)
source("PAX6_Coexpression_Analysis.R")
```

The script stops if the sample count, sample IDs, manifest metadata, expression columns or numerical benchmark checks fail. It does not silently remove samples or substitute reference results for calculated results.

## Generated outputs

### Figures

- `Figures/Figure5_PAX6_Coexpression_Scatterplots.png`
- `Figures/Figure5_PAX6_Coexpression_Scatterplots.tiff`
- `Figures/Figure6_FiveGene_Correlation_Heatmap.png`
- `Figures/Figure6_FiveGene_Correlation_Heatmap.tiff`

Both TIFF figures are saved at 600 dpi with LZW compression. Figure 5 contains four PAX6-centred scatterplots and no linear fitted regression lines. Figure 6 displays the full symmetric 5 x 5 correlation matrix.

### Tables

- `Tables/Table4_Coexpression.xlsx` — Times New Roman workbook with manuscript table, full-precision results, Kendall tau-b results, leave-one-out results, input QC and analysis notes.
- `Tables/Table4_Coexpression_Manuscript.csv` — manuscript-formatted Table 4.
- `Tables/Table4_Coexpression_FullPrecision.csv` — full-precision Spearman results.

### Additional results and validation

- `Results/Spearman_All10_FullPrecision.csv`
- `Results/Kendall_TauB_Sensitivity.csv`
- `Results/LeaveOneOut_Sensitivity_Summary.csv`
- `Results/LeaveOneOut_Sensitivity_Details.csv`
- `Results/Input_QC_Summary.csv`
- `Results/R_sessionInfo.txt`

The leave-one-out detail file contains omission indices, not sample identifiers.

## Reproducibility and provenance

The script records aggregate input checks and the R session information. The accompanying manifest identifies the GDC workflow metadata, sample type and source file identifiers. The script validates the manifest-to-expression sample mapping but does not regenerate the log2(TPM + 1) expression matrix from the underlying GDC STAR count files. To reproduce that upstream preprocessing step independently, retain and document the corresponding GDC data retrieval and TPM extraction procedure.
