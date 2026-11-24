# PAX6 Cervical Cancer — Five-Gene Co-expression Analysis

Reproducible R workflow for analysing the transcriptomic co-expression of **PAX6 and four selected developmental regulators** in The Cancer Genome Atlas cervical squamous cell carcinoma and endocervical adenocarcinoma cohort (TCGA-CESC).

This repository contains the expression matrix, aligned sample manifest, one integrated R analysis script, and documentation required to reproduce the five-gene co-expression analysis and its associated statistical outputs.

## Research objective

The objective of this analysis was to investigate pairwise transcriptomic co-expression among PAX6 and four selected developmental regulators:

- **PAX6** — Paired Box 6
- **SOX1** — SRY-Box Transcription Factor 1
- **SOX2** — SRY-Box Transcription Factor 2
- **NEUROG2** — Neurogenin 2
- **SIX3** — SIX Homeobox 3

The analysis examines whether PAX6 expression is coordinated with these selected genes across primary cervical tumour samples.

The findings are interpreted as transcriptomic associations and do not establish direct transcriptional regulation, physical protein interactions or causality.

## Study cohort and data source

The analysis uses 304 primary tumour samples from the TCGA-CESC project.

| Characteristic | Description |
|---|---|
| Dataset | TCGA-CESC |
| Cohort size | 304 primary tumour samples |
| Unique patients | 304 |
| Sample type | Primary Tumor |
| Data workflow | STAR - Counts |
| Workflow version | `5d8c131bbff59fb0c969217fc1d44e6d1503cd1f` |
| Expression scale | log2(TPM + 1) |
| Genes analysed | PAX6, SOX1, SOX2, NEUROG2, SIX3 |

The expression matrix and aligned GDC manifest contain the same 304 sample identifiers. Each patient is represented by one primary tumour sample.

The manifest provides sample-level provenance, including GDC file identifiers, file names, patient identifiers, sample types, workflow information, file sizes and MD5 checksums.

Data source: [NCI Genomic Data Commons — TCGA-CESC](https://portal.gdc.cancer.gov/projects/TCGA-CESC)

## Repository structure

```text
PAX6-CervicalCancer-MultiOmics/
│
├── README.md
│
├── PAX6_Coexpression_Analysis.R
│
└── data/
    ├── TCGA_304_FiveGene_log2TPM.csv
    └── TCGA_304_Aligned_Manifest.csv
```

The R script automatically creates the following output directories when executed:

```text
Results/
Figures/
Tables/
```

These directories contain the generated statistical results, publication-quality figures and formatted tables.

## Input data

Two CSV files are required.

### 1. Five-gene expression matrix

**Filename:** `TCGA_304_FiveGene_log2TPM.csv`

This file contains the expression values of the five selected genes across 304 primary tumour samples.

Required columns:

| Column | Description |
|---|---|
| `sample_id` | TCGA sample barcode |
| `PAX6` | PAX6 expression |
| `SOX1` | SOX1 expression |
| `SOX2` | SOX2 expression |
| `NEUROG2` | NEUROG2 expression |
| `SIX3` | SIX3 expression |

The expression values are already transformed using:

\[
\log_2(\mathrm{TPM}+1)
\]

**The R script uses the supplied transformed values directly. It does not apply a second logarithmic transformation.**

### 2. Aligned sample manifest

**Filename:** `TCGA_304_Aligned_Manifest.csv`

The manifest contains the corresponding GDC metadata for the 304 expression samples.

Relevant fields include:

- `file_id`
- `file_name`
- `patient_id`
- `sample_id`
- `sample_type`
- `workflow`
- `workflow_version`
- `file_size_bytes`
- `md5sum`

The script verifies sample identity, matching between expression and manifest records, unique patient representation, primary tumour sample type and STAR - Counts workflow consistency.

Both input files are included to support reproduction of the downstream co-expression analysis.

**Data provenance limitation:** The repository contains the prepared log2(TPM + 1) expression matrix and its aligned GDC manifest. The R script does not independently retrieve the original GDC STAR count files or regenerate the upstream TPM transformation. Reproduction of that upstream preprocessing step requires the corresponding GDC source files and processing procedure.

## Statistical methodology

### 1. Pairwise correlation analysis

Two-sided Spearman rank correlation was used to evaluate the relationships among all five genes.

The five-gene analysis produces ten unique pairwise comparisons:

\[
\frac{5(5-1)}{2}=10
\]

For each comparison, the analysis calculates:

- Spearman correlation coefficient (ρ).
- Raw two-sided p-value.
- 95% percentile-bootstrap confidence interval.
- Benjamini–Hochberg-adjusted q-value.

### 2. Multiple-testing correction

The Benjamini–Hochberg false discovery rate procedure was applied across all ten unique pairwise comparisons.

Statistical significance was defined as:

\[
q < 0.05
\]

The mirrored cells in the correlation heatmap represent the same pairwise association and are not counted as separate statistical tests.

### 3. Bootstrap confidence intervals

Pointwise 95% percentile-bootstrap confidence intervals were estimated using 5,000 non-parametric resamples for each unique gene pair.

A fixed random seed is specified in the R script to support reproducibility.

These are pointwise confidence intervals and are not adjusted for multiplicity.

### 4. Kendall's tau-b sensitivity analysis

Kendall's tau-b correlation was calculated for all ten gene pairs as a sensitivity analysis.

Benjamini–Hochberg correction was independently applied across the ten Kendall correlation tests.

This analysis provides a tie-aware assessment of the observed rank-based associations.

### 5. Leave-one-out sensitivity analysis

Leave-one-out analysis was performed for:

- PAX6–SOX2
- PAX6–SIX3

Each of the 304 samples was omitted individually, and the corresponding Spearman correlation coefficient was recalculated.

This analysis evaluates whether the direction and magnitude of the two principal PAX6-centred associations are strongly influenced by individual samples.

## Principal findings

Three unique gene-pair associations remained statistically significant after Benjamini–Hochberg correction.

| Gene pair | N | Spearman ρ | BH-adjusted q-value |
|---|---:|---:|---:|
| PAX6–SOX1 | 304 | 0.011 | 0.845 |
| **PAX6–SOX2** | **304** | **0.296** | **1.41 × 10⁻⁶** |
| PAX6–NEUROG2 | 304 | 0.055 | 0.499 |
| **PAX6–SIX3** | **304** | **0.274** | **6.41 × 10⁻⁶** |
| SOX1–SOX2 | 304 | 0.038 | 0.561 |
| SOX1–NEUROG2 | 304 | 0.054 | 0.499 |
| SOX1–SIX3 | 304 | 0.081 | 0.316 |
| SOX2–NEUROG2 | 304 | 0.125 | 0.074 |
| **SOX2–SIX3** | **304** | **0.160** | **0.0172** |
| NEUROG2–SIX3 | 304 | 0.046 | 0.530 |

The three statistically significant associations were:

1. PAX6–SOX2
2. PAX6–SIX3
3. SOX2–SIX3

Although SOX2–NEUROG2 showed a nominal p-value below 0.05, its BH-adjusted q-value was 0.074. It was therefore not considered statistically significant after multiple-testing correction.

### Sensitivity analysis findings

Kendall's tau-b supported the two principal PAX6-centred associations:

| Gene pair | Kendall's tau-b | BH-adjusted q-value |
|---|---:|---:|
| PAX6–SOX2 | 0.198 | 2.63 × 10⁻⁶ |
| PAX6–SIX3 | 0.190 | 7.57 × 10⁻⁶ |

Leave-one-out analysis showed that both associations remained positive across all 304 sample omissions.

| Gene pair | Minimum ρ | Maximum ρ |
|---|---:|---:|
| PAX6–SOX2 | 0.290 | 0.308 |
| PAX6–SIX3 | 0.267 | 0.283 |

Together, these findings support modest positive co-expression associations between PAX6 and SOX2 and between PAX6 and SIX3 in the TCGA-CESC primary tumour cohort.

These observations should not be interpreted as evidence of direct transcriptional regulation or causal biological interactions.

## Generated figures and tables

The integrated R script generates the following outputs.

### Figure 5 — PAX6-centred co-expression scatterplots

Four scatterplots represent the relationships between PAX6 and:

- SOX1
- SOX2
- NEUROG2
- SIX3

Each panel reports the Spearman correlation coefficient, pointwise 95% bootstrap confidence interval and BH-adjusted q-value across the ten unique comparisons.

No linear regression lines are displayed.

Output files:

```text
Figures/
├── Figure5_PAX6_Coexpression_Scatterplots.png
└── Figure5_PAX6_Coexpression_Scatterplots.tiff
```

### Figure 6 — Five-gene correlation heatmap

A symmetric 5 × 5 heatmap displays all pairwise Spearman correlation coefficients.

- All 25 matrix cells are displayed.
- Diagonal cells represent self-correlation (ρ = 1.00).
- Mirrored cells represent identical pairwise associations.
- Asterisks indicate BH-adjusted q < 0.05.

Output files:

```text
Figures/
├── Figure6_FiveGene_Correlation_Heatmap.png
└── Figure6_FiveGene_Correlation_Heatmap.tiff
```

The TIFF figures are exported at 600 dpi with LZW compression.

### Table 4 — Complete pairwise co-expression results

Table 4 contains all ten unique gene-pair comparisons and reports:

- Sample size.
- Spearman correlation coefficient.
- 95% bootstrap confidence interval.
- Raw p-value.
- BH-adjusted q-value.

The Excel workbook is formatted using Times New Roman.

Output files:

```text
Tables/
├── Table4_Coexpression.xlsx
├── Table4_Coexpression_Manuscript.csv
└── Table4_Coexpression_FullPrecision.csv
```

## Software and dependencies

The analysis was conducted using R version 4.6.1.

Required R packages:

- `ggplot2`
- `openxlsx`

Install the dependencies using:

```r
install.packages(c("ggplot2", "openxlsx"))
```

Base R functions are used for correlation calculations, multiple-testing correction, bootstrap resampling and sensitivity analyses.

## Reproducing the analysis

### Step 1. Clone the repository

```bash
git clone https://github.com/atarskushwah/PAX6-CervicalCancer-MultiOmics.git
```

### Step 2. Set up the input directory

Ensure the repository contains:

```text
data/
├── TCGA_304_FiveGene_log2TPM.csv
└── TCGA_304_Aligned_Manifest.csv
```

### Step 3. Open R or RStudio at the repository root

### Step 4. Install required packages

```r
install.packages(c("ggplot2", "openxlsx"))
```

### Step 5. Execute the analysis

```r
source("PAX6_Coexpression_Analysis.R")
```

Alternatively, run from the terminal:

```bash
Rscript PAX6_Coexpression_Analysis.R
```

The script validates the input files and stops if required columns, sample counts, sample identifiers, metadata or statistical benchmark checks fail.

### Step 6. Retrieve generated outputs

After successful execution, the following folders will be created:

```text
Results/
Figures/
Tables/
```

The analysis also exports R session information to support software-environment documentation.

## Data availability and responsible use

The expression matrix and aligned sample manifest are derived from the publicly accessible TCGA-CESC project through the NCI Genomic Data Commons.

The repository provides the prepared five-gene expression matrix and its corresponding sample-level manifest.

The original underlying GDC source data should be referenced using the supplied file identifiers and manifest metadata.

TCGA data must be used in accordance with applicable GDC data-use policies. Sample-level identifiers must not be used for re-identification.

No patient-level clinical peripheral-blood qRT-PCR data, private laboratory records or confidential institutional information are included in this repository.

## Reproducibility limitations

This repository is designed to reproduce the downstream five-gene co-expression analysis using the supplied expression matrix and aligned manifest.

It does not independently reproduce the upstream GDC download, STAR count processing or TPM extraction procedure.

Reproducing the full upstream data-processing workflow requires the corresponding original GDC files and documented preprocessing steps.

Bootstrap confidence intervals are pointwise percentile intervals and may differ from previously reported intervals if the random seed, resampling implementation or software environment changes.

## Scientific interpretation

The analysis identifies positive transcriptomic associations among selected developmental regulators in primary cervical tumour samples.

The strongest PAX6-centred associations were observed with SOX2 and SIX3. These findings provide evidence of coordinated expression patterns within the analysed TCGA-CESC cohort.

However, correlation does not establish direct transcriptional regulation, physical interaction or causality. Functional experimental validation is necessary to determine the biological significance of these relationships.

---

**Repository:** [PAX6-CervicalCancer-MultiOmics](https://github.com/atarskushwah/PAX6-CervicalCancer-MultiOmics)

**Primary analysis:** `PAX6_Coexpression_Analysis.R`

**Dataset:** TCGA-CESC, 304 primary tumour samples
