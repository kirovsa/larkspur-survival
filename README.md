# Larkspur Survival — OS/PFS Survival Analysis

Cox proportional hazards analysis of predefined gene signatures in TCGA cohorts,
using cBioPortal data and Quarto reporting.

---

## Project Overview

This pipeline:

1. Reads gene signature definitions from `metadata/signatures.yaml`
2. Reads cohort definitions (cBioPortal study IDs) from `metadata/cohorts.yaml`
3. Loads clinical + mRNA expression data from the `cbioportal/` directory
4. Computes per-sample **signature scores** (z-scored mean expression)
5. Stratifies patients into **High / Low** groups at the cohort median
6. Fits **Cox proportional hazards** models for OS and PFS
7. Tests the PH assumption (Schoenfeld residuals)
8. Generates **Kaplan–Meier curves** and **forest plots**
9. Renders a self-contained HTML report via Quarto

---

## Repository Structure

```
larkspur-survival/
├── analysis/
│   └── survival_analysis.qmd    # Main Quarto report
├── cbioportal/
│   ├── README.md                # Data format docs & download instructions
│   └── brca_tcga_pub/           # Example synthetic TCGA-BRCA data
│       ├── data_clinical_patient.txt
│       ├── data_clinical_sample.txt
│       ├── data_mrna_seq_v2_rsem.txt
│       └── meta_study.txt
├── metadata/
│   ├── cohorts.yaml             # Cohort definitions (study_id → folder path)
│   └── signatures.yaml          # Gene signature definitions
├── R/
│   └── utils.R                  # Data loading, scoring, Cox PH, plotting helpers
├── _quarto.yml                  # Quarto project configuration
├── renv.lock                    # R package lockfile (reproducible environment)
└── .Rprofile                    # R project settings
```

---

## Quick Start

### 1. Prerequisites

- **R ≥ 4.3** — https://cran.r-project.org/
- **Quarto ≥ 1.4** — https://quarto.org/docs/get-started/

### 2. Install R packages

```r
install.packages("renv")
renv::restore()   # installs packages from renv.lock
```

### 3. Render the report (example / synthetic data)

```bash
cd analysis
quarto render survival_analysis.qmd
```

The report is written to `analysis/survival_analysis.html`.

### 4. Use real TCGA data

1. Download a TCGA study from cBioPortal (https://www.cbioportal.org/) — click
   **Download → Bulk Data** on any study page.
2. Extract the archive into `cbioportal/<study_id>/` (e.g. `cbioportal/luad_tcga/`).
3. Add a new entry to `metadata/cohorts.yaml`:

   ```yaml
   - study_id: luad_tcga
     label: "Lung Adenocarcinoma (TCGA)"
     cancer_type: LUAD
     clinical_file: data_clinical_patient.txt
     expression_file: data_mrna_seq_v2_rsem.txt
     endpoints: [OS, PFS]
   ```

4. Re-render: `quarto render analysis/survival_analysis.qmd`

---

## Configuration Reference

### `metadata/cohorts.yaml`

| Key | Description |
|-----|-------------|
| `study_id` | Folder name under `cbioportal/` |
| `label` | Human-readable name for the report |
| `cancer_type` | Short acronym (BRCA, LUAD, …) |
| `clinical_file` | Filename of the clinical patient file |
| `expression_file` | Filename of the mRNA expression file |
| `endpoints` | List of `OS` and/or `PFS` |

### `metadata/signatures.yaml`

| Key | Description |
|-----|-------------|
| `name` | Unique identifier (used in filenames) |
| `label` | Human-readable label |
| `genes` | Hugo gene symbol list |
| `scoring` | `mean` (default), `median`, or `pc1` |
| `direction` | `positive` (high = bad) or `negative` (high = good) |

---

## Methods Summary

| Step | Details |
|------|---------|
| Signature scoring | Mean of per-gene z-scores across all samples in cohort |
| Stratification | Median split → High / Low |
| Model | Unadjusted Cox PH (`coxph`, `survival` package) |
| PH test | Schoenfeld residuals (`cox.zph`) |
| Plots | KM curves (`survminer`), forest plots (`ggplot2`) |
| Report | Quarto HTML with collapsible code |

---

## Output

The rendered report (`analysis/survival_analysis.html`) contains:

- Configuration summary
- Data overview table
- Per-cohort sections with:
  - Cox PH coefficient table
  - PH assumption test result
  - Kaplan–Meier curves
- Cross-cohort summary table (HR, 95 % CI, *p*-value)
- Forest plots per cohort × endpoint
- Session info / software versions

---

## License

MIT
