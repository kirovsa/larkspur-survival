# cBioPortal Data Directory

This directory holds TCGA cohort data downloaded from [cBioPortal](https://www.cbioportal.org/).

## Directory Structure

```
cbioportal/
└── <study_id>/                        # e.g. brca_tcga_pub, luad_tcga
    ├── meta_study.txt                 # Study metadata
    ├── data_clinical_patient.txt      # Patient-level clinical data
    ├── data_clinical_sample.txt       # Sample-level clinical data
    └── data_mrna_seq_v2_rsem.txt      # mRNA expression (log2 RSEM)
```

## Included Example Data

`brca_tcga_pub/` — Synthetic data (150 patients) mimicking TCGA Breast Cancer structure.
Used for demonstration and testing. Replace with real data downloaded from cBioPortal.

## Downloading Real TCGA Data from cBioPortal

1. Visit https://www.cbioportal.org/
2. Select the desired TCGA study (e.g., *Breast Invasive Carcinoma (TCGA, PanCancer Atlas)*)
3. Click **Download** → **Bulk Data** (or use the API)
4. Extract the archive and place the study folder here, matching the `study_id` in `metadata/cohorts.yaml`

### Required Files per Study

| File | Description | Required Columns |
|------|-------------|-----------------|
| `data_clinical_patient.txt` | Patient clinical data | `PATIENT_ID`, `OS_STATUS`, `OS_MONTHS`, `PFS_STATUS`, `PFS_MONTHS` |
| `data_mrna_seq_v2_rsem.txt` | mRNA expression | `Hugo_Symbol`, `Entrez_Gene_Id`, one column per sample |

### cBioPortal File Format

Clinical files have **4 comment lines** (starting with `#`) before the column header:

```
#Display Name Row
#Column ID Row
#Data Type Row      (STRING / NUMBER)
#Priority Row       (1 = shown by default)
PATIENT_ID   OS_STATUS   OS_MONTHS   ...
data...
```

### Survival Status Coding

| Column | Event Value | Censored Value |
|--------|------------|----------------|
| `OS_STATUS` | `1:DECEASED` | `0:LIVING` |
| `PFS_STATUS` | `1:PROGRESSED` | `0:NOT PROGRESSED` |
