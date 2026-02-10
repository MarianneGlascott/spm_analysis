# Manuscript 4 — SPM impacts on kelp zoospore motility

[![Reproducible R](https://img.shields.io/badge/Reproducible-renv-blue)](https://rstudio.github.io/renv/)
[![GitHub Repo](https://img.shields.io/badge/GitHub-spm_analysis-black)](https://github.com/MarianneGlascott/spm_analysis)

## Project overview

This repository contains the **fully reproducible analysis pipeline** for *Manuscript 4*, examining how **suspended particulate matter (SPM)** influences **kelp zoospore motility** under controlled experimental conditions.

The project is designed to support:

* transparent, file-based workflows
* unbalanced experimental designs
* explicit quality control decisions
* reviewer-ready reproducibility

All analytical decisions are logged, and **no result is considered valid unless written to disk**.

---

## Experimental design (fixed at analysis stage)

* **Response (primary):** `motility_ratio`
* **Experimental unit:** well
* **Fixed effects:**

  * Site × Species (crossed)
  * Season (categorical)
* **Random effects:**

  * `(1 | culture / experiment)`
* **Design notes:**

  * Unbalanced designs are expected
  * Estimated marginal means used where appropriate

No exclusion thresholds are applied without explicit documentation.

---

## Reproducibility principles

This project follows strict reproducibility rules:

* `renv` for package version control
* `here::here()` for path safety
* explicit script boundaries (setup → import → QC → EDA → models → figures)
* no silent filtering or row dropping
* all outputs written to disk with informative filenames

---

## Repository structure

```
spm_analysis/
├── scripts/              # Analysis scripts (numbered, single-purpose)
│   ├── 00_project_setup.R
│   ├── 01_import_clean.R
│   ├── 02_qc_report.R
│   ├── 03_eda.R
│   ├── 04_models_motility.R
│   ├── 05_models_settlement.R
│   ├── 06_figures_pub.R
│   ├── 07_tables_pub.R
│   └── 08_writeup_helpers.R
│
├── R/                    # Shared helpers and style definitions
│   ├── helpers.R
│   └── kelp_style.R
│
├── data_raw/             # Raw input data (not tracked)
├── data_processed/       # Cleaned datasets (not tracked)
├── outputs/              # Logs, summaries, QC reports (not tracked)
├── figures/              # Generated figures (not tracked)
├── tables/               # Generated tables (not tracked)
├── models/               # Saved model objects (.rds) (not tracked)
├── renv.lock             # Locked R package versions
├── .gitignore
└── README.md
```

---

## How to run the analysis

### 1. Clone the repository

```bash
git clone https://github.com/MarianneGlascott/spm_analysis.git
cd spm_analysis
```

### 2. Restore the R environment

```r
renv::restore()
```

### 3. Add raw data

Place the agreed raw dataset(s) into:

```
data_raw/
```

(No raw data are stored in the repository.)

### 4. Run scripts in order

Scripts are designed to be run **sequentially**:

```r
source("scripts/00_project_setup.R")
source("scripts/01_import_clean.R")
source("scripts/02_qc_report.R")
# ... continue in numeric order
```

Each script:

* performs one clearly defined task
* writes all outputs to disk
* logs decisions and warnings

---

## Data governance

* Raw data remain local and version-controlled externally
* Processed data and outputs are regenerated on demand
* This repository contains **no sensitive or proprietary data**

---

## Citation and status

This repository supports PhD thesis chapter and associated manuscript(s). It is under active development and should not be cited without contacting the author.

**Author:** Marianne Glascott
**Affiliation:** University of Sussex
**Manuscript:** *SPM impacts on kelp zoospore motility*

---

## Notes for reviewers and collaborators

* All modelling choices are explicitly encoded in scripts
* No hidden preprocessing steps exist
* Quality control rules are defined *before* exclusion
* Random effects structure is fixed by experimental design

Questions or issues should be raised via GitHub Issues.
