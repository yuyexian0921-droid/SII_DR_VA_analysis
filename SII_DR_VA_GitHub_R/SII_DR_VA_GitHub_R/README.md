# SII–DR–Visual Decline: Reproducible R analysis framework

Research study: Association Between Systemic Immune-Inflammation Index and Three-Year Visual Decline in Patients with Diabetic Retinopathy.

## What this repository is
A cleaned, modular **reconstruction** of the older exploratory analysis script. **Not** an assurance that newly recomputed coefficients equal the submitted manuscript. Outputs that do not reproduce manuscript results require resolution before treating the repository as the paper's exact source code.

## Setup
1. Install `readxl`, `dplyr`, `tidyr`, `mice`, `survival`, `broom`, `ggplot2`, `rms`.
2. Put the institution-approved private Excel dataset at `data/final_dataset_3y_logmar2_1.xlsx` (NEVER commit it).
3. From repository root, use `source('run_all.R')`.
4. To estimate ΔC-index, run `source('scripts/07_delta_cindex_bootstrap.R')`. Set `Sys.setenv(SII_BOOTSTRAP_B='1000')` to control resamples.

## Module organization
| R file | Responsibility |
|---|---|
| `00_config.R` | Paths, models, helpers |
| `01_prepare_data.R` | Read Excel, normalize variables, stop if cohort totals differ |
| `02_multiple_imputation.R` | Create 20 MICE completed datasets (new controlled run) |
| `03_cox_models.R` | M1–M4 cluster-robust Cox, pooled MI Wald D1 |
| `04_ph_and_rcs.R` | PH diagnostic and 4-knot RCS, on first completion only |
| `05_subgroups_S4.R` | Eight subgroup analyses, observed counts and imputed-stratum sizes |
| `06_sensitivity_S5.R` | Recomputable sensitivities; identifies missing source inputs |
| `07_delta_cindex_bootstrap.R` | Patient-level paired bootstrap; average across completed MI data |
| `08_km_overall.R` | 1/2/3-year overall KM from recorded time-to-event data |
| `09_imputation_diagnostics.R` | MICE convergence and observed/imputed density plots |

## Mandatory study-specific audit before public claims
- **Original script has m=10 and m=20 calls**. Author confirms all reported analyses used 20, but the exact original imputation specifications must be recovered to claim exact reproducibility.
- **Eye-specific SII versus patient-level imputation**: the manuscript says patient-level systemic covariates were imputed at patient level. This reconstruction currently uses eye-record-level MICE to preserve differing eye-index lab values. Confirm the original two-level method and amend either code or manuscript to what was actually done.
- **RCS**: old exploratory code used 3 and 4 knots. This cleaned code uses 4 knots at 5/35/65/95 percentiles, matching manuscript; reported spline statistics must be validated against original.
- **HbA1c subgroups**: older code once used a median split. This project uses `<7%` vs `>=7%` as published.
- **Subgroup counts**: when HbA1c/eGFR values are imputed, their strata can vary across completed datasets. Report MI counts with range or clearly mark original observed counts; do not silently equate them.
- **Bootstrap**: original code resampled eyes and used one imputation. Clean script resamples patient clusters and uses all 20 completed datasets, but is conditional on imputations (does not redo MICE inside each resample). Its CI is methodologically **not identical** to a fully nested bootstrap-plus-MICE CI. Do not assert that this script yielded any existing manuscript number unless its output is independently verified.
- **Sensitivity analyses**: baseline BCVA, lab dates for 30-day window, and visit-level repeat-confirmation are unavailable in the present workbook; code does not fabricate them. For the >=0.3 endpoint the `days` column is reused; verify if the real outcome has a different event time before running as paper-level analysis.
- **PH diagnostic and RCS** are currently on the first completed dataset, not pooled; do not describe these as MI-pooled analyses. Fill in exactly the procedure actually used in the final manuscript.
- Eye-level records are private and excluded via `.gitignore`. Remove identifiable data, paths, institution-specific raw values from output before posting to GitHub.
