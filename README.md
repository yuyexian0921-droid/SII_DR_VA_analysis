# SII_DR_VA_analysis

R code for the statistical analyses evaluating the association between the **systemic immune-inflammation index (SII)** and **3-year visual decline** in patients with diabetic retinopathy.

## Repository structure

```text
SII_DR_VA_analysis/
├── 01_data_preparation.R
├── 02_descriptive_analysis.R
├── 03_multiple_imputation.R
├── 04_primary_cox_models.R
├── 05_model_performance.R
├── 06_nonlinear_analysis.R
├── 07_subgroup_analysis.R
├── 08_sensitivity_analysis.R
├── .gitignore
└── README.md
```

## Analysis workflow

| Script | Purpose |
|---|---|
| `01_data_preparation.R` | Imports the analytic dataset, derives age, AAR, SII/100, and comorbidity variables, and prepares the final analysis dataset. |
| `02_descriptive_analysis.R` | Assesses missingness and generates baseline descriptive statistics. |
| `03_multiple_imputation.R` | Performs multiple imputation using `mice` (10 imputations, 20 iterations). |
| `04_primary_cox_models.R` | Runs univariable Cox regression, sequential multivariable Cox models, model-comparison tests, and PH-assumption checks. |
| `05_model_performance.R` | Compares Harrell's C-index for the clinical model with and without SII and bootstraps the C-index difference. |
| `06_nonlinear_analysis.R` | Performs SII tertile and restricted cubic spline analyses. |
| `07_subgroup_analysis.R` | Performs prespecified subgroup and interaction analyses. |
| `08_sensitivity_analysis.R` | Performs sensitivity analyses under alternative restrictions and outcome definitions. |

## Main modeling convention

SII is rescaled as:

```r
SII100 = SII / 100
```

Therefore, hazard ratios for `SII100` represent the association with a **100-unit increase in SII**.

Analyses are conducted at the eye level. Robust standard errors in the Cox models account for within-patient correlation using:

```r
cluster(ID)
```

## Data

The patient-level analytic dataset is **not included in this public repository**.

By default, `01_data_preparation.R` expects:

```text
data/final_dataset_3y_logmar2_1.xlsx
```

You can instead point to a local file by setting the `SII_DATA_FILE` environment variable before running the scripts.

Example:

```r
Sys.setenv(SII_DATA_FILE = "D:/your_private_folder/final_dataset_3y_logmar2_1.xlsx")
source("01_data_preparation.R")
```

## Required R packages

The analysis uses the following packages:

```r
readxl
dplyr
tidyr
lubridate
tableone
mice
survival
broom
rms
ggplot2
Hmisc
boot
```

Install missing packages with, for example:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "lubridate", "tableone",
  "mice", "survival", "broom", "rms", "ggplot2", "Hmisc", "boot"
))
```

## Running the analysis

Run the scripts from the repository root in numerical order:

```r
source("01_data_preparation.R")
source("02_descriptive_analysis.R")
source("03_multiple_imputation.R")
source("04_primary_cox_models.R")
source("05_model_performance.R")
source("06_nonlinear_analysis.R")
source("07_subgroup_analysis.R")
source("08_sensitivity_analysis.R")
```

The later scripts contain checks that automatically source the required earlier script if its main analysis object is not already available.

## Reproducibility note

Random seeds are specified for multiple imputation and bootstrap analyses. Results can still vary across R/package versions, so reporting the R session information used for the final manuscript is recommended.
