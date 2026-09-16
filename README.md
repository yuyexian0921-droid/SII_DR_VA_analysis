
# SII_DR_VA_analysis_code

Analysis workflow for:
Systemic Immune-Inflammation Index and Risk of Three-Year Visual Decline in Patients with Diabetic Retinopathy:
A Real-World Longitudinal Cohort Study

Workflow:
01_data_preparation.R
02_descriptive_analysis.R
03_multiple_imputation.R
04_primary_cox_models.R
05_model_performance.R
06_nonlinear_analysis.R
07_subgroup_analysis.R
08_sensitivity_analysis.R

Notes:
- Patient-level data are not included.
- The unit of analysis is eye-level.
- Correlation between two eyes from the same patient is handled using cluster(ID).
- Missing covariates are handled using multiple imputation by chained equations (mice).
- Results are pooled using Rubin's rules.
