# ==============================================================================
# 02_descriptive_analysis.R
# Missing-data assessment and baseline characteristics
# ==============================================================================

if (!exists("analysis_data")) {
  source("01_data_preparation.R")
}

required_packages <- c("dplyr", "tidyr", "tableone")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(dplyr)
library(tidyr)
library(tableone)

# ------------------------------------------------------------------------------
# Missing-data assessment
# ------------------------------------------------------------------------------

missing_table <- analysis_data %>%
  summarise(
    across(
      everything(),
      ~ mean(is.na(.x)) * 100
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_percent"
  ) %>%
  arrange(desc(missing_percent))

print(missing_table)

# ------------------------------------------------------------------------------
# Baseline characteristics (Table 1)
# ------------------------------------------------------------------------------

baseline_vars <- c(
  "age",
  "sex",
  "HbA1c",
  "DM_duration",
  "eGFR",
  "SII100",
  "AAR",
  "Hb",
  "RDW_CV",
  "GGT",
  "comorbidity_score",
  "HTN",
  "DR_grading",
  "DME",
  "IOP",
  "CAT",
  "GLAU",
  "ERM"
)

baseline_vars <- intersect(baseline_vars, names(analysis_data))

table1 <- CreateTableOne(
  vars = baseline_vars,
  strata = "outcome_3y_logmar2",
  data = analysis_data,
  test = TRUE
)

print(
  table1,
  showAllLevels = TRUE
)
