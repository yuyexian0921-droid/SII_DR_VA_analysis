
# ============================================================
# 02 Descriptive analysis
# ============================================================
# Generate baseline characteristics and missing data summary.
# ============================================================

library(tableone)
library(mice)

# Missing data pattern
md.pattern(data)

# Example Table 1
vars <- c(
  "age","sex","HbA1c","eGFR",
  "DR_grading","DME","CAT","GLAU","SII100"
)

CreateTableOne(
  vars = vars,
  strata = "outcome_3y_logmar2",
  data = data,
  factorVars = c("sex","DR_grading","DME","CAT","GLAU")
)
