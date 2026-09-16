
# ============================================================
# 03 Multiple imputation
# ============================================================
# Missing values are assumed to be missing at random.
# Five imputed datasets are generated using MICE.
# ============================================================

library(mice)

imp_vars <- c(
  "days",
  "outcome_3y_logmar2",
  "outcome_logmar3",
  "age",
  "sex",
  "HbA1c",
  "eGFR",
  "DR_grading",
  "DME",
  "CAT",
  "GLAU",
  "SII100",
  "ID"
)

imp <- mice(
  data[, imp_vars],
  m = 5,
  seed = 2026,
  maxit = 20
)

summary(imp)
