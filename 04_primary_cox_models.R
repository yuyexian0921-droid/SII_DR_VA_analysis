
# ============================================================
# 04 Primary Cox regression models
# ============================================================
# Eye-level Cox models with robust variance estimation.
# cluster(ID) accounts for correlation between two eyes.
#
# Model 4:
# age + sex + HbA1c + eGFR +
# DR severity + DME + CAT + GLAU + SII100
# ============================================================

library(survival)
library(mice)
library(broom)

fit_model4 <- with(
  imp,
  coxph(
    Surv(days, outcome_3y_logmar2) ~
      age +
      sex +
      HbA1c +
      eGFR +
      DR_grading +
      DME +
      CAT +
      GLAU +
      SII100 +
      cluster(ID)
  )
)

model4_result <- summary(
  pool(fit_model4),
  conf.int = TRUE,
  exponentiate = TRUE
)

model4_result
