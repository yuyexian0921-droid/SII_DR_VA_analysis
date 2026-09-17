# ==============================================================================
# 05_model_performance.R
# C-index comparison of the clinical model with and without SII
# ==============================================================================

if (!exists("imp")) {
  source("03_multiple_imputation.R")
}

required_packages <- c("survival", "mice", "Hmisc", "boot")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(survival)
library(mice)
library(Hmisc)
library(boot)

complete_data_1 <- complete(imp, 1)

# ------------------------------------------------------------------------------
# C-index for Model 3 and Model 4
# ------------------------------------------------------------------------------

cox3_complete <- coxph(
  Surv(days, outcome_3y_logmar2) ~
    age + sex + HbA1c + eGFR +
    DR_grading + DME + CAT + GLAU,
  data = complete_data_1,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

cox4_complete <- coxph(
  Surv(days, outcome_3y_logmar2) ~
    age + sex + HbA1c + eGFR +
    DR_grading + DME + CAT + GLAU + SII100,
  data = complete_data_1,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

c_index_model3 <- Hmisc::rcorr.cens(
  -predict(cox3_complete),
  cox3_complete$y
)["C Index"]

c_index_model4 <- Hmisc::rcorr.cens(
  -predict(cox4_complete),
  cox4_complete$y
)["C Index"]

cat("Model 3 C-index:", c_index_model3, "\n")
cat("Model 4 C-index:", c_index_model4, "\n")
cat("Difference:", c_index_model4 - c_index_model3, "\n")

# ------------------------------------------------------------------------------
# Bootstrap 95% CI for the C-index difference
# ------------------------------------------------------------------------------

cindex_diff <- function(data, index) {
  d <- data[index, , drop = FALSE]

  fit3 <- coxph(
    Surv(days, outcome_3y_logmar2) ~
      age + sex + HbA1c + eGFR +
      DR_grading + DME + CAT + GLAU,
    data = d,
    x = TRUE,
    y = TRUE
  )

  fit4 <- coxph(
    Surv(days, outcome_3y_logmar2) ~
      age + sex + HbA1c + eGFR +
      DR_grading + DME + CAT + GLAU + SII100,
    data = d,
    x = TRUE,
    y = TRUE
  )

  c3 <- Hmisc::rcorr.cens(-predict(fit3), fit3$y)["C Index"]
  c4 <- Hmisc::rcorr.cens(-predict(fit4), fit4$y)["C Index"]

  as.numeric(c4 - c3)
}

set.seed(2026)
boot_result <- boot::boot(
  data = complete_data_1,
  statistic = cindex_diff,
  R = 1000
)

cindex_difference_ci <- boot::boot.ci(
  boot_result,
  type = "perc"
)

print(cindex_difference_ci)
