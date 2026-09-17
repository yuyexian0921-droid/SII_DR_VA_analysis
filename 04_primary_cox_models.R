# ==============================================================================
# 04_primary_cox_models.R
# Univariable and sequential multivariable Cox regression models
# ==============================================================================

if (!exists("imp")) {
  source("03_multiple_imputation.R")
}

required_packages <- c("survival", "mice", "dplyr", "broom")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(survival)
library(mice)
library(dplyr)
library(broom)

# ------------------------------------------------------------------------------
# Univariable Cox regression
# ------------------------------------------------------------------------------

univariable_vars <- c(
  "age",
  "sex",
  "HbA1c",
  "eGFR",
  "DR_grading",
  "DME",
  "CAT",
  "GLAU",
  "SII100"
)

run_univariable_cox <- function(var, data = analysis_data) {
  formula <- as.formula(
    paste0("Surv(days, outcome_3y_logmar2) ~ ", var)
  )

  fit <- coxph(
    formula = formula,
    data = data,
    cluster = ID
  )

  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(Variable = var, .before = 1)
}

univariable_results <- bind_rows(
  lapply(univariable_vars, run_univariable_cox)
) %>%
  mutate(
    HR_CI = sprintf("%.3f (%.3f-%.3f)", estimate, conf.low, conf.high),
    P = format.pval(p.value, digits = 3, eps = 0.001)
  )

print(univariable_results)

# ------------------------------------------------------------------------------
# Sequential multivariable Cox models
# ------------------------------------------------------------------------------

model1 <- with(
  imp,
  coxph(
    Surv(days, outcome_3y_logmar2) ~
      age +
      sex +
      cluster(ID)
  )
)

model2 <- with(
  imp,
  coxph(
    Surv(days, outcome_3y_logmar2) ~
      age +
      sex +
      HbA1c +
      eGFR +
      cluster(ID)
  )
)

model3 <- with(
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
      cluster(ID)
  )
)

model4 <- with(
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

model1_pool <- pool(model1)
model2_pool <- pool(model2)
model3_pool <- pool(model3)
model4_pool <- pool(model4)

extract_pooled_hr <- function(pooled_model, model_name) {
  summary(
    pooled_model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    mutate(Model = model_name, .before = 1)
}

HR_table <- bind_rows(
  extract_pooled_hr(model1_pool, "Model 1"),
  extract_pooled_hr(model2_pool, "Model 2"),
  extract_pooled_hr(model3_pool, "Model 3"),
  extract_pooled_hr(model4_pool, "Model 4")
) %>%
  select(Model, term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    P = p.value
  )

print(HR_table)

# ------------------------------------------------------------------------------
# Nested model comparison
# ------------------------------------------------------------------------------

comparison_12 <- D1(model2, model1)
comparison_23 <- D1(model3, model2)
comparison_34 <- D1(model4, model3)

print(comparison_12)
print(comparison_23)
print(comparison_34)

# ------------------------------------------------------------------------------
# Proportional-hazards assumption
# ------------------------------------------------------------------------------
# cox.zph() is applied to one completed dataset because it cannot be directly
# pooled across multiple imputations.

complete_data_1 <- complete(imp, 1)

cox_model4_complete <- coxph(
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
    cluster(ID),
  data = complete_data_1,
  x = TRUE,
  y = TRUE
)

ph_test <- cox.zph(cox_model4_complete)
print(ph_test)

# Optional PH diagnostic plots
ph_terms <- c(
  "age", "sex", "HbA1c", "eGFR", "DR_grading",
  "DME", "CAT", "GLAU", "SII100"
)

ph_titles <- c(
  "Age", "Sex", "HbA1c", "eGFR", "DR severity",
  "DME", "Cataract", "Glaucoma", "SII (per 100 units)"
)

available_terms <- intersect(ph_terms, rownames(ph_test$table))

if (interactive() && length(available_terms) > 0) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

  for (v in available_terms) {
    idx <- match(v, ph_terms)
    plot(
      ph_test[v],
      main = ph_titles[idx],
      xlab = "Time",
      ylab = "Scaled Schoenfeld residuals"
    )
    abline(h = 0, lty = 2)
  }
}
