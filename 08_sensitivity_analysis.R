# ==============================================================================
# 08_sensitivity_analysis.R
# Sensitivity analyses
# ==============================================================================

if (!exists("imp")) {
  source("03_multiple_imputation.R")
}

required_packages <- c("mice", "survival", "dplyr", "broom")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(mice)
library(survival)
library(dplyr)
library(broom)

# ------------------------------------------------------------------------------
# Convert imputed datasets to a list
# ------------------------------------------------------------------------------

imp_list <- lapply(
  seq_len(imp$m),
  function(i) complete(imp, i)
)

base_covariates <- c(
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

build_cox_formula <- function(
    outcome = "outcome_3y_logmar2",
    remove_covariates = character(0)) {

  covars <- setdiff(base_covariates, remove_covariates)

  as.formula(
    paste0(
      "Surv(days, ", outcome, ") ~ ",
      paste(c(covars, "cluster(ID)"), collapse = " + ")
    )
  )
}

run_sensitivity_cox <- function(
    data_list,
    outcome = "outcome_3y_logmar2",
    remove_covariates = character(0)) {

  formula <- build_cox_formula(
    outcome = outcome,
    remove_covariates = remove_covariates
  )

  fits <- lapply(data_list, function(d) {
    coxph(formula, data = d)
  })

  pooled <- pool(as.mira(fits))

  summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term == "SII100") %>%
    transmute(
      HR = estimate,
      Lower95 = conf.low,
      Upper95 = conf.high,
      P = p.value,
      HR_CI = sprintf("%.3f (%.3f-%.3f)", estimate, conf.low, conf.high),
      P_value = format.pval(p.value, digits = 3, eps = 0.001)
    )
}

# ------------------------------------------------------------------------------
# Sensitivity subsets
# ------------------------------------------------------------------------------

DME_list <- lapply(
  imp_list,
  function(d) d[as.character(d$DME) == "0", , drop = FALSE]
)

CAT_list <- lapply(
  imp_list,
  function(d) d[as.character(d$CAT) == "0", , drop = FALSE]
)

GLAU_list <- lapply(
  imp_list,
  function(d) d[as.character(d$GLAU) == "0", , drop = FALSE]
)

NPDR_list <- lapply(
  imp_list,
  function(d) d[as.character(d$DR_grading) == "NPDR", , drop = FALSE]
)

FU180_list <- lapply(
  imp_list,
  function(d) d[d$days >= 180, , drop = FALSE]
)

# ------------------------------------------------------------------------------
# Main and sensitivity analyses
# ------------------------------------------------------------------------------

main_result <- run_sensitivity_cox(imp_list)

DME_result <- run_sensitivity_cox(
  DME_list,
  remove_covariates = "DME"
)

CAT_result <- run_sensitivity_cox(
  CAT_list,
  remove_covariates = "CAT"
)

GLAU_result <- run_sensitivity_cox(
  GLAU_list,
  remove_covariates = "GLAU"
)

FU180_result <- run_sensitivity_cox(FU180_list)

NPDR_result <- run_sensitivity_cox(
  NPDR_list,
  remove_covariates = "DR_grading"
)

# Stricter visual-decline definition, if present in the dataset
logmar03_result <- NULL
if ("outcome_logmar3" %in% names(analysis_data)) {
  logmar03_result <- run_sensitivity_cox(
    imp_list,
    outcome = "outcome_logmar3"
  )
}

# ------------------------------------------------------------------------------
# Complete-case analysis
# ------------------------------------------------------------------------------

complete_case_data <- analysis_data %>%
  filter(
    complete.cases(
      days,
      outcome_3y_logmar2,
      age,
      sex,
      HbA1c,
      eGFR,
      DR_grading,
      DME,
      CAT,
      GLAU,
      SII100
    )
  )

complete_formula <- build_cox_formula()

complete_result <- broom::tidy(
  coxph(complete_formula, data = complete_case_data),
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term == "SII100") %>%
  transmute(
    HR = estimate,
    Lower95 = conf.low,
    Upper95 = conf.high,
    P = p.value,
    HR_CI = sprintf("%.3f (%.3f-%.3f)", estimate, conf.low, conf.high),
    P_value = format.pval(p.value, digits = 3, eps = 0.001)
  )

# ------------------------------------------------------------------------------
# Final sensitivity table
# ------------------------------------------------------------------------------

Sensitivity_table <- bind_rows(
  data.frame(Analysis = "Main analysis", main_result),
  data.frame(Analysis = "Excluding baseline DME", DME_result),
  data.frame(Analysis = "Excluding baseline cataract", CAT_result),
  data.frame(Analysis = "Excluding baseline glaucoma", GLAU_result),
  data.frame(Analysis = "Follow-up >=180 days", FU180_result),
  data.frame(Analysis = "Complete-case analysis", complete_result),
  data.frame(Analysis = "Restricting to NPDR eyes", NPDR_result)
)

if (!is.null(logmar03_result)) {
  Sensitivity_table <- bind_rows(
    Sensitivity_table,
    data.frame(
      Analysis = "Visual decline defined as logMAR increase >=0.3",
      logmar03_result
    )
  )
}

print(Sensitivity_table)
