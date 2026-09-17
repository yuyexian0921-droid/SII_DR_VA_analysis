# ==============================================================================
# 07_subgroup_analysis.R
# Subgroup and interaction analyses
# ==============================================================================

if (!exists("imp")) {
  source("03_multiple_imputation.R")
}

required_packages <- c("mice", "survival", "dplyr")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(mice)
library(survival)
library(dplyr)

# ------------------------------------------------------------------------------
# Create subgroup variables in all imputed datasets
# ------------------------------------------------------------------------------

imp_long <- complete(
  imp,
  action = "long",
  include = TRUE
) %>%
  group_by(.imp) %>%
  mutate(
    Age_group = if_else(age < 60, "<60", ">=60"),
    HbA1c_group = if_else(
      HbA1c < median(HbA1c, na.rm = TRUE),
      "Low",
      "High"
    ),
    eGFR_group = if_else(eGFR < 60, "<60", ">=60"),
    DR_stage = if_else(DR_grading == "PDR", "PDR", "NPDR")
  ) %>%
  ungroup()

imp_sub <- as.mids(imp_long)

base_covariates <- c(
  "age",
  "sex",
  "HbA1c",
  "eGFR",
  "DR_grading",
  "DME",
  "CAT",
  "GLAU"
)

# ------------------------------------------------------------------------------
# Subgroup-specific effect of SII
# ------------------------------------------------------------------------------
# A factor used to define a subgroup is removed from the adjustment set if it is
# constant within that subgroup (e.g., sex within male/female strata).

run_subgroup <- function(
    imp_data,
    subgroup_var,
    subgroup_value,
    remove_covariates = character(0)) {

  completed_sets <- lapply(seq_len(imp_data$m), function(i) {
    d <- complete(imp_data, i)
    d[d[[subgroup_var]] == subgroup_value, , drop = FALSE]
  })

  covars <- setdiff(base_covariates, remove_covariates)
  rhs <- paste(c(covars, "SII100", "cluster(ID)"), collapse = " + ")
  formula <- as.formula(
    paste0("Surv(days, outcome_3y_logmar2) ~ ", rhs)
  )

  fits <- lapply(completed_sets, function(d) {
    coxph(formula, data = d)
  })

  pooled <- pool(as.mira(fits))

  summary(
    pooled,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term == "SII100")
}

subgroup_spec <- list(
  DR_NPDR = list(var = "DR_stage", value = "NPDR", remove = "DR_grading"),
  DR_PDR = list(var = "DR_stage", value = "PDR", remove = "DR_grading"),
  DME_0 = list(var = "DME", value = "0", remove = "DME"),
  DME_1 = list(var = "DME", value = "1", remove = "DME"),
  CAT_0 = list(var = "CAT", value = "0", remove = "CAT"),
  CAT_1 = list(var = "CAT", value = "1", remove = "CAT"),
  GLAU_0 = list(var = "GLAU", value = "0", remove = "GLAU"),
  GLAU_1 = list(var = "GLAU", value = "1", remove = "GLAU"),
  Sex_F = list(var = "sex", value = "F", remove = "sex"),
  Sex_M = list(var = "sex", value = "M", remove = "sex"),
  Age_low = list(var = "Age_group", value = "<60", remove = character(0)),
  Age_high = list(var = "Age_group", value = ">=60", remove = character(0)),
  HbA1c_low = list(var = "HbA1c_group", value = "Low", remove = character(0)),
  HbA1c_high = list(var = "HbA1c_group", value = "High", remove = character(0)),
  eGFR_low = list(var = "eGFR_group", value = "<60", remove = character(0)),
  eGFR_high = list(var = "eGFR_group", value = ">=60", remove = character(0))
)

subgroup_results <- lapply(names(subgroup_spec), function(name) {
  spec <- subgroup_spec[[name]]

  run_subgroup(
    imp_data = imp_sub,
    subgroup_var = spec$var,
    subgroup_value = spec$value,
    remove_covariates = spec$remove
  ) %>%
    mutate(Subgroup = name, .before = 1)
})

subgroup_table <- bind_rows(subgroup_results) %>%
  select(Subgroup, estimate, conf.low, conf.high, p.value) %>%
  rename(
    HR = estimate,
    Lower95 = conf.low,
    Upper95 = conf.high,
    P = p.value
  )

print(subgroup_table)

# ------------------------------------------------------------------------------
# Interaction analyses
# ------------------------------------------------------------------------------

run_interaction_model <- function(imp_data, interaction_var, covariates) {
  rhs <- paste(
    c(covariates, paste0("SII100 * ", interaction_var), "cluster(ID)"),
    collapse = " + "
  )

  formula <- as.formula(
    paste0("Surv(days, outcome_3y_logmar2) ~ ", rhs)
  )

  fit <- with(
    imp_data,
    coxph(formula)
  )

  summary(
    pool(fit),
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(grepl("SII100:", term))
}

interaction_definitions <- list(
  `DR severity` = list(
    var = "DR_stage",
    covars = setdiff(base_covariates, "DR_grading")
  ),
  DME = list(
    var = "DME",
    covars = setdiff(base_covariates, "DME")
  ),
  Cataract = list(
    var = "CAT",
    covars = setdiff(base_covariates, "CAT")
  ),
  Glaucoma = list(
    var = "GLAU",
    covars = setdiff(base_covariates, "GLAU")
  ),
  Sex = list(
    var = "sex",
    covars = setdiff(base_covariates, "sex")
  ),
  Age = list(
    var = "Age_group",
    covars = setdiff(base_covariates, "age")
  ),
  HbA1c = list(
    var = "HbA1c_group",
    covars = setdiff(base_covariates, "HbA1c")
  ),
  eGFR = list(
    var = "eGFR_group",
    covars = setdiff(base_covariates, "eGFR")
  )
)

interaction_table_final <- bind_rows(
  lapply(names(interaction_definitions), function(label) {
    definition <- interaction_definitions[[label]]

    run_interaction_model(
      imp_data = imp_sub,
      interaction_var = definition$var,
      covariates = definition$covars
    ) %>%
      mutate(Subgroup = label, .before = 1)
  })
) %>%
  select(Subgroup, term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    Interaction = term,
    HR = estimate,
    Lower95 = conf.low,
    Upper95 = conf.high,
    P_interaction = p.value
  )

print(interaction_table_final)
