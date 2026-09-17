# ==============================================================================
# 06_nonlinear_analysis.R
# SII tertile and restricted cubic spline analyses
# ==============================================================================

if (!exists("imp")) {
  source("03_multiple_imputation.R")
}

required_packages <- c("dplyr", "mice", "survival", "rms", "ggplot2")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(dplyr)
library(mice)
library(survival)
library(rms)
library(ggplot2)

# ------------------------------------------------------------------------------
# SII tertiles
# ------------------------------------------------------------------------------

sii_breaks <- quantile(
  analysis_data$SII100,
  probs = c(0, 1 / 3, 2 / 3, 1),
  na.rm = TRUE,
  names = FALSE
)

sii_breaks <- unique(sii_breaks)

if (length(sii_breaks) == 4) {
  analysis_data_tertile <- analysis_data %>%
    mutate(
      SII_tertile = cut(
        SII100,
        breaks = sii_breaks,
        include.lowest = TRUE,
        labels = c("Low", "Medium", "High")
      ),
      SII_tertile = factor(
        SII_tertile,
        levels = c("Low", "Medium", "High")
      )
    )

  imp_tertile <- mice(
    analysis_data_tertile,
    m = 10,
    maxit = 20,
    seed = 2026,
    printFlag = FALSE
  )

  sii_tertile_model <- with(
    imp_tertile,
    coxph(
      Surv(days, outcome_3y_logmar2) ~
        age + sex + HbA1c + eGFR +
        DR_grading + DME + CAT + GLAU +
        SII_tertile + cluster(ID)
    )
  )

  sii_tertile_results <- summary(
    pool(sii_tertile_model),
    conf.int = TRUE,
    exponentiate = TRUE
  )

  print(sii_tertile_results)
} else {
  warning("SII tertiles could not be created because quantile cut points were duplicated.")
}

# ------------------------------------------------------------------------------
# Restricted cubic spline (RCS)
# ------------------------------------------------------------------------------
# RCS is fitted on complete cases for variables required in this model.

rcs_data <- analysis_data %>%
  select(
    ID,
    days,
    outcome_3y_logmar2,
    eGFR,
    sex,
    HbA1c,
    age,
    SII100,
    DR_grading,
    DME,
    CAT,
    GLAU
  ) %>%
  filter(complete.cases(.))

rcs_dd <- datadist(rcs_data)
rcs_dd$limits["Adjust to", "SII100"] <- median(
  rcs_data$SII100,
  na.rm = TRUE
)
options(datadist = "rcs_dd")

rcs_model <- cph(
  Surv(days, outcome_3y_logmar2) ~
    eGFR +
    sex +
    HbA1c +
    age +
    rcs(SII100, 4) +
    DR_grading +
    DME +
    CAT +
    GLAU +
    cluster(ID),
  data = rcs_data,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

rcs_pred <- as.data.frame(
  Predict(
    rcs_model,
    SII100,
    fun = exp,
    ref.zero = TRUE,
    conf.int = 0.95
  )
)

rcs_plot <- ggplot(rcs_pred, aes(x = SII100, y = yhat)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "grey80",
    alpha = 0.8
  ) +
  geom_line(
    linewidth = 1.1,
    color = "#2C7FB8"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    limits = c(0, 22),
    breaks = seq(0, 20, 5)
  ) +
  scale_y_continuous(
    limits = c(0.8, 2.0),
    breaks = seq(0.8, 2.0, 0.2)
  ) +
  labs(
    x = "SII (per 100-unit increase)",
    y = "Hazard ratio (95% CI)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

print(rcs_plot)
print(anova(rcs_model))
