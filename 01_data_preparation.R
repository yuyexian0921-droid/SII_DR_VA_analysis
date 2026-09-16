
# ============================================================
# 01 Data preparation
# ============================================================
# Purpose:
# Import raw analytic dataset, define variables, and prepare
# the dataset for statistical analyses.
#
# Required variables:
# ID, days, outcome_3y_logmar2, outcome_logmar3,
# age, sex, HbA1c, eGFR, DR_grading, DME, CAT, GLAU, SII100
# ============================================================

library(readxl)
library(dplyr)

data <- read_excel(
  "D:/phd/SII_DR_VA/final_dataset_3y_logmar2_1.xlsx"
)

# Convert categorical variables to factors
data <- data %>%
  mutate(
    sex = factor(sex),
    DR_grading = factor(DR_grading),
    DME = factor(DME),
    CAT = factor(CAT),
    GLAU = factor(GLAU)
  )

# Check data structure
str(data)
summary(data)
