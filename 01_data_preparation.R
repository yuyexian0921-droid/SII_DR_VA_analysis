# ==============================================================================
# 01_data_preparation.R
# Data import and preparation
# ==============================================================================

required_packages <- c("readxl", "dplyr", "lubridate")

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}))

library(readxl)
library(dplyr)
library(lubridate)

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------
# Keep patient-level data outside the public repository. By default, the script
# looks for the file below. Alternatively, set the environment variable
# SII_DATA_FILE to a local file path before running the analysis.

DATA_FILE <- Sys.getenv(
  "SII_DATA_FILE",
  unset = "data/final_dataset_3y_logmar2_1.xlsx"
)
DATA_SHEET <- 1

if (!file.exists(DATA_FILE)) {
  stop(
    paste0(
      "Data file not found: ", DATA_FILE, "\n",
      "Place the analytic dataset in data/ or set SII_DATA_FILE to its path."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------

data_raw <- read_excel(
  path = DATA_FILE,
  sheet = DATA_SHEET
)

cat("Raw data dimensions:", nrow(data_raw), "rows x", ncol(data_raw), "columns\n")

# ------------------------------------------------------------------------------
# Derived variables and recoding
# ------------------------------------------------------------------------------

data <- data_raw %>%
  mutate(
    index_date = as.Date(index_date),
    index_year = year(index_date),
    age = index_year - birth_year,

    # ALT-to-AST ratio
    AAR = if_else(!is.na(AST) & AST != 0, ALT / AST, NA_real_),

    # SII modeled per 100-unit increase
    SII100 = SII / 100,

    # Comorbidity burden
    comorbidity_score = HTN + HLD + CAD + CeVD + DKD + DFU,
    comorbidity_group = case_when(
      comorbidity_score == 0 ~ "0",
      comorbidity_score == 1 ~ "1",
      comorbidity_score >= 2 ~ ">=2",
      TRUE ~ NA_character_
    ),

    # Sex recoding
    sex = recode(
      as.character(sex),
      "男" = "M",
      "女" = "F",
      .default = as.character(sex)
    ),

    # Factor variables
    sex = factor(sex),
    DR_grading = factor(DR_grading),
    DME = factor(DME),
    anti_VEGF = factor(anti_VEGF),
    PRP = factor(PRP),
    PPV = factor(PPV),
    CAT = factor(CAT),
    GLAU = factor(GLAU),
    ERM = factor(ERM),
    HTN = factor(HTN),
    comorbidity_group = factor(
      comorbidity_group,
      levels = c("0", "1", ">=2")
    )
  )

cat("Number of eyes:", nrow(data), "\n")
cat("Number of patients:", n_distinct(data$ID), "\n")

# ------------------------------------------------------------------------------
# Final analysis dataset
# ------------------------------------------------------------------------------

analysis_data <- data %>%
  select(
    ID,
    days,
    outcome_3y_logmar2,
    any_of("outcome_logmar3"),

    age,
    sex,

    HbA1c,
    DM_duration,
    eGFR,
    SII,
    SII100,
    ALT,
    AST,
    AAR,
    Hb,
    RDW_CV,
    GGT,
    comorbidity_score,
    comorbidity_group,
    HTN,

    DR_grading,
    DME,
    IOP,

    anti_VEGF,
    PRP,
    PPV,

    CAT,
    GLAU,
    ERM
  )

cat("Analysis dataset prepared successfully.\n")
