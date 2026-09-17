# ==============================================================================
# 03_multiple_imputation.R
# Multiple imputation of missing covariate data
# ==============================================================================

if (!exists("analysis_data")) {
  source("01_data_preparation.R")
}

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Package 'mice' is required but not installed.", call. = FALSE)
}

library(mice)

set.seed(20260915)

# Ten imputed datasets and 20 iterations, matching the main analysis.
imp <- mice(
  analysis_data,
  m = 10,
  maxit = 20,
  seed = 20260915,
  printFlag = FALSE
)

print(imp)
cat("Multiple imputation completed.\n")
