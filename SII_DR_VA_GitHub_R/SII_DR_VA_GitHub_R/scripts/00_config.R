# Project-wide configuration. Run from repository root.
required <- c('readxl','dplyr','tidyr','mice','survival','broom','ggplot2','rms')
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_pkgs)) stop('Install packages: ',paste(missing_pkgs,collapse=', '))
INPUT <- Sys.getenv('SII_INPUT', unset='data/final_dataset_3y_logmar2_1.xlsx')
OUT <- 'output'
dir.create(OUT, showWarnings=FALSE,recursive=TRUE)
set.seed(20260915)
MODEL3 <- 'Surv(days, outcome_3y_logmar2) ~ age + sex + HbA1c + eGFR + DR_grading + DME + CAT + GLAU + cluster(ID)'
MODEL4 <- 'Surv(days, outcome_3y_logmar2) ~ age + sex + HbA1c + eGFR + DR_grading + DME + CAT + GLAU + SII100 + cluster(ID)'
compact <- function(x) ifelse(is.na(x),NA_character_,trimws(as.character(x)))
get_completed <- function() { x<-readRDS(file.path(OUT,'imputations.rds'));lapply(seq_len(x$m),function(i) mice::complete(x,i)) }
model_formula <- function(x) stats::as.formula(x)
write_csv <- function(x,name) utils::write.csv(x,file.path(OUT,name),row.names=FALSE,fileEncoding='UTF-8')
