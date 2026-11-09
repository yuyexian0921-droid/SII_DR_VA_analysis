source('scripts/00_config.R')
d <- readRDS(file.path(OUT,'data_prepared.rds'))
# The original large script has conflicting m=10/m=20 calls.
# This module specifies ONE controlled m=20 run.
mi_vars <- c('ID','days','outcome_3y_logmar2','age','sex','HbA1c','eGFR',
             'DR_grading','DME','CAT','GLAU','SII100')
x <- d[mi_vars]
init <- mice::mice(x,maxit=0,printFlag=FALSE)
methods <- init$method; methods[] <- ''
for(v in c('HbA1c','eGFR','SII100')) methods[v] <- 'pmm'
pred <- init$predictorMatrix;pred[,] <- 0
# Use all relevant observed covariates and outcome for missing variables;
# patient identifier must NOT be treated as a continuous predictor.
for(v in c('HbA1c','eGFR','SII100')) pred[v,setdiff(mi_vars,c('ID',v))] <- 1
# Imputed predictor models include outcome and follow-up, but ID/eye is excluded.
# This is an eye-record-level MICE implementation. Verify against the actual
# patient-level imputation used for published estimates before claiming exact reproduction.
imp <- mice::mice(x,m=20,maxit=20,method=methods,predictorMatrix=pred,
                  seed=20260915,printFlag=FALSE)
saveRDS(imp,file.path(OUT,'imputations.rds'))
# Keep eye-specific rows and event data stable across all completions.
stopifnot(imp$m==20L)
message('20 MICE datasets saved. Verify original imputation-level methodology before publication.')
