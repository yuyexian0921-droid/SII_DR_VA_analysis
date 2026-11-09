source('scripts/00_config.R')
imp <- readRDS(file.path(OUT,'imputations.rds'))
form <- c(
 M1='Surv(days,outcome_3y_logmar2)~age+sex+cluster(ID)',
 M2='Surv(days,outcome_3y_logmar2)~age+sex+HbA1c+eGFR+cluster(ID)',
 M3=MODEL3,M4=MODEL4)
fits <- lapply(form,function(f) with(imp,survival::coxph(model_formula(f))))
pooled <- lapply(fits,mice::pool)
extract <- function(i) {
 x <- summary(pooled[[i]],conf.int=TRUE,exponentiate=TRUE)
 x$model <- names(form)[i]; x
}
res <- dplyr::bind_rows(lapply(seq_along(pooled),extract))
write_csv(res,'cox_models_M1_M4.csv')
# Nested-model pooled Wald tests (mice::D1), not conventional likelihood-ratio tests.
comp <- list(M2_vs_M1=mice::D1(fits$M2,fits$M1),
             M3_vs_M2=mice::D1(fits$M3,fits$M2),
             M4_vs_M3=mice::D1(fits$M4,fits$M3))
capture.output(comp,file=file.path(OUT,'pooled_wald_D1.txt'))
saveRDS(fits,file.path(OUT,'cox_fits.rds'))
message('Pooled Cox results and D1 Wald tests saved')
