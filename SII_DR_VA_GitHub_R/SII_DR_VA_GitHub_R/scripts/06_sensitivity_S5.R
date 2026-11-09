source('scripts/00_config.R')
raw<-readRDS(file.path(OUT,'data_prepared.rds'))
completed<-get_completed()
cov<-c('age','sex','HbA1c','eGFR','DR_grading','DME','CAT','GLAU','SII100')
run_one<-function(label,filter_fun=function(d) rep(TRUE,nrow(d)),outcome='outcome_3y_logmar2',exclude=NULL,add=NULL,use_mi=TRUE){
 ds<-if(use_mi) completed else list(raw)
 fits<-lapply(ds,function(d){
   d<-d[filter_fun(d),,drop=FALSE]
   cv<-c(setdiff(cov,exclude),add)
   f<-as.formula(paste0('Surv(days,',outcome,')~',paste(cv,collapse='+'),'+cluster(ID)'))
   survival::coxph(f,data=d)
 })
 s<-if(use_mi) summary(mice::pool(fits),conf.int=TRUE,exponentiate=TRUE) else {
   sm<-summary(fits[[1]])
   data.frame(term=rownames(sm$coefficients),estimate=exp(sm$coefficients[,'coef']),
              `2.5 %`=exp(sm$conf.int[,'lower .95']),`97.5 %`=exp(sm$conf.int[,'upper .95']),
              p.value=sm$coefficients[,'Pr(>|z|)'],check.names=FALSE)
 }
 e<-s[s$term=='SII100',,drop=FALSE]
 data.frame(Analysis=label,HR=e$estimate,lower=e$`2.5 %`,upper=e$`97.5 %`,P=e$p.value,
            n_eyes=sum(filter_fun(raw)),events=sum(raw$outcome_3y_logmar2[filter_fun(raw)]))
}
rows<-list(run_one('Main analysis'),
 run_one('Excluding baseline DME',function(d) d$DME=='0',exclude='DME'),
 run_one('Excluding baseline glaucoma',function(d) d$GLAU=='0',exclude='GLAU'),
 run_one('Follow-up >=180 days',function(d) d$days>=180))
# Sensitivity analyses that require columns absent from baseline source are NOT invented.
# Alternative 0.3 endpoint has its own event times; check the column for 0.3 endpoint.
if(all(c('outcome_3y_logmar3') %in% names(raw))){
 # Completed MI frames omit logmar3; merge by row order only if row IDs/order intact.
 completed03<-Map(function(d,i){d$outcome_3y_logmar3<-raw$outcome_3y_logmar3;d},completed,seq_along(completed))
 old<-completed;completed<-completed03
 rows[[length(rows)+1]]<-run_one('Visual decline >=0.3',outcome='outcome_3y_logmar3')
 completed<-old
}
cc <- stats::complete.cases(raw[c('age','sex','HbA1c','eGFR','DR_grading','DME','CAT','GLAU','SII100','days','outcome_3y_logmar2')])
rows[[length(rows)+1]]<-run_one('Complete case',function(d) cc,use_mi=FALSE)
write_csv(dplyr::bind_rows(rows),'Supplementary_Table_S5_recomputed_part.csv')
# Additional requested sensitivities require separate genuine source data:
# baseline_BCVA, 30-day SII match flag, and repeated BCVA encounter timeline.
writeLines(c('Not produced from this workbook alone:',
 'Baseline BCVA adjustment: baseline BCVA absent from analysis workbook',
 '30-day SII window: no laboratory test dates/window flag',
 'Repeated confirmation: original visit-level BCVA sequences absent',
 'Worse-eye: baseline BCVA needed to select worse eye',
 'Check 0.3 endpoint event-time column before treating baseline days as valid.'),
 file.path(OUT,'sensitivity_missing_inputs.txt'))
