source('scripts/00_config.R')
imp<-readRDS(file.path(OUT,'imputations.rds'))
raw<-readRDS(file.path(OUT,'data_prepared.rds'))
completed<-get_completed()
add_groups<-function(d){
 d$Age_group<-factor(ifelse(d$age<60,'<60','>=60'),levels=c('<60','>=60'))
 d$HbA1c_group<-factor(ifelse(d$HbA1c<7,'<7','>=7'),levels=c('<7','>=7'))
 d$eGFR_group<-factor(ifelse(d$eGFR<60,'<60','>=60'),levels=c('<60','>=60'))
 d
}
completed<-lapply(completed,add_groups)
raw<-add_groups(raw)
spec<-list('DR severity'='DR_grading',DME='DME',Cataract='CAT',Glaucoma='GLAU',
           Sex='sex',Age='Age_group',HbA1c='HbA1c_group',eGFR='eGFR_group')
# For subgroup-specific Cox, remove the stratification factor from covariates.
# This prevents singular terms; remain adjusted for the other main-model variables.
base_covars<-c('age','sex','HbA1c','eGFR','DR_grading','DME','CAT','GLAU','SII100')
remove_covars<-c('DR severity'='DR_grading',DME='DME',Cataract='CAT',Glaucoma='GLAU',
                 Sex='sex',Age='age',HbA1c='HbA1c',eGFR='eGFR')
fit_levels<-function(group_name,var,level) {
 fs<-lapply(completed,function(d){
  z<-d[!is.na(d[[var]]) & as.character(d[[var]])==level,,drop=FALSE]
  cv<-setdiff(base_covars,remove_covars[[group_name]])
  f<-as.formula(paste0('Surv(days,outcome_3y_logmar2)~',paste(cv,collapse='+'),'+cluster(ID)'))
  survival::coxph(f,data=z)
 })
 pooled<-summary(mice::pool(fs),conf.int=TRUE,exponentiate=TRUE)
 effect<-pooled[pooled$term=='SII100',,drop=FALSE]
 sizes<-vapply(completed,function(d)sum(!is.na(d[[var]]) & as.character(d[[var]])==level),numeric(1))
 events<-vapply(completed,function(d)sum(d$outcome_3y_logmar2[!is.na(d[[var]]) & as.character(d[[var]])==level]),numeric(1))
 observed<-raw[!is.na(raw[[var]]) & as.character(raw[[var]])==level,,drop=FALSE]
 data.frame(Subgroup=group_name,Level=level,
  observed_eyes=nrow(observed),observed_events=sum(observed$outcome_3y_logmar2),
  MI_eyes_median=median(sizes),MI_eyes_min=min(sizes),MI_eyes_max=max(sizes),
  MI_events_median=median(events),MI_events_min=min(events),MI_events_max=max(events),
  HR=effect$estimate,Lower95=effect$`2.5 %`,Upper95=effect$`97.5 %`,P=effect$p.value)
}
result<-list()
for(g in names(spec))for(l in levels(factor(completed[[1]][[spec[[g]]]]))){
 message('Subgroup ',g,': ',l)
 result[[length(result)+1]]<-tryCatch(fit_levels(g,spec[[g]],l),error=function(e){
  warning(g,' ',l,' failed: ',conditionMessage(e));NULL})
}
out<-dplyr::bind_rows(result)
write_csv(out,'Supplementary_Table_S4_with_counts_MI.csv')
# HbA1c/eGFR observed counts DO NOT sum to 1725; use MI counts explicitly when
# reporting subgroup-adjusted HR on imputed datasets.
# Interaction tests: categorical split interacts with SII100, while continuous covariate
# adjusted in original model should not simultaneously duplicate its group derivative.
for(g in names(spec)){
 var<-spec[[g]]
 cv<-setdiff(base_covars,remove_covars[[g]])
 f<-as.formula(paste0('Surv(days,outcome_3y_logmar2)~',paste(cv,collapse='+'),
                      '+',var,'+SII100:',var,'+cluster(ID)'))
 fits<-lapply(completed,function(d)survival::coxph(f,data=d))
 term<-grep('SII100:',names(stats::coef(fits[[1]])),value=TRUE)
 if(length(term)!=1) {warning('Interaction requires joint Wald test: ',g);next}
 sumtab<-summary(mice::pool(fits),conf.int=TRUE)
 row<-sumtab[sumtab$term==term,]
 write_csv(data.frame(group=g,term=term,p_interaction=row$p.value),paste0('interaction_',gsub(' ','_',g),'.csv'))
}
