source('scripts/00_config.R')
# Conditional-on-completed-imputations cluster bootstrap, NOT nested imputation.
# This code reports an MI-averaged point estimate and percentile interval conditional
# on existing 20 completed datasets. Does not redo MICE within bootstrap replicates.
completed<-get_completed()
bootstrap_reps <- as.integer(Sys.getenv('SII_BOOTSTRAP_B',unset='1000'))
# concordance can take risk score (reverse=TRUE), censoring-aware Harrell C
c_index <- function(d,formula) {
 fit<-survival::coxph(model_formula(formula),data=d)
 lp<-stats::predict(fit,type='lp')
 as.numeric(survival::concordance(survival::Surv(d$days,d$outcome_3y_logmar2)~lp,reverse=TRUE)$concordance)
}
delta_one<-function(d) c_index(d,MODEL4)-c_index(d,MODEL3)
point<-mean(vapply(completed,delta_one,numeric(1)))
ids<-unique(completed[[1]]$ID)
set.seed(20261109)
replicates<-rep(NA_real_,bootstrap_reps)
for(b in seq_len(bootstrap_reps)){
 selected<-sample(ids,length(ids),replace=TRUE)
 # Reassign sampled cluster IDs to ensure duplicate patients are independent
 indices<-unlist(lapply(seq_along(selected),function(j) which(completed[[1]]$ID==selected[j])),use.names=FALSE)
 counts<-vapply(selected,function(i) sum(completed[[1]]$ID==i),integer(1))
 new_id<-rep(seq_along(selected),counts)
 replicates[b]<-tryCatch(mean(vapply(completed,function(df){
  z<-df[indices,,drop=FALSE];z$ID<-as.character(new_id);delta_one(z)
 },numeric(1))),error=function(e) NA_real_)
 if(b%%100==0) message('Bootstrap ',b,'/',bootstrap_reps)
}
valid<-replicates[is.finite(replicates)]
if(length(valid)<.9*bootstrap_reps) warning('Too many failed resamples; investigate results')
ci<-stats::quantile(valid,c(.025,.975),na.rm=TRUE)
write_csv(data.frame(method='20 MI completed datasets; patient cluster bootstrap conditional on imputation',
  B=bootstrap_reps,successful=length(valid),delta=point,lower=ci[1],upper=ci[2]),
  'Delta_Cindex_95CI_conditional_MI.csv')
write_csv(data.frame(replicate=seq_along(replicates),delta=replicates),'Delta_Cindex_bootstrap_replicates.csv')
message(sprintf('Delta C-index %.4f [%.4f, %.4f]',point,ci[1],ci[2]))
