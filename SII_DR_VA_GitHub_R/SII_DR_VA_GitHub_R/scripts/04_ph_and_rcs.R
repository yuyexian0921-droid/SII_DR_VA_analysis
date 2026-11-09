source('scripts/00_config.R')
imp<-readRDS(file.path(OUT,'imputations.rds'))
d<-mice::complete(imp,1)
# PH diagnostic on representative completion. NOT automatically pooled across 20 imputations.
fit<-survival::coxph(model_formula(MODEL4),data=d,x=TRUE)
ph<-survival::cox.zph(fit)
ph_table<-as.data.frame(ph$table);ph_table$term<-rownames(ph_table)
write_csv(ph_table,'ph_Schoenfeld_imputation1.csv')
pdf(file.path(OUT,'ph_Schoenfeld_imputation1.pdf'),width=9,height=8)
plot(ph);dev.off()
# Four SII knots AT the manuscript's stated 5th, 35th, 65th and 95th percentiles.
# This is a check on one completed dataset; manuscript RCS results must be validated
# against the final original analysis (not automatically pooled MI).
knots<-stats::quantile(d$SII100,c(.05,.35,.65,.95),na.rm=TRUE)
dd<-rms::datadist(d); dd$limits['Adjust to','SII100']<-stats::median(d$SII100)
options(datadist='dd')
rcs_fit<-rms::cph(survival::Surv(days,outcome_3y_logmar2)~
   age+sex+HbA1c+eGFR+DR_grading+DME+CAT+GLAU+
   rms::rcs(SII100,knots)+survival::cluster(ID),data=d,x=TRUE,y=TRUE,surv=TRUE)
capture.output(rms::anova(rcs_fit),file=file.path(OUT,'rcs_anova_imputation1.txt'))
pr<-as.data.frame(rms::Predict(rcs_fit,SII100,fun=exp,ref.zero=TRUE,conf.int=.95))
write_csv(pr,'rcs_points_imputation1.csv')
p<-ggplot2::ggplot(pr,ggplot2::aes(x=SII100,y=yhat))+ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=.2)+ggplot2::geom_line()+ggplot2::geom_hline(yintercept=1,linetype='dashed')+ggplot2::labs(x='SII (per 100 units)',y='HR (95% CI)')+ggplot2::theme_classic()
ggplot2::ggsave(file.path(OUT,'rcs_imputation1.pdf'),p,width=7,height=5)
message('RCS/PH are labeled single-imputation diagnostics, NOT manuscript-confirmed pooled results')
