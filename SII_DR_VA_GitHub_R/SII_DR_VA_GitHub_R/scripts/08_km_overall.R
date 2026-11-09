source('scripts/00_config.R')
d<-readRDS(file.path(OUT,'data_prepared.rds'))
# This is a Kaplan-Meier estimate of the observed time-to-documented-event data.
# Do not call it a 30-day landmark analysis, or assume original event histories
# can be reconstructed from the eye-level event/censor summaries alone.
fit<-survival::survfit(survival::Surv(days,outcome_3y_logmar2)~1,data=d,conf.type='log-log')
s<-summary(fit,times=c(365,730,1095),extend=TRUE)
out<-data.frame(day=c(365,730,1095),eyes_at_risk=s$n.risk,
 cumulative_event_probability=1-s$surv,lower95=1-s$upper,upper95=1-s$lower)
write_csv(out,'KM_1_2_3_year.csv')
pdf(file.path(OUT,'KM_overall.pdf'),width=7,height=5)
plot(fit,xlim=c(0,1095),xlab='Days since index',ylab='Event-free probability',conf.int=TRUE)
dev.off()
