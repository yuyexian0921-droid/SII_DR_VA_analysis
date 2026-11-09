source('scripts/00_config.R')
imp<-readRDS(file.path(OUT,'imputations.rds'))
pdf(file.path(OUT,'mice_trace.pdf'),width=10,height=8)
plot(imp)
dev.off()
long<-mice::complete(imp,action='long',include=TRUE)
vars<-c('SII100','HbA1c','eGFR')
stack<-dplyr::bind_rows(lapply(vars,function(v){
 imputed_index<-which(is.na(imp$data[[v]]))
 obs<-data.frame(variable=v,value=imp$data[[v]][!is.na(imp$data[[v]])],source='Observed',imp=0)
 z<-long[long$.imp>0 & as.integer(long$.id) %in% imputed_index,c('.imp',v)]
 imp_vals<-data.frame(variable=v,value=z[[v]],source='Imputed',imp=z$.imp)
 dplyr::bind_rows(obs,imp_vals)
}))
p<-ggplot2::ggplot(stack,ggplot2::aes(x=value,color=source,group=interaction(source,imp)))+
 ggplot2::geom_density(alpha=.4)+ggplot2::facet_wrap(~variable,scales='free',nrow=1)+ggplot2::theme_bw()
ggplot2::ggsave(file.path(OUT,'Figure_S5_observed_vs_imputed.pdf'),p,width=11,height=4)
