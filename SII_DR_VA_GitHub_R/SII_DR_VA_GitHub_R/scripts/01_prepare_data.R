source('scripts/00_config.R')
if(!file.exists(INPUT)) stop('Input workbook not found: ',INPUT)
d <- as.data.frame(readxl::read_excel(INPUT,sheet=1))
need <- c('ID','birth_year','sex','index_date','days','outcome_3y_logmar2',
          'DR_grading','DME','HbA1c','eGFR','SII','CAT','GLAU')
stopifnot(all(need %in% names(d)))
# Excel serial dates / date or string representations.
to_date <- function(x) {
 if(inherits(x,'Date')) return(x)
 if(inherits(x,'POSIXt')) return(as.Date(x))
 if(is.numeric(x)) return(as.Date(x,origin='1899-12-30'))
 as.Date(as.character(x))
}
d$index_date <- to_date(d$index_date)
d$age <- as.integer(format(d$index_date,'%Y')) - as.numeric(d$birth_year)
d$ID <- as.character(d$ID)
d$sex <- factor(dplyr::case_when(compact(d$sex) %in% c('男','M','Male')~'Male',
                                   compact(d$sex) %in% c('女','F','Female')~'Female',TRUE~NA_character_),
                levels=c('Female','Male'))
d$DR_grading <- factor(d$DR_grading,levels=c('NPDR','PDR'))
for(v in c('DME','CAT','GLAU')) d[[v]] <- factor(as.character(d[[v]]),levels=c('0','1'))
for(v in c('HbA1c','eGFR','SII','days','age')) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
d$SII100 <- d$SII/100
d$outcome_3y_logmar2 <- as.integer(d$outcome_3y_logmar2)
if('outcome_3y_logmar3' %in% names(d)) d$outcome_3y_logmar3 <- as.integer(d$outcome_3y_logmar3)
stopifnot(nrow(d)==1725L, dplyr::n_distinct(d$ID)==914L,
          sum(d$outcome_3y_logmar2,na.rm=TRUE)==748L,
          !anyNA(d[c('ID','age','sex','days','outcome_3y_logmar2','DR_grading','DME','CAT','GLAU')]),
          all(d$days>=0),all(d$outcome_3y_logmar2 %in% 0:1))
saveRDS(d,file.path(OUT,'data_prepared.rds'))
missing <- data.frame(variable=names(d), n_missing=vapply(d,function(x)sum(is.na(x)),numeric(1)))
missing$percent <- 100*missing$n_missing/nrow(d)
write_csv(missing,'missingness.csv')
write_csv(data.frame(n_patients=dplyr::n_distinct(d$ID),n_eyes=nrow(d),events=sum(d$outcome_3y_logmar2)),'cohort_counts.csv')
message('Prepared ',nrow(d),' eyes / ',dplyr::n_distinct(d$ID),' patients')
