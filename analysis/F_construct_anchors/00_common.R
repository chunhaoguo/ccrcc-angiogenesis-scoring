# Author / creator / last modified by: Chunhao Guo
# post hoc follow-up (specified after modules B2 and D)
options(stringsAsFactors=FALSE, warn=1)
Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1)
.libPaths(c('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B1/Rlib', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(1);set.seed(20261113)
ROOT <- '/LOCAL_PROJECT_ROOT'
A <- file.path(ROOT,'84_JCO_PO_Upgrade_20261113/05_analysis')
OUT <- file.path(A,'F_construct_anchors');setwd(OUT)
SCR <- '/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/F'
LABEL <- 'post hoc follow-up (specified after modules B2 and D)'
rd <- function(f) fread(f,na.strings=c('NA',''))
wr <- function(d,f) {d<-as.data.table(copy(d));d[,analysis_label:=LABEL];fwrite(d,file.path(OUT,'results',f),na='NA')}
source(file.path(A,'shared/operators.R'))
cohorts <- c('JAVELIN','IMmotion151','CheckMate025'); sigs<-c('Angio36','Angio16');ops<-paste0('O',1:7)
load_scores <- function(){
 ss<-list()
 for(co in cohorts[1:2]){
  q<-rd(file.path(A,'B1_operator_panel_trials/results',paste0(co,'_scores_local_only.csv')))
  ss[[co]]<-q[signature%in%sigs,.(cohort=trial,ID,signature,operator,score)]
 }
 for(co in c('CheckMate025','CPTAC')){
  nm<-if(co=='CheckMate025')'CM025'else co
  q<-rd(file.path(A,'B2_independent_cohorts/results',paste0(nm,'_scores.csv')))
  q<-melt(q,id.vars=c('signature','sample','patient'),measure.vars=ops,variable.name='operator',value.name='score')
  ss[[co]]<-q[,.(cohort=co,ID=patient,signature,operator=as.character(operator),score)]
 }
 rbindlist(ss)
}
interval<-function(x) if(any(is.finite(x)))as.numeric(quantile(x,c(.025,.975),na.rm=TRUE,type=7)) else c(NA_real_,NA_real_)
# Centred two-sided bootstrap tail test; plus-one correction gives finite Monte Carlo resolution.
bp<-function(b,theta){b<-b[is.finite(b)];if(!length(b)||!is.finite(theta))return(NA_real_);(1+sum(abs(b-theta)>=abs(theta)))/(length(b)+1)}
