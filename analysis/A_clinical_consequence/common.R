# Author: Chunhao Guo
# Locked Module A; one serial R process, no new packages.
options(stringsAsFactors=FALSE, width=160)
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages(library(survival))
set.seed(20261113)
ROOT <- '/LOCAL_PROJECT_ROOT'
OUT <- file.path(ROOT,'84_JCO_PO_Upgrade_20261113/05_analysis/A_clinical_consequence')
for (d in c('results','logs','plots')) dir.create(file.path(OUT,d),recursive=TRUE,showWarnings=FALSE)
setwd(OUT) # Any automatic graphics device must also remain in the authorised output tree.
options(device=function(...)grDevices::pdf(file.path(OUT,'logs/default_device_diagnostic.pdf'),...))
path37 <- file.path(ROOT,'37_A_clinical_workflow_validation_20260917/local_only')
path48 <- file.path(ROOT,'48_Polan_ccRCC_GenomeMedicine_20260922')
path67 <- file.path(ROOT,'67_Mathematical_Extensions_20261008')
inputs <- readRDS(file.path(path37,'INPUTS.rds'))
folds <- readRDS(file.path(path37,'FOLD_ASSIGNMENTS.rds'))
saved_held <- read.csv(gzfile(file.path(path48,'local_only/HELDOUT_PATIENT_RANKS.csv.gz')))
saved_full <- read.csv(gzfile(file.path(path48,'local_only/PATIENT_RANKS.csv.gz')))
cutoffs <- read.csv(file.path(path48,'new_analysis/aggregate/THRESHOLD_MANIFEST.csv'))
saved_summary <- read.csv(file.path(ROOT,'70_npj_Visual_Narrative_Final_20261011/figures/source_data/sc/new/HELDOUT_REPEAT_SUMMARY.csv'))
policy <- lapply(names(inputs),function(tr) {
 p<-read.csv(file.path(path67,'local_only',paste0('B_patient_policy_',tr,'.csv')))
 p[p$tau==12,]
}); names(policy)<-names(inputs)
wc <- function(x,f) write.csv(x,file.path(OUT,'results',f),row.names=FALSE,na='NA')
transition <- function(g0,g1) ifelse(is.na(g0)|is.na(g1),'Unavailable',
 ifelse(g0,ifelse(g1,'Stable high','High to low'),ifelse(g1,'Low to high','Stable low')))
apply_reference <- function(M,r) {
 if(any(!is.finite(r$gene_sd)|r$gene_sd<=1e-8)) return(rep(NA_real_,ncol(M)))
 raw<-colMeans(sweep(sweep(M,1,r$gene_mean,'-'),1,r$gene_sd,'/'))
 (raw-r$center)/r$unit
}
full_scores <- function(M,method) {
 if(method=='I0') raw<-colMeans(M) else {
  ss<-apply(M,1,sd);if(any(!is.finite(ss)|ss<=1e-8))return(rep(NA_real_,ncol(M)))
  raw<-colMeans(sweep(sweep(M,1,rowMeans(M),'-'),1,ss,'/'))
 }; (raw-mean(raw))/sd(raw)
}
held_frame <- function(tr,tar,re,prob=.5) {
 obj<-inputs[[tr]]; d<-obj$d; M<-obj$M[[tar]]; n<-nrow(d)
 z0<-z1<-c0<-c1<-rep(NA_real_,n)
 refs<-readRDS(file.path(path37,paste0('CV_REFERENCES__',tr,'__',tar,'.rds')))
 for(fo in 1:5){r<-refs[[paste(re,fo,sep='__')]];it<-which(folds[[tr]][[re]]!=fo);iv<-which(folds[[tr]][[re]]==fo)
  stopifnot(setequal(r$training_ID,d$ID[it]),setequal(r$test_ID,d$ID[iv]),!length(intersect(r$training_ID,r$test_ID)))
  s0<-apply_reference(M,r$versions$I0);s1<-apply_reference(M,r$versions$I1)
  z0[iv]<-s0[iv];z1[iv]<-s1[iv];c0[iv]<-quantile(s0[it],prob,type=7);c1[iv]<-quantile(s1[it],prob,type=7)
 }
 data.frame(trial=tr,target=tar,partition=re,threshold_probability=prob,ID=d$ID,fold=folds[[tr]][[re]],
  score_I0=z0,score_I1=z1,threshold_I0=c0,threshold_I1=c1,
  signed_distance_I0=z0-c0,signed_distance_I1=z1-c1,absolute_distance_I0=abs(z0-c0),absolute_distance_I1=abs(z1-c1),
  high_I0=z0>=c0,high_I1=z1>=c1,transition=transition(z0>=c0,z1>=c1))
}
