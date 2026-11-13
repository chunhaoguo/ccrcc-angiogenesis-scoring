# Author / creator / last modified by: Chunhao Guo
# post hoc follow-up (specified after modules B2 and D)
source('scripts/00_common.R')
f1<-rbind(rd('results/F1_GENOTYPE_ALL.csv'),rd('results/F1_SPECIMEN_SENSITIVITY_ALL.csv'))
f3<-rd('results/F3_GROUP_ANCHORS_ALL.csv');checks<-list()
add<-function(check,error=0,tolerance=1e-12){checks[[length(checks)+1]]<<-data.table(check,error,tolerance,pass=is.finite(error)&&error<=tolerance)}
for(k in seq_len(nrow(f1))){
 r<-f1[k];if(!is.finite(r$estimate))next
 z<-readRDS(file.path(SCR,paste0('F1_',r$cohort,'_',r$signature,'_',r$anchor,'_',r$mode,'.rds')))
 b<-z$bootstrap[,r$operator];delta<-b-z$bootstrap[,'O1']
 add(paste('F1 stored bootstrap CI',r$cohort,r$signature,r$anchor,r$operator,r$mode),max(abs(c(interval(b)-c(r$lower,r$upper),interval(delta)-c(r$difference_lower,r$difference_upper)))))
 add(paste('F1 stored bootstrap p',r$cohort,r$signature,r$anchor,r$operator,r$mode),abs(bp(delta,r$difference_vs_O1)-r$p_difference))
}
f1[,q_check:=p.adjust(p_difference,'BH',n=21),by=.(cohort,signature,mode)]
add('F1 all BH family values',max(abs(f1$q_check-f1$q_BH),na.rm=TRUE))
for(k in seq_len(nrow(f3))){
 r<-f3[k];gr<-gsub('[^A-Za-z0-9]','',r$gene_group)
 if(r$cohort=='CPTAC'){
  z<-readRDS(file.path(SCR,paste0('F3_CPTAC_',gr,'.rds')));j<-match(r$operator,c('O1','O2'));b<-z$bootstrap[,j];delta<-b-z$bootstrap[,1]
 }else{
  z<-readRDS(file.path(SCR,paste0('F3_',r$cohort,'_Angio36 fixed-group decomposition_',r$anchor,'_main_',gr,'.rds')))
  b<-z$bootstrap[,paste0(r$gene_group,'__',r$operator)];delta<-b-z$bootstrap[,paste0(r$gene_group,'__O1')]
 }
 add(paste('F3 stored paired CI',r$cohort,r$anchor,r$gene_group,r$operator),max(abs(c(interval(b)-c(r$lower,r$upper),interval(delta)-c(r$difference_lower,r$difference_upper)))))
 add(paste('F3 stored bootstrap p',r$cohort,r$anchor,r$gene_group,r$operator),abs(bp(delta,r$difference_vs_O1)-r$p_difference))
}
f3[,q_check:=1];f3[operator=='O2',q_check:=p.adjust(p_difference,'BH',n=4),by=.(cohort,anchor)]
add('F3 separate descriptive BH families',max(abs(f3$q_check-f3$q_BH)))
# Direct group-score formula, independent of shared wrapper functions.
ob<-readRDS(file.path(A,'D_biology/results/qualified_inputs.rds'));groups<-rd('results/F3_FIXED_GENE_GROUPS.csv');saved<-rd('results/F3_GROUP_SCORES_LOCAL_ONLY.csv')
for(co in c(cohorts,'CPTAC')){
 if(co=='CPTAC'){q<-readRDS('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2/CPTAC.rds');X<-q$x[groups$gene,,drop=FALSE];ids<-q$patient;rm(q);gc()}else{
  ids<-ob$ids[[co]];a<-ob$mat[[co]]$Angio16;b<-ob$mat[[co]]$Angio36;X<-rbind(a,b[!rownames(b)%in%rownames(a),,drop=FALSE])
 }
 for(gr in unique(groups$group)){
  Y<-X[groups[group==gr,gene],,drop=FALSE];a<-colMeans(Y);b<-colMeans((Y-rowMeans(Y))/apply(Y,1,sd));a<-(a-mean(a))/sd(a);b<-(b-mean(b))/sd(b)
  s<-saved[cohort==co&gene_group==gr];s<-s[match(ids,ID)]
  add(paste('direct complete-group formula',co,gr),max(abs(a-s$O1),abs(b-s$O2)))
 }
 rm(X);gc()
}
qa<-rbindlist(checks);wr(qa,'BOOTSTRAP_VALIDATION.csv');stopifnot(all(qa$pass));cat('PASS',nrow(qa),'bootstrap/CI/BH/direct-group checks\n')
