# Author: Chunhao Guo
source('/LOCAL_PROJECT_ROOT/84_JCO_PO_Upgrade_20261113/05_analysis/A_clinical_consequence/scripts/common.R')
checks<-list(); frames<-list(); add<-function(tr,test,observed,expected,tolerance=1e-8){
 checks[[length(checks)+1]]<<-data.frame(cohort=tr,check=test,observed=observed,expected=expected,tolerance=tolerance,pass=isTRUE(abs(observed-expected)<=tolerance))
}
membership<-read.csv(file.path(ROOT,'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv'))
for(tr in names(inputs)) {
 d<-inputs[[tr]]$d;n<-nrow(d);add(tr,'RNA-linked patients',n,if(tr=='JAVELIN')726 else 823)
 add(tr,'PFS events',sum(d$event),if(tr=='JAVELIN')358 else 529)
 stopifnot(!anyDuplicated(d$ID),all(d$time>0),all(d$event%in%0:1),all(d$trt%in%0:1))
 for(tar in c('Angio36','Angio16')) {
  M<-inputs[[tr]]$M[[tar]];add(tr,paste(tar,'complete gene count'),nrow(M),if(tar=='Angio36')36 else 16)
  add(tr,paste(tar,'final manuscript gene-list mismatch'),length(setdiff(membership$gene[membership[[tar]]],rownames(M)))+length(setdiff(rownames(M),membership$gene[membership[[tar]]])),0)
  for(re in 1:5) {
   fr<-held_frame(tr,tar,re)
   sh<-saved_held[saved_held$trial==tr&saved_held$comparison==paste('Algorithm',tar)&saved_held$rep==re,]
   sh<-sh[match(fr$ID,sh$ID),];stopifnot(!anyNA(sh$ID),!anyDuplicated(sh$ID))
   add(tr,paste(tar,'partition',re,'patient transition mismatch'),sum(fr$transition!=sh$transition),0)
   add(tr,paste(tar,'partition',re,'fold mismatch'),sum(fr$fold!=sh$fold),0)
   cm<-cutoffs[cutoffs$trial==tr&cutoffs$comparison==paste('Algorithm',tar)&cutoffs$rep==re,]
   add(tr,paste(tar,'partition',re,'threshold max error'),max(abs(fr$threshold_I0-cm$cutoff0[match(fr$fold,cm$fold)]),abs(fr$threshold_I1-cm$cutoff1[match(fr$fold,cm$fold)])),0)
   ss<-saved_summary[saved_summary$trial==tr&saved_summary$comparison==paste('Algorithm',tar)&saved_summary$rep==re,]
   add(tr,paste(tar,'partition',re,'reclassified count'),sum(fr$high_I0!=fr$high_I1),ss$switches)
   frames[[paste(tr,tar,re,sep='__')]]<-fr
  }
  sf<-saved_full[saved_full$trial==tr&saved_full$comparison==paste('Algorithm',tar),]
  sf<-sf[match(d$ID,sf$ID),];z0<-full_scores(M,'I0');z1<-full_scores(M,'I1');trans<-transition(z0>=median(z0),z1>=median(z1))
  add(tr,paste(tar,'full-cohort patient transition mismatch'),sum(sf$transition!=trans),0)
  if(tar=='Angio36') add(tr,'Angio36 full-cohort median disagreement',sum(sf$switch),if(tr=='JAVELIN')48 else 80)
  p<-policy[[tr]];p<-p[match(d$ID,p$ID),];stopifnot(!anyNA(p$ID),!anyDuplicated(p$ID))
  add(tr,paste(tar,'policy arm mismatch'),sum(p$arm!=d$trt),0)
  add(tr,paste(tar,'policy fold mismatch'),sum(p$fold!=folds[[tr]][[1]]),0)
  nd<-sum(p[[paste0('d_',tar,'_I0')]]!=p[[paste0('d_',tar,'_I1')]])
  expected<-if(tar=='Angio36')if(tr=='JAVELIN')40 else 76 else if(tr=='JAVELIN')9 else 15
  add(tr,paste(tar,'12-month policy disagreement'),nd,expected)
 }
 v<-sapply(1:5,function(re)100*mean(with(frames[[paste(tr,'Angio36',re,sep='__')]],high_I0!=high_I1)))
 add(tr,'Angio36 held-out minimum rounded percent',round(min(v),2),if(tr=='JAVELIN')5.79 else 9.72)
 add(tr,'Angio36 held-out maximum rounded percent',round(max(v),2),if(tr=='JAVELIN')7.44 else 10.33)
}
tab<-do.call(rbind,checks);wc(tab,'reproduction_checks.csv');print(tab)
if(!all(tab$pass))stop('Reproduction failed: no new estimates permitted; inspect reproduction_checks.csv')
saveRDS(frames,file.path(OUT,'results/verified_heldout_frames.rds'))
writeLines('PASS: all manuscript-number and patient-level reproduction checks passed before any new estimate.',file.path(OUT,'logs/REPRODUCTION_PASS.txt'))
cat('All reproduction checks PASS; new estimation may begin.\n')
