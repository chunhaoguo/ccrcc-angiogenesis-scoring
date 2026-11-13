# Author / creator / last modified by: Chunhao Guo
source('scripts/00_common.R')
scores<-load_scores();stopifnot(!anyDuplicated(scores[,.(cohort,ID,signature,operator)]))
gt<-rd(file.path(A,'D_biology/results/GENOTYPE_PATIENT_JOIN.csv'))[gene%in%c('VHL','PBRM1','BAP1')]
stopifnot(!anyDuplicated(gt[,.(cohort,ID,gene)]))
qa<-list(); add<-function(check,cohort,signature,observed,expected,tolerance=1e-10){qa[[length(qa)+1]]<<-data.table(check,cohort,signature,observed,expected,tolerance,pass=is.finite(observed)&&abs(observed-expected)<=tolerance)}
for(co in c(cohorts,'CPTAC'))for(sig in sigs){
 d<-scores[cohort==co&signature==sig];S<-dcast(d,ID~operator,value.var='score');m<-as.matrix(S[,..ops])
 if(co%in%cohorts[1:2]){
  old<-rd(file.path(A,'B1_operator_panel_trials/results/pairwise_metrics.csv'))[trial==co&signature==sig]
  for(k in seq_len(nrow(old))){z<-old[k];if(is.finite(z$disagreement_n))add(paste('all saved pairwise disagreement',z$operator1,z$operator2),co,sig,sum((S[[z$operator1]]>=median(S[[z$operator1]]))!=(S[[z$operator2]]>=median(S[[z$operator2]]))),z$disagreement_n,0)}
 }else{
  nm<-if(co=='CheckMate025')'CM025'else co
  old<-rd(file.path(A,'B2_independent_cohorts/results',paste0(nm,'_pairwise_metrics.csv')))[signature==sig]
  for(k in seq_len(nrow(old))){z<-old[k];if(is.finite(z$disagreement_percent))add(paste('all saved pairwise disagreement',z$operator_a,z$operator_b),co,sig,100*mean((S[[z$operator_a]]>=median(S[[z$operator_a]]))!=(S[[z$operator_b]]>=median(S[[z$operator_b]]))),z$disagreement_percent)}
 }
 if(co%in%cohorts){
  old<-rd(file.path(A,'D_biology/results/PATIENT_SCORES_RANKS.csv'))[cohort==co&signature==sig];old<-old[match(S$ID,ID)]
  add('module D O1 score max absolute error',co,sig,max(abs(S$O1-old$I0)),0,1e-8)
  add('module D O2 score max absolute error',co,sig,max(abs(S$O2-old$I1)),0,1e-8)
 }
}
for(co in cohorts[1:2]){
 s<-dcast(scores[cohort==co&signature=='Angio36'],ID~operator,value.var='score')
 add('required headline disagreement count',co,'Angio36',sum((s$O1>=median(s$O1))!=(s$O2>=median(s$O2))),if(co=='JAVELIN')48 else 80,0)
 add('required headline denominator',co,'Angio36',nrow(s),if(co=='JAVELIN')726 else 823,0)
}
counts<-gt[,. (n_altered=sum(status=='Altered'),n_source_negative=sum(status=='No_reported_alteration'),n_unknown=sum(status=='Unknown_not_assayed'),RNA_n=.N),by=.(cohort,gene)]
audit<-rd(file.path(A,'D_biology/results/GENOTYPE_JOIN_AUDIT.csv'))
for(k in seq_len(nrow(counts))){z<-counts[k];e<-audit[cohort==z$cohort&gene==z$gene];for(v in c('n_altered','n_source_negative','RNA_n'))add(paste('module D genotype',z$gene,v),z$cohort,'all',z[[v]],e[[v]],0);add(paste('module D genotype',z$gene,'unknown'),z$cohort,'all',z$n_unknown,e$n_unknown_not_assayed,0)}
for(v in list(c('JAVELIN','VHL',383,308,35),c('IMmotion151','PBRM1',328,374,121),c('CheckMate025','PBRM1',63,112,75))){z<-counts[cohort==v[1]&gene==v[2]];for(j in 1:3)add(paste('required genotype count',v[2],c('altered','negative','unknown')[j]),v[1],'all',unlist(z[,.(n_altered,n_source_negative,n_unknown)])[j],as.numeric(v[j+2]),0)}
anchor<-rd(file.path(A,'B2_independent_cohorts/results/CPTAC_protein_proxy.csv'))
pa<-rd(file.path(A,'B2_independent_cohorts/results/CPTAC_protein_anchor.csv'))
for(sig in sigs)for(op in ops){s<-scores[cohort=='CPTAC'&signature==sig&operator==op];v<-s$score[match(anchor$case,s$ID)];old<-pa[signature==sig&operator==op];if(all(is.finite(v)))add(paste('protein point-estimate identity',op),'CPTAC',sig,cor(v,anchor$protein_anchor,method='spearman'),old$spearman)}
o1<-pa[signature=='Angio36'&operator=='O1',spearman];add('required protein correlation rounded 0.740','CPTAC','Angio36',round(o1,3),.740,0)
add('protein complete cases','CPTAC','both',nrow(anchor),78,0)
q<-rbindlist(qa);wr(q,'REPRODUCTION_CHECKS.csv');wr(counts,'GENOTYPE_COUNTS.csv')
stopifnot(all(q$pass))
saveRDS(list(scores=scores,gt=gt),file.path(SCR,'verified_inputs.rds'))
writeLines(capture.output(sessionInfo()),'sessionInfo.txt');cat('PASS:',nrow(q),'reproduction checks before F computation\n')
