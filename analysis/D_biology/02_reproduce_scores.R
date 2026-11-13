# Author: Chunhao Guo
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(splines);library(jsonlite)})
setDTthreads(1);set.seed(20261113)
P<-normalizePath(commandArgs(TRUE)[1]);O<-file.path(P,'84_JCO_PO_Upgrade_20261113/05_analysis/D_biology');R<-file.path(O,'results');wr<-function(x,n)fwrite(x,file.path(R,n),na='NA')
ob<-readRDS(file.path(P,'48_Polan_ccRCC_GenomeMedicine_20260922/local_only/POLAN_INPUTS.rds'))
legacy<-readRDS(file.path(P,'30_ccRCC_A_GenomeMedicine_final_v1_0/local_only/IMPLEMENTATION_SCORES.rds'))
mem<-fread(file.path(P,'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv'));genes<-lapply(c('Angio16','Angio36'),function(t)mem[get(t)==TRUE,gene]);names(genes)<-c('Angio16','Angio36')
pr<-function(x)100*(rank(x,ties.method='average')-.5)/length(x)
score<-function(X,tr=seq_len(ncol(X))){mu<-rowMeans(X[,tr,drop=FALSE]);ss<-apply(X[,tr,drop=FALSE],1,sd);stopifnot(all(is.finite(ss)),all(ss>1e-8));r0<-colMeans(X);r1<-colMeans((X-mu)/ss);w0<-rep(1/(nrow(X)*sd(r0[tr])),nrow(X));w1<-1/(nrow(X)*ss*sd(r1[tr]));list(I0=(r0-mean(r0[tr]))/sd(r0[tr]),I1=(r1-mean(r1[tr]))/sd(r1[tr]),w0=w0,w1=w1,mu=mu,sd=ss)}
mat<-lapply(ob$inputs,function(x)x$M[c('Angio16','Angio36')]);ids<-lapply(ob$inputs,function(x)x$d$ID)
cb<-file.path(P,'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized');cm<-fread(file.path(cb,'RNA_LINKED_CLINICAL.csv'))[Cohort=='CM-025'];cx<-fread(file.path(cb,'TARGET_GENE_EXPRESSION.csv'));xx<-as.matrix(cx[,-1]);rownames(xx)<-cx[[1]];mat$CheckMate025<-lapply(genes,function(g)xx[g,cm$RNA_ID,drop=FALSE]);ids$CheckMate025<-cm$SUBJID
checks<-list();patients<-list();foldlog<-list()
cv<-fread(file.path(P,'48_Polan_ccRCC_GenomeMedicine_20260922/local_only/HELDOUT_PATIENT_RANKS.csv.gz'))[rep==1&comparison%in%paste('Algorithm',names(genes))]
for(co in names(mat))for(t in names(genes)){
 X<-mat[[co]][[t]];stopifnot(setequal(rownames(X),genes[[t]]),length(genes[[t]])==nrow(X),all(is.finite(X)))
 mat[[co]][[t]]<-X;s<-score(X);r0<-pr(s$I0);r1<-pr(s$I1)
 if(co!='CheckMate025'){
  ls<-legacy[[paste(co,t,sep='__')]]$z;err<-max(abs(s$I0-ls$I0),abs(s$I1-ls$I1));stopifnot(err<1e-8)
  old<-as.data.table(ob$ranks)[trial==co&comparison==paste('Algorithm',t)];old<-old[match(ids[[co]],ID)];er<-max(abs(pr(ls$I0)-old$r0),abs(pr(ls$I1)-old$r1));stopifnot(er<1e-8)
  # Use the existing saved scores and ranks for the actual downstream analysis.
  s$I0<-ls$I0;s$I1<-ls$I1;r0<-old$r0;r1<-old$r1
  dd<-cv[trial==co&comparison==paste('Algorithm',t)];dd<-dd[match(ids[[co]],ID)];held<-dd$switch;direction<-as.character(dd$transition);fold<-dd$fold;classification<-'nominated saved training-frozen partition 1'
 }else{
  err<-er<-NA_real_;fp<-file.path(P,'84_JCO_PO_Upgrade_20261113/05_analysis/shared/CM025_folds.csv')
  if(file.exists(fp)){
   f<-fread(fp); idcol<-intersect(c('ID','SUBJID','patient_id','patient'),names(f))[1]; fc<-intersect(c('fold','Fold'),names(f))[1];stopifnot(!is.na(idcol),!is.na(fc));fold<-f[[fc]][match(ids[[co]],f[[idcol]])];stopifnot(!anyNA(fold),setequal(unique(fold),1:5));h0<-h1<-rep(NA,ncol(X))
   for(k in 1:5){tr<-which(fold!=k);ee<-which(fold==k);z<-score(X,tr);h0[ee]<-z$I0[ee]>=median(z$I0[tr]);h1[ee]<-z$I1[ee]>=median(z$I1[tr])}
   held<-h0!=h1;direction<-ifelse(h0,ifelse(h1,'Stable high','High to low'),ifelse(h1,'Low to high','Stable low'));classification<-'B2 shared treatment-stratified training-frozen folds'
  }else{fold<-rep(NA_integer_,ncol(X));h0<-s$I0>=median(s$I0);h1<-s$I1>=median(s$I1);held<-h0!=h1;direction<-ifelse(h0,ifelse(h1,'Stable high','High to low'),ifelse(h1,'Low to high','Stable low'));classification<-'descriptive full-cohort median groups; B2 folds absent at run'}
 }
 h0<-s$I0>=median(s$I0);h1<-s$I1>=median(s$I1)
 d<-data.table(cohort=co,signature=t,ID=ids[[co]],I0=s$I0,I1=s$I1,r0=r0,r1=r1,delta_r=r1-r0,m_r=(r0+r1)/2,abs_delta=abs(r1-r0),switch_full=h0!=h1,transition_full=ifelse(h0,ifelse(h1,'Stable high','High to low'),ifelse(h1,'Low to high','Stable low')),switch=held,transition=direction,fold=fold,classification=classification)
 patients[[paste(co,t)]]<-d
 checks[[length(checks)+1]]<-data.table(check='saved scores and ranks / full-list recomputation',cohort=co,signature=t,n=nrow(d),estimate=median(d$abs_delta),expected=if(t=='Angio36'&&co%in%c('JAVELIN','IMmotion151'))c(JAVELIN=3.58,IMmotion151=6.08)[co]else NA_real_,max_score_error=err,max_rank_error=er,status=if(co=='CheckMate025')'computed'else'PASS')
 foldlog[[length(foldlog)+1]]<-data.table(cohort=co,signature=t,classification=classification,n=nrow(d),reclassified=sum(held),full_cohort_reclassified=sum(d$switch_full))
}
pa<-rbindlist(patients);qc<-rbindlist(checks);stopifnot(all(round(qc[signature=='Angio36'&cohort%in%c('JAVELIN','IMmotion151'),estimate],2)==c(3.58,6.08)))
# Reproduce point clinical coefficients before carrying forward the saved bootstrap intervals.
phdir<-file.path(P,'70_npj_Visual_Narrative_Final_20261011/figures/source_data/sc/new');oldc<-fread(file.path(phdir,'PRIMARY_PHENOTYPE_COEFFICIENTS.csv'));oldt<-fread(file.path(phdir,'PRIMARY_PHENOTYPE_TESTS.csv'));oldsec<-fread(file.path(P,'48_Polan_ccRCC_GenomeMedicine_20260922/new_analysis/aggregate/SECONDARY_PHENOTYPE_TESTS.csv'))
repcl<-list()
for(t in names(genes)){
 d<-pa[cohort=='IMmotion151'&signature==t];dd<-cbind(ob$meta,as.data.frame(d[,.(delta_r,m_r)]));f<-lm(delta_r~ns(m_r,3)+IMDC+sarcomatoid+specimen+age+sex,dd)
 for(k in oldc[target==t,coefficient]){z<-oldc[target==t&coefficient==k];b<-coef(f)[k];stopifnot(abs(b-z$estimate)<1e-8);repcl[[length(repcl)+1]]<-data.table(check='IMmotion151 clinical point coefficient',cohort='IMmotion151',signature=t,coefficient=k,n=nobs(f),estimate=b,expected=z$estimate,max_score_error=abs(b-z$estimate),status='PASS')}
}
wr(rbindlist(list(qc,rbindlist(repcl)),fill=TRUE),'REPRODUCTION_CHECKS.csv');wr(oldc,'CARRIED_FORWARD_CLINICAL_COEFFICIENTS.csv');wr(oldt,'CARRIED_FORWARD_CLINICAL_TESTS.csv');wr(oldsec,'CARRIED_FORWARD_CLINICAL_SECONDARY_TESTS.csv');wr(pa,'PATIENT_SCORES_RANKS.csv');wr(rbindlist(foldlog),'CLASSIFICATION_PROVENANCE.csv')
saveRDS(list(mat=mat,ids=ids,patients=pa,genes=genes,meta=ob$meta,folds=ob$folds),file.path(R,'qualified_inputs.rds'))
writeLines(capture.output(sessionInfo()),file.path(O,'sessionInfo.txt'));print(qc);print(rbindlist(repcl)[signature=='Angio36'])
