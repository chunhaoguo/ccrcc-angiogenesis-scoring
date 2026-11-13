# Author / creator / last modified by: Chunhao Guo
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival)})
setDTthreads(1);set.seed(20261113)
ROOT<-'/LOCAL_PROJECT_ROOT'
WORK<-file.path(ROOT,'84_JCO_PO_Upgrade_20261113');OUT<-file.path(WORK,'05_analysis/B2_independent_cohorts')
SCR<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2'
shared<-file.path(WORK,'05_analysis/shared/operators.R');stopifnot(file.exists(shared));operator_hash<-digest::digest(file=shared,algo='sha256');source(shared);stopifnot(identical(operator_hash,digest::digest(file=shared,algo='sha256')))
snapshot<-file.path(OUT,'scripts/operators_used.R')
if(!file.exists(snapshot))file.copy(shared,snapshot)else stopifnot(identical(operator_hash,digest::digest(file=snapshot,algo='sha256')))
sets<-fromJSON(file.path(OUT,'results/gene_sets.json'));sets<-sets[c('Angio16','Angio36')]
cohort<-commandArgs(TRUE)[1];z<-readRDS(file.path(SCR,paste0(cohort,'.rds')));X<-z$x;patient<-z$patient
wr<-function(d,n)fwrite(d,file.path(OUT,'results',paste0(cohort,'_',n,'.csv')))
ci<-function(x)if(any(is.finite(x)))unname(quantile(x[is.finite(x)],c(.025,.975),type=7))else c(NA_real_,NA_real_)
exact<-function(k,n)if(n>0)as.numeric(binom.test(k,n)$conf.int)else c(NA_real_,NA_real_)
metrics<-function(S,pid,tag='primary'){
 n<-nrow(S);ops<-colnames(S);av<-vapply(seq_len(ncol(S)),function(k)all(is.finite(S[,k])),logical(1));valid<-S[,av,drop=FALSE]
 calls<-sweep(valid,2,apply(valid,2,median),'>=');non<-rowSums(calls)>0&rowSums(calls)<ncol(calls)
 grp<-split(seq_len(n),pid);B<-2000L;bootdis<-array(NA_real_,dim=c(B,ncol(valid),ncol(valid)));bootnon<-numeric(B)
 set.seed(20261113)
 for(b in seq_len(B)){
  ids<-sample.int(length(grp),length(grp),replace=TRUE);ii<-unlist(grp[ids],use.names=FALSE);v<-valid[ii,,drop=FALSE]
  cc<-sweep(v,2,apply(v,2,median),'>=');nb<-nrow(cc);ones<-colSums(cc);inter<-crossprod(cc)
  bootdis[b,,]<-(outer(ones,ones,'+')-2*inter)*100/nb
  bootnon[b]<-mean(rowSums(cc)>0&rowSums(cc)<ncol(cc))*100
 }
 pairs<-list()
 for(a in seq_along(ops))for(b in seq_along(ops)){
  if(av[a]&&av[b]){
   ia<-match(ops[a],colnames(valid));ib<-match(ops[b],colnames(valid));dif<-abs(rank(S[,a])-rank(S[,b]))*100/n
   conf<-ci(bootdis[,ia,ib]);value<-mean(calls[,ia]!=calls[,ib])*100
   rho<-suppressWarnings(cor(S[,a],S[,b],method='spearman'));med<-median(dif);p90<-unname(quantile(dif,.9))
  }else{conf<-c(NA_real_,NA_real_);value<-rho<-med<-p90<-NA_real_}
  pairs[[length(pairs)+1]]<-data.table(label=tag,cohort,operator_a=ops[a],operator_b=ops[b],n_total=n,n_patients=length(grp),n_complete=if(av[a]&&av[b])n else 0,
   disagreement_percent=value,lower=conf[1],upper=conf[2],spearman=rho,median_absolute_rank_displacement_pp=med,p90_absolute_rank_displacement_pp=p90,
   bootstrap='2000 fixed-score paired resamples; medians recomputed; patient-clustered for regions',operators_recomputed=FALSE)
 }
 conf<-ci(bootnon)
 list(pair=rbindlist(pairs),non=data.table(label=tag,cohort,n_total=n,n_patients=length(grp),available_operators=paste(ops[av],collapse=';'),n_available_operators=sum(av),nonunanimous_n=if(sum(av)>=2)sum(non)else NA_integer_,nonunanimous_percent=if(sum(av)>=2)mean(non)*100 else NA_real_,lower=if(sum(av)>=2)conf[1]else NA_real_,upper=if(sum(av)>=2)conf[2]else NA_real_),
  calls=data.table(sample=rownames(S),patient=pid,as.data.table(calls),nonunanimous=non))
}
scoreall<-function(x,tag){
 statuses<-coverage<-scores<-pairs<-non<-calls<-list()
 for(sig in names(sets)){
  g<-sets[[sig]];missing<-setdiff(g,rownames(x));coverage[[sig]]<-data.table(cohort,signature=sig,required=length(g),measured=sum(g%in%rownames(x)),missing=paste(missing,collapse=';'),complete=!length(missing))
  S<-matrix(NA_real_,nrow=ncol(x),ncol=7,dimnames=list(colnames(x),paste0('O',1:7)))
  for(op in colnames(S)){
   t<-proc.time()[3];r<-get(op)(x,g,whole_transcriptome=TRUE);S[,op]<-r$scores
   statuses[[length(statuses)+1]]<-data.table(cohort,signature=sig,operator=op,n_total=ncol(x),n_available=sum(r$available),reason=r$reason,seconds=proc.time()[3]-t,implementation=normalizePath(shared),implementation_sha256=operator_hash,label=tag)
   cat(cohort,tag,sig,op,'available',sum(r$available),'of',ncol(x),'seconds',proc.time()[3]-t,'\n');flush.console();gc()
  }
  scores[[sig]]<-S;m<-metrics(S,if(tag=='sensitivity')sub('[-_].*$','',colnames(x)) else patient,tag)
  m$pair[,signature:=sig];m$non[,signature:=sig];m$calls[,signature:=sig];pairs[[sig]]<-m$pair;non[[sig]]<-m$non;calls[[sig]]<-m$calls
 }
 list(scores=scores,status=rbindlist(statuses),coverage=rbindlist(coverage),pair=rbindlist(pairs),non=rbindlist(non),calls=rbindlist(calls,fill=TRUE))
}
r<-scoreall(X,'primary');wr(r$status,'operator_status');wr(r$coverage,'signature_coverage');wr(r$pair,'pairwise_metrics');wr(r$non,'nonunanimous');wr(r$calls,'research_groups')
qc<-rbindlist(lapply(names(sets),function(s){g<-sets[[s]];m<-X[g,,drop=FALSE];data.table(signature=s,gene=g,reference_SD=apply(m,1,sd),reference_unscaled_MAD=apply(m,1,mad,constant=1),minimum=apply(m,1,min),fraction_at_minimum=apply(m,1,function(v)mean(v==min(v))))}));wr(qc,'signature_reference_QC')
long<-rbindlist(lapply(names(r$scores),function(s)data.table(signature=s,sample=colnames(X),patient,as.data.table(r$scores[[s]]))))
wr(long,'scores');saveRDS(r$scores,file.path(OUT,'results',paste0(cohort,'_scores.rds')))
if(cohort=='TRACERx'){
 count<-table(patient);multi<-names(count)[count>=2];regional<-variances<-list()
 for(sig in names(sets))for(op in paste0('O',1:7)){
  v<-r$scores[[sig]][,op];available<-all(is.finite(v));g<-if(available)v>=median(v) else rep(NA,length(v));flags<-vapply(multi,function(p){u<-g[patient==p];if(anyNA(u))NA else any(u)&any(!u)},logical(1));k<-sum(flags,na.rm=TRUE);n<-sum(!is.na(flags));conf<-exact(k,n)
  regional[[length(regional)+1]]<-data.table(label='exploratory',cohort,signature=sig,operator=op,n_patients_total=length(count),n_multiregion=length(multi),n_available_patients=n,discordant_patients=if(n)k else NA,discordance_percent=if(n)100*k/n else NA,lower=100*conf[1],upper=100*conf[2])
  if(available){fit<-lme4::lmer(v~1+(1|patient),REML=TRUE);vc<-as.data.frame(lme4::VarCorr(fit));between<-vc$vcov[vc$grp=='patient'];within<-sigma(fit)^2;icc<-between/(between+within)}else between<-within<-icc<-NA_real_
  variances[[length(variances)+1]]<-data.table(label='exploratory',cohort,signature=sig,operator=op,n_regions=length(v),n_patients=length(count),between_patient_variance=between,within_patient_variance=within,ICC=icc,method='REML random-intercept model; all tumour regions; patient grouping')
 }
 wr(rbindlist(regional),'within_patient_discordance');wr(rbindlist(variances),'variance_ICC')
 set.seed(20261113);sel<-vapply(split(seq_along(patient),patient),function(ii)ii[sample.int(length(ii),1)],integer(1));wr(data.table(patient=patient[sel],sample=colnames(X)[sel]),'one_region_selection')
 one<-scoreall(X[,sel,drop=FALSE],'sensitivity');wr(one$status,'one_region_operator_status');wr(one$pair,'one_region_pairwise_metrics');wr(one$non,'one_region_nonunanimous')
}
if(cohort=='CPTAC'){
 readmat<-function(f){d<-fread(file.path(ROOT,'43_evidence_linked_exploration_20260918/protein_inputs',f));m<-as.matrix(d[,-1]);rownames(m)<-d[[1]];m}
 p<-fread(file.path(ROOT,'43_evidence_linked_exploration_20260918/protein_inputs/CASE_PAIRING.csv'))[eligible==TRUE];stopifnot(nrow(p)==78,all(p$case%in%colnames(X)))
 protein<-readmat('PROTEIN_MVD4.csv')[,p$case,drop=FALSE];anchor<-colMeans(t(scale(t(protein))));ix<-match(p$case,colnames(X));anchors<-list()
 set.seed(20261113);boots<-replicate(2000,sample.int(78,78,replace=TRUE),simplify=FALSE)
 for(sig in names(sets))for(op in paste0('O',1:7)){
  v<-r$scores[[sig]][ix,op];base<-r$scores[[sig]][ix,'O1'];ok<-all(is.finite(v));est<-if(ok)cor(v,anchor,method='spearman')else NA_real_
  bb<-if(ok)vapply(boots,function(ii)cor(v[ii],anchor[ii],method='spearman'),numeric(1))else NA_real_
  bd<-if(ok)vapply(boots,function(ii)cor(v[ii],anchor[ii],method='spearman')-cor(base[ii],anchor[ii],method='spearman'),numeric(1))else NA_real_
  conf<-ci(bb);dc<-ci(bd);anchors[[length(anchors)+1]]<-data.table(label='exploratory',cohort,signature=sig,operator=op,n=78,spearman=est,lower=conf[1],upper=conf[2],difference_vs_O1=if(ok)est-cor(base,anchor,method='spearman')else NA_real_,difference_lower=dc[1],difference_upper=dc[2],valid_bootstrap=sum(is.finite(bb)),score_reference_cases=ncol(X),protein_reference_cases=78)
 }
 wr(rbindlist(anchors),'protein_anchor');wr(data.table(case=p$case,protein_anchor=anchor),'protein_proxy')
}
writeLines(capture.output(sessionInfo()),file.path(OUT,paste0('logs/',cohort,'_sessionInfo.txt')))
