# Author / creator / last modified by: Chunhao Guo
# post hoc follow-up (specified after modules B2 and D)
source('scripts/00_common.R')
stopifnot(all(rd('results/REPRODUCTION_CHECKS.csv')$pass))
v<-readRDS(file.path(SCR,'verified_inputs.rds'));scores<-v$scores;gt<-v$gt
ob<-readRDS(file.path(A,'D_biology/results/qualified_inputs.rds'))
spec<-list(IMmotion151=as.data.table(ob$meta)[,.(ID,specimen=as.character(specimen))],
 CheckMate025=rd(file.path(A,'B2_independent_cohorts/results/CM025_sample_manifest.csv'))[,.(ID=patient,specimen)])
spec$CheckMate025[is.na(specimen)|!nzchar(trimws(specimen)),specimen:=NA_character_]
wr(rbindlist(lapply(cohorts,function(co){s<-spec[[co]];data.table(cohort=co,n_total=length(ob$ids[[co]]),n_specimen_observed=if(is.null(s))0L else sum(!is.na(s$specimen)),reason=if(is.null(s))'no qualified origin covariate in reused B1/D inputs'else'qualified existing source field; observed complete cases only')})),'SPECIMEN_AVAILABILITY.csv')

pooled_sd<-function(Y,g){
 n1<-sum(g);n0<-sum(!g);if(n1<2||n0<2)return(rep(NA_real_,ncol(Y)))
 sqrt(((n1-1)*apply(Y[g,,drop=FALSE],2,var)+(n0-1)*apply(Y[!g,,drop=FALSE],2,var))/(n1+n0-2))
}
effect<-function(Y,g,X=NULL){
 den<-pooled_sd(Y,g)
 if(is.null(X))b<-colMeans(Y[g,,drop=FALSE])-colMeans(Y[!g,,drop=FALSE])
 else {f<-lm.fit(X,Y);if(f$rank<ncol(X))return(rep(NA_real_,ncol(Y)));b<-f$coefficients[2,]}
 as.numeric(b/den)
}
run_genotype<-function(S,co,sig,gene,mode='main',group=NA_character_,module='F1'){
 gname<-gene;d<-gt[cohort==co&gene==gname]
 d<-d[match(rownames(S),ID)];stopifnot(!anyNA(d$ID),!anyDuplicated(d$ID))
 nunknown<-sum(d$status=='Unknown_not_assayed');keep<-d$status%in%c('Altered','No_reported_alteration')
 origin<-NULL;reason<-'';missing_origin<-0L
 if(mode=='specimen_adjusted'){
  s<-spec[[co]]
  if(is.null(s))reason<-'not run: no qualified specimen-origin covariate'
  else{origin<-s$specimen[match(d$ID,s$ID)];missing_origin<-sum(keep&is.na(origin));keep<-keep&!is.na(origin)}
 }
 dd<-d[keep];Y<-S[keep,,drop=FALSE];altered<-dd$status=='Altered';n1<-sum(altered);n0<-sum(!altered)
 good<-apply(Y,2,function(x)all(is.finite(x)))
 theta<-rep(NA_real_,ncol(S));names(theta)<-colnames(S);B<-2000L
 bs<-matrix(NA_real_,B,ncol(S),dimnames=list(NULL,colnames(S)))
 X<-if(mode=='specimen_adjusted'&&!nzchar(reason))model.matrix(~ altered+factor(origin[keep]))else NULL
 if(!nzchar(reason)&&n1>=2&&n0>=2&&any(good)){
  yy<-Y[,good,drop=FALSE];theta[good]<-effect(yy,altered,X)
  set.seed(20261113)
  for(b in seq_len(B)){
   ii<-sample.int(nrow(Y),replace=TRUE)
   bs[b,good]<-effect(yy[ii,,drop=FALSE],altered[ii],if(is.null(X))NULL else X[ii,,drop=FALSE])
  }
 }
 base<-if(module=='F1')'O1'else paste0(group,'__O1');baseix<-match(base,colnames(S));stopifnot(!is.na(baseix))
 res<-rbindlist(lapply(seq_len(ncol(S)),function(j){
  ci<-interval(bs[,j]);delta<-theta[j]-theta[baseix];db<-bs[,j]-bs[,baseix];dc<-interval(db)
  why<-if(nzchar(reason))reason else if(!good[j])'unavailable: inherited complete-list operator failure'else if(!is.finite(theta[j]))'not estimable'else'estimable'
  data.table(module,anchor=gene,cohort=co,signature=sig,gene_group=group,operator=if(module=='F1')colnames(S)[j]else sub('^.*__','',colnames(S)[j]),
   estimand=if(mode=='main')'Cohen d altered minus source-negative'else'specimen-adjusted altered coefficient / pooled within-genotype SD',
   role=if(mode=='main')'main estimand; exploratory'else'sensitivity; exploratory',mode,n_total=nrow(S),n=nrow(Y),n_altered=n1,n_source_negative=n0,n_unknown_excluded=nunknown,n_missing_specimen_excluded=missing_origin,
   estimate=theta[j],lower=ci[1],upper=ci[2],difference_vs_O1=delta,difference_lower=dc[1],difference_upper=dc[2],p_difference=if(j==baseix&&is.finite(delta))1 else bp(db,delta),
   valid_bootstrap=sum(is.finite(bs[,j])),valid_difference_bootstrap=sum(is.finite(db)),status=why,unit='standardised mean difference',bootstrap='2000 paired patient draws; fixed scores; type-7 percentile CI')
 }))
 saveRDS(list(estimate=theta,bootstrap=bs,IDs=dd$ID),file.path(SCR,paste0(module,'_',co,'_',sig,'_',gene,'_',mode,if(!is.na(group))paste0('_',gsub('[^A-Za-z0-9]','',group))else'','.rds')))
 res
}
F1<-list()
for(co in cohorts)for(sig in sigs){
 s<-dcast(scores[cohort==co&signature==sig],ID~operator,value.var='score');S<-as.matrix(s[,..ops]);rownames(S)<-s$ID
 for(gene in c('VHL','PBRM1','BAP1'))for(mode in c('main','specimen_adjusted')){
  cat('F1',co,sig,gene,mode,'\n');flush.console();F1[[length(F1)+1]]<-run_genotype(S,co,sig,gene,mode)
 }
}
f1<-rbindlist(F1)
f1[,q_BH:=p.adjust(p_difference,method='BH',n=21),by=.(cohort,signature,mode)]
f1[,BH_family:='trial x signature x mode; 7 operators x 3 genes; 21 planned slots including self/unavailable']
wr(f1[mode=='main'],'F1_GENOTYPE_ALL.csv');wr(f1[mode=='specimen_adjusted'],'F1_SPECIMEN_SENSITIVITY_ALL.csv')

# F2: exact carry-forward, no new full-signature protein resampling.
p<-rd(file.path(A,'B2_independent_cohorts/results/CPTAC_protein_anchor.csv'))
wr(p,'F2_PROTEIN_CARRIED_FORWARD.csv')
f2<-p[,.(module='F2',anchor='Vascular protein proxy',cohort,signature,gene_group=NA_character_,operator,
 estimand='Spearman rho with four-protein proxy',role='carried forward; exploratory',mode='main',n_total=score_reference_cases,n,
 estimate=spearman,lower,upper,difference_vs_O1,difference_lower,difference_upper,p_difference=NA_real_,q_BH=NA_real_,
 valid_bootstrap,valid_difference_bootstrap=valid_bootstrap,status=ifelse(is.finite(spearman),'carried forward unchanged','unavailable: inherited complete-list operator failure'),unit='Spearman rho',
 q_status='not available in B2 source; no recomputation under F2',source=file.path(A,'B2_independent_cohorts/results/CPTAC_protein_anchor.csv'))]
wr(rbindlist(list(f1[mode=='main'],f2),fill=TRUE),'F1_F2_COMBINED_ANCHORS.csv')

# Fixed manuscript group definitions, reused from module D, cross-checked against membership.
w<-rd(file.path(A,'D_biology/results/CARRIED_FORWARD_ALL_GENE_WEIGHTS.csv'))
groups<-unique(w[cohort=='JAVELIN',.(gene,group)])
mem<-rd(file.path(ROOT,'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv'))
groups<-merge(groups,mem[,.(gene,Angio16,Angio36)],by='gene',all.x=TRUE)
groupnames<-c('Shared 6','Angio16-only','Angio36-only FGF/FGFR','Angio36-only other')
stopifnot(setequal(groups$group,groupnames),nrow(groups)==46,!anyNA(groups$Angio36),!anyDuplicated(groups$gene))
stopifnot(all(groups[group=='Shared 6',Angio16&Angio36]),all(groups[group=='Angio16-only',Angio16&!Angio36]),all(groups[grepl('Angio36-only',group),!Angio16&Angio36]))
stopifnot(identical(as.integer(table(factor(groups$group,levels=groupnames))),c(6L,10L,26L,4L)))
wr(groups,'F3_FIXED_GENE_GROUPS.csv')
allF3<-allweights<-allgroupscores<-list();protein<-rd(file.path(A,'B2_independent_cohorts/results/CPTAC_protein_proxy.csv'))
for(co in c(cohorts,'CPTAC')){
 if(co=='CPTAC'){
  z<-readRDS('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2/CPTAC.rds');X<-z$x[groups$gene,,drop=FALSE];ids<-z$patient;rm(z);gc()
 }else{
  ids<-ob$ids[[co]];m1<-ob$mat[[co]]$Angio16;m2<-ob$mat[[co]]$Angio36
  X<-rbind(m1,m2[!rownames(m2)%in%rownames(m1),,drop=FALSE]);X<-X[groups$gene,,drop=FALSE]
 }
 stopifnot(ncol(X)==length(ids),!anyDuplicated(ids),all(is.finite(X)),setequal(rownames(X),groups$gene))
 colnames(X)<-ids
 # Group O1/O2 references use the full RNA cohort, never the anchor-selected subset.
 for(gr in groupnames){
  g<-groups[group==gr,gene];S<-matrix(NA_real_,ncol(X),2,dimnames=list(ids,paste0(gr,'__',c('O1','O2'))))
  for(j in 1:2)S[,j]<-get(paste0('O',j))(X,g)$scores
  allgroupscores[[length(allgroupscores)+1]]<-data.table(cohort=co,ID=ids,gene_group=gr,O1=S[,1],O2=S[,2],reference_n=ncol(X))
  if(co!='CPTAC')for(gene in c('VHL','PBRM1')){
   cat('F3',co,gr,gene,'\n');flush.console();allF3[[length(allF3)+1]]<-run_genotype(S,co,'Angio36 fixed-group decomposition',gene,group=gr,module='F3')
  }else{
   ix<-match(protein$case,ids);stopifnot(!anyNA(ix));Y<-S[ix,,drop=FALSE];a<-protein$protein_anchor
   est<-apply(Y,2,function(x)cor(x,a,method='spearman'));B<-2000L;bs<-matrix(NA_real_,B,2)
   set.seed(20261113);for(b in 1:B){ii<-sample.int(length(a),replace=TRUE);bs[b,]<-apply(Y[ii,,drop=FALSE],2,function(x)cor(x,a[ii],method='spearman'))}
   db<-bs-bs[,1];for(j in 1:2){ci<-interval(bs[,j]);dc<-interval(db[,j]);delta<-est[j]-est[1]
    allF3[[length(allF3)+1]]<-data.table(module='F3',anchor='Vascular protein proxy',cohort=co,signature='Angio36 fixed-group decomposition',gene_group=gr,operator=paste0('O',j),estimand='Spearman rho with four-protein proxy',role='exploratory gene-group anchor',mode='main',n_total=ncol(X),n=nrow(Y),estimate=est[j],lower=ci[1],upper=ci[2],difference_vs_O1=delta,difference_lower=dc[1],difference_upper=dc[2],p_difference=if(j==1)1 else bp(db[,j],delta),valid_bootstrap=sum(is.finite(bs[,j])),valid_difference_bootstrap=sum(is.finite(db[,j])),status='estimable',unit='Spearman rho',bootstrap='2000 paired case draws; fixed 109-case score reference; type-7 percentile CI')
   };saveRDS(list(estimate=est,bootstrap=bs,IDs=protein$case),file.path(SCR,paste0('F3_CPTAC_',gsub('[^A-Za-z0-9]','',gr),'.rds')))
  }
 }
 # Full-signature coefficients explain arithmetic upweighting, separately from group-standardized scores.
 m<-X[mem[Angio36==TRUE,gene],,drop=FALSE];sdgene<-apply(m,1,sd);r0<-colMeans(m);r2<-colMeans(sweep(sweep(m,1,rowMeans(m),'-'),1,sdgene,'/'))
 W0<-1/(nrow(m)*sd(r0));W2<-1/(nrow(m)*sdgene*sd(r2))
 ww<-copy(groups);ww[,`:=`(cohort=co,reference_n=ncol(X),reference_sd=apply(X,1,sd)[match(gene,rownames(X))],weight_full_Angio36_O1=ifelse(Angio36,W0,NA_real_),weight_full_Angio36_O2=W2[match(gene,names(W2))])]
 ww[,weight_ratio_full_Angio36_O2_O1:=weight_full_Angio36_O2/weight_full_Angio36_O1];allweights[[co]]<-ww
 rm(X);gc()
}
f3<-rbindlist(allF3,fill=TRUE)
f3[,q_BH:=NA_real_]
f3[operator=='O2',q_BH:=p.adjust(p_difference,'BH',n=4),by=.(cohort,anchor)]
f3[operator=='O1',q_BH:=1]
f3[,BH_family:='descriptive; cohort x anchor; four O2-minus-O1 group contrasts; separate from F1']
wr(f3,'F3_GROUP_ANCHORS_ALL.csv');ww<-rbindlist(allweights);wr(ww,'F3_ALL_GENE_WEIGHTS.csv')
wg<-ww[,.(n_genes=.N,n_Angio36=sum(Angio36),median_gene_SD=median(reference_sd),sum_full36_O1_weight=if(any(Angio36))sum(weight_full_Angio36_O1,na.rm=TRUE)else NA_real_,sum_full36_O2_weight=if(any(Angio36))sum(weight_full_Angio36_O2,na.rm=TRUE)else NA_real_,median_full36_weight_ratio=if(any(Angio36))median(weight_ratio_full_Angio36_O2_O1,na.rm=TRUE)else NA_real_),by=.(cohort,group)]
wg[,sum_full36_weight_ratio:=sum_full36_O2_weight/sum_full36_O1_weight];wr(wg,'F3_GROUP_WEIGHT_SUMMARY.csv')
wr(rbindlist(allgroupscores),'F3_GROUP_SCORES_LOCAL_ONLY.csv')
writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
cat('COMPLETE F1 main',nrow(f1[mode=='main']),'sensitivity',nrow(f1[mode=='specimen_adjusted']),'F2',nrow(f2),'F3',nrow(f3),'\n')
