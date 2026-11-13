# Author: Chunhao Guo
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(splines);library(sandwich);library(jsonlite)})
setDTthreads(1);set.seed(20261113);P<-normalizePath(commandArgs(TRUE)[1]);O<-file.path(P,'84_JCO_PO_Upgrade_20261113/05_analysis/D_biology');R<-file.path(O,'results');wr<-function(x,n)fwrite(x,file.path(R,n),na='NA');rd<-function(n)fread(file.path(R,n),na.strings=c('NA',''))
ob<-readRDS(file.path(R,'qualified_inputs.rds'));pa<-copy(ob$patients);gs<-c('VHL','PBRM1','BAP1','SETD2','KDM5C')
j<-rd('source_JAVELIN_genotype.csv');jw<-rd('source_JAVELIN_WES_eligibility.csv');im<-rd('source_IMmotion151_genotype_cluster.csv');cm<-rd('source_CheckMate_WES_crosswalk.csv');maf<-rd('source_CheckMate_target_MAF.csv');mafid<-rd('source_CheckMate_MAF_assayed_ids.csv')$MAF_Tumor_ID
stopifnot(!anyDuplicated(j$ID),!anyDuplicated(im$PATIENT_ID),!anyDuplicated(cm$SUBJID))
excluded<-c('Silent','Intron','IGR',"3'UTR","5'UTR","3'Flank","5'Flank",'RNA','lincRNA');coding<-maf[!Variant_Classification%in%excluded]
allg<-list();audit<-list()
for(co in unique(pa$cohort)){
 ids<-unique(pa[cohort==co,ID]);dd<-if(co=='JAVELIN')j[match(ids,ID)]else if(co=='IMmotion151')im[match(ids,PATIENT_ID)]else cm[match(ids,SUBJID)]
 for(g in gs){
  raw<-if(co=='IMmotion151')dd[[paste0('FMI_',g)]]else dd[[g]]
  if(co=='JAVELIN'){ass<-raw%in%c(0,1);pos<-raw==1;assay<-'WES';def<-'S21 source binary 1 mutated / 0 not mutated; mutation class/copy-number scope is source-defined';spec<-file.path(P,'09_fulltext/downloads/CCRCC-00889/CCRCC-00889_Motzer_2020_Nature_Medicine_JAVELIN_biomarker_tables.xlsx#S21_Gene_mutation_status');reason<-ifelse(ass,'qualified S21 explicit flag','no S21 assay/call')}
  else if(co=='IMmotion151'){ass<-!is.na(dd$FMI_SAMPLE_ID)&nzchar(dd$FMI_SAMPLE_ID);pos<-!is.na(raw)&nzchar(raw);stopifnot(!any(!ass&pos));assay<-'FoundationOne targeted FMI panel';def<-'Any reported FMI alteration string positive; assayed FMI_SAMPLE_ID with empty gene field = no reported alteration; no FMI ID = unknown';spec<-file.path(P,'57_GenomeMedicine_Final_20260927/analysis/source_documents/NIHMS1733314-supplement-Table_S1.xlsx#TableS1');reason<-ifelse(ass,'FMI_SAMPLE_ID present','no FMI_SAMPLE_ID')}
  else{ass<-!is.na(dd$MAF_Tumor_ID)&dd$MAF_Tumor_ID%in%mafid;pos<-dd$MAF_Tumor_ID%in%coding[Hugo_Symbol==g,Tumor_Sample_Barcode];assay<-'WES somatic MAF';def<-'S2 MAF coding nonsynonymous somatic variant; exclude silent/intronic/UTR/flank/RNA; S1 tumor-MAF ID must be present in S2; not a CNA call';spec<-file.path(P,'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized/source/Braun2020_supplement.xlsx#S2_WES_Data+S1_Clinical_and_Immune_Data');reason<-ifelse(ass,'S1 MAF tumor ID verified in S2','no qualified S1-to-S2 MAF ID')}
  status<-ifelse(!ass,'Unknown_not_assayed',ifelse(pos,'Altered','No_reported_alteration'))
  allg[[length(allg)+1]]<-data.table(cohort=co,ID=ids,gene=g,status=status,qualified=ass,source_value=as.character(raw),assay_eligibility_reason=reason)
  audit[[length(audit)+1]]<-data.table(cohort=co,gene=g,RNA_n=length(ids),n_qualified=sum(ass),n_altered=sum(status=='Altered'),n_source_negative=sum(status=='No_reported_alteration'),n_unknown_not_assayed=sum(!ass),assay=assay,source=spec,status_definition=def,WES_summary_observed_if_JAVELIN=if(co=='JAVELIN')sum(is.finite(jw$total_reads[match(ids,jw$ID)]))else NA_integer_)
 }
}
gt<-rbindlist(allg);wr(gt,'GENOTYPE_PATIENT_JOIN.csv');wr(rbindlist(audit),'GENOTYPE_JOIN_AUDIT.csv');wr(coding[, .N,by=.(Hugo_Symbol,Variant_Classification)],'CHECKMATE_QUALIFYING_VARIANT_CLASSES.csv')
# Explicit source-flag comparison demonstrates why S1 WT does not mean no nonsynonymous MAF variant.
cmcheck<-list();for(g in gs){z<-cm[!is.na(MAF_Tumor_ID)&MAF_Tumor_ID%in%mafid];z[,maf_coding:=MAF_Tumor_ID%in%coding[Hugo_Symbol==g,Tumor_Sample_Barcode]];z[,source_flag:=as.character(get(g))];cmcheck[[g]]<-z[,.N,by=.(source_flag,maf_coding)][,gene:=g]};wr(rbindlist(cmcheck),'CHECKMATE_MAF_VS_SOURCE_TRUNCATING_FLAGS.csv')
dat<-merge(pa,gt,by=c('cohort','ID'),allow.cartesian=TRUE);dat[,status:=factor(status,levels=c('No_reported_alteration','Altered','Unknown_not_assayed'))]
# Models fit all score-available patients, with unknown genotype a separate factor level.
fitcat<-function(d,y,term='status'){
 d<-copy(d);d[[term]]<-droplevels(d[[term]]);counts<-table(d[[term]]);wanted<-if(term=='status')c('Altered','Unknown_not_assayed')else levels(d[[term]])[-1]
 out<-list();for(level in wanted){
  cn<-paste0(term,level);ref<-levels(d[[term]])[1];n1<-sum(as.character(d[[term]])==level);n0<-sum(as.character(d[[term]])==ref);why<-'';b<-lo<-hi<-p<-NA_real_
  if(n1<5||n0<5)why<-'not estimable: fewer than 5 in comparison category'
  else if(y=='switch'&&any(table(factor(d[[term]][d[[term]]%in%c(ref,level)],levels=c(ref,level)),factor(d$switch[d[[term]]%in%c(ref,level)],levels=c(FALSE,TRUE)))<5))why<-'not estimable: sparse events/non-events (<5)'
  else{
   fo<-as.formula(paste(y,if(y=='switch')'~ns(abs(r0-50),3)+'else'~ns(m_r,3)+',term));f<-if(y=='switch')try(glm(fo,d,family=binomial()),silent=TRUE)else try(lm(fo,d),silent=TRUE)
   if(inherits(f,'try-error')||!cn%in%names(coef(f))||!is.finite(coef(f)[cn]))why<-'not estimable: model rank/fit failure'
   else{V<-try(vcovHC(f,type='HC3'),silent=TRUE);if(inherits(V,'try-error')||!is.finite(V[cn,cn])||V[cn,cn]<=0)why<-'not estimable: covariance failure'else{b<-coef(f)[cn];se<-sqrt(V[cn,cn]);crit<-qt(.975,df.residual(f));lo<-b-crit*se;hi<-b+crit*se;p<-2*pt(-abs(b/se),df.residual(f));if(y=='switch'){b<-exp(b);lo<-exp(lo);hi<-exp(hi)}}}
  }
  out[[length(out)+1]]<-data.table(contrast=paste(level,'vs',ref),n=nrow(d),n_reference=n0,n_category=n1,estimate=b,lower=lo,upper=hi,p=p,unit=if(y=='switch')'odds ratio'else'percentile points',status=if(nzchar(why))why else'estimable')
 };rbindlist(out)
}
models<-list()
for(co in unique(dat$cohort))for(t in c('Angio36','Angio16'))for(g in gs)for(y in c('delta_r','abs_delta','switch')){
 d<-dat[cohort==co&signature==t&gene==g];a<-fitcat(d,y);a[,`:=`(cohort=co,signature=t,gene=g,endpoint=y,label=if(t=='Angio36'&&g%in%gs[1:3]&&y=='delta_r')'primary'else'exploratory',classification=unique(d$classification))];models[[length(models)+1]]<-a
}
res<-rbindlist(models);primary<-res[label=='primary'&grepl('^Altered vs',contrast)];stopifnot(nrow(primary)==9);primary[,`:=`(q_BH=p.adjust(p,'BH',n=9),family='3 genes x 3 qualified cohorts; 9 tests')];wr(primary,'D_PRIMARY_GENOTYPE.csv')
secondary<-res[!(label=='primary'&grepl('^Altered vs',contrast))];secondary[,label:='exploratory'];secondary[,q_BH_secondary:=p.adjust(p,'BH')];secondary[,family:='all remaining genotype contrasts/endpoints; exploratory, separate from primary'];wr(secondary,'D_SECONDARY_GENOTYPE_ALL_ENDPOINTS.csv')
counts<-dat[,{
 k<-sum(switch);ci<-binom.test(k,.N)$conf.int;list(n=.N,mean_signed_pp=mean(delta_r),median_signed_pp=median(delta_r),mean_absolute_pp=mean(abs_delta),reclassified=k,reclassification_fraction=k/.N,reclassification_lower=ci[1],reclassification_upper=ci[2])
},by=.(cohort,signature,gene,status)];wr(counts,'GENOTYPE_CATEGORY_DESCRIPTIVES.csv')
# Independent DNA purity: CM025 ABSOLUTE source estimates, unit = 10 percentage-point increase.
pc<-cm[match(ob$ids$CheckMate025,SUBJID)];pur<-list()
for(t in c('Angio36','Angio16'))for(y in c('delta_r','abs_delta','switch')){
 d<-copy(pa[cohort=='CheckMate025'&signature==t]);d[,purity_10pp:=as.numeric(pc$Purity)/.1];d<-d[is.finite(purity_10pp)];fo<-as.formula(paste(y,if(y=='switch')'~ns(abs(r0-50),3)+purity_10pp'else'~ns(m_r,3)+purity_10pp'));f<-if(y=='switch')glm(fo,d,family=binomial())else lm(fo,d);V<-vcovHC(f,type='HC3');b<-coef(f)['purity_10pp'];se<-sqrt(V['purity_10pp','purity_10pp']);ci<-b+c(-1,1)*qt(.975,df.residual(f))*se;p<-2*pt(-abs(b/se),df.residual(f));if(y=='switch'){b<-exp(b);ci<-exp(ci)}
 pur[[length(pur)+1]]<-data.table(label='exploratory',cohort='CheckMate025',signature=t,endpoint=y,n=nrow(d),missing=250-nrow(d),estimate=b,lower=ci[1],upper=ci[2],p=p,unit=if(y=='switch')'odds ratio per 10pp DNA purity'else'rank percentile points per 10pp DNA purity')
};pu<-rbindlist(pur);pu[,q_BH_secondary:=p.adjust(p,'BH')];wr(pu,'DNA_PURITY_ASSOCIATIONS.csv');wr(data.table(cohort=c('JAVELIN','IMmotion151','CheckMate025'),n_RNA=c(726,823,250),n_DNA_purity=c(0,0,sum(is.finite(as.numeric(pc$Purity)))),reason=c('No qualified independent DNA-derived purity field; H&E is a different measurement','No qualified independent DNA-derived purity field; xCell is same-RNA derived','S1 ABSOLUTE DNA purity source field')),'DNA_PURITY_AVAILABILITY.csv')
# Fixed source clusters, no re-fitting or new subtype claims.
jcl<-rd('source_JAVELIN_cluster.csv');clusters<-list(JAVELIN=data.table(ID=jcl$ID,cluster=jcl$TCGA_cluster),IMmotion151=data.table(ID=im$PATIENT_ID,cluster=as.character(im$NMF_GROUP)))
cldesc<-cltests<-clcoef<-list()
cramer<-function(tb){chi<-suppressWarnings(chisq.test(tb,correct=FALSE)$statistic);sqrt(as.numeric(chi)/(sum(tb)*min(nrow(tb)-1,ncol(tb)-1)))}
for(co in names(clusters))for(t in c('Angio36','Angio16')){
 d<-merge(pa[cohort==co&signature==t],clusters[[co]],by='ID',all.x=TRUE);d[is.na(cluster),cluster:='Unknown'];d[,cluster:=factor(cluster)];tb<-table(d$cluster,d$switch);v<-cramer(tb);bs<-replicate(2000,{ii<-sample.int(nrow(d),replace=TRUE);cramer(table(d$cluster[ii],factor(d$switch[ii],levels=c(FALSE,TRUE))))});ci<-quantile(bs,c(.025,.975),na.rm=TRUE);fip<-fisher.test(tb,simulate.p.value=TRUE,B=20000)$p.value
 cltests[[length(cltests)+1]]<-data.table(label='exploratory',cohort=co,signature=t,n=nrow(d),estimand='fixed-cluster composition: reclassified vs stable',estimate=v,lower=ci[1],upper=ci[2],p=fip,unit="Cramer's V",method='patient bootstrap CI 2000; Fisher Monte Carlo 20000')
 a<-fitcat(d,'delta_r','cluster');a[,`:=`(label='exploratory',cohort=co,signature=t)];clcoef[[length(clcoef)+1]]<-a
 for(k in levels(d$cluster)){
  z<-d[cluster==k];n<-nrow(z);sw<-sum(z$switch);bci<-binom.test(sw,n)$conf.int
  for(y in c('delta_r','abs_delta')){x<-z[[y]];eci<-mean(x)+c(-1,1)*qt(.975,n-1)*sd(x)/sqrt(n);cldesc[[length(cldesc)+1]]<-data.table(label='exploratory',cohort=co,signature=t,source_cluster=k,n=n,endpoint=y,estimate=mean(x),lower=eci[1],upper=eci[2],median=median(x),reclassified=sw,stable=n-sw,reclassification_rate=sw/n,reclassification_lower=bci[1],reclassification_upper=bci[2],fraction_of_all_reclassified=sw/sum(d$switch),fraction_of_all_stable=(n-sw)/sum(!d$switch))}
 }
}
wr(rbindlist(cldesc),'FIXED_CLUSTER_DISPLACEMENT_COMPOSITION.csv');wr(rbindlist(cltests),'FIXED_CLUSTER_COMPOSITION_TESTS.csv');wr(rbindlist(clcoef),'FIXED_CLUSTER_ADJUSTED_SIGNED_CONTRASTS.csv')
print(primary[,.(cohort,gene,n,n_reference,n_category,estimate,lower,upper,p,q_BH)]);print(rbindlist(audit)[gene=='VHL',.(cohort,RNA_n,n_qualified,n_unknown_not_assayed)]);print(pu)
