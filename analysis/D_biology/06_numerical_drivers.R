# Author: Chunhao Guo
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(1);set.seed(20261113);P<-normalizePath(commandArgs(TRUE)[1]);O<-file.path(P,'84_JCO_PO_Upgrade_20261113/05_analysis/D_biology');R<-file.path(O,'results');wr<-function(x,n)fwrite(x,file.path(R,n),na='NA');ob<-readRDS(file.path(R,'qualified_inputs.rds'));mat<-ob$mat;ids<-ob$ids
A<-file.path(P,'60_TME_Definition_Operator_Exploration_20260928/results');oldw<-fread(file.path(A,'ALL_GENE_WEIGHTS.csv'));oldv<-fread(file.path(A,'ALL_GENE_VARIANCE_CONTRIBUTIONS.csv'));oldg<-fread(file.path(A,'GROUP_VARIANCE_CONTRIBUTIONS.csv'));wr(oldw,'CARRIED_FORWARD_ALL_GENE_WEIGHTS.csv');wr(oldv,'CARRIED_FORWARD_ALL_GENE_VARIANCE_CONTRIBUTIONS.csv');wr(oldg,'CARRIED_FORWARD_GROUP_VARIANCE_CONTRIBUTIONS.csv')
cpt<-fread(file.path(P,'60_TME_Definition_Operator_Exploration_20260928/local_only/AUDIT_INPUT_CPTAC.csv'));mat$CPTAC<-lapply(names(ob$genes),function(t){z<-as.matrix(cpt[,paste0(t,'__',ob$genes[[t]]),with=FALSE]);z<-t(z);rownames(z)<-ob$genes[[t]];z});names(mat$CPTAC)<-names(ob$genes);ids$CPTAC<-cpt$ID
trp<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/Documents/生信/ccRCC_independent_bioinformatics_2026/05_data/reference/TRACERx_tx100_rna_2024/data/processed/tumour_tpm.rds';tr<-readRDS(trp);stopifnot(all(is.finite(tr)),all(tr>=0));mat$TRACERx_regions<-lapply(ob$genes,function(g){alias<-ifelse(g=='VEGFD','FIGF',g);z<-log2(1+tr[alias,,drop=FALSE]);rownames(z)<-g;z});ids$TRACERx_regions<-colnames(tr);rm(tr);gc(FALSE)
tc<-fread(file.path(R,'TCGA_PRIMARY_TARGET_EXPRESSION.csv'));tcx<-as.matrix(tc[,-1]);rownames(tcx)<-tc[[1]];mat$TCGA_KIRC_primary<-lapply(ob$genes,function(g)tcx[g,,drop=FALSE]);ids$TCGA_KIRC_primary<-colnames(tcx)
score<-function(X,tr=seq_len(ncol(X))){mu<-rowMeans(X[,tr,drop=FALSE]);ss<-apply(X[,tr,drop=FALSE],1,sd);r0<-colMeans(X);r1<-if(all(is.finite(ss)&ss>1e-8))colMeans((X-mu)/ss)else rep(NA_real_,ncol(X));w0<-rep(1/(nrow(X)*sd(r0[tr])),nrow(X));w1<-1/(nrow(X)*ss*sd(r1[tr]));list(I0=(r0-mean(r0[tr]))/sd(r0[tr]),I1=(r1-mean(r1[tr]))/sd(r1[tr]),w0=w0,w1=w1,mu=mu,sd=ss)}
fl<-qc<-atts<-sums<-patientdec<-relations<-availability<-list()
for(co in names(mat))for(t in names(ob$genes)){
 X<-mat[[co]][[t]];s<-score(X);g<-rownames(X);arch<-oldw[cohort==co&signature==t];valid<-all(is.finite(s$I1));err<-NA_real_
 if(nrow(arch)){
  arch<-arch[match(g,gene)];err<-max(abs(arch$weight_I0-s$w0),abs(arch$weight_I1-s$w1));stopifnot(err<1e-8)
  # Reuse saved weight quantities; only floor statistics and switcher attribution are new.
  s$w0<-arch$weight_I0;s$w1<-arch$weight_I1;group<-arch$group
 }else{group<-ifelse(g%in%intersect(ob$genes$Angio16,ob$genes$Angio36),'Shared 6',ifelse(t=='Angio16','Angio16-only',ifelse(grepl('^FGF',g),'Angio36-only FGF/FGFR','Angio36-only other')))}
 minimum<-apply(X,1,min);floor_fraction<-rowMeans(X==minimum);mode<-apply(X,1,function(z){v<-unique(z);v[which.max(tabulate(match(z,v)))]});zero<-rowMeans(X==0);ratio<-s$w1/s$w0
 fl[[length(fl)+1]]<-data.table(label='exploratory',interpretation='arithmetic; expression floor is not a validated assay LOD',cohort=co,signature=t,gene=g,group=group,n=ncol(X),reference_sd=s$sd,released_minimum=minimum,fraction_at_released_minimum=floor_fraction,fraction_exact_zero=zero,modal_value=mode,fraction_at_modal_value=rowMeans(X==mode),weight_I0=s$w0,weight_I1=s$w1,weight_I1_relative_I0=ratio,weight_source=if(nrow(arch))'existing module60 ALL_GENE_WEIGHTS.csv'else'new complete-list reference calculation')
 relations[[length(relations)+1]]<-data.table(label='exploratory',cohort=co,signature=t,n_genes=length(g),spearman_floor_weight_gain=if(valid)cor(floor_fraction,ratio,method='spearman')else NA_real_,median_floor_fraction_FGF=if(t=='Angio36')median(floor_fraction[grepl('^FGF',g)])else NA_real_,n_FGF_at_floor_ge90pct=sum(grepl('^FGF',g)&floor_fraction>=.9),note='finite gene-set arithmetic descriptor; no independent-gene inference or biological activation claim')
 availability[[length(availability)+1]]<-data.table(cohort=co,signature=t,n=ncol(X),I1_available=valid,failed_genes=paste(g[!is.finite(s$sd)|s$sd<=1e-8],collapse='|'),unit=if(co=='TRACERx_regions')'regions, not patients'else'patients/cases')
 if(!valid)next
 # Both full-reference and nominated frozen-reference calculations; complete gene list preserved.
 for(modeid in if(co%in%names(ob$mat))c('full_cohort','heldout_nominated')else'full_cohort'){
  n<-ncol(X);C<-matrix(NA_real_,nrow(X),n,dimnames=list(g,ids[[co]]));S0<-S1<-c0<-c1<-rep(NA_real_,n);floorct<-X==minimum
  if(modeid=='full_cohort'){
   C<-(X-s$mu)*(s$w1-s$w0);S0<-s$I0;S1<-s$I1;c0[]<-median(S0);c1[]<-median(S1)
  }else{
   d<-ob$patients[cohort==co&signature==t];fold<-d$fold
   if(all(is.na(fold))){C<-(X-s$mu)*(s$w1-s$w0);S0<-s$I0;S1<-s$I1;c0[]<-median(S0);c1[]<-median(S1)}else for(k in sort(unique(fold))){tr<-which(fold!=k);ee<-which(fold==k);z<-score(X,tr);C[,ee]<-(X[,ee,drop=FALSE]-z$mu)*(z$w1-z$w0);S0[ee]<-z$I0[ee];S1[ee]<-z$I1[ee];c0[ee]<-median(z$I0[tr]);c1[ee]<-median(z$I1[tr])}
  }
  errrec<-max(abs(colSums(C)-(S1-S0)));stopifnot(errrec<1e-8);sw<-(S0>=c0)!=(S1>=c1);direction<-ifelse(S0>=c0,'High to low','Low to high')
  if(co%in%names(ob$mat)){d<-ob$patients[cohort==co&signature==t];if(modeid=='heldout_nominated')stopifnot(identical(as.logical(sw),d$switch))else stopifnot(identical(as.logical(sw),d$switch_full))}
  qc[[length(qc)+1]]<-data.table(cohort=co,signature=t,reference=modeid,n=n,max_archived_weight_error=err,max_reconstruction_error=errrec,n_reclassified=sum(sw))
  for(dir in c('Low to high','High to low')){
   ii<-which(sw&direction==dir);if(!length(ii))next
   cc<-C[,ii,drop=FALSE];den<-sum(cc);absden<-sum(abs(cc));mc<-rowMeans(cc);se<-if(length(ii)>1)apply(cc,1,sd)/sqrt(length(ii))else rep(NA_real_,nrow(X));crit<-qt(.975,length(ii)-1)
   a<-data.table(label='exploratory',interpretation='exact arithmetic contribution to standardised score difference; not a biological or rank-displacement share',cohort=co,signature=t,reference=modeid,direction=dir,gene=g,group=group,n_reclassified=length(ii),mean_signed_contribution=mc,lower_conditional_mean=mc-crit*se,upper_conditional_mean=mc+crit*se,mean_absolute_contribution=rowMeans(abs(cc)),share_signed_total=rowSums(cc)/den,share_absolute_total=rowSums(abs(cc))/absden,fraction_at_full_reference_floor=rowMeans(floorct[,ii,drop=FALSE]))
   sums[[length(sums)+1]]<-a
   # All per-gene per-patient contributions retained, not only largest drivers.
   atts[[length(atts)+1]]<-data.table(label='exploratory_arithmetic',cohort=co,signature=t,reference=modeid,direction=dir,ID=rep(ids[[co]][ii],each=nrow(X)),gene=rep(g,length(ii)),group=rep(group,length(ii)),released_expression=as.vector(X[,ii,drop=FALSE]),at_full_reference_floor=as.vector(floorct[,ii,drop=FALSE]),contribution=as.vector(cc))
  }
  patientdec[[length(patientdec)+1]]<-data.table(cohort=co,signature=t,reference=modeid,ID=ids[[co]][sw],direction=direction[sw],I0=S0[sw],I1=S1[sw],sum_gene_contributions=colSums(C)[sw],threshold_I0=c0[sw],threshold_I1=c1[sw],threshold_contribution=(c0-c1)[sw],margin_difference=(S1-c1-S0+c0)[sw])
 }
}
wr(rbindlist(fl),'GENE_FLOOR_SD_WEIGHT.csv');wr(rbindlist(relations),'FLOOR_WEIGHT_RELATIONSHIP.csv');wr(rbindlist(availability),'DRIVER_SCORE_AVAILABILITY.csv');wr(rbindlist(qc),'ATTRIBUTION_REPRODUCTION_QA.csv');su<-rbindlist(sums);wr(su,'RECLASSIFIED_GENE_ATTRIBUTION_SUMMARY.csv');wr(rbindlist(atts),'RECLASSIFIED_PATIENT_GENE_ATTRIBUTIONS.csv.gz');wr(rbindlist(patientdec),'RECLASSIFIED_PATIENT_DECOMPOSITION.csv')
gr<-su[,.(n_reclassified=unique(n_reclassified),mean_signed_contribution=sum(mean_signed_contribution),mean_absolute_contribution=sum(mean_absolute_contribution),share_signed_total=sum(share_signed_total),share_absolute_total=sum(share_absolute_total)),by=.(label,interpretation,cohort,signature,reference,direction,group)];wr(gr,'RECLASSIFIED_GENE_GROUP_ATTRIBUTION.csv');top<-su[signature=='Angio36'][order(-share_absolute_total),head(.SD,5),by=.(cohort,reference,direction)];wr(top,'TOP_FIVE_ARITHMETIC_DRIVERS.csv');print(rbindlist(relations)[signature=='Angio36']);print(rbindlist(qc));writeLines(capture.output(sessionInfo()),file.path(O,'logs/06_driver_sessionInfo.txt'))
