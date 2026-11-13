# Author / creator / last modified by: Chunhao Guo
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(GSVA)})
setDTthreads(1);set.seed(20261113)
ROOT<-'/LOCAL_PROJECT_ROOT';WORK<-file.path(ROOT,'84_JCO_PO_Upgrade_20261113');OUT<-file.path(WORK,'05_analysis/B2_independent_cohorts')
source(file.path(WORK,'05_analysis/shared/operators.R'))
sources<-c(file.path(WORK,'05_analysis/ANALYSIS_PLAN_LOCKED_20261113.md'),file.path(WORK,'05_analysis/shared/operators.R'),file.path(WORK,'05_analysis/shared/CM025_folds.csv'),
 file.path(ROOT,'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv'),
 file.path(ROOT,'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized',c('CM025_CLINICAL.csv','FEATURE_SETS.json','CLINICAL_SOURCE_BRIDGE.json','SOURCE_PROCESSING_PASSAGES.txt')),
 file.path(ROOT,'43_evidence_linked_exploration_20260918',c('results/CPTAC_GDC_FILE_METADATA.json','protein_inputs/CASE_PAIRING.csv','protein_inputs/PROTEIN_MVD4.csv','protein_inputs/TARGET_GENE_DEFINITIONS.csv')),
 file.path(ROOT,'70_npj_Visual_Narrative_Final_20261011/figures/source_data/evidence/ALL_SCORE_MARKER_CORRELATIONS.csv'),
 '/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/Documents/生信/ccRCC_independent_bioinformatics_2026/05_data/reference/TRACERx_tx100_rna_2024/analysis/scripts/04_consensustme_vs_histo.R')
fwrite(data.table(path=sources,bytes=file.info(sources)$size,sha256=vapply(sources,function(p)digest::digest(file=p,algo='sha256'),character(1))),file.path(OUT,'results/contract_source_hashes.csv'))
x<-matrix(rnorm(100*24),100,24,dimnames=list(paste0('g',1:100),paste0('s',1:24)));x[1:16,]<-x[1:16,]+rep(seq(-2,2,length.out=24),each=16);g<-paste0('g',1:16)
std<-function(v)as.numeric(scale(v));m<-x[g,,drop=FALSE];z<-t(scale(t(m)));rob<-sweep(sweep(m,1,apply(m,1,median),'-'),1,1.4826*apply(m,1,mad,constant=1),'/')
pc<-prcomp(t(z),center=FALSE,scale.=FALSE,rank.=1);v<-as.numeric(t(z)%*%pc$rotation[,1]);v<-v*sign(cor(v,colMeans(m)))
direct<-list(O1=std(colMeans(m)),O2=std(colMeans(z)),O3=std(colMeans(rob)),
 O4=apply(x,2,function(v)(mean(rank(v,ties.method='min')[1:16])-8.5)/(100-16)-.5),
 O5=as.numeric(gsva(ssgseaParam(x,list(signature=g),minSize=16,maxSize=16,alpha=.25,normalize=FALSE),verbose=FALSE,BPPARAM=BiocParallel::SerialParam())),
 O6=as.numeric(gsva(gsvaParam(x,list(signature=g),minSize=16,maxSize=16,kcdf='Gaussian',tau=1,filterRows=FALSE),verbose=FALSE,BPPARAM=BiocParallel::SerialParam())),O7=std(v))
parity<-rbindlist(lapply(names(direct),function(o){got<-get(o)(x,g,whole_transcriptome=TRUE)$scores;data.table(operator=o,max_absolute_difference=max(abs(got-direct[[o]])),comparison='shared operators.R versus direct locked-definition calculation; same GSVA API for O5/O6')}))
stopifnot(all(parity$max_absolute_difference<1e-10));fwrite(parity,file.path(OUT,'results/operator_synthetic_parity.csv'))
# Exact centred singscore toy endpoints and tie handling.
toy<-matrix(c(1:5,5:1,rep(1,5)),5,3,dimnames=list(letters[1:5],c('bottom','top','tied')))
toy_result<-O4(toy,c('a','b'),whole_transcriptome=TRUE)$scores;stopifnot(max(abs(toy_result-c(-.5,.5,-2/3)))<1e-12)
fwrite(data.table(sample=names(toy_result),expected=c(-.5,.5,-2/3),observed=toy_result,ties_method='min; singscore rankGenes default'),file.path(OUT,'results/singscore_toy_check.csv'))
checks<-list();check<-function(name,value,detail){checks[[length(checks)+1]]<<-data.table(check=name,pass=isTRUE(value),detail)}
repro<-fread(file.path(OUT,'results/CM025_I0_reproduction.csv'));check('I0 manuscript interaction reproduction',all(repro$pass_rounded_manuscript),'both signatures, n250, 222 PFS events')
fold<-fread(file.path(WORK,'05_analysis/shared/CM025_folds.csv'));check('fold balance',all(fold[,.N,by=.(treatment,fold)]$N%in%c(24,26)),'120 nivolumab and 130 everolimus, five balanced folds')
for(co in c('CM025','TRACERx','TCGA','CPTAC')){
 p<-fread(file.path(OUT,'results',paste0(co,'_pairwise_metrics.csv')));st<-fread(file.path(OUT,'results',paste0(co,'_operator_status.csv')));cov<-fread(file.path(OUT,'results',paste0(co,'_signature_coverage.csv')))
 check(paste(co,'full 7x7 matrices'),nrow(p)==98,'49 entries per signature; unavailable entries retained')
 check(paste(co,'symmetry'),isTRUE(all.equal(p$disagreement_percent,p$disagreement_percent[order(p$signature,p$operator_b,p$operator_a)],check.attributes=FALSE)),'includes NAs')
 check(paste(co,'complete signatures'),all(cov$complete),'no omitted or zero-filled signature genes')
 check(paste(co,'all operator status rows'),nrow(st)==14,'unavailable operators explicitly reported')
 check(paste(co,'same shared source hash'),length(unique(st$implementation_sha256))==1,unique(st$implementation_sha256))
 check(paste(co,'final shared source hash'),all(st$implementation_sha256==digest::digest(file=file.path(WORK,'05_analysis/shared/operators.R'),algo='sha256')),'all cohorts use the final unchanged shared implementation')
 xx<-readRDS(file.path('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2',paste0(co,'.rds')))$x
 ss<-readRDS(file.path(OUT,'results',paste0(co,'_scores.rds')));sets<-jsonlite::fromJSON(file.path(OUT,'results/gene_sets.json'))
 for(s in c('Angio16','Angio36'))for(op in c('O1','O2')){
  expected<-get(op)(xx,sets[[s]])$scores;check(paste(co,s,op,'direct replay'),max(abs(expected-ss[[s]][,op]),na.rm=TRUE)<1e-10,'replay final full-cohort reference from immutable input')
 }
 rm(xx);gc()
}
tr<-fread(file.path(OUT,'results/TRACERx_sample_manifest.csv'));check('TRACERx corrected lineage',sum(tr$include)==230&&uniqueN(tr[include==TRUE]$patient)==79&&!anyDuplicated(tr[include==TRUE,.(patient,region)]),'231 source measurements, 230 distinct retained region keys, 79 patients')
sel<-fread(file.path(OUT,'results/TRACERx_one_region_selection.csv'));tx<-readRDS('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2/TRACERx.rds')$x[,sel$sample,drop=FALSE]
gs<-jsonlite::fromJSON(file.path(OUT,'results/gene_sets.json'));tqc<-rbindlist(lapply(c('Angio16','Angio36'),function(s){g<-gs[[s]];m<-tx[g,,drop=FALSE];data.table(signature=s,gene=g,reference_SD=apply(m,1,sd),reference_unscaled_MAD=apply(m,1,mad,constant=1),minimum=apply(m,1,min),fraction_at_minimum=apply(m,1,function(v)mean(v==min(v))))}));fwrite(tqc,file.path(OUT,'results/TRACERx_one_region_signature_reference_QC.csv'))
fr<-fread(file.path(OUT,'results/CM025_frozen_scores_transitions.csv'));check('CM025 frozen sample rows',nrow(fr)==500&&all(fr[,.N,by=signature]$N==250),'one held-out row per patient and signature')
pa<-fread(file.path(OUT,'results/CPTAC_protein_anchor.csv'));check('CPTAC protein anchor',all(pa$n==78)&&nrow(pa)==14,'four-protein proxy; all operators/signatures, including unavailable')
old_anchor<-fread(file.path(ROOT,'70_npj_Visual_Narrative_Final_20261011/figures/source_data/evidence/ALL_SCORE_MARKER_CORRELATIONS.csv'))
check('CPTAC existing protein anchor reproduction',abs(pa[signature=='Angio36'&operator=='O1',spearman]-old_anchor[cohort=='CPTAC'&variant=='primary'&marker=='Vascular_protein'&score=='A36_I0',rho])<1e-12,'same 78-case O1 association as existing source-qualified vascular proxy; folder70 ALL_SCORE_MARKER_CORRELATIONS.csv')
fwrite(rbindlist(checks),file.path(OUT,'results/validation_checks.csv'));print(parity);print(rbindlist(checks));stopifnot(all(rbindlist(checks)$pass))
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/validation_sessionInfo.txt'))
