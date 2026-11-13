# Author / creator / last modified by: Chunhao Guo
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(1);set.seed(20261113)
ROOT<-'/LOCAL_PROJECT_ROOT'
WORK<-file.path(ROOT,'84_JCO_PO_Upgrade_20261113');OUT<-file.path(WORK,'05_analysis/B2_independent_cohorts')
SCR<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2'
ext<-file.path(ROOT,'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized')
member<-fread(file.path(ROOT,'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv'))
sets<-list(Angio16=member[Angio16==TRUE,gene],Angio36=member[Angio36==TRUE,gene])
oldsets<-fromJSON(file.path(ext,'FEATURE_SETS.json'));sets$APM12<-oldsets$APM12
write_json(sets,file.path(OUT,'results/gene_sets.json'),pretty=TRUE,auto_unbox=FALSE)
cl<-fread(file.path(ext,'CM025_CLINICAL.csv'));stopifnot(nrow(cl)==250,!anyDuplicated(cl$SUBJID),!anyDuplicated(cl$RNA_ID))
foldpath<-file.path(WORK,'05_analysis/shared/CM025_folds.csv')
if(!file.exists(foldpath)){
 folds<-cl[,.(patient=SUBJID,RNA_ID,treatment=Arm)];folds[,fold:=NA_integer_]
 for(a in sort(unique(folds$treatment))){ii<-which(folds$treatment==a);folds$fold[ii]<-sample(rep(1:5,length.out=length(ii)))}
 folds[,seed:=20261113L];fwrite(folds,foldpath)
}else folds<-fread(foldpath)
stopifnot(nrow(folds)==250,setequal(folds$patient,cl$SUBJID),all(folds$seed==20261113))
fwrite(folds[,.N,by=.(treatment,fold)],file.path(OUT,'results/CM025_fold_counts.csv'))
cat('FOLDS SAVED BEFORE ANY OUTCOME MODEL:',foldpath,'\n');flush.console()
fwrite(cl[,.(sample=RNA_ID,patient=SUBJID,include=TRUE,arm=Arm,specimen=Tumor_Sample_Primary_or_Metastasis,
 days_collection_to_treatment=Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy,SampleType)],file.path(OUT,'results/CM025_sample_manifest.csv'))
sources<-list();note<-function(cohort,path,scale){sources[[length(sources)+1]]<<-data.table(cohort,path=normalizePath(path),scale,bytes=file.info(path)$size,sha256=digest::digest(file=path,algo='sha256'))}
savecohort<-function(n,x,patient,scale,path){stopifnot(is.matrix(x),is.numeric(x),all(is.finite(x)),!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)),length(patient)==ncol(x));saveRDS(list(x=x,patient=patient,scale=scale,source=path),file.path(SCR,paste0(n,'.rds')),compress=FALSE);cat(n,dim(x),'saved\n');flush.console()}
# Source-released log2 upper-quartile-normalized TPM, ComBat corrected and +21;
# do NOT log again. Stream XLSX XML with 00_extract_cm025.py before running this.
wb<-file.path(ext,'source/Braun2020_supplement.xlsx');note('CM025',wb,'released log2 UQ-normalized TPM; ComBat; +21; no additional transform')
dd<-fread(file.path(SCR,'CM025_complete_expression.csv'));stopifnot(nrow(dd)==43893,ncol(dd)==251)
gene<-dd[[1]];x<-as.matrix(dd[,-1]);rownames(x)<-gene;x<-x[,cl$RNA_ID,drop=FALSE]
rownames(x)<-make.unique(gene)
alias<-c(CXCR7='ACKR3',FIGF='VEGFD');for(a in names(alias))if(!alias[[a]]%in%rownames(x)&&a%in%rownames(x))rownames(x)[rownames(x)==a]<-alias[[a]]
fwrite(data.table(original=gene,analysis=rownames(x)),file.path(OUT,'results/CM025_gene_mapping.csv'))
savecohort('CM025',x,cl$SUBJID,'released log2 UQ-normalized TPM; ComBat; +21; no additional transform',wb)
old<-fread(file.path(ext,'TARGET_GENE_EXPRESSION.csv'));ox<-as.matrix(old[,-1]);rownames(ox)<-old[[1]]
eq<-intersect(rownames(ox),rownames(x));fwrite(data.table(genes=length(eq),samples=ncol(x),max_abs_difference=max(abs(x[eq,]-ox[eq,colnames(x)]))),file.path(OUT,'results/CM025_source_expression_check.csv'))
rm(dd,x,ox);gc()
tp<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/Documents/生信/ccRCC_independent_bioinformatics_2026/05_data/reference/TRACERx_tx100_rna_2024/data/processed/tumour_tpm.rds'
note('TRACERx',tp,'log2(TPM+1) applied here');x<-readRDS(tp);pid<-sub('[-_].*$','',colnames(x));region<-sub('^[^_-]+[-_]','',colnames(x))
keep<-colnames(x)!='K207-R4'
stopifnot(all(grepl('^K[0-9]{3}$',pid)),all(grepl('^(R|LN|THR|M|T|CB)[0-9]*$',region)),!anyDuplicated(paste(pid[keep],region[keep])))
fwrite(data.table(sample=colnames(x),patient=pid,region,include=keep,source_tumour_matrix=TRUE,reason=ifelse(keep,'source tumour matrix; unique normalized region key','source-author excluded duplicate normalized K207_R4 key')),file.path(OUT,'results/TRACERx_sample_manifest.csv'))
x<-x[,keep,drop=FALSE];pid<-pid[keep]
fwrite(data.table(patient=pid)[,.N,by=patient],file.path(OUT,'results/TRACERx_patient_region_counts.csv'))
for(a in names(alias))if(!alias[[a]]%in%rownames(x)&&a%in%rownames(x))rownames(x)[rownames(x)==a]<-alias[[a]]
savecohort('TRACERx',log2(x+1),pid,'log2(TPM+1)',tp);rm(x);gc()
tc<-'/LOCAL_DATA_ROOT/Yang_RCC_metabolic_immune_plasticity_2026/02_data_inventory/raw/tcga_kirc/TCGA.KIRC.HiSeqV2.gz'
note('TCGA',tc,'released log2(normalized RSEM+1); no additional transform');dd<-fread(tc);gene<-dd[[1]];x<-as.matrix(dd[,-1]);rownames(x)<-make.unique(gene)
for(a in names(alias))if(!alias[[a]]%in%rownames(x)&&a%in%rownames(x))rownames(x)[rownames(x)==a]<-alias[[a]]
samp<-colnames(x);pid<-substr(samp,1,12);type<-substr(samp,14,15);keep<-type=='01'&!duplicated(ifelse(type=='01',pid,samp))
fwrite(data.table(sample=samp,patient=pid,sample_type=type,include=keep,reason=ifelse(keep,'one primary tumour per patient','nonprimary or repeated patient')),file.path(OUT,'results/TCGA_sample_manifest.csv'))
fwrite(data.table(original=gene,analysis=rownames(x)),file.path(OUT,'results/TCGA_gene_mapping.csv'))
savecohort('TCGA',x[,keep,drop=FALSE],pid[keep],'released log2(normalized RSEM+1)',tc);rm(dd,x);gc()
meta<-fromJSON(file.path(ROOT,'43_evidence_linked_exploration_20260918/results/CPTAC_GDC_FILE_METADATA.json'),simplifyVector=FALSE)$data$hits
files<-list.files('/LOCAL_EXTERNAL_DRIVE/生信_public_data/cptac_ccrcc_gdc_rna_star_counts_v21',pattern='tsv$',full.names=TRUE);names(files)<-sub('__.*','',basename(files))
mp<-unique(fread(file.path(ROOT,'43_evidence_linked_exploration_20260918/protein_inputs/TARGET_GENE_DEFINITIONS.csv'))[,.(ensembl_gene_id,approved_symbol)])
mp[approved_symbol=='COXFA4L2',approved_symbol:='NDUFA4L2'] # same ENSG00000185633; locked membership name
xx<-list();manifest<-list();gm<-NULL
for(i in seq_along(meta)){
 h<-meta[[i]];stopifnot(length(h$cases)==1);cas<-h$cases[[1]];ss<-cas$samples
 stopifnot(all(vapply(ss,function(s)s$sample_type=='Primary Tumor',logical(1))))
 f<-files[[h$file_id]];dd<-fread(f,skip=1,select=c('gene_id','gene_name','tpm_unstranded'));dd<-dd[grepl('^ENSG',gene_id)]
 ids<-sub('\\..*$','',dd$gene_id);gn<-dd$gene_id
 stopifnot(!anyDuplicated(gn),!anyDuplicated(ids[ids%in%mp$ensembl_gene_id]));j<-match(mp$ensembl_gene_id,ids);stopifnot(!anyNA(j));gn[j]<-mp$approved_symbol
 if(is.null(gm)){gm<-data.table(gene_id=dd$gene_id,source_symbol=dd$gene_name,analysis=gn);xx<-matrix(NA_real_,nrow=nrow(dd),ncol=length(meta),dimnames=list(gn,rep('',length(meta))))}else stopifnot(identical(gn,rownames(xx)))
 xx[,i]<-dd$tpm_unstranded;colnames(xx)[i]<-cas$submitter_id
 manifest[[i]]<-data.table(sample=cas$submitter_id,patient=cas$submitter_id,file=f,file_id=h$file_id,RNA_aliquot=h$associated_entities[[1]]$entity_submitter_id,
 source_sample_labels=paste(vapply(ss,function(s)s$submitter_id,character(1)),collapse=';'),n_source_labels=length(ss),include=TRUE,sha256=digest::digest(file=f,algo='sha256'))
}
fwrite(gm,file.path(OUT,'results/CPTAC_gene_mapping.csv'));fwrite(rbindlist(manifest),file.path(OUT,'results/CPTAC_sample_manifest.csv'))
old<-fread(file.path(ROOT,'43_evidence_linked_exploration_20260918/protein_inputs/RNA_TARGET_TPM.csv'));ox<-as.matrix(old[,-1]);rownames(ox)<-old[[1]]
rownames(ox)[rownames(ox)=='COXFA4L2']<-'NDUFA4L2'
fwrite(data.table(genes=nrow(ox),cases=ncol(ox),max_abs_difference=max(abs(xx[rownames(ox),colnames(ox)]-ox))),file.path(OUT,'results/CPTAC_source_expression_check.csv'))
savecohort('CPTAC',log2(xx+1),colnames(xx),'log2(GDC STAR tpm_unstranded+1)',dirname(files[1]));rm(xx);gc()
fwrite(rbindlist(sources),file.path(OUT,'results/input_provenance.csv'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/preparation_sessionInfo.txt'))
