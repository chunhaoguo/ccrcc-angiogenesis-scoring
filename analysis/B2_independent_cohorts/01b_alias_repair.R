# Author / creator / last modified by: Chunhao Guo
# One-time repair of a matrix label to the immutable locked membership name.
# The full prepare script now performs this mapping directly.
suppressPackageStartupMessages(library(data.table));setDTthreads(1)
ROOT<-'/LOCAL_PROJECT_ROOT';OUT<-file.path(ROOT,'84_JCO_PO_Upgrade_20261113/05_analysis/B2_independent_cohorts');SCR<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2'
mp<-fread(file.path(ROOT,'43_evidence_linked_exploration_20260918/protein_inputs/TARGET_GENE_DEFINITIONS.csv'))
stopifnot(mp[printed_gene=='NDUFA4L2',unique(ensembl_gene_id)]=='ENSG00000185633')
z<-readRDS(file.path(SCR,'CPTAC.rds'));rn<-rownames(z$x);stopifnot('COXFA4L2'%in%rn,!'NDUFA4L2'%in%rn);rn[rn=='COXFA4L2']<-'NDUFA4L2';rownames(z$x)<-rn;saveRDS(z,file.path(SCR,'CPTAC.rds'),compress=FALSE)
gm<-fread(file.path(OUT,'results/CPTAC_gene_mapping.csv'));gm[analysis=='COXFA4L2',analysis:='NDUFA4L2'];fwrite(gm,file.path(OUT,'results/CPTAC_gene_mapping.csv'))
old<-fread(file.path(ROOT,'43_evidence_linked_exploration_20260918/protein_inputs/RNA_TARGET_TPM.csv'));ox<-as.matrix(old[,-1]);rownames(ox)<-old[[1]];rownames(ox)[rownames(ox)=='COXFA4L2']<-'NDUFA4L2'
err<-max(abs((2^z$x[rownames(ox),colnames(ox)]-1)-ox));stopifnot(err<1e-8)
fwrite(data.table(locked_symbol='NDUFA4L2',alternate_symbol='COXFA4L2',ensembl_gene_id='ENSG00000185633',max_TPM_roundtrip_difference=err),file.path(OUT,'results/CPTAC_alias_check.csv'))
