# Author / creator / last modified by: Chunhao Guo
.libPaths(c('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B1/Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(survival)})
P<-normalizePath('ccRCC_signature_transportability_benchmark_2026');W<-file.path(P,'84_JCO_PO_Upgrade_20261113');O<-file.path(W,'05_analysis/B1_operator_panel_trials');S<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B1';source(file.path(W,'05_analysis/shared/operators.R'));setDTthreads(1);set.seed(20261113)
tr<-commandArgs(TRUE)[1];old<-readRDS(file.path(P,'37_A_clinical_workflow_validation_20260917/local_only/INPUTS.rds'))[[tr]];d<-old$d
prepared_path<-file.path(S,paste0(tr,'_prepared.rds'))
if(file.exists(prepared_path)){pp<-readRDS(prepared_path);x<-pp$x;d<-pp$d;sets<-pp$sets;scale_desc<-pp$input_scale;rm(pp);gc()}else{
if(tr=='JAVELIN'){path<-file.path(P,'27_A_incremental_background_revision_v1_0/local_only/JAVELIN_EXPRESSION.csv.gz');e<-fread(path);x<-as.matrix(e[,d$ID,with=FALSE]);rownames(x)<-e$symbol;scale_desc<-'released log2 TPM (floored), no second transformation';keys<-d$ID}else{path<-file.path(P,'21_shared_datasets/IMmotion151_20201106/raw/IMmotion151.expression.data.TPM.anon.20201106(1).csv');e<-fread(path);x<-as.matrix(e[,d$RNA_ID,with=FALSE]);rn<-e$symbol;bad<-duplicated(rn)|duplicated(rn,fromLast=TRUE);if(any(bad)){dup<-data.frame(symbol=rn[bad],gene_id=e[[1]][bad],qualified_id=paste(rn[bad],e[[1]][bad],sep='|'));write.csv(dup,file.path(O,'results/IMmotion151_duplicate_symbol_resolution.csv'),row.names=FALSE);rn[bad]<-dup$qualified_id};rownames(x)<-rn;scale_desc<-'released log2(TPM+1), no second transformation';keys<-d$RNA_ID}
if(tr=='JAVELIN'){rn<-rownames(x);bad<-duplicated(rn)|duplicated(rn,fromLast=TRUE);if(any(bad)){dup<-data.frame(symbol=rn[bad],source_row=which(bad)+2,qualified_id=paste0(rn[bad],'|S13row',which(bad)+2));write.csv(dup,file.path(O,'results/JAVELIN_duplicate_symbol_resolution.csv'),row.names=FALSE);rn[bad]<-dup$qualified_id;rownames(x)<-rn}}
rm(e);gc();stopifnot(!anyDuplicated(rownames(x)),all(is.finite(x)),ncol(x)==nrow(d))
# Preserve source symbols; only frozen one-to-one aliases from existing mapping accepted.
mem<-read.csv(file.path(P,'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv'));sets<-list(Angio16=mem$gene[mem$Angio16],Angio36=mem$gene[mem$Angio36]);immune<-read.csv(file.path(O,'results/immune_gene_sources.csv'));sets<-c(sets,split(immune$gene,immune$signature));map<-fread(file.path(P,'24_A_GenomeMedicine_Revision_20260915/source_data/GENE_RESOLUTION.csv'));mappings<-list()
for(nm in names(sets)){g<-sets[[nm]];resolved<-vapply(g,function(z){if(sum(rownames(x)==z)==1)return(z);a<-unique(na.omit(map$resolved_javelin_symbol[map$printed_source_symbol==z]));a<-a[a%in%rownames(x)];if(length(a)==1)return(a);z},character(1));mappings[[nm]]<-data.frame(trial=tr,signature=nm,source_gene=g,matrix_gene=resolved,present=resolved%in%rownames(x));sets[[nm]]<-resolved}
write.csv(do.call(rbind,mappings),file.path(O,'results',paste0(tr,'_gene_coverage.csv')),row.names=FALSE)
# Verify the complete released matrix exactly reconstructs each old input subset.
for(sig in c('Angio16','Angio36'))stopifnot(max(abs(unname(x[sets[[sig]],,drop=FALSE])-unname(old$M[[sig]][mappings[[sig]]$source_gene,,drop=FALSE])))<1e-8)
write.csv(data.frame(trial=tr,genes=nrow(x),samples=ncol(x),unique_genes=length(unique(rownames(x))),whole_transcriptome=TRUE,input_scale=scale_desc,source_path=path,nonfinite=sum(!is.finite(x)),constant_rows=sum(apply(x,1,function(z)diff(range(z))==0))),file.path(O,'results',paste0(tr,'_matrix_qualification.csv')),row.names=FALSE)
# RDS intermediate >50 MB stays on external scratch.
saveRDS(list(x=x,d=d,sets=sets,input_scale=scale_desc),file.path(S,paste0(tr,'_prepared.rds')),compress=FALSE)
}
gm_canonical<-fread(file.path(O,'results',paste0(tr,'_gene_coverage.csv')))
for(sig in c('Angio16','Angio36')){source_order<-rownames(old$M[[sig]]);mm<-gm_canonical[signature==sig];stopifnot(setequal(source_order,mm$source_gene));sets[[sig]]<-mm$matrix_gene[match(source_order,mm$source_gene)]}
out<-list();avail<-list()
for(sig in names(sets)){ops<-if(grepl('^Angio',sig))paste0('O',1:7)else c('O1','O2','O5');out[[sig]]<-list();for(op in ops){cat(format(Sys.time()),tr,sig,op,'start\n');flush.console();t0<-proc.time()[3];b<-get(op)(x,sets[[sig]],whole_transcriptome=TRUE);out[[sig]][[op]]<-b;avail[[length(avail)+1]]<-data.frame(trial=tr,signature=sig,operator=op,n=ncol(x),available_n=sum(b$available),reason=b$reason,seconds=proc.time()[3]-t0);cat(tr,sig,op,'available',sum(b$available),b$reason,'seconds',proc.time()[3]-t0,'\n');flush.console();gc();saveRDS(out,file.path(S,paste0(tr,'_scores_progress.rds')))} }
saveRDS(list(scores=out,d=d,sets=sets),file.path(S,paste0(tr,'_scores.rds')))
write.csv(do.call(rbind,avail),file.path(O,'results',paste0(tr,'_operator_availability.csv')),row.names=FALSE)
# Export compact participant results locally (restricted trial inputs are not redistributed).
rr<-list();for(sig in names(out))for(op in names(out[[sig]])){b<-out[[sig]][[op]];rr[[length(rr)+1]]<-data.frame(trial=tr,ID=d$ID,signature=sig,operator=op,score=b$scores,available=b$available,research_high=if(all(b$available))b$scores>=median(b$scores)else NA)};fwrite(rbindlist(rr),file.path(O,'results',paste0(tr,'_scores_local_only.csv')))
writeLines(capture.output(sessionInfo()),file.path(O,paste0('sessionInfo_',tr,'.txt')))
