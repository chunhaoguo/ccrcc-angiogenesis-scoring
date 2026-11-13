# Author, creator and last modified by: Chunhao Guo
# Research scorer; not a clinical decision tool. No fitting on incoming cases.
read_reference <- function(path, manifest_path) {
  if (!requireNamespace('jsonlite',quietly=TRUE) || !requireNamespace('digest',quietly=TRUE))
    stop('DEPENDENCY_MISSING: jsonlite and digest required')
  manifest <- jsonlite::read_json(manifest_path,simplifyVector=FALSE)
  hit <- Filter(function(z) basename(z$file)==basename(path),manifest$parameters)
  if(length(hit)!=1L)stop('REFERENCE_NOT_IN_MANIFEST')
  check_reference_hash(path,hit[[1]]$sha256)
  r <- jsonlite::read_json(path,simplifyVector=TRUE)
  validate_reference(r)
  r
}
check_reference_hash <- function(path,expected_sha256) {
  if(digest::digest(file=path,algo='sha256')!=expected_sha256)stop('REFERENCE_HASH_MISMATCH')
  invisible(TRUE)
}
validate_reference <- function(r) {
  n<-length(r$genes)
  if(!r$method %in% c('I0','I1') || !r$target %in% c('Angio16','Angio36'))stop('UNSUPPORTED_IMPLEMENTATION')
  expected<-if(r$target=='Angio16')16L else 36L
  if(n!=expected || anyDuplicated(r$genes) || any(!nzchar(r$genes)))stop('INVALID_GENE_LIST')
  for(k in c('gene_mean','gene_sd','weights'))if(length(r[[k]])!=n || any(!is.finite(r[[k]])))stop('INVALID_REFERENCE_PARAMETERS')
  if(length(r$unit)!=1L || !is.finite(r$unit) || r$unit<=1e-8)stop('UNAVAILABLE_FINAL_SCALE')
  if(any(!is.finite(c(r$center,r$intercept))))stop('INVALID_REFERENCE_PARAMETERS')
  if(any(r$gene_sd<=1e-8))stop('UNAVAILABLE_GENE_REFERENCE_SD')
  w<-1/(n*r$gene_sd*r$unit)
  b<-(-mean(r$gene_mean/r$gene_sd)-r$center)/r$unit
  if(max(abs(w-r$weights))>1e-10 || abs(b-r$intercept)>1e-10)stop('INCONSISTENT_REFERENCE_FORMULA')
  if(r$method=='I0' && (any(r$gene_mean!=0)||any(r$gene_sd!=1)))stop('INVALID_MEAN_EXPRESSION_REFERENCE')
  if(!is.null(r$research_threshold) && (length(r$research_threshold)!=1L || !is.finite(r$research_threshold)))stop('INVALID_THRESHOLD')
  invisible(TRUE)
}
classify_score <- function(score, threshold=NULL) {
  if(is.null(threshold))return(ifelse(is.finite(score),'unavailable: threshold not frozen','unavailable'))
  if(length(threshold)!=1L || !is.finite(threshold))stop('INVALID_THRESHOLD')
  ifelse(!is.finite(score),'unavailable',ifelse(score>=threshold,'high','low'))
}
score_batch <- function(expression, reference, declaration) {
  r<-reference
  validate_reference(r)
  for(k in c('platform','input_scale','gene_identifiers'))
    if(is.null(declaration[[k]]) || !identical(declaration[[k]],r[[k]]))stop(paste0('DECLARATION_MISMATCH:',k))
  # A declaration is necessary but does not certify numeric scale or assay processing.
  x<-as.data.frame(expression,check.names=FALSE,stringsAsFactors=FALSE)
  if(ncol(x)<2L || nrow(x)<1L)stop('INVALID_MATRIX')
  genes<-x[[1]];samples<-names(x)[-1]
  if(!is.character(genes)||anyNA(genes)||any(!nzchar(genes))||anyDuplicated(genes))stop('INVALID_OR_DUPLICATE_GENE')
  if(anyNA(samples)||any(!nzchar(samples))||anyDuplicated(samples))stop('INVALID_SAMPLE_NAMES')
  if(!all(vapply(x[-1],is.numeric,logical(1))))stop('NONNUMERIC_INPUT')
  missing<-setdiff(r$genes,genes)
  if(length(missing))return(data.frame(sample=samples,score=NA_real_,research_group='unavailable',status='unavailable',reason=paste0('MISSING_REQUIRED_GENE:',paste(missing,collapse=';'))))
  m<-as.matrix(x[match(r$genes,genes),-1,drop=FALSE]);ok<-apply(m,2,function(v)all(is.finite(v)))
  z<-rep(NA_real_,length(samples));z[ok]<-as.numeric(crossprod(r$weights,m[,ok,drop=FALSE]))+r$intercept
  ok<-ok & is.finite(z)
  z[!ok]<-NA_real_
  data.frame(sample=samples,score=z,research_group=classify_score(z,r$research_threshold),
    status=ifelse(ok,'available','unavailable'),reason=ifelse(ok,'','NONFINITE_REQUIRED_INPUT_OR_OUTPUT'))
}
select_reference <- function(package_dir,target,method,platform,design='primary') {
  if(!target %in% c('Angio16','Angio36') || !method %in% c('I0','I1') || !platform %in% c('RNAseq','Agilent') || !design %in% c('primary','strict11'))stop('UNSUPPORTED_IMPLEMENTATION')
  if(target=='Angio36')stop(if(platform=='Agilent')'UNAVAILABLE_FROZEN_PARAMETERS: complete Angio36 array lacks FGF16' else 'UNAVAILABLE_FROZEN_PARAMETERS: Angio36 RNAseq parameters not present in declared source assets')
  p<-file.path(package_dir,'parameters',paste0(paste(design,target,method,platform,sep='__'),'.json'))
  read_reference(p,file.path(package_dir,'parameter_manifest.json'))
}
score_single <- function(gene_values,reference,declaration,sample_name='SYNTH_SINGLE') {
  if(!is.numeric(gene_values)||is.null(names(gene_values)))stop('NAMED_NUMERIC_VECTOR_REQUIRED')
  d<-data.frame(gene=names(gene_values),value=as.numeric(gene_values),check.names=FALSE)
  names(d)[2]<-sample_name;score_batch(d,reference,declaration)
}
