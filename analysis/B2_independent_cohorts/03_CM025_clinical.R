# Author / creator / last modified by: Chunhao Guo
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(survival)})
setDTthreads(1);set.seed(20261113)
ROOT<-'/LOCAL_PROJECT_ROOT'
WORK<-file.path(ROOT,'84_JCO_PO_Upgrade_20261113');OUT<-file.path(WORK,'05_analysis/B2_independent_cohorts')
SCR<-'/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2'
source(file.path(WORK,'05_analysis/shared/operators.R'))
sets<-fromJSON(file.path(OUT,'results/gene_sets.json'));ext<-file.path(ROOT,'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized')
cl<-fread(file.path(ext,'CM025_CLINICAL.csv'));z<-readRDS(file.path(SCR,'CM025.rds'));X<-z$x[,cl$RNA_ID,drop=FALSE]
foldpath<-file.path(WORK,'05_analysis/shared/CM025_folds.csv');stopifnot(file.exists(foldpath));folds<-fread(foldpath);fold<-folds$fold[match(cl$SUBJID,folds$patient)]
stopifnot(nrow(cl)==250,!anyNA(fold),sum(cl$PFS_CNSR)==222,sum(cl$Arm=='NIVOLUMAB')==120,all(sets$APM12%in%rownames(X)))
apm<-O2(X,sets$APM12)$scores
d<-data.frame(patient=cl$SUBJID,trt=as.integer(cl$Arm=='NIVOLUMAB'),i=apm)
fitint<-function(v,endpoint){
 dd<-d;dd$v<-as.numeric(scale(v));dd$time<-cl[[endpoint]];dd$event<-cl[[paste0(endpoint,'_CNSR')]]
 if(any(!is.finite(dd$v)))return(data.table(n_total=250,n=0,events=NA,estimate=NA,lower=NA,upper=NA,p=NA,reason='operator unavailable'))
 ok<-complete.cases(dd)&is.finite(dd$time)&dd$time>0;dd<-dd[ok,];fit<-coxph(Surv(time,event)~v+i+trt+i:trt+v:trt,data=dd,ties='efron');s<-summary(fit)$coefficients['v:trt',];b<-s[1];se<-s[3]
 data.table(n_total=250,n=nrow(dd),events=sum(dd$event),estimate=exp(b),lower=exp(b-qnorm(.975)*se),upper=exp(b+qnorm(.975)*se),p=s[5],reason='available',PH_interaction_p=cox.zph(fit)$table['v:trt',3])
}
if(commandArgs(TRUE)[1]=='reproduce'){
 out<-rbindlist(lapply(c('Angio16','Angio36'),function(s){v<-fitint(O1(X,sets[[s]])$scores,'PFS');v[,signature:=s];v}))
 expected<-data.table(signature=c('Angio16','Angio36'),expected_HR=c(1.137,1.071),expected_lower=c(.878,.823),expected_upper=c(1.473,1.393))
 out<-merge(out,expected,by='signature');out[,pass_rounded_manuscript:=round(estimate,3)==expected_HR&round(lower,3)==expected_lower&round(upper,3)==expected_upper]
 out[,`:=`(model='Surv(PFS,PFS_CNSR) ~ score_per_SD + gene_standardised_APM12 + NIVOLUMAB + APM12:NIVOLUMAB + score:NIVOLUMAB; Efron',fold_file_sha256=digest::digest(file=foldpath,algo='sha256'),model_start=as.character(Sys.time()),fold_file_modified=as.character(file.info(foldpath)$mtime))]
 fwrite(out,file.path(OUT,'results/CM025_I0_reproduction.csv'));print(out);stopifnot(all(out$pass_rounded_manuscript));quit(save='no')
}
rep<-fread(file.path(OUT,'results/CM025_I0_reproduction.csv'));stopifnot(all(rep$pass_rounded_manuscript))
scores<-readRDS(file.path(OUT,'results/CM025_scores.rds'));ints<-list()
for(s in c('Angio16','Angio36'))for(ep in c('PFS','OS'))for(op in paste0('O',1:7)){
 row<-fitint(scores[[s]][,op],ep);row[,`:=`(signature=s,operator=op,endpoint=ep,label='exploratory',endpoint_role=if(ep=='PFS')'primary endpoint in B2-CM025 clinical' else 'secondary endpoint')];ints[[length(ints)+1]]<-row
}
ints<-rbindlist(ints,fill=TRUE);ints[,BH_q:=p.adjust(p,'BH'),by=.(signature,endpoint)];fwrite(ints,file.path(OUT,'results/CM025_interactions.csv'))
frozen<-list();frsumm<-list()
for(s in c('Angio16','Angio36')){
 dd<-data.table(patient=cl$SUBJID,sample=cl$RNA_ID,fold,arm=cl$Arm)
 for(op in c('O1','O2')){
  dd[,(op):=NA_real_];dd[,(paste0(op,'_threshold')):=NA_real_]
  for(f in 1:5){train<-which(fold!=f);test<-which(fold==f);score<-get(op)(X,sets[[s]],ref_idx=train)$scores;cut<-if(all(is.finite(score[train])))median(score[train])else NA_real_;dd[test,(op):=score[test]];dd[test,(paste0(op,'_threshold')):=cut]}
  dd[,(paste0(op,'_distance')):=get(op)-get(paste0(op,'_threshold'))];dd[,(paste0(op,'_high')):=get(paste0(op,'_distance'))>=0]
 }
 dd[,signature:=s];dd[,transition:=fifelse(is.na(O1_high)|is.na(O2_high),'unavailable',fifelse(!O1_high&!O2_high,'stable low',fifelse(!O1_high&O2_high,'low to high',fifelse(O1_high&!O2_high,'high to low','stable high'))))]
 frozen[[s]]<-dd
 k<-sum(dd$O1_high!=dd$O2_high,na.rm=TRUE);nn<-sum(complete.cases(dd[,.(O1_high,O2_high)]));conf<-if(nn)binom.test(k,nn)$conf.int else c(NA,NA)
 frsumm[[s]]<-data.table(signature=s,label='exploratory',n_total=250,n_available=nn,n_unavailable=250-nn,reclassified_n=k,reclassified_per_100=100*k/250,lower=100*conf[1],upper=100*conf[2],CI_denominator=nn,interval='exact binomial conditional on nominated frozen fitted references')
}
frozen<-rbindlist(frozen);fwrite(frozen,file.path(OUT,'results/CM025_frozen_scores_transitions.csv'));fwrite(rbindlist(frsumm),file.path(OUT,'results/CM025_frozen_reclassification.csv'))
kmrow<-function(q,endpoint){
 time<-q[[endpoint]];event<-q[[paste0(endpoint,'_CNSR')]];ok<-is.finite(time)&time>0&!is.na(event);time<-time[ok];event<-event[ok]
 if(!length(time))return(data.table(n_total=nrow(q),n=0,events=0))
 f<-survfit(Surv(time,event)~1,conf.type='log-log');tab<-summary(f)$table
 zero_tail<-isTRUE(tail(f$surv,1)==0);at<-summary(f,times=12,extend=zero_tail)
 supported<-max(time)>=12||zero_tail
 tt<-c(0,f$time[f$time<12],12);ss<-c(1,f$surv[f$time<12]);rmst<-if(supported)sum(diff(tt)*ss)else NA_real_
 data.table(n_total=nrow(q),n=length(time),events=sum(event),median_months=unname(tab['median']),median_lower=unname(tab['0.95LCL']),median_upper=unname(tab['0.95UCL']),
 survival_12m=if(length(at$surv))at$surv else NA_real_,survival_12m_lower=if(length(at$lower))at$lower else NA_real_,survival_12m_upper=if(length(at$upper))at$upper else NA_real_,
 n_at_risk_12m=if(length(at$n.risk))at$n.risk else 0,RMST_12m=rmst,RMST_supported_to_12m=supported,max_followup_months=max(time),KM_zero_tail=zero_tail)
}
outcomes<-orr<-list();groups<-c('stable low','low to high','high to low','stable high','unavailable')
for(s in c('Angio16','Angio36'))for(a in sort(unique(cl$Arm)))for(group in groups){
 ids<-frozen[signature==s&arm==a&transition==group,patient];q<-cl[SUBJID%in%ids]
 for(ep in c('PFS','OS')){row<-kmrow(q,ep);row[,`:=`(signature=s,arm=a,transition=group,endpoint=ep,label='exploratory')];outcomes[[length(outcomes)+1]]<-row}
 ne<-sum(q$ORR=='NE'|is.na(q$ORR));eval<-sum(q$ORR%in%c('CRPR','SD','PD'));resp<-sum(q$ORR=='CRPR',na.rm=TRUE);conf<-if(eval)binom.test(resp,eval)$conf.int else c(NA,NA)
 orr[[length(orr)+1]]<-data.table(signature=s,arm=a,transition=group,label='exploratory',n_total=nrow(q),response_evaluable=eval,not_evaluable=ne,responders=resp,ORR=if(eval)resp/eval else NA_real_,lower=conf[1],upper=conf[2])
}
fwrite(rbindlist(outcomes,fill=TRUE),file.path(OUT,'results/CM025_transition_survival.csv'));fwrite(rbindlist(orr),file.path(OUT,'results/CM025_transition_ORR.csv'))
fwrite(cl[,.(arm=Arm,PFS,PFS_event=PFS_CNSR,OS,OS_event=OS_CNSR,ORR)][,.(n=.N,PFS_events=sum(PFS_event),OS_events=sum(OS_event),ORR_evaluable=sum(ORR%in%c('CRPR','PD','SD')),responders=sum(ORR=='CRPR')),by=arm],file.path(OUT,'results/CM025_endpoint_denominators.csv'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
