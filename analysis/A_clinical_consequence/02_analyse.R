# Author: Chunhao Guo
source('/LOCAL_PROJECT_ROOT/84_JCO_PO_Upgrade_20261113/05_analysis/A_clinical_consequence/scripts/common.R')
stopifnot(file.exists(file.path(OUT,'logs/REPRODUCTION_PASS.txt')),
 all(read.csv(file.path(OUT,'results/reproduction_checks.csv'))$pass))
frames<-readRDS(file.path(OUT,'results/verified_heldout_frames.rds'))
groups<-c('Stable low','Low to high','High to low','Stable high')
arms<-function(tr,a)ifelse(a==0,'sunitinib',if(tr=='JAVELIN')'avelumab + axitinib' else 'atezolizumab + bevacizumab')

# Verify clinical response/time/event/arm against the retained IMmotion151 source.
cl<-read.csv(file.path(ROOT,'21_shared_datasets/IMmotion151_20201106/raw/IMmotion151_clinical_anon_20201106(1).csv'),skip=3,na.strings=c('NA',''))
d<-inputs$IMmotion151$d;cl<-cl[match(d$ID,cl$PATIENT_ID),]
stopifnot(!anyNA(cl$PATIENT_ID),identical(cl$PATIENT_ID,d$ID),identical(cl$RNASEQ_SAMPLE_ID,d$RNA_ID),
 max(abs(cl$PFS_MONTHS-d$time))<1e-7,all(as.integer(cl$PFS_CENSOR)==d$event),
 all((cl$ARM!='sunitinib')==d$trt),all(cl$OBJECTIVE_RESPONSE==d$response))
wc(data.frame(cohort='IMmotion151',n=nrow(d),RNA_key_match=TRUE,patient_key_match=TRUE,
 PFS_time_match=TRUE,PFS_event_match=TRUE,arm_match=TRUE,response_match=TRUE),'clinical_join_audit.csv')

kmstats<-function(dd,tau=12) {
 ans<-data.frame(n=nrow(dd),events=sum(dd$event),n_at_risk_12=sum(dd$time>=tau),
  median_PFS=NA_real_,median_low=NA_real_,median_high=NA_real_,median_status='not estimable',
  median_CI_status='not estimable',PFS12=NA_real_,PFS12_low=NA_real_,PFS12_high=NA_real_,
  PFS12_status='not estimable',RMST12=NA_real_,RMST12_SE=NA_real_,RMST12_status='not estimable')
 if(!nrow(dd))return(ans)
 fit<-survfit(Surv(time,event)~1,data=dd,conf.type='log-log')
 tab<-summary(fit,rmean=tau)$table
 ans$median_PFS<-unname(tab['median']);ans$median_low<-unname(tab['0.95LCL']);ans$median_high<-unname(tab['0.95UCL'])
 if(is.finite(ans$median_PFS))ans$median_status<-'estimable'
 if(all(is.finite(c(ans$median_low,ans$median_high))))ans$median_CI_status<-'estimable'
 support<-max(dd$time)>=tau || tail(fit$surv,1)==0
 if(support){
  ss<-summary(fit,times=tau,extend=TRUE)
  ans$PFS12<-ss$surv[1];ans$PFS12_low<-ss$lower[1];ans$PFS12_high<-ss$upper[1]
  if(ans$PFS12==0){ans$PFS12_low<-0;ans$PFS12_high<-0}
  ans$PFS12_status<-'estimable';ans$RMST12<-unname(tab['rmean']);ans$RMST12_SE<-unname(tab['se(rmean)']);ans$RMST12_status<-'estimable'
 }
 ans
}
rmst<-function(dd,tau=12){
 if(!nrow(dd))return(NA_real_)
 ff<-survfit(Surv(time,event)~1,data=dd)
 if(max(dd$time)<tau && tail(ff$surv,1)>0)return(NA_real_)
 ii<-which(ff$time<tau)
 sum(diff(c(0,ff$time[ii],tau))*c(1,ff$surv[ii]))
}
exact<-function(k,n)if(n>0)unname(binom.test(k,n)$conf.int) else c(NA_real_,NA_real_)
cox_est<-function(dd,var){
 ans<-data.frame(HR=NA_real_,HR_low=NA_real_,HR_high=NA_real_,HR_status='not estimable',warning='')
 if(!nrow(dd)||length(unique(dd[[var]]))<2||sum(dd$event)==0)return(ans)
 warnings<-character()
 fit<-tryCatch(withCallingHandlers(coxph(as.formula(paste('Surv(time,event)~',var)),data=dd,ties='efron'),
 warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error')){ans$warning<-conditionMessage(fit);return(ans)}
 b<-unname(coef(fit)[1]);se<-sqrt(vcov(fit)[1,1]);est<-exp(c(b,b-1.95996398454*se,b+1.95996398454*se))
 ans$warning<-paste(warnings,collapse='; ')
 if(length(warnings)||!all(is.finite(est))||se<=0)return(ans)
 ans[1,1:3]<-est;ans$HR_status<-'estimable';ans
}
attachclinical<-function(fr){
 d<-inputs[[fr$trial[1]]]$d;ix<-match(fr$ID,d$ID);stopifnot(!anyNA(ix),!anyDuplicated(fr$ID))
 fr$time<-d$time[ix];fr$event<-d$event[ix];fr$trt<-d$trt[ix]
 fr$response<-if('response'%in%names(d))d$response[ix] else NA_character_
 fr$arm<-arms(fr$trial[1],fr$trt);fr
}
outcomes<-response<-reclassification<-treatment<-interactions<-distance<-list()
run_groups<-function(fr,label,design){
 fr<-attachclinical(fr);tr<-fr$trial[1];tar<-fr$target[1];re<-fr$partition[1]
 base<-data.frame(module='A',label=label,overall_status='exploratory secondary analysis',cohort=tr,target=tar,
  design=design,partition=re,threshold_probability=fr$threshold_probability[1])
 for(a in 0:1)for(g in groups){
  dd<-fr[fr$trt==a&fr$transition==g,];b<-cbind(base,arm=arms(tr,a),trt=a,transition=g)
  outcomes[[length(outcomes)+1]]<<-cbind(b,kmstats(dd))
  if(tr=='IMmotion151'){
   cats<-table(factor(dd$response,levels=c('CR','PR','SD','PD','NE')));nr<-sum(dd$response%in%c('CR','PR'));n<-nrow(dd);ne<-sum(dd$response=='NE'|is.na(dd$response));ev<-n-ne
   ci<-exact(nr,n);cie<-exact(nr,ev)
   response[[length(response)+1]]<<-cbind(b,n=n,CR=unname(cats['CR']),PR=unname(cats['PR']),SD=unname(cats['SD']),PD=unname(cats['PD']),NE=ne,
    response_evaluable=ev,responders=nr,ORR=if(n>0)nr/n else NA_real_,ORR_low=ci[1],ORR_high=ci[2],
    ORR_status=if(n>0)'estimable' else 'not estimable',ORR_evaluable=if(ev>0)nr/ev else NA_real_,ORR_evaluable_low=cie[1],ORR_evaluable_high=cie[2])
  }
 }
 # Entire cohort and by randomised arm; unavailable stays in denominator.
 for(a in c(-1,0,1)){
  dd<-if(a==-1)fr else fr[fr$trt==a,];un<-sum(dd$transition=='Unavailable');nc<-sum(dd$high_I0!=dd$high_I1,na.rm=TRUE);ci<-exact(nc,nrow(dd))*100
  reclassification[[length(reclassification)+1]]<<-cbind(base,arm=if(a==-1)'All randomised arms' else arms(tr,a),n=nrow(dd),unavailable=un,reclassified=nc,
   per100=100*nc/nrow(dd),low=ci[1],high=ci[2],CI_method='Clopper-Pearson conditional on saved classifications')
 }
 fr$discordant<-as.integer(fr$high_I0!=fr$high_I1)
 base_treat<-base;base_treat$classification_definition_label<-label
 base_treat$label<-if(label=='primary')'exploratory 2' else 'sensitivity'
 for(g in 0:1){dd<-fr[!is.na(fr$discordant)&fr$discordant==g,];ct<-cox_est(dd,'trt')
  treatment[[length(treatment)+1]]<<-cbind(base_treat,
   classification_group=if(g==1)'Discordant / reclassified' else 'Concordant / stable',n=nrow(dd),events=sum(dd$event),n_experimental=sum(dd$trt==1),n_sunitinib=sum(dd$trt==0),ct)
 }
 warn<-character();ff<-tryCatch(withCallingHandlers(coxph(Surv(time,event)~trt*discordant,data=fr,ties='efron'),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 p<-b<-se<-NA_real_;status<-'not estimable'
 if(!inherits(ff,'error')&&'trt:discordant'%in%names(coef(ff))){b<-unname(coef(ff)['trt:discordant']);se<-sqrt(vcov(ff)['trt:discordant','trt:discordant']);if(!length(warn)&&is.finite(b)&&is.finite(se)&&se>0){p<-2*pnorm(-abs(b/se));status<-'estimable'}}
 interactions[[length(interactions)+1]]<<-cbind(base_treat,n=nrow(fr),
  interaction_log_HR=b,SE=se,ratio_of_HRs=if(status=='estimable')exp(b) else NA_real_,
  low=if(status=='estimable')exp(b-1.95996398454*se) else NA_real_,high=if(status=='estimable')exp(b+1.95996398454*se) else NA_real_,
  p_interaction=p,status=status,warning=if(inherits(ff,'error'))conditionMessage(ff) else paste(warn,collapse='; '))
 distance[[length(distance)+1]]<<-cbind(label=label,design=design,fr)
 invisible(fr)
}
for(tr in names(inputs))for(tar in c('Angio36','Angio16'))for(re in 1:5){
 fr<-frames[[paste(tr,tar,re,sep='__')]]
 run_groups(fr,if(tar=='Angio36'&&re==1)'primary' else 'sensitivity','training-frozen held-out median')
}
for(tr in names(inputs))for(tar in c('Angio36','Angio16')){
 M<-inputs[[tr]]$M[[tar]];d<-inputs[[tr]]$d;z0<-full_scores(M,'I0');z1<-full_scores(M,'I1');c0<-median(z0);c1<-median(z1)
 fr<-data.frame(trial=tr,target=tar,partition=0,threshold_probability=.5,ID=d$ID,fold=NA_integer_,score_I0=z0,score_I1=z1,threshold_I0=c0,threshold_I1=c1,
  signed_distance_I0=z0-c0,signed_distance_I1=z1-c1,absolute_distance_I0=abs(z0-c0),absolute_distance_I1=abs(z1-c1),high_I0=z0>=c0,high_I1=z1>=c1,transition=transition(z0>=c0,z1>=c1))
 run_groups(fr,'sensitivity','full-cohort median descriptive comparison')
}
wc(do.call(rbind,outcomes),'transition_PFS_outcomes.csv');wc(do.call(rbind,response),'transition_ORR_IMmotion151.csv')
wc(do.call(rbind,reclassification),'median_reclassification.csv');wc(do.call(rbind,treatment),'classification_group_treatment_HR.csv');wc(do.call(rbind,interactions),'classification_treatment_interaction.csv')
wc(do.call(rbind,distance),'patient_score_threshold_distances_LOCAL_ONLY.csv')

# Expression-only threshold sensitivity; no outcome fitting for extra cutoffs.
threshold_rows<-list();threshold_patient<-list()
for(tr in names(inputs))for(tar in c('Angio36','Angio16'))for(re in 1:5)for(pr in c(.25,1/3,.5,2/3,.75)){
 fr<-held_frame(tr,tar,re,pr);n<-nrow(fr);nc<-sum(fr$high_I0!=fr$high_I1,na.rm=TRUE);ci<-exact(nc,n)*100
 threshold_rows[[length(threshold_rows)+1]]<-data.frame(label='sensitivity',overall_status='exploratory secondary analysis',cohort=tr,target=tar,partition=re,
  threshold_probability=pr,n=n,unavailable=sum(fr$transition=='Unavailable'),reclassified=nc,per100=100*nc/n,low=ci[1],high=ci[2])
 for(fo in 1:5){v<-fr[fr$fold==fo,];threshold_patient[[length(threshold_patient)+1]]<-data.frame(cohort=tr,target=tar,partition=re,fold=fo,
  threshold_probability=pr,training_n=n-nrow(v),test_n=nrow(v),threshold_I0=v$threshold_I0[1],threshold_I1=v$threshold_I1[1])}
}
wc(do.call(rbind,threshold_rows),'threshold_sensitivity_reclassification.csv');wc(do.call(rbind,threshold_patient),'training_threshold_manifest.csv')
cat('Primary and classification sensitivities complete.\n');flush.console()

# Exploratory 1: fixed saved policy-discordant subgroup, both signatures.
policy_results<-policy_draws<-risk_tables<-policy_patient<-list()
set.seed(20261113)
for(tr in names(inputs))for(tar in c('Angio36','Angio16')){
 d<-inputs[[tr]]$d;p<-policy[[tr]];p<-p[match(d$ID,p$ID),]
 a0<-p[[paste0('d_',tar,'_I0')]];a1<-p[[paste0('d_',tar,'_I1')]];ix<-which(a0!=a1)
 dd<-d[ix,c('ID','time','event','trt')];dd$received_I1<-as.integer(dd$trt==a1[ix]);dd$received_I0<-as.integer(dd$trt==a0[ix])
 stopifnot(all(dd$received_I0+dd$received_I1==1))
 ct<-cox_est(dd,'received_I1');r0<-rmst(dd[dd$received_I1==0,]);r1<-rmst(dd[dd$received_I1==1,]);delta<-r1-r0
 draw<-rep(NA_real_,2000);why<-rep('estimable',2000)
 for(b in 1:2000){ii<-sample.int(nrow(dd),nrow(dd),replace=TRUE);db<-dd[ii,]
  if(length(unique(db$received_I1))<2){why[b]<-'absent agreement group';next}
  draw[b]<-rmst(db[db$received_I1==1,])-rmst(db[db$received_I1==0,])
  if(!is.finite(draw[b]))why[b]<-'insufficient 12-month support'
 }
 ci<-if(sum(is.finite(draw))>0)unname(quantile(draw,c(.025,.975),na.rm=TRUE,type=7)) else c(NA_real_,NA_real_)
 policy_results[[length(policy_results)+1]]<-cbind(data.frame(label='exploratory 1',overall_status='exploratory secondary analysis',cohort=tr,target=tar,
  partition=1,n=nrow(dd),n_received_I0=sum(dd$received_I0),n_received_I1=sum(dd$received_I1),events_received_I0=sum(dd$event[dd$received_I0==1]),events_received_I1=sum(dd$event[dd$received_I1==1])),ct,
  RMST12_received_I0=r0,RMST12_received_I1=r1,RMST12_difference_I1_minus_I0=delta,RMST12_low=ci[1],RMST12_high=ci[2],
  RMST_status=if(is.finite(delta))'estimable' else 'not estimable',bootstrap_requested=2000,bootstrap_estimable=sum(is.finite(draw)),bootstrap_undefined=sum(!is.finite(draw)),
  bootstrap_absent_group=sum(why=='absent agreement group'),bootstrap_insufficient_support=sum(why=='insufficient 12-month support'))
 policy_draws[[length(policy_draws)+1]]<-data.frame(cohort=tr,target=tar,draw=1:2000,RMST12_difference=draw,status=why)
 policy_patient[[length(policy_patient)+1]]<-cbind(cohort=tr,target=tar,dd,recommended_I0=a0[ix],recommended_I1=a1[ix])
 for(g in 0:1)for(t in seq(0,24,3))risk_tables[[length(risk_tables)+1]]<-data.frame(label='exploratory 1',cohort=tr,target=tar,received_rule=if(g==1)'I1 gene-standardised score' else 'I0 mean-expression score',month=t,n_at_risk=sum(dd$received_I1==g&dd$time>=t))
 png(file.path(OUT,'plots',paste0('KM_policy_discordant_',tr,'_',tar,'.png')),width=1400,height=1000,res=140,type='quartz')
 layout(matrix(1:2,ncol=1),heights=c(3,1));par(mar=c(4,4,4,2))
 fit<-survfit(Surv(time,event)~received_I1,data=dd,conf.type='log-log')
 plot(fit,col=c('#2563A6','#C04B3B'),lwd=2,conf.int=FALSE,xlim=c(0,24),xlab='Months',ylab='Progression-free survival',
  main=paste('Exploratory 1:',tr,tar,'policy-discordant subgroup'))
 abline(v=12,lty=3,col='grey50');legend('topright',c(paste0('Received I0 rule arm (n=',sum(dd$received_I0),')'),paste0('Received I1 rule arm (n=',sum(dd$received_I1),')')),col=c('#2563A6','#C04B3B'),lwd=2,bty='n',cex=.9)
 par(mar=c(2,4,0,2));plot(c(0,24),c(-.5,2.5),type='n',axes=FALSE,xlab='',ylab='');axis(1,at=seq(0,24,3))
 for(g in 0:1){text(0,2-g,if(g==0)'I0 at risk' else 'I1 at risk',adj=1,xpd=NA,cex=.8);text(seq(0,24,3),1.5-g,sapply(seq(0,24,3),function(t)sum(dd$received_I1==g&dd$time>=t)),cex=.85)}
 dev.off();cat('Exploratory 1 complete:',tr,tar,'\n');flush.console()
}
wc(do.call(rbind,policy_results),'policy_discordant_PFS_contrasts.csv');wc(do.call(rbind,policy_draws),'policy_RMST_bootstrap_draws.csv')
wc(do.call(rbind,risk_tables),'policy_KM_numbers_at_risk.csv');wc(do.call(rbind,policy_patient),'policy_discordant_patients_LOCAL_ONLY.csv')
wc(expand.grid(cohort=names(inputs),target=c('Angio36','Angio16'))|>transform(label='exploratory 2',estimand='Model-predicted 12-month benefit margins, reclassified versus stable',status=if(file.exists(file.path(OUT,'logs/POLICY_MARGIN_RECOVERY_PASS.txt')))'completed by exact original-estimator reconstruction' else 'pending exact reconstruction',reason='Policy predictions not retained; recover only with original model and verification gates in 02b_recover_policy_margins.R'),'benefit_margin_input_status.csv')

# Quick primary transition-group bars by randomised arm.
primary<-do.call(rbind,outcomes);primary<-primary[primary$label=='primary',]
png(file.path(OUT,'plots/primary_transition_group_counts.png'),width=1500,height=800,res=140,type='quartz');par(mfrow=c(1,2),mar=c(6,4,4,1))
for(tr in names(inputs)){v<-primary[primary$cohort==tr,];mat<-matrix(v$n,nrow=4,ncol=2,dimnames=list(groups,arms(tr,0:1)))
 bp<-barplot(t(mat),beside=TRUE,col=c('#2563A6','#C04B3B'),names.arg=c('Stable\nlow','Low to\nhigh','High to\nlow','Stable\nhigh'),ylab='Patients',ylim=c(0,max(mat)*1.2),main=paste('A-primary:',tr),las=1)
 text(bp,t(mat)+max(mat)*.035,labels=t(mat),cex=.8);legend('topright',legend=colnames(mat),fill=c('#2563A6','#C04B3B'),bty='n',cex=.7)
};dev.off()
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
cat('Module A estimation completed.\n')
