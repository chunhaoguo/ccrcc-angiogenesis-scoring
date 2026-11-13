# Author: Chunhao Guo
# Exact R reconstruction of discarded predictions from folder67 core.py/run_policy.py.
source('/LOCAL_PROJECT_ROOT/84_JCO_PO_Upgrade_20261113/05_analysis/A_clinical_consequence/scripts/common.R')
stopifnot(all(read.csv(file.path(OUT,'results/reproduction_checks.csv'))$pass))
lse<-function(x){m<-max(x);m+log(sum(exp(x-m)))}
newton<-function(fun,init,maxit=80,tol=1e-7){
 b<-init
 for(k in 0:(maxit-1)){
  o<-fun(b);if(max(abs(o$g))<tol)return(list(b=b,iterations=k,gradient=max(abs(o$g))))
  step<-solve(o$H,o$g);alpha<-1;accepted<-FALSE
  for(j in 1:40){nb<-b-alpha*step;nval<-fun(nb)$value
   if(is.finite(nval)&&nval<=o$value-1e-4*alpha*sum(o$g*step)+1e-12*(1+abs(o$value))){b<-nb;accepted<-TRUE;break}
   alpha<-alpha*.5
  }
  if(!accepted){if(max(abs(o$g))<1e-5)return(list(b=b,iterations=k,gradient=max(abs(o$g))));stop('Original Newton line-search criterion failed')}
 };stop('Original Newton iteration bound exceeded')
}
pem<-function(X,time,event){
 p<-ncol(X);ng<-2;groups<-c(1,1,2,2);ends<-c(3,6,9,12);starts<-c(0,3,6,9)
 exposure<-sapply(1:4,function(j)pmax(0,pmin(time,ends[j])-starts[j]));events<-matrix(0,nrow(X),4)
 ev<-which(event==1&time<=12);jj<-as.integer(cut(time[ev],c(0,3,6,9,12),include.lowest=TRUE,right=TRUE));events[cbind(ev,jj)]<-1
 D<-colSums(events);suff<-t(events)%*%X;logexp<-log(exposure)
 obj<-function(flat){
  b<-matrix(flat,ng,p,byrow=TRUE);value<-.5*sum(b*b);grad<-b;H<-diag(ng*p)
  for(j in 1:4){if(D[j]==0)next;g<-groups[j];eta<-as.vector(X%*%b[g,]);wl<-logexp[,j]+eta;ls<-lse(wl);w<-exp(wl-ls);mu<-as.vector(crossprod(w,X))
   value<-value+D[j]*ls-sum(suff[j,]*b[g,]);grad[g,]<-grad[g,]+D[j]*mu-suff[j,]
   cov<-crossprod(X,X*w)-tcrossprod(mu);sl<-((g-1)*p+1):(g*p);H[sl,sl]<-H[sl,sl]+D[j]*cov
  };list(value=value,g=as.vector(t(grad)),H=H)
 }
 fit<-newton(obj,rep(0,ng*p));beta<-matrix(fit$b,ng,p,byrow=TRUE);base<-rep(-Inf,4)
 for(j in 1:4)if(D[j]>0)base[j]<-log(D[j])-lse(logexp[,j]+as.vector(X%*%beta[groups[j],]))
 predict<-function(V)exp(sweep(V%*%t(beta[groups,,drop=FALSE]),2,base,'+'))
 list(predict=predict,iterations=fit$iterations,gradient=fit$gradient,events=sum(event*(time<=12)),zero_event_intervals=sum(D==0))
}
integrate12<-function(rates){tail<-rep(0,nrow(rates));for(j in 4:1){r<-rates[,j];integ<-ifelse(r<=1e-12,3,-expm1(-r*3)/r);tail<-integ+exp(-r*3)*tail};tail}
score_train<-function(M,it,method){
 raw<-if(method=='I0')colMeans(M) else {s<-apply(M[,it,drop=FALSE],1,sd);stopifnot(all(is.finite(s)),all(s>1e-8));colMeans(sweep(sweep(M,1,rowMeans(M[,it,drop=FALSE]),'-'),1,s,'/'))}
 (raw-mean(raw[it]))/sd(raw[it])
}
qa<-list();fits<-list();predictions<-list()
add<-function(cohort,check,observed,limit,pass){qa[[length(qa)+1]]<<-data.frame(cohort=cohort,check=check,observed=observed,limit=limit,pass=pass)}
for(tr in names(inputs)){
 d<-inputs[[tr]]$d;M<-inputs[[tr]]$M;fo<-folds[[tr]][[1]];p<-policy[[tr]];p<-p[match(d$ID,p$ID),]
 orig<-read.csv(file.path(path67,'local_only',paste0(tr,'_clinical.csv')));orig<-orig[match(d$ID,orig$ID),]
 add(tr,'Original policy clinical/fold input max error',max(abs(as.matrix(orig[,c('time','event','trt','fold')])-cbind(d$time,d$event,d$trt,fo))),1e-8,max(abs(as.matrix(orig[,c('time','event','trt','fold')])-cbind(d$time,d$event,d$trt,fo)))<1e-8)
 for(tar in names(M)){
  mm<-read.csv(file.path(path67,'local_only',paste0(tr,'_',tar,'.csv')),check.names=FALSE);mm<-mm[match(d$ID,mm$ID),]
  err<-max(abs(as.matrix(mm[,-1])-t(M[[tar]])));add(tr,paste(tar,'matrix input max error'),err,1e-10,err<1e-10)
 }
 archived<-read.csv(file.path(path67,'qa',paste0('B_fit_log_',tr,'.csv')))
 independent<-read.csv(file.path(path67,'local_only',paste0('B_independent_check_inputs_',tr,'.csv')))
 out<-list()
 for(f in 1:5){it<-which(fo!=f);iv<-which(fo==f);a<-score_train(M$APM12,it,'I1');zs<-list()
  for(tar in c('Angio16','Angio36'))for(me in c('I0','I1'))zs[[paste(tar,me,sep='_')]]<-score_train(M[[tar]],it,me)
  nuisance<-cbind(zs$Angio16_I0,zs$Angio36_I0,a)
  for(arm in 0:1){ix<-it[d$trt[it]==arm];fit<-pem(nuisance[ix,,drop=FALSE],d$time[ix],d$event[ix]);pred<-fit$predict(nuisance[iv,,drop=FALSE]);di<-independent[independent$ID%in%d$ID[iv],];j<-match(di$ID,d$ID[iv]);expected<-as.matrix(di[,paste0('q',arm,'_',0:3)])
   err<-max(abs(pred[j,,drop=FALSE]-expected));add(tr,paste('fold',f,'arm',arm,'saved nuisance hazard max error'),err,1e-9,err<1e-9)
  }
  for(tar in c('Angio16','Angio36'))for(me in c('I0','I1')){
   key<-paste(tar,me,sep='_');V<-cbind(zs[[key]],a);values<-matrix(NA_real_,length(iv),2)
   for(arm in 0:1){ix<-it[d$trt[it]==arm];fit<-pem(V[ix,,drop=FALSE],d$time[ix],d$event[ix]);values[,arm+1]<-integrate12(fit$predict(V[iv,,drop=FALSE]))
    old<-archived[archived$fold==f&archived$component==key&archived$arm==arm,]
    pass<-nrow(old)==1&&fit$gradient<1e-5&&fit$events==old$events&&fit$zero_event_intervals==old$zero_event_intervals
    add(tr,paste(f,key,arm,'archived-fit event/support and gradient gate'),fit$gradient,1e-5,pass)
    fits[[length(fits)+1]]<-data.frame(cohort=tr,fold=f,target=tar,method=me,arm=arm,iterations=fit$iterations,archived_iterations=old$iterations,gradient=fit$gradient,archived_gradient=old$gradient,events=fit$events,zero_event_intervals=fit$zero_event_intervals)
   }
   recommended<-as.integer(values[,2]>=values[,1]);nd<-sum(recommended!=p[[paste0('d_',key)]][iv]);add(tr,paste(f,key,'saved recommendation mismatches'),nd,0,nd==0)
   out[[length(out)+1]]<-data.frame(cohort=tr,target=tar,method=me,ID=d$ID[iv],fold=f,RMST12_sunitinib=values[,1],RMST12_experimental=values[,2],benefit_margin=values[,2]-values[,1],recommended_arm=recommended)
  }
 };predictions[[tr]]<-do.call(rbind,out)
}
q<-do.call(rbind,qa);wc(q,'policy_margin_reconstruction_checks.csv');wc(do.call(rbind,fits),'policy_margin_reconstruction_fits.csv');print(q)
if(!all(q$pass))stop('Margin reconstruction gate failed: do not report affected estimates')
pred<-do.call(rbind,predictions)
pred$classification_group<-rep(NA_character_,nrow(pred))
held<-readRDS(file.path(OUT,'results/verified_heldout_frames.rds'))
summ<-list()
for(tr in names(inputs))for(tar in c('Angio36','Angio16'))for(me in c('I0','I1')){
 v<-pred[pred$cohort==tr&pred$target==tar&pred$method==me,];fr<-held[[paste(tr,tar,1,sep='__')]];j<-match(v$ID,fr$ID);v$classification_group<-ifelse(fr$high_I0[j]!=fr$high_I1[j],'Reclassified','Stable')
 pred$classification_group[match(paste(v$cohort,v$target,v$method,v$ID),paste(pred$cohort,pred$target,pred$method,pred$ID))]<-v$classification_group
 for(g in c('Stable','Reclassified')){dd<-v[v$classification_group==g,];z<-dd$benefit_margin
  summ[[length(summ)+1]]<-data.frame(label='exploratory 2',overall_status='exploratory secondary analysis',cohort=tr,target=tar,method=me,partition=1,classification_group=g,n=nrow(dd),
   mean_benefit_margin=mean(z),median_benefit_margin=median(z),q25_benefit_margin=unname(quantile(z,.25)),q75_benefit_margin=unname(quantile(z,.75)),
   min_benefit_margin=min(z),max_benefit_margin=max(z),median_absolute_margin=median(abs(z)),n_recommend_experimental=sum(dd$recommended_arm==1),
   status='recovered exactly from original estimator; all saved recommendation and numerical gates pass',CI='not specified; descriptive distribution of fixed model predictions')
 }
}
wc(pred,'policy_model_predictions_reconstructed_LOCAL_ONLY.csv');wc(do.call(rbind,summ),'predicted_benefit_margin_summary.csv')
status<-expand.grid(cohort=names(inputs),target=c('Angio36','Angio16'))|>transform(label='exploratory 2',estimand='Model-predicted 12-month benefit margins, reclassified versus stable',status='completed by exact original-estimator reconstruction',reason='All saved recommendation, nuisance numerical parity and archived-fit checks pass; see policy_margin_reconstruction_checks.csv')
wc(status,'benefit_margin_input_status.csv')
# Preserve the initial assessment filename while explicitly superseding its status.
if(file.exists(file.path(OUT,'results/benefit_margin_NOT_RUN.csv')))wc(status,'benefit_margin_NOT_RUN.csv')
writeLines('PASS: original-estimator prediction recovery; all saved decisions, independent nuisance predictions, and archived-fit gates reproduced.',file.path(OUT,'logs/POLICY_MARGIN_RECOVERY_PASS.txt'))
cat('Benefit margin reconstruction passed all gates.\n')
