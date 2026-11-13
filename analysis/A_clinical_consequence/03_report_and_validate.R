# Author: Chunhao Guo
source('/LOCAL_PROJECT_ROOT/84_JCO_PO_Upgrade_20261113/05_analysis/A_clinical_consequence/scripts/common.R')
rd<-function(f)read.csv(file.path(OUT,'results',f))
oo<-rd('transition_PFS_outcomes.csv');orr<-rd('transition_ORR_IMmotion151.csv');rc<-rd('median_reclassification.csv')
pc<-rd('policy_discordant_PFS_contrasts.csv');tx<-rd('classification_group_treatment_HR.csv');it<-rd('classification_treatment_interaction.csv')
th<-rd('threshold_sensitivity_reclassification.csv');ch<-rd('reproduction_checks.csv');draws<-rd('policy_RMST_bootstrap_draws.csv')
dist<-rd('patient_score_threshold_distances_LOCAL_ONLY.csv');risk<-rd('policy_KM_numbers_at_risk.csv')
margin<-rd('predicted_benefit_margin_summary.csv');mchecks<-rd('policy_margin_reconstruction_checks.csv');mpred<-rd('policy_model_predictions_reconstructed_LOCAL_ONLY.csv')
fmt<-function(x,d=2)ifelse(is.finite(x),formatC(x,format='f',digits=d),'not estimable')
ci<-function(lo,hi,d=2)paste0('[',fmt(lo,d),', ',fmt(hi,d),']')
mdtable<-function(df){df[]<-lapply(df,as.character);c(paste0('| ',paste(names(df),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(df)),collapse=' | '),' |'),apply(df,1,function(x)paste0('| ',paste(x,collapse=' | '),' |')))}
headlines<-list();add<-function(estimand,label,cohort,n,estimate,interval,src){headlines[[length(headlines)+1]]<<-data.frame(estimand=estimand,label=label,cohort=cohort,n=n,estimate=estimate,`95% CI`=interval,`source file`=paste0('results/',src),check.names=FALSE)}
for(tr in names(inputs)){
 z<-rc[rc$cohort==tr&rc$label=='primary'&rc$arm=='All randomised arms',]
 add('Angio36 reclassified patients per 100, nominated partition','primary',tr,z$n,fmt(z$per100),ci(z$low,z$high),'median_reclassification.csv')
}
for(j in seq_len(nrow(pc))){z<-pc[j,];role<-if(z$target=='Angio16')' (secondary signature)' else ''
 add(paste0(z$target,' policy-discordant PFS HR, received I1 vs I0 rule arm',role),'exploratory 1',z$cohort,z$n,fmt(z$HR),ci(z$HR_low,z$HR_high),'policy_discordant_PFS_contrasts.csv')
 add(paste0(z$target,' policy-discordant 12-month RMST difference, I1 minus I0 (months)',role),'exploratory 1',z$cohort,z$n,fmt(z$RMST12_difference_I1_minus_I0),ci(z$RMST12_low,z$RMST12_high),'policy_discordant_PFS_contrasts.csv')
}
for(tr in names(inputs)){
 z<-it[it$cohort==tr&it$label=='exploratory 2',]
 add('Angio36 treatment HR ratio, discordant / concordant classification','exploratory 2',tr,z$n,fmt(z$ratio_of_HRs,3),ci(z$low,z$high,3),'classification_treatment_interaction.csv')
 z<-rc[rc$cohort==tr&rc$target=='Angio16'&rc$partition==1&rc$arm=='All randomised arms',]
 add('Angio16 reclassified patients per 100, nominated partition','sensitivity',tr,z$n,fmt(z$per100),ci(z$low,z$high),'median_reclassification.csv')
}
h<-do.call(rbind,headlines);wc(h,'headline_table.csv')

# Independent numerical and completeness checks; no new outcome estimands.
qa<-list();check<-function(name,pass,details=''){qa[[length(qa)+1]]<<-data.frame(check=name,pass=isTRUE(pass),details=details)}
check('Reproduction checks',all(ch$pass),paste(nrow(ch),'checks'))
check('All 24 median designs have 8 transition-by-arm rows',nrow(oo)==192)
check('All 12 IMmotion151 designs have 8 ORR rows',nrow(orr)==96)
check('Exactly 16 A-primary PFS rows',sum(oo$label=='primary')==16)
check('Exactly 8 A-primary ORR rows',sum(orr$label=='primary')==8)
check('All threshold boundaries/signatures/partitions reported',nrow(th)==100)
check('All 4 policy subgroup contrasts reported',nrow(pc)==4)
check('Exactly 2000 attempted bootstrap draws per contrast',all(table(draws$cohort,draws$target)==2000))
check('All 48 classification-group treatment HR rows reported',nrow(tx)==48)
check('All 24 classification-by-treatment interaction rows reported',nrow(it)==24)
check('All 18588 score-distance rows retained',nrow(dist)==18588)
check('All policy-margin recovery gates pass',all(mchecks$pass),paste(nrow(mchecks),'checks'))
check('All 6196 recovered policy-prediction rows retained',nrow(mpred)==6196&&!anyNA(mpred$classification_group))
check('All 16 predicted-margin summaries retained',nrow(margin)==16)
check('No duplicate transition-by-arm keys',!anyDuplicated(oo[,c('cohort','target','design','partition','trt','transition')]))
for(tr in names(inputs))for(tar in c('Angio36','Angio16'))for(re in 0:5){
 z<-oo[oo$cohort==tr&oo$target==tar&oo$partition==re,]
 check(paste(tr,tar,'partition',re,'PFS counts/events conserved'),sum(z$n)==nrow(inputs[[tr]]$d)&&sum(z$events)==sum(inputs[[tr]]$d$event))
 if(tr=='IMmotion151'){
  v<-orr[orr$cohort==tr&orr$target==tar&orr$partition==re,]
  check(paste(tr,tar,'partition',re,'response categories conserved'),sum(v$n)==823&&sum(v$CR)==32&&sum(v$PR)==257&&sum(v$SD)==318&&sum(v$PD)==151&&sum(v$NE)==65)
 }
}
for(j in seq_len(nrow(pc))){z<-pc[j,];dd<-draws[draws$cohort==z$cohort&draws$target==z$target,]
 check(paste(z$cohort,z$target,'bootstrap accounting'),sum(is.finite(dd$RMST12_difference))==z$bootstrap_estimable)
 v<-dist[dist$trial==z$cohort&dist$target==z$target&dist$partition==1&dist$design=='training-frozen held-out median',]
 # Cross-check manual KM integration against survival's restricted-mean implementation.
 for(a in 0:1)for(g in c('Stable low','Low to high','High to low','Stable high')){
  d<-v[v$trt==a&v$transition==g,];if(nrow(d)&&max(d$time)>=12){ff<-survfit(Surv(time,event)~1,data=d);idx<-ff$time<12
   manual<-sum(diff(c(0,ff$time[idx],12))*c(1,ff$surv[idx]));sv<-unname(summary(ff,rmean=12)$table['rmean'])
   check(paste(z$cohort,z$target,a,g,'RMST integration parity'),abs(manual-sv)<1e-10)
  }
 }
}
q<-do.call(rbind,qa);wc(q,'validation_checks.csv');stopifnot(all(q$pass));print(q)

lines<-c(mdtable(h),'',
 '# Module A — clinical consequence of implementation-dependent reclassification','',
 '**Author and last modified by:** Chunhao Guo. Date: 2026-11-13.','',
 '**Status:** A-primary, A-exploratory 1 and 2, and A-sensitivity completed. The model-predicted benefit margins were absent from saved policy CSVs but recovered with an exact R reconstruction of the original estimator, gated against saved decisions and numerical predictions (input-route deviation recorded before reconstruction). Every result is a post hoc exploratory secondary analysis overall; “primary” refers only to the locked upgrade estimand, not confirmatory evidence.','',
 '**Definitions:** I0 = mean-expression score; I1 = gene-standardised score; research high = score at or above the training median (R quantile type 7; ties high). Transition direction is I0 → I1. Experimental arms are avelumab + axitinib in JAVELIN and atezolizumab + bevacizumab in IMmotion151; comparator is sunitinib. Time is months. Patient is the analysis unit.','',
 '## Reproduction checks completed before any new estimate','',
 paste0('**PASS: ',nrow(ch),' numerical, gene-list, threshold, fold, and patient-transition checks.** Source: `results/reproduction_checks.csv`; execution: `logs/01_reproduce.log`. The saved transition of every patient was reproduced under all five partitions for both signatures before new outcome estimation. Thresholds reproduced within 1e-8 using the saved fold-specific reference parameters. No branch proceeded on unverified inputs.'),'',
 '| Check | JAVELIN | IMmotion151 | Source |','| --- | --- | --- | --- |',
 '| Angio36 held-out counts, partitions 1–5 | 42, 54, 51, 52, 49 / 726 | 85, 80, 82, 83, 81 / 823 | `reproduction_checks.csv`; saved folder48 patient ranks and folder70 held-out summary |',
 '| Angio36 held-out rate range | 5.79–7.44 per 100 | 9.72–10.33 per 100 | Same sources |',
 '| Angio36 full-cohort median disagreement | 48/726 | 80/823 | Folder48 saved full-cohort patient ranks, verified against complete-score recalculation |',
 '| Angio36 saved 12-month policy disagreement | 40/726 | 76/823 | Folder67 `B_patient_policy_*.csv`, tau=12 |',
 '| Angio16 saved 12-month policy disagreement | 9/726 | 15/823 | Same policy files |','',
 'All manuscript headline counts match. PFS totals are 358 events in 726 JAVELIN patients and 529 in 823 IMmotion151 patients. IMmotion151 patient/RNA keys, PFS time/event, randomised arm and response also match the original retained clinical CSV exactly: `results/clinical_join_audit.csv`.','',
 '## A-primary — Angio36 nominated-partition observed outcomes','',
 'The four transitions retain every RNA-linked patient in each randomised arm. Median PFS intervals invert pointwise log-log KM bands; PFS at 12 months uses log-log 95% CI. RMST is the KM area through 12 months. A non-estimable bound or horizon is printed explicitly. No curve is extended beyond available follow-up while survival remains positive. Source for every cell below: `results/transition_PFS_outcomes.csv`, `label=primary`.','')
z<-oo[oo$label=='primary',]
tab<-data.frame(cohort=z$cohort,arm=z$arm,transition=z$transition,n=z$n,events=z$events,`at risk 12m`=z$n_at_risk_12,
 `median PFS [95% CI], months`=paste(fmt(z$median_PFS),ci(z$median_low,z$median_high)),
 `12-month PFS [95% CI], %`=paste(fmt(100*z$PFS12,1),ci(100*z$PFS12_low,100*z$PFS12_high,1)),`12-month RMST, months`=fmt(z$RMST12),check.names=FALSE)
lines<-c(lines,mdtable(tab),'',
 'Reclassified patients per 100 with exact conditional binomial intervals, both overall and by arm, are in `results/median_reclassification.csv`. These intervals condition on the saved classification rules and do not include reference-fitting uncertainty.','',
 '### A-primary — IMmotion151 objective response','',
 'ORR is documented CR/PR divided by all patients in that transition-by-arm group; NE remains in the denominator. Exact Clopper–Pearson 95% CI. Evaluable-only ORR and exact intervals are retained as denominator information in `results/transition_ORR_IMmotion151.csv`; no response is imputed. Source for every cell below: that CSV, `label=primary`.','')
z<-orr[orr$label=='primary',]
lines<-c(lines,mdtable(data.frame(arm=z$arm,transition=z$transition,n=z$n,`CR/PR`=paste(z$CR,z$PR,sep='/'),SD=z$SD,PD=z$PD,NE=z$NE,
 `ORR [95% CI], %`=paste(fmt(100*z$ORR,1),ci(100*z$ORR_low,100*z$ORR_high,1)),check.names=FALSE)),'',
 '### A-primary — distance from frozen thresholds','',
 'Each patient’s score, frozen method-specific threshold, signed distance and absolute distance under I0 and I1 are saved in `results/patient_score_threshold_distances_LOCAL_ONLY.csv`. Distance is in the corresponding training score SD units, not a calibrated clinical probability; the two implementations have separate training scales. The file also retains PFS/arm, partition and transition for audit, is local participant-level material, and must not be redistributed without the data holders’ permission. The following summaries use only the nominated Angio36 held-out partition.','')
ds<-list()
for(tr in names(inputs))for(a in 0:1)for(g in c('Stable low','Low to high','High to low','Stable high')){
 v<-dist[dist$trial==tr&dist$target=='Angio36'&dist$partition==1&dist$trt==a&dist$transition==g&dist$design=='training-frozen held-out median',]
 ds[[length(ds)+1]]<-data.frame(label='primary',cohort=tr,arm=if(a==0)'sunitinib' else if(tr=='JAVELIN')'avelumab + axitinib' else 'atezolizumab + bevacizumab',transition=g,n=nrow(v),
  median_absolute_distance_I0=median(v$absolute_distance_I0),median_absolute_distance_I1=median(v$absolute_distance_I1))
}
ds<-do.call(rbind,ds);wc(ds,'primary_distance_summary.csv');ds[,6:7]<-lapply(ds[,6:7],fmt)
lines<-c(lines,mdtable(ds),'',
 '## A-exploratory 1 — randomised comparison within saved policy-discordant patients','',
 'HR compares patients who received the randomised arm recommended by the I1 rule with those who received the arm recommended by I0; lower HR favours the I1-agreement group. RMST difference is I1 minus I0, so positive differences favour I1 agreement. Each discordant patient agrees with exactly one rule. The subgroup is defined by the saved cross-fitted learned hypothetical policies; it is distinct from median-threshold reclassification. Angio16 is the secondary signature. Source: `results/policy_discordant_PFS_contrasts.csv`.','')
tab<-data.frame(cohort=pc$cohort,signature=pc$target,n=pc$n,`n I0/I1 agreement`=paste(pc$n_received_I0,pc$n_received_I1,sep='/'),
 `events I0/I1 agreement`=paste(pc$events_received_I0,pc$events_received_I1,sep='/'),`HR [95% CI]`=paste(fmt(pc$HR),ci(pc$HR_low,pc$HR_high)),
 `RMST I0/I1, months`=paste(fmt(pc$RMST12_received_I0),fmt(pc$RMST12_received_I1),sep='/'),`RMST difference [95% percentile interval], months`=paste(fmt(pc$RMST12_difference_I1_minus_I0),ci(pc$RMST12_low,pc$RMST12_high)),
 `estimable / 2000 bootstrap draws`=paste0(pc$bootstrap_estimable,'/2000'),check.names=FALSE)
lines<-c(lines,mdtable(tab),'',
 '**Bootstrap qualification:** all 2,000 ordinary patient draws per contrast were attempted, with seed 20261113 and the saved rules held fixed. Percentile intervals use estimable draws only; they are conditional intervals when draws fail, without a nominal-coverage claim. No replacement draws or tail extrapolation were used. JAVELIN Angio36: 157 undefined draws, all insufficient 12-month support; JAVELIN Angio16: 507 undefined (19 absent agreement group, 488 insufficient support); IMmotion151 Angio36: 0 undefined; IMmotion151 Angio16: 7 undefined (6 absent group, 1 insufficient support). All draws and reasons are in `results/policy_RMST_bootstrap_draws.csv`. The particularly incomplete JAVELIN Angio16 bootstrap cannot be treated as a stable unconditional 95% CI.','',
 'The observed Angio36 HRs point towards worse PFS in the I1-agreement groups in both trials, but all four HR intervals include 1, and all four RMST intervals include 0. These small-subgroup, post hoc comparisons cannot identify individual treatment benefit, prove superiority of either implementation, or describe patients as incorrectly treated. Saved policies depend on outcome-fitted cross-fitted models; this is not an independently validated clinical rule.','',
 'Quick KM diagnostic PNGs are in `plots/KM_policy_discordant_{JAVELIN,IMmotion151}_{Angio36,Angio16}.png`; exact numbers at risk (0–24 months, every 3 months) are in `results/policy_KM_numbers_at_risk.csv`. Curves are displayed through 24 months, while the prespecified RMST horizon remains 12 months.','',
 '## A-exploratory 2 — treatment HR within classification-discordance groups','',
 'Experimental versus sunitinib PFS HR is estimated separately in concordant/stable and discordant/reclassified median-classification groups. The one-df interaction is the Wald test from `Surv(time,event) ~ treatment * discordance`, Efron ties; ratio of HRs = discordant treatment HR / concordant treatment HR. These unadjusted comparisons are exploratory, with no multiplicity-adjusted confirmatory claim. Source: `results/classification_group_treatment_HR.csv` and `results/classification_treatment_interaction.csv`.','')
tt<-list()
for(j in seq_len(nrow(it))){z<-it[j,];v<-tx[tx$cohort==z$cohort&tx$target==z$target&tx$partition==z$partition,];co<-v[v$classification_group=='Concordant / stable',];di<-v[v$classification_group=='Discordant / reclassified',]
 tt[[length(tt)+1]]<-data.frame(label=z$label,cohort=z$cohort,signature=z$target,partition=if(z$partition==0)'full-cohort' else as.character(z$partition),
  `n stable / reclassified`=paste(co$n,di$n,sep='/'),`stable HR [95% CI]`=paste(fmt(co$HR),ci(co$HR_low,co$HR_high)),
  `reclassified HR [95% CI]`=paste(fmt(di$HR),ci(di$HR_low,di$HR_high)),`interaction HR ratio [95% CI]`=paste(fmt(z$ratio_of_HRs),ci(z$low,z$high)),`interaction P`=fmt(z$p_interaction,4),check.names=FALSE)
}
lines<-c(lines,mdtable(do.call(rbind,tt)),'',
 '### A-exploratory 2 — recovered model-predicted 12-month benefit margins','',
 '**Completed with documented input-route deviation.** The original policy code discards arm-specific predictions after deriving each recommendation; the saved `phi`/`v` fields are augmented policy-value contributions, not individual predictions. The original deterministic estimator was reconstructed in R using the exact retained inputs and nominated training folds: arm-specific piecewise-exponential PFS models, intervals 0/3/6/9/12 months, early/late slopes for complete signature score and APM12, L2 penalty=1, original Newton convergence/line-search rules, and analytic 12-month integration. No model/endpoint/cut-off was changed. Neither augmented contributions nor a newly selected model was substituted.','',
 paste0('**PASS: ',nrow(mchecks),' recovery checks** in `results/policy_margin_reconstruction_checks.csv`: zero mismatches for every saved recommendation under both signatures/operators and both trials; exact retained clinical/matrix inputs; archived-fit events and zero-event intervals; original gradient criteria; nuisance hazard predictions reproduced within 1e-9 against the saved independent-check patients. The maximum nuisance-hazard discrepancy was ',format(max(mchecks$observed[grepl('nuisance hazard',mchecks$check)]),scientific=TRUE,digits=3),'. Per-fit diagnostics: `results/policy_margin_reconstruction_fits.csv`. These checks qualify deterministic prediction recovery, not new clinical validation.'),'',
 'Benefit margin = model-predicted experimental-arm 12-month RMST minus sunitinib RMST (months). Classification groups are nominated-partition median reclassified versus stable patients, not the policy-discordant subgroup. The table describes fixed prediction distributions; the locked plan does not specify a margin CI, so IQR is reported explicitly rather than relabelled as a confidence interval. Source: `results/predicted_benefit_margin_summary.csv`; all arm-specific predictions: `results/policy_model_predictions_reconstructed_LOCAL_ONLY.csv`. The historical `benefit_margin_NOT_RUN.csv` filename is retained only as an audit record and its rows now explicitly state completed by reconstruction; current availability is `benefit_margin_input_status.csv`.','',
 mdtable(data.frame(label=margin$label,cohort=margin$cohort,signature=margin$target,method=margin$method,group=margin$classification_group,n=margin$n,
 `mean margin, months`=fmt(margin$mean_benefit_margin),`median margin [IQR], months`=paste(fmt(margin$median_benefit_margin),ci(margin$q25_benefit_margin,margin$q75_benefit_margin)),
 `median absolute margin, months`=fmt(margin$median_absolute_margin),`n recommending experimental`=margin$n_recommend_experimental,check.names=FALSE)),'',
 '## A-sensitivity — all partitions, Angio16 and alternative fixed research boundaries','',
 'The full transition-by-arm PFS table includes all other four saved partitions, Angio16 under all five partitions, and both signatures’ full-cohort median descriptive comparisons (176 sensitivity rows in `results/transition_PFS_outcomes.csv`). IMmotion151 response is retained for all corresponding sensitivity designs (88 sensitivity rows in `results/transition_ORR_IMmotion151.csv`). All treatment HRs and interaction tests for these designs are printed above, including favourable and unfavourable results. Partitions reuse the same patients and are dependent; they are not five independent replications.','',
 'Median-threshold reclassification totals (every design; exact conditional 95% intervals) follow. Source: `results/median_reclassification.csv`; arm-specific rates remain in that source CSV.','')
z<-rc[rc$arm=='All randomised arms',]
lines<-c(lines,mdtable(data.frame(label=z$label,cohort=z$cohort,signature=z$target,partition=ifelse(z$partition==0,'full-cohort',z$partition),n=z$n,reclassified=z$reclassified,
 `per 100 [95% CI]`=paste(fmt(z$per100),ci(z$low,z$high)),unavailable=z$unavailable,check.names=FALSE)),'',
 'Each prespecified training boundary is evaluated as a binary research high/low split: lower and upper tertile (1/3, 2/3), and quartile boundaries (1/4, 1/2, 3/4). Both complete signatures and all five saved partitions are reported below; source: `results/threshold_sensitivity_reclassification.csv`. Training-fold thresholds, test sizes and training sizes are in `results/training_threshold_manifest.csv`. No outcomes enter threshold construction. The median column overlaps the main held-out median analysis by design.','')
tw<-list()
for(tr in names(inputs))for(tar in c('Angio36','Angio16'))for(re in 1:5){v<-th[th$cohort==tr&th$target==tar&th$partition==re,];v<-v[order(v$threshold_probability),]
 tw[[length(tw)+1]]<-data.frame(label='sensitivity',cohort=tr,signature=tar,partition=re,n=v$n[1],`1/4: n; per100 [CI]`=paste0(v$reclassified[1],'; ',fmt(v$per100[1]),' ',ci(v$low[1],v$high[1])),
 `1/3: n; per100 [CI]`=paste0(v$reclassified[2],'; ',fmt(v$per100[2]),' ',ci(v$low[2],v$high[2])),`1/2: n; per100 [CI]`=paste0(v$reclassified[3],'; ',fmt(v$per100[3]),' ',ci(v$low[3],v$high[3])),
 `2/3: n; per100 [CI]`=paste0(v$reclassified[4],'; ',fmt(v$per100[4]),' ',ci(v$low[4],v$high[4])),`3/4: n; per100 [CI]`=paste0(v$reclassified[5],'; ',fmt(v$per100[5]),' ',ci(v$low[5],v$high[5])),check.names=FALSE)
}
lines<-c(lines,mdtable(do.call(rbind,tw)),'',
 '## Methods, provenance and reproducibility','',
 '- Run order and commands: `scripts/README.md`; scoring and statistical code: `scripts/common.R`, `01_reproduce.R`, `02_analyse.R`, `02b_recover_policy_margins.R`, `03_report_and_validate.R`. No manuscript, retained analysis input/source code or locked analysis plan was edited. Exception: the failed initial graphics device overwrote the pre-existing incidental `/LOCAL_DATA_ROOT/Rplots.pdf`; previous PDF contents were not recoverable from the current file and have not been restored. See DEVIATIONS.md. Subsequent automatic graphics output is confined to this module.',
 '- Runtime: `/usr/local/bin/Rscript`, R 4.6.0; base R and survival only, no survminer and no package installation. `sessionInfo.txt` records versions. New resampling seed is 20261113; original saved partitions and references were reused.',
 '- Input matrices are complete Angio16/Angio36 lists from folder37. JAVELIN uses released log2 TPM without a second transform; IMmotion151 uses released log2(TPM+1) without a second transform. Final aggregate standardisation is learned in training only for held-out analyses. Complete gene membership was checked against final folder83 membership.',
 '- Input paths, resolved symlink targets, byte sizes and SHA-256 hashes: `results/input_manifest.csv`; final output hashes: `results/output_manifest.csv`. The locked plan’s current SHA-256 is recorded with the input manifest for immutable provenance.',
 '- Validation: `results/validation_checks.csv` checks counts, response-category conservation, completeness and independent KM integration parity. Explicit missing/not-estimable rows are retained in all output tables. Numeric CSVs use NA plus status fields; readable report tables print “not estimable”.',
 '- Resource evidence: `logs/02_resource_usage.log` and `logs/02b_resource_usage.log` record elapsed time and maximum resident memory for the serial analyses. One R process, one BLAS thread, no workers; all intermediates fit well below 50 MB, so no external scratch or project-local library was needed. The unavailable Cairo/X11 renderer was replaced with the native Quartz bitmap renderer; no statistical method changed (see DEVIATIONS.md).',
 '- Diagnostic plots are unpolished PNGs; `plots/primary_transition_group_counts.png` shows the complete primary transition-by-arm counts. Participant-level outputs ending LOCAL_ONLY and the verified-frame RDS remain governed by the original data conditions.',
 '- Claim-evidence gate: PASS for reproducibility, prespecified patient-unit estimands and complete reporting; exploratory for all clinical interpretation. Independent validation and individual-benefit/clinical-deployment claims are not established. Main reviewer risks are small discordant groups, incomplete 12-month support, cross-fitted outcome-derived policies, and multiple dependent sensitivities.','',
 '## Anything not run and why','',
 '1. **No specified available-input Module A component remains unrun.** The initially missing benefit predictions were recovered with the original estimator and passed all gates; this input-route deviation is documented in DEVIATIONS.md.',
 '2. **JAVELIN OS and objective response; IMmotion151 OS:** not available in retained source schemas, as declared in the locked plan; no substitutes used.',
 '3. **Non-estimable primary 12-month PFS/RMST and median bounds, and undefined bootstrap draws:** retained as non-estimable because of follow-up/event support, not compute failure. Exact affected groups and attempted draws are in the outcome/bootstrapping CSVs. No requested analysis was replaced with downsampling, a proxy or a lighter method.',
 '4. Policy contrasts and predicted-margin summaries use the nominated-partition saved policy definitions; other saved partitions supply the classification/outcome sensitivities. No additional learned policies were invented.','',
 '## What this does and does not show','',
 'Implementation changes research high/low assignments for 42/726 JAVELIN and 85/823 IMmotion151 patients on the nominated held-out partition, with explicit observed outcomes and response denominators by randomised arm.',
 'The policy-discordant groups are different model-defined groups (40 and 76 Angio36 patients); their observed PFS contrasts are imprecise and do not identify an individually preferable treatment.',
 'All prespecified signatures, saved partitions and training-quantile sensitivities are retained, including results that differ in direction or nominal significance.',
 'Small groups and insufficient follow-up produce genuine non-estimability; conditional bootstrap intervals do not establish unconditional 95% coverage.',
 'Discarded benefit predictions were recovered from the original deterministic estimator after numerical/decision reproduction gates; no augmented contribution is relabelled as an individual prediction.',
 'This post hoc secondary analysis supports an implementation-reporting conclusion, not prospective biomarker validation or a treatment recommendation.')
writeLines(lines,file.path(OUT,'RESULTS.md'))
cat('RESULTS.md written; all validation checks passed.\n')
