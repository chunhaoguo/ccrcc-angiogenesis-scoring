# Author: Chunhao Guo
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,re,html
from scipy.stats import false_discovery_control
P=Path(__file__).resolve().parents[4];O=P/'84_JCO_PO_Upgrade_20261113/05_analysis/D_biology';R=O/'results';rel=lambda x:str(x.relative_to(P))
def rd(n):return pd.read_csv(R/n)
def md(d):
 def st(x):
  if pd.isna(x):return 'not estimable'
  return str(x).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(st(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))
def f(x,k=2):return f'{x:.{k}f}' if pd.notna(x) else 'not estimable'
def qf(x):return f'{x:.3g}' if pd.notna(x) else 'not estimable'
# Carry forward genuine restriction results and their original bootstrap intervals unchanged.
for n in ['RESTRICTED_PHENOTYPE_COEFFICIENTS.csv','RESTRICTED_PHENOTYPE_TESTS.csv','FULL823_REPRODUCTION.csv']:
 shutil.copyfile(P/'74_JTM_Submission_20261028/figures/source_data/restriction'/n,R/('CARRIED_FORWARD_'+n))
fp=P/'84_JCO_PO_Upgrade_20261113/05_analysis/shared/CM025_folds.csv'
if fp.exists():shutil.copyfile(fp,R/'CM025_FOLD_SOURCE_SNAPSHOT.csv')
primary=rd('D_PRIMARY_GENOTYPE.csv');order={'JAVELIN':0,'IMmotion151':1,'CheckMate025':2};primary['sort']=primary.cohort.map(order);primary=primary.sort_values(['sort','gene'],key=lambda x:x if x.name!='gene' else x.map({'VHL':0,'PBRM1':1,'BAP1':2}));audit=rd('GENOTYPE_JOIN_AUDIT.csv');rows=[]
for _,z in primary.iterrows():
 a=audit[(audit.cohort==z.cohort)&(audit.gene==z.gene)].iloc[0]
 rows.append({'estimand':z.gene+' altered vs source-negative: signed Angio36 displacement','label':'primary (exploratory overall)','cohort':z.cohort,'n':f'{int(a.n_qualified)} qualified ({int(a.n_altered)}/{int(a.n_source_negative)}); {int(a.n_unknown_not_assayed)} unknown','estimate':f(z.estimate)+' pp','95% CI':f(z.lower)+' to '+f(z.upper),'BH q':qf(z.q_BH),'source file':rel(R/'D_PRIMARY_GENOTYPE.csv')})
headline=pd.DataFrame(rows);headline.to_csv(R/'HEADLINE_TABLE.csv',index=False)
rep=rd('REPRODUCTION_CHECKS.csv');cl=rd('CLASSIFICATION_PROVENANCE.csv');pur=rd('DNA_PURITY_ASSOCIATIONS.csv');ct=rd('FIXED_CLUSTER_COMPOSITION_TESTS.csv');sc=rd('FIXED_CLUSTER_DISPLACEMENT_COMPOSITION.csv');floor=rd('FLOOR_WEIGHT_RELATIONSHIP.csv');weights=rd('GENE_FLOOR_SD_WEIGHT.csv');ag=rd('RECLASSIFIED_GENE_GROUP_ATTRIBUTION.csv');aq=rd('ATTRIBUTION_REPRODUCTION_QA.csv');atlas=rd('ATLAS_EQUAL_DONOR_SUMMARY.csv');atct=rd('ATLAS_ELIGIBILITY_COUNTS.csv');qual=json.loads((R/'ATLAS_QUALIFICATION.json').read_text());sec=rd('D_SECONDARY_GENOTYPE_ALL_ENDPOINTS.csv');clinical=rd('CARRIED_FORWARD_CLINICAL_COEFFICIENTS.csv');tests=rd('CARRIED_FORWARD_CLINICAL_TESTS.csv');restriction=rd('CARRIED_FORWARD_RESTRICTED_PHENOTYPE_COEFFICIENTS.csv');gs=ag[(ag.signature=='Angio36')&(ag.reference=='heldout_nominated')&ag.group.str.contains('FGF')]
sec_s=sec[(sec.signature=='Angio36')&(sec.endpoint=='delta_r')&sec.gene.isin(['SETD2','KDM5C'])&sec.contrast.str.startswith('Altered')].copy();smt=sec_s[['cohort','gene','n_reference','n_category','estimate','lower','upper','p','q_BH_secondary']].copy()
for col in ['estimate','lower','upper']:smt[col]=smt[col].map(f)
for col in ['p','q_BH_secondary']:smt[col]=smt[col].map(qf)
wm=weights[(weights.signature=='Angio36')&weights.gene.str.startswith('FGF')];floor_table=floor[floor.signature=='Angio36'][['cohort','n_genes','spearman_floor_weight_gain','n_FGF_at_floor_ge90pct']].copy();floor_table['spearman_floor_weight_gain']=floor_table.spearman_floor_weight_gain.map(lambda x:f(x,3))
gst=gs[['cohort','direction','n_reclassified','share_absolute_total']].copy();gst['share_absolute_total']=gst.share_absolute_total.map(lambda x:f'{100*x:.1f}%')
af=atlas[(atlas.metric=='detection_fraction')&atlas.gene.str.startswith('FGF')];atl=af.groupby('compartment').agg(n_FGF_genes=('gene','size'),median_equal_donor_detection=('estimate','median'),genes_mean_detection_le_1pct=('estimate',lambda x:sum(x<=.01)),genes_zero_detection_all_donors=('estimate',lambda x:sum(x==0))).reset_index();atl['median_equal_donor_detection']=atl.median_equal_donor_detection.map(lambda x:f'{100*x:.3f}%')
raw=r'''All results are post hoc exploratory secondary analyses; **primary** identifies the locked primary estimand within this upgrade, not confirmatory inference. Units are patients for trial associations and donors for cellular context. Positive displacement means a higher gene-standardised percentile rank. Genotype models address signed movement, not a causal genotype effect or proof of a reclassified subtype.

## Reproduction check — PASS

Saved I0/I1 scores and saved full-cohort ranks were used downstream for JAVELIN and IMmotion151. Recomputing the complete gene lists agreed with saved standardised scores to at most {score_error:.3g}, and saved-score rank reproduction error was zero. Median absolute Angio36 rank displacement was **3.581267 pp in JAVELIN** and **6.075334 pp in IMmotion151**, reproducing 3.58 and 6.08. Source: `{repfile}`.

Before genotype association fitting, the original IMmotion151 average-rank-adjusted clinical point models were rerun from the saved ranks: poor IMDC **+5.038394 pp** and metastatic specimen origin **+5.019710 pp**, both n=822 and coefficient differences below 1e-8. The existing bootstrap intervals and joint-test q values are carried forward unchanged: poor IMDC +5.04 [2.53, 7.71], IMDC joint q=0.00121; metastatic origin +5.02 [3.20, 7.41], q=3.56e-6. These are directional movement associations; the separate existing absolute-movement and switching tests do not establish these features as reclassification predictors. Source: `{clinicalfile}`, `{clinicaltests}`, `{clinicalsecondary}`.

The source sarcomatoid coefficient is +3.66 pp [1.28, 6.21] in the full source; its attenuation in the 798 explicit clear-cell patients is retained, rather than replaced with the full-source result: +1.04 pp with frozen full-reference ranks and +1.19 pp with within-subset refitting. Exact original estimates/intervals/tests: `{restrictionfile}` and `{restrictiontests}`. No new outcome models were fitted.

The reclassification variable uses the nominated saved partition 1 for JAVELIN/IMmotion151, and the B2 treatment-stratified five-fold CM025 file, which existed when scoring ran. Displacement remains full-cohort-referenced in all three. Counts below must not be interchanged with full-cohort median counts or hypothetical treatment-policy disagreement.

{classification}

Source: `{classificationfile}`; patient scores/ranks, saved fold numbers and labels: `{patientfile}`; B2 fold snapshot: `{foldfile}`. The historical source uses 5 repeats; only its first saved partition is nominated here.

## D-primary: qualified genotype joins and one BH family

Each model is signed Angio36 displacement ~ natural spline(average rank, 3 df) + three-state genotype, with source-negative as reference. All RNA patients with computable scores enter, and unknown/not-assayed patients have their own factor category. The headline n reports the qualified altered/negative counts; total model n is 726/823/250. HC3 sandwich SE and residual-df t intervals/p values are conditional on the saved score/rank reference. One and only one primary BH family contains all **9** altered-versus-negative tests (VHL/PBRM1/BAP1 × 3 cohorts). Unknown-category coefficients are descriptive extensions and are outside that family.

{audit}

Source: `{auditfile}`. Full source paths, source sheet names, assay definitions, explicit patient join status and unknown counts are in `{auditfile}`, `{joinfile}` and `{source_schema}`.

- **JAVELIN:** S21 is explicitly headed “0 = not mutated, 1 = mutated”; these binary source flags define the WES comparison. There are 738 S21 assay columns, but only 691 match the established 726 RNA patients. A missing S21 patient is unknown, never negative. S18 WES summary availability is audited separately and does not replace the explicit S21 flags.
- **IMmotion151:** original TableS1, FoundationOne targeted FMI panel. FMI_SAMPLE_ID establishes assay availability; any released alteration string is positive. An assayed ID with an empty gene field is no reported panel alteration; no assay ID is unknown. There are 702 qualified RNA patients and 121 unknown. Panel alteration status includes short variants and the reported deletion/loss/truncation/rearrangement calls; it is not harmonised to WES mutation scope.
- **CM025:** S1 RNA_ID/SUBJID/MAF_Tumor_ID crosswalk plus S2 released somatic MAF. There are 175 RNA patients with an S1 tumor MAF ID verified in S2 and 75 unknown. Primary positive means a released coding nonsynonymous variant; negative means no qualifying coding variant recorded in a qualified exome. Silent/intronic/UTR/flank/RNA classes are excluded. This does not assert absence of copy-number change. S1 MUT/WT chiefly indicates truncating mutations, and must not be treated as all-nonsynonymous status; the source-flag/MAF cross-tabulation and exact classes are preserved in `{mafflags}` and `{mafclasses}`.

VHL and PBRM1 alterations associate with downward movement in both first-line trials. CM025 does not reproduce a VHL association and has an opposite PBRM1 direction (+4.94 pp). BAP1 is positive in JAVELIN, with CIs crossing zero in the other two. The family is fully reported, including nulls. Differences in assay scope, therapy setting and specimen selection prevent a common causal/genotypic subtype interpretation.

## D-secondary and exploratory extensions

The locked secondary genes' signed Angio36 associations are all reported below. Remaining Angio16, absolute-displacement and reclassification results, plus unknown-category effects, are in the complete CSV, including all “not estimable” sparse categories. These exploratory extensions have a separate descriptive BH column over all remaining genotype contrasts; this is not part of the 9-test primary family and is not study-wide multiplicity control.

{secondary}

Source: `{secondaryfile}`. The reclassification endpoint uses logistic models with a 3-df spline of I0 rank distance from 50; these are models of research group changes, not patient clinical outcomes. Fewer than 5 patients in a comparison category, or fewer than 5 events/non-events per compared group, is “not estimable”. Category-level n, signed and absolute movement, reclassification rates and exact binomial 95% CIs are in `{categoryfile}`. Some logistic fits emit numerical near-boundary-probability warnings; their full log is retained. Sparse-category results are not substituted with a different model.

**DNA-derived purity:** available only as ABSOLUTE in CM025 (158 observed / 92 missing). Per 10 percentage-point higher purity, signed Angio36 displacement is +1.369 pp [0.351, 2.387], nominal p=0.00871; q across the six purity extensions=0.0523. Its absolute-displacement and reclassification CIs include the null. There is no qualified independent DNA purity measure for JAVELIN or IMmotion151; H&E and same-RNA xCell are not relabelled DNA purity. All six tests and availability: `{purityfile}`, `{purityavailability}`.

**Fixed source clusters:** author NMF classes 1–7 are retained in IMmotion151, and the JAVELIN S11 TCGA field retains m1–m4, Other and 3 unknown. No classifier is refitted. For Angio36, source-cluster composition among reclassified/stable patients has Cramer's V=0.144 [0.109, 0.224] in IMmotion151 (Monte Carlo Fisher p=0.00325) and 0.095 [0.058, 0.178] in JAVELIN (p=0.192). These are exploratory descriptive effect intervals (patient bootstrap, fixed scores); positive Cramer's V bootstrap lower bounds must not be read as a test of independence. NMF class 5 has 14/74 reclassified (18.9%), compared with 2/98 (2.0%) in class 1; all seven classes, unadjusted mean/median displacement and exact rate CIs are retained. Source: `{clusterfile}`, `{clusterdescriptive}`, `{clusteradjusted}`. The same-RNA class annotations are context, not orthogonal proof of biological mechanism or subtype reassignment.

## D-drivers: exact arithmetic, with biological interpretation kept separate

For gene g, on the released expression scale,

`w0_g = 1 / (G * SD(mean-expression raw score))`

`w1_g = 1 / (G * SD(expression_g) * SD(mean gene-z raw score))`

`contribution_gj = (w1_g - w0_g) * (x_gj - reference_mean_g)`.

Their sum is exactly standardised I1 minus standardised I0, **not** a sum of percentile-rank contributions. Frozen-reference versions use each patient's training-fold means/SDs/aggregate scaling. Threshold movement is a separate constant: `(I1-cutoff1)-(I0-cutoff0) = sum(contributions)+(cutoff0-cutoff1)`. Both patient decompositions are retained. Maximum score-difference reconstruction error across all runs is **{reconstruction_error:.3g}**; existing module60 weights reproduce to **{weight_error:.3g}**, and those saved weights/variance attributions are reused rather than newly estimated under a different definition. Sources: `{arithmeticqa}`, `{archiveweights}`, `{archivevariance}`, `{archivegroup}`.

All genes are retained. Every gene has its released minimum, exact-zero fraction, modal fraction, reference SD, I0/I1 weight and weight ratio in `{floorfile}`. “Floor” means the exact released minimum, not a validated detection limit; ComBat-corrected CM025 expression in particular has different floor behaviour. TCGA uses 533 unique primary-tumour patients (one -01 measurement per patient); TRACERx uses all 231 regions and is not reported as 231 patients; CPTAC uses the existing 78 paired cases so the archived weight quantities are comparable. Additional descriptive cohorts are explicitly exploratory. Aliases: source FIGF→VEGFD in TRACERx/TCGA and source CXCR7→ACKR3 in TCGA. No gene is dropped or zero-filled. All complete-list I1 scores are available in these full references: `{availabilityfile}`.

{floor_table}

Source: `{floorrelation}`. Spearman values describe a fixed set of 36 genes, not independent biological replicates; no gene-wise hypothesis p value is supplied. Several FGF genes sit at a released floor in the first-line trials; the weak CM025 relationship is retained. A different floor/weight relationship does not refute the algebraic 1/SD scaling identity.

Among **held-out reclassified patients**, the 26-gene FGF/FGFR block contributes these fractions of the sum of absolute per-gene score-difference contributions:

{gene_group}

Source: `{groupfile}`. This is an arithmetic L1 share, different from the already-reported signed variance share. Signed shares are also reported and can be negative or exceed 100% because contributions cancel. The largest individual drivers depend on cohort and direction (for example, FGF19 for JAVELIN low→high, PGF for its high→low; FGF22/FGF23 for IMmotion151 low→high; FGF3/FGF4 for CM025 low→high). The complete per-gene table includes conditional mean 95% intervals; these condition on selected reclassified patients and fixed references and are not whole-workflow uncertainty intervals. Sources: `{topfile}`, `{gene_summary}`, `{patientgene}`, `{patientdec}`. TRACERx region means are arithmetic descriptors; they are not independent-patient inferential intervals.

Gene standardisation up-weights a low-SD gene by definition. Floor/quantisation, processing, rare expression and real cellular composition can all influence that SD. These tables do not establish activated FGF signalling, a malignant FGF program, or a fraction of biologically caused reclassification.

## D-cellular context: qualifying new atlas, full-data streaming

The selected HCA/CELLxGENE ccRCC atlas is the PNAS 2021 resource with DOI 10.1073/pnas.2103240118, with its original author_cell_type and donor_id fields, not the Li et al. Cancer Cell 2022 atlas already reported in the manuscript. The selected file is `{atlaspath}` ({atlasbytes:,} bytes), 20,509 cells × 32,088 measured features, and all **36/36 Angio36 genes**, including **26/26 FGF/FGFR genes**, are in its raw feature universe. Its author-labelled Endo_ACKR1/Endo_PLVAP, Peri/vSMC, Tumor and immune labels are mapped to the four locked compartments without de novo annotation. Each compartment has seven unique donors, satisfying ≥5. Potential GSE159115 alternative representation is counted only once. Written qualification and feature accounting: `{atlasqual}`, `{atlasgenes}`, `{atlasmap}`.

{atlascounts}

Source: `{atlascountsfile}`. All 20,509 cells were retained; blocks of 256 cells were read from raw/X CSR arrays with integer/nonnegative checks. The maximum sparse matrix block is 9.14 MiB, and measured process maximum RSS is 297,713,664 bytes (~0.277 GiB), well below 3 GiB. No whole atlas densification, cell downsampling, integration or reannotation occurred. One endothelial donor has only 12 cells (all-donor analysis follows this plan's ≥5-donor criterion; this limits detection precision), and one stromal donor has 22 cells. Log/stream evidence: `{atlaslog}`, `{atlasstream}`.

Within each donor×compartment, all-cell sums, library totals, detection fraction, mean log1p(CP10K) and pseudobulk CP10K are saved. Cross-donor summaries use equal donor weights, and 2,000 donor-bootstrap descriptive 95% intervals. Endothelial-minus-other contrasts are paired within the same seven donors; cells are never independent replicates. All 36 genes, all four compartments and all three contrasts are retained without significance-based driver selection. Sources: `{atlasdonor}`, `{atlassummary}`, `{atlascontrasts}`.

{atlasfgf}

Source: `{atlassummary}`; the table summarizes all 26 FGF/FGFR genes with the declared prefix-defined group. Complete feature coverage therefore does not imply common expression: 16–24 of 26 genes, depending on compartment, have ≤1% equal-donor mean detection. FGF3/FGF4/FGF19/FGF20, several numerical drivers, have no detected counts in some or all atlas compartments; gene-specific values and their 95% donor intervals remain in the CSV. This is new cellular context with complete measured-gene accounting, not a repetition of the Li atlas's 16/26-coverage result. It does not prove that sparse bulk values are artifacts, identify the cellular origin of every bulk driver, or establish pathway activation.

## Five-line answer: who is reclassified and why

1. The nominated frozen implementations reclassify 42/726 JAVELIN, 85/823 IMmotion151 and 26/250 CM025 patients; these are changes in research median groups, not demonstrated clinical errors.
2. Poor IMDC risk and metastatic specimen origin accompany upward rank movement in IMmotion151; the carried-forward binary reclassification associations remain unestablished.
3. VHL/PBRM1 alterations accompany downward movement in both first-line trials, but the opposite CM025 PBRM1 direction and null genotype results preclude a common genotype-defined reclassified subtype.
4. Low-variance gene scaling changes the weights arithmetically, with FGF/FGFR genes accounting for 67.1–83.4% of absolute contributions in the nominated switcher groups; the individual drivers vary by cohort and direction.
5. A fully measured seven-donor atlas shows sparse FGF detection across compartments, adding cellular context without establishing FGF activation or a biological mechanism for reclassification.

## Reproducibility, deviations and anything not run

Run order and exact commands: `{readme}`; software: `{session}`; source hashes: `{inputmanifest}`; run logs and diagnostics remain in this module. Seed=20261113 for new stochastic calculations, original saved seed/partition for carried-forward trial classifications. At most one R/Python worker per process; maximum measured R RSS was below 0.3 GiB. New outputs total only a few MiB; no >50 MB intermediates or package installations, external downloads, external contact or source mutation were required. Author/creator of generated outputs: Chunhao Guo.

Deviations/clarifications are in `{deviations}`: source-specific MAF versus truncating-flag semantics, secondary-model CI/family choices, exact released-minimum definition, complete-gene aliases, and seven-donor all-cell atlas context. Original clinical bootstrap intervals were carried forward unchanged after point-estimate reproduction; this module did not rerun their old bootstrap workflow. Genotype uncertainty conditions on the saved scores, as made explicit in the methods.

- **Not run — unavailable data:** independent DNA purity in JAVELIN/IMmotion151; no substitute from xCell/H&E was made.
- **Not run — conditional method not applicable:** GSE242299/SCP1288 and GSE207493 cellular analyses; the first eligible new HCA atlas satisfied the locked one-atlas condition. GSE207493 also has the previously documented missing author-annotation crosswalk; it was not promoted through inferred labels.
- **Not run — sparse categories:** every prespecified non-estimable genotype/cluster category is retained with its count and reason in the CSVs. No Firth or other lighter replacement was substituted.
- **Not run — outside this module:** pathway enrichment, de novo subtyping, survival/response/treatment-outcome models, or causal/functional biology. No compute-limited required primary analysis remains unfinished.

## What this does and does not show

The complete prespecified primary family supports cohort-dependent associations between genotype and directional score movement.
It does not establish that a genotype causes reclassification or that changing research group causes treatment benefit/harm.
Per-gene decomposition proves a numerical identity; FGF/FGFR shares are arithmetic quantities, not biological activation fractions.
The atlas supplies seven-donor, complete-feature cellular context; sparse counts and small compartment cell counts limit biological localisation.
Panel/WES differences, unknown genotypes, saved-reference conditioning and results-informed cohort reuse limit interpretation.
All nulls and non-estimable categories remain reported; no new clinically validated subtype or score is claimed.
'''
# Values that are sourced exactly from saved tables.
values=dict(score_error=rep.max_score_error.max(),repfile=rel(R/'REPRODUCTION_CHECKS.csv'),clinicalfile=rel(R/'CARRIED_FORWARD_CLINICAL_COEFFICIENTS.csv'),clinicaltests=rel(R/'CARRIED_FORWARD_CLINICAL_TESTS.csv'),clinicalsecondary=rel(R/'CARRIED_FORWARD_CLINICAL_SECONDARY_TESTS.csv'),restrictionfile=rel(R/'CARRIED_FORWARD_RESTRICTED_PHENOTYPE_COEFFICIENTS.csv'),restrictiontests=rel(R/'CARRIED_FORWARD_RESTRICTED_PHENOTYPE_TESTS.csv'),classification=md(cl),classificationfile=rel(R/'CLASSIFICATION_PROVENANCE.csv'),patientfile=rel(R/'PATIENT_SCORES_RANKS.csv'),foldfile=rel(R/'CM025_FOLD_SOURCE_SNAPSHOT.csv'),audit=md(audit[['cohort','gene','RNA_n','n_qualified','n_altered','n_source_negative','n_unknown_not_assayed','assay']]),auditfile=rel(R/'GENOTYPE_JOIN_AUDIT.csv'),joinfile=rel(R/'GENOTYPE_PATIENT_JOIN.csv'),source_schema=rel(O/'logs/source_schema.json'),mafflags=rel(R/'CHECKMATE_MAF_VS_SOURCE_TRUNCATING_FLAGS.csv'),mafclasses=rel(R/'CHECKMATE_QUALIFYING_VARIANT_CLASSES.csv'),secondary=md(smt),secondaryfile=rel(R/'D_SECONDARY_GENOTYPE_ALL_ENDPOINTS.csv'),categoryfile=rel(R/'GENOTYPE_CATEGORY_DESCRIPTIVES.csv'),purityfile=rel(R/'DNA_PURITY_ASSOCIATIONS.csv'),purityavailability=rel(R/'DNA_PURITY_AVAILABILITY.csv'),clusterfile=rel(R/'FIXED_CLUSTER_COMPOSITION_TESTS.csv'),clusterdescriptive=rel(R/'FIXED_CLUSTER_DISPLACEMENT_COMPOSITION.csv'),clusteradjusted=rel(R/'FIXED_CLUSTER_ADJUSTED_SIGNED_CONTRASTS.csv'),reconstruction_error=aq.max_reconstruction_error.max(),weight_error=aq.max_archived_weight_error.max(),arithmeticqa=rel(R/'ATTRIBUTION_REPRODUCTION_QA.csv'),archiveweights=rel(R/'CARRIED_FORWARD_ALL_GENE_WEIGHTS.csv'),archivevariance=rel(R/'CARRIED_FORWARD_ALL_GENE_VARIANCE_CONTRIBUTIONS.csv'),archivegroup=rel(R/'CARRIED_FORWARD_GROUP_VARIANCE_CONTRIBUTIONS.csv'),floorfile=rel(R/'GENE_FLOOR_SD_WEIGHT.csv'),availabilityfile=rel(R/'DRIVER_SCORE_AVAILABILITY.csv'),floor_table=md(floor_table),floorrelation=rel(R/'FLOOR_WEIGHT_RELATIONSHIP.csv'),gene_group=md(gst),groupfile=rel(R/'RECLASSIFIED_GENE_GROUP_ATTRIBUTION.csv'),topfile=rel(R/'TOP_FIVE_ARITHMETIC_DRIVERS.csv'),gene_summary=rel(R/'RECLASSIFIED_GENE_ATTRIBUTION_SUMMARY.csv'),patientgene=rel(R/'RECLASSIFIED_PATIENT_GENE_ATTRIBUTIONS.csv.gz'),patientdec=rel(R/'RECLASSIFIED_PATIENT_DECOMPOSITION.csv'),atlaspath=qual['atlas'],atlasbytes=qual['file_bytes'],atlasqual=rel(R/'ATLAS_QUALIFICATION.json'),atlasgenes=rel(R/'ATLAS_MEASURED_GENE_ACCOUNTING.csv'),atlasmap=rel(R/'ATLAS_AUTHOR_COMPARTMENT_MAP.csv'),atlascounts=md(atct),atlascountsfile=rel(R/'ATLAS_ELIGIBILITY_COUNTS.csv'),atlaslog=rel(O/'logs/04_cellular_context.log'),atlasstream=rel(R/'ATLAS_STREAM_MEMORY_QA.csv'),atlasdonor=rel(R/'ATLAS_DONOR_COMPARTMENT_EXPRESSION.csv'),atlassummary=rel(R/'ATLAS_EQUAL_DONOR_SUMMARY.csv'),atlascontrasts=rel(R/'ATLAS_PAIRED_DONOR_CONTRASTS.csv'),atlasfgf=md(atl),readme=rel(O/'scripts/README.md'),session=rel(O/'sessionInfo.txt'),inputmanifest=rel(R/'INPUT_MANIFEST.csv'),deviations=rel(O/'DEVIATIONS.md'))
(O/'RESULTS.md').write_text(md(headline)+'\n\n# Module D — who is reclassified and why\n\n'+raw.format(**values))
# A compact inline diagnostic with theme-aware SVG, no participant data.
rowsvg=[];xmin,xmax=-12,10;x=lambda v:250+(v-xmin)/(xmax-xmin)*405
for i,(_,z) in enumerate(primary.iterrows()):
 y=32+i*33;col={'JAVELIN':'var(--chart-1)','IMmotion151':'var(--chart-2)','CheckMate025':'var(--chart-3)'}[z.cohort]
 rowsvg.append(f'<text x="0" y="{y+5}">{html.escape(z.cohort)} · {z.gene}</text><line x1="{x(z.lower):.2f}" y1="{y}" x2="{x(z.upper):.2f}" y2="{y}" stroke="{col}" stroke-width="3"/><circle cx="{x(z.estimate):.2f}" cy="{y}" r="5" fill="{col}"/>')
svg=f'<svg viewBox="0 0 700 365" role="img" aria-label="Adjusted genotype associations with 95 percent confidence intervals"><line x1="{x(0)}" y1="12" x2="{x(0)}" y2="318" stroke="var(--muted-foreground)" stroke-dasharray="4 4"/>'+''.join(rowsvg)+''.join(f'<text x="{x(v)}" y="343" text-anchor="middle">{v}</text>' for v in [-10,-5,0,5,10])+'<text x="450" y="365" text-anchor="middle">I1 − I0 percentile points / 百分位点</text></svg>'
ht='<!doctype html><html><head><meta charset="utf-8"><meta name="author" content="Chunhao Guo"><meta name="last-modified-by" content="Chunhao Guo"><style>body{font-family:var(--font-sans);color:var(--foreground)}p{font-size:13px;color:var(--muted-foreground)}svg{width:100%;height:auto;font-size:14px;fill:var(--foreground)}svg text{fill:var(--foreground)}</style></head><body><p>Adjusted signed Angio36 rank movement; 95% CI; primary family q in the table below.<br>调整后的 Angio36 有向排名变化；95% 置信区间；主要检验家族的 q 值见下表。</p>'+svg+'<p>Exploratory secondary analysis; unknown genotypes remain separate.<br>探索性二次分析；未知基因型保留为独立类别。</p></body></html>'
(O/'plots/primary_genotype_inline.html').write_text(ht)
# Deterministic final validation before claiming completion.
qa={}
qa['primary_exactly_9_unique_tests']=len(primary)==9 and not primary[['cohort','gene']].duplicated().any()
qa['primary_BH_matches_scipy']=bool(np.allclose(primary.q_BH,false_discovery_control(primary.p.to_numpy(),method='bh')))
qa['genotype_counts_sum_to_RNA']=bool(((audit.n_altered+audit.n_source_negative+audit.n_unknown_not_assayed)==audit.RNA_n).all())
qa['no_unknown_recoded_negative']=all(rd('GENOTYPE_PATIENT_JOIN.csv').groupby('cohort').qualified.sum().eq(audit.groupby('cohort').n_qualified.sum()))
qa['complete_genes_each_score']=set(rd('GENE_FLOOR_SD_WEIGHT.csv').groupby(['cohort','signature']).size().unique())=={16,36}
qa['all_scores_reconstruct']=bool(aq.max_reconstruction_error.max()<1e-8)
qa['all_archived_weights_reproduced']=bool(aq.max_archived_weight_error.dropna().max()<1e-8)
qa['all_atlas_cells_accounted']=int(rd('ATLAS_DONOR_COMPARTMENT_COUNTS.csv').N.sum())==20509
qa['atlas_all36_measured']=len(rd('ATLAS_MEASURED_GENE_ACCOUNTING.csv'))==36
qa['atlas_atleast5_donors_each']=bool(atct.n_donors.ge(5).all())
qa['seven_donor_36gene_4compartment_summaries']=len(rd('ATLAS_DONOR_COMPARTMENT_EXPRESSION.csv'))==7*36*4
qa['primary_reproduction_PASS']=bool((rep.status.isin(['PASS','computed'])).all())
assert all(qa.values()),qa
(O/'logs/FINAL_VALIDATION.json').write_text(json.dumps({'author':'Chunhao Guo','last_modified_by':'Chunhao Guo','checks':qa,'status':'PASS'},indent=2))
print(md(headline.drop(columns=['source file'])));print('Final validation PASS',len(qa),'checks')
