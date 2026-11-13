#!/usr/bin/env python3
# Author / creator / last modified by: Chunhao Guo
# post hoc follow-up (specified after modules B2 and D)
from pathlib import Path
import csv, hashlib, html, json, math
import numpy as np
import pandas as pd

O=Path(__file__).resolve().parents[1]
A=O.parent
ROOT=A.parents[1]
LABEL='post hoc follow-up (specified after modules B2 and D)'
R=O/'results'
read=lambda p: pd.read_csv(p,keep_default_na=True)
f1=read(R/'F1_GENOTYPE_ALL.csv');sens=read(R/'F1_SPECIMEN_SENSITIVITY_ALL.csv')
f2=read(R/'F2_PROTEIN_CARRIED_FORWARD.csv');f3=read(R/'F3_GROUP_ANCHORS_ALL.csv')
comb=read(R/'F1_F2_COMBINED_ANCHORS.csv');groups=read(R/'F3_FIXED_GENE_GROUPS.csv')
checks=[]
def check(name, passed, detail=''):
    checks.append(dict(check=name,pass_check=bool(passed),detail=str(detail),analysis_label=LABEL))
def wr(df,name):
    df=df.copy();df['analysis_label']=LABEL;df.to_csv(R/name,index=False,na_rep='NA')
check('F1 complete planned rows',len(f1)==126)
check('F1 sensitivity complete planned rows',len(sens)==126)
check('F2 complete planned rows',len(f2)==14)
check('F3 complete planned rows',len(f3)==56)
check('all F1 keys unique',not f1.duplicated(['cohort','signature','anchor','operator']).any())
check('all F3 keys unique',not f3.duplicated(['cohort','anchor','gene_group','operator']).any())
check('required reproduction checks',read(R/'REPRODUCTION_CHECKS.csv').pass_check.all() if 'pass_check' in read(R/'REPRODUCTION_CHECKS.csv') else read(R/'REPRODUCTION_CHECKS.csv')['pass'].all())
orig=read(A/'B2_independent_cohorts/results/CPTAC_protein_anchor.csv')
check('F2 every source value carried unchanged',orig.equals(f2[orig.columns]))
check('F2 no fabricated p/q',comb.loc[comb.module=='F2',['p_difference','q_BH']].isna().all().all())
gt=read(A/'D_biology/results/GENOTYPE_PATIENT_JOIN.csv')
scores=[]
for co in ['JAVELIN','IMmotion151']:
    s=read(A/f'B1_operator_panel_trials/results/{co}_scores_local_only.csv')
    s['cohort']=co;scores.append(s[['cohort','ID','signature','operator','score']])
s=read(A/'B2_independent_cohorts/results/CM025_scores.csv')
s=s.melt(id_vars=['patient','signature'],value_vars=[f'O{x}' for x in range(1,8)],var_name='operator',value_name='score').rename(columns={'patient':'ID'})
s['cohort']='CheckMate025';scores.append(s)
scores=pd.concat(scores,ignore_index=True)
def cohen_d(score,altered):
    x=score[altered];y=score[~altered]
    pooled=np.sqrt(((len(x)-1)*np.var(x,ddof=1)+(len(y)-1)*np.var(y,ddof=1))/(len(x)+len(y)-2))
    return (np.mean(x)-np.mean(y))/pooled
errors=[]
for r in f1.itertuples():
    ss=scores[(scores.cohort==r.cohort)&(scores.signature==r.signature)&(scores.operator==r.operator)]
    gg=gt[(gt.cohort==r.cohort)&(gt.gene==r.anchor)]
    d=ss.merge(gg,on=['cohort','ID'],validate='one_to_one')
    d=d[d.status!='Unknown_not_assayed']
    check(f'F1 counts {r.cohort}/{r.signature}/{r.anchor}/{r.operator}',len(d)==r.n and (d.status=='Altered').sum()==r.n_altered)
    if math.isfinite(r.estimate):
        est=cohen_d(d.score.to_numpy(),(d.status=='Altered').to_numpy())
        errors.append(abs(est-r.estimate))
check('independent numpy F1 pooled-SD point estimates',max(errors)<1e-12,max(errors))
gs=read(R/'F3_GROUP_SCORES_LOCAL_ONLY.csv');errors=[]
for r in f3[f3.cohort!='CPTAC'].itertuples():
    d=gs[(gs.cohort==r.cohort)&(gs.gene_group==r.gene_group)].merge(gt[(gt.cohort==r.cohort)&(gt.gene==r.anchor)],on=['cohort','ID'],validate='one_to_one')
    d=d[d.status!='Unknown_not_assayed']
    est=cohen_d(d[r.operator].to_numpy(),(d.status=='Altered').to_numpy());errors.append(abs(est-r.estimate))
check('independent numpy F3 genotype point estimates',max(errors)<1e-12,max(errors))
anchor=read(A/'B2_independent_cohorts/results/CPTAC_protein_proxy.csv')
errors=[]
for r in f3[f3.cohort=='CPTAC'].itertuples():
    d=gs[(gs.cohort=='CPTAC')&(gs.gene_group==r.gene_group)].merge(anchor,left_on='ID',right_on='case',validate='one_to_one')
    est=d[r.operator].rank().corr(d.protein_anchor.rank());errors.append(abs(est-r.estimate))
check('independent pandas F3 protein point estimates',max(errors)<1e-12,max(errors))
for name,df in [('F1',f1),('F1sensitivity',sens),('F3',f3)]:
    x=df[df.estimate.notna()];check(name+' every estimable result has 2000 bootstrap draws',(x.valid_bootstrap==2000).all())
    check(name+' all CI endpoints ordered',(x.lower<=x.upper).all() and (x.difference_lower<=x.difference_upper).all())
    x=df[df.operator=='O1'];x=x[x.estimate.notna()]
    check(name+' all self differences/intervals exactly zero',np.all(x[['difference_vs_O1','difference_lower','difference_upper']].to_numpy()==0))
    check(name+' q values in bounds',df.q_BH.dropna().between(0,1).all())
check('F1 O3 failure in Angio36 all trials',f1[(f1.signature=='Angio36')&(f1.operator=='O3')].estimate.isna().all())
check('JAVELIN specimen sensitivity explicitly not run',sens[sens.cohort=='JAVELIN'].status.str.startswith('not run').all())
check('all source genotype labels reused',set(gt.status)=={'Altered','No_reported_alteration','Unknown_not_assayed'})
check('all fixed groups reported',set(f3.gene_group)==set(groups.group) and len(groups)==46)
for p in R.glob('*.csv'):
    d=read(p);check('required label '+p.name,'analysis_label' in d and (d.analysis_label==LABEL).all())
wr(pd.DataFrame(checks),'VALIDATION_CHECKS.csv')
assert all(x['pass_check'] for x in checks),[x for x in checks if not x['pass_check']]

def fmt(x,d=3):
    return 'NA' if pd.isna(x) else f'{x:.{d}f}'
def ci(row,lo='lower',hi='upper'):
    return 'NA' if pd.isna(row[lo]) else f'[{fmt(row[lo])}, {fmt(row[hi])}]'
def tab(df,extra=False):
    heads=['anchor','cohort','signature','operator','estimate','95% CI','difference vs O1 [95% CI]; BH q','n','source CSV']
    if extra:heads.insert(3,'gene group')
    lines=['| '+' | '.join(heads)+' |','| '+' | '.join(['---']*len(heads))+' |']
    for _,r in df.iterrows():
        src='F3_GROUP_ANCHORS_ALL.csv' if r.module=='F3' else 'F1_F2_COMBINED_ANCHORS.csv'
        cells=[r.anchor,r.cohort,r.signature,r.operator,fmt(r.estimate),ci(r),f'{fmt(r.difference_vs_O1)} {ci(r,"difference_lower","difference_upper")}; q={fmt(r.q_BH)}',str(int(r.n)),f'`results/{src}`']
        if extra:cells.insert(3,r.gene_group)
        lines.append('| '+' | '.join(map(str,cells))+' |')
    return '\n'.join(lines)

headline=comb[(comb.signature=='Angio36')&comb.operator.isin(['O1','O2'])].copy()
wr(headline,'HEADLINE_TABLE.csv')
count=read(R/'GENOTYPE_COUNTS.csv')
cnt='| cohort | gene | altered | source-negative | unknown excluded |\n|---|---|---|---|---|\n'+'\n'.join(f'| {r.cohort} | {r.gene} | {r.n_altered} | {r.n_source_negative} | {r.n_unknown} |' for r in count.itertuples())
f3compact=f3.copy()
statements=[
 'Scoring operators track the genotype and protein anchors with different strengths and sometimes opposite differences from O1.',
 'For Angio36, O2 reduces the VHL and PBRM1 standardised contrasts in both first-line trials; the CheckMate 025 PBRM1 difference is positive (+0.151 [0.058, 0.249], BH q=0.052).',
 'In CPTAC, Angio36 O2 has lower vascular-proxy correlation than O1 (0.647 versus 0.740; difference -0.092 [-0.176, -0.014]), whereas O7 has higher correlation (0.847; difference +0.107 [0.037, 0.201]); F2 q is unavailable in the source.',
 'Shared 6 has strong CPTAC protein correlations (O1 0.870; O2 0.875), while the FGF/FGFR group is weaker (0.153; 0.239); all four fixed groups and all trial genotype contrasts are reported.',
 'Gene-wise standardisation upweights full-Angio36 FGF/FGFR coefficients by 2.00-, 7.43-, 1.52- and 6.82-fold in JAVELIN, IMmotion151, CheckMate 025 and CPTAC, respectively; these arithmetic ratios do not establish biological activation or an optimal operator.'
]
report=f'''# Module F — {LABEL}

{tab(headline)}

Author / creator / last modified by: **Chunhao Guo**. All results are **exploratory**. The headline shows the prespecified O1/O2 comparison for Angio36; complete panels below retain every operator, gene, signature and group. Genotype estimates are Cohen's d (altered minus source-negative), protein estimates are Spearman rho; these units must not be pooled. All relative paths resolve within `{O}`. Every row carries `{LABEL}` in its source CSV.

## Reproduction gate: PASS before new F computation

All **344** reproduction checks passed: Angio36 O1/O2 disagreement **48/726** and **80/823**, all available saved pairwise disagreement cells in four cohorts, all 13 available protein point estimates, D's O1/O2 trial scores and all nine VHL/PBRM1/BAP1 genotype count sets. CPTAC Angio36 O1 rho **0.739855**, rounding to **0.740**. The three required genotype gates are JAVELIN VHL **383/308/35**, IMmotion151 PBRM1 **328/374/121**, and CheckMate 025 PBRM1 **63/112/75** (altered/source-negative/unknown). See `results/REPRODUCTION_CHECKS.csv`, `logs/01_reproduction_gate.log`, `results/GENOTYPE_COUNTS.csv`.

No full-signature scores needed regeneration. B1/B2 saved scores are loaded unchanged; F2 correlations, percentile intervals and paired differences are carried forward without recomputation. Protein point-estimate recalculation occurred only within the requested reproduction gate, not as a new F2 analysis.

## Five-line plain statement

'''+ '\n'.join(f'{i+1}. {s}' for i,s in enumerate(statements))+f'''

## Methods and missingness

Exact choices were recorded before F associations in `METHODS_DECISIONS.md`. F1-main uses unadjusted pooled-SD standardised mean differences, with 2,000 unstratified paired patient resamples and fixed scores. CI: percentile, type 7. Each operator-versus-O1 difference uses identical resampled patients. The centred two-sided bootstrap tail p value uses a plus-one correction (minimum 1/2001); BH retains 21 planned slots per trial x signature, including O1 self-comparisons and unavailable tests. These q values are exploratory and conditional on fixed scores, and percentile CIs need not agree exactly with centred bootstrap p values. No reference-fitting uncertainty is included. All F1 O1-O7 x 3 genes x 2 signatures x 3 trials = **126** rows are retained, including nine unavailable Angio36 O3 rows (zero required-gene MAD under B1/B2's rule). F2 adds **14** rows, including unavailable CPTAC Angio36 O3.

Qualified genotype joins and categories are **exactly** module D's `GENOTYPE_PATIENT_JOIN.csv`; their source-defined assay scopes differ (JAVELIN binary WES flags, IMmotion151 FMI alteration strings and assay IDs, CM025 coding nonsynonymous MAF calls). Source-negative never means unassayed or universally mutation-free. No genotype status is redefined. Unknown patients are excluded from the F contrasts and counted below; full RNA-reference scores remain unchanged.

{cnt}

### Specimen-origin sensitivity

`results/F1_SPECIMEN_SENSITIVITY_ALL.csv` contains all **126** planned rows. The sensitivity is the altered coefficient from `score ~ altered + specimen origin`, divided by the pooled within-genotype score SD in the same complete cases; coefficient and denominator are recomputed within each patient bootstrap. IMmotion151 uses D's frozen `qualified_inputs.rds$meta$specimen`; CM025 uses B2's saved specimen manifest, with two missing origins excluded where relevant. Exact genotype-specific exclusions are in the CSV. JAVELIN lacks a qualified origin covariate in reused B1/D inputs (also recorded as unavailable in module E's Table 1 source); its 42 sensitivity rows are explicitly not run. Families are separate from main F1 families. This adjustment does not address all clinical, purity or assay confounding.

### F2 carry-forward and q availability

`results/F2_PROTEIN_CARRIED_FORWARD.csv` preserves every B2 source value and adds only the mandatory follow-up label. All 78 paired cases and the **109-case** RNA score reference are preserved. The proxy is the existing average of standardised CDH5, ENG, PCDH12 and SPARCL1 protein values. B2 did not provide protein p values or bootstrap draws, so **F2 p and BH q are NA**; they are neither inferred from intervals nor recomputed. This is explicitly marked in `results/F1_F2_COMBINED_ANCHORS.csv`.

O1 = mean expression; O2 = gene-standardised mean; O3 = robust gene-standardised mean; O4 = singscore; O5 = ssGSEA (alpha 0.25, normalization off); O6 = Gaussian-kernel GSVA; O7 = expression-oriented PC1. Definitions and source are read-only `../shared/operators.R`; no scoring operator is redefined.

## Complete F1/F2 panel — {LABEL}

{tab(comb)}

## F3 fixed gene-group anchors — {LABEL}

The four groups are reused from D's carried-forward module60 gene weights and cross-checked against the final membership CSV: Shared 6 (6 genes), Angio16-only (10; a fixed comparator outside Angio36), Angio36-only FGF/FGFR (26), and Angio36-only other (4). Exact membership: `results/F3_FIXED_GENE_GROUPS.csv`. O1/O2 are applied separately to each complete group, with full-cohort references; they are not group contributions summed to the complete-signature standardised score. F3 trial inputs reuse D's qualified expression subsets; CPTAC uses B2's 109-case full RNA reference and existing 78-case protein proxy. No genes are dropped, imputed or relabelled. The group bootstrap retains source input patient order and seed 20261113 per group/anchor; F1's score pivots use sorted IDs.

`results/F3_GROUP_ANCHORS_ALL.csv` includes all **56** rows: 4 groups x 2 operators x (2 genotype anchors x 3 trials + 1 CPTAC protein anchor). The F3 q family is descriptive and separate: four O2-minus-O1 group contrasts per cohort x anchor; O1 self q=1. BAP1 is required in F1 but is not an F3 anchor in the addendum. None of the group results is an operator recommendation.

{tab(f3compact,True)}

### Arithmetic upweighting

Every one of the 46 union genes is retained for each of four cohorts in `results/F3_ALL_GENE_WEIGHTS.csv` (**184 rows**); Angio16-only genes have NA full-Angio36 weights because they are not members. The complete Angio36 O1 coefficient is `1/(36*SD(raw mean))`; O2 is `1/(36*gene SD*SD(raw mean gene-z))`. Group sum ratios, SDs and median ratios are in `results/F3_GROUP_WEIGHT_SUMMARY.csv` (**16 rows**). These ratios are coefficient arithmetic on the released expression scale, not variance shares or independent-gene inference. CPTAC weights use the B2 109-case reference, distinct from D's historical 78-case paired-only weight reference. Shared gene ratios are not uniformly reduced across cohorts: CM025 Shared 6 sum ratio is 1.76. This preserves a material opposite numerical pattern.

## Diagnostics, validation and runtime

`plots/F1_Angio36_O2_minus_O1.png` and `plots/F3_CPTAC_group_anchors.png` are original-design, single-panel R diagnostics with matching PDFs, not final submission figures or published-panel reproductions. The first plots all nine Angio36 O2-minus-O1 genotype differences; the second plots all four CPTAC groups under both operators. Source CSVs: `results/PLOT_SOURCE_F1.csv`, `results/PLOT_SOURCE_F3.csv`. `plots/anchor_explorer.html` contains aggregate-only estimates and complete operator/group access without network requests.

`results/VALIDATION_CHECKS.csv` independently verifies all 117 available F1 pooled-SD point estimates with NumPy, all 48 F3 genotype estimates with NumPy, all eight group-protein estimates with pandas tied ranks, source carry-forward identity, planned row counts, missingness, all 2,000 bootstrap replicates per estimable result, self differences, interval order and labels. `logs/05_bootstrap_validation.log` and `results/BOOTSTRAP_VALIDATION.csv` report 520 passing stored-bootstrap, paired-interval, p-value, BH and direct-group-score checks. Input SHA256 hashes: `results/INPUT_MANIFEST.csv`; output hashes: `results/OUTPUT_MANIFEST.csv`. R/software versions: `sessionInfo.txt` and `logs/03_plot_sessionInfo.txt`.

F1-F3 compute completed in **36.88 seconds wall time**, maximum resident set **282,198,016 bytes** (about 269 MiB), peak memory footprint **256,459,640 bytes**, zero swaps (`logs/02_run_F1_F3.log`, macOS `/usr/bin/time -l`). One R process and one numerical thread were used. No compute escalation, downsampling, alternative scoring method, installation, outcome model or external transfer occurred. Scratch inputs/bootstrap arrays: `/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/F/`; shared library reuse: `/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B1/Rlib`.

## Deviations and anything not run

- No deviation from requested score operators, genotype definitions, fixed groups, full-cohort references, or 2,000 resamples; unspecified statistical details were frozen before associations.
- JAVELIN specimen adjustment: not run, missing qualified data field; all not-run rows retained. This is not a compute limitation.
- Angio36 O3: unavailable in all four cohorts under the inherited zero-MAD rule; no substitute method or omitted gene.
- F2 p/q: unavailable in the carried source; carry-forward-only instruction respected. Original estimates and intervals are complete.
- No outcomes or operator selection were run, as forbidden by the addendum. No important F1/F3 main estimand remains unrun for computing reasons.
- Initial diagnostic Cairo export lacked X11 libraries; the existing native R PDF device with Helvetica succeeded, and the unchanged-data diagnostics passed visual, PDF font-size and collision checks. See `DEVIATIONS.md`. Final manuscript figure design was outside this module; diagnostic render checks do not certify publication readiness.

## What this does and does not show

The independently measured anchors make operator-dependent construct tracking observable across genotype and protein assays.
The first-line genotype attenuation is not universal: CM025 PBRM1 has the opposite difference and several intervals include zero.
Gene-group anchors separate strong shared vascular signal from less consistent FGF/FGFR tracking, while coefficient ratios explain arithmetic upweighting only.
Assay scope, specimen origin, tissue composition and clinical setting remain competing explanations; these observational associations are not causal mechanisms.
All analyses were specified after B2/D were seen, are exploratory, use fixed-score uncertainty and are not untouched external validation or evidence for an optimal operator.
'''
(O/'RESULTS.md').write_text(report)

# Aggregate-only self-contained explorer; R remains the source of scientific plotted marks in PDFs/PNGs.
rows=pd.concat([comb,f3],ignore_index=True).replace({np.nan:None}).to_dict(orient='records')
payload=json.dumps(rows,allow_nan=False)
page='''<!doctype html><html><head><meta charset="utf-8"><meta name="author" content="Chunhao Guo"><style>
:root{color-scheme:light dark}body{margin:0;font-family:var(--font-sans);color:var(--foreground)}select{padding:7px;background:var(--secondary);color:var(--foreground);border:1px solid var(--border);border-radius:var(--radius);max-width:100%}.controls{display:flex;gap:8px;flex-wrap:wrap}p{font-size:12px;color:var(--muted-foreground)}svg{display:block;width:100%;height:330px}table{font-size:11px;border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}td,th{text-align:left;border-bottom:1px solid var(--border);padding:7px 3px} .scroll{overflow:auto}text{fill:var(--foreground);font-family:var(--font-sans)}label{font-size:12px}#note{min-height:28px}
</style></head><body><p>post hoc follow-up (specified after modules B2 and D)<br>事后随访分析（在模块 B2 和 D 完成后制定）</p>
<div class="controls"><label>Panel / 分析 <select id="panel"><option value="full">Full signature / 完整特征集</option><option value="group">Fixed gene groups / 固定基因组别</option></select></label><label>Cohort / 队列 <select id="co"></select></label><label>Anchor / 锚点 <select id="anchor"></select></label><label id="siglabel">Signature / 特征集 <select id="sig"><option>Angio36</option><option>Angio16</option></select></label><label>View / 视图 <select id="view"><option value="estimate">Anchor association / 锚点关联</option><option value="difference">Difference vs O1 / 相较 O1 的差异</option></select></label></div><p id="note"></p><svg viewBox="0 0 710 330" aria-label="Aggregate anchor estimates and percentile 95% confidence intervals"></svg><div class="scroll"><table><thead><tr><th>Operator / 算子</th><th>Estimate / 估计</th><th>95% CI / 区间</th><th>Δ vs O1 / 差异 (95% CI)</th><th>BH q</th></tr></thead><tbody></tbody></table></div><p>Fixed-score patient/case bootstrap; 2,000 draws; exploratory.<br>固定评分的患者自助抽样，共 2,000 次；结果为探索性。</p><script>
const data=PAYLOAD;const $=id=>document.getElementById(id);const cohorts=['JAVELIN','IMmotion151','CheckMate025','CPTAC'];$('co').innerHTML=cohorts.map(x=>`<option>${x}</option>`).join('');
function anchors(){const co=$('co').value,group=$('panel').value==='group';$('anchor').innerHTML=(co==='CPTAC'?['Vascular protein proxy']:group?['VHL','PBRM1']:['VHL','PBRM1','BAP1']).map(x=>`<option>${x}</option>`).join('');$('siglabel').style.display=group?'none':'block';draw()}
const f=x=>x===null||x===undefined?'NA':Number(x).toFixed(3);const interval=(a,b)=>a===null?'NA':`[${f(a)}, ${f(b)}]`;
function draw(){const group=$('panel').value==='group';const rows=data.filter(r=>(group?r.module==='F3':r.module!=='F3')&&r.cohort===$('co').value&&r.anchor===$('anchor').value&&(group||r.signature===$('sig').value));let val=$('view').value;const lo=val==='estimate'?'lower':'difference_lower',hi=val==='estimate'?'upper':'difference_upper',v=val==='estimate'?'estimate':'difference_vs_O1';let nums=rows.flatMap(r=>[r[lo],r[hi]]).filter(x=>x!==null);let min=Math.min(0,...nums)-.08,max=Math.max(.01,...nums)+.08;const left=220,right=685;const x=n=>left+(n-min)/(max-min)*(right-left),y=i=>35+i*30;const stroke='var(--chart-1)';let svg=`<line x1="${x(0)}" x2="${x(0)}" y1="12" y2="275" stroke="var(--border)" stroke-dasharray="4 3"/>`;rows.forEach((r,i)=>{let label=group?r.gene_group+' '+r.operator:r.operator;svg+=`<text x="210" y="${y(i)+4}" text-anchor="end" font-size="11">${label}</text>`;if(r[v]!==null){svg+=`<line x1="${x(r[lo])}" x2="${x(r[hi])}" y1="${y(i)}" y2="${y(i)}" stroke="${stroke}" stroke-width="2"/><circle cx="${x(r[v])}" cy="${y(i)}" r="4" fill="${stroke}"/>`}else svg+=`<text x="${left}" y="${y(i)+4}" font-size="11">Unavailable / 不可用</text>`});for(let i=0;i<=4;i++){let n=min+(max-min)*i/4;svg+=`<text x="${x(n)}" y="300" text-anchor="middle" font-size="11">${n.toFixed(2)}</text>`}document.querySelector('svg').innerHTML=svg;document.querySelector('tbody').innerHTML=rows.map(r=>`<tr><td>${group?r.gene_group+' ':''}${r.operator}</td><td>${f(r.estimate)}</td><td>${interval(r.lower,r.upper)}</td><td>${f(r.difference_vs_O1)} ${interval(r.difference_lower,r.difference_upper)}</td><td>${f(r.q_BH)}</td></tr>`).join('');$('note').textContent=$('co').value==='CPTAC'?'Spearman rho; 78 paired cases; F2 q unavailable in source. / Spearman 相关；78 例配对；F2 源表无 q 值。':'Altered minus source-negative, pooled-SD units; unknown genotype excluded. / 改变组减源数据阴性组，单位为合并标准差；排除未知基因型。'}
['co','panel'].forEach(id=>$(id).onchange=anchors);['anchor','sig','view'].forEach(id=>$(id).onchange=draw);anchors();
</script></body></html>'''.replace('PAYLOAD',payload)
(O/'plots/anchor_explorer.html').write_text(page)

inputs=[A/'ANALYSIS_PLAN_LOCKED_20261113.md',A/'ANALYSIS_PLAN_ADDENDUM_F_20261113.md',A/'shared/operators.R',A/'D_biology/results/qualified_inputs.rds',A/'D_biology/results/GENOTYPE_PATIENT_JOIN.csv',A/'D_biology/results/GENOTYPE_JOIN_AUDIT.csv',A/'D_biology/results/PATIENT_SCORES_RANKS.csv',A/'D_biology/results/CARRIED_FORWARD_ALL_GENE_WEIGHTS.csv',ROOT/'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv']
inputs += [A/f'B1_operator_panel_trials/results/{co}_scores_local_only.csv' for co in ['JAVELIN','IMmotion151']]
inputs += [A/'B1_operator_panel_trials/results/pairwise_metrics.csv']
inputs += [A/f'B2_independent_cohorts/results/{nm}.csv' for nm in ['CM025_scores','CPTAC_scores','CM025_pairwise_metrics','CPTAC_pairwise_metrics','CPTAC_protein_proxy','CPTAC_protein_anchor','CM025_sample_manifest']]
inputs += [Path('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2/CPTAC.rds')]
digest=lambda p:hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
wr(pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in inputs]),'INPUT_MANIFEST.csv')
metadata={'author':'Chunhao Guo','creator':'Chunhao Guo','last_modified_by':'Chunhao Guo','analysis_label':LABEL,'claim_status':'exploratory','seed':20261113,'bootstrap_draws':2000,'input_paths_read_only':True,'validation_checks':len(checks),'validation_passed':True}
(R/'OUTPUT_METADATA.json').write_text(json.dumps(metadata,indent=2)+'\n')
wr(pd.DataFrame([dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.csv' and p.suffix not in ['.log']]),'OUTPUT_MANIFEST.csv')
print(json.dumps(metadata,indent=2))
