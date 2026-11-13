# Author / last modified by: Chunhao Guo
from pathlib import Path
import csv, json, re, hashlib
import openpyxl
P=Path('ccRCC_signature_transportability_benchmark_2026').resolve();O=P/'84_JCO_PO_Upgrade_20261113/05_analysis/B1_operator_panel_trials'
p=P/'09_fulltext/downloads/CCRCC-00889/CCRCC-00889_Motzer_2020_Nature_Medicine_JAVELIN_biomarker_tables.xlsx'
w=openpyxl.load_workbook(p,read_only=True,data_only=True);sets={};loc={}
for i,r in enumerate(w['S15_Pathway_genesets'].iter_rows(values_only=True),1):
 if r[0] in ['B9991003_Javelin_Renal_101_genes26','MCDERMOTT_TEFF','MCDERMOTT_MYELOID_INFLAMMATION']:
  sets[r[0]]=[str(x) for x in r[2:] if x is not None];loc[r[0]]=f'{p}#S15_Pathway_genesets!C{i}:last_nonempty'
textp=P/'09_fulltext/downloads/CCRCC-01007/CCRCC-01007_McDermott_2018_Nature_Medicine_IMmotion150_fulltext.txt'
t=textp.read_text();line=next(x for x in t.splitlines() if x.startswith('Gene signatures were defined as follows:'))
teff=re.search(r'Teff\d*: (.*?);',line).group(1);myeloid=re.search(r'myeloid inflammation[^:]*: (.*?)\.',line).group(1)
def parse(s):return [x.strip().replace('IL-6','IL6') for x in s.replace(', and ',', ').replace(' and ',', ').split(',')]
a=parse(teff);b=parse(myeloid);assert a==sets['MCDERMOTT_TEFF'];assert b==sets['MCDERMOTT_MYELOID_INFLAMMATION'];assert len(sets['B9991003_Javelin_Renal_101_genes26'])==26
rows=[]
for name,key in [('Teff5','MCDERMOTT_TEFF'),('Myeloid6','MCDERMOTT_MYELOID_INFLAMMATION'),('Immuno26','B9991003_Javelin_Renal_101_genes26')]:
 for g in sets[key]:rows.append(dict(signature=name,gene=g,source_location=loc[key],crosscheck_source=str(textp)+'#line 193' if name!='Immuno26' else 'JAVELIN S15 row 2 (source-extracted)',normalization_note='IL-6 printed in paper maps to IL6 in source workbook' if g=='IL6' else 'exact source symbol'))
with (O/'results/immune_gene_sources.csv').open('w') as f:
 c=csv.DictWriter(f,fieldnames=list(rows[0]));c.writeheader();c.writerows(rows)
s=w['S16_Pathway_scores'];it=s.iter_rows(min_row=2,values_only=True);h=next(it);keys=['row_id']+list(sets);ix=[h.index(k) for k in keys]
with (O/'results/javelin_source_score_anchors.csv').open('w') as f:
 c=csv.writer(f);c.writerow(keys)
 for r in it:c.writerow([r[i] for i in ix])
w.close()
with (O/'results/source_hashes.json').open('w') as f:json.dump({str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [p,textp,P/'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv',P/'37_A_clinical_workflow_validation_20260917/local_only/INPUTS.rds',P/'37_A_clinical_workflow_validation_20260917/local_only/FOLD_ASSIGNMENTS.rds']},f,indent=2)
print(json.dumps(sets,indent=2));print('Source binding and S15/paper crosscheck PASS')
