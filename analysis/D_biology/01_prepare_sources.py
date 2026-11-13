# Author: Chunhao Guo
from pathlib import Path
import pandas as pd, numpy as np, openpyxl, json, hashlib, time
P=Path(__file__).resolve().parents[4]; O=P/'84_JCO_PO_Upgrade_20261113/05_analysis/D_biology'; R=O/'results'
genes=['VHL','PBRM1','BAP1','SETD2','KDM5C']; sources=[]
def reg(p,role):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 sources.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest(),role=role))
jp=P/'09_fulltext/downloads/CCRCC-00889/CCRCC-00889_Motzer_2020_Nature_Medicine_JAVELIN_biomarker_tables.xlsx'
w=openpyxl.load_workbook(jp,read_only=True,data_only=True); it=w['S21_Gene_mutation_status'].iter_rows(values_only=True); note=next(it)[0]; head=next(it); ids=head[1:]; found={}
for row in it:
 if row[0] in genes:found[row[0]]=row[1:]
 if len(found)==5:break
assert len(set(ids))==len(ids) and set(found)==set(genes)
jg=pd.DataFrame(found,index=ids).rename_axis('ID').reset_index(); jg.to_csv(R/'source_JAVELIN_genotype.csv',index=False)
jc=pd.read_excel(jp,sheet_name='S11_Clinical_data',header=1); jc[['ID','TCGA_cluster']].to_csv(R/'source_JAVELIN_cluster.csv',index=False)
jw=pd.read_excel(jp,sheet_name='S18_WES_summary_stats',header=1); jw.to_csv(R/'source_JAVELIN_WES_eligibility.csv',index=False);w.close();reg(jp,'S21 binary source mutation flags; S18 WES qualification; S11 fixed TCGA clusters')
ip=P/'57_GenomeMedicine_Final_20260927/analysis/source_documents/NIHMS1733314-supplement-Table_S1.xlsx'
i=pd.read_excel(ip,sheet_name='TableS1',header=5); keep=['PATIENT_ID','RNASEQ_SAMPLE_ID','FMI_SAMPLE_ID','NMF_GROUP']+[x for x in i.columns if str(x).startswith('FMI_') and any(str(x)=='FMI_'+g for g in genes)]
i[keep].to_csv(R/'source_IMmotion151_genotype_cluster.csv',index=False)
w=openpyxl.load_workbook(ip,read_only=True,data_only=True); notes=[]
for row in w['TableS1'].iter_rows(values_only=True):
 for x in row:
  if isinstance(x,str) and (('FMI' in x or 'mutation' in x.lower() or 'alteration' in x.lower() or 'Foundation' in x) and len(x)>25):notes.append(x)
w.close();reg(ip,'Original TableS1 FMI targeted-panel calls, assay IDs and seven fixed NMF classes')
bp=P/'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized/source/Braun2020_supplement.xlsx'
c=pd.read_excel(bp,sheet_name='S1_Clinical_and_Immune_Data',header=1); cols=['SUBJID','Cohort','MAF_Tumor_ID','MAF_Normal_ID','RNA_ID','Purity','TMB_Counts']+genes
c[cols].to_csv(R/'source_CheckMate_WES_crosswalk.csv',index=False)
w=openpyxl.load_workbook(bp,read_only=True,data_only=True);it=w['S2_WES_Data'].iter_rows(values_only=True);bn=next(it)[0]; hh=next(it); hh=list(hh); selected=[];variant_counts={};mafids=set()
for row in it:
 d=dict(zip(hh,row));mid=d['Tumor_Sample_Barcode'];mafids.add(mid); cl=d['Variant_Classification'];variant_counts[cl]=variant_counts.get(cl,0)+1
 if d['Hugo_Symbol'] in genes:selected.append({k:d[k] for k in ['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','Cohort']})
w.close();pd.DataFrame(selected).to_csv(R/'source_CheckMate_target_MAF.csv',index=False);pd.DataFrame({'MAF_Tumor_ID':sorted(mafids)}).to_csv(R/'source_CheckMate_MAF_assayed_ids.csv',index=False);reg(bp,'S1 assay/RNA crosswalk and source flags; S2 released somatic MAF including silent variants')
(O/'logs/source_schema.json').write_text(json.dumps({'JAVELIN_S21_definition':note,'IMmotion151_notes':list(dict.fromkeys(notes)),'IMmotion151_FMI_values':{k:i[k].value_counts(dropna=False).to_dict() for k in keep if k.startswith('FMI_') and k!='FMI_SAMPLE_ID'},'CheckMate_S2_definition':bn,'CheckMate_variant_classes':variant_counts},indent=2,default=str))
(O/'results/source_manifest.json').write_text(json.dumps(sources,indent=2));print('Source extraction complete; IMmotion151 fields:',keep);print('FMI values:',{k:i[k].value_counts(dropna=False).to_dict() for k in keep if k.startswith('FMI_') and k!='FMI_SAMPLE_ID'});print('IMmotion151 notes:',list(dict.fromkeys(notes)))
