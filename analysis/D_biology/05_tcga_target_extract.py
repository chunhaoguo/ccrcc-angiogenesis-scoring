# Author: Chunhao Guo
from pathlib import Path
import pandas as pd,json
P=Path(__file__).resolve().parents[4];R=P/'84_JCO_PO_Upgrade_20261113/05_analysis/D_biology/results'
g=pd.read_csv(P/'83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv');want=set(g.loc[g.Angio16 | g.Angio36,"gene"]);p=Path('/LOCAL_DATA_ROOT/Yang_RCC_metabolic_immune_plasticity_2026/02_data_inventory/raw/tcga_kirc/TCGA.KIRC.HiSeqV2.gz');rows=[]
for d in pd.read_csv(p,sep='\t',index_col=0,chunksize=512):
 d=d.rename(index={'CXCR7':'ACKR3','FIGF':'VEGFD'});rows.append(d[d.index.isin(want)])
x=pd.concat(rows);assert not x.index.duplicated().any();assert want.issubset(x.index)
cols=sorted([c for c in x.columns if c.endswith('-01')]);patients=[c[:12] for c in cols];assert len(set(patients))==len(cols)
x[cols].to_csv(R/'TCGA_PRIMARY_TARGET_EXPRESSION.csv');(R/'TCGA_SELECTION.json').write_text(json.dumps(dict(source=str(p),n_measurements=len(x.columns),n_primary_patients=len(cols),units='released log2(normalized RSEM counts+1), no second transform',selection='one -01 primary tumour measurement per 12-character patient ID; uniqueness verified; source scale previously fingerprinted in folder42/code/prepare.py'),indent=2));print(len(cols),'unique primary patients')
