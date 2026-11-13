"""Author / creator / last modified by: Chunhao Guo.
Stream every gene and all 250 CM025 columns; no expression filtering or transforms.
Standard library only; XML elements are cleared after each row.
"""
from pathlib import Path
import csv, zipfile, xml.etree.ElementTree as ET
ROOT=Path('/LOCAL_PROJECT_ROOT')
EXT=ROOT/'20_harmonized_four_cancer_audit_v1_0/phase5_manuscript_v1_0/external_randomized'
SCR=Path('/LOCAL_EXTERNAL_DRIVE/MacBook_Project_Storage/JCO_PO_Upgrade_scratch_20261113/B2')
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with (EXT/'CM025_CLINICAL.csv').open() as f:
    ids=[r['RNA_ID'] for r in csv.DictReader(f)]
with zipfile.ZipFile(EXT/'source/Braun2020_supplement.xlsx') as z:
    ss=ET.fromstring(z.read('xl/sharedStrings.xml'))
    strings=[''.join(s.itertext()) for s in ss]
    del ss
    with z.open('xl/worksheets/sheet5.xml') as source, (SCR/'CM025_complete_expression.csv').open('w') as out:
        writer=csv.writer(out); header=None; n=0
        for _,row in ET.iterparse(source,events=('end',)):
            if row.tag.rsplit('}',1)[-1]!='row':continue
            vals={}
            for cell in row:
                val=cell.find('m:v',ns)
                if val is None:continue
                text=val.text
                if cell.get('t')=='s':text=strings[int(text)]
                col=''.join(c for c in cell.get('r') if c.isalpha())
                vals[col]=text
            if row.get('r')=='2':
                header={v:k for k,v in vals.items()}
                assert set(ids)<=header.keys() and 'gene_name' in header
                cols=[header['gene_name']]+[header[i] for i in ids]
                writer.writerow(['gene_name']+ids)
            elif header is not None and vals:
                assert all(c in vals for c in cols),row.get('r')
                writer.writerow([vals[c] for c in cols]);n+=1
            row.clear()
        assert n==43893,n
        print('Exact full expression extracted:',n,'genes x',len(ids),'CM025 patients')
        (SCR/'CM025_extraction_complete.txt').write_text(f'{n} genes x {len(ids)} samples\nAuthor: Chunhao Guo\n')
