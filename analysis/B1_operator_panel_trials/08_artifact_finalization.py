# Author / creator / last modified by: Chunhao Guo
from pathlib import Path
import hashlib,json
from PIL import Image,PngImagePlugin
P=Path('ccRCC_signature_transportability_benchmark_2026').resolve();O=P/'84_JCO_PO_Upgrade_20261113/05_analysis/B1_operator_panel_trials';p=O/'results/diagnostic_headline_contrasts.png';im=Image.open(p);meta=PngImagePlugin.PngInfo()
for k in ['Author','Creator','LastModifiedBy']:meta.add_text(k,'Chunhao Guo')
meta.add_text('Description','Diagnostic plot of the locked B1-primary contrasts; fixed-score patient-bootstrap confidence intervals.');im.save(p,pnginfo=meta)
assert all(Image.open(p).info[k]=='Chunhao Guo'for k in ['Author','Creator','LastModifiedBy'])
(O/'results/diagnostic-headline-contrasts.html').write_text('<!doctype html><html><head><meta name="author" content="Chunhao Guo"><meta name="last-modified-by" content="Chunhao Guo"></head><body style="margin:0"><img style="width:100%;height:auto" src="'+str(p)+'" alt="B1 primary operator disagreement and fixed-score bootstrap confidence intervals. B1 主分析评分方法分歧及固定评分自举置信区间。"></body></html>')
files=[]
for p in sorted(O.rglob('*')):
 if p.is_file()and p.name!='OUTPUT_MANIFEST.json':files.append(dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(O/'OUTPUT_MANIFEST.json').write_text(json.dumps(dict(author='Chunhao Guo',creator='Chunhao Guo',last_modified_by='Chunhao Guo',files=files),indent=2));print('Attribution verified; manifest written',len(files),'files')
