#!/usr/bin/env python3
"""Author: Chunhao Guo. Final immutable-input and delivery metadata checks."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import struct
import subprocess
import sys
import zlib

ROOT = Path('/LOCAL_PROJECT_ROOT')
OUT = ROOT / '84_JCO_PO_Upgrade_20261113/05_analysis/A_clinical_consequence'
WORK = ROOT / '84_JCO_PO_Upgrade_20261113'
AUTHOR = 'Chunhao Guo'

def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()

def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

paths = [
    WORK / '00_context/SHARED_CONTEXT.md',
    WORK / '05_analysis/ANALYSIS_PLAN_LOCKED_20261113.md',
    WORK / '05_analysis/ANALYSIS_PLAN_LOCK.sha256',
    WORK / '02_resource_inventory/UPGRADE_OPTIONS_RANKED.md',
    WORK / '02_resource_inventory/DATA_INVENTORY.md',
    WORK / '04_figure_audit/FIGURE_PROVENANCE.md',
    ROOT / '83_JTM_Final_CrossReferences_20261102/Additional_File_2_Gene_Membership.csv',
    ROOT / '37_A_clinical_workflow_validation_20260917/local_only/INPUTS.rds',
    ROOT / '37_A_clinical_workflow_validation_20260917/local_only/FOLD_ASSIGNMENTS.rds',
    ROOT / '37_A_clinical_workflow_validation_20260917/code/common.R',
    ROOT / '37_A_clinical_workflow_validation_20260917/code/prepare.R',
    ROOT / '48_Polan_ccRCC_GenomeMedicine_20260922/local_only/HELDOUT_PATIENT_RANKS.csv.gz',
    ROOT / '48_Polan_ccRCC_GenomeMedicine_20260922/local_only/PATIENT_RANKS.csv.gz',
    ROOT / '48_Polan_ccRCC_GenomeMedicine_20260922/new_analysis/aggregate/THRESHOLD_MANIFEST.csv',
    ROOT / '48_Polan_ccRCC_GenomeMedicine_20260922/code/polan_analysis.R',
    ROOT / '70_npj_Visual_Narrative_Final_20261011/figures/source_data/sc/new/HELDOUT_REPEAT_SUMMARY.csv',
    ROOT / '21_shared_datasets/IMmotion151_20201106/raw/IMmotion151_clinical_anon_20201106(1).csv',
    ROOT / '67_Mathematical_Extensions_20261008/code/core.py',
    ROOT / '67_Mathematical_Extensions_20261008/code/run_policy.py',
    ROOT / '67_Mathematical_Extensions_20261008/aggregate/B_policy_summary.csv',
]
for cohort in ['JAVELIN', 'IMmotion151']:
    for signature in ['Angio16', 'Angio36']:
        paths.append(ROOT / f'37_A_clinical_workflow_validation_20260917/local_only/CV_REFERENCES__{cohort}__{signature}.rds')
    for filename in [f'B_patient_policy_{cohort}.csv', f'B_independent_check_inputs_{cohort}.csv', f'{cohort}_clinical.csv']:
        paths.append(ROOT / '67_Mathematical_Extensions_20261008/local_only' / filename)
    for signature in ['Angio16', 'Angio36', 'APM12']:
        paths.append(ROOT / f'67_Mathematical_Extensions_20261008/local_only/{cohort}_{signature}.csv')
    paths.append(ROOT / f'67_Mathematical_Extensions_20261008/qa/B_fit_log_{cohort}.csv')

manifest = [dict(path=str(p), resolved_path=str(p.resolve(strict=True)), bytes=p.stat().st_size, sha256=sha(p)) for p in paths]
write_csv(OUT / 'results/input_manifest.csv', manifest)
lock = (WORK / '05_analysis/ANALYSIS_PLAN_LOCK.sha256').read_text().split()[0]
actual = sha(WORK / '05_analysis/ANALYSIS_PLAN_LOCKED_20261113.md')
assert actual == lock, 'Locked plan checksum changed'

pngs = sorted((OUT / 'plots').glob('*.png'))
assert len(pngs) == 5, f'Expected five PNG diagnostics; found {len(pngs)}'
if importlib.util.find_spec('lxml') is not None:
    subprocess.run([sys.executable, '/LOCAL_USER_PATH', *map(str, pngs)], check=True)
else:
    # Same lossless tEXt-chunk operation as the supplied finalizer, without its
    # unrelated Office/XML dependency. Preserve the compressed image bytes.
    def chunk(kind, value):
        return struct.pack('>I', len(value)) + kind + value + struct.pack('>I', zlib.crc32(kind + value) & 0xffffffff)
    for path in pngs:
        data = path.read_bytes()
        assert data[:8] == b'\x89PNG\r\n\x1a\n'
        result, pos, attributed = bytearray(data[:8]), 8, False
        while pos < len(data):
            length = struct.unpack('>I', data[pos:pos+4])[0]
            kind, value = data[pos+4:pos+8], data[pos+8:pos+8+length]
            if kind == b'IDAT' and not attributed:
                for key in [b'Author', b'Creator']:
                    result.extend(chunk(b'tEXt', key + b'\0' + AUTHOR.encode('ascii')))
                attributed = True
            if not (kind == b'tEXt' and value.split(b'\0', 1)[0] in [b'Author', b'Creator']):
                result.extend(data[pos:pos+12+length])
            pos += length + 12
        assert attributed
        path.write_bytes(result)
    print('Supplied finalizer requires unavailable lxml in this Python; used the identical lossless PNG chunk operation with the standard library. No installation.')
metadata = []
for path in pngs:
    data = path.read_bytes()
    assert data[:8] == b'\x89PNG\r\n\x1a\n'
    width, height = struct.unpack('>II', data[16:24])
    pos, fields = 8, {}
    while pos < len(data):
        length = struct.unpack('>I', data[pos:pos+4])[0]
        kind, value = data[pos+4:pos+8], data[pos+8:pos+8+length]
        if kind == b'tEXt':
            key, text = value.split(b'\0', 1)
            fields[key.decode('latin1')] = text.decode('latin1')
        pos += length + 12
    assert fields.get('Author') == fields.get('Creator') == AUTHOR
    assert width >= 1400 and height >= 800
    metadata.append(dict(file=str(path.relative_to(OUT)), width=width, height=height, author=fields['Author'], creator=fields['Creator'], verified=True))
write_csv(OUT / 'results/file_metadata_checks.csv', metadata)

files = sorted(p for p in OUT.rglob('*') if p.is_file())
assert all(p.stat().st_size < 50 * 1024 * 1024 for p in files), 'Large intermediate needs external scratch'
rows = [dict(file=str(p.relative_to(OUT)), bytes=p.stat().st_size, sha256=sha(p)) for p in files if p.name != 'output_manifest.csv']
write_csv(OUT / 'results/output_manifest.csv', rows)
print(f'PASS: locked plan SHA-256 unchanged; {len(paths)} input files hashed; five PNG metadata checks passed; {len(rows)} output hashes saved.')
print(f'Total output size: {sum(p.stat().st_size for p in files):,} bytes; largest file: {max(p.stat().st_size for p in files):,} bytes.')
