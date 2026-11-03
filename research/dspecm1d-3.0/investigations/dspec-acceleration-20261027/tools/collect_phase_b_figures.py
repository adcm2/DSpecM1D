#!/usr/bin/env python3
"""Collect completed Phase B figure-job artifacts while later jobs continue."""
import argparse
import base64
import hashlib
import io
import json
import sys
import tarfile
from pathlib import Path

A=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(A.parent/'tunya-profile-20261027'))
import ssh_transport as st

p=argparse.ArgumentParser();p.add_argument('degree',type=int);p.add_argument('destination',type=Path);p.add_argument('--plot-v2',action='store_true');a=p.parse_args()
assert a.degree in [384,1536]
root='/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/phase-b-figures-control/continuation/phase-b-figures-01'
if a.plot_v2:root=root.replace('phase-b-figures-control/continuation/phase-b-figures-01','phase-b-plot-control/continuation/phase-b-plot-02')
source='''import pathlib,json,hashlib,base64,io,tarfile
root=pathlib.Path(ROOT)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def job(name):
 matches=[j for j in json.loads((root/'state.json').read_text())['jobs'] if j['id']==name]
 assert len(matches)==1
 j=matches[0]
 assert j['status']=='completed' and j['returncode']==0 and j['inputs_unchanged'] and j['cleanup']['verified_gone'] and not j.get('failure_reasons')
 return j
names=[f'section-l{DEGREE}',f'slices-l{DEGREE}',f'validate-l{DEGREE}',f'render-l{DEGREE}']
records={n:job(n) for n in names};files=[]
for name,j in records.items():
 for rel,digest in {**j['outputs_sha256'],j['log']:j['log_sha256']}.items():
  path=root/name/rel;assert sha(path)==digest,str(path)
  if path.suffix in ['.json','.csv','.log','.png','.pdf']:
   assert path.stat().st_size<30000000
   files.append(path)
assert records=={n:job(n) for n in names}
receipt=json.dumps({'remote_task':str(root),'degree':DEGREE,'verified_completed_jobs':records},indent=2).encode()
buf=io.BytesIO()
with tarfile.open(fileobj=buf,mode='w:gz') as tf:
 for path in files:tf.add(path,arcname=str(path.relative_to(root)),recursive=False)
 ti=tarfile.TarInfo('collection-receipt.json');ti.size=len(receipt);tf.addfile(ti,io.BytesIO(receipt))
print(base64.b64encode(buf.getvalue()).decode())
'''.replace('ROOT',repr(root)).replace('DEGREE',str(a.degree))
if a.plot_v2:source=source.replace("names=[f'section-l{"+str(a.degree)+"}',f'slices-l{"+str(a.degree)+"}',f'validate-l{"+str(a.degree)+"}',f'render-l{"+str(a.degree)+"}']", "names=['render-v2-l"+str(a.degree)+"']")
blob=base64.b64decode(st.python(source))
a.destination.mkdir(parents=True)
with (a.destination/'snapshot.tar.gz').open('xb') as f:f.write(blob)
with tarfile.open(fileobj=io.BytesIO(blob)) as tf:
    for member in tf.getmembers():
        rel=Path(member.name)
        assert member.isfile() and not rel.is_absolute() and '..' not in rel.parts
        out=a.destination/rel;out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('xb') as f:f.write(tf.extractfile(member).read())
print(json.dumps({'destination':str(a.destination),'archive_sha256':hashlib.sha256(blob).hexdigest(),'archive_bytes':len(blob)}))
