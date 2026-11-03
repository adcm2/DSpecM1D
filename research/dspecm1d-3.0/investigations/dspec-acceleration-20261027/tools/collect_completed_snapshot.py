#!/usr/bin/env python3
"""Collect hash-verified small artifacts; retain original large binaries remotely."""
import argparse,base64,hashlib,io,json,sys,tarfile
from pathlib import Path
A=Path(__file__).resolve().parents[1];sys.path.insert(0,str(A.parent/'tunya-profile-20261027'));import ssh_transport as st
p=argparse.ArgumentParser();p.add_argument('remote_task');p.add_argument('name');args=p.parse_args();dest=A/'receipts'/args.name;dest.mkdir()
source='''import pathlib,json,hashlib,io,base64,tarfile
p=pathlib.Path(TASK);before=(p/'state.json').read_bytes();d=json.loads(before);assert d['status']=='completed'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
files=[p/'state.json'];verified={}
for j in d['jobs']:
 assert j['status']=='completed' and j['returncode']==0 and j['inputs_unchanged'] and j['cleanup']['verified_gone'] and not j.get('failure_reasons')
 for rel,digest in {**j['outputs_sha256'],j['log']:j['log_sha256']}.items():
  x=p/j['id']/rel;assert sha(x)==digest,str(x);verified[str(x.relative_to(p))]=digest
  if x.suffix in ('.json','.csv','.log') and x.stat().st_size<20000000:files.append(x)
assert (p/'state.json').read_bytes()==before
buf=io.BytesIO()
with tarfile.open(fileobj=buf,mode='w:gz') as tf:
 for x in files:tf.add(x,arcname=str(x.relative_to(p)),recursive=False)
 value=json.dumps({'remote_task':str(p),'verified_artifact_sha256':verified,'state_sha256':hashlib.sha256(before).hexdigest()},indent=2).encode();ti=tarfile.TarInfo('snapshot-verification.json');ti.size=len(value);tf.addfile(ti,io.BytesIO(value))
print(base64.b64encode(buf.getvalue()).decode())
'''.replace('TASK',repr(args.remote_task))
blob=base64.b64decode(st.python(source));archive=dest/'snapshot.tar.gz';archive.write_bytes(blob)
with tarfile.open(fileobj=io.BytesIO(blob)) as tf:
 for member in tf.getmembers():
  rel=Path(member.name);assert member.isfile() and not rel.is_absolute() and '..' not in rel.parts
  target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True)
  with target.open('xb') as f:f.write(tf.extractfile(member).read())
print(json.dumps({'destination':str(dest),'archive_sha256':hashlib.sha256(blob).hexdigest(),'archive_bytes':len(blob)}))
