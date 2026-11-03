#!/usr/bin/env python3
"""Read sealed producer artifacts over SSH; write only new local evidence."""
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('transport', ROOT.parents[1] / 'tunya-profile-20261027/ssh_transport.py')
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)
REMOTE = r'''
from pathlib import Path
import hashlib, io, json, sys, tarfile
root=Path('/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/dspec-science-study2-v1-control/continuation/dspec-science-study2-v1-run-v1')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
state=root/'state.json'
assert sha(state)=='9cd09f3e4bd0979a39a240cfa8d2a1a2bf631750ab50a8869a55db49cdfad704'
s=json.loads(state.read_text()); assert s['status']=='completed' and len(s['jobs'])==11
files={'state.json':state}
for job in s['jobs']:
    assert job['status']=='completed' and job['returncode']==0 and job['inputs_unchanged'] and job['cleanup']['verified_gone']
    for name in job['outputs']:
        path=root/job['id']/name; assert sha(path)==job['outputs_sha256'][name]
        files[str(path.relative_to(root))]=path
    path=root/job['id']/job['log']; assert sha(path)==job['log_sha256']
    files[str(path.relative_to(root))]=path
manifest={name:{'sha256':sha(p),'bytes':p.stat().st_size} for name,p in files.items()}
data=io.BytesIO()
with tarfile.open(fileobj=data,mode='w:gz') as tar:
    for name,p in sorted(files.items()): tar.add(p,arcname=name,recursive=False)
    raw=(json.dumps(manifest,sort_keys=True,indent=2)+'\n').encode()
    item=tarfile.TarInfo('manifest.json');item.size=len(raw);tar.addfile(item,io.BytesIO(raw))
sys.stdout.buffer.write(data.getvalue())
'''

def main():
    archive = ROOT / 'receipts/study2-producer-raw-v1.tar.gz'
    target = ROOT / 'evidence/producer-raw-v1'
    assert not archive.exists() and not target.exists()
    result = subprocess.run(['ssh', *transport.OPTIONS, '-T', 'tunya', 'python3 -'],
                            input=REMOTE.encode(), capture_output=True, check=True)
    payload = result.stdout
    with archive.open('xb') as stream: stream.write(payload)
    with tarfile.open(fileobj=io.BytesIO(payload), mode='r:gz') as tar:
        members = tar.getmembers()
        assert all(m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in members)
        assert len({m.name for m in members}) == len(members)
        manifest=json.load(tar.extractfile('manifest.json'))
        assert {m.name for m in members} == set(manifest) | {'manifest.json'}
        target.mkdir(parents=True)
        for member in members:
            raw=tar.extractfile(member).read()
            if member.name!='manifest.json':
                info=manifest[member.name]
                assert len(raw)==info['bytes'] and hashlib.sha256(raw).hexdigest()==info['sha256']
            path=target/member.name; path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as stream: stream.write(raw)
    receipt={'archive_sha256':hashlib.sha256(payload).hexdigest(),'archive_bytes':len(payload),
             'verified_files':len(manifest),'producer_state_sha256':manifest['state.json']['sha256'],
             'remote_mutations':False,'local_root':str(target)}
    with (ROOT/'receipts/study2-producer-raw-v1-retrieval.json').open('x') as stream:
        json.dump(receipt,stream,indent=2);stream.write('\n')
    print(json.dumps(receipt))

if __name__=='__main__': main()
