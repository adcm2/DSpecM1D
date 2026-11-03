#!/usr/bin/env python3
"""Read the finalized angular-v1 producer outputs over the established SSH transport."""
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
TRANSPORT = ROOT.parents[1] / 'tunya-profile-20261027'
BASE = '/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya'
RUN_ROOT = BASE + '/dspec-science-angular-v1-control/continuation/dspec-science-angular-v1-run-v1'
SNAPSHOT = BASE + '/scientific-inputs/dspec-science-angular-v1-check/data/producer-state-finalized.json'
PLAN = ROOT / 'configs/angular-run-v1-plan.json'
TERMINAL_PLAN = ROOT / 'configs/angular-terminal-check-v2-plan.json'
FINALIZATION = ROOT / 'receipts/angular-run-v1-finalization.json'
ARCHIVE = ROOT / 'receipts/angular-producer-raw-v1.tar.gz'
TARGET = ROOT / 'evidence/producer-raw-v1'
RETRIEVAL = ROOT / 'receipts/angular-producer-raw-v1-retrieval.json'


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def load_transport():
    spec = importlib.util.spec_from_file_location('ssh_transport', TRANSPORT / 'ssh_transport.py')
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load the established Tunya SSH transport')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_payload():
    plan = json.loads(PLAN.read_text(encoding='utf-8'))
    final = json.loads(FINALIZATION.read_text(encoding='utf-8'))
    plan_sha = sha(PLAN)
    terminal_sha = sha(TERMINAL_PLAN)
    state_sha = final.get('producer_state_sha256')
    if (final.get('status') != 'completed_and_frozen' or
            final.get('producer_task') != RUN_ROOT.rsplit('/', 1)[-1] or
            final.get('plan_sha256') != plan_sha or
            final.get('terminal_plan_sha256') != terminal_sha or
            final.get('snapshot_sha256') != state_sha or
            final.get('snapshot_path') != SNAPSHOT or
            final.get('job_count') != 8):
        raise ValueError('producer finalization receipt does not bind the exact angular-v1 plan/snapshot')
    job_by_id = {job['id']: job for job in plan.get('jobs', [])}
    if len(job_by_id) != 8:
        raise ValueError('angular-v1 plan job inventory is not eight jobs')
    expected = {}
    for item in final.get('outputs', []):
        path = item.get('path')
        if path in expected or not path or not item.get('sha256') or item.get('bytes', 0) < 0:
            raise ValueError('finalization output inventory has duplicate or invalid records')
        expected[path] = {'sha256': item['sha256'], 'bytes': item['bytes']}
    expected_paths = set()
    for job in plan['jobs']:
        task = RUN_ROOT + '/' + job['id']
        expected_paths.add(task + '/job.log')
        expected_paths.update(task + '/' + output for output in job['outputs'])
    if set(expected) != expected_paths or len(expected) != len(expected_paths):
        raise ValueError('finalization output/log inventory differs from the exact producer plan')
    return {'run_root': RUN_ROOT, 'snapshot_path': SNAPSHOT,
            'snapshot_sha256': state_sha, 'plan_sha256': plan_sha,
            'jobs': {jid: {'outputs': job['outputs']} for jid, job in job_by_id.items()},
            'expected_files': expected}


REMOTE = r'''
from pathlib import Path
import hashlib, io, json, sys, tarfile
CFG=json.loads(PAYLOAD)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
root=Path(CFG['run_root']); state_path=root/'state.json'
snapshot=Path(CFG['snapshot_path'])
raw=state_path.read_bytes(); assert sha(state_path)==CFG['snapshot_sha256']
assert snapshot.is_file() and snapshot.read_bytes()==raw and sha(snapshot)==CFG['snapshot_sha256']
state=json.loads(raw); assert state.get('status')=='completed' and state.get('plan_sha256')==CFG['plan_sha256']
byid={j.get('id'):j for j in state.get('jobs',[])}
assert set(byid)==set(CFG['jobs']) and len(byid)==8
files={'state.json':state_path,'producer-state-finalized.json':snapshot}
for jid,job in CFG['jobs'].items():
    rec=byid[jid]
    assert rec.get('status')=='completed' and rec.get('returncode')==0
    assert not rec.get('failure_reasons') and rec.get('inputs_unchanged') is True
    assert rec.get('cleanup',{}).get('verified_gone') is True
    assert set(rec.get('outputs',[]))==set(job['outputs'])
    for rel in job['outputs']:
        p=root/jid/rel; expected=rec.get('outputs_sha256',{}).get(rel)
        assert p.is_file() and expected and sha(p)==expected
        files[str(p.relative_to(root))]=p
    logrel=rec.get('log','job.log'); log=root/jid/logrel
    assert log.is_file() and sha(log)==rec.get('log_sha256')
    files[str(log.relative_to(root))]=log
manifest={name:{'sha256':sha(p),'bytes':p.stat().st_size} for name,p in files.items()}
for remote,expected in CFG['expected_files'].items():
    rel=str(Path(remote).relative_to(root)); assert rel in manifest
    assert manifest[rel]==expected,(rel,manifest[rel],expected)
assert set(manifest)-{'state.json','producer-state-finalized.json'}=={
    str(Path(path).relative_to(root)) for path in CFG['expected_files']}
data=io.BytesIO()
with tarfile.open(fileobj=data,mode='w:gz') as tar:
    for name,p in sorted(files.items()): tar.add(p,arcname=name,recursive=False)
    raw_manifest=(json.dumps(manifest,sort_keys=True,indent=2)+'\n').encode()
    info=tarfile.TarInfo('manifest.json');info.size=len(raw_manifest);tar.addfile(info,io.BytesIO(raw_manifest))
sys.stdout.buffer.write(data.getvalue())
'''


def main():
    if ARCHIVE.exists() or TARGET.exists() or RETRIEVAL.exists():
        raise FileExistsError('raw retrieval outputs already exist; preserve them and choose a new version')
    payload = make_payload()
    source = 'import json\nPAYLOAD=' + repr(json.dumps(payload, sort_keys=True)) + '\n' + REMOTE
    compile(source, '<angular-v1-read-only-producer-retrieval>', 'exec')
    transport = load_transport()
    result = subprocess.run(['ssh', *transport.OPTIONS, '-T', 'tunya', 'python3 -'],
                            input=source.encode(), capture_output=True, check=True)
    raw_archive = result.stdout
    archive_sha = sha_bytes(raw_archive)
    final = json.loads(FINALIZATION.read_text(encoding='utf-8'))
    if final.get('status') != 'completed_and_frozen' or final.get('snapshot_sha256') != payload['snapshot_sha256']:
        raise ValueError('finalization receipt changed during read-only retrieval')
    with tarfile.open(fileobj=io.BytesIO(raw_archive), mode='r:gz') as tar:
        members = tar.getmembers()
        if any(not m.isfile() or Path(m.name).is_absolute() or '..' in Path(m.name).parts for m in members):
            raise ValueError('retrieved archive contains unsafe/non-file members')
        names = [m.name for m in members]
        if len(names) != len(set(names)) or 'manifest.json' not in names:
            raise ValueError('retrieved archive has duplicate members or lacks its manifest')
        manifest = json.load(tar.extractfile('manifest.json'))
        if set(names) != set(manifest) | {'manifest.json'}:
            raise ValueError('retrieved archive member set differs from manifest')
        for member in members:
            content = tar.extractfile(member).read()
            if member.name != 'manifest.json':
                record = manifest[member.name]
                if len(content) != record['bytes'] or sha_bytes(content) != record['sha256']:
                    raise ValueError('retrieved member hash/size mismatch: ' + member.name)
        # Check every finalization-bound producer file, not only archive-internal hashes.
        for path, expected in payload['expected_files'].items():
            name = str(Path(path).relative_to(Path(payload['run_root'])))
            if manifest.get(name) != expected:
                raise ValueError('retrieved output differs from finalized receipt: ' + name)
        if manifest['state.json']['sha256'] != payload['snapshot_sha256'] or \
                manifest['producer-state-finalized.json']['sha256'] != payload['snapshot_sha256']:
            raise ValueError('retrieved producer state/snapshot hash differs from finalization receipt')
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        with ARCHIVE.open('xb') as stream:
            stream.write(raw_archive)
        TARGET.mkdir(parents=True)
        for member in members:
            dest = TARGET / member.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            with dest.open('xb') as stream:
                stream.write(tar.extractfile(member).read())
    receipt = {'archive_sha256': archive_sha, 'archive_bytes': len(raw_archive),
               'verified_files': len(manifest), 'producer_state_sha256': manifest['state.json']['sha256'],
               'finalization_receipt_sha256': sha(FINALIZATION), 'run_plan_sha256': sha(PLAN),
               'remote_mutations': False, 'local_root': str(TARGET)}
    with RETRIEVAL.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
