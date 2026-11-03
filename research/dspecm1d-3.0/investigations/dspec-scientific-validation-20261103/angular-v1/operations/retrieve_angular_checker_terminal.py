#!/usr/bin/env python3
"""Read-only retrieval of angular-v1 standalone checker terminal evidence."""
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
PLAN = ROOT / 'configs/angular-terminal-check-v2-plan.json'
LAUNCH = ROOT / 'receipts/angular-terminal-check-v1-launch.json'
FINALIZATION = ROOT / 'receipts/angular-run-v1-finalization.json'
ARCHIVE = ROOT / 'receipts/angular-checker-terminal-v2.tar.gz'
TARGET = ROOT / 'evidence/checker-terminal-v2'
RETRIEVAL = ROOT / 'receipts/angular-checker-terminal-v2-retrieval.json'
TASK = 'dspec-science-angular-v1-check-v1'
RUN_ROOT = BASE + '/dspec-science-angular-v1-check-control/continuation/' + TASK
STAGE = BASE + '/scientific-inputs/dspec-science-angular-v1-check'


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def load_transport():
    spec = importlib.util.spec_from_file_location('ssh_transport', TRANSPORT / 'ssh_transport.py')
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load established Tunya SSH transport')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REMOTE = r'''
from pathlib import Path
import hashlib, io, json, sys, tarfile
CFG=json.loads(PAYLOAD)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
task=Path(CFG['task_root']); stage=Path(CFG['stage']); state_path=task/'state.json'
assert state_path.is_file(), 'checker terminal state is absent'
state_raw=state_path.read_bytes(); state=json.loads(state_raw)
assert state.get('plan_sha256')==CFG['plan_sha256']
assert state.get('task_id') in (None,CFG['task_id'])
jobs={j.get('id'):j for j in state.get('jobs',[])}
assert set(jobs)=={CFG['job_id']},('checker job inventory',sorted(jobs))
rec=jobs[CFG['job_id']]
files={'state.json':state_path}
snapshot=Path(CFG['snapshot_path'])
assert snapshot.is_file() and sha(snapshot)==CFG['snapshot_sha256']
files['producer-state-finalized.json']=snapshot
jobdir=task/CFG['job_id']
tasklog=jobdir/(rec.get('log') or 'job.log')
assert tasklog.is_file() and sha(tasklog)==rec.get('log_sha256'), 'checker log differs from finalized state'
files[str(tasklog.relative_to(task))]=tasklog
outputs=rec.get('outputs',[])
for rel in outputs:
    p=jobdir/rel
    assert p.is_file() and sha(p)==rec.get('outputs_sha256',{}).get(rel),('checker output differs from finalized state',rel)
    files[str(p.relative_to(task))]=p
for rel in CFG['extra_stage_files']:
    p=stage/rel
    assert p.is_file(),('missing executed staged file',str(p))
    files['executed-stage/'+rel]=p
launch=Path(CFG['launch_path'])
assert launch.is_file(), 'checker launch receipt is absent'
files['receipts/checker-launch.json']=launch
manifest={name:{'sha256':sha(p),'bytes':p.stat().st_size} for name,p in files.items()}
data=io.BytesIO()
with tarfile.open(fileobj=data,mode='w:gz') as tar:
    for name,p in sorted(files.items()): tar.add(p,arcname=name,recursive=False)
    mraw=(json.dumps(manifest,sort_keys=True,indent=2)+'\n').encode()
    info=tarfile.TarInfo('manifest.json');info.size=len(mraw);tar.addfile(info,io.BytesIO(mraw))
sys.stdout.buffer.write(data.getvalue())
'''


def main():
    if ARCHIVE.exists() or TARGET.exists() or RETRIEVAL.exists():
        raise FileExistsError('checker retrieval outputs exist; preserve them and choose a new version')
    plan=json.loads(PLAN.read_text(encoding='utf-8'))
    launch=json.loads(LAUNCH.read_text(encoding='utf-8'))
    final=json.loads(FINALIZATION.read_text(encoding='utf-8'))
    job=plan['jobs'][0]
    if (len(plan.get('jobs',[]))!=1 or job.get('id')!=TASK or
            sha(PLAN)!=launch.get('plan_sha256') or launch.get('task_id')!=TASK or
            launch.get('control_root')!=plan['campaign_root'] or
            final.get('status')!='completed_and_frozen' or
            final.get('snapshot_sha256')!=launch.get('snapshot_sha256')):
        raise ValueError('checker launch/plan/producer finalization binding mismatch')
    stage=STAGE
    command=job['command']
    staged=[]
    for flag in ('--helpers',):
        if flag in command:
            staged.append(Path(command[command.index(flag)+1]).relative_to(stage).as_posix())
    staged.extend(['tools/check_angular_v1.py','configs/angular-terminal-check-v2-plan.json'])
    if len(staged)!=3 or len(set(staged))!=3:
        raise ValueError('unexpected checker/helper staged source inventory')
    payload={'task_root':RUN_ROOT,'stage':STAGE,'task_id':TASK,'job_id':TASK,
             'plan_sha256':sha(PLAN),'launch_path':launch['log'].rsplit('/',1)[0]+'/'+TASK+'.launch.json',
             'snapshot_path':final['snapshot_path'],'snapshot_sha256':final['snapshot_sha256'],
             'extra_stage_files':staged}
    source='import json\nPAYLOAD='+repr(json.dumps(payload,sort_keys=True))+'\n'+REMOTE
    compile(source,'<angular-v1-checker-terminal-readonly>','exec')
    transport=load_transport()
    result=subprocess.run(['ssh',*transport.OPTIONS,'-T','tunya','python3','-'],
                          input=source.encode(),capture_output=True,check=True)
    raw=result.stdout
    with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as tar:
        members=tar.getmembers(); names=[m.name for m in members]
        if any(not m.isfile() or Path(m.name).is_absolute() or '..' in Path(m.name).parts for m in members):
            raise ValueError('unsafe/non-file terminal archive member')
        if len(names)!=len(set(names)) or 'manifest.json' not in names:
            raise ValueError('duplicate archive member or missing manifest')
        manifest=json.load(tar.extractfile('manifest.json'))
        if set(names)!=set(manifest)|{'manifest.json'}: raise ValueError('terminal archive/manifest inventory mismatch')
        for m in members:
            if m.name=='manifest.json': continue
            body=tar.extractfile(m).read(); item=manifest[m.name]
            if len(body)!=item['bytes'] or sha_bytes(body)!=item['sha256']:
                raise ValueError('terminal member hash/size mismatch: '+m.name)
        ARCHIVE.parent.mkdir(parents=True,exist_ok=True)
        with ARCHIVE.open('xb') as f: f.write(raw)
        TARGET.mkdir(parents=True)
        for m in members:
            dest=TARGET/m.name; dest.parent.mkdir(parents=True,exist_ok=True)
            with dest.open('xb') as f: f.write(tar.extractfile(m).read())
    task_state=json.loads((TARGET/'state.json').read_text())
    jobrec=next(item for item in task_state['jobs'] if item.get('id')==TASK)
    receipt={'archive_sha256':sha(ARCHIVE),'archive_bytes':ARCHIVE.stat().st_size,
             'verified_files':len(manifest),'checker_plan_sha256':sha(PLAN),
             'checker_launch_receipt_sha256':sha(LAUNCH),
             'producer_finalization_sha256':sha(FINALIZATION),
             'producer_state_sha256':final['snapshot_sha256'],
             'checker_status':jobrec.get('status'),'checker_returncode':jobrec.get('returncode'),
             'checker_inputs_unchanged':jobrec.get('inputs_unchanged'),
             'checker_cleanup':jobrec.get('cleanup'),
             'remote_mutations':False,'local_root':str(TARGET)}
    with RETRIEVAL.open('x',encoding='utf-8') as f:
        json.dump(receipt,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(receipt,sort_keys=True))


if __name__=='__main__': main()
