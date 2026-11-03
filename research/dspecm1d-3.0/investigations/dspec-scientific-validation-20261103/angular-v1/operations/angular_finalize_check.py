#!/usr/bin/env python3
"""Freeze the completed angular study producer state and launch its separate checker."""

from pathlib import Path
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
TRANSPORT = Path(__file__).resolve().parents[3] / "tunya-profile-20261027"
RUN_PLAN = ROOT / "configs/angular-run-v1-plan.json"
CHECK_PLAN = ROOT / "configs/angular-terminal-check-v2-plan.json"
RECEIPTS = ROOT / "receipts"
STAGER = ROOT / "operations/angular_run_stage.py"
FINALIZE_RECEIPT = RECEIPTS / "angular-run-v1-finalization.json"
CHECK_LAUNCH_RECEIPT = RECEIPTS / "angular-terminal-check-v1-launch.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resource_module():
    spec=importlib.util.spec_from_file_location('angular_resource',ROOT/'operations/resource_admission.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    assert sha(ROOT/'operations/resource_admission.py')=='47311038ba6555e2df4331dc2f8f3f6f9360fb120b90d1db28ecfd694c161339'
    return module


def load_payload() -> dict:
    spec=importlib.util.spec_from_file_location('angular_runner',STAGER)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    stage=module.load_payload()
    return {**stage,'run_plan':json.loads(RUN_PLAN.read_text()),
            'check_plan':json.loads(CHECK_PLAN.read_text()),
            'finalizer_sha256':sha(Path(__file__))}


REMOTE = r'''
from pathlib import Path
import hashlib, importlib.util, json, os, subprocess, time
CFG=json.loads(PAYLOAD)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def tree_hash(path):
    spec=importlib.util.spec_from_file_location('build_tool',CFG['build_tool'])
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod.tree_hash(Path(path))
def check_tools():
    for p,h in CFG['tools'].items(): assert sha(p)==h,('tool hash mismatch',p)
def quiet():
    active=[]; roots=(CFG['producer_stage'],CFG['checker_stage'],CFG['run_control'],CFG['check_control'])
    names={'g++','gcc','c++','cc1plus','phase-d-driver','physical-driver','doc-driver','phase-c2-driver'}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit(): continue
        try:
            if p.stat().st_uid!=os.getuid(): continue
            raw=(p/'cmdline').read_bytes().replace(bytes([0]),b' ').decode(errors='replace')
            if not raw: continue
            if Path(raw.split()[0]).name in names or any(r in raw for r in roots) and 'supervise_campaign.py' in raw:
                active.append({'pid':int(p.name),'command':raw[:500]})
        except (FileNotFoundError,ProcessLookupError,PermissionError): pass
    assert not active,('competing process/task',active)
    return active
def verify_build():
    assert Path(CFG['build']).is_dir() and Path(CFG['build_control']).is_dir()
    assert set(CFG['build_evidence']['command_receipts'])=={'dspec-timer-disabled'}
    for key in ('dspec-timer-disabled',):
        variant=CFG['build_evidence']['command_receipts'][key]
        root=Path(CFG['build'])/key; binary=root/'phase-d-driver'; cmd=root/'command.json'
        assert binary.is_file() and cmd.is_file()
        assert sha(binary)==variant['binary_sha256']==variant['binary_sha256_observed']
        rec=json.loads(cmd.read_text())
        assert rec.get('returncode')==0 and rec.get('binary_sha256')==variant['binary_sha256']
        assert rec.get('command_sha256')==variant['command_sha256']
        assert rec.get('tracked_inputs_unchanged') is True
        assert rec.get('tracked_input_sha256')==variant['tracked_input_sha256']==rec.get('tracked_input_sha256_after')
    expected=CFG['build_stage_receipt']['external_dependencies']
    observed={p:tree_hash(p) for p in expected}
    assert observed==expected,'external dependency trees changed after reviewed build'
    for p,h in CFG['build_stage_receipt']['tools'].items(): assert sha(p)==h,p
    return {'binaries':{k:v['binary_sha256'] for k,v in CFG['build_evidence']['command_receipts'].items()},
            'external_dependencies':observed}
def verify_inputs(plan):
    missing=[]
    for job in plan['jobs']:
        for name in job.get('inputs',[]):
            if '{JOB_DIR}' not in name and not Path(name).is_file(): missing.append(name)
    assert not missing,('declared inputs missing',missing)
    return len({name for job in plan['jobs'] for name in job.get('inputs',[]) if '{JOB_DIR}' not in name})
def finalize():
    check_tools(); quiet(); build=verify_build()
    run=CFG['run_plan']; term=CFG['check_plan']
    assert sha(CFG['producer_stage']+'/configs/angular-run-v1-plan.json')==CFG['run_plan_sha256']
    assert sha(CFG['checker_stage']+'/configs/angular-terminal-check-v2-plan.json')==CFG['check_plan_sha256']
    assert run['campaign_root']==CFG['run_control'] and term['campaign_root']==CFG['check_control']
    assert len(run['jobs'])==8 and run['thread_cap']==20
    task=Path(CFG['run_control'])/'continuation'/CFG['run_task']
    state_path=task/'state.json'
    assert state_path.is_file(),'producer state is absent'
    raw=state_path.read_bytes(); state=json.loads(raw)
    assert state.get('status')=='completed' and state.get('plan_sha256')==CFG['run_plan_sha256']
    by_id={j.get('id'):j for j in state.get('jobs',[])}
    expected_jobs={j['id']:j for j in run['jobs']}
    assert set(by_id)==set(expected_jobs) and len(by_id)==8
    output_records=[]
    for job_id, job in expected_jobs.items():
        rec=by_id[job_id]
        assert rec.get('status')=='completed' and rec.get('returncode')==0
        assert not rec.get('failure_reasons') and rec.get('inputs_unchanged') is True
        assert rec.get('cleanup',{}).get('verified_gone') is True
        assert set(rec.get('outputs',[]))==set(job['outputs'])
        out_hashes=rec.get('outputs_sha256',{})
        jobroot=task/job_id
        log=jobroot/rec.get('log','job.log')
        assert log.is_file() and sha(log)==rec.get('log_sha256'),(job_id,'log hash')
        output_records.append({'job':job_id,'path':str(log),'sha256':sha(log),'bytes':log.stat().st_size})
        for rel in job['outputs']:
            path=jobroot/rel; expected=out_hashes.get(rel)
            assert path.is_file() and expected and sha(path)==expected,(job_id,rel)
            output_records.append({'job':job_id,'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
    check_stage=Path(CFG['checker_stage']); snapshot=check_stage/'data/producer-state-finalized.json'
    assert not snapshot.exists(),'immutable snapshot already exists'
    snapshot.parent.mkdir(parents=True,exist_ok=True)
    with snapshot.open('xb') as stream: stream.write(raw)
    assert sha(snapshot)==hashlib.sha256(raw).hexdigest()
    assert state_path.read_bytes()==raw,'producer state changed during finalization'
    count=verify_inputs(term)
    receipt={'status':'completed_and_frozen','producer_task':CFG['run_task'],
        'plan_sha256':CFG['run_plan_sha256'],'producer_state_sha256':sha(state_path),
        'snapshot_path':str(snapshot),'snapshot_sha256':sha(snapshot),'snapshot_bytes':snapshot.stat().st_size,
        'terminal_plan_sha256':CFG['check_plan_sha256'],'terminal_unique_inputs':count,
        'job_count':len(by_id),'verified_outputs_and_logs':len(output_records),
        'outputs':output_records,'build':build,'active_processes':[],'finalized_unix':time.time()}
    return receipt
def launch_checker():
    check_tools(); quiet(); build=verify_build()
    assert FINALIZE_SHA256==sha(CFG['checker_stage']+'/data/producer-state-finalized.json')
    plan_path=Path(CFG['checker_stage'])/'configs/angular-terminal-check-v2-plan.json'
    assert sha(plan_path)==CFG['check_plan_sha256']
    plan=json.loads(plan_path.read_text())
    assert plan['campaign_root']==CFG['check_control'] and plan['thread_cap']==4 and len(plan['jobs'])==1
    unique_inputs=verify_inputs(plan); assert unique_inputs==241
    control=Path(CFG['check_control']); assert not control.exists(),('checker control not fresh',str(control))
    control.mkdir(parents=False,exist_ok=False); continuation=control/'continuation'; continuation.mkdir(exist_ok=False)
    taskroot=continuation/CFG['check_task']; log=continuation/(CFG['check_task']+'.nohup.log')
    launch_record=continuation/(CFG['check_task']+'.launch.json')
    assert not taskroot.exists() and not log.exists() and not launch_record.exists()
    with log.open('xb'): pass
    resource_receipt=resource_launch(['/usr/bin/nohup','/usr/bin/python3',CFG['supervisor'],
        str(plan_path),'--task-id',CFG['check_task']],CFG['checker_stage'],log,'check')
    receipt={'pid':resource_receipt['live']['main_pid'],'resource':resource_receipt,'task_id':CFG['check_task'],'plan_sha256':sha(plan_path),'snapshot_sha256':FINALIZE_SHA256,
        'log':str(log),'stage':CFG['checker_stage'],'control_root':CFG['check_control'],'thread_cap':4,'memory_gib_per_job':8,
        'unique_inputs':unique_inputs,'build_binary_hashes':build['binaries'],'launched_unix':time.time()}
    with launch_record.open('x') as f: json.dump(receipt,f,indent=2);f.write('\n')
    return receipt
if ACTION=='finalize': RESULT=finalize()
elif ACTION=='launch-checker': RESULT=launch_checker()
else: raise ValueError(ACTION)
print(json.dumps(RESULT,sort_keys=True))
'''


def run(action: str) -> dict:
    spec = importlib.util.spec_from_file_location("ssh_transport", TRANSPORT / "ssh_transport.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load the established Tunya SSH transport")
    transport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(transport)
    review=(ROOT/'reviews/angular-finalize-operator-v1-independent-review-20261103.txt').read_text()
    assert 'PASS' in review and sha(Path(__file__)) in review,'current finalizer not independently sealed'
    payload = load_payload()
    payload["action"] = action
    if action == "launch-checker":
        final = json.loads(FINALIZE_RECEIPT.read_text())
        if final.get("status") != "completed_and_frozen":
            raise ValueError("producer snapshot has not been finalized")
        payload["finalize_sha256"] = final["snapshot_sha256"]
    source = "import json\nPAYLOAD=" + repr(json.dumps(payload)) + "\nACTION=" + repr(action) + "\n"
    if action == "launch-checker":
        source += "FINALIZE_SHA256=" + repr(payload["finalize_sha256"]) + "\n"
    source += resource_module().REMOTE_RESOURCE_CODE + REMOTE
    compile(source, f"<angular-{action}>", "exec")
    raw = transport.python(source)
    result = json.loads(raw)
    result["remote_script_sha256"] = hashlib.sha256(source.encode()).hexdigest()
    result["operator_sha256"] = sha(Path(__file__))
    path = FINALIZE_RECEIPT if action == "finalize" else CHECK_LAUNCH_RECEIPT
    with path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("finalize", "launch-checker"):
        raise SystemExit("usage: angular_finalize_check.py finalize|launch-checker")
    print(json.dumps(run(sys.argv[1]), indent=2, sort_keys=True))
