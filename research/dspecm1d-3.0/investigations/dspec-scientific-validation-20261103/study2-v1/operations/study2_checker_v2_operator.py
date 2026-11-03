#!/usr/bin/env python3
"""Stage and launch only the reviewed Study 2 checker-v2 retry."""

from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
TRANSPORT = Path(__file__).resolve().parents[3] / "tunya-profile-20261027"
BASE = "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya"
STAGE = BASE + "/scientific-inputs/dspec-science-study2-v1-check-v2"
CONTROL = BASE + "/dspec-science-study2-v1-check-control-v2"
TASK = "dspec-science-study2-v1-check-v2"
PLAN = ROOT / "configs/study2-terminal-check-v2-plan.json"
ARCHIVE = ROOT / "receipts/study2-checker-v2-allowlist.tar.gz"
MANIFEST = ROOT / "receipts/study2-checker-v2-allowlist.json"
RUN_OPERATOR = ROOT / "operations/study2_run_stage.py"
FINALIZATION = ROOT / "receipts/study2-run-v1-finalization.json"
FAILED_LAUNCH = ROOT / "receipts/study2-terminal-check-v1-launch.json"
LOCAL_RECEIPT = ROOT / "receipts/study2-terminal-check-v2-launch.json"
PLAN_SHA = "041eff1e6f75c5e464a88f23d7348c699a096d9219d5f93fff7020fdb1192a0e"
CHECKER_SHA = "d138e8869c349ce8510774e71e425bee0841c002330110775c42d5550ea12b3b"
ARCHIVE_SHA = "eaa7622f8f3857af20fbcff3f55246e9f893932f43ee43bd8778659160e3a9f9"
MANIFEST_SHA = "813b8ae915fd1e8a820824d999415fb7c0214e4c8621b8e23b93e02ea836f14a"
SNAPSHOT_SHA = "9cd09f3e4bd0979a39a240cfa8d2a1a2bf631750ab50a8869a55db49cdfad704"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload() -> dict:
    runner_spec = importlib.util.spec_from_file_location("study2_runner", RUN_OPERATOR)
    if runner_spec is None or runner_spec.loader is None:
        raise RuntimeError("cannot load reviewed Study 2 run operator")
    runner = importlib.util.module_from_spec(runner_spec)
    runner_spec.loader.exec_module(runner)
    config = runner.load_payload()
    assert sha(PLAN) == PLAN_SHA
    assert sha(ARCHIVE) == ARCHIVE_SHA and sha(MANIFEST) == MANIFEST_SHA
    assert sha(ROOT / "tools/check_study2.py") == CHECKER_SHA
    final = json.loads(FINALIZATION.read_text())
    assert final["snapshot_sha256"] == SNAPSHOT_SHA
    manifest = json.loads(MANIFEST.read_text())
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        members = tar.getmembers()
        expected = {item["path"]: item for item in manifest["members"]}
        assert len(members) == len(expected) == 3
        assert {item.name for item in members} == set(expected)
        for item in members:
            record = expected[item.name]
            data = tar.extractfile(item).read()
            assert item.isfile() and item.size == record["bytes"]
            assert hashlib.sha256(data).hexdigest() == record["sha256"]
    return {
        "base": BASE,
        "stage": STAGE,
        "control": CONTROL,
        "task": TASK,
        "plan_sha256": PLAN_SHA,
        "archive_sha256": ARCHIVE_SHA,
        "manifest_sha256": MANIFEST_SHA,
        "manifest": manifest,
        "snapshot_path": config["checker_stage"] + "/data/producer-state-finalized.json",
        "snapshot_sha256": SNAPSHOT_SHA,
        "producer_control": config["run_control"],
        "producer_task": config["run_task"],
        "producer_plan_sha256": config["run_plan_sha256"],
        "failed_control": config["check_control"],
        "failed_task": config["check_task"],
        "build": config["build"],
        "build_control": config["build_control"],
        "build_evidence": config["build_evidence"],
        "build_stage_receipt": config["build_stage_receipt"],
        "build_tool": config["build_tool"],
        "tools": config["tools"],
        "supervisor": config["supervisor"],
        "plan": json.loads(PLAN.read_text()),
        "operator_sha256": sha(Path(__file__)),
    }


REMOTE = r'''
from pathlib import Path
import hashlib, importlib.util, json, os, subprocess, tarfile, time
C=json.loads(PAYLOAD)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tree_hash(p):
    spec=importlib.util.spec_from_file_location('build_tool',C['build_tool'])
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod.tree_hash(Path(p))
def check_tools():
    for p,h in C['tools'].items(): assert sha(p)==h,('tool hash mismatch',p)
def quiet():
    active=[]
    roots=(C['stage'],C['control'],C['producer_control'],C['failed_control'])
    names={'g++','gcc','c++','cc1plus','phase-d-driver','physical-driver','doc-driver','phase-c2-driver'}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit(): continue
        try:
            if p.stat().st_uid!=os.getuid(): continue
            raw=(p/'cmdline').read_bytes().replace(bytes([0]),b' ').decode(errors='replace')
            if not raw: continue
            if Path(raw.split()[0]).name in names or ('supervise_campaign.py' in raw and any(r in raw for r in roots)):
                active.append({'pid':int(p.name),'command':raw[:400]})
        except (FileNotFoundError,ProcessLookupError,PermissionError): pass
    assert not active,('competing process/task',active)
    return active
def verify_build():
    assert Path(C['build']).is_dir() and Path(C['build_control']).is_dir()
    for key in ('dspec-timer-disabled','dspec-timer-enabled'):
        rec=C['build_evidence']['command_receipts'][key]; root=Path(C['build'])/key
        binary=root/'phase-d-driver'; cmd=root/'command.json'
        assert binary.is_file() and cmd.is_file() and sha(binary)==rec['binary_sha256']==rec['binary_sha256_observed']
        command=json.loads(cmd.read_text())
        assert command.get('returncode')==0 and command.get('command_sha256')==rec['command_sha256']
        assert command.get('tracked_input_sha256')==rec['tracked_input_sha256']==command.get('tracked_input_sha256_after')
    expected=C['build_stage_receipt']['external_dependencies']
    observed={p:tree_hash(p) for p in expected}; assert observed==expected
    return {'binaries':{k:v['binary_sha256'] for k,v in C['build_evidence']['command_receipts'].items()},'external_dependencies':observed}
def verify_producer():
    control=Path(C['producer_control'])
    task=control/'continuation'/C['producer_task']; state_path=task/'state.json'; state_raw=state_path.read_bytes(); state=json.loads(state_raw)
    assert state.get('status')=='completed' and state.get('plan_sha256')==C['producer_plan_sha256']
    assert len(state.get('jobs',[]))==11
    for job in state['jobs']:
        assert job.get('status')=='completed' and job.get('returncode')==0
        assert not job.get('failure_reasons') and job.get('inputs_unchanged') is True
        assert job.get('cleanup',{}).get('verified_gone') is True
    snap=Path(C['snapshot_path']); assert sha(state_path)==C['snapshot_sha256']==sha(snap)
    failed=Path(C['failed_control'])/'continuation'/C['failed_task']/'state.json'
    failure=json.loads(failed.read_text()); assert failure.get('status')=='failed'
    fr=failure['jobs'][0]; assert fr.get('returncode')==1 and fr.get('cleanup',{}).get('verified_gone') is True
    return {'producer_state_sha256':sha(state_path),'snapshot_sha256':sha(snap),'job_count':len(state['jobs']),
            'preserved_failed_checker_rc':fr['returncode']}
def init_stage():
    check_tools(); quiet(); build=verify_build(); producer=verify_producer()
    stage=Path(C['stage']); control=Path(C['control'])
    assert not stage.exists() and not control.exists()
    stage.mkdir(parents=False,exist_ok=False)
    return {'status':'fresh_stage_created','stage':str(stage),'control_absent':not control.exists(),
            'build':build,'producer':producer,'active_processes':[]}
def extract():
    check_tools(); quiet(); verify_build(); verify_producer()
    stage=Path(C['stage']); assert stage.is_dir() and not Path(C['control']).exists()
    archive=stage/'study2-checker-v2-allowlist.tar.gz'; manifest_path=stage/'study2-checker-v2-allowlist.json'
    assert archive.is_file() and manifest_path.is_file()
    assert sha(archive)==C['archive_sha256'] and sha(manifest_path)==C['manifest_sha256']
    manifest=json.loads(manifest_path.read_text()); assert manifest==C['manifest']
    expected={r['path']:r for r in manifest['members']}
    with tarfile.open(archive,'r:gz') as tar:
        members=tar.getmembers(); assert len(members)==len(expected)==3
        assert {m.name for m in members}==set(expected)
        for m in members:
            assert m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts
            data=tar.extractfile(m).read(); rec=expected[m.name]
            assert len(data)==rec['bytes'] and sha_bytes(data)==rec['sha256']
            target=stage/m.name; target.parent.mkdir(parents=True,exist_ok=True)
            with target.open('xb') as out: out.write(data)
            assert sha(target)==rec['sha256']
    plan=stage/'configs/study2-terminal-check-v2-plan.json'
    assert sha(plan)==C['plan_sha256']
    parsed=json.loads(plan.read_text()); job=parsed['jobs'][0]
    assert parsed['campaign_root']==C['control'] and parsed['thread_cap']==4 and len(parsed['jobs'])==1
    assert job['id']==C['task'] and job['cwd']==str(stage)
    assert job['command'][1]==str(stage/'tools/check_study2.py')
    missing=[p for p in job['inputs'] if not Path(p).is_file()]
    assert not missing,('checker inputs missing after stage',missing)
    return {'status':'staged','members_verified':len(expected),'unique_inputs':len(set(job['inputs'])),
            'plan_sha256':sha(plan),'snapshot_sha256':sha(C['snapshot_path'])}
def sha_bytes(data): return hashlib.sha256(data).hexdigest()
def launch():
    check_tools(); quiet(); build=verify_build(); producer=verify_producer()
    stage=Path(C['stage']); plan=stage/'configs/study2-terminal-check-v2-plan.json'
    assert sha(plan)==C['plan_sha256']
    parsed=json.loads(plan.read_text()); job=parsed['jobs'][0]
    assert parsed['campaign_root']==C['control'] and len(parsed['jobs'])==1 and parsed['thread_cap']==4
    assert len(set(job['inputs']))==252 and not [p for p in job['inputs'] if not Path(p).is_file()]
    for rec in C['manifest']['members']:
        assert sha(stage/rec['path'])==rec['sha256']
    control=Path(C['control']); assert not control.exists(); control.mkdir(parents=False,exist_ok=False)
    continuation=control/'continuation'; continuation.mkdir(exist_ok=False)
    taskroot=continuation/C['task']; log=continuation/(C['task']+'.nohup.log'); record=continuation/(C['task']+'.launch.json')
    assert not taskroot.exists() and not log.exists() and not record.exists()
    env=os.environ.copy(); env.update(LD_LIBRARY_PATH='/opt/gcc-13.2.0/lib64',LD_PRELOAD='',PYTHONDONTWRITEBYTECODE='1')
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','GOTO_NUM_THREADS'): env[k]='1'
    supervisor=C['supervisor']
    with log.open('xb') as out:
        proc=subprocess.Popen(['/usr/bin/nohup','/usr/bin/python3',supervisor,str(plan),'--task-id',C['task']],
            cwd=stage,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
    receipt={'pid':proc.pid,'task_id':C['task'],'plan_sha256':sha(plan),'stage':str(stage),'control_root':str(control),
        'thread_cap':4,'memory_gib_per_job':8,'unique_inputs':len(set(job['inputs'])),'producer_state_sha256':producer['producer_state_sha256'],
        'build_binary_hashes':build['binaries'],'operator_payload_sha256':C['operator_payload_sha256'],'launched_unix':time.time()}
    with record.open('x') as f: json.dump(receipt,f,indent=2,sort_keys=True); f.write('\n')
    return receipt
if ACTION=='init': RESULT=init_stage()
elif ACTION=='extract': RESULT=extract()
elif ACTION=='launch': RESULT=launch()
else: raise ValueError(ACTION)
print(json.dumps(RESULT,sort_keys=True))
'''


def remote_run(action: str) -> dict:
    spec = importlib.util.spec_from_file_location("ssh_transport", TRANSPORT / "ssh_transport.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load established Tunya transport")
    transport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(transport)
    config = payload()
    config["action"] = action
    config["operator_payload_sha256"] = sha(Path(__file__))
    source = "import json\nPAYLOAD=" + repr(json.dumps(config)) + "\nACTION=" + repr(action) + "\n" + REMOTE
    compile(source, f"<study2-checker-v2-{action}>", "exec")
    return json.loads(transport.python(source))


def stage_and_launch() -> dict:
    transport_spec = importlib.util.spec_from_file_location("ssh_transport", TRANSPORT / "ssh_transport.py")
    if transport_spec is None or transport_spec.loader is None:
        raise RuntimeError("cannot load established Tunya transport")
    transport = importlib.util.module_from_spec(transport_spec)
    transport_spec.loader.exec_module(transport)
    initial = remote_run("init")
    assert initial["status"] == "fresh_stage_created"
    remote_payload = STAGE + "/study2-checker-v2-allowlist.tar.gz"
    remote_manifest = STAGE + "/study2-checker-v2-allowlist.json"
    transport.send([ARCHIVE, MANIFEST], STAGE)
    extracted = remote_run("extract")
    assert extracted["status"] == "staged" and extracted["members_verified"] == 3
    launched = remote_run("launch")
    launched["stage_receipt"] = initial
    launched["extract_receipt"] = extracted
    launched["operator_sha256"] = sha(Path(__file__))
    with LOCAL_RECEIPT.open("x", encoding="utf-8") as stream:
        json.dump(launched, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return launched


if __name__ == "__main__":
    print(json.dumps(stage_and_launch(), indent=2, sort_keys=True))
