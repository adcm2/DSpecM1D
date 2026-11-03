#!/usr/bin/env python3
"""Exclusive, hash-checked Study 1 staging and detached producer launch."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
ACC = ROOT.parent / 'dspec-acceleration-20261027'
TRANSPORT = ROOT.parent / 'tunya-profile-20261027'
BASE = '/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya'
PRODUCER_STAGE = BASE + '/scientific-inputs/dspec-science-study1-v1-validation'
CHECKER_STAGE = BASE + '/scientific-inputs/dspec-science-study1-v1-check'
BUILD_STAGE = BASE + '/scientific-inputs/dspec-science-study1-v1'
BUILD = '/space/adcm2/research/builds/dspec-science-study1-v1'
BUILD_CONTROL = BASE + '/dspec-science-study1-v1-build-control'
RUN_CONTROL = BASE + '/dspec-science-study1-v1-control'
CHECK_CONTROL = BASE + '/dspec-science-study1-v1-check-control'
RUN_TASK = 'dspec-science-study1-v1-run-v1'
CHECK_TASK = 'dspec-science-study1-v1-check-v1'
RUN_PLAN = ROOT / 'configs/study1-run-v1-plan.json'
CHECK_PLAN = ROOT / 'configs/study1-terminal-check-v1-plan.json'
BUILD_PLAN = ROOT / 'configs/study1-build-v1-plan.json'
RECEIPTS = ROOT / 'receipts'
RUN_ARCHIVE = RECEIPTS / 'study1-run-v1-validation-allowlist.tar.gz'
RUN_MANIFEST = RECEIPTS / 'study1-run-v1-validation-allowlist.json'
CHECK_ARCHIVE = RECEIPTS / 'study1-terminal-check-v1-allowlist.tar.gz'
CHECK_MANIFEST = RECEIPTS / 'study1-terminal-check-v1-allowlist.json'
SUPERVISOR = BASE + '/continuation-inputs/v6/tools/supervise_campaign.py'
BUILD_TOOL = BASE + '/continuation-inputs/v3/tools/build_campaign.py'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_payload():
    run_manifest = json.loads(RUN_MANIFEST.read_text())
    check_manifest = json.loads(CHECK_MANIFEST.read_text())
    build_stage = json.loads((RECEIPTS / 'study1-build-v1-stage.json').read_text())
    build_evidence = json.loads((RECEIPTS / 'study1-build-v1-evidence.json').read_text())
    assert run_manifest['archive_sha256'] == sha(RUN_ARCHIVE)
    assert check_manifest['archive_sha256'] == sha(CHECK_ARCHIVE)
    for manifest, archive in ((run_manifest, RUN_ARCHIVE), (check_manifest, CHECK_ARCHIVE)):
        with __import__('tarfile').open(archive) as tar:
            members = tar.getmembers()
            assert len(members) == manifest['file_count']
            assert {m.name for m in members} == set(manifest['files'])
            for m in members:
                assert m.isfile() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts
                assert m.size == manifest['files'][m.name]['bytes']
                data = tar.extractfile(m).read()
                assert hashlib.sha256(data).hexdigest() == manifest['files'][m.name]['sha256']
    assert sha(RUN_PLAN) == 'a1bf9d07ed56799537db489432c42f57b58287d24d5ee537eb94e91c20219e62'
    assert sha(CHECK_PLAN) == 'd4f3cc3c203dddb69d43d1c78d65e31af1a66309f0c79c0223f4ccf478971201'
    assert sha(BUILD_PLAN) == '31eb960a56cb3cadb1f11e90e864848b8d1f99c2ed90e85c4ac3c6c2a9fc9ef2'
    transport_root = TRANSPORT
    tools = {
        SUPERVISOR: sha(transport_root / 'tools-v4/supervise_campaign.py'),
        BUILD_TOOL: sha(transport_root / 'tools/build_campaign.py'),
        BASE + '/continuation-inputs/v3/build_receipt.py': sha(transport_root / 'build_receipt.py'),
    }
    return {
        'base': BASE, 'producer_stage': PRODUCER_STAGE, 'checker_stage': CHECKER_STAGE,
        'build_stage': BUILD_STAGE, 'build': BUILD, 'build_control': BUILD_CONTROL,
        'run_control': RUN_CONTROL, 'check_control': CHECK_CONTROL,
        'run_task': RUN_TASK, 'check_task': CHECK_TASK,
        'supervisor': SUPERVISOR, 'build_tool': BUILD_TOOL,
        'tools': tools,
        'run_archive_sha256': sha(RUN_ARCHIVE), 'run_manifest_sha256': sha(RUN_MANIFEST),
        'run_archive_name': RUN_ARCHIVE.name, 'run_manifest_name': RUN_MANIFEST.name,
        'run_manifest': run_manifest,
        'check_archive_sha256': sha(CHECK_ARCHIVE), 'check_manifest_sha256': sha(CHECK_MANIFEST),
        'check_archive_name': CHECK_ARCHIVE.name, 'check_manifest_name': CHECK_MANIFEST.name,
        'check_manifest': check_manifest,
        'build_evidence': build_evidence, 'build_stage_receipt': build_stage,
        'run_plan_sha256': sha(RUN_PLAN), 'check_plan_sha256': sha(CHECK_PLAN),
        'build_plan_sha256': sha(BUILD_PLAN),
        'operator_sha256': sha(__file__),
    }


REMOTE = r'''
from pathlib import Path
import hashlib, importlib.util, json, os, subprocess, tarfile, time
CFG=json.loads(PAYLOAD)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def tree_hash(path):
    spec=importlib.util.spec_from_file_location('frozen_build_tool',CFG['build_tool'])
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod.tree_hash(Path(path))
def check_tools():
    for p,h in CFG['tools'].items(): assert sha(p)==h,('tool hash changed',p)
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
    assert not active,('competing Study 1 task',active)
    return active
def verify_build():
    assert Path(CFG['build']).is_dir() and Path(CFG['build_control']).is_dir()
    for key in ('dspec-timer-disabled','dspec-timer-enabled'):
        variant=CFG['build_evidence']['command_receipts'][key]
        root=Path(CFG['build'])/key
        binary=root/'phase-d-driver'; cmd=root/'command.json'
        assert binary.is_file() and cmd.is_file()
        assert sha(binary)==variant['binary_sha256']==variant['binary_sha256_observed']
        record=json.loads(cmd.read_text())
        assert record.get('returncode')==0 and record.get('binary_sha256')==variant['binary_sha256']
        assert record.get('command_sha256')==variant['command_sha256']
        assert record.get('tracked_inputs_unchanged') is True
        assert record.get('tracked_input_sha256')==variant['tracked_input_sha256']
        assert record.get('tracked_input_sha256_after')==variant['tracked_input_sha256']
    expected=CFG['build_stage_receipt']['external_dependencies']
    observed={p:tree_hash(p) for p in expected}
    assert observed==expected,'external dependency tree changed since reviewed Study 1 build'
    for p,h in CFG['build_stage_receipt']['tools'].items(): assert sha(p)==h,p
    return {'binaries':{k:v['binary_sha256'] for k,v in CFG['build_evidence']['command_receipts'].items()},
            'external_dependencies':observed}
def verify_plan_inputs(plan, allowed_missing=()):
    missing=set()
    for j in plan['jobs']:
        for name in j.get('inputs',[]):
            if '{JOB_DIR}' in name: continue
            if not Path(name).is_file(): missing.add(name)
    assert set(missing)==set(allowed_missing),('missing input set differs',missing,allowed_missing)
    return sorted(missing)
def extract_exact(stage, archive_name, manifest_name, manifest):
    stage=Path(stage)
    assert {p.name for p in stage.iterdir()}=={archive_name,manifest_name}
    assert sha(stage/archive_name)==manifest['archive_sha256']
    with tarfile.open(stage/archive_name) as tar:
        ms=tar.getmembers(); assert len(ms)==manifest['file_count'] and {m.name for m in ms}==set(manifest['files'])
        for m in ms:
            rel=Path(m.name); assert m.isfile() and not rel.is_absolute() and '..' not in rel.parts
            assert m.size==manifest['files'][m.name]['bytes']
            dest=stage/rel;dest.parent.mkdir(parents=True,exist_ok=True)
            with dest.open('xb') as out: out.write(tar.extractfile(m).read())
    for rel,item in manifest['files'].items(): assert sha(stage/rel)==item['sha256'],rel
def report_preflight():
    check_tools();quiet();build=verify_build()
    roots=[CFG[k] for k in ('producer_stage','checker_stage','run_control','check_control')]
    assert all(not Path(p).exists() for p in roots),('nonfresh roots',[(p,Path(p).exists()) for p in roots])
    # Before staging, validate the reviewed producer-plan object and classify
    # only the exact producer-stage package members as expected missing files.
    plan=CFG['run_plan_object']
    assert plan['thread_cap']==20 and len(plan['jobs'])==10
    generated=[]
    for j in plan['jobs']:
        for out in j['outputs']: generated.append(j['id']+'/'+out)
    assert len(generated)==len(set(generated))
    stage=Path(CFG['producer_stage'])
    referenced_stage_paths=set()
    for j in plan['jobs']:
        for name in j.get('inputs',[]):
            path=Path(name)
            if path.is_relative_to(stage):
                rel=str(path.relative_to(stage))
                assert rel in CFG['run_manifest']['files'],('producer command input absent from allowlist',rel)
                referenced_stage_paths.add(name)
    missing=verify_plan_inputs(plan,allowed_missing=referenced_stage_paths)
    affinity=sorted(os.sched_getaffinity(0))
    mem={}
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'): mem['available_bytes']=int(line.split()[1])*1024
    assert mem.get('available_bytes',0)>=16*(1024**3),('insufficient available memory for one bounded job',mem)
    assert set(range(16,24)).issubset(affinity),('planned taskset CPUs unavailable',affinity)
    for job in plan['jobs']:
        cmd=job['command']
        assert 'OMP_NUM_THREADS=1' in cmd and 'OMP_THREAD_LIMIT=8' in cmd
        if 'yspec' in cmd:
            assert cmd[:3]==['/usr/bin/env','-u','DSPEC_FREQUENCY_WORKERS']
        else:
            assert 'DSPEC_FREQUENCY_WORKERS=8' in cmd
    fs=os.statvfs(CFG['base'])
    return {'status':'pass','fresh_roots':roots,'binary_dependencies':build,'quiet_processes':[],
            'affinity':affinity,'meminfo':mem,'free_bytes':fs.f_bavail*fs.f_frsize,
            'run_plan_sha256':CFG['run_plan_sha256'],'terminal_plan_sha256':CFG['check_plan_sha256'],
            'producer_inputs_expected_from_allowlist':missing,'declared_outputs':len(generated)}
def stage():
    check_tools();quiet();build=verify_build()
    for name,archive,manifest in (
        ('producer_stage',CFG['run_archive_name'],CFG['run_manifest_name']),
        ('checker_stage',CFG['check_archive_name'],CFG['check_manifest_name'])):
        path=Path(CFG[name]);assert path.is_dir(),('stage root missing',str(path))
        assert {p.name for p in path.iterdir()}=={archive,manifest},('unexpected staged payload',str(path))
    # Payload files are copied by the caller into each exclusively created stage.
    extract_exact(CFG['producer_stage'],CFG['run_archive_name'],CFG['run_manifest_name'],CFG['run_manifest'])
    extract_exact(CFG['checker_stage'],CFG['check_archive_name'],CFG['check_manifest_name'],CFG['check_manifest'])
    run=json.loads((Path(CFG['producer_stage'])/'configs/study1-run-v1-plan.json').read_text())
    term=json.loads((Path(CFG['checker_stage'])/'configs/study1-terminal-check-v1-plan.json').read_text())
    assert sha(Path(CFG['producer_stage'])/'configs/study1-run-v1-plan.json')==CFG['run_plan_sha256']
    assert sha(Path(CFG['checker_stage'])/'configs/study1-terminal-check-v1-plan.json')==CFG['check_plan_sha256']
    assert run['campaign_root']==CFG['run_control'] and term['campaign_root']==CFG['check_control']
    missing=verify_plan_inputs(run)
    outputs={(j['id'],o) for j in run['jobs'] for o in j['outputs']}
    runroot=Path(CFG['run_control'])/'continuation'/CFG['run_task']
    snapshot=Path(CFG['checker_stage'])/'data/producer-state-finalized.json'
    allowed_missing={str(snapshot)}
    for j in run['jobs']:
        for o in j['outputs']: allowed_missing.add(str(runroot/j['id']/o))
        allowed_missing.add(str(runroot/j['id']/'job.log'))
    term_missing=verify_plan_inputs(term,allowed_missing=allowed_missing)
    return {'status':'staged','producer_member_count':len(CFG['run_manifest']['files']),
            'checker_member_count':len(CFG['check_manifest']['files']),'build':build,
            'plan_hashes':{'producer':CFG['run_plan_sha256'],'terminal':CFG['check_plan_sha256']},
            'producer_declared_inputs_missing':missing,'checker_declared_generated_inputs_missing':term_missing}
def launch():
    check_tools();quiet();build=verify_build()
    stage=Path(CFG['producer_stage']);plan=stage/'configs/study1-run-v1-plan.json'
    assert sha(plan)==CFG['run_plan_sha256']
    parsed=json.loads(plan.read_text());assert parsed['campaign_root']==CFG['run_control'] and parsed['thread_cap']==20 and len(parsed['jobs'])==10
    control=Path(CFG['run_control'])
    continuation=control/'continuation'
    assert control.is_dir() and {p.name for p in control.iterdir()}=={'continuation'}
    assert continuation.is_dir() and not any(continuation.iterdir())
    missing=verify_plan_inputs(parsed);assert not missing,missing
    for m in (CFG['run_manifest'],CFG['check_manifest']):
        # Recheck staged members against the exact locally reviewed hash manifests.
        target=Path(CFG['producer_stage'] if m is CFG['run_manifest'] else CFG['checker_stage'])
        for rel,item in m['files'].items(): assert sha(target/rel)==item['sha256'],rel
    assert continuation.is_dir() and not any(continuation.iterdir()),'partial launch directory is not empty'
    taskroot=continuation/CFG['run_task']
    log=continuation/(CFG['run_task']+'.nohup.log')
    launch_record=continuation/(CFG['run_task']+'.launch.json')
    assert not taskroot.exists() and not log.exists() and not launch_record.exists()
    env=os.environ.copy();env.update(LD_LIBRARY_PATH='/opt/gcc-13.2.0/lib64',LD_PRELOAD='',PYTHONDONTWRITEBYTECODE='1')
    for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','GOTO_NUM_THREADS'):
        env[k]='1'
    with log.open('xb') as out:
        proc=subprocess.Popen(['/usr/bin/nohup','/usr/bin/python3',CFG['supervisor'],str(plan),'--task-id',CFG['run_task']],
            cwd=stage,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
    receipt={'pid':proc.pid,'task_id':CFG['run_task'],'plan_sha256':sha(plan),'log':str(log),
        'stage':str(stage),'control_root':CFG['run_control'],'thread_cap':20,'memory_gib_per_job':16,
        'job_count':10,'run_archive_sha256':CFG['run_archive_sha256'],'run_manifest_sha256':CFG['run_manifest_sha256'],
        'checker_archive_sha256':CFG['check_archive_sha256'],'checker_manifest_sha256':CFG['check_manifest_sha256'],
        'build_binary_hashes':build['binaries'],'external_dependencies':build['external_dependencies'],
        'operator_sha256':CFG['operator_sha256'],'launched_unix':time.time()}
    with launch_record.open('x') as f: json.dump(receipt,f,indent=2);f.write('\n')
    return receipt
if ACTION=='preflight': RESULT=report_preflight()
elif ACTION=='stage': RESULT=stage()
elif ACTION=='launch': RESULT=launch()
else: raise ValueError(ACTION)
print(json.dumps(RESULT,sort_keys=True))
'''


def run(action):
    spec = importlib.util.spec_from_file_location('ssh_transport', TRANSPORT / 'ssh_transport.py')
    transport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(transport)
    payload = load_payload()
    payload['run_plan'] = str(RUN_PLAN)
    payload['run_plan_object'] = json.loads(RUN_PLAN.read_text())
    payload['action'] = action
    source = 'import json\nPAYLOAD=' + repr(json.dumps(payload)) + '\nACTION=' + repr(action) + '\n' + REMOTE
    compile(source, '<study1-remote-' + action + '>', 'exec')
    if action == 'stage':
        # The remote roots are created exclusively before any bytes are transferred.
        create = ('from pathlib import Path\n'
                  'Path(' + repr(PRODUCER_STAGE) + ').mkdir(parents=True,exist_ok=False)\n'
                  'Path(' + repr(CHECKER_STAGE) + ').mkdir(parents=True,exist_ok=False)\n')
        compile(create, '<study1-stage-root-creation>', 'exec')
        transport.python(create)
        transport.send([RUN_ARCHIVE, RUN_MANIFEST], PRODUCER_STAGE + '/')
        transport.send([CHECK_ARCHIVE, CHECK_MANIFEST], CHECKER_STAGE + '/')
    raw = transport.python(source)
    result = json.loads(raw)
    result['remote_script_sha256'] = hashlib.sha256(source.encode()).hexdigest()
    result['operator_sha256'] = sha(__file__)
    output = RECEIPTS / ('study1-run-v1-' + action + '.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    return result


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('preflight', 'stage', 'launch'):
        raise SystemExit('usage: study1_run_stage.py preflight|stage|launch')
    print(json.dumps(run(sys.argv[1]), indent=2, sort_keys=True))
