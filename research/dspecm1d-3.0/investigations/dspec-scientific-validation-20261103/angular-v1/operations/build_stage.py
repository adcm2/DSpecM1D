#!/usr/bin/env python3
"""Minimal, exclusive staging of the independently reviewed angular-convergence build.

Network actions use the established SSH transport; invoke them with the normal
external-access approval mechanism. This does not change SSH configuration.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import sys
import tarfile
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
INVESTIGATIONS = ROOT.parent.parent
BASE = '/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya'
STAGE = BASE + '/scientific-inputs/dspec-science-angular-v1'
BUILD = '/space/adcm2/research/builds/dspec-science-angular-v1'
CONTROL = BASE + '/dspec-science-angular-v1-build-control'
TASK = 'dspec-science-angular-v1-build'
RECEIPTS = ROOT / 'receipts'
ARCHIVE = 'angular-build-v1-allowlist.tar.gz'
MANIFEST = 'angular-build-v1-allowlist.json'
SUPERVISOR = BASE + '/continuation-inputs/v6/tools/supervise_campaign.py'
BUILD_TOOL = BASE + '/continuation-inputs/v3/tools/build_campaign.py'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')


def log(message):
    with (ROOT / 'operations/implementation_status.md').open('a') as stream:
        stream.write('\n2026-11-03: ' + message + '\n')


def common_payload():
    resource_spec=importlib.util.spec_from_file_location('angular_resource',
        ROOT/'operations/resource_admission.py')
    resource=importlib.util.module_from_spec(resource_spec)
    resource_spec.loader.exec_module(resource)
    old = json.loads((INVESTIGATIONS / 'dspec-acceleration-20261027/receipts/'
                     'phase-e-stage2-modal-v1-rebuild3-build-evidence-20261103.json').read_text())
    tracked = old['command_receipts']['dspec-timer-disabled']['tracked_input_sha256']
    deps = {p: h for p, h in tracked.items()
            if '/acceleration-inputs/phase-e-stage2-modal-v1-rebuild3/' not in p}
    assert len(deps) == 12, 'unexpected frozen external dependency inventory'
    transport = INVESTIGATIONS / 'tunya-profile-20261027'
    tools = {
        SUPERVISOR: sha(transport / 'tools-v4/supervise_campaign.py'),
        BUILD_TOOL: sha(transport / 'tools/build_campaign.py'),
        BASE + '/continuation-inputs/v3/build_receipt.py': sha(transport / 'build_receipt.py'),
    }
    return dict(stage=STAGE, build=BUILD, control=CONTROL, task=TASK,
                supervisor=SUPERVISOR, build_tool=BUILD_TOOL, tools=tools,
                dependencies=deps, operator_sha256=sha(__file__),
                **resource.resource_payload())


def package(review_name):
    review = ROOT / review_name
    review_text = review.read_text()
    required = ['candidate/dspec/phase_d_driver.cpp', 'candidate/dspec/phase_d_controls.hpp',
                'candidate/dspec/phase_d_batched_reconstruction.hpp', 'configs/angular-build-v1.json',
                'configs/angular-build-v1-plan.json']
    assert 'PASS' in review_text
    for rel in required:
        assert sha(ROOT / rel) in review_text, 'review does not seal current ' + rel
    candidate = (ROOT / 'candidate/dspec').resolve()
    queue = [candidate / 'phase_d_driver.cpp']
    closure = set()
    while queue:
        path = queue.pop()
        if path in closure:
            continue
        assert path.is_relative_to(candidate)
        closure.add(path)
        for name in re.findall(r'^\s*#\s*include\s*"([^"]+)"', path.read_text(), re.M):
            matches = [(path.parent / name).resolve(), (candidate / name).resolve()]
            found = next((p for p in matches if p.is_file() and p.is_relative_to(candidate)), None)
            if found is not None:
                queue.append(found)
    files = closure | {ROOT / name for name in required[1:]}
    files |= {candidate / 'case/kernels-small.json', candidate / 'case/prem.200', review}
    files = sorted(files, key=lambda p: str(p.relative_to(ROOT)))
    members = {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size}
               for p in files}
    with tarfile.open(RECEIPTS / ARCHIVE, 'x:gz') as archive:
        for path in files:
            archive.add(path, arcname=str(path.relative_to(ROOT)), recursive=False)
    assert members == {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size}
                       for p in files}, 'source changed while packaging'
    manifest = dict(schema=1, purpose='minimal reviewed angular-convergence solver build',
                    intended_remote_stage_root=STAGE, intended_build_root=BUILD,
                    file_count=len(files), files=members,
                    archive_sha256=sha(RECEIPTS / ARCHIVE), review=review_name,
                    review_sha256=sha(review), operator_sha256=sha(__file__),
                    resource_helper_sha256=common_payload()['resource_helper_sha256'],
                    resource_admission_receipt_sha256=common_payload()['resource_admission_receipt_sha256'])
    write_json(RECEIPTS / MANIFEST, manifest)
    log('Prepared exact reviewed build allowlist, ' + str(len(files)) + ' files; archive SHA-256 '
        + manifest['archive_sha256'] + '. No remote action.')
    return manifest


REMOTE_COMMON = '''
from pathlib import Path
import hashlib, importlib.util, json, os, subprocess, time
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def quiet():
    active=[]
    names={'g++','gcc','c++','cc1plus','phase-d-driver','physical-driver','doc-driver','phase-c2-driver'}
    for path in Path('/proc').iterdir():
        if not path.name.isdigit(): continue
        try:
            if path.stat().st_uid != os.getuid(): continue
            args=(path/'cmdline').read_bytes().replace(bytes([0]),b' ').decode(errors='replace')
            if not args: continue
            first=Path(args.split()[0]).name
            if first in names or ('supervise_campaign.py' in args and CFG['stage'].split('/scientific-inputs')[0] in args):
                active.append({'pid':int(path.name),'command':args[:1000]})
        except (FileNotFoundError,ProcessLookupError,PermissionError): continue
    assert not active, ('competing campaign/build process',active)
    return active
def check_external():
    for path,digest in CFG['tools'].items(): assert sha(path)==digest,path
    spec=importlib.util.spec_from_file_location('frozen_build_tool',CFG['build_tool'])
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    observed={path:mod.tree_hash(Path(path)) for path in CFG['dependencies']}
    assert observed==CFG['dependencies'],'shared dependency changed since final modal build'
    return observed
'''


def remote(action):
    spec = importlib.util.spec_from_file_location('ssh_transport',
              INVESTIGATIONS / 'tunya-profile-20261027/ssh_transport.py')
    transport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(transport)
    payload = common_payload()
    resource_spec=importlib.util.spec_from_file_location('angular_resource',
        ROOT/'operations/resource_admission.py')
    resource=importlib.util.module_from_spec(resource_spec)
    resource_spec.loader.exec_module(resource)
    if action != 'preflight':
        manifest = json.loads((RECEIPTS / MANIFEST).read_text())
        assert manifest['operator_sha256']==payload['operator_sha256'],'build operator changed after packaging'
        assert manifest['resource_helper_sha256']==payload['resource_helper_sha256'],'resource helper changed after packaging'
        assert manifest['resource_admission_receipt_sha256']==payload['resource_admission_receipt_sha256']
        assert sha(RECEIPTS / ARCHIVE) == manifest['archive_sha256']
        for rel, entry in manifest['files'].items():
            assert sha(ROOT / rel) == entry['sha256'], 'local reviewed input changed: ' + rel
        payload.update(manifest=manifest, manifest_sha256=sha(RECEIPTS / MANIFEST),
                       archive_name=ARCHIVE, manifest_name=MANIFEST)
    prefix = 'CFG=json.loads(' + repr(json.dumps(payload)) + ')\n'
    # Imports precede payload decoding; compile every remote script before SSH.
    prefix = 'import json\n' + prefix
    if action == 'preflight':
        body = '''
for name in ['stage','build','control']: assert not Path(CFG[name]).exists(),CFG[name]
deps=check_external(); quiet()
admission=resource_admit()
free=os.statvfs(Path(CFG['stage']).parent.parent)
print(json.dumps({'status':'pass','fresh_roots':{k:CFG[k] for k in ['stage','build','control']},
 'external_dependencies':deps,'tools':CFG['tools'],'free_bytes':free.f_bavail*free.f_frsize,
 'resource_admission':admission}))
'''
    elif action == 'stage':
        create = 'from pathlib import Path\nPath(' + repr(STAGE) + ').mkdir(parents=True,exist_ok=False)\n'
        compile(create, '<remote-stage-mkdir>', 'exec')
        body = '''
import tarfile
stage=Path(CFG['stage']); manifest=CFG['manifest']
assert {p.name for p in stage.iterdir()}=={CFG['archive_name'],CFG['manifest_name']}
assert sha(stage/CFG['archive_name'])==manifest['archive_sha256']
assert sha(stage/CFG['manifest_name'])==CFG['manifest_sha256']
with tarfile.open(stage/CFG['archive_name']) as archive:
    members=archive.getmembers()
    assert len(members)==len(manifest['files'])
    assert {m.name for m in members}==set(manifest['files'])
    for member in members:
        rel=Path(member.name)
        assert member.isfile() and not rel.is_absolute() and '..' not in rel.parts
        assert member.size==manifest['files'][member.name]['bytes']
        dest=stage/rel; dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('xb') as out: out.write(archive.extractfile(member).read())
for rel,item in manifest['files'].items(): assert sha(stage/rel)==item['sha256'],rel
deps=check_external(); quiet()
assert not Path(CFG['build']).exists() and not Path(CFG['control']).exists()
print(json.dumps({'status':'staged','file_count':len(manifest['files']),
 'archive_sha256':manifest['archive_sha256'],'manifest_sha256':CFG['manifest_sha256'],
 'external_dependencies':deps,'tools':CFG['tools']}))
'''
    elif action == 'launch':
        body = '''
stage=Path(CFG['stage']); manifest=CFG['manifest']
for rel,item in manifest['files'].items(): assert sha(stage/rel)==item['sha256'],rel
check_external(); quiet()
assert not Path(CFG['build']).exists() and not Path(CFG['control']).exists()
plan=stage/'configs/angular-build-v1-plan.json'; parsed=json.loads(plan.read_text())
assert parsed['campaign_root']==CFG['control'] and parsed['thread_cap']==8
assert len(parsed['jobs'])==1 and parsed['jobs'][0]['id']==CFG['task']
assert parsed['jobs'][0]['cwd']==str(stage) and parsed['jobs'][0]['memory_gib']==32
for job in parsed['jobs']:
    for path in job['inputs']: assert Path(path).exists(),path
control=Path(CFG['control']); (control/'continuation').mkdir(parents=True,exist_ok=False)
log=control/'continuation'/(CFG['task']+'.nohup.log')
with log.open('xb'): pass
resource_receipt=resource_launch(['/usr/bin/nohup','/usr/bin/python3',CFG['supervisor'],
    str(plan),'--task-id',CFG['task']],stage,log,'build')
receipt={'pid':resource_receipt['live']['main_pid'],'task_id':CFG['task'],'stage':str(stage),'control_root':str(control),
 'log':str(log),'plan_sha256':sha(plan),'archive_sha256':manifest['archive_sha256'],
 'manifest_sha256':CFG['manifest_sha256'],'operator_sha256':CFG['operator_sha256'],
 'thread_cap':8,'memory_gib':32,'launched_unix':time.time(),'resource':resource_receipt}
with (control/'continuation'/(CFG['task']+'.launch.json')).open('x') as out: json.dump(receipt,out,indent=2);out.write('\\n')
print(json.dumps(receipt))
'''
    else:
        raise ValueError(action)
    source = prefix + REMOTE_COMMON + resource.REMOTE_RESOURCE_CODE + body
    compile(source, '<reviewed-remote-' + action + '>', 'exec')
    if action == 'stage':
        transport.python(create)
        transport.send([RECEIPTS / ARCHIVE, RECEIPTS / MANIFEST], STAGE + '/')
    raw = transport.python(source)
    result = json.loads(raw)
    result['remote_script_sha256'] = hashlib.sha256(source.encode()).hexdigest()
    receipt_name=('angular-build-v1-resource-preflight.json' if action=='preflight'
                  else 'angular-build-v1-' + action + '.json')
    write_json(RECEIPTS / receipt_name, result)
    log('Build operation `' + action + '` completed; receipt saved. '
        + ('Detached build launched.' if action == 'launch' else 'No solver job launched by this operation.'))
    return result


if __name__ == '__main__':
    RECEIPTS.mkdir(exist_ok=True)
    try:
        if len(sys.argv) >= 3 and sys.argv[1] == 'package':
            result = package(sys.argv[2])
        elif len(sys.argv) == 2:
            result = remote(sys.argv[1])
        else:
            raise ValueError('usage: build_stage.py package REVIEW | preflight | stage | launch')
        print(json.dumps(result, indent=2))
    except Exception:
        failure = {'operation':sys.argv[1:], 'time_unix':time.time(),
                   'traceback':traceback.format_exc(), 'operator_sha256':sha(__file__)}
        write_json(RECEIPTS / ('build-operation-failure-' + str(time.time_ns()) + '.json'), failure)
        raise
