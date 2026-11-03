#!/usr/bin/env python3
"""Stage one reviewed Phase E plan and launch its sequential Tunya queue."""
from pathlib import Path
import argparse, hashlib, json, sys, tarfile

A = Path(__file__).resolve().parents[1]
N = A.parent / "tunya-profile-20261027"
sys.path.insert(0, str(N))
import ssh_transport as st

parser = argparse.ArgumentParser()
parser.add_argument("stage_dir")
parser.add_argument("plan_name")
parser.add_argument("task")
args = parser.parse_args()
C = "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya"
S = C + "/acceleration-inputs/" + args.stage_dir
files = [A / "configs" / args.plan_name]
files += sorted((A / "tools").glob("check*.py"))
files += sorted((A / "tools").glob("test_check*.py"))
files += sorted((A / "review").glob("*.md"))
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {str(path.relative_to(A)): digest(path) for path in files}
archive = A / f"{args.stage_dir}-inputs.tar.gz"
with tarfile.open(archive, "x:gz") as stream:
    for path in files:
        stream.add(path, arcname=str(path.relative_to(A)), recursive=False)
assert manifest == {str(path.relative_to(A)): digest(path) for path in files}
archive_record = {"files": manifest, "archive_sha256": digest(archive)}
with (A / f"{args.stage_dir}-inputs.json").open("x") as stream:
    json.dump(archive_record, stream, indent=2)

st.python("from pathlib import Path\nPath(" + repr(S) + ").mkdir(parents=True)")
st.send([archive], S + "/")
tool_hashes = {
    "continuation-inputs/v6/tools/supervise_campaign.py": digest(N / "tools-v4/supervise_campaign.py"),
    "continuation-inputs/v3/build_receipt.py": digest(N / "build_receipt.py"),
    "continuation-inputs/v3/tools/build_campaign.py": digest(N / "tools/build_campaign.py"),
}
remote_source = '''import pathlib,hashlib,json,tarfile,subprocess,os
S=pathlib.Path(SPATH);C=pathlib.Path(CPATH);M=MANIFEST
archive=S/ARCHIVE
assert hashlib.sha256(archive.read_bytes()).hexdigest()==M['archive_sha256']
with tarfile.open(archive) as tf:
 for member in tf.getmembers():
  rel=pathlib.Path(member.name);assert member.isfile() and not rel.is_absolute() and '..' not in rel.parts
  dst=S/rel;dst.parent.mkdir(parents=True,exist_ok=True)
  with dst.open('xb') as f:f.write(tf.extractfile(member).read())
for rel,digest in M['files'].items():assert hashlib.sha256((S/rel).read_bytes()).hexdigest()==digest,rel
for name,digest in TOOLS.items():assert hashlib.sha256((C/name).read_bytes()).hexdigest()==digest,name
supervisor=C/'continuation-inputs/v6/tools/supervise_campaign.py'
plan=S/'configs'/PLANFILE;parsed=json.loads(plan.read_text());control=pathlib.Path(parsed['campaign_root'])
control.joinpath('continuation').mkdir(parents=True,exist_ok=True)
produced=set()
for job in parsed['jobs']:
 for name in job.get('inputs',[]):
  path=pathlib.Path(name)
  assert path.exists() or str(path) in produced,name
 for name in job.get('outputs',[]):
  rel=pathlib.Path(name);assert not rel.is_absolute() and '..' not in rel.parts
  produced.add(str(control/'continuation'/TASK/job['id']/rel))
task=TASK;log=control/'continuation'/f'{task}.nohup.log'
assert not (control/'continuation'/task).exists()
env=os.environ.copy();env.update(LD_LIBRARY_PATH='/opt/gcc-13.2.0/lib64',LD_PRELOAD='',PYTHONDONTWRITEBYTECODE='1')
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
with log.open('xb') as f:
 proc=subprocess.Popen(['/usr/bin/nohup','/usr/bin/python3',str(supervisor),str(plan),'--task-id',task],cwd=S,env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
receipt={'pid':proc.pid,'task_id':task,'control_root':str(control),'log':str(log),'plan_sha256':hashlib.sha256(plan.read_bytes()).hexdigest(),'staged_manifest':M,'thread_cap':parsed['thread_cap']}
with (control/'continuation'/f'{task}.launch.json').open('x') as f:json.dump(receipt,f,indent=2)
print(json.dumps(receipt))
'''
for key, value in {"SPATH": S, "CPATH": C, "MANIFEST": archive_record,
                   "ARCHIVE": archive.name, "TOOLS": tool_hashes,
                   "PLANFILE": args.plan_name, "TASK": args.task}.items():
    remote_source = remote_source.replace(key, repr(value))
raw = st.python(remote_source)
with (A / "receipts" / f"{args.task}.launch.json").open("x") as stream:
    stream.write(raw)
receipt = json.loads(raw)
print({key: receipt[key] for key in ("pid", "task_id", "control_root", "plan_sha256", "thread_cap")})
print(st.python(
    "import subprocess,pathlib\nprint(subprocess.check_output(['ps','-o','pid,ppid,sid,nlwp,etime,stat,args','-p',"
    + repr(str(receipt["pid"]))
    + "],text=True))\np=pathlib.Path("
    + repr(receipt["control_root"] + "/continuation/" + receipt["task_id"] + "/state.json")
    + ");print(p.read_text()[:300] if p.exists() else 'state pending')"
))
