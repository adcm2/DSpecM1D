#!/usr/bin/env python3
"""Phase B operator: exclusive hashed staging with sequential in-plan dependencies."""
from pathlib import Path
import argparse,hashlib,json,sys,tarfile
A=Path(__file__).resolve().parents[1]
N=A.parent/'tunya-profile-20261027'
sys.path.insert(0,str(N));import ssh_transport as st
p=argparse.ArgumentParser();p.add_argument('stage_dir');p.add_argument('plan_name');p.add_argument('task');p.add_argument('--frozen-trace');p.add_argument('--reference-components');p.add_argument('--reference-products');args=p.parse_args()
C='/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya';S=C+'/acceleration-inputs/'+args.stage_dir
files=[A/'configs'/args.plan_name]+sorted((A/'tools').glob('check*.py'))+sorted((A/'tools').glob('test_check*.py'))+sorted((A/'review').glob('*.md'))
if args.frozen_trace: files.append(Path(args.frozen_trace).resolve())
if args.reference_components: files.append(Path(args.reference_components).resolve())
if args.reference_products: files.append(Path(args.reference_products).resolve())
hashfile=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
manifest={str(x.relative_to(A)):hashfile(x) for x in files};archive=A/f'{args.stage_dir}-inputs.tar.gz'
with tarfile.open(archive,'x:gz') as tf:
 for x in files:tf.add(x,arcname=str(x.relative_to(A)),recursive=False)
assert manifest=={str(x.relative_to(A)):hashfile(x) for x in files}
m={'files':manifest,'archive_sha256':hashfile(archive)}
with (A/f'{args.stage_dir}-inputs.json').open('x') as f:json.dump(m,f,indent=2)
st.python('from pathlib import Path\nPath('+repr(S)+').mkdir(parents=True)')
st.send([archive],S+'/')
checks={'continuation-inputs/v6/tools/supervise_campaign.py':hashfile(N/'tools-v4/supervise_campaign.py'),'continuation-inputs/v3/build_receipt.py':hashfile(N/'build_receipt.py'),'continuation-inputs/v3/tools/build_campaign.py':hashfile(N/'tools/build_campaign.py')}
source='''import pathlib,hashlib,json,tarfile,subprocess,os
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
plan=S/'configs'/PLANFILE;control=pathlib.Path(json.loads(plan.read_text())['campaign_root']);(control/'continuation').mkdir(parents=True,exist_ok=True)
parsed=json.loads(plan.read_text())
produced=set()
for external in parsed.get('external_campaigns',[]):
 ep=pathlib.Path(external['plan_path']);assert hashlib.sha256(ep.read_bytes()).hexdigest()==external['plan_sha256']
 ed=json.loads(ep.read_text());er=pathlib.Path(external['task_path'])
 for ej in ed['jobs']:
  for name in ej.get('outputs',[]):
   rel=pathlib.Path(name);assert not rel.is_absolute() and '..' not in rel.parts
   produced.add(str(er/ej['id']/rel))
for job in parsed['jobs']:
 for name in job.get('inputs',[]):
  path=pathlib.Path(name);assert path.exists() or str(path) in produced,name
  if path.name=='physical-driver':
   command=json.loads((path.parent/'command.json').read_text());assert command['returncode']==0 and command['tracked_inputs_unchanged']
   assert hashlib.sha256(path.read_bytes()).hexdigest()==command['binary_sha256']
 for name in job.get('outputs',[]):
  rel=pathlib.Path(name);assert not rel.is_absolute() and '..' not in rel.parts
  produced.add(str(control/'continuation'/TASK/job['id']/rel))
task=TASK;log=control/'continuation'/f'{task}.nohup.log'
assert not (control/'continuation'/task).exists()
env=os.environ.copy();env.update(LD_LIBRARY_PATH='/opt/gcc-13.2.0/lib64',LD_PRELOAD='',PYTHONDONTWRITEBYTECODE='1')
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
with log.open('xb') as f:
 proc=subprocess.Popen(['/usr/bin/nohup','/usr/bin/python3',str(supervisor),str(plan),'--task-id',task],cwd=S,env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
receipt={'pid':proc.pid,'task_id':task,'control_root':str(control),'log':str(log),'plan_sha256':hashlib.sha256(plan.read_bytes()).hexdigest(),'staged_manifest':M,'thread_cap':json.loads(plan.read_text())['thread_cap']}
with (control/'continuation'/f'{task}.launch.json').open('x') as f:json.dump(receipt,f,indent=2)
print(json.dumps(receipt))
'''
for key,value in {'SPATH':S,'CPATH':C,'MANIFEST':m,'ARCHIVE':archive.name,'TOOLS':checks,'PLANFILE':args.plan_name,'TASK':args.task}.items():source=source.replace(key,repr(value))
raw=st.python(source)
with (A/'receipts'/f'{args.task}.launch.json').open('x') as f:f.write(raw)
d=json.loads(raw);print({k:d[k] for k in ('pid','task_id','control_root','plan_sha256','thread_cap')})
print(st.python('import subprocess,pathlib\nprint(subprocess.check_output(["ps","-o","pid,ppid,sid,nlwp,etime,stat,args","-p",'+repr(str(d['pid']))+'],text=True))\np=pathlib.Path('+repr(d['control_root']+'/continuation/'+d['task_id']+'/state.json')+');print(p.read_text()[:300] if p.exists() else "state pending")'))
