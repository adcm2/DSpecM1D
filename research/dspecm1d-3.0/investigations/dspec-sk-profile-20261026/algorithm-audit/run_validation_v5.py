from pathlib import Path
import json,subprocess,os,resource,hashlib
from validate_arrays import compare
R=Path(__file__).resolve().parent;root=R/'final-validation-v5';root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
exe=R/'fast/fast_driver-v5';sha=hashlib.sha256(exe.read_bytes()).hexdigest();records=[]
for prev in json.loads((R/'final-validation-v4/records.json').read_text()):
 out=root/prev['name'];out.mkdir();data=out/'data';data.mkdir()
 cmd=prev['commands'][1].copy();cmd[0]=str(exe);cmd[3]=str(data)
 e=env|{'DSFAST_BATCH_RESPONSES':str(prev['batch']),'DSFAST_CONTACTS':str(prev['contacts']),'DSFAST_NATIVE_DERIVATIVE':prev['derivative'],'DSFAST_EXPORT_PACKET':'0'}
 with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:p=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=e,preexec_fn=guard,timeout=210)
 if p.returncode:raise RuntimeError(str(out)+(out/'stderr.log').read_text()[-3000:])
 ref=Path(prev['commands'][0][3]);assert (ref/'coordinates.bin').read_bytes()==(data/'coordinates.bin').read_bytes()
 c=compare(data/'arrays_candidate.bin',ref/'arrays_reference.bin');meta=json.loads((data/'result.json').read_text())
 r={'name':prev['name'],'command':cmd,'reference':str(ref),'candidate_sha256':sha,'configuration':{k:prev[k] for k in ['degree','bins','angles','radial','objectives','batch','contacts','derivative','backend']},'comparison':c,'metrics':meta,'data':str(data)};records.append(r)
 (out/'comparison.json').write_text(json.dumps(c,indent=2)+'\n');(root/'records.json').write_text(json.dumps(records,indent=2)+'\n')
 print(prev['name'],c['passed'],meta['candidate_s'],'cache',meta['tensor_grid_hits'],meta['tensor_plan_hits'],meta['tensor_basis_hits'],flush=True)
assert hashlib.sha256(exe.read_bytes()).hexdigest()==sha
print('Frozen v5 unchanged',flush=True)
