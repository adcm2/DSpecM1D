from pathlib import Path
import json,subprocess,os,resource,gzip,shutil,numpy as np
R=Path(__file__).resolve().parent;root=R/'field-validation-v4';root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(120,120));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
records=[]
for backend in [0]:
 for bridge in [0,1]:
  out=root/f'backend{backend}-bridge{bridge}';out.mkdir();cmd=[str(R/'bin/fields-validation-v4'),str(out/'contacts.csv'),str(backend)]
  with (out/'fields.csv').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
   p=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=env|{'DSFAST_RSTRAIN_BRIDGE':str(bridge)},preexec_fn=guard,timeout=150)
  if p.returncode:raise RuntimeError((out/'stderr.log').read_text()[-3000:])
  pack=[json.loads(x[len('PACKERR '):]) for x in (out/'stderr.log').read_text().splitlines() if x.startswith('PACKERR ')]
  maxima=[max(x['relative_norms'][i] for x in pack) for i in range(3)]
  for name in ['fields.csv','contacts.csv']:
   src=out/name
   with src.open('rb') as f,gzip.open(str(src)+'.gz','wb',compresslevel=3) as g:shutil.copyfileobj(f,g)
   src.unlink()
  r={'backend':backend,'bridge':bridge,'command':cmd,'packet_groups':len(pack),'max_internal_packed_norms_displacement_forwardstrain_adjointstrain':maxima,'internal_packed_passed':all(a<=b for a,b in zip(maxima,[1e-12,1e-12,1e-9])),'output':str(out)}
  records.append(r);(root/'records.json').write_text(json.dumps(records,indent=2)+'\n');print(backend,bridge,maxima,flush=True)
print('Full-field arrays retained for frozen-parent comparison',flush=True)
