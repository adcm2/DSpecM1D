from pathlib import Path
import json,subprocess,os,resource,time
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(120,120));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
root=R/'frozen-reference-cases';root.mkdir(exist_ok=True)
cases=[(3,1,3,'direct'),(3,1,3,'native'),(12,3,24,'direct'),(12,3,1352,'direct'),(12,3,24,'native'),(12,12,192,'direct'),(12,48,192,'direct'),(6,3,1352,'direct')]
manifest=[]
for l,nf,a,rad in cases:
 out=root/f'l{l}-f{nf}-a{a}-{rad}'
 if out.exists():raise RuntimeError('Refusing reference overwrite '+str(out))
 out.mkdir();data=out/'data';data.mkdir()
 cmd=list(map(str,[R/'bin/frozen-reference',repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',R/'sample-grids'/f'l{l}-a{a}.csv',data,l,nf,1,'openblas',rad]))
 with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
  p=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=env,preexec_fn=guard,timeout=150)
 if p.returncode:raise RuntimeError((out/'stderr.log').read_text())
 meta=json.loads((data/'result.json').read_text());manifest.append({'degree':l,'bins':nf,'angles':a,'radial':rad,'command':cmd,'data':str(data),'metrics':meta})
 (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(out.name,meta,flush=True)
