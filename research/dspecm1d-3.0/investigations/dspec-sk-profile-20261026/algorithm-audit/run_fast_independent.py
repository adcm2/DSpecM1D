from pathlib import Path
import json,os,subprocess,resource,hashlib,time
from validate_arrays import compare
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
root=R/'fast-vs-frozen-01';root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
exe=R/'fast/fast_driver';exe_sha=hashlib.sha256(exe.read_bytes()).hexdigest();manifest=[]
for ref in json.loads((R/'frozen-reference-cases/manifest.json').read_text()):
 for bridge in [0,1]:
  name=f'l{ref["degree"]}-f{ref["bins"]}-a{ref["angles"]}-{ref["radial"]}-bridge{bridge}'
  out=root/name;out.mkdir();data=out/'data';data.mkdir()
  cmd=list(map(str,[exe,repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',R/'sample-grids'/f'l{ref["degree"]}-a{ref["angles"]}.csv',data,ref['degree'],ref['bins'],1,'openblas',ref['radial'],'equilibrated']))
  with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
   p=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=env|{'DSFAST_RSTRAIN_BRIDGE':str(bridge)},preexec_fn=guard,timeout=210)
  if p.returncode:raise RuntimeError(str(out)+' '+(out/'stderr.log').read_text()[-2000:])
  refdata=Path(ref['data']);assert (data/'coordinates.bin').read_bytes()==(refdata/'coordinates.bin').read_bytes(),name
  result=compare(data/'arrays_candidate.bin',refdata/'arrays_reference.bin');meta=json.loads((data/'result.json').read_text())
  r={'case':name,'reference':ref,'command':cmd,'bridge':bridge,'candidate':str(data),'candidate_exe_sha_at_start':exe_sha,'metrics':meta,'frozen_comparison':result};manifest.append(r)
  (out/'independent_comparison.json').write_text(json.dumps(result,indent=2)+'\n');(root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
  print(name,'passed',result['passed'],'time',meta.get('candidate_s'),flush=True)
assert hashlib.sha256(exe.read_bytes()).hexdigest()==exe_sha,'Candidate executable changed during validation; results are preliminary'
print('Allcases',all(x['frozen_comparison']['passed'] for x in manifest),flush=True)
