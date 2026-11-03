from pathlib import Path
import json,os,subprocess,resource,hashlib
from validate_arrays import compare
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
root=R/'final-validation-v4';root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
exe=R/'fast/fast_driver-v4';sha=hashlib.sha256(exe.read_bytes()).hexdigest();records=[]
cases=[(3,1,192,'direct',1,0,0,'physical','tensor-real'),(3,3,192,'direct',2,1,1,'physical','tensor-real'),(3,3,192,'direct',3,1,1,'physical','tensor-real'),(6,3,1352,'direct',1,0,1,'physical','tensor-real'),(12,3,1352,'direct',1,0,1,'physical','tensor-real'),(12,12,192,'direct',1,0,1,'physical','tensor-real'),(12,48,192,'direct',1,0,1,'physical','tensor-real'),(12,3,1352,'direct',2,1,1,'physical','tensor-real'),(12,3,24,'direct',1,0,1,'physical','openblas'),(12,3,24,'direct',1,0,1,'physical','scalar'),(12,3,24,'direct',1,0,1,'physical','eigen'),(3,1,192,'native',1,0,1,'physical','tensor-real'),(12,3,24,'native',1,0,1,'physical','tensor-real'),(12,3,24,'native',1,0,1,'reference','tensor-real'),(12,3,24,'direct',1,0,1,'physical','tensor')]
for l,nf,a,rad,nobj,batch,contacts,deriv,backend in cases:
 name=f'l{l}-f{nf}-a{a}-{rad}-o{nobj}-batch{batch}-c{contacts}-{deriv}-{backend}';out=root/name;out.mkdir();ref=out/'reference';ref.mkdir();cand=out/'candidate';cand.mkdir()
 e=env|{'DSFAST_BATCH_RESPONSES':str(batch),'DSFAST_CONTACTS':str(contacts),'DSFAST_NATIVE_DERIVATIVE':deriv,'DSFAST_EXPORT_PACKET':'0'}
 base=list(map(str,[repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',R/'sample-grids'/f'l{l}-a{a}.csv']))
 cmds=[[str(R/'bin/frozen-reference-v2')]+base+list(map(str,[ref,l,nf,nobj,backend,rad])),[str(exe)]+base+list(map(str,[cand,l,nf,nobj,backend,rad,'equilibrated']))]
 for i,cmd in enumerate(cmds):
  with (out/f'stdout{i}.log').open('w') as stdout,(out/f'stderr{i}.log').open('w') as stderr:
   p=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=e,preexec_fn=guard,timeout=210)
  if p.returncode:raise RuntimeError(str(out)+' '+(out/f'stderr{i}.log').read_text()[-2000:])
 assert (ref/'coordinates.bin').read_bytes()==(cand/'coordinates.bin').read_bytes()
 c=compare(cand/'arrays_candidate.bin',ref/'arrays_reference.bin');meta=json.loads((cand/'result.json').read_text())
 rec={'name':name,'degree':l,'bins':nf,'angles':a,'radial':rad,'objectives':nobj,'batch':batch,'contacts':contacts,'derivative':deriv,'backend':backend,'commands':cmds,'candidate_sha256':sha,'comparison':c,'candidate_metrics':meta,'reference_metrics':json.loads((ref/'result.json').read_text())}
 records.append(rec);(out/'comparison.json').write_text(json.dumps(c,indent=2)+'\n');(root/'records.json').write_text(json.dumps(records,indent=2)+'\n')
 failed=[(r['objective'],k,v['relative_L2']) for r in c['objectives'] for k,v in r['checks'].items() if not v['passed']]
 print(name,c['passed'],failed,meta['candidate_s'],flush=True)
assert hashlib.sha256(exe.read_bytes()).hexdigest()==sha
print('Frozen binary unchanged',flush=True)
