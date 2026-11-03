from pathlib import Path
import json,subprocess,os,resource
from validate_arrays import compare
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D');root=R/'forward-cache-gates-v6';root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
cases=[(3,192,3,0,1,'tensor'),(12,1352,2,0,1,'tensor'),(12,1352,3,1,0,'tensor'),(12,1352,3,1,0,'tensor-real')];records=[]
for l,a,nobj,batch,contacts,backend in cases:
 name=f'l{l}-a{a}-o{nobj}-b{batch}-c{contacts}-{backend}';out=root/name;out.mkdir();ref=out/'reference';ref.mkdir()
 base=list(map(str,[repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',R/'sample-grids'/f'l{l}-a{a}.csv']))
 e=env|{'DSFAST_CONTACTS':str(contacts),'DSFAST_BATCH_RESPONSES':str(batch),'DSFAST_NATIVE_DERIVATIVE':'physical','DSFAST_FORWARD_FIELDS':'1','DSFAST_EXPORT_PACKET':'0'}
 cmd=[str(R/'bin/frozen-reference-v2')]+base+list(map(str,[ref,l,3,nobj,backend,'direct']))
 with (out/'reference.stdout').open('w') as so,(out/'reference.stderr').open('w') as se:subprocess.run(cmd,stdout=so,stderr=se,env=e,preexec_fn=guard,timeout=210,check=True)
 for cap in [0,134217728]:
  data=out/f'cap{cap}';data.mkdir();cmd=[str(R/'fast/fast_driver-v6')]+base+list(map(str,[data,l,3,nobj,backend,'direct','equilibrated']))
  with (out/f'cap{cap}.stdout').open('w') as so,(out/f'cap{cap}.stderr').open('w') as se:subprocess.run(cmd,stdout=so,stderr=se,env=e|{'DSFAST_FORWARD_CACHE_BYTES':str(cap)},preexec_fn=guard,timeout=210,check=True)
  assert (data/'coordinates.bin').read_bytes()==(ref/'coordinates.bin').read_bytes();c=compare(data/'arrays_candidate.bin',ref/'arrays_reference.bin');meta=json.loads((data/'result.json').read_text());rec={'case':name,'cap_bytes':cap,'command':cmd,'comparison':c,'metrics':meta,'data':str(data)};records.append(rec);(root/'records.json').write_text(json.dumps(records,indent=2)+'\n');print(name,cap,c['passed'],'hits',meta['forward_tensor_cache_hits'],'skips',meta['forward_tensor_cache_admission_skips'],flush=True)
assert all(x['comparison']['passed'] for x in records)
