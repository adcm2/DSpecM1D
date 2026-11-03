from pathlib import Path
import json,subprocess,os,resource,statistics,hashlib
from validate_arrays import compare
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D');root=R/'final-ablation-v5';root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
configs=[('v4-real',R/'fast/fast_driver-v4','tensor-real',0,0,False,3),('v5-real',R/'fast/fast_driver-v5','tensor-real',0,0,False,3),('v5-complex',R/'fast/fast_driver-v5','tensor',0,0,False,3),('packed-scalar',R/'fast/fast_driver-v5','scalar',0,0,False,3),('packed-eigen',R/'fast/fast_driver-v5','eigen',0,0,False,3),('packed-openblas',R/'fast/fast_driver-v5','openblas',0,0,False,3),('minimal-basis',R/'fast/fast_driver-v5','tensor-real',1,0,False,3),('rstrain-bridge',R/'fast/fast_driver-v5','tensor-real',0,1,False,3),('strict-native-build',R/'bin/fast-v5-native','tensor-real',0,0,False,3),('three-trace-f3',R/'fast/fast_driver-v5','tensor-real',0,0,True,3),('three-trace-f48',R/'fast/fast_driver-v5','tensor-real',0,0,True,48)]
records=[]
for rep in range(3):
 for name,exe,backend,batch,bridge,three,nf in (configs if rep%2==0 else list(reversed(configs))):
  out=root/f'{name}-r{rep}';out.mkdir();data=out/'data';data.mkdir();sample=R/'sample-grids'/('l12-three-traces-full-angles.csv' if three else 'l12-a1352.csv')
  cmd=list(map(str,[exe,repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',sample,data,12,nf,1,backend,'direct','equilibrated']))
  e=env|{'DSFAST_CONTACTS':'0','DSFAST_BATCH_RESPONSES':str(batch),'DSFAST_RSTRAIN_BRIDGE':str(bridge),'DSFAST_NATIVE_DERIVATIVE':'physical','DSFAST_EXPORT_PACKET':'0'}
  with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:p=subprocess.run(cmd,stdout=stdout,stderr=stderr,env=e,preexec_fn=guard,timeout=210)
  if p.returncode:raise RuntimeError(name+(out/'stderr.log').read_text()[-2000:])
  meta=json.loads((data/'result.json').read_text());gate=None
  if not three:
   ref=R/'final-validation-v4/l12-f3-a1352-direct-o1-batch0-c1-physical-tensor-real/reference/arrays_reference.bin'
   gate=compare(data/'arrays_candidate.bin',ref)
   for obj in gate['objectives']:
    obj['checks']={k:v for k,v in obj['checks'].items() if not k.startswith('contact_')};obj['passed']=all(x['passed'] for x in obj['checks'].values())
   gate['passed']=all(x['passed'] for x in gate['objectives'])
  rec={'variant':name,'repeat':rep,'command':cmd,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'metrics':meta,'volume_only_frozen_comparison':gate,'data':str(data)};records.append(rec);(root/'records.json').write_text(json.dumps(records,indent=2)+'\n');print(name,rep,meta['candidate_s'],gate['passed'] if gate else 'same-engine only',flush=True)
summary=[]
for name,*_ in configs:
 rows=[x for x in records if x['variant']==name]
 summary.append({'variant':name,'processes':len(rows),'median_api_s':statistics.median(sum(x['metrics']['candidate_s']) for x in rows),'median_sampling_s':statistics.median(sum(y['packed_sampling'] for y in x['metrics']['candidate_stage_s']) for x in rows),'min_api_s':min(sum(x['metrics']['candidate_s']) for x in rows),'max_api_s':max(sum(x['metrics']['candidate_s']) for x in rows),'passed_volume_frozen':all(x['volume_only_frozen_comparison']['passed'] for x in rows) if rows[0]['volume_only_frozen_comparison'] else None})
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
