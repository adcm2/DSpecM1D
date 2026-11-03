from pathlib import Path
import json,os,shutil,subprocess,sys,importlib.util,resource,hashlib
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D');here=repo/'benchmarks/yspecpp/kernels/3d';out=R/'born-v5';out.mkdir(exist_ok=False)
for name in ['dspec','yspecpp','oracle']:(out/name).mkdir()
shutil.copyfile(R.parent/'phase1/born-01/fixture.json',out/'fixture.json');shutil.copyfile(R.parent/'phase1/born-01/common_homogeneous.deck',out/'common_homogeneous.deck')
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
env|={'DSFAST_CONTACTS':'1','DSFAST_NATIVE_DERIVATIVE':'physical','DSFAST_BATCH_RESPONSES':'0','DSFAST_EXPORT_PACKET':'0'}
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
cmds=[[R/'bin/born-candidate',out/'common_homogeneous.deck',out/'dspec',.05,8],[repo/'build/phase-d-build-cache/build-yspecpp/phase_d_adapter',out/'common_homogeneous.deck',out/'yspecpp',out/'dspec/radial_nodes.csv',out/'dspec/angles.csv'],[sys.executable,here/'cartesian_ti_oracle.py','--baseline',out/'dspec','--output',out/'oracle','--fixture',out/'fixture.json']]
for i,cmd in enumerate(cmds):
 with (out/f'run-{i}.log').open('w') as f:subprocess.run(list(map(str,cmd)),stdout=f,stderr=f,env=env,preexec_fn=guard,timeout=210,check=True)
 print('Born stage',i,'done',flush=True)
spec=importlib.util.spec_from_file_location('phase_d',here/'run_phase_d.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);metrics=m.baseline_metrics(out)
(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');(out/'commands.json').write_text(json.dumps([list(map(str,x)) for x in cmds],indent=2)+'\n')
(out/'provenance.json').write_text(json.dumps({'prototype':'v5 tensor-real, shared existing equilibrated factor, batchoff, contactson, physical-node compatibility derivative','candidate_executable_sha256':hashlib.sha256((R/'bin/born-candidate').read_bytes()).hexdigest(),'runtime_threads':1,'fresh_cartesian_oracle':True,'no_finite_difference_sweep':True,'all_baseline_gates_pass':metrics['all_baseline_gates_pass']},indent=2)+'\n')
print('Fresh optimized Born gates',metrics['all_baseline_gates_pass'],flush=True)
assert metrics['all_baseline_gates_pass']
