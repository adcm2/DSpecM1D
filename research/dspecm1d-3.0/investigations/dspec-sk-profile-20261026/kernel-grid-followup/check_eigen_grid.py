from pathlib import Path
import os,json,subprocess,resource,hashlib,importlib.util
F=Path(__file__).resolve().parent;A=F.parent/'algorithm-audit';M=F.parent/'math-profile'
spec=importlib.util.spec_from_file_location('arrays',A/'validate_arrays.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:env[k]='1'
env.update(DSFAST_CONTACTS='0',DSFAST_BATCH_RESPONSES='0',DSFAST_FORWARD_FIELDS='0',DSFAST_NATIVE_DERIVATIVE='physical',DSFAST_EXPORT_PACKET='0')
def guard():
 os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
root=F/'evidence/regression-v2';root.mkdir()
records=[]
for l in [3,12]:
 dirs={}
 for mode,exe in [('eigen',F/'eigen-candidate/eigen_candidate_driver'),('frozen_cblas',A/'fast/fast_driver-v6')]:
  out=root/f'l{l}-{mode}';out.mkdir();dirs[mode]=out
  cmd=list(map(str,[exe,'/home/adcm2/Documents/Research/Projects/DSpecM1D/tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',F/'evidence'/f'gsh-kernel-l{l}.csv',out,l,3,1,'tensor','direct','equilibrated']))
  p=subprocess.run(cmd,capture_output=True,text=True,env=env,preexec_fn=guard,timeout=200);(out.parent/f'{out.name}-stdout.log').write_text(p.stdout);(out.parent/f'{out.name}-stderr.log').write_text(p.stderr);assert p.returncode==0,p.stderr[-3000:]
  result=json.loads((out/'result.json').read_text());assert result['max_kernel_l2_relative_error']<=1e-12 and result['max_kernel_pointwise_scaled_error']<=1e-10 and result['max_field_relative_error']<=1e-12
  records.append({'degree':l,'mode':mode,'command':cmd,'result':result,'output_sha256':hashlib.sha256((out/'arrays_candidate.bin').read_bytes()).hexdigest()})
 gate=v.compare(dirs['eigen']/'arrays_candidate.bin',dirs['frozen_cblas']/'arrays_candidate.bin');assert gate['passed'];records.append({'degree':l,'comparison':'Eigen internal vs frozen explicit CBLAS at identical native kernel-band grid','gate':gate})
(F/'evidence/eigen-grid-gates.json').write_text(json.dumps(records,indent=2)+'\n')
print('Two new-grid Eigen/production and Eigen/frozen-CBLAS regression cases pass; not a physical benchmark.')
