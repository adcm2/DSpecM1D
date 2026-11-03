from pathlib import Path
import os,sys,json,time,hashlib,subprocess,resource,gzip,shutil,importlib.util
R=Path(__file__).resolve().parent;A=R.parent/'algorithm-audit';repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
spec=importlib.util.spec_from_file_location('arrays',A/'validate_arrays.py');arr=importlib.util.module_from_spec(spec);spec.loader.exec_module(arr)
root=R/(sys.argv[1] if len(sys.argv)>1 else 'runs-01');root.mkdir(exist_ok=False)
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
cpu=min(os.sched_getaffinity(0))
def guard():
 os.sched_setaffinity(0,{cpu});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
exes={'dspec':{'native':R/'dspec/fast_driver-profile-disabled-v4','profile':R/'dspec/fast_driver-profiled-v4','original':A/'fast/fast_driver-v6'},'yspec':{'native':R/'yspec/build/v2/native/yspec-math-profile','profile':R/'yspec/build/v2/instrumented/yspec-math-profile','original':A/'bin/yspec-audit'}}
exe_sha={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for dd in exes.values() for p in dd.values()}
cases=[(3,3,1,0),(6,3,1,0),(12,3,1,0),(12,3,2,0),(12,12,1,0),(12,48,1,0)];records=[]
if '--pilot' in sys.argv:cases=[(12,3,1,0)]
if '--contacts' in sys.argv:cases=[(12,3,1,1)]
for rep in range(1 if '--pilot' in sys.argv else 3):
 for l,nf,m,contacts in (cases if rep%2==0 else list(reversed(cases))):
  a=192 if l==3 else 1352
  for solver in (['dspec','yspec'] if rep%2==0 else ['yspec','dspec']):
   if '--contacts' in sys.argv and solver!='dspec':continue
   for mode in (['original','native','profile'] if rep==0 else ['profile','native'] if rep%2 else ['native','profile']):
    out=root/f'l{l}-f{nf}-m{m}-c{contacts}-r{rep}-{solver}-{mode}';out.mkdir();data=out/'data';data.mkdir()
    exe=exes[solver][mode];sample=A/'sample-grids'/f'l{l}-a{a}.csv';model=repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck'
    if solver=='dspec':cmd=list(map(str,[exe,model,sample,data,l,nf,m,'tensor','direct','equilibrated']))
    else:cmd=list(map(str,[exe,model,sample,data,l,m,1,a]))
    e=env|{'SK_BINS':str(nf),'DSFAST_CONTACTS':str(contacts),'DSFAST_BATCH_RESPONSES':str(int(m>1)),'DSFAST_FORWARD_FIELDS':str(int(m>1)),'DSFAST_NATIVE_DERIVATIVE':'physical','DSFAST_EXPORT_PACKET':'0'}
    start=time.monotonic();peak=0
    with (out/'stdout.log').open('w') as so,(out/'stderr.log').open('w') as se:
     p=subprocess.Popen(cmd,stdout=so,stderr=se,env=e,preexec_fn=guard)
     while p.poll() is None:
      if time.monotonic()-start>210:p.kill();raise RuntimeError('timeout '+str(out))
      try:
       if Path(os.readlink(f'/proc/{p.pid}/exe'))==exe.resolve():
        v=next(x for x in Path(f'/proc/{p.pid}/status').read_text().splitlines() if x.startswith('VmHWM:'));peak=max(peak,int(v.split()[1])*1024)
      except (FileNotFoundError,ProcessLookupError,PermissionError,StopIteration):pass
      time.sleep(.01)
    if p.returncode:raise RuntimeError(str(out)+(out/'stderr.log').read_text()[-3000:])
    stdout=(out/'stdout.log').read_text();meta=json.loads((data/'result.json').read_text()) if solver=='dspec' else json.loads(stdout.splitlines()[-1])
    profiles=[json.loads(x[len('MATHPROFILE '):]) for x in (out/'stderr.log').read_text().splitlines() if x.startswith('MATHPROFILE ')]
    for pro in profiles:
     ss=pro['stages'];total=sum(x['exclusive_s'] for x in ss.values());roots=['candidate_api_including_local_cleanup','api.ComputeSampledKernels','frontend.model_sampling_configuration_and_input_setup','frontend.sampled_kernel_csv_export','setup_all']
     recognized=[name for name in roots if name in ss]
     assert len(recognized)==1,('Unknown accounting root',out,list(ss))
     pro['accounting_root']=recognized[0];pro['exclusive_sum_s']=total;pro['closure_error_s']=total-ss[recognized[0]]['inclusive_s'];assert abs(pro['closure_error_s'])<1e-8,(out,pro)
    if mode=='profile':assert profiles,'No profile output '+str(out)
    else:assert not profiles,'Native contains active timers '+str(out)
    target=data/('arrays_candidate.bin' if solver=='dspec' else 'kernels.csv');h=hashlib.sha256(target.read_bytes()).hexdigest()
    if solver=='yspec':
     with target.open('rb') as f,gzip.open(str(target)+'.gz','wb',compresslevel=3) as g:shutil.copyfileobj(f,g)
     target.unlink()
    native_prev=next((x for x in records if (x['degree'],x['bins'],x['objectives'],x['solver'],x['mode'])==(l,nf,m,solver,'original')),None)
    gate=None
    if native_prev:
     if solver=='dspec':gate=arr.compare(target,Path(native_prev['data'])/'arrays_candidate.bin')
     else:gate={'byte_identical':h==native_prev['output_sha256'],'passed':h==native_prev['output_sha256']}
     assert gate['passed'],('Instrumentation equivalence failure',out,gate)
    rec={'degree':l,'bins':nf,'objectives':m,'contacts':contacts,'repeat':rep,'solver':solver,'mode':mode,'command':cmd,'env':{k:e[k] for k in ['SK_BINS','DSFAST_CONTACTS','DSFAST_BATCH_RESPONSES','DSFAST_FORWARD_FIELDS','DSFAST_NATIVE_DERIVATIVE']},'metrics':meta,'profiles':profiles,'peak_current_executable_hwm_bytes':peak,'wall_s':time.monotonic()-start,'output_sha256':h,'original_equivalence':gate,'data':str(data)};records.append(rec)
    (root/'records.json').write_text(json.dumps(records,indent=2)+'\n');print(l,nf,m,solver,mode,meta.get('candidate_s',meta.get('api_s')),gate['passed'] if gate else 'oracle',flush=True)
for p,sha in exe_sha.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha,'binary changed: '+p
(root/'provenance.json').write_text(json.dumps({'executables':exe_sha,'runtime_threads':1,'cpu':cpu,'all_equivalence_gates_pass':all(x['original_equivalence']['passed'] for x in records if x['original_equivalence'])},indent=2)+'\n')
