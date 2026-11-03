from pathlib import Path
import os,subprocess,time,json,resource,hashlib,statistics,sys,gzip,shutil
from validate_arrays import compare
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
root=R/(sys.argv[1] if len(sys.argv)>1 else 'final-timings-v5');root.mkdir(exist_ok=False)
fast=R/'fast/fast_driver-v5';old=R/'bin/frozen-reference-v2';ys=R/'bin/yspec-audit';exes=[fast,old,ys]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in exes}
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
cpu=min(os.sched_getaffinity(0))
def guard():
 os.sched_setaffinity(0,{cpu});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
cases=[(3,3,1),(6,3,1),(12,3,1),(12,3,2),(12,12,1),(12,48,1)]
records=[]
for rep in range(3):
 for l,nf,nobj in (cases if rep%2==0 else list(reversed(cases))):
  a=192 if l==3 else 1352
  for kind in (['old','fast','yspec'] if rep%2==0 else ['yspec','fast','old']):
   name=f'l{l}-f{nf}-o{nobj}-r{rep}-{kind}';out=root/name;out.mkdir();data=out/'data';data.mkdir()
   model=repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck';sample=R/'sample-grids'/f'l{l}-a{a}.csv'
   e=env|{'DSFAST_CONTACTS':'0','DSFAST_BATCH_RESPONSES':str(int(nobj>1)),'DSFAST_NATIVE_DERIVATIVE':'physical','DSFAST_EXPORT_PACKET':'0','SK_BINS':str(nf)}
   if kind=='yspec':cmd=list(map(str,[ys,model,sample,data,l,nobj,1,a]));exe=ys
   elif kind=='old':cmd=list(map(str,[old,model,sample,data,l,nf,nobj,'tensor-real','direct']));exe=old
   else:cmd=list(map(str,[fast,model,sample,data,l,nf,nobj,'tensor-real','direct','equilibrated']));exe=fast
   peak=0;samples=0;start=time.monotonic()
   with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
    p=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,env=e,preexec_fn=guard)
    while p.poll() is None:
     if time.monotonic()-start>210:p.kill();raise RuntimeError('Timeout '+name)
     try:
      if Path(os.readlink(f'/proc/{p.pid}/exe'))==exe.resolve():
       status=Path(f'/proc/{p.pid}/status').read_text();v=next(x for x in status.splitlines() if x.startswith('VmHWM:'));peak=max(peak,int(v.split()[1])*1024);samples+=1
     except (FileNotFoundError,ProcessLookupError,StopIteration,PermissionError):pass
     time.sleep(.01)
   if p.returncode:raise RuntimeError(name+' '+(out/'stderr.log').read_text()[-2000:])
   meta=json.loads((out/'stdout.log').read_text().splitlines()[-1]) if kind=='yspec' else json.loads((data/'result.json').read_text())
   seconds=meta['api_s'][0] if kind=='yspec' else sum(meta['api_s'] if kind=='old' else meta['candidate_s'])
   if kind=='yspec':
    csv=data/'kernels.csv';h=hashlib.sha256(csv.read_bytes()).hexdigest()
    with csv.open('rb') as f,gzip.open(str(csv)+'.gz','wb',compresslevel=3) as g:shutil.copyfileobj(f,g)
    csv.unlink()
   else:h=hashlib.sha256((data/('arrays_reference.bin' if kind=='old' else 'arrays_candidate.bin')).read_bytes()).hexdigest()
   rec={'name':name,'degree':l,'bins':nf,'objectives':nobj,'repeat':rep,'kind':kind,'seconds':seconds,'peak_process_hwm_bytes':peak,'hwm_samples':samples,'cpu':cpu,'command':cmd,'environment_controls':{k:e[k] for k in ['DSFAST_CONTACTS','DSFAST_BATCH_RESPONSES','DSFAST_EXPORT_PACKET','DSFAST_NATIVE_DERIVATIVE','SK_BINS']},'metrics':meta,'output_sha256':h,'data':str(data),'wall_s':time.monotonic()-start}
   records.append(rec);(root/'records.json').write_text(json.dumps(records,indent=2)+'\n');print(name,seconds,peak,flush=True)
summary=[]
for l,nf,nobj in cases:
 for kind in ['old','fast','yspec']:
  rows=[x for x in records if (x['degree'],x['bins'],x['objectives'],x['kind'])==(l,nf,nobj,kind)];assert len(rows)==3
  summary.append({'degree':l,'bins':nf,'objectives':nobj,'kind':kind,'processes':len(rows),'median_s':statistics.median(x['seconds'] for x in rows),'min_s':min(x['seconds'] for x in rows),'max_s':max(x['seconds'] for x in rows),'peak_process_hwm_bytes':max(x['peak_process_hwm_bytes'] for x in rows),'deterministic':len({x['output_sha256'] for x in rows})==1})
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
gates=[]
for row in records:
 if row['kind']!='fast':continue
 ref=next(x for x in records if (x['degree'],x['bins'],x['objectives'],x['repeat'],x['kind'])==(row['degree'],row['bins'],row['objectives'],row['repeat'],'old'))
 cand=Path(row['data']);olddata=Path(ref['data']);assert (cand/'coordinates.bin').read_bytes()==(olddata/'coordinates.bin').read_bytes()
 gates.append({'name':row['name'],**compare(cand/'arrays_candidate.bin',olddata/'arrays_reference.bin')})
(root/'frozen_comparisons.json').write_text(json.dumps(gates,indent=2)+'\n')
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha for p,sha in hashes.items())
(root/'provenance.json').write_text(json.dumps({'executables':hashes,'runtime_cpu':cpu,'runtime_threads':1,'setup_and_output_excluded_from_API_time':True,'contacts_included':False,'memory':'10ms sampled current executable /proc VmHWM; includes all process phases and internal same-engine reference in fast driver; not aggregate or hard cap','all_frozen_array_comparisons_pass':all(x['passed'] for x in gates)},indent=2)+'\n')
print('Gates',all(x['passed'] for x in gates),flush=True)
