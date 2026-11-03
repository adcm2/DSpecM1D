from pathlib import Path
import subprocess,json,os,time,gzip,hashlib,resource
ROOT=Path(__file__).resolve().parent
REPO=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
MODEL=REPO/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck'
ENV=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:ENV[k]='1'
CPU=min(os.sched_getaffinity(0))
def guard():
 os.sched_setaffinity(0,{CPU});resource.setrlimit(resource.RLIMIT_CPU,(120,120));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
records=[]
for degree in [3,12]:
 for measurements in [1,2]:
  for kind in ['native','profile']:
   out=ROOT/'current-profile-02'/f'l{degree}-m{measurements}-{kind}'
   out.mkdir(parents=True,exist_ok=False)
   samples=REPO/('build/kernel-examples-phase12/comparison-small-03/samples.csv' if degree==3 else 'build/kernel-examples-phase12/comparison-illustration-12-02/samples.csv')
   exe=ROOT/'bin/profile-current' if kind=='profile' else ROOT.parent/'phase1/attempt-02/bin/after-timing'
   data=out/'data';data.mkdir()
   command=list(map(str,[exe,MODEL,samples,data,degree,measurements,1]))
   start=time.monotonic()
   with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
    proc=subprocess.run(command,stdout=stdout,stderr=stderr,env=ENV,preexec_fn=guard,timeout=150)
   if proc.returncode:raise RuntimeError(str(out)+': '+(out/'stderr.log').read_text()[-2000:])
   metrics=json.loads((out/'stdout.log').read_text().splitlines()[-1])
   profile=[json.loads(x[len('AUDITPROFILE '):]) for x in (out/'stderr.log').read_text().splitlines() if x.startswith('AUDITPROFILE ')]
   csv=data/'kernels.csv';sha=hashlib.sha256(csv.read_bytes()).hexdigest()
   with csv.open('rb') as src,gzip.open(str(csv)+'.gz','wb',compresslevel=3) as dst:
    import shutil;shutil.copyfileobj(src,dst)
   csv.unlink()
   r={'degree':degree,'measurements':measurements,'kind':kind,'cpu':CPU,'command':command,'metrics':metrics,'profile':profile,'csv_sha256':sha,'wall_s':time.monotonic()-start}
   (out/'record.json').write_text(json.dumps(r,indent=2)+'\n');records.append(r)
   (ROOT/'current_profile_records.json').write_text(json.dumps(records,indent=2)+'\n')
   print(degree,measurements,kind,metrics['calculation_s']-metrics['setup_s'],profile,flush=True)
for degree in [3,12]:
 for m in [1,2]:
  rr=[r for r in records if r['degree']==degree and r['measurements']==m]
  assert rr[0]['csv_sha256']==rr[1]['csv_sha256'],(degree,m)
print('All instrumented outputs identical to native.',flush=True)
