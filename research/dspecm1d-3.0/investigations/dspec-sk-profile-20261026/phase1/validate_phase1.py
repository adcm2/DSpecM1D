"""Parent validation of approved Phase 1 against frozen accepted binaries."""
from pathlib import Path
import csv,gzip,hashlib,importlib.util,json,os,shutil,statistics,subprocess,tempfile,time
BASE=Path(__file__).resolve().parent
ROOT=BASE/os.environ.get('DSPEC_PHASE1_ATTEMPT','attempt-01')
ROOT.mkdir(exist_ok=True)
(ROOT/'bin').mkdir(exist_ok=True)
for name in ['before-fields.csv.gz','before-contacts.csv.gz']:
 if not (ROOT/name).exists():shutil.copyfile(BASE/name,ROOT/name)
for name in ['before-timing','before-fields']:
 if not (ROOT/'bin'/name).exists():shutil.copy2(BASE/'bin'/name,ROOT/'bin'/name)
REPO=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
ENV=os.environ.copy()
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): ENV[k]='1'
ENV['OMP_DYNAMIC']='FALSE';ENV['MPLCONFIGDIR']=str(ROOT/'mpl-cache')
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def run(command,out):
 out=Path(out)
 if out.exists():raise RuntimeError(f'Refusing nonempty/reused output: {out}')
 out.mkdir(parents=True)
 start=time.monotonic();peak_rss=0;peak_hwm=0
 with out.with_suffix('.stdout.log').open('w') as stdout,out.with_suffix('.stderr.log').open('w') as stderr:
  p=subprocess.Popen(list(map(str,command)),stdout=stdout,stderr=stderr,env=ENV)
  while p.poll() is None:
   try:
    lines=Path(f'/proc/{p.pid}/status').read_text().splitlines()
    for line in lines:
     if line.startswith('VmRSS:'):peak_rss=max(peak_rss,int(line.split()[1])*1024)
     if line.startswith('VmHWM:'):peak_hwm=max(peak_hwm,int(line.split()[1])*1024)
   except (OSError,ValueError):pass
   if time.monotonic()-start>180:p.kill();raise RuntimeError('Bounded 180 s timeout')
   time.sleep(.01)
 if p.returncode:raise RuntimeError((command,out.with_suffix('.stderr.log').read_text()[-4000:]))
 return {'command':list(map(str,command)),'wall_s':time.monotonic()-start,'peak_sampled_rss_bytes':peak_rss,'peak_sampled_hwm_bytes':peak_hwm}
def compress(p):
 with p.open('rb') as f,gzip.open(str(p)+'.gz','wb',compresslevel=1) as out:shutil.copyfileobj(f,out)
 p.unlink()
def build():
 commands=json.loads((BASE.parent/'build_commands.json').read_text())
 command=commands['current'];command[command.index(str(BASE.parent/'timing_current.cpp'))]=str(BASE.parent/'timing_current.cpp');command[command.index('-o')+1]=str(ROOT/'bin/after-timing')
 fs=json.loads((REPO/'build/kernel-reconstruction-pole-frames-20261026/optimized_field_command.json').read_text())
 fs[fs.index('/tmp/dspec-pole-frames-20261026/fields.cpp')]=str(REPO/'build/kernel-reconstruction-pole-frames-20261026/fields.cpp');fs[fs.index('-o')+1]=str(ROOT/'bin/after-fields')
 (ROOT/'build_commands.json').write_text(json.dumps({'timing':command,'fields':fs},indent=2)+'\n')
 for name,cmd in [('timing',command),('fields',fs)]:
  with (ROOT/f'build-{name}.log').open('w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=f)
  print('Built',name,flush=True)
def timings():
 records=[]
 for degree,folder in [(3,'small-05'),(12,'comparison-illustration-12-02')]:
  samples=REPO/'build/kernel-examples-phase12'/folder/'samples.csv'
  for measurements in (1,2):
   for repeat,labels in enumerate([('before','after','yspec'),('yspec','after','before'),('after','before','yspec')]):
    for label in labels:
     exe=(REPO/'build/phase-d-build-cache/build-yspecpp/yspec_kernel_timing') if label=='yspec' else ROOT/'bin'/f'{label}-timing'
     out=ROOT/'timings'/f'l{degree}-m{measurements}-r{repeat}-{label}'
     cmd=[exe,REPO/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',samples,out,degree,measurements,1]
     r=run(cmd,out);meta=json.loads(out.with_suffix('.stdout.log').read_text().strip().splitlines()[-1])
     r.update({'solver':label,'degree':degree,'measurements':measurements,'repeat':repeat,'metrics':meta,'api_s':meta['first_workload_api_s'] if label=='yspec' else sum(meta['api_s']),'csv_sha256':digest(out/'kernels.csv'),'output':str(out)})
     compress(out/'kernels.csv');records.append(r)
     (ROOT/'timing_records.json').write_text(json.dumps(records,indent=2)+'\n');print(out.name,round(r['api_s'],6),flush=True)
 summary=[]
 for degree in (3,12):
  for measurements in (1,2):
   row={'degree':degree,'measurements':measurements}
   for label in ('before','after','yspec'):
    cases=[r for r in records if r['degree']==degree and r['measurements']==measurements and r['solver']==label];t=[r['api_s'] for r in cases]
    row[label]={'median_s':statistics.median(t),'min_s':min(t),'max_s':max(t),'peak_sampled_rss_bytes':max(r['peak_sampled_rss_bytes'] for r in cases),'peak_sampled_hwm_bytes':max(r['peak_sampled_hwm_bytes'] for r in cases)}
   row['speedup']=row['before']['median_s']/row['after']['median_s'];row['after_over_yspec']=row['after']['median_s']/row['yspec']['median_s'];summary.append(row)
 (ROOT/'timing_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 target=next(r for r in summary if r['degree']==12 and r['measurements']==1)
 print('Performance gate',target['speedup']>=3,'speedup',target['speedup'],flush=True)
def compare_numeric(x,y,point_limit=1e-10):
 import numpy as np
 assert np.isfinite(x).all() and np.isfinite(y).all()
 l2=float(np.linalg.norm(x-y)/max(np.linalg.norm(x),np.linalg.norm(y),1e-300));peak=max(float(np.max(abs(x))),float(np.max(abs(y))),1e-300)
 point=float(np.max(abs(x-y)/np.maximum(np.maximum(abs(x),abs(y)),.01*peak)))
 return {'relative_L2':l2,'L2_limit':1e-12,'max_pointwise_scaled':point,'pointwise_limit':point_limit,'passed':l2<=1e-12 and point<=point_limit}
def compare_kernels():
 import numpy as np
 spec=importlib.util.spec_from_file_location('runtime',REPO/'benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py');runtime=importlib.util.module_from_spec(spec);spec.loader.exec_module(runtime)
 records=json.loads((ROOT/'timing_records.json').read_text());comparisons=[]
 for degree,folder in [(3,'small-05'),(12,'comparison-illustration-12-02')]:
  for m in (1,2):
   cases={s:next(r for r in records if r['degree']==degree and r['measurements']==m and r['solver']==s and r['repeat']==0) for s in ('before','after','yspec')}
   def rows(label):
    with gzip.open(Path(cases[label]['output'])/'kernels.csv.gz','rt') as f:return list(csv.DictReader(f))
   a,b=rows('before'),rows('after');assert len(a)==len(b)
   checks=[]
   for obj in ('J_Z','J_N')[:m]:
    for mod in 'ACFLN':
     aa=[r for r in a if r['objective']==obj and r['modulus']==mod];bb=[r for r in b if r['objective']==obj and r['modulus']==mod]
     assert len(aa)==len(bb)>0;value='kernel_per_Pa_m3'
     assert all(all(x[k]==y[k] for k in x if k!=value) for x,y in zip(aa,bb))
     checks.append({'objective':obj,'modulus':mod,**compare_numeric(np.array([float(x[value]) for x in aa]),np.array([float(x[value]) for x in bb]))})
   with tempfile.TemporaryDirectory(prefix='dspec-phase1-kernels-') as t:
    paths={}
    for s in ('after','yspec'):
     paths[s]=Path(t)/f'{s}.csv'
     with gzip.open(Path(cases[s]['output'])/'kernels.csv.gz','rb') as f,paths[s].open('wb') as out:shutil.copyfileobj(f,out)
    gate=runtime.accuracy_gate(paths['after'],paths['yspec'],REPO/'build/kernel-examples-phase12'/folder/'samples.csv',1e-3,m,{'objective_scale_S':cases['after']['metrics']['objective_scale_S']})
   result={'degree':degree,'measurements':m,'old_new':checks,'yspec':gate};comparisons.append(result)
   (ROOT/'kernel_comparison.json').write_text(json.dumps(comparisons,indent=2)+'\n');print('Kernels',degree,m,all(r['passed'] for r in checks),flush=True)
 assert all(c['passed'] for r in comparisons for c in r['old_new'])
def fields():
 import numpy as np
 out=ROOT/'field-run';cmd=[ROOT/'bin/after-fields',out/'contacts.csv'];record=run(cmd,out)
 shutil.copyfile(out.with_suffix('.stdout.log'),out/'fields.csv');compress(out/'fields.csv');compress(out/'contacts.csv')
 arrays=[]
 for name in ['before-fields.csv.gz','field-run/fields.csv.gz']:
  with gzip.open(ROOT/name,'rt') as f:arrays.append(np.loadtxt(f,delimiter=','))
 a,b=arrays;assert a.shape==b.shape and np.array_equal(a[:,:5],b[:,:5]);assert np.isfinite(a).all() and np.isfinite(b).all()
 x,y=a[:,5::2]+1j*a[:,6::2],b[:,5::2]+1j*b[:,6::2]
 groups={'forward_raw_displacement':(0,3),'observable_SI_displacement':(3,6),'adjoint_dual_displacement':(6,9),'forward_raw_strain':(9,15),'adjoint_dual_strain':(15,21),'forward_rscaled_strain':(21,27),'adjoint_rscaled_strain':(27,33)};checks=[]
 for case,freq in np.unique(a[:,:2],axis=0):
  mask=(a[:,0]==case)&(a[:,1]==freq)
  for key,(lo,hi) in groups.items():
   aa,bb=x[mask,lo:hi],y[mask,lo:hi];norm=float(np.linalg.norm(aa-bb)/max(np.linalg.norm(aa),np.linalg.norm(bb),1e-30));limit=1e-9 if key.startswith('adjoint') and 'strain' in key else 1e-12
   checks.append({'case':int(case),'frequency':int(freq),'field':key,'relative_norm':norm,'limit':limit,'passed':norm<=limit})
 contacts=[];ca=[]
 for name in ['before-contacts.csv.gz','field-run/contacts.csv.gz']:
  with gzip.open(ROOT/name,'rt') as f:ca.append(np.loadtxt(f,delimiter=','))
 aa,bb=ca;assert aa.shape==bb.shape and np.array_equal(aa[:,:3],bb[:,:3])
 for case,freq in np.unique(aa[:,:2],axis=0):
  mask=(aa[:,0]==case)&(aa[:,1]==freq)
  for col,mod in enumerate('ACFLN',start=3):contacts.append({'case':int(case),'frequency':int(freq),'modulus':mod,**compare_numeric(aa[mask,col],bb[mask,col],1e-12)})
 result={'run':record,'field_checks':checks,'contacts':contacts,'contacts_byte_numeric_identical':bool(np.array_equal(aa,bb)),'passed':all(r['passed'] for r in checks+contacts)}
 (ROOT/'field_contact_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print('Fields/contact',result['passed'],flush=True)
 assert result['passed']

def plane():
 """Old/new equivalence outside the regular shell; no YSpec validity claim."""
 import numpy as np
 samples=REPO/'build/kernels-plot-redesign/plane-run-07/points.csv';model=REPO/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck'
 results={}
 for label in ('before','after'):
  out=ROOT/f'plane-{label}';r=run([ROOT/'bin'/f'{label}-timing',model,samples,out,12,1,1],out)
  with (out/'kernels.csv').open() as f:results[label]=list(csv.DictReader(f))
  compress(out/'kernels.csv')
 checks=[]
 for mod in 'ACFLN':
  a=[r for r in results['before'] if r['modulus']==mod];b=[r for r in results['after'] if r['modulus']==mod];value='kernel_per_Pa_m3'
  assert len(a)==len(b)>0 and all(all(x[k]==y[k] for k in x if k!=value) for x,y in zip(a,b))
  checks.append({'modulus':mod,**compare_numeric(np.array([float(x[value]) for x in a]),np.array([float(x[value]) for x in b]))})
 (ROOT/'plane_comparison.json').write_text(json.dumps({'scope':'Old/new equality only, not receiver-region validation against YSpec.','checks':checks,'passed':all(r['passed'] for r in checks)},indent=2)+'\n')
 print('Full plane',all(r['passed'] for r in checks),flush=True);assert all(r['passed'] for r in checks)

if __name__=='__main__':
 import sys
 for command in sys.argv[1:]:globals()[command]()
