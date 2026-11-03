from pathlib import Path
import importlib.util,json,gzip,tempfile,shutil
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
spec=importlib.util.spec_from_file_location('runtime',repo/'benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py');runtime=importlib.util.module_from_spec(spec);spec.loader.exec_module(runtime)
yspec=json.loads((R/'yspec-audit-01/records.json').read_text());native=json.loads((R/'current_profile_records.json').read_text());results=[]
for l,m in [(3,1),(12,1),(12,2)]:
 y=next(x for x in yspec if x['degree']==l and x['measurements']==m and x['frequency_bins']==3 and x['repeat']==0)
 d=next(x for x in native if x['degree']==l and x['measurements']==m and x['kind']=='native')
 csvd=Path(d['command'][3])/'kernels.csv.gz';csvy=Path(y['output'])/'kernels.csv.gz'
 samples=R/'sample-grids'/f'l{l}-a{192 if l==3 else 1352}.csv'
 # Gate requires deltaX/dr/dOmega columns from the original full tutorial table.
 samples=repo/'build/kernel-examples-phase12'/('small-05' if l==3 else 'comparison-illustration-12-02')/'samples.csv'
 with tempfile.TemporaryDirectory(prefix='dspec-audit-gate-') as tmp:
  a=Path(tmp)/'dspec.csv';b=Path(tmp)/'yspec.csv'
  for source,dest in [(csvd,a),(csvy,b)]:
   with gzip.open(source,'rb') as src,dest.open('wb') as out:shutil.copyfileobj(src,out)
  result=runtime.accuracy_gate(a,b,samples,1e-3,m,{'objective_scale_S':d['metrics']['objective_scale_S']})
 results.append({'degree':l,'objectives':m,'yspec':result})
 print(l,m,'gates passed',flush=True)
(R/'yspec_frozen_gates.json').write_text(json.dumps(results,indent=2)+'\n')
