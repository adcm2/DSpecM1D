from pathlib import Path
import os,subprocess,resource,json,time,gzip,shutil,hashlib,statistics
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
baseenv=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:baseenv[k]='1'
CPU=min(os.sched_getaffinity(0))
def guard():
 os.sched_setaffinity(0,{CPU});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
records=[]
root=R/'cache-ablation-01';root.mkdir(exist_ok=False)
for repeat in range(3):
 cases=[(4096,32),(1352,32),(2704,32),(75712,32),(4096,128)]
 if repeat%2:cases.reverse()
 for block,mb in cases:
  out=root/f'b{block}-mb{mb}-r{repeat}';out.mkdir();data=out/'data';data.mkdir()
  samples=repo/'build/kernel-examples-phase12/comparison-illustration-12-02/samples.csv'
  cmd=list(map(str,[R/'bin/block-ablation',repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',samples,data,12,1,1]))
  env=baseenv|{'SK_BLOCK':str(block),'SK_PLAN_MB':str(mb)}
  with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
   p=subprocess.run(cmd,env=env,preexec_fn=guard,stdout=stdout,stderr=stderr,timeout=210)
  if p.returncode:raise RuntimeError((out/'stderr.log').read_text())
  m=json.loads((out/'stdout.log').read_text().splitlines()[-1]);csv=data/'kernels.csv';sha=hashlib.sha256(csv.read_bytes()).hexdigest()
  with csv.open('rb') as src,gzip.open(str(csv)+'.gz','wb',compresslevel=3) as dst:shutil.copyfileobj(src,dst)
  csv.unlink()
  rec={'block':block,'plan_budget_mb':mb,'repeat':repeat,'cpu':CPU,'command':cmd,'env':{'SK_BLOCK':str(block),'SK_PLAN_MB':str(mb)},'metrics':m,'sha256':sha}
  (out/'record.json').write_text(json.dumps(rec,indent=2)+'\n');records.append(rec)
  (root/'records.json').write_text(json.dumps(records,indent=2)+'\n')
  print(block,mb,repeat,m['api_s'][0],flush=True)
assert len({r['sha256'] for r in records})==1
summary=[{'block':b,'plan_budget_mb':mb,'median_s':statistics.median(x['metrics']['api_s'][0] for x in records if x['block']==b and x['plan_budget_mb']==mb),'peak_rss_bytes':max(x['metrics']['peak_rss_bytes'] for x in records if x['block']==b and x['plan_budget_mb']==mb)}for b,mb in [(4096,32),(1352,32),(2704,32),(75712,32),(4096,128)]]
(root/'summary.json').write_text(json.dumps({'cases':summary,'all_arrays_byte_identical':True},indent=2)+'\n')
print(summary,flush=True)
