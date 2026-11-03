from pathlib import Path
import os,subprocess,resource,json,time,gzip,shutil,hashlib,statistics
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
ENV=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:ENV[k]='1'
CPU=min(os.sched_getaffinity(0))
def guard():
 os.sched_setaffinity(0,{CPU});resource.setrlimit(resource.RLIMIT_CPU,(180,180));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
root=R/'yspec-audit-01';root.mkdir(exist_ok=False);records=[]
cases=[(3,3,1,0),(12,3,1,0),(12,3,2,0),(3,3,1,1),(12,3,2,1),(12,3,1,1),(12,3,1,2),(12,3,2,2),(3,3,1,2),(6,3,1,0),(12,12,1,0),(12,48,1,0),(3,12,1,0),(3,48,1,0)]
for l,nf,m,rep in cases:
 out=root/f'l{l}-f{nf}-m{m}-r{rep}';out.mkdir();data=out/'data';data.mkdir()
 a=192 if l==3 else 1352
 samples=R/'sample-grids'/f'l{l}-a{a}.csv'
 cmd=list(map(str,[R/'bin/yspec-audit',repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck',samples,data,l,m,1,a]))
 env=ENV|{'SK_BINS':str(nf)}
 with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
  p=subprocess.run(cmd,env=env,preexec_fn=guard,stdout=stdout,stderr=stderr,timeout=210)
 if p.returncode:raise RuntimeError((out/'stderr.log').read_text())
 meta=json.loads((out/'stdout.log').read_text().splitlines()[-1]);csv=data/'kernels.csv';sha=hashlib.sha256(csv.read_bytes()).hexdigest()
 with csv.open('rb') as src,gzip.open(str(csv)+'.gz','wb',compresslevel=3) as dst:shutil.copyfileobj(src,dst)
 csv.unlink()
 rec={'degree':l,'frequency_bins':nf,'measurements':m,'repeat':rep,'cpu':CPU,'command':cmd,'env':{'SK_BINS':str(nf)},'metrics':meta,'sha256':sha,'output':str(data)}
 (out/'record.json').write_text(json.dumps(rec,indent=2)+'\n');records.append(rec)
 (root/'records.json').write_text(json.dumps(records,indent=2)+'\n')
 print(l,nf,m,rep,meta['api_s'][0],flush=True)
summary=[]
for l,nf,m in dict.fromkeys((x['degree'],x['frequency_bins'],x['measurements']) for x in records):
 group=[x for x in records if (x['degree'],x['frequency_bins'],x['measurements'])==(l,nf,m)]
 assert len({x['sha256'] for x in group})==1
 summary.append({'degree':l,'frequency_bins':nf,'measurements':m,'processes':len(group),'median_s':statistics.median(x['metrics']['api_s'][0] for x in group),'min_s':min(x['metrics']['api_s'][0] for x in group),'max_s':max(x['metrics']['api_s'][0] for x in group),'peak_rss_bytes':max(x['metrics']['peak_rss_bytes'] for x in group),'max_solved_frequency_hz':nf/2560,'deterministic':True})
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
