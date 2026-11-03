"""Single-thread physical case runner; raw evidence retained on failure."""
from pathlib import Path
import argparse,hashlib,json,os,resource,subprocess,time
def sha256_file(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024**2),b''):h.update(chunk)
 return h.hexdigest()
p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--timeout',type=int,default=3600);p.add_argument('--memory-gb',type=int,default=8);p.add_argument('--cpu-index',type=int,default=0);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
root=Path(__file__).resolve().parent/'runs';root.mkdir(exist_ok=True);out=root/a.name;out.mkdir(exist_ok=False)
cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
cmd=[s.replace('{OUTPUT}',str(out/'data')) for s in cmd]
if not cmd:raise SystemExit('missing command')
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
available=sorted(os.sched_getaffinity(0))
if a.cpu_index<0 or a.cpu_index>=len(available):raise SystemExit('cpu-index outside allowed affinity')
cpu=available[a.cpu_index]
def guard():
 os.sched_setaffinity(0,{cpu});resource.setrlimit(resource.RLIMIT_CPU,(a.timeout,a.timeout));resource.setrlimit(resource.RLIMIT_AS,(a.memory_gb*1024**3,a.memory_gb*1024**3))
record={'command':cmd,'cpu':cpu,'threads':1,'timeout_s':a.timeout,'address_space_limit_gb':a.memory_gb,'executable_sha256':sha256_file(cmd[0]),'start_unix':time.time(),'cpu_index':a.cpu_index}
record['runner_sha256']=sha256_file(__file__)
record['underflow_diagnostic_environment']={k:env[k] for k in ['DSPEC_PROFILE_FTZ','LD_PRELOAD'] if k in env}
if env.get('LD_PRELOAD'):
 record['preload_sha256']={s:sha256_file(s) for s in env['LD_PRELOAD'].replace(':',' ').split() if Path(s).is_file()}
(out/'command.json').write_text(json.dumps(record,indent=2)+'\n');start=time.monotonic();peak=0;timedout=False
with (out/'stdout.log').open('w') as so,(out/'stderr.log').open('w') as se:
 proc=subprocess.Popen(cmd,stdout=so,stderr=se,env=env,preexec_fn=guard)
 while proc.poll() is None:
  if time.monotonic()-start>a.timeout+30:proc.kill();timedout=True
  try:
   status=Path(f'/proc/{proc.pid}/status').read_text().splitlines()
   peak=max(peak,int(next(s for s in status if s.startswith('VmHWM:')).split()[1])*1024)
  except (FileNotFoundError,ProcessLookupError,PermissionError,StopIteration):pass
  time.sleep(.1)
record.update(returncode=proc.returncode,wall_s=time.monotonic()-start,peak_rss_bytes=peak,timedout=timedout)
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
record.update(child_user_cpu_s=usage.ru_utime,child_system_cpu_s=usage.ru_stime,
              peak_rss_bytes=max(peak,usage.ru_maxrss*1024))
record['profiles']=[json.loads(s[12:]) for s in (out/'stderr.log').read_text().splitlines() if s.startswith('MATHPROFILE ')]
record['output_files']={str(x.relative_to(out)):sha256_file(x) for x in (out/'data').rglob('*') if x.is_file()} if (out/'data').is_dir() else {}
(out/'run.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:record[k] for k in ['returncode','wall_s','peak_rss_bytes','timedout','child_user_cpu_s','child_system_cpu_s']}),flush=True)
raise SystemExit(proc.returncode or 0)
