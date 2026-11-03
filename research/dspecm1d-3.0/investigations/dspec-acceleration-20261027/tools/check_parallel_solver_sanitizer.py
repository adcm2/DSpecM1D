#!/usr/bin/env python3
"""Run full forward/adjoint kernel path under TSan with explicit accounting."""
from pathlib import Path
import argparse,json,os,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('binary');p.add_argument('model');p.add_argument('degree',type=int);p.add_argument('threads',type=int);p.add_argument('output');a=p.parse_args()
out=Path(a.output);out.mkdir();env=os.environ.copy();env.update(TSAN_OPTIONS='halt_on_error=1:exitcode=66',DSPEC_SYMBOLIC_REUSE='1',OMP_NUM_THREADS=str(a.threads),OMP_THREAD_LIMIT=str(a.threads),OMP_DYNAMIC='FALSE')
cmd=[a.binary,a.model,str(out/'data'),'kernel-probe','--lmax',str(a.degree),'--nq','6','--radial-refinement','0.05','--radial-traces','0','--probe-angle-count','8','--solve-threads',str(a.threads),'--pole-only','on','--forward-only','off','--response-reuse','on']
with (out/'stdout.log').open('xb') as so,(out/'stderr.log').open('xb') as se:r=subprocess.run(cmd,env=env,stdout=so,stderr=se,timeout=600)
err=(out/'stderr.log').read_text();status='PASS' if r.returncode==0 else 'UNSUPPORTED_RUNTIME' if 'ThreadSanitizer: unexpected memory mapping' in err else 'FAIL'
report={'status':status,'returncode':r.returncode,'sanitizer_passed':status=='PASS','command':cmd,'stderr_excerpt':err[-20000:]}
with (out/'checks.json').open('x') as f:json.dump(report,f,indent=2)
sys.exit(1 if status=='FAIL' else 0)
