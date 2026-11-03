#!/usr/bin/env python3
"""Diagnostic controls for TSan/OpenMP interoperability, not a passing race gate."""
from pathlib import Path
import hashlib,json,os,subprocess,sys
out=Path(sys.argv[1]);out.mkdir();binary=sys.argv[2];model=sys.argv[3]
source='''#include <vector>\n#include <numeric>\n#include <iostream>\nint main(){std::vector<int> x(4096);\n#pragma omp parallel for num_threads(3)\nfor(int i=0;i<4096;++i)x[i]=i;\nstd::cout<<std::accumulate(x.begin(),x.end(),0LL)<<"\\n";\n}\n'''
(out/'probe.cpp').write_text(source)
cmd=['/opt/gcc-13.2.0/bin/g++','-O1','-g','-std=c++20','-fopenmp','-fsanitize=thread',str(out/'probe.cpp'),'-o',str(out/'probe.bin')]
r=subprocess.run(cmd,capture_output=True,text=True);(out/'compile.log').write_text(r.stdout+r.stderr)
if r.returncode:raise RuntimeError('probe compile failed')
report={'diagnostic_completed':True,'not_a_passing_race_gate':True,'compile_command':cmd,'probe_source_sha256':hashlib.sha256(source.encode()).hexdigest(),'probe_binary_sha256':hashlib.sha256((out/'probe.bin').read_bytes()).hexdigest(),'cases':[]}
for name,setup,solve in [('minimal-openmp',1,1),('serial-solver-parallel-setup',3,1),('serial-all',1,1),('parallel-solver-serial-setup',1,3)]:
 env=os.environ.copy();env.update(TSAN_OPTIONS='halt_on_error=1:exitcode=66',DSPEC_SYMBOLIC_REUSE='1',OMP_NUM_THREADS=str(setup),OMP_THREAD_LIMIT='3',OMP_DYNAMIC='FALSE')
 command=[str(out/'probe.bin')] if name=='minimal-openmp' else [binary,model,str(out/name),'kernel-probe','--lmax','3','--nq','6','--radial-refinement','0.05','--radial-traces','0','--probe-angle-count','8','--solve-threads',str(solve),'--pole-only','on','--forward-only','off','--response-reuse','on']
 with (out/(name+'.stdout')).open('xb') as so,(out/(name+'.stderr')).open('xb') as se:r=subprocess.run(command,env=env,stdout=so,stderr=se,timeout=300)
 err=(out/(name+'.stderr')).read_text();report['cases'].append({'name':name,'command':command,'returncode':r.returncode,'race_reported':'WARNING: ThreadSanitizer: data race' in err,'unsupported_mapping':'unexpected memory mapping' in err,'stderr_excerpt':err[-16000:]})
with (out/'controls.json').open('x') as f:json.dump(report,f,indent=2)
