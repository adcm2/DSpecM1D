"""Fresh approved bounded baseline/Born checks; no finite-difference sweep."""
from pathlib import Path
import hashlib,importlib.util,json,os,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parent;REPO=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D');HERE=REPO/'benchmarks/yspecpp/kernels/3d'
ENV=os.environ.copy()
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):ENV[k]='1'
ENV['MPLCONFIGDIR']=str(ROOT/'mpl-cache')
OUT=ROOT/sys.argv[1]
if OUT.exists() and any(OUT.iterdir()):raise RuntimeError('Refusing nonempty evidence output')
OUT.mkdir(exist_ok=True)
records=[]
def run(cmd,log,cwd=None):
 start=time.monotonic()
 with (OUT/log).open('w') as f:p=subprocess.run(list(map(str,cmd)),cwd=cwd,env=ENV,stdout=f,stderr=subprocess.STDOUT,timeout=300)
 record={'command':list(map(str,cmd)),'cwd':str(cwd) if cwd else os.getcwd(),'wall_s':time.monotonic()-start,'exit_status':p.returncode};records.append(record);(OUT/'commands.json').write_text(json.dumps(records,indent=2)+'\n')
 if p.returncode:raise RuntimeError((log,(OUT/log).read_text()[-5000:]))
 print(log,record['wall_s'],flush=True)
spec=importlib.util.spec_from_file_location('phase_d',HERE/'run_phase_d.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
shutil.copyfile(HERE/'fixture.json',OUT/'fixture.json');module.write_deck(OUT/'common_homogeneous.deck')
base=json.loads((ROOT.parent/'build_commands.json').read_text())['current']
base[base.index(str(ROOT.parent/'timing_current.cpp'))]=str(HERE/'dspec_3d_adapter.cpp');base[base.index('-o')+1]=str(ROOT/'bin/after-adapter')
run(base,'build-dspec.log')
for name in ['dspec','yspecpp','oracle']:(OUT/name).mkdir()
model=OUT/'common_homogeneous.deck'
run([ROOT/'bin/after-adapter',model,OUT/'dspec',.05,8],'dspec.log')
yexe=REPO/'build/phase-d-build-cache/build-yspecpp/phase_d_adapter'
run([yexe,model,OUT/'yspecpp',OUT/'dspec/radial_nodes.csv',OUT/'dspec/angles.csv'],'yspec.log')
run([sys.executable,HERE/'cartesian_ti_oracle.py','--baseline',OUT/'dspec','--output',OUT/'oracle','--fixture',OUT/'fixture.json'],'oracle.log')
metrics=module.baseline_metrics(OUT);(OUT/'baseline_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source_paths=[REPO/'DSpecM1D/src/Windowed3DFields.h',REPO/'DSpecM1D/src/Windowed3DKernels.h',REPO/'DSpecM1D/src/SensitivityKernels.h',HERE/'dspec_3d_adapter.cpp',HERE/'yspec_3d_adapter.cpp',HERE/'cartesian_ti_oracle.py',HERE/'fixture.json']
provenance={'kind':'Fresh bounded baseline and Born calculation; no finite-differences. Same frozen fixture, single numerical thread. Existing native YSpec adapter reused; freshly compiled DSpec adapter and freshly recalculated Cartesian oracle.','runtime_threads':{k:ENV[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS')},'sources':{str(p):sha(p) for p in source_paths},'executables':{str(p):sha(p) for p in [ROOT/'bin/after-adapter',yexe]},'all_baseline_gates_pass':metrics['all_baseline_gates_pass']}
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');print('BOUNDED BORN',metrics['all_baseline_gates_pass'],flush=True)
assert metrics['all_baseline_gates_pass']
