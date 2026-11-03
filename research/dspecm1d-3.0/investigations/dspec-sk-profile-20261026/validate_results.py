from pathlib import Path
import json,gzip,shutil,tempfile,importlib.util,hashlib,subprocess
root=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
spec=importlib.util.spec_from_file_location('runtime',repo/'benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py')
runtime=importlib.util.module_from_spec(spec);spec.loader.exec_module(runtime)
a=json.loads((root/'records.json').read_text()); results=[]
for degree,folder in [(3,'small-05'),(12,'comparison-illustration-12-02')]:
 for m in (1,2):
    cases={label:next(r for r in a if r['degree']==degree and r['measurements']==m and r['solver']==label and r['repeat']==0 and not r['variant']) for label in ['dspec','yspec','dspec_timed','yspec_timed']}
    instrument={label:cases[label]['csv_sha256']==cases[label+'_timed']['csv_sha256'] for label in ['dspec','yspec']}
    with tempfile.TemporaryDirectory(prefix='dspec-sk-accuracy-') as t:
        files={}
        for label in ['dspec','yspec']:
            files[label]=Path(t)/f'{label}.csv'
            with gzip.open(Path(cases[label]['output'])/'kernels.csv.gz','rb') as f,files[label].open('wb') as g:shutil.copyfileobj(f,g)
        gate=runtime.accuracy_gate(files['dspec'],files['yspec'],repo/'build/kernel-examples-phase12'/folder/'samples.csv',1e-3,m,{'objective_scale_S':cases['dspec']['metrics']['objective_scale_S']})
    results.append({'degree':degree,'measurements':m,'instrumentation_byte_identical':instrument,'gate':gate})
    (root/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
    print(degree,m,instrument,'profile_max',max(r['relative_difference'] for r in gate['profiles']),'integral_max',max(r['relative_difference'] for r in gate['integrals']),flush=True)
# Record current source identity and verify unchanged tracked dirty inputs.
provenance={'scope':'Investigation only; native runs alternate solvers, three independent processes per case, one thread. Timed outputs saved losslessly. No convergence sweep.', 'repos':{}}
for source in [repo,repo.parent/'YSpecpp']:
    provenance['repos'][source.name]=runtime.git_info(source)
provenance['compiler']=subprocess.check_output(['/usr/bin/g++','--version'],text=True)
provenance['uname']=subprocess.check_output(['uname','-a'],text=True)
provenance['executables']={str(p):runtime.sha256(p) for p in (root/'bin').iterdir() if p.is_file()}
provenance['yspec_native_executable']={str(repo/'build/phase-d-build-cache/build-yspecpp/yspec_kernel_timing'):runtime.sha256(repo/'build/phase-d-build-cache/build-yspecpp/yspec_kernel_timing')}
provenance['inputs']={str(repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck'):runtime.sha256(repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck')}
for folder in ['small-05','comparison-illustration-12-02']:
 p=repo/'build/kernel-examples-phase12'/folder/'samples.csv';provenance['inputs'][str(p)]=runtime.sha256(p)
provenance['source_hashes_initial']=json.loads((root/'source_hashes.json').read_text())
provenance['tracked_profile_inputs_unchanged']=all(runtime.sha256(Path(p))==h for p,h in provenance['source_hashes_initial'].items())
provenance['thread_environment']={k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']}
(root/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
