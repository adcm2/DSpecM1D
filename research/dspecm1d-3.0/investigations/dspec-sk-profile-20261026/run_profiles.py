from pathlib import Path
import os,json,subprocess,time,hashlib,gzip,shutil,statistics,csv
root=Path(__file__).resolve().parent
repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
model=repo/'tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck'
yspec=repo/'build/phase-d-build-cache/build-yspecpp/yspec_kernel_timing'
env=os.environ.copy()
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[k]='1'
env['OMP_DYNAMIC']='FALSE'
runs=root/'runs';runs.mkdir(exist_ok=True)
records=json.loads((root/'records.json').read_text()) if (root/'records.json').exists() else []
def run(label,degree,measurements,repeat,samples,suffix='',angles=None):
    for previous in records:
        if (previous['solver'],previous['degree'],previous['measurements'],previous['repeat'],previous['variant'])==(label,degree,measurements,repeat,suffix):return previous
    output=runs/f'l{degree}-m{measurements}-r{repeat}-{label}{suffix}'
    if output.exists():output=output.with_name(output.name+'-retry1')
    output.mkdir()
    exe={'dspec':root/'bin/current','dspec_timed':root/'bin/instrumented','yspec':yspec,'yspec_timed':root/'bin/yspec_instrumented'}[label]
    cmd=[str(exe),str(model),str(samples),str(output),str(degree),str(measurements),'1']
    if angles is not None:cmd.append(str(angles))
    start=time.monotonic()
    p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=env,timeout=120)
    output.with_suffix('.stdout.log').write_text(p.stdout);output.with_suffix('.stderr.log').write_text(p.stderr)
    if p.returncode:raise RuntimeError((cmd,p.stderr[-4000:]))
    meta=json.loads(p.stdout.strip().splitlines()[-1])
    profile=[json.loads(v.split('SKPROFILE ',1)[1]) for v in p.stderr.splitlines() if v.startswith('SKPROFILE ')]
    t=meta['first_workload_api_s'] if label.startswith('yspec') else sum(meta['api_s'])
    f=output/'kernels.csv';digest=hashlib.sha256(f.read_bytes()).hexdigest()
    with f.open('rb') as src,gzip.open(str(f)+'.gz','wb',compresslevel=1) as out:shutil.copyfileobj(src,out)
    f.unlink() # Lossless compression of this newly generated evidence only.
    r={'solver':label,'degree':degree,'measurements':measurements,'repeat':repeat,'variant':suffix,'api_seconds':t,'wall_seconds':time.monotonic()-start,'metrics':meta,'profile':profile,'csv_sha256':digest,'command':cmd,'output':str(output)}
    records.append(r);(root/'records.json').write_text(json.dumps(records,indent=2)+'\n')
    print(output.name,round(t,6),flush=True)
    return r
for degree,folder in [(3,'small-05'),(12,'comparison-illustration-12-02')]:
    samples=repo/'build/kernel-examples-phase12'/folder/'samples.csv'
    for measurements in (1,2):
        for repeat in range(3):
            for label in (('dspec','yspec') if repeat%2==0 else ('yspec','dspec')):
                run(label,degree,measurements,repeat,samples)
        for label in ('dspec_timed','yspec_timed'):
            run(label,degree,measurements,0,samples)
# Fixed lmax12 diagnostic samples change output size, not numerical resolution.
samples=repo/'build/kernel-examples-phase12/comparison-illustration-12-02/samples.csv'
with samples.open() as f:
    reader=csv.DictReader(f);fields=reader.fieldnames
    a=[row for row in reader if row['modulus']=='A' and row['objective']=='J_Z']
first_radius=a[0]['radius_m']; theta=a[0]['theta_rad']; phi=a[0]['phi_rad']
three={(r['theta_rad'],r['phi_rad']) for r in a[:3]}
for variant,selected,angles in [('one_radius',[r for r in a if r['radius_m']==first_radius],1352),('one_angle',[r for r in a if r['theta_rad']==theta and r['phi_rad']==phi],1),('one_point',[a[0]],1),('three_angles',[r for r in a if (r['theta_rad'],r['phi_rad']) in three],3),('three_points',a[:3],3)]:
    p=root/f'{variant}.csv'
    with p.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(selected)
    for label in ('dspec','yspec','dspec_timed','yspec_timed'):
        if label.startswith('yspec') and angles<3:continue
        run(label,12,1,0,p,'-'+variant,angles if label.startswith('yspec') else None)
summary=[]
for degree in (3,12):
 for measurements in (1,2):
    row={'degree':degree,'measurements':measurements}
    for solver in ('dspec','yspec'):
        times=[r['api_seconds'] for r in records if r['solver']==solver and r['degree']==degree and r['measurements']==measurements and not r['variant']]
        row[solver]={'median_s':statistics.median(times),'min_s':min(times),'max_s':max(times)}
    row['dspec_over_yspec']=row['dspec']['median_s']/row['yspec']['median_s'];summary.append(row)
(root/'timing_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
