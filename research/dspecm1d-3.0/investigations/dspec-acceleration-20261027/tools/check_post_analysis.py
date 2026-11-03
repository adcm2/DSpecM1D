import json,pathlib,subprocess,sys,hashlib
sys.path.insert(0,sys.argv[1]);import check_phase_a_summary as summary
output=pathlib.Path(sys.argv[2]);output.mkdir()
roots=json.loads(sys.argv[3]);results={}
for name,root in roots.items():
 value=summary.summarize(pathlib.Path(root))
 with (output/(name+'.json')).open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False)
 results[name]={'status':value['status'],'state_sha256':value['state_sha256'],'job_count':len(value['jobs'])}

import check_pair
repeat_root=pathlib.Path(sys.argv[7]);repeats={}
for degree in (128,384):
 for version in ('v4','v5'):
  name=f'paired-{version}-l{degree}';receipt=check_pair.verified_job(repeat_root,name);job=receipt[0];folder=repeat_root/name
  repeats[name]={'process_wall_s':job['elapsed_seconds'],'peak_threads':job.get('peak_observed_threads'),'peak_rss_gib_userspace':job.get('peak_observed_rss_gib_userspace'),'result':check_pair.load(folder/'data/result.json'),'profile':check_pair.profile(folder/job['log']),'verified_artifact_sha256':receipt[1]}
  after=check_pair.verified_job(repeat_root,name)
  assert receipt[0:2]==after[0:2],name
with (output/'paired_repeats.json').open('x') as stream:json.dump(repeats,stream,indent=2,allow_nan=False)

resolution=output/'resolution.json'
command=[sys.executable,str(pathlib.Path(sys.argv[1])/'check_global_resolution.py'),'--l128-root',sys.argv[4],'--l128-job','global-generated-w1-l128','--l384-root',sys.argv[4],'--l384-job','global-generated-w1-l384','--layout-root',sys.argv[5],'--theta128-job','reference-l128','--theta384-job','reference-l384','--output',str(resolution)]
run=subprocess.run(command,check=False)
old=pathlib.Path(sys.argv[6]);data={'queue_statuses':results,'resolution_returncode':run.returncode,'resolution_status':json.loads(resolution.read_text())['status'] if resolution.exists() else 'missing','prior_yspec_summary':{'path':str(old),'sha256':hashlib.sha256(old.read_bytes()).hexdigest()},'claim_limits':['Bookkeeping completion is not scientific acceptance.','Original-global 2e-13 gate remains failed; internal thread/storage parity is diagnostic.','No figures or 179-side global volume generated in this phase.','Resolution comparison uses approximate angular interpolation, not a convergence certificate.']}
with (output/'summary.json').open('x') as stream:json.dump(data,stream,indent=2,allow_nan=False)
