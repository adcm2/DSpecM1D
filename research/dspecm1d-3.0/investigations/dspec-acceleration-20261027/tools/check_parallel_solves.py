#!/usr/bin/env python3
"""Verify scientific bytes and separate wall/worker profiling accounting."""
from pathlib import Path
import argparse,json,math,subprocess,sys
import check_pair as cp

def parallel_profile(path):
 records=[json.loads(x.split('MATHPROFILE ',1)[1]) for x in path.read_text().splitlines() if 'MATHPROFILE ' in x]
 if not records:return None
 if len(records)!=1:raise ValueError('duplicate profile record')
 r=records[0];cp.finite(r)
 stages=r['stages'];total=stages['physical_case_total']['inclusive_s']
 ordinary={n:v for n,v in stages.items() if not n.startswith('parallel_worker.')}
 exclusive=sum(v['exclusive_s'] for v in ordinary.values())
 if abs(exclusive-total)>max(1e-8,total*2e-12):raise ValueError('parent wall-time closure failed')
 for v in stages.values():
  if any(v[k]<0 for k in ('inclusive_s','exclusive_s','calls')):raise ValueError('negative profile value')
 r['parallel_worker_durations_are_summed_elapsed_not_wall']=True
 return r

def compare(a):
 cp.profile=parallel_profile
 if not a.spool:
  r=cp.assess(Path(a.reference_root),a.reference_job,Path(a.candidate_root),a.candidate_job,True)
 else:
  roots=[Path(a.reference_root),Path(a.candidate_root)];ids=[a.reference_job,a.candidate_job]
  rec=[cp.verified_job(p,i) for p,i in zip(roots,ids)];ds=[p/i/'data' for p,i in zip(roots,ids)]
  values=[cp.load(p/'result.json') for p in ds]
  for key in cp.PHYSICAL:
   if values[0][key]!=values[1][key]:raise ValueError('physical mismatch '+key)
  for v in values:
   if not(v['response_record_spool_used'] and v['energy_m2_s']>0 and v['energy_adjoint_relative_error_vs_2E']<1e-10):raise ValueError('invalid spool energy')
  exact={}
  for name in ['trace.csv','spectrum.csv','radial_nodes.csv','pole_response_spool.bin','kernel.bin']:
   for x in rec:
    if 'data/'+name not in x[1]:raise ValueError('unreceipted file '+name)
   exact[name]=cp.sha(ds[0]/name)==cp.sha(ds[1]/name)
  r={'status':'PASS' if all(exact.values()) else 'FAIL','byte_exact_files':exact,'profiles':[parallel_profile(p/i/x[0]['log']) for p,i,x in zip(roots,ids,rec)],'telemetry':[x[0] for x in rec]}
  if values[0]['energy_m2_s']!=values[1]['energy_m2_s']:raise ValueError('energy mismatch')
  for p,i in zip(roots,ids):cp.verified_job(p,i)
 r['reference_job']=a.reference_job;r['candidate_job']=a.candidate_job
 return r

def main():
 p=argparse.ArgumentParser()
 for n in ['reference-root','reference-job','candidate-root','candidate-job','output']:p.add_argument('--'+n,required=True)
 p.add_argument('--spool',action='store_true');a=p.parse_args()
 try:r=compare(a)
 except Exception as e:r={'status':'FAIL','error':str(e),'type':type(e).__name__}
 with Path(a.output).open('x') as f:json.dump(r,f,indent=2,allow_nan=False)
 return 0 if r['status']=='PASS' else 1
if __name__=='__main__':sys.exit(main())
