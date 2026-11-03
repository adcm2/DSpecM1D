from pathlib import Path
import json,gzip,numpy as np
from validate_arrays import metric
R=Path(__file__).resolve().parent
refdir=R.parent/'phase1/attempt-02/field-run'
def load(p):
 with gzip.open(p,'rt') as f:return np.loadtxt(f,delimiter=',')
a=load(refdir/'fields.csv.gz');ca=load(refdir/'contacts.csv.gz')
groups={'forward_raw_displacement':(0,3),'observable_SI_displacement':(3,6),'adjoint_dual_displacement':(6,9),'forward_raw_strain':(9,15),'adjoint_dual_strain':(15,21),'forward_rscaled_strain':(21,27),'adjoint_rscaled_strain':(27,33)}
results=[]
for rec in json.loads((R/'field-validation-v5-tensor/records.json').read_text()):
 out=Path(rec['output']);b=load(out/'fields.csv.gz');cb=load(out/'contacts.csv.gz');assert a.shape==b.shape and np.array_equal(a[:,:5],b[:,:5]);assert ca.shape==cb.shape and np.array_equal(ca[:,:3],cb[:,:3])
 x,y=a[:,5::2]+1j*a[:,6::2],b[:,5::2]+1j*b[:,6::2];checks=[];contacts=[]
 for case,freq in np.unique(a[:,:2],axis=0):
  mask=(a[:,0]==case)&(a[:,1]==freq)
  for key,(lo,hi) in groups.items():
   xx,yy=x[mask,lo:hi],y[mask,lo:hi]
   err=float(np.linalg.norm(xx-yy)/max(np.linalg.norm(xx),np.linalg.norm(yy),1e-30));limit=1e-9 if key.startswith('adjoint') and 'strain' in key else 1e-12
   checks.append({'case':int(case),'frequency':int(freq),'field':key,'relative_norm':err,'limit':limit,'passed':err<=limit})
  cmask=(ca[:,0]==case)&(ca[:,1]==freq)
  for i,mod in enumerate('ACFLN',start=3):contacts.append({'case':int(case),'frequency':int(freq),'modulus':mod,**metric(ca[cmask,i],cb[cmask,i],1e-12,1e-12)})
 result={'backend':rec['backend'],'bridge':rec['bridge'],'field_checks':checks,'contacts':contacts,'packed_checks':rec['internal_packed_passed'],'passed':all(x['passed'] for x in checks+contacts) and rec['internal_packed_passed']}
 results.append(result)
 print(rec['backend'],rec['bridge'],result['passed'],'maxfield',max(x['relative_norm'] for x in checks),'maxcontactL2',max(x['relative_L2'] for x in contacts),flush=True)
(R/'fields_v5_tensor_vs_frozen.json').write_text(json.dumps(results,indent=2)+'\n')
summary={'comparisons':len(results),'field_groups':sum(len(x['field_checks']) for x in results),'contact_groups':sum(len(x['contacts']) for x in results),'passed':all(x['passed'] for x in results),'max_field_norm':max(c['relative_norm'] for x in results for c in x['field_checks']),'max_contact_L2':max(c['relative_L2'] for x in results for c in x['contacts']),'max_contact_pointwise':max(c['max_scaled_pointwise'] for x in results for c in x['contacts'])}
(R/'fields_v5_tensor_vs_frozen_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
raise SystemExit(0 if summary['passed'] else 1)
