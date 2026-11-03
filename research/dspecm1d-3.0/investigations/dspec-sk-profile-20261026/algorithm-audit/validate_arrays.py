"""Independent parent parser and frozen numerical metrics for actual outputs."""
from pathlib import Path
import struct,json,numpy as np

def evaluations(path):
 data=Path(path).read_bytes();offset=0;out=[]
 def scalar(code):
  nonlocal offset
  result=struct.unpack_from('<'+code,data,offset)[0];offset+=struct.calcsize('<'+code);return result
 def matrix():
  nonlocal offset
  r=scalar('I');c=scalar('I');a=np.frombuffer(data,dtype='<c16',count=r*c,offset=offset).reshape(r,c).copy();offset+=16*r*c;return a
 def real_array(n,cols):
  nonlocal offset
  a=np.frombuffer(data,dtype='<f8',count=n*cols,offset=offset).reshape(n,cols).copy();offset+=8*n*cols;return a
 while offset<len(data):
  raw=matrix();processed=matrix();kernel=real_array(scalar('Q'),6);contacts=real_array(scalar('Q'),5);obj=scalar('d')
  for a in [raw,processed,kernel,contacts,np.array(obj)]:assert np.isfinite(a).all(),str(path)
  out.append({'raw':raw,'processed':processed,'kernels':kernel,'contacts':contacts,'objective':obj})
 assert offset==len(data);return out

def metric(a,b,l2limit=1e-12,pointlimit=1e-10):
 assert a.shape==b.shape
 den=max(float(np.linalg.norm(a)),float(np.linalg.norm(b)),1e-300)
 l2=float(np.linalg.norm(a-b)/den)
 if a.size:
  peak=max(float(np.max(abs(a))),float(np.max(abs(b))),1e-300)
  point=float(np.max(abs(a-b)/np.maximum(np.maximum(abs(a),abs(b)),.01*peak)))
 else:point=0.0
 return {'relative_L2':l2,'max_scaled_pointwise':point,'L2_limit':l2limit,'pointwise_limit':pointlimit if np.isfinite(pointlimit) else None,'passed':l2<=l2limit and point<=pointlimit}

def compare(candidate,reference):
 a=evaluations(candidate);b=evaluations(reference);assert len(a)==len(b)>0
 result=[]
 for i,(x,y) in enumerate(zip(a,b)):
  checks={name:metric(x[name],y[name],1e-12,float('inf')) for name in ['raw','processed']}
  for j,name in enumerate(['A','C','F','L','N','kappa']):checks['kernel_'+name]=metric(x['kernels'][:,j],y['kernels'][:,j])
  for j,name in enumerate(['A','C','F','L','N']):checks['contact_'+name]=metric(x['contacts'][:,j],y['contacts'][:,j],1e-12,1e-12)
  checks['objective']=metric(np.array(x['objective']),np.array(y['objective']),1e-12,float('inf'))
  result.append({'objective':i,'checks':checks,'passed':all(c['passed'] for c in checks.values())})
 return {'objectives':result,'passed':all(r['passed'] for r in result)}
if __name__=='__main__':
 import sys
 r=compare(sys.argv[1],sys.argv[2]);print(json.dumps(r,indent=2));raise SystemExit(0 if r['passed'] else 1)
