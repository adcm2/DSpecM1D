"""Recompute physical energy and refinement differences from saved SI traces."""
from pathlib import Path
import argparse,json,numpy as np
p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('reference',nargs='?');p.add_argument('--output');a=p.parse_args()
def load(name):
 p=Path(name);p=p/'trace.csv' if p.is_dir() else p
 d=np.genfromtxt(p,delimiter=',',names=True)
 t=d['time_s'];u=d['vertical_m'];w=d['window_weight']
 if not np.all(np.isfinite(u)) or not np.all(np.isfinite(w)):raise RuntimeError('non-finite trace')
 if len(t)<2:raise RuntimeError('short trace')
 dt=float(t[1]-t[0]);expect=np.where((t>=1000)&(t<=1200),.5*(1-np.cos(2*np.pi*(t-1000)/200)),0)
 if not np.allclose(w,expect,atol=1e-12,rtol=1e-12):raise RuntimeError('window_weight must mean Hann, not sqrt(Hann)')
 return t,u,w,dt
x,u,w,dt=load(a.candidate);e=.5*dt*float(np.sum(w*u*u));r={'candidate':a.candidate,'energy_m2_s':e,'dt_s':dt,'samples':len(x),'window_peak_m':float(np.max(np.abs(u[(x>=1000)&(x<=1200)])))}
if a.reference:
 y,v,z,ds=load(a.reference);use=(x>=1000)&(x<=1200)
 vr=np.interp(x[use],y,v);ur=u[use];wr=w[use]
 er=.5*ds*float(np.sum(z*v*v));den=float(np.linalg.norm(np.sqrt(wr)*vr))
 if e<=0 or er<=0 or den<=0:raise RuntimeError('noninformative waveform or energy')
 r.update(reference=a.reference,reference_energy_m2_s=er,energy_relative=abs(e-er)/abs(er),window_weighted_relative_L2=float(np.linalg.norm(np.sqrt(wr)*(ur-vr))/den),window_unweighted_relative_L2=float(np.linalg.norm(ur-vr)/np.linalg.norm(vr)),window_peak_scaled_max=float(np.max(np.abs(ur-vr))/np.max(np.abs(vr))))
 r['passed_1e_3']=bool(r['energy_relative']<=1e-3 and r['window_weighted_relative_L2']<=1e-3)
s=json.dumps(r,indent=2)+'\n'
if a.output:Path(a.output).write_text(s)
print(s,end='')
