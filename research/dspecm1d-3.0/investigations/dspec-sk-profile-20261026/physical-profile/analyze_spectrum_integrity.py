"""Independent NumPy inverse FFT of SI spectra for adapter normalization checks."""
from pathlib import Path
import argparse,json,numpy as np
p=argparse.ArgumentParser();p.add_argument('directory');p.add_argument('--filtered',action='store_true');p.add_argument('--nt',type=int,default=1024);p.add_argument('--epsilon',type=float,default=5/3600);p.add_argument('--output');a=p.parse_args()
r=Path(a.directory);s=np.genfromtxt(r/'spectrum.csv',delimiter=',',names=True);t=np.genfromtxt(r/'trace.csv',delimiter=',',names=True)
z=np.zeros(a.nt//2+1,dtype=complex);k=np.array(s['bin'],dtype=int);v=s['real']+1j*s['imag'];f=s['frequency_hz'];df=float(np.diff(f)[0]);dt=1/(a.nt*df)
if not a.filtered:
 b=np.where(f<.008,1,np.where(f<.01,.5*(1+np.cos(np.pi*(f-.008)/.002)),0));v=v*b
z[k]=v;u=np.fft.irfft(z,n=a.nt)*a.nt*df*np.exp(a.epsilon*np.arange(a.nt)*dt)
ix=np.rint(t['time_s']/dt).astype(int);sel=(t['time_s']>=1000)&(t['time_s']<=1200);ref=u[ix[sel]];cand=t['vertical_m'][sel]
r={'spectrum_directory':str(r),'frequency_spacing_hz':df,'dt_s':dt,'window_trace_spectrum_relative_L2':float(np.linalg.norm(cand-ref)/np.linalg.norm(ref)),'independent_energy_m2_s':float(.5*dt*np.sum(t['window_weight'][sel]*ref**2)),'adapter_energy_m2_s':float(.5*dt*np.sum(t['window_weight'][sel]*cand**2))}
r['passed']=r['window_trace_spectrum_relative_L2']<1e-10
j=json.dumps(r,indent=2)+'\n'
if a.output:Path(a.output).write_text(j)
print(j,end='')
