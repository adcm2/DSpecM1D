#!/usr/bin/env python3
"""Compare completed sections and saved waveforms while slices continue."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.signal import resample
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from check_phase_b import save, sha

p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('baseline_root',type=Path);p.add_argument('geometry',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output.mkdir();degrees=[384,1536];channels=['A','C','F','L','N','kappa']
geometry=json.loads(a.geometry.read_text());receiver=geometry['receiver_along_deg']
data={};raw={};metadata={};angles={};radial={};physical={};traces={}
for degree in degrees:
    folder=a.root/f'section-l{degree}/data';m=json.loads((folder/'section.json').read_text());metadata[degree]=m
    angle=np.genfromtxt(folder/'section_angles.csv',delimiter=',',names=True);angles[degree]=angle
    radial[degree]=np.genfromtxt(folder/'section_radial_sides.csv',delimiter=',',names=True)
    values=np.memmap(folder/'section.bin',dtype='f8',mode='r').reshape(179,m['angles'],6)
    assert np.isfinite(values).all() and m['angles']==4*degree+1
    source=int(np.argmin(np.abs(angle['along_deg'])));assert angle['along_deg'][source]==0
    data[degree]=np.roll(values,-source,axis=1);raw[degree]=values
    base=a.baseline_root/f'global-generated-t32-l{degree}/data'
    physical[degree]=json.loads((base/'result.json').read_text())
    traces[degree]=np.genfromtxt(base/'trace.csv',delimiter=',',names=True)
assert np.array_equal(radial[384],radial[1536])
fine=data[1536];n=fine.shape[1];along=360*np.arange(n)/n
dist=lambda point:np.abs((along-point+180)%360-180)
mask=(dist(0)>=5)&(dist(receiver)>=5)
depth=(6371000-radial[1536]['radius_m'])/1000
interior=(depth>=100)&(depth<=1500)
rows=[]
for c,name in enumerate(channels):
    coarse=resample(data[384][:,:,c],n,axis=1)
    x=coarse[:,mask];y=fine[:,:,c][:,mask]
    def errors(x,y):
        norm=float(np.linalg.norm(y));scale=float(np.max(np.abs(y)));err=x-y
        return {'relative_l2_vs_fine':float(np.linalg.norm(err)/norm) if norm else 0,
                'max_difference_scaled_by_fine_max':float(np.max(np.abs(err))/scale) if scale else 0}
    rows.append({'channel':name,'all_native_sides':errors(x,y),'mantle100_1500km':errors(x[interior],y[interior])})
energies={d:physical[d]['energy_m2_s'] for d in degrees}
t=traces[1536]['time_s'];assert np.array_equal(t,traces[384]['time_s'])
visible=t<=3600;win=(t>=1000)&(t<=1200)
u0=traces[384]['vertical_m'];u1=traces[1536]['vertical_m'];w=traces[1536]['window_weight']
wave_l2=float(np.linalg.norm((u0-u1)[visible])/np.linalg.norm(u1[visible]))
window_l2=float(np.linalg.norm((u0-u1)[win]*np.sqrt(w[win]))/np.linalg.norm(u1[win]*np.sqrt(w[win])))

# One normalization per channel shared between both degrees. This prevents
# independent per-image normalization from hiding amplitude differences.
scales=[]
for c in range(5):
    candidates=[]
    for d in degrees:
        angle=angles[d]['along_deg'];keep=(angle>=-10)&(angle<=receiver+10)&(np.abs(angle)>=5)&(np.abs((angle-receiver+180)%360-180)>=5)
        candidates.append(float(np.max(np.abs(raw[d][interior,:,c][:,keep]))))
    scales.append(max(candidates))
fig,axes=plt.subplots(2,5,figsize=(18,9),sharex=True,sharey=True)
for row,d in enumerate(degrees):
    angle=angles[d]['along_deg'];keep=(angle>=-10)&(angle<=receiver+10);indices=np.flatnonzero(keep);delta=360/len(angle)
    edges=np.r_[angle[indices]-.5*delta,angle[indices[-1]]+.5*delta]
    endpoint=(np.abs(angle)<5)|(np.abs((angle-receiver+180)%360-180)<5)
    for c,ax in enumerate(axes[row]):
        for element in np.unique(radial[d]['element']).astype(int):
            selected=np.flatnonzero(radial[d]['element']==element);r=radial[d]['radius_m'][selected]/1000
            if r[-1]<4871:continue
            re=np.r_[r[-1]-radial[d]['width_m'][selected[-1]]/1000,(r[:-1]+r[1:])/2,r[-1]]
            values=np.array(raw[d][selected,:,c][:,indices],copy=True)/(scales[c] or 1)
            values[:,endpoint[indices]]=np.nan
            ax.pcolormesh(edges,6371-re,values,cmap='RdBu_r',vmin=-1,vmax=1,shading='flat',rasterized=True)
        ax.set(xlim=(-10,receiver+10),ylim=(1500,0));ax.grid(alpha=.15)
        ax.scatter([0],[50],marker='*',c='black',s=35);ax.scatter([receiver],[3],marker='v',c='black',s=25)
        if row==0:ax.set_title(channels[c])
        if c==0:ax.set_ylabel(f'L={d}\nDepth (km)')
        if row==1:ax.set_xlabel('Along great circle (°)')
fig.subplots_adjust(left=.055,right=.89,bottom=.16,top=.88,wspace=.15,hspace=.15)
cax=fig.add_axes([.92,.25,.012,.5]);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1,1),cmap='RdBu_r'),cax=cax,label='Shared normalization per channel; clipped at ±1')
fig.suptitle('DSpecM1D · L=384 versus L=1536 · source–receiver mantle sections\nVertical-displacement energy, 1000–1200 s · 0–10 mHz')
fig.text(.5,.055,'Same scale for both degrees, per channel: max |K| at100–1500km depth outside5° endpoint masks. Shallower peaks clipped.\nRaw data unmasked and unclipped; native SEM interface sides retained. Exploratory; convergence not established.',ha='center',fontsize=9)
for ext in ['png','pdf']:fig.savefig(a.output/f'section-comparison.{ext}',dpi=160)
plt.close(fig)

fig,axes=plt.subplots(2,1,figsize=(11,7))
for d,color in [(384,'#b35a19'),(1536,'#174f77')]:
    u=traces[d]['vertical_m']*1e6
    axes[0].plot(t[visible],u[visible],label=f'L={d}',c=color,lw=.9)
    near=(t>=950)&(t<=1250);axes[1].plot(t[near],u[near],label=f'L={d}',c=color,lw=1.2)
for ax in axes:
    ax.axvspan(1000,1200,color='grey',alpha=.12);ax.set_ylabel('Vertical displacement (µm)');ax.legend(frameon=False);ax.grid(alpha=.2)
axes[0].set_xlim(0,3600);axes[1].set_xlim(950,1250);axes[1].set_xlabel('Time after earthquake (s)')
fig.suptitle(f'Waveform comparison · energy difference {abs(energies[1536]-energies[384])/abs(energies[1536]):.2%} relative to L=1536\nSame filter, damping and Hann measurement window')
fig.tight_layout(rect=[0,0,1,.92])
for ext in ['png','pdf']:fig.savefig(a.output/f'waveform-comparison.{ext}',dpi=160)
plt.close(fig)
save(a.output/'comparison.json',{'status':'DESCRIPTIVE_SECTION_COMPARISON;SLICES_PENDING_AT_LAUNCH','degrees':degrees,'channels':rows,
     'mask':'outside5degree source/receiver wedges; same179 native radial sides',
     'alignment':'periodic Fourier resampling of L384 section onto L1536 angular samples; no radial interpolation',
     'energy_m2_s':energies,'energy_difference_relative_to_fine':abs(energies[1536]-energies[384])/abs(energies[1536]),
     'waveform_relative_l2_0_3600s':wave_l2,'hann_weighted_window_relative_l2':window_l2,
     'comparison_figure_shared_raw_scales':dict(zip(channels[:5],scales)),
     'section_source_sha256':{d:sha(a.root/f'section-l{d}/data/section.bin') for d in degrees},
     'figure_sha256':{f.name:sha(f) for f in a.output.iterdir() if f.suffix in ['.png','.pdf']},
     'limits':['Two degrees do not establish convergence.','Unweighted sample norms are descriptive; not volume-integrated Sobolev norms.',
               'Source-shell contacts and receiver enrichment remain omitted.']})
