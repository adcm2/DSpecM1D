#!/usr/bin/env python3
"""Render reproducible exploratory kernel sections/slices from bounded outputs."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, SymLogNorm
from check_phase_b import sha, save

CHANNELS=['A','C','F','L','N','kappa']
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'figure.facecolor':'white',
                     'savefig.facecolor':'white','pdf.fonttype':42})

def save_figure(fig, out, name):
    fig.savefig(out/(name+'.png'),dpi=160)
    fig.savefig(out/(name+'.pdf'),dpi=160)
    plt.close(fig)

def angular_distance(lat, lon, lat0, lon0):
    a=np.deg2rad(lat);a0=np.deg2rad(lat0)
    return np.rad2deg(np.arccos(np.clip(np.sin(a)*np.sin(a0)+np.cos(a)*np.cos(a0)*np.cos(np.deg2rad(lon-lon0)),-1,1)))

def render(section, slices, baseline, geometry, checks, output):
    assert json.loads(checks.read_text())['status']=='PASS'
    output.mkdir()
    sm=json.loads((section/'section.json').read_text());dm=json.loads((slices/'slices.bin.json').read_text())
    physical=json.loads((baseline/'result.json').read_text());geom=json.loads(geometry.read_text())
    degree=sm['degree'];assert degree==dm['degree']==physical['degree']
    angles=np.genfromtxt(section/'section_angles.csv',delimiter=',',names=True)
    radial=np.genfromtxt(section/'section_radial_sides.csv',delimiter=',',names=True)
    data=np.memmap(section/'section.bin',dtype='f8',mode='r').reshape(179,len(angles),6)
    delta=360/len(angles);edges=np.deg2rad(np.r_[angles['along_deg']-.5*delta,angles['along_deg'][-1]+.5*delta])
    recv=geom['receiver_along_deg']
    def along_dist(point):return np.abs((angles['along_deg']-point+180)%360-180)
    mask=(along_dist(0)<5)|(along_dist(recv)<5)
    scales={}
    norm=SymLogNorm(linthresh=1e-3,linscale=.5,vmin=-1,vmax=1,base=10)
    fig,axes=plt.subplots(2,3,subplot_kw={'projection':'polar'},figsize=(15,10.7))
    for c,ax in enumerate(axes.flat):
        scale=float(np.max(np.abs(data[:,~mask,c])))
        scales[CHANNELS[c]]=scale
        for element in np.unique(radial['element']).astype(int):
            rows=np.flatnonzero(radial['element']==element)
            radii=radial['radius_m'][rows]/1000
            lower=radii[-1]-radial['width_m'][rows[-1]]/1000
            redges=np.r_[lower,(radii[:-1]+radii[1:])/2,radii[-1]]
            values=np.array(data[rows,:,c],copy=True)/(scale if scale else 1)
            values[:,mask]=np.nan
            im=ax.pcolormesh(edges,redges,values,cmap='RdBu_r',norm=norm,shading='flat',rasterized=True)
        ax.set_theta_zero_location('N');ax.set_theta_direction(-1)
        ax.set_ylim(0,6371);ax.set_yticks([1221.5,3480,5701,6371]);ax.set_yticklabels(['ICB','CMB','670 km','surface'],fontsize=7)
        ax.grid(alpha=.3);ax.set_title(CHANNELS[c])
        ax.scatter([0],[6321],marker='*',s=65,c='black',zorder=5)
        ax.scatter([np.deg2rad(recv)],[6368],marker='v',s=30,c='black',zorder=5)
    cax=fig.add_axes([.91,.25,.015,.5]);fig.colorbar(im,cax=cax,label='K / max |K|; symmetric log, linear below 10⁻³')
    fig.suptitle(f'DSpecM1D · L={degree} · source–receiver great-circle section\nVertical-displacement energy, 1000–1200 s · 0–10 mHz',y=.98)
    fig.text(.5,.015,'Exploratory; degree convergence not established. Black star: source; triangle: receiver.\n5° angular display masks; raw data retained. Native SEM interface sides separate. Source-shell contacts omitted.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.89,bottom=.08,wspace=.32,hspace=.3,right=.86)
    save_figure(fig,output,f'section-l{degree}')

    zoom=(angles['along_deg']>=-10)&(angles['along_deg']<=recv+10)
    zi=np.flatnonzero(zoom)
    zedges=np.r_[angles['along_deg'][zi]-.5*delta,angles['along_deg'][zi[-1]]+.5*delta]
    zscales={}
    fig,axes=plt.subplots(1,5,figsize=(18,5.5),sharex=True,sharey=True)
    interior=(radial['radius_m']>=4871000)&(radial['radius_m']<=6271000)
    for c,ax in enumerate(axes):
        zscale=float(np.max(np.abs(data[interior,:,c][:,zoom&~mask])))
        zscales[CHANNELS[c]]=zscale
        for element in np.unique(radial['element']).astype(int):
            rows=np.flatnonzero(radial['element']==element)
            radii=radial['radius_m'][rows]/1000
            if radii[-1]<4871:continue
            lower=radii[-1]-radial['width_m'][rows[-1]]/1000
            redges=np.r_[lower,(radii[:-1]+radii[1:])/2,radii[-1]]
            depths=6371-redges
            v=np.array(data[rows,:,c][:,zi],copy=True)/(zscale or 1)
            v[:,mask[zi]]=np.nan
            ax.pcolormesh(zedges,depths,v,cmap='RdBu_r',vmin=-1,vmax=1,shading='flat',rasterized=True)
        ax.set(xlim=(-10,recv+10),ylim=(1500,0),xlabel='Along great circle (°)',title=CHANNELS[c])
        ax.axvline(0,color='black',lw=.6,ls=':');ax.axvline(recv,color='black',lw=.6,ls=':');ax.grid(alpha=.15)
        ax.scatter([0],[50],marker='*',s=40,c='black',zorder=4);ax.scatter([recv],[3],marker='v',s=25,c='black',zorder=4)
    axes[0].set_ylabel('Depth (km)')
    fig.suptitle(f'DSpecM1D · L={degree} · source–receiver mantle section\nVertical-displacement energy, 1000–1200 s · 0–10 mHz',y=.98)
    fig.subplots_adjust(left=.05,right=.89,bottom=.24,top=.8,wspace=.15)
    cax=fig.add_axes([.92,.3,.012,.45]);fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1,1),cmap='RdBu_r'),cax=cax,label='K / interior max |K|; clipped at ±1')
    fig.text(.5,.07,'Exploratory; not degree-converged. Linear scale uses max |K| at100–1500km depth outside5° endpoint masks; shallower peaks clipped.\nRaw unmasked, unclipped data retained. Each SEM element shown separately. Source-shell contacts omitted.',ha='center',fontsize=9)
    save_figure(fig,output,f'section-zoom-l{degree}')

    # The raw slice uses the complete bandwidth grid. Overview rendering takes
    # native angular samples at a stride set by plot pixel dimensions. This
    # changes only the displayed raster, never the saved numerical output.
    values=np.memmap(slices/'slices.bin',dtype='f8',mode='r').reshape(3,dm['n_theta'],dm['n_phi'],6)
    nt,np_=dm['n_theta'],dm['n_phi']
    # Existing native metadata has exactly the same Gauss-Legendre theta rule.
    q=np.genfromtxt(baseline/'angular_grid.csv',delimiter=',',names=True) if (baseline/'angular_grid.csv').exists() else None
    if q is None:
        from scipy.special import roots_legendre
        nodes=roots_legendre(nt)[0]
        theta=np.arccos(-nodes)
    else:
        theta=q['theta_rad']
    assert len(theta)==nt and np.all(np.diff(theta)>0)
    st=max(1,int(np.ceil(nt/900)));sp=max(1,int(np.ceil(np_/1800)))
    ti=np.arange(0,nt,st);pi=np.arange(0,np_,sp)
    lat=90-np.rad2deg(theta[ti]);lon=360*pi/np_
    lon[lon>=180]-=360;order=np.argsort(lon);lon=lon[order]
    lon2,lat2=np.meshgrid(lon,lat)
    endpoint=(angular_distance(lat2,lon2,10,20)<5)|(angular_distance(lat2,lon2,45,75)<5)
    displayed=np.array(values[:,ti,:,:][:,:,pi,:][:,:,order,:],copy=True)
    scale_by_channel=[]
    for c in range(5):scale_by_channel.append(float(np.max(np.abs(displayed[:,:,:,c][:,~endpoint]))))
    fig,axes=plt.subplots(3,5,figsize=(19,10.5),sharex=True,sharey=True)
    actual_depths=[]
    for row,metadata in enumerate(dm['native_side_metadata']):
        depth=(6371000-metadata['radius_m'])/1000;actual_depths.append(depth)
        for c in range(5):
            ax=axes[row,c];v=displayed[row,:,:,c]/(scale_by_channel[c] or 1)
            v=np.where(endpoint,np.nan,v)
            ax.pcolormesh(lon,lat,v,shading='nearest',cmap='RdBu_r',norm=norm,rasterized=True)
            ax.scatter([20],[10],marker='*',s=40,c='black');ax.scatter([75],[45],marker='v',s=20,c='black')
            ax.set_xlim(-180,180);ax.set_ylim(-90,90);ax.set_xticks([-180,-90,0,90,180]);ax.set_yticks([-90,-45,0,45,90]);ax.grid(alpha=.15)
            if row==0:ax.set_title(CHANNELS[c])
            if c==0:ax.set_ylabel(f'{depth:.2f} km depth\nLatitude (°)')
            if row==2:ax.set_xlabel('Longitude (°)')
    fig.suptitle(f'DSpecM1D · L={degree} · native-depth slices\nVertical-displacement energy, 1000–1200 s · 0–10 mHz',y=.98)
    fig.subplots_adjust(left=.055,right=.88,top=.89,bottom=.12,wspace=.13,hspace=.22)
    cax=fig.add_axes([.91,.25,.01,.5]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='RdBu_r'),cax=cax,label='K / displayed max |K|; symmetric log; shared depths')
    fig.text(.5,.035,f'Exploratory; not degree-converged. Solid channels A,C,F,L,N; κ=0 in these mantle slices.\n5° source/receiver display masks; full-resolution unmasked data retained. Overview stride θ={st}, φ={sp}; source-shell contacts omitted.',ha='center',fontsize=9)
    save_figure(fig,output,f'depth-slices-l{degree}')

    trace=np.genfromtxt(baseline/'trace.csv',delimiter=',',names=True)
    fig,axes=plt.subplots(2,1,figsize=(10,7),gridspec_kw={'height_ratios':[2,1]})
    t=trace['time_s'];u=trace['vertical_m'];w=trace['window_weight']
    visible=t<=3600
    axes[0].plot(t[visible],u[visible]*1e6,lw=.8,color='#174f77');axes[0].axvspan(1000,1200,color='orange',alpha=.2,label='1000–1200 s window')
    axes[0].set_xlim(0,3600);axes[0].set_ylabel('Vertical displacement (µm)');axes[0].legend(frameon=False);axes[0].grid(alpha=.2)
    f=np.linspace(0,.012,1000);band=np.where(f<.008,1,np.where(f<.010,.5*(1+np.cos(np.pi*(f-.008)/.002)),0))
    axes[1].plot(f*1000,band,color='#174f77');axes[1].set(xlabel='Frequency (mHz)',ylabel='Filter amplitude',xlim=(0,12),ylim=(-.05,1.05));axes[1].grid(alpha=.2)
    energy=float(.5*physical['dt_s']*np.sum(w*u*u))
    assert abs(energy/physical['energy_m2_s']-1)<1e-12
    fig.suptitle(f'DSpecM1D · L={degree} · vertical-displacement energy\nE = {energy:.8e} m² s; Hann weight over 1000–1200 s')
    fig.tight_layout(rect=[0,.015,1,.92]);save_figure(fig,output,f'waveform-filter-l{degree}')
    files=sorted(p for p in output.iterdir() if p.suffix in ['.png','.pdf'])
    save(output/'figure-manifest.json',{'degree':degree,'status':'EXPLORATORY_FIGURES',
         'full_global_volume':False,'section_scales_raw':scales,'section_zoom_scales_raw':zscales,'section_overview_color_norm':'SymLogNorm;linthresh1e-3;normalized per channel','section_zoom_color_norm':'linear;interior100-1500km max;clipped+-1','slice_scales_raw':dict(zip(CHANNELS[:5],scale_by_channel)),
         'slice_color_norm':'SymLogNorm;linthresh1e-3;shared per channel across depths','slice_depths_km':actual_depths,'slice_overview_native_stride':[st,sp],
         'raw_numerical_grid':'complete kernel-bandwidth grid for three slices;4L+1 section samples;179native sides',
         'energy_recomputed_from_saved_trace':energy,'endpoint_display_mask_deg':5,
         'files_sha256':{p.name:sha(p) for p in files},'validation_receipt_sha256':sha(checks),
         'normalization':'per channel; section across all unmasked native points; slices shared across3depths using displayed native samples',
         'scientific_limits':['Not degree-converged.','No source-shell contact or receiver/Sobolev terms.','No new cross-code comparison.']})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('paths',nargs=6,type=Path);a=p.parse_args();render(*a.paths)
