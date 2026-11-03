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

def render_section(section,baseline,geometry,output):
    output.mkdir()
    sm=json.loads((section/'section.json').read_text())
    physical=json.loads((baseline/'result.json').read_text());geom=json.loads(geometry.read_text())
    degree=sm['degree'];assert degree==physical['degree']
    assert sm['radial_sides']==179 and sm['angles']==4*degree+1 and not sm['full_global_volume']
    angles=np.genfromtxt(section/'section_angles.csv',delimiter=',',names=True)
    radial=np.genfromtxt(section/'section_radial_sides.csv',delimiter=',',names=True)
    data=np.memmap(section/'section.bin',dtype='f8',mode='r').reshape(179,len(angles),6)
    assert np.isfinite(data).all()
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

    save(output/'section-figure-manifest.json',{'degree':degree,'section_only':True,'scales':scales,'zoom_scales':zscales,'raw_section_sha256':sha(section/'section.bin'),'files_sha256':{p.name:sha(p) for p in output.iterdir() if p.suffix in ['.png','.pdf']}})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('paths',nargs=4,type=Path);a=p.parse_args();render_section(*a.paths)
