#!/usr/bin/env python3
"""Validate and render matched documentation-case native-side kernels."""
from pathlib import Path
import argparse,json,hashlib,sys
import numpy as np
import check_pair as cp

def data(root,name):
 job,hashes,state=cp.verified_job(root,name);p=root/name/'data';r=json.loads((p/'result.json').read_text())
 assert r['attenuation'] and abs(r['time_shift_adjoint_identity']-1)<1e-10 and r['gravity_max_abs_nd_error']<1e-10
 return p,r,job

def table(path):return np.atleast_1d(np.genfromtxt(path,names=True,delimiter=','))
def dump(path,value):
 with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False)
def parity(root,first,second,out):
 p,r,j=data(root,first);q,s,k=data(root,second)
 exact={name:cp.sha(p/name)==cp.sha(q/name) for name in ['trace.csv','spectrum.csv','components.csv','nodes.csv','along.bin','along-angles.csv']}
 assert all(exact.values()),exact
 dump(out,{'status':'PASS','byte_exact_files':exact,'jobs':[first,second]})
def fd(root,base,plus,minus,out):
 p,r,_=data(root,base);q,pr,_=data(root,plus);s,mr,_=data(root,minus)
 assert pr['epsilon_perturbation']==-mr['epsilon_perturbation']>0
 t=table(p/'trace.csv');tp=table(q/'trace.csv');tm=table(s/'trace.csv');dt=t['time_s'][1]-t['time_s'][0]
 finite=dt*np.dot(t['adjoint_source'],(tp['velocity_m_s']-tm['velocity_m_s'])/(2*pr['epsilon_perturbation']))
 predicted=r['integrated_profile_kernel'];error=abs(finite-predicted)/max(abs(finite),abs(predicted))
 report={'status':'PASS' if error<2e-4 else 'FAIL','finite_difference_seconds':float(finite),'integrated_kernel_seconds':predicted,'relative_error':float(error),'limit':2e-4,'perturbation':'ln(vs) smooth sin^2 bump, radius5300-5800km; fixedvp/rho/effectiveQ; source/receiver support excluded','jobs':[base,plus,minus]}
 dump(out,report);assert report['status']=='PASS',report

def pair(root,ds,ys):
 p,r,j=data(root,ds);q,s,k=data(root,ys)
 for key in ['degree','nq','refinement','attenuation','frequency_count','first_bin','last_bin','fft_samples','fft_period_s','native_side_count','along_angles','across_angles']:assert r[key]==s[key],key
 assert cp.sha(p/'nodes.csv')==cp.sha(q/'nodes.csv')
 nodes=table(p/'nodes.csv');tr=table(p/'trace.csv');ty=table(q/'trace.csv');assert np.array_equal(tr['time_s'],ty['time_s']) and np.array_equal(tr['window'],ty['window'])
 w=tr['window'];wave=np.sqrt(np.sum(w*(tr['velocity_m_s']-ty['velocity_m_s'])**2)/np.sum(w*ty['velocity_m_s']**2))
 report={'dspec_job':ds,'yspec_job':ys,'degree':r['degree'],'waveform_window_relative_l2':float(wave),'timings':{'DSpec':r,'YSpec':s},'resources':{'DSpec':{key:j[key] for key in ['peak_observed_threads','peak_observed_rss_gib_userspace']},'YSpec':{key:k[key] for key in ['peak_observed_threads','peak_observed_rss_gib_userspace']}},'sections':{}}
 arrays={}
 for section in ['along','across']:
  a=table(p/(section+'-angles.csv'));b=table(q/(section+'-angles.csv'));assert np.array_equal(a,b)
  n=len(nodes);m=len(a);x=np.fromfile(p/(section+'.bin'),np.float64).reshape(n,m);y=np.fromfile(q/(section+'.bin'),np.float64).reshape(n,m)
  assert np.isfinite(x).all() and np.isfinite(y).all()
  radial=(nodes['radius_m']<=6371e3-150e3)&(nodes['solid']==1)
  loc=np.ones(m,dtype=bool)
  coordinate=a['along_deg'] if section=='along' else a['across_deg']
  if section=='along':loc=(coordinate>5)&(coordinate<55)
  mask=radial[:,None]&loc[None,:]
  difference=x-y;weights=np.broadcast_to((nodes['radius_m']*nodes['dr_m'])[:,None],x.shape)
  relative=np.sqrt(np.sum(weights[mask]*difference[mask]**2)/np.sum(weights[mask]*y[mask]**2))
  limit=float(np.percentile(np.abs(np.concatenate((x[mask],y[mask]))),99))
  diff_limit=float(np.percentile(np.abs(difference[mask]),99))
  norm=float(np.linalg.norm(y[mask]));corr=float(np.dot(x[mask],y[mask])/(np.linalg.norm(x[mask])*norm))
  report['sections'][section]={'native_sides':n,'angles':m,'bulk_area_weighted_relative_l2':float(relative),'bulk_point_relative_l2':float(np.linalg.norm(difference[mask])/norm),'bulk_cosine_similarity':corr,'shared_limit':limit,'difference_own_limit':diff_limit,'max_abs':{'DSpec':float(np.max(np.abs(x))),'YSpec':float(np.max(np.abs(y))),'difference':float(np.max(np.abs(difference)))},'saturated_fraction':{label:float(np.mean(np.abs(v)>limit)) for label,v in [('DSpec',x),('YSpec',y),('difference',difference)]}}
  arrays[section]=(nodes,a,x,y,difference)
 return report,arrays,(tr,ty)

def plot_pair(root,ds,ys,out):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 from matplotlib.cm import ScalarMappable
 out.mkdir();report,arrays,traces=pair(root,ds,ys)
 for section,(nodes,angles,x,y,diff) in arrays.items():
  coordinate=angles['along_deg'] if section=='along' else angles['across_deg'];center=30 if section=='along' else 0
  theta=np.radians(coordinate-center);lim=report['sections'][section]['shared_limit']
  for own in [False,True]:
   fig,axes=plt.subplots(1,3,figsize=(16,6.3) if section=='along' else (6.5,8.5),layout=None)
   fig.subplots_adjust(top=.84,bottom=.18,wspace=.18)
   for ax,values,title in zip(axes,[x,y,diff],['DSpec','YSpec','DSpec − YSpec']):
    scale=report['sections'][section]['difference_own_limit'] if own and title.startswith('DSpec −') else lim
    norm=Normalize(-scale,scale)
    for element in np.unique(nodes['element']):
     ids=np.flatnonzero(nodes['element']==element)
     if len(ids)<2:continue
     radius=nodes['radius_m'][ids]/1000
     X=radius[:,None]*np.sin(theta)[None,:];Y=radius[:,None]*np.cos(theta)[None,:]
     ax.pcolormesh(X,Y,values[ids],cmap='RdBu',norm=norm,shading='gouraud',rasterized=True)
    dense=np.linspace(theta.min(),theta.max(),400)
    for radius in [3000,3480,5701,5961,6371]:ax.plot(radius*np.sin(dense),radius*np.cos(dense),color='.52',lw=.4)
    for edge in [theta.min(),theta.max()]:ax.plot(np.array([3000,6371])*np.sin(edge),np.array([3000,6371])*np.cos(edge),color='.45',lw=.5)
    if section=='along':
     for angle,radius,marker,color in [(theta.min(),6371-29.5,'*','#ffd400'),(theta.max(),6371-3,'v','white')]:ax.plot(radius*np.sin(angle),radius*np.cos(angle),marker=marker,ms=8,mfc=color,mec='black',mew=.65,clip_on=False)
    ax.set_aspect('equal');ax.set_xlim(6400*np.sin(theta.min())-50,6400*np.sin(theta.max())+50);ax.set_ylim(3000*np.cos(max(abs(theta.min()),abs(theta.max())))-60,6400);ax.axis('off');ax.set_title(title,fontsize=13,pad=16)
    stat=report['sections'][section];ax.text(.01,-.015,f"limit ±{scale:.3g}\nmax |sample| {np.max(np.abs(values)):.3g}\nclipped {np.mean(np.abs(values)>scale):.2%}",transform=ax.transAxes,fontsize=7,va='top')
    if own:
     fig.colorbar(ScalarMappable(norm=norm,cmap='RdBu'),ax=ax,fraction=.045,pad=.035,shrink=.72,label=r'$K_{\ln v_s}$ [s m$^{-3}$]')
   if not own:fig.colorbar(ScalarMappable(norm=Normalize(-lim,lim),cmap='RdBu'),ax=list(axes),fraction=.025,pad=.025,shrink=.78,label=r'$K_{\ln v_s}$ [s m$^{-3}$]')
   err=report['sections'][section]['bulk_area_weighted_relative_l2']
   fig.suptitle(f"{'Along the path' if section=='along' else 'Across the path at its midpoint'} · YSpec documentation shear packet · L={report['degree']}\nT velocity, 1000–1200 s, 2–10 mHz · bulk section difference {100*err:.3g}%",fontsize=11,y=.98)
   fig.text(.5,.025,('Separate labelled residual scale. ' if own else 'All three panels share one colour scale. ')+"Limits: combined 99th percentile below 150 km depth; along path excludes 5° near endpoints.\nColours interpolate within native SEM elements; clipping affects display only. Radius ≥ 3000 km. No degree-convergence claim.",ha='center',fontsize=7)
   name=section+('-residual-scale' if own else '-common-scale')
   for suffix in ['png','pdf']:fig.savefig(out/(name+'.'+suffix),dpi=200,bbox_inches='tight')
   plt.close(fig)
 tr,ty=traces;fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
 for axis,field,label in zip(axes,['velocity_m_s','adjoint_source','velocity_m_s'],['T velocity [m/s]','Traveltime adjoint source','DSpec − YSpec [m/s]']):
  if axis is axes[-1]:axis.plot(tr['time_s'],tr[field]-ty[field],lw=.8,color='black')
  else:axis.plot(tr['time_s'],tr[field],lw=.9,label='DSpec');axis.plot(ty['time_s'],ty[field],lw=.8,ls='--',label='YSpec');axis.legend(fontsize=8)
  axis.axvspan(1000,1200,color='.8',alpha=.3);axis.set_ylabel(label);axis.grid(alpha=.2)
 axes[-1].set_xlabel('Time after origin [s]');axes[-1].set_xlim(800,1400);fig.suptitle(f"Documentation case · waveform window relative L2 {report['waveform_window_relative_l2']:.3%}");fig.tight_layout()
 for suffix in ['png','pdf']:fig.savefig(out/('waveform-comparison.'+suffix),dpi=180,bbox_inches='tight')
 plt.close(fig);report['status']='COMPARISON';report['colour_scaling']='shared DSpec/YSpec combined99percentile; primary difference shares same scale; additional explicitly labelled residual-scale version';dump(out/'comparison.json',report)
 hashes={p.name:cp.sha(p) for p in out.iterdir() if p.is_file()};dump(out/'figure-manifest.json',{'files_sha256':hashes,'source_jobs':[ds,ys]})

def main():
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['parity','fd','pilot','plot']);p.add_argument('root',type=Path);p.add_argument('names',nargs='+');a=p.parse_args()
 if a.mode=='parity':parity(a.root,*a.names)
 elif a.mode=='fd':fd(a.root,*a.names)
 elif a.mode=='pilot':
  ds,ys,out=a.names;r,_,_=pair(a.root,ds,ys);ok=r['waveform_window_relative_l2']<.1 and all(x['bulk_area_weighted_relative_l2']<.2 for x in r['sections'].values());r['status']='PASS' if ok else 'FAIL';r['pilot_limits']={'waveform_window_relative_l2':.1,'bulk_area_weighted_relative_l2':.2};dump(out,r);assert ok,r
 else:ds,ys,out=a.names;plot_pair(a.root,ds,ys,Path(out))
if __name__=='__main__':main()
