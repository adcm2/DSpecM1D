#!/usr/bin/env python3
"""Render signed percentage residuals from immutable matched section kernels."""
from pathlib import Path
import argparse,json
import numpy as np
import check_doc_wedges as original
import check_pair as cp

def render(root,ds,ys,out):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 from matplotlib.cm import ScalarMappable
 report,arrays,_=original.pair(root,ds,ys);out.mkdir()
 report['percentage_formula']='100*(DSpec-YSpec)/max(abs(DSpec))'
 report['denominator_scope']='Each section independently; every original unclipped native-side sample, radius>=3000km, including shallow and endpoint samples.'
 report['percentage_sections']={}
 for section,(nodes,angles,x,y,diff) in arrays.items():
  denominator=float(np.max(np.abs(x)));assert denominator>0
  percent=100*diff/denominator;assert np.isfinite(percent).all()
  maximum=float(np.max(np.abs(percent)));assert np.isclose(maximum,100*np.max(np.abs(diff))/denominator,rtol=1e-14)
  shared=report['sections'][section]['shared_limit']
  bulk=100*report['sections'][section]['difference_own_limit']/denominator
  assert maximum>0 and bulk>0
  report['percentage_sections'][section]={'dspec_max_abs_denominator_s_m3':denominator,'full_range_percent_limit':maximum,'bulk_zoom_percent_limit':bulk,'bulk_zoom_clipped_fraction':float(np.mean(np.abs(percent)>bulk)),'full_range_clipped_fraction':float(np.mean(np.abs(percent)>maximum))}
  coordinate=angles['along_deg'] if section=='along' else angles['across_deg'];center=30 if section=='along' else 0;theta=np.radians(coordinate-center)
  for zoom in [False,True]:
   diff_limit=bulk if zoom else maximum
   fig=plt.figure(figsize=(17,6.3) if section=='along' else (9,9))
   grid=fig.add_gridspec(1,5,width_ratios=[1,1,.045,1,.045],left=.025,right=.96,bottom=.21,top=.80,wspace=.30)
   axes=[fig.add_subplot(grid[0,i]) for i in [0,1,3]]
   for ax,values,title,scale in zip(axes,[x,y,percent],['DSpec','YSpec','Difference [% of max |DSpec|]'],[shared,shared,diff_limit]):
    norm=Normalize(-scale,scale)
    for element in np.unique(nodes['element']):
     ids=np.flatnonzero(nodes['element']==element)
     if len(ids)<2:continue
     radius=nodes['radius_m'][ids]/1000
     ax.pcolormesh(radius[:,None]*np.sin(theta)[None,:],radius[:,None]*np.cos(theta)[None,:],values[ids],cmap='RdBu',norm=norm,shading='gouraud',rasterized=True)
    dense=np.linspace(theta.min(),theta.max(),400)
    for radius in [3000,3480,5701,5961,6371]:ax.plot(radius*np.sin(dense),radius*np.cos(dense),color='.52',lw=.4)
    for edge in [theta.min(),theta.max()]:ax.plot(np.array([3000,6371])*np.sin(edge),np.array([3000,6371])*np.cos(edge),color='.45',lw=.5)
    if section=='along':
     for angle,radius,marker,color in [(theta.min(),6371-29.5,'*','#ffd400'),(theta.max(),6371-3,'v','white')]:ax.plot(radius*np.sin(angle),radius*np.cos(angle),marker=marker,ms=8,mfc=color,mec='black',mew=.65,clip_on=False)
    ax.set_aspect('equal');ax.set_xlim(6400*np.sin(theta.min())-50,6400*np.sin(theta.max())+50);ax.set_ylim(3000*np.cos(max(abs(theta.min()),abs(theta.max())))-60,6400);ax.axis('off');ax.set_title(title,fontsize=11,pad=15)
    suffix='%' if ax is axes[-1] else ''
    ax.text(.01,-.025,f"limit ±{scale:.3g}{suffix}\nmax |sample| {np.max(np.abs(values)):.3g}{suffix}\nclipped {np.mean(np.abs(values)>scale):.2%}",transform=ax.transAxes,fontsize=7,va='top')
   fig.colorbar(ScalarMappable(norm=Normalize(-shared,shared),cmap='RdBu'),cax=fig.add_subplot(grid[0,2]),label=r'$K_{\ln v_s}$ [s m$^{-3}$]')
   fig.colorbar(ScalarMappable(norm=Normalize(-diff_limit,diff_limit),cmap='RdBu'),cax=fig.add_subplot(grid[0,4]),label=r'$100(K_D-K_Y)/\max|K_D|$ [%]')
   fig.suptitle(f"{'Along the path' if section=='along' else 'Across the path at its midpoint'} · YSpec documentation case · L={report['degree']}\nT velocity, 1000–1200 s, 2–10 mHz · {'bulk residual colour zoom' if zoom else 'full residual range'}",fontsize=11,y=.97)
   detail=('Difference colour limit: 99th percentile below 150 km depth'+('; >5° from endpoints.' if section=='along' else '.')+' Larger values are clipped for display.' if zoom else 'Difference colour limit covers every residual sample; no residual clipping.')
   fig.text(.5,.065,f"Signed difference = 100 × (DSpec − YSpec) / max |DSpec|. Denominator = {denominator:.6g} s m⁻³ over the full, unclipped section.\n"+detail+" DSpec and YSpec retain their shared colour scale.\nColours interpolate within native SEM elements. Radius ≥ 3000 km. No degree-convergence claim.",ha='center',fontsize=7)
   name=section+('-percent-bulk-zoom' if zoom else '-percent-full-range')
   for suffix in ['png','pdf']:fig.savefig(out/(name+'.'+suffix),dpi=200,bbox_inches='tight')
   plt.close(fig)
 original.dump(out/'percentage-comparison.json',report)
 original.dump(out/'figure-manifest.json',{'files_sha256':{p.name:cp.sha(p) for p in out.iterdir() if p.is_file()},'source_jobs':[ds,ys],'source_root':str(root)})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('ds');p.add_argument('ys');p.add_argument('out',type=Path);a=p.parse_args();render(a.root,a.ds,a.ys,a.out)
