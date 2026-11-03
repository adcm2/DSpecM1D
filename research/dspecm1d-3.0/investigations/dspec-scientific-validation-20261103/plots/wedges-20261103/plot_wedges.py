#!/usr/bin/env python3
"""Render saved great-circle sections; no solver calls or evidence edits."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, Normalize
from matplotlib.ticker import FuncFormatter
from matplotlib.cm import ScalarMappable

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
RAW=ROOT/"study2-v1/evidence/producer-raw-v1"
REF=ROOT/"study2-v1/reference/study1-yspec-tight-l200"
inputs={}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def record(p):
    inputs[str(p.relative_to(ROOT))]={"sha256":digest(p),"bytes":p.stat().st_size}
def rows(p):
    record(p)
    with p.open() as f:return list(csv.DictReader(f))
def load(d):
    nodes=rows(d/"nodes.csv");angles=rows(d/"along-angles.csv")
    record(d/"along.bin")
    data=np.fromfile(d/"along.bin",dtype="<f8").reshape(len(nodes),len(angles))
    assert np.isfinite(data).all()
    return nodes,angles,data
# Bind producer outputs to the immutable manifest, not just their current path.
manifest=json.loads((RAW/"manifest.json").read_text())
assert digest(RAW/"state.json")=="9cd09f3e4bd0979a39a240cfa8d2a1a2bf631750ab50a8869a55db49cdfad704"
for name,entry in manifest.items():
    assert (RAW/name).stat().st_size==entry["bytes"] and digest(RAW/name)==entry["sha256"]
cases={}
for L,dname,yroot in [(200,"l200-dspec-sparse-bridge",REF),(300,"l300-dspec",RAW/"l300-yspec-tight/data")]:
    nodes,angles,D=load(RAW/dname/"data");yn,ya,Y=load(yroot)
    assert nodes==yn
    target=np.array([float(x["along_deg"]) for x in angles])
    source=np.array([float(x["along_deg"]) for x in ya])
    indices=[int(np.argmin(abs(source-x))) for x in target]
    assert np.max(abs(source[indices]-target))<1e-12
    Y=Y[:,indices]  # Exact saved coordinates, no interpolation between backend grids.
    cases[L]={"nodes":nodes,"angles":target,"D":D,"Y":Y,"diff":D-Y}
    denominator=float(np.max(abs(D)))
    assert denominator>0
    cases[L]["percent_denominator"]=denominator
    cases[L]["percent_diff"]=100*(D-Y)/denominator
assert cases[200]["nodes"]==cases[300]["nodes"]
assert np.array_equal(cases[200]["angles"],cases[300]["angles"])
scale=1e-15
vmax=max(np.max(abs(c[k]))/scale for c in cases.values() for k in ("D","Y"))
dmax=max(np.max(abs(c["percent_diff"])) for c in cases.values())
norm=SymLogNorm(linthresh=.01,linscale=.7,vmin=-vmax,vmax=vmax,base=10)
dlimit=float(np.ceil(dmax/10)*10)
dnorm=Normalize(vmin=-dlimit,vmax=dlimit)
percent_formatter=FuncFormatter(lambda x,pos: f"{x:g}%")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"svg.fonttype":"none"})
arc=np.linspace(-np.pi/6,np.pi/6,500)
def wedge(ax,c,values,color_norm,value_scale=scale):
    nodes=c["nodes"];theta=np.radians(c["angles"]-30)
    elems={}
    for i,n in enumerate(nodes):elems.setdefault(n["element"],[]).append(i)
    for indices in elems.values():
        r=np.array([float(nodes[i]["radius_m"])/1000 for i in indices])
        if len(r)<2:continue
        assert np.all(np.diff(r)>0)
        x=r[:,None]*np.sin(theta)[None,:];z=r[:,None]*np.cos(theta)[None,:]
        ax.pcolormesh(x,z,values[indices]/value_scale,shading="gouraud",cmap="RdBu_r",
                      norm=color_norm,rasterized=True)
    for r in (6371,6346.6,6336,6151,5971,5701,3480,3022.953):
        ax.plot(r*np.sin(arc),r*np.cos(arc),color=".35",lw=.45,alpha=.65)
    for t in (-np.pi/6,np.pi/6):
        ax.plot([3022.953*np.sin(t),6371*np.sin(t)],[3022.953*np.cos(t),6371*np.cos(t)],color=".3",lw=.6)
    ax.scatter([6341.5*np.sin(-np.pi/6)],[6341.5*np.cos(-np.pi/6)],marker="*",s=95,
               facecolor="#ffd046",edgecolor="black",linewidth=.7,zorder=5,clip_on=False)
    ax.scatter([6368*np.sin(np.pi/6)],[6368*np.cos(np.pi/6)],marker="v",s=55,
               facecolor="#44b8c2",edgecolor="black",linewidth=.7,zorder=5,clip_on=False)
    ax.text(-3185,5635,"Source",ha="right",fontsize=9)
    ax.text(3185,5635,"Receiver",ha="left",fontsize=9)
    ax.text(0,6420,"Surface",ha="center",fontsize=9)
    ax.text(0,3550,"CMB",ha="center",fontsize=8,color=".3")
    ax.text(0,2870,"r = 3023 km",ha="center",fontsize=8,color=".3")
    ax.set_xlim(-3650,3650);ax.set_ylim(2500,6580);ax.set_aspect("equal");ax.axis("off")
stats={}
for L,c in cases.items():
    fig,axes=plt.subplots(1,3,figsize=(16,4.65),layout="constrained")
    for ax,key,title in zip(axes,("D","Y","diff"),("DSpec","YSpec","DSpec - YSpec (%)")):
        color_norm=dnorm if key=="diff" else norm
        values=c["percent_diff"] if key=="diff" else c[key]
        wedge(ax,c,values,color_norm,1 if key=="diff" else scale);ax.set_title(title,fontsize=11)
        cb=fig.colorbar(ScalarMappable(norm=color_norm,cmap="RdBu_r"),ax=ax,
                       orientation="horizontal",fraction=.05,pad=.015,shrink=.88)
        if key=="diff":
            cb.set_ticks(np.linspace(-dlimit,dlimit,5))
            cb.ax.xaxis.set_major_formatter(percent_formatter)
            cb.set_label("Percentage difference (linear)",fontsize=12)
        else:
            cb.set_label(r"$K_{\ln V_s}$  ($10^{-15}$ s m$^{-3}$)",fontsize=9)
        cb.ax.tick_params(labelsize=12 if key=="diff" else 8)
        stats[f"L{L}-{key}"]={"minimum":float(np.min(c[key])),"maximum":float(np.max(c[key])),
                             "nonfinite":0,"display_clipped_values":0}
        if key=="diff":
            stats[f"L{L}-{key}"].update(percent_denominator=c["percent_denominator"],
                minimum_percent=float(np.min(values)),maximum_percent=float(np.max(values)))
    fig.suptitle(f"L{L}: along-path great-circle kernel wedge (0-60 degrees)",fontsize=14)
    fig.supxlabel("Difference normalized by the maximum |DSpec| over each section; shared LINEAR percentage scale.\n"
        "All saved values retained. Native source-shell traces are not matched: the percentage difference includes that known mismatch.",
        fontsize=9)
    fig.savefig(OUT/f"kernel-wedge-L{L}-percent-linear.png",dpi=160)
    fig.savefig(OUT/f"kernel-wedge-L{L}-percent-linear.svg")
    plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(16,4.65),layout="constrained")
c=cases[200]
assert np.max(abs(cases[300]["D"]-c["D"]))/scale <= vmax
cutoff_percent=100*(cases[300]["D"]-c["D"])/c["percent_denominator"]
cutoff_max=float(np.max(abs(cutoff_percent)))
cutoff_limit=float(np.ceil(cutoff_max/50)*50)
cutoff_norm=Normalize(vmin=-cutoff_limit,vmax=cutoff_limit)
for i,(ax,values,title) in enumerate(zip(axes,(c["D"],cases[300]["D"],cutoff_percent),
                           ("DSpec L200","DSpec L300","DSpec L300 - L200 (%)"))):
    color_norm=cutoff_norm if i==2 else norm
    wedge(ax,c,values,color_norm,1 if i==2 else scale);ax.set_title(title,fontsize=11)
    cb=fig.colorbar(ScalarMappable(norm=color_norm,cmap="RdBu_r"),ax=ax,orientation="horizontal",
                   fraction=.05,pad=.015,shrink=.88)
    if i==2:
        cb.set_ticks(np.linspace(-cutoff_limit,cutoff_limit,5))
        cb.ax.xaxis.set_major_formatter(percent_formatter)
        cb.set_label("Percentage difference (linear)",fontsize=12)
    else:
        cb.set_label(r"$K_{\ln V_s}$  ($10^{-15}$ s m$^{-3}$)",fontsize=9)
    cb.ax.tick_params(labelsize=12 if i==2 else 8)
fig.suptitle("Changing the physical angular cutoff changes the kernel",fontsize=14)
fig.supxlabel("Difference normalized by the maximum |DSpec L200| over the section; LINEAR percentage scale, no clipping.\n"
    "Saved sparse section: 68 angles. Interpolation within each SEM element is for display; no cross-interface averaging.",
    fontsize=9)
fig.savefig(OUT/"kernel-wedge-cutoff-change-percent-linear.png",dpi=160)
fig.savefig(OUT/"kernel-wedge-cutoff-change-percent-linear.svg")
plt.close(fig)
# Full-width difference panels keep the percentage axis readable in chat.
for filename,c,values,color_norm,title,note in [
    (f"kernel-difference-L{L}-linear",c,c["percent_diff"],dnorm,
     f"L{L}: DSpec - YSpec",
     "100 (DSpec - YSpec) / max |DSpec| over this section\nRaw source-shell trace mismatch retained.")
    for L,c in cases.items()
]+[("kernel-difference-cutoff-linear",cases[200],cutoff_percent,cutoff_norm,
    "DSpec: L300 - L200",
    "100 (L300 - L200) / max |L200| over this section\nSame saved coordinates at both cutoffs.")]:
    fig,ax=plt.subplots(figsize=(7.4,5.7),layout="constrained")
    wedge(ax,c,values,color_norm,value_scale=1)
    fig.suptitle(title,fontsize=17)
    cb=fig.colorbar(ScalarMappable(norm=color_norm,cmap="RdBu_r"),ax=ax,
        orientation="horizontal",fraction=.09,pad=.025,shrink=.94)
    cb.set_ticks(np.linspace(color_norm.vmin,color_norm.vmax,5))
    cb.ax.xaxis.set_major_formatter(percent_formatter)
    cb.ax.tick_params(labelsize=15,pad=7)
    cb.set_label("Percentage difference — linear scale",fontsize=15,labelpad=9)
    fig.supxlabel(note,fontsize=11)
    fig.savefig(OUT/f"{filename}.png",dpi=160)
    fig.savefig(OUT/f"{filename}.svg")
    plt.close(fig)
    zoom_limit=10.0 if "cutoff" in filename else .01
    zoom_norm=Normalize(vmin=-zoom_limit,vmax=zoom_limit)
    fig,ax=plt.subplots(figsize=(7.4,5.7),layout="constrained")
    wedge(ax,c,values,zoom_norm,value_scale=1)
    fig.suptitle(title,fontsize=17)
    cb=fig.colorbar(ScalarMappable(norm=zoom_norm,cmap="RdBu_r"),ax=ax,
        orientation="horizontal",fraction=.09,pad=.025,shrink=.94,extend="both")
    cb.set_ticks(np.linspace(-zoom_limit,zoom_limit,5))
    cb.ax.xaxis.set_major_formatter(percent_formatter)
    cb.ax.tick_params(labelsize=15,pad=7)
    cb.set_label("Percentage difference — linear scale",fontsize=15,labelpad=9)
    peak=float(np.max(abs(values)))
    fig.supxlabel(f"Normalized by the section maximum of the reference kernel.\n"
        f"Colours saturate beyond ±{zoom_limit:g}%; full-data peak = {peak:.3g}%.",fontsize=11)
    fig.savefig(OUT/f"{filename}-zoom.png",dpi=160)
    fig.savefig(OUT/f"{filename}-zoom.svg")
    plt.close(fig)
# Decompose the exact 89-probe metric; this is not a spatial volume integral.
def probes(job):
    return {p["probe_id"]:p for p in rows(RAW/job/"data/fixed-probe-kernel.csv")}
a,b=probes("l200-dspec-sparse-bridge"),probes("l300-dspec")
depths={}
for pid,p in a.items():
    shell=p["layer"]=="9" and p["radial_trace"] in ("inner","outer")
    if (p["material_side"]=="solid" and p["location_class"]=="bulk" and
        min(float(p["source_distance_3d_m"]),float(p["receiver_distance_3d_m"]))>=100000 and not shell):
        depth=round((6371000-float(p["radius_m"]))/1000,6)
        u,v=float(p["kernel_vs_s_m-3"]),float(b[pid]["kernel_vs_s_m-3"])
        depths.setdefault(depth,[]).append((u,v))
assert sum(map(len,depths.values()))==89
err2=sum((v-u)**2 for pairs in depths.values() for u,v in pairs)
breakdown={"metric":"same 89 far-solid-interior fixed probes; unweighted discrete vectors, not volume integrals",
    "relative_change_percent":100*np.sqrt(err2/sum(u*u for pairs in depths.values() for u,v in pairs)),
    "depth_km":[{"depth":depth,"count":len(pairs),
        "relative_change_percent":100*np.sqrt(sum((v-u)**2 for u,v in pairs)/sum(u*u for u,v in pairs)),
        "fraction_of_total_squared_change_percent":100*sum((v-u)**2 for u,v in pairs)/err2}
        for depth,pairs in sorted(depths.items())],
    "24_to30km_fraction_of_total_squared_change_percent":
        100*sum((v-u)**2 for depth,pairs in depths.items() if 24<=depth<=30 for u,v in pairs)/err2}
(OUT/"cutoff-depth-breakdown.json").write_text(json.dumps(breakdown,indent=2)+"\n")
output={"inputs":inputs,"plot_semantics":{
    "data":"saved along-path great-circle native sections; physical solve degrees 200/300",
    "angular_sampling":"68 saved angles at both cutoffs; not enough to certify full angular bandwidth",
    "radial_extent_m":[3022952.9061119016,6371000],"full_radial_operator_retained":True,
    "interpolation":"Gouraud display interpolation within SEM elements; interfaces and duplicate traces never averaged",
    "L200_yspec_mapping":"exact matching subset of saved 135 angles to DSpec's 68; tolerance 1e-12 degrees",
    "color":"symmetric-log, no clipping, shared physical kernel range across L200/L300",
    "raw_difference":"includes unmatched source-shell radial trace convention; not a matched-trace error everywhere",
    "percentage_difference":"100*(D-Y)/max_section(abs(D)), separately at each cutoff; no pointwise division or denominator floor",
    "percentage_color_scale":{"type":"linear","backend_limits_percent":[-dlimit,dlimit],
                              "cutoff_limits_percent":[-cutoff_limit,cutoff_limit]},
    "standalone_zoom_scales":{"type":"linear","backend_limits_percent":[-.01,.01],
        "cutoff_limits_percent":[-10,10],"saturation":"both colorbar extensions and full-data peak labelled; values not removed"},
    "cutoff_percentage_difference":"100*(D300-D200)/max_section(abs(D200))",
    "source_contact_terms_included":False},
    "stats":stats,"outputs":{p.name:{"sha256":digest(p),"bytes":p.stat().st_size}
        for p in sorted(OUT.iterdir()) if p.suffix in (".png",".svg")}}
(OUT/"provenance-percent-linear.json").write_text(json.dumps(output,indent=2)+"\n")
print(json.dumps({"color_limit_kernel_1e15":vmax,"color_limit_diff_percent":dmax,
                  "cutoff_percent_max_abs":cutoff_max,"stats":stats},indent=2))
