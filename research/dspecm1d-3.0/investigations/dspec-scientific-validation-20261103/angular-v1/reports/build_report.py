#!/usr/bin/env python3
"""Reproduce the report/figures from actual sealed angular-study evidence."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, SymLogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FuncFormatter

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
RAW=ROOT/"evidence/producer-raw-v1"
TERM=ROOT/"evidence/checker-terminal-v2"
REPORT=TERM/"dspec-science-angular-v1-check-v1/data/angular-report.json"
METRICS=ROOT/"reviews/angular-independent-raw-metrics-20261103.json"
SUPPLEMENT=ROOT/"reviews/angular-independent-supplementary-metrics-20261103.json"
CLOSURE=ROOT/"resource-coordination/final-closure-observation.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(REPORT)=="c9a162ad9dd4d05b86a0e7b8515e427f10ba18ee7d23510e9c6af695c789892c"
assert sha(METRICS)=="0e03770ade5be6a3da38e17e4a6434d4896506def848909ae1c41814eeee36f5"
assert sha(SUPPLEMENT)=="cd5764633dddbb08967e4864a3555a5ef3c0f324f8ca759dc0ddf7364d19dfc8"
assert sha(CLOSURE)=="809d4924f536086733c5cafe8259215b162c24d3f135096af5dcacaae2134d79"
closure=json.loads(CLOSURE.read_text())
assert not closure["live_dspecm1d_services"] and closure["slice"]["TasksCurrent"]=="0"
for folder in (RAW,TERM):
    manifest=json.loads((folder/"manifest.json").read_text())
    for name,entry in manifest.items():
        p=folder/name
        assert p.stat().st_size==entry["bytes"] and sha(p)==entry["sha256"]
m=json.loads(METRICS.read_text());s=json.loads(SUPPLEMENT.read_text());report=json.loads(REPORT.read_text())
assert report["state_sha256"]==m["producer_state_sha256"]=="35ac18a3b42d13c3ae2f524004bda3a8fcbffe32414bfb7e34ed19b127ef603f"
assert not any(v["certified"] for levels in report["acceptance"].values() for v in levels.values())
assert not any(report["radial_control_pass"].values())
def pct(value):return 100*value
def f(value):return f"{value:.6f}"
def worst_depth(d):
    return max((max(v["relative_l2"],v["max_relative"]),k) for k,v in d["depth_strata"].items())
def metric_row(label,d):
    return f'{label} & {f(pct(d["window_weighted"]["relative_l2"]))} & {f(pct(d["normalization_relative"]))} & {f(pct(d["P"]["relative_l2"]))} & {f(pct(d["kernel_far"]["relative_l2"]))} '+r"\\"
degree_rows="\n".join(metric_row(f'{dict(dspec="DSpec", yspec="YSpec")[code]} {a} to {b}',m["degree"][code][f"{a}-to-{b}"])
    for code in ("dspec","yspec") for a,b in ((300,400),(400,600),(600,800)))
radial_rows="\n".join(f'NQ{a} to NQ{b} & {f(pct(m["radial"][key]["window_weighted"]["relative_l2"]))} & '
    f'{f(pct(m["radial"][key]["kernel_far"]["relative_l2"]))} & {f(pct(worst_depth(m["radial"][key])[0]))} '
    +r"\\"
    for a,b,key in ((6,8,"nq6-to-nq8"),(6,10,"nq6-to-nq10"),(8,10,"nq8-to-nq10")))
cross_rows="\n".join(f'L{L}, NQ6 & {f(pct(m["cross_backend_nongating"][str(L)]["window_weighted"]["relative_l2"]))} & '
    f'{f(pct(m["cross_backend_nongating"][str(L)]["kernel_far"]["relative_l2"]))} '+r"\\"
    for L in (400,600,800))
order_rows="\n".join(f'NQ{n} & {f(pct(s["matched_800_D_order_vs_Y_diagnostic"][str(n)]["relative_l2"]))} '+r"\\"
    for n in (6,8,10))
depth_rows="\n".join(f'{depth} & {v["count"]} & {f(pct(v["relative_l2"]))} & {f(pct(v["max_relative"]))} '+r"\\"
    for depth,v in m["degree"]["dspec"]["600-to-800"]["depth_strata"].items())
plt.rcParams.update({"font.family":"serif","font.size":10,"axes.spines.top":False,
                     "axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
fig,axes=plt.subplots(1,2,figsize=(7.1,3.15),constrained_layout=True)
for ax,key,title in zip(axes,("window_weighted","kernel_far"),("Tapered measurement-window velocity","89 far-solid-interior kernel probes")):
    for code,color in (("dspec","#244f75"),("yspec","#a12432")):
        values=[pct(m["degree"][code][f"{a}-to-{b}"][key]["relative_l2"]) for a,b in ((300,400),(400,600),(600,800))]
        ax.plot(range(3),values,"o-",color=color,label="DSpec NQ6" if code=="dspec" else "YSpec",ms=4,lw=1)
    ax.axhline(.1,color=".35",ls="--",lw=.8,label="0.1% target")
    ax.set_yscale("log");ax.set_xticks(range(3),["300 to 400","400 to 600","600 to 800"],rotation=18)
    ax.set_xlabel("Physical degree cutoff");ax.set_ylabel("Higher-reference vector change (%)")
    ax.set_title(title,fontsize=9);ax.legend(frameon=False,fontsize=7)
fig.savefig(HERE/"figures/degree-convergence.pdf",bbox_inches="tight")
fig.savefig(HERE/"figures/degree-convergence.png",dpi=170,bbox_inches="tight")
plt.close(fig)

# Native great-circle sections: paired at the exact SAME nodes/angles, NQ6.
def rows(p):
    with p.open() as file:return list(csv.DictReader(file))
cases={}
for L in (400,600,800):
    d=RAW/f"angular-l{L}-dspec-nq6/data";y=RAW/f"angular-l{L}-yspec-nq6/data"
    nodes,angles=rows(d/"nodes.csv"),rows(d/"along-angles.csv")
    assert nodes==rows(y/"nodes.csv") and angles==rows(y/"along-angles.csv")
    D=np.fromfile(d/"along.bin",dtype="<f8").reshape(len(nodes),len(angles))
    Y=np.fromfile(y/"along.bin",dtype="<f8").reshape(len(nodes),len(angles))
    assert np.isfinite(D).all() and np.isfinite(Y).all()
    cases[L]={"nodes":nodes,"angles":np.radians([float(v["along_deg"])-30 for v in angles]),
              "D":D,"Y":Y,"percent":100*(D-Y)/np.max(abs(D))}
vmax=max(np.max(abs(c[k]))/1e-15 for c in cases.values() for k in ("D","Y"))
knorm=SymLogNorm(linthresh=.01,linscale=.7,vmin=-vmax,vmax=vmax)
pnorm=Normalize(vmin=-.1,vmax=.1)
arc=np.linspace(-np.pi/6,np.pi/6,500)
def wedge(ax,c,values,norm,unit_scale):
    groups={}
    for i,node in enumerate(c["nodes"]):groups.setdefault(node["element"],[]).append(i)
    for indices in groups.values():
        radii=np.array([float(c["nodes"][i]["radius_m"])/1000 for i in indices])
        if len(radii)<2:continue
        ax.pcolormesh(radii[:,None]*np.sin(c["angles"]),radii[:,None]*np.cos(c["angles"]),
            values[indices]/unit_scale,norm=norm,cmap="RdBu_r",shading="gouraud",rasterized=True)
    for radius in (6371,6346.6,6336,6151,5971,5701,3480,3022.953):
        ax.plot(radius*np.sin(arc),radius*np.cos(arc),color=".4",lw=.45,alpha=.7)
    for theta in (-np.pi/6,np.pi/6):
        ax.plot([3022.953*np.sin(theta),6371*np.sin(theta)],
                [3022.953*np.cos(theta),6371*np.cos(theta)],color=".35",lw=.5)
    ax.scatter([-3170.75],[6341.5*np.cos(np.pi/6)],marker="*",s=95,c="#ffd046",edgecolor="black",zorder=5)
    ax.scatter([3184],[6368*np.cos(np.pi/6)],marker="v",s=55,c="#44b8c2",edgecolor="black",zorder=5)
    ax.text(-3185,5635,"Source",ha="right",fontsize=9)
    ax.text(3185,5635,"Receiver",ha="left",fontsize=9)
    ax.text(0,6420,"Surface",ha="center",fontsize=9)
    ax.text(0,3550,"CMB",ha="center",fontsize=9,color=".35")
    ax.set_xlim(-3700,3700);ax.set_ylim(2500,6580);ax.set_aspect("equal");ax.axis("off")
percent_formatter=FuncFormatter(lambda x,pos:f"{x:g}%")
figure_stats={}
for L,c in cases.items():
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),layout="constrained")
    for ax,key,title in zip(axes,("D","Y","percent"),("DSpec NQ6","YSpec","DSpec - YSpec (%)")):
        norm=pnorm if key=="percent" else knorm
        wedge(ax,c,c[key],norm,1 if key=="percent" else 1e-15);ax.set_title(title,fontsize=13)
        cb=fig.colorbar(ScalarMappable(norm=norm,cmap="RdBu_r"),ax=ax,orientation="horizontal",
            fraction=.06,pad=.02,shrink=.9,extend="both" if key=="percent" else "neither")
        if key=="percent":
            cb.set_ticks([-.1,-.05,0,.05,.1]);cb.ax.xaxis.set_major_formatter(percent_formatter)
            cb.set_label("Percentage difference (linear)",fontsize=12);cb.ax.tick_params(labelsize=12)
        else:
            cb.set_ticks([-100,-1,0,1,100])
            cb.set_label(r"$K_{\ln V_s}$ ($10^{-15}$ s m$^{-3}$; symmetric log)",fontsize=10)
    peak=float(np.max(abs(c["percent"])))
    fig.suptitle(f"L{L}: native along-path great-circle kernel wedges",fontsize=15)
    fig.supxlabel(f"Percent = 100 (DSpec - YSpec) / max section |DSpec|; display saturates beyond ±0.1%, raw peak {peak:.3g}%.\n"
        "Raw difference includes radial discretization differences and unmatched source-shell traces. Saved 68-angle grid is for display.",fontsize=10)
    fig.savefig(HERE/f"figures/kernel-comparison-L{L}.png",dpi=160)
    fig.savefig(HERE/f"figures/kernel-comparison-L{L}.svg")
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.4,5.7),layout="constrained")
    wedge(ax,c,c["percent"],pnorm,1);fig.suptitle(f"L{L}: DSpec NQ6 - YSpec",fontsize=17)
    cb=fig.colorbar(ScalarMappable(norm=pnorm,cmap="RdBu_r"),ax=ax,orientation="horizontal",
        fraction=.09,pad=.02,shrink=.94,extend="both")
    cb.set_ticks([-.1,-.05,0,.05,.1]);cb.ax.xaxis.set_major_formatter(percent_formatter)
    cb.ax.tick_params(labelsize=15,pad=6);cb.set_label("Percentage difference — linear scale",fontsize=15,labelpad=9)
    fig.supxlabel(f"Normalized by max section |DSpec|; colours saturate beyond ±0.1%.\n"
                  f"Raw peak {peak:.3g}%; source-shell traces are not fully paired.",fontsize=11)
    fig.savefig(HERE/f"figures/kernel-difference-L{L}-linear.png",dpi=160)
    plt.close(fig)
    figure_stats[str(L)]={"raw_difference_peak_percent_of_D_section_max":peak,
        "linear_display_limit_percent":.1,"saturated_sample_percent":100*float(np.mean(abs(c["percent"])>.1))}

peak_rss=max(v["peak_rss_gib_userspace"] for v in m["resources"].values())
summary={"checker_report_sha256":sha(REPORT),"independent_metrics_sha256":sha(METRICS),
    "supplement_sha256":sha(SUPPLEMENT),"resource_closure_sha256":sha(CLOSURE),"certified_cutoff":None,
    "conditional_extension_triggered":False,"peak_numeric_rss_gib_userspace":peak_rss,
    "figure_stats":figure_stats,"scope":"2/4/8/10mHz;1000..1200s window;40s taper;one fixture only"}
(HERE/"report-data.json").write_text(json.dumps(summary,indent=2)+"\n")
text=r'''\documentclass[11pt,a4paper]{article}
\usepackage{alex-report}
\ReportSetup{DSpecM1D scientific validation}{Angular convergence}{3 November 2026}
\hypersetup{pdftitle={DSpecM1D angular convergence at L400, L600 and L800}}
\begin{document}
\ReportTitle{Angular convergence at L400, L600 and L800}{No tested cutoff establishes the requested 0.1\% target}
\textbf{Result.} The eight isolated calculations and their independent checks completed successfully, but the requested numerical target was not demonstrated. L600 to L800 changes the tapered measurement-window velocity by 0.251811\% and the selected far-solid-interior kernel vector by 137.884259\%. YSpec shows essentially the same degree dependence. The L800 radial checks also fail the complete predeclared criterion. L800 is therefore not a converged reference merely because it is the highest tested level.

\section{The cutoff is specific to the filter and window}
The fixture is attenuated anisotropic PREM with full gravity and the full radial operator, source depth 29.5 km, solid-side receiver depth 3 km, separation 60 degrees and azimuth 30 degrees. The observable is transverse-velocity traveltime with the four-corner 2/4/8/10 mHz filter, 1000--1200 s measurement window and 40 s taper. The stored FFT grid has 35 bins (7--41), of which 32 have nonzero unit-adjoint support; bins 7, 8 and 41 are skipped placeholders, not evaluated zero unit fields. All forward bins are retained.

These results certify neither a different frequency band nor a different time window. The user's illustrative [10,50] and [10,100] bands were not substituted for this fixture. Each cutoff uses its own waveform-derived adjoint weights $P_L$.

\section{Observed degree changes}
\begin{center}\small\begin{tabular}{lrrrr}\toprule
Change & Window $u$ (\%) & $|\Delta N|$ (\%) & $P$ (\%) & Far kernel (\%)\\\midrule
DEGREE
\bottomrule\end{tabular}\end{center}
The higher cutoff is the reference: vector changes are $100\|a-b\|_2/\|b\|_2$. The velocity gate uses tapered samples $w u$ in the measurement window. $N$ is the measurement-rate normalization and $P$ the complex adjoint-weight vector over the stored band. The far kernel is an unweighted vector of 89 fixed probes in solid bulk material at least 100 km from both source and receiver, excluding the source shell. These are discrete sample diagnostics, not volume or $H^1$ norms.

L400 versus L800 changes DSpec's window velocity by 0.128427\% and far kernel by 78.362329\%. Neither the successive increments nor this direct comparison passes. The nonmonotone sequence does not support extrapolating a sufficient cutoff.

\newpage
\begin{figure}[H]\centering\includegraphics[width=0.95\linewidth]{figures/degree-convergence.pdf}
\caption{Observed higher-reference changes. The DSpec NQ6 radial control fails at L800, so its changes are not a pure angular-error estimate. YSpec independently exhibits the same large cutoff dependence. The dashed line is the 0.1\% operational target.}
\end{figure}
\begin{figure}[H]\centering\includegraphics[width=\linewidth]{figures/kernel-comparison-L400.png}
\caption{L400: DSpec, YSpec and linear percentage difference. Physical kernel panels share the symmetric-log scale; percentage colours saturate at $\pm0.1\%$, with raw peak 8.30\%.}
\end{figure}
\begin{figure}[H]\centering\includegraphics[width=\linewidth]{figures/kernel-comparison-L800.png}
\caption{L800: the same scales and conventions. Raw percentage peak is 3.65\%. Displayed differences include radial error and unmatched source-shell traces; these section extrema differ from the required 89-probe norm. The coarse saved angles do not resolve high-L spatial bandwidth.}
\end{figure}

\newpage
\section{Radial accuracy must be revisited at higher degree}
All angular-ladder DSpec cases use NQ6 and a fixed normalized mesh parameter. Its value is $h=0.02035382429229197$; it is not a length in metres. L800 was repeated at NQ8 and NQ10 with exactly the same element boundaries.
\begin{center}\small\begin{tabular}{lrrr}\toprule
L800 change & Window $u$ (\%) & Far kernel (\%) & Worst depth metric (\%)\\\midrule
RADIAL
\bottomrule\end{tabular}\end{center}
The last column is the maximum over each required depth cohort's relative vector norm and maximum absolute difference divided by that cohort's reference maximum. All three comparisons fail the complete criterion. NQ8 versus NQ10 passes the aggregate far-vector test (0.038370\%) but fails the 50 km depth stratum: vector change 0.187392\%, maximum-scaled change 0.185817\%. An aggregate norm alone would hide that failure.

Thus the radial stability established at L200 did not transfer to NQ6 at L800. This does not establish that NQ10 itself is inaccurate; agreement with a lower order is not an exact error estimate.

\section{Cross-code agreement and its limits}
\begin{center}\small\begin{tabular}{lrr}\toprule
Matched calculation & Window $u$ (\%) & Far kernel (\%)\\\midrule
CROSS
\bottomrule\end{tabular}\end{center}
These differences use YSpec as reference and each code's own measurement weights. The earlier close NQ6 agreement at L100--L300 must not be extrapolated to higher cutoffs. The L800 kernel comparison improves strongly with radial polynomial order:
\begin{center}\small\begin{tabular}{lr}\toprule
DSpec order versus YSpec at L800 & Far kernel difference (\%)\\\midrule
ORDER
\bottomrule\end{tabular}\end{center}
NQ10 agrees with YSpec very closely on this vector. That is encouraging evidence about radial accuracy at fixed L800; it does not establish angular convergence, or replace the predeclared depth-by-depth self-convergence checks. YSpec's SEM-like output mesh is sampling only, not its ODE-solver discretization; its ODE tolerance remains $10^{-7}$.

\newpage
\section{Where the kernel still changes}
The L600-to-L800 DSpec NQ6 kernel change is not confined to exact source-shell samples:
\begin{center}\small\begin{tabular}{rrrr}\toprule
Depth (km) & Probe count & Vector change (\%) & Max-scaled change (\%)\\\midrule
DEPTH
\bottomrule\end{tabular}\end{center}
The shallow groups are especially unstable, but even the 300 and 500 km groups exceed 0.1\%. Exact source-shell traces and material-interface/near-contact cohorts are retained separately as diagnostics; they are not required bounded-pointwise convergence tests.

Holding the L800 adjoint weights fixed changes the far-vector L600-to-L800 result only from 137.884259\% to 137.835946\%. The large observed difference therefore cannot be explained simply by the small change in waveform-derived weighting. The counterfactual is explicitly $P_{800}C_L$; it is not a native lower-cutoff kernel or an independent YSpec unit-field calculation. Products include cross terms, so their cutoff differences are not additive isolated-degree kernel contributions.

All saved native sections are finite. This establishes neither a bounded point-source/receiver limit nor divergence. No new $H^1$ formulation, regularization, smoothing, receiver averaging or contact correction was introduced. The current API does not support an independent high-L three-dimensional action validation within the bounded study; that gap was documented before results. Study 1's low-degree finite-difference action test retains only its original scope.

\section{Acceptance and stopping point}
Before results, acceptance required the chosen cutoff versus L800 and two successive increments to pass 0.1\% for the tapered waveform, $N$, $P$, and every required kernel-cohort vector/max-scaled measure. Zero reference cases are explicit, with no arbitrary denominator floors. DSpec additionally required the L800 radial controls to pass. YSpec's angular result is reported separately.

No candidate passed. L600 also has only one higher increment, and L800 has no higher reference. The conditional two-job NQ10 extension required \emph{all} NQ8-to-NQ10 radial measures to pass; its trigger failed, so no extra jobs were run. No calculation above L800 was started.

The next scientific stage should first secure a radial reference at L800 by modest $h$ refinement or higher order, then investigate frequency-resolved high-degree field tails and their shallow spatial behavior before selecting a further angular ladder. That is a proposal, not executed work. Changing the frequency filter or measurement window requires its own convergence check.

\newpage
\section{Execution, memory coordination and evidence}
Eight sequential physical jobs completed with unchanged declared inputs, verified cleanup, and 132 data outputs plus eight logs hash-checked. The finalized producer snapshot begins \texttt{35ac18a3}; the separately supervised checker returned zero with all 241 inputs unchanged and cleanup verified. Its report begins \texttt{c9a162ad}; operational success is distinct from its failed accuracy criteria.

The scientific source changes were limited to the isolated bounded high-degree controls, explicit workspace policy and metadata. The radial/modal/angular mathematics were retained. One frozen native binary was built on Tunya; its SHA-256 begins \texttt{d7ffb6cb}. Luna high implemented and Sol 6.1 high reviewed the plan, source, build, input closures, operators, executed evidence and numerical findings independently.

The user requested a combined active-memory ceiling of 900 GB on Tunya. Coordination with CCDSM and DSpecM3D Implementation preserved Pilot's separate 8 GB reservation and the existing CCDSM limits, and admitted DSpecM1D to its own 40 GB/12-task/zero-swap/eight-core slice. One 32 GiB build, 16 GiB numerical campaign or 8 GiB checker ran at a time. All descendants were inside kernel-owned services; the original identity-aware user-space supervisor also remained active.

The same shared active-RAM guardian was extended only to discover the new slice/prefix. Its active-memory definition excludes reclaimable Linux memory and reported evictable ARC, with protective thresholds at 825/875 GB below the 900 GB boundary. The guardian's preserved old epoch is distinct from its new observed-peak epoch. The maximum numerical resident memory observed by the user-space supervisor was RSS GiB; this is not kernel aggregate accounting or a machine-wide kernel memory cap.

At completion, all three DSpecM1D services were inactive and our slice contained zero tasks. Its kernel-accounted peak was 7.287 GB, with zero memory-limit or OOM events. Shared active RAM was 81.949 GB and its peak observed in the new guardian epoch was 90.397 GB. These are sampled active-memory estimates under the stated definition, not cache-inclusive occupancy. CCDSM and DSpecM3D Implementation were notified to release our logical 40 GB reservation; their other work and guardian were left running.

The staging checks initially caught two incorrectly declared external build paths and stopped before numerical execution. The old plans/packages and stopped attempt were preserved; a reviewed recovery added the corrected terminal plan without overwriting the original. Local retrieval-path/receipt-writer failures were also retained. No scientific calculation was repeated to repair metadata or retrieval.

Figures are reproducible from hash-bound raw files. The saved wedges have 68 along-path and 10 cross-path angles, separately from the solved L; interpolation is for display and does not resolve the full high-L angular bandwidth. Linear percentage difference panels use the full displayed DSpec-section maximum, whereas the tables use the named cohort norms/maxima. Display saturation and unmatched source-shell traces are labelled.

All changes/results remain isolated. No production integration, new optimization, broader band/window campaign or unapproved high-degree extension followed.

\medskip\noindent\small Full hashes and supporting source are retained in \texttt{reports/report-data.json}, the delivery manifest, the independent AFTER-results review, and the producer/checker raw evidence directories. Scientific conclusions apply only to the stated fixture and sampled quantities.
\end{document}
'''
for key,value in {"DEGREE":degree_rows,"RADIAL":radial_rows,"CROSS":cross_rows,
                  "ORDER":order_rows,"DEPTH":depth_rows,"RSS":f"{peak_rss:.4f}"}.items():
    text=text.replace(key,value)
(HERE/"angular-convergence.tex").write_text(text)
print(json.dumps(summary,indent=2))
