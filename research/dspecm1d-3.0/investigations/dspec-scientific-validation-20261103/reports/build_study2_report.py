#!/usr/bin/env python3
"""Produce a finite-convergence report from sealed actual evidence only."""
from pathlib import Path
import csv
import hashlib
import json
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
S=ROOT/"study2-v1"
RAW=S/"evidence/producer-raw-v1"
REPORT=S/"evidence/checker-terminal-v2/v2/task/study2-report.json"
METRICS=S/"reviews/study2-independent-raw-metrics-20261103.json"
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(REPORT)=="aedda47c7cfe830021f93e08320c7922b617c6b31c341cc8382fa29e60fe1e61"
report=json.loads(REPORT.read_text())
metrics=json.loads(METRICS.read_text())
assert report["study2_state_sha256"]==metrics["source_state_sha256"]==sha(RAW/"state.json")
manifest=json.loads((RAW/"manifest.json").read_text())
for name, item in manifest.items():
    p=RAW/name
    assert p.stat().st_size==item["bytes"] and sha(p)==item["sha256"]
c=metrics["comparisons"]
def val(name,kind):
    return 100*(c[name]["waveform"]["velocity_m_s"]["relative_change"] if kind=="wave"
        else c[name]["kernel_subsets"][kind]["relative_change"])
def row(label,name):
    return f'{label} & {val(name,"wave"):.6f} & {val(name,"far_solid_bulk_100km"):.6f} & {100*abs(c[name]["measurement_norm_relative_change"]):.6f} '+r"\\"
radial="\n".join(row(label,key) for label,key in [
    (r"$h_0\to0.75h_0$","h_base_to075"),(r"$0.75h_0\to0.5625h_0$","h075_to05625"),
    ("NQ6 to NQ8","nq6_to8"),("NQ8 to NQ10","nq8_to10")])
degree="\n".join(row(label,key) for label,key in [
    ("DSpec L100 to L150","D100_to150"),("DSpec L150 to L200","D150_to200"),
    ("DSpec L200 to L300","D200_to300"),("YSpec L100 to L150","Y100_to150"),
    ("YSpec L150 to L200","Y150_to200"),("YSpec L200 to L300","Y200_to300")])
matched="\n".join(row(f"L{l}",f"DY{l}") for l in (100,150,200,300))
# An independent window-only diagnostic: same nonzero-window sample indices, no weighting.
def trace(job):
    return list(csv.DictReader((RAW/job/"data/trace.csv").open()))
a,b=trace("l200-dspec-sparse-bridge"),trace("l300-dspec")
ab=[(float(x["velocity_m_s"]),float(y["velocity_m_s"])) for x,y in zip(a,b) if float(x["window"])>0]
window_change=100*math.sqrt(sum((x-y)**2 for x,y in ab)/sum(x*x for x,y in ab))
assert len(ab)>0 and window_change>0.1
plt.rcParams.update({"font.family":"serif","font.size":10,"axes.spines.top":False,
    "axes.spines.right":False,"pdf.fonttype":42})
fig,axes=plt.subplots(1,2,figsize=(7,3.2),constrained_layout=True)
for ax,kind,title in zip(axes,("wave","far_solid_bulk_100km"),("Whole velocity trace","89 far-solid-interior probes")):
    for code,color in (("D","#244f75"),("Y","#a12432")):
        ys=[val(f"{code}{x}_to{y}",kind) for x,y in ((100,150),(150,200),(200,300))]
        ax.plot(range(3),ys,"o-",color=color,label="DSpec" if code=="D" else "YSpec",lw=1.2,ms=4)
    ax.axhline(.1,color=".35",ls="--",lw=.8,label="0.1% target")
    ax.set_xticks(range(3),["100 to 150","150 to 200","200 to 300"],rotation=18)
    ax.set_yscale("log");ax.set_ylabel("Relative vector change (%)")
    ax.set_xlabel("Physical degree cutoff L");ax.set_title(title,fontsize=10)
    ax.legend(frameon=False,fontsize=8)
fig.savefig(HERE/"figures/study2-degree-controls.pdf",bbox_inches="tight")
fig.savefig(HERE/"figures/study2-degree-controls.png",dpi=150,bbox_inches="tight")
plt.close(fig)
summary={"checker_report_sha256":sha(REPORT),"independent_metrics_sha256":sha(METRICS),
    "producer_manifest_sha256":sha(RAW/"manifest.json"),
    "window_nonzero_sample_count":len(ab),"L200_to300_window_velocity_relative_change_percent":window_change,
    "interpretation":"radial observed changes below 0.1%; degree ladder does not establish convergence"}
(HERE/"study2-report-data.json").write_text(json.dumps(summary,indent=2)+"\n")
text=r'''\documentclass[11pt,a4paper]{article}
\usepackage{alex-report}
\ReportSetup{DSpecM1D scientific validation}{Study 2}{3 November 2026}
\hypersetup{pdftitle={DSpecM1D study 2: finite convergence controls}}
\begin{document}
\ReportTitle{Finite convergence controls}{Radial resolution is stable; angular cutoff is not yet stabilized}
\textbf{Finding.} The observed radial-refinement and NQ changes are below 0.1\% for the waveform, measurement normalization and selected far-solid-interior kernel vector. The physical angular cutoff is a different limitation: L200 to L300 changes the whole velocity trace by 4.16\% and the selected far-interior kernel vector by 165.90\%. YSpec shows essentially the same degree dependence. Close agreement between codes at fixed L does not establish convergence to the continuum problem.

\section{Case, controls and meaning of the numbers}
This is the same attenuated anisotropic PREM transverse-velocity traveltime fixture as Study 1: source depth 29.5 km; solid-side receiver at the 3 km solid--fluid boundary; 60-degree separation; 1000--1200 s window with 40 s taper; existing 2/4/8/10 mHz filter. All 35 bins, full gravity and full radial systems are retained. The kernel parameter is $\ln V_s$ at fixed $V_p$, density and effective $Q$.

Eleven new sequential solver jobs were run on Tunya in isolated directories. DSpec used eight frequency workers and YSpec eight outer workers, each with one inner library thread. Physical degree cutoffs were L100, L150, L200 and L300. YSpec's radial ODE tolerance was fixed at $10^{-7}$. Each unperturbed case used its own waveform-derived adjoint weights. The L200 reference comes from the accepted first study.

All percentages below are $100\|b-a\|_2/\|a\|_2$ for the named discrete sample vectors, except the normalization column, which is $100|N_b/N_a-1|$. The far-interior vector contains 89 fixed probes in solid material, away from material interfaces and at least 100 km from both source and receiver; source-shell traces are excluded. These are sample-vector diagnostics, not three-dimensional $L^2$ or $H^1$ norms or universal pointwise error bounds.

\section{Radial resolution and polynomial/quadrature order}
\begin{center}\small\begin{tabular}{lrrr}\toprule
Change at L200 & Velocity (\%) & Far kernel (\%) & $|\Delta N|$ (\%)\\\midrule
RADIAL
\bottomrule\end{tabular}\end{center}
The normalized mesh-step parameter is $h_0=0.02035382429229197$, not a length in metres. The mesh has 59, 75 and 95 elements on the three $h$ levels; maximum physical element widths are 125.472, 94.104 and 72.855 km. NQ6, NQ8 and NQ10 use exactly the same 59-element boundaries at $h_0$. Thus increasing NQ is not confounded with a mesh change. NQ6 versus NQ10 changes the selected far-interior vector by 0.007464\%; NQ8 versus NQ10 is smaller still.

These controls support the current radial discretization for this finite-L fixture. They do not transfer the 0.1\% statement to an untested angular cutoff, every point or a singular-limit norm. YSpec's sampled radial mesh is an output grid; it is not an independent YSpec ODE refinement control.

\newpage
\section{Physical angular cutoff: unresolved convergence}
\begin{center}\small\begin{tabular}{lrrr}\toprule
Change & Velocity (\%) & Far kernel (\%) & $|\Delta N|$ (\%)\\\midrule
DEGREE
\bottomrule\end{tabular}\end{center}
\begin{figure}[H]\centering
\includegraphics[width=\linewidth]{figures/study2-degree-controls.pdf}
\caption{Actual adjacent cutoff changes with each case's own waveform-derived adjoint weights. Curves for the two codes nearly overlap. The dashed line is a 0.1\% observed-change target, not an established error estimator. Fixed probes are reused; their vector norm has no volume weights.}
\end{figure}
The L200-to-L300 DSpec velocity change within the nonzero measurement-window samples alone is WINDOW\%, so the discrepancy is not confined to times outside the measurement. The measurement-rate normalization changes by $-1.43944$\%, and the adjoint-source vector changes by 0.794524\%. The selected far-interior kernel changes by 165.90\%, solid material-interface samples by 146.81\%, and source-shell inner/outer samples by 193.53\%. The latter remain contact-limited regular-volume samples.

The far-interior L200 vector norm is $8.71906\times10^{-15}$ and the L300 norm is $2.00077\times10^{-14}$; their difference norm is $1.44648\times10^{-14}$ in the sampled kernel units. This is a material change, not merely a percentage computed at one near-zero point. It does not prove divergence. The nonmonotone finite ladder does not establish an asymptotic rate or a sufficient final L. No L greater than 300 was run.

\newpage
\section{Agreement between the two implementations}
\begin{center}\small\begin{tabular}{lrrr}\toprule
Matched cutoff & Velocity (\%) & Far kernel (\%) & $|\Delta N|$ (\%)\\\midrule
MATCHED
\bottomrule\end{tabular}\end{center}
The denominator here is DSpec's vector norm, and both implementations use their own waveform-derived weights. These small discrepancies establish agreement for the matched finite-cutoff calculations on the supported sample cohort. They do not establish that the common finite-cutoff answer is accurate to 0.1\% against a converged angular solution.

Source-shell traces unavailable through YSpec's public interface are not included in the matched far-interior comparison. The Study 1 source-radius trace mismatch remains a separate interpretation issue. The receiver functional is unchanged on the solid side of the physical solid--fluid interface. No smoothing, new $H^1$ gradient, source regularization, receiver averaging or singular-contact correction was implemented.

\section{Finite samples and the source/receiver question}
All saved native section samples are finite. This cannot establish bounded continuum values near ideal point sources or receivers. Neither this study nor Study 1 proves or disproves a singular limit, global $H^1$ regularity or a valid pointwise product at contact. The source-shell regular-volume samples and any omitted distributional contact contributions must not be interchanged.

The output sections deliberately use a common sparse grid of 68 along-path and 10 cross-path angles while the actual solved degree changes up to L300. Sparse section extrema remain unmasked in the saved evidence, but they are sampled extrema, not global angular maxima. The 168 fixed probes provide exact matched physical positions for the control comparisons. The historical vertical-displacement energy example was not rerun and is not certified by this traveltime study.

\section{Validation record and stopping point}
All 11 jobs completed cleanly, with unchanged declared inputs and verified cleanup. The root retrieval rehashed 197 outputs/logs plus the exact finalized state. The standalone checker passed after one reviewed provenance-only correction: the absolute model path in \texttt{parameters.txt} necessarily differs between the isolated Study 1 and Study 2 folders. Both model hashes and every non-path byte were checked. The failed first checker state and log are retained; no numerical job was rerun.

The checker status PASS means that the planned finite controls and evidence checks passed. It does not mean the angular solution converged. Sol 6.1 independently recomputed the trends from raw outputs, in addition to its plan, source, build and terminal reviews. Resource peaks were 9 threads and approximately 1.82 GiB under the accepted user-space watchdog; these are not kernel aggregate accounting.

The supported conclusion is radial stability and close cross-code agreement at fixed L, with unresolved angular-cutoff convergence. Study 1's supported smooth-perturbation directional derivative check remains valid within its own low-degree scope. All changes and results remain isolated. Production integration, higher-cutoff work and further optimization have not been started.

\medskip\noindent\small Evidence: \texttt{study2-v1/evidence/producer-raw-v1/}; standalone checker report and independent raw metrics are retained alongside the editable report generator. Frozen producer state SHA-256 begins \texttt{9cd09f3e}; checker report begins \texttt{aedda47c}. The delivery manifest records full hashes.
\end{document}
'''
text=text.replace("RADIAL",radial).replace("DEGREE",degree).replace("MATCHED",matched).replace("WINDOW",f"{window_change:.6f}")
(HERE/"study2.tex").write_text(text)
print(json.dumps(summary,indent=2))
