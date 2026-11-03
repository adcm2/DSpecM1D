#!/usr/bin/env python3
"""Check the Phase C direct/common L3 kernel against frozen Stage 2 FD traces."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path

def read_json(path):
    obj=json.loads(Path(path).read_text())
    def finite(v):
        if isinstance(v,bool) or v is None or isinstance(v,str): return
        if isinstance(v,(int,float)) and not math.isfinite(float(v)): raise ValueError(f"nonfinite value in {path}")
        if isinstance(v,list):
            for x in v: finite(x)
        if isinstance(v,dict):
            for x in v.values(): finite(x)
    finite(obj); return obj

def rows(path):
    with Path(path).open(newline="") as f: return list(csv.DictReader(f))

def require(ok,msg):
    if not ok: raise ValueError(msg)

def exact_csv_column(a,b,key):
    require(len(a)==len(b),f"row count differs for {key}")
    for i,(x,y) in enumerate(zip(a,b)):
        require(key in x and key in y,f"missing {key} in baseline comparison")
        require(float(x[key])==float(y[key]),f"exact grid mismatch row {i} column {key}")

def relative_l2_csv(a,b,key,limit=1e-11):
    require(len(a)==len(b),f"row count differs for {key}")
    err=0.0;norm=0.0
    for i,(x,y) in enumerate(zip(a,b)):
        require(key in x and key in y,f"missing {key} in baseline comparison")
        xv=float(x[key]);yv=float(y[key])
        require(math.isfinite(xv) and math.isfinite(yv),f"nonfinite {key} at row {i}")
        err+=(xv-yv)**2;norm+=yv*yv
    rel=math.sqrt(err/norm) if norm else (0.0 if err==0 else math.inf)
    require(rel<=limit,f"baseline {key} relative L2 mismatch {rel:.17g} exceeds {limit}")
    return rel

def compare_csv(candidate,reference,grid_columns,physical_columns):
    a=rows(candidate);b=rows(reference)
    require(len(a)==len(b),f"row count differs: {candidate} vs {reference}")
    for key in grid_columns: exact_csv_column(a,b,key)
    return {key:relative_l2_csv(a,b,key) for key in physical_columns}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidate",type=Path,required=True)
    ap.add_argument("--base",type=Path,required=True)
    ap.add_argument("--plus",type=Path,required=True)
    ap.add_argument("--minus",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    c=read_json(a.candidate/"result.json");b=read_json(a.base/"result.json")
    p=read_json(a.plus/"result.json");m=read_json(a.minus/"result.json")
    for name,d in (("C direct/common",c),("base",b),("plus",p),("minus",m)):
        require(d.get("degree")==3 and d.get("nq")==8 and d.get("refinement")==0.05,
                f"{name} is not the frozen L3/NQ8/h=.05 audit")
        require(d.get("frequency_count")==35 and d.get("first_bin")==7 and d.get("last_bin")==41,
                f"{name} frequency grid changed")
        require(d.get("fft_samples")==2048 and d.get("fft_period_s")==4096,
                f"{name} sample grid changed")
        require(d.get("native_side_count")==157 and d.get("along_angles")==91 and d.get("across_angles")==0,
                f"{name} audit inventory changed")
        require(abs(float(d.get("gravity_max_abs_nd_error",math.inf)))<=1e-10,
                f"{name} gravity check failed")
        require(abs(float(d.get("time_shift_adjoint_identity",0))-1)<=1e-10,
                f"{name} time-shift identity failed")
    require(c.get("backend")=="dspec" and c.get("phase_c_direct_adjoint") is True and
            c.get("phase_c_combined_reconstruction") is True and
            c.get("phase_c_validate_reconstruction") is True,
            "C audit did not use the validated direct-adjoint/common-reconstruction route")
    eps=float(p.get("epsilon_perturbation",0))
    require(eps>0 and float(m.get("epsilon_perturbation",0))==-eps and
            float(b.get("epsilon_perturbation",math.nan))==0 and
            float(c.get("epsilon_perturbation",math.nan))==0,
            "frozen central-difference pair is not matched")
    trace_parity=compare_csv(a.candidate/"trace.csv",a.base/"trace.csv",
                ["time_s","window"],["velocity_m_s","rate_m_s2","adjoint_source"])
    spectrum_parity=compare_csv(a.candidate/"spectrum.csv",a.base/"spectrum.csv",
                ["bin","frequency_hz"],["real_m","imag_m"])
    component_parity=compare_csv(a.candidate/"components.csv",a.base/"components.csv",
                ["bin"],["up_re","up_im","north_re","north_im","east_re","east_im"])
    base=rows(a.base/"trace.csv");plus=rows(a.plus/"trace.csv");minus=rows(a.minus/"trace.csv")
    require(len(base)==len(plus)==len(minus)==1201,"frozen F8 saved-trace sample count changed")
    dt=2.0
    for i,(x,y,z) in enumerate(zip(base,plus,minus)):
        expected=2.0*i
        require(float(x["time_s"])==expected and float(y["time_s"])==expected and float(z["time_s"])==expected,
                f"frozen saved-trace time grid changed at row {i}")
    fd=sum(dt*float(x["adjoint_source"])*(float(y["velocity_m_s"])-float(z["velocity_m_s"]))/(2*eps)
           for x,y,z in zip(base,plus,minus))
    kernel=float(c["integrated_profile_kernel"])
    rel=abs(fd-kernel)/max(abs(fd),abs(kernel),1e-300)
    require(rel<=2e-4,f"C direct/common kernel failed existing 2e-4 FD gate: {rel:.17g}")
    report={"status":"pass","route":"phase C direct measurement adjoint + combined common reconstruction",
            "degree":3,"nq":8,"refinement":0.05,"frequency_count":35,"radial_sides":157,
            "angular_points":91,"frozen_epsilon":eps,"finite_difference":fd,
            "integrated_common_kernel":kernel,"fd_relative_error":rel,"fd_limit":2e-4,
            "baseline_trace_spectrum_components_relative_l2": {"trace":trace_parity,"spectrum":spectrum_parity,"components":component_parity},
            "gravity_max_abs_nd_error":float(c["gravity_max_abs_nd_error"]),
            "time_shift_adjoint_identity":float(c["time_shift_adjoint_identity"])}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__": main()
