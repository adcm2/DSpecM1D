#!/usr/bin/env python3
"""Validate the sealed angular-v1 degree ladder and its radial controls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import struct
from pathlib import Path
from typing import Any


JOB_IDS = (
    "angular-l400-dspec-nq6", "angular-l400-yspec-nq6",
    "angular-l600-dspec-nq6", "angular-l600-yspec-nq6",
    "angular-l800-dspec-nq6", "angular-l800-yspec-nq6",
    "angular-l800-dspec-nq8", "angular-l800-dspec-nq10",
)
DEGREES = (200, 300, 400, 600, 800)
TARGET = 0.001
FILTER_MHZ = [2, 4, 8, 10]
WINDOW = (1000.0, 1200.0)
DEPTH_STRATA = {10: 10, 24: 10, 29: 10, 30: 10, 50: 10, 100: 11, 300: 14, 500: 14}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def finite(value: Any, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} is non-finite")
    return number


def load_helper(path: Path):
    spec = importlib.util.spec_from_file_location("angular_study_helpers", path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot import the hash-bound Study 1 helper module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_manifest(root: Path, path: Path) -> dict[str, Any]:
    manifest = read_json(path)
    checked = 0
    for group, items in manifest.get("references", {}).items():
        for entry in items:
            f = root / entry["path"]
            if not f.is_file() or f.stat().st_size != entry["bytes"] or sha(f) != entry["sha256"]:
                raise ValueError(f"reference file mismatch in {group}: {f}")
            checked += 1
    return {"manifest_sha256": sha(path), "verified_files": checked,
            "study2_after_review_sha256": manifest.get("study2_after_review_sha256")}


def verify_state(run_root: Path, plan_path: Path, state_path: Path) -> dict[str, dict[str, Any]]:
    state = read_json(state_path)
    if state.get("status") != "completed" or state.get("plan_sha256") != sha(plan_path):
        raise ValueError("producer state is not the completed snapshot of the exact angular-v1 plan")
    by_id = {j.get("id"): j for j in state.get("jobs", [])}
    if set(by_id) != set(JOB_IDS) or len(by_id) != len(JOB_IDS):
        raise ValueError("producer job inventory differs from the exact eight-job queue")
    for job_id in JOB_IDS:
        entry = by_id[job_id]
        if (entry.get("status") != "completed" or entry.get("returncode") != 0 or
                entry.get("failure_reasons") or entry.get("inputs_unchanged") is not True or
                not entry.get("cleanup", {}).get("verified_gone")):
            raise ValueError(f"{job_id} lacks clean terminal completion evidence")
        task_dir = run_root / job_id
        log = task_dir / entry["log"]
        if not log.is_file() or sha(log) != entry.get("log_sha256"):
            raise ValueError(f"{job_id} producer log hash mismatch")
        for rel in entry.get("outputs", []):
            out = task_dir / rel
            expected = entry.get("outputs_sha256", {}).get(rel)
            if not out.is_file() or not expected or sha(out) != expected:
                raise ValueError(f"{job_id} output missing or hash mismatch: {rel}")
    return by_id


def rel_l2(reference: list[complex | float], candidate: list[complex | float], label: str) -> dict[str, Any]:
    if len(reference) != len(candidate):
        raise ValueError(f"{label} vector lengths differ")
    ref = math.sqrt(sum(abs(x) ** 2 for x in reference))
    err = math.sqrt(sum(abs(y - x) ** 2 for x, y in zip(reference, candidate)))
    ratio = err / ref if ref else (0.0 if err == 0 else None)
    maxref = max((abs(x) for x in reference), default=0.0)
    maxerr = max((abs(y - x) for x, y in zip(reference, candidate)), default=0.0)
    maxratio = maxerr / maxref if maxref else (0.0 if maxerr == 0 else None)
    if ref == 0 and err != 0:
        raise ValueError(f"{label} differs from an exact-zero reference")
    if maxref == 0 and maxerr != 0:
        raise ValueError(f"{label} max differs from an exact-zero reference")
    return {"count": len(reference), "reference_l2": ref, "difference_l2": err,
            "relative_l2": ratio, "reference_max_abs": maxref,
            "maximum_absolute_difference": maxerr, "max_relative_to_reference_max": maxratio,
            "zero_reference_class": ref == 0.0}


def ratio_pass(metric: dict[str, Any], target: float, label: str) -> bool:
    value = metric.get("relative_l2")
    if value is None:
        return metric["zero_reference_class"] and metric["difference_l2"] == 0.0
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{label} has invalid relative L2")
    return value <= target


def compare_scalar(a: float, b: float, label: str) -> dict[str, Any]:
    a, b = finite(a, label + " candidate"), finite(b, label + " reference")
    diff = abs(a - b)
    ratio = diff / abs(b) if b != 0 else (0.0 if diff == 0 else None)
    if b == 0 and diff != 0:
        raise ValueError(f"{label} differs from an exact-zero reference")
    return {"candidate": a, "reference": b, "absolute_difference": diff,
            "relative_difference": ratio, "zero_reference_class": b == 0}


def check_trace(data: Path, label: str) -> dict[str, Any]:
    trace = rows(data / "trace.csv")
    required = {"time_s", "velocity_m_s", "rate_m_s2", "window", "adjoint_source"}
    if not trace or set(trace[0]) != required:
        raise ValueError(f"{label} trace schema differs")
    times, u, w = [], [], []
    for row in trace:
        t = finite(row["time_s"], label + " time")
        v = finite(row["velocity_m_s"], label + " velocity")
        win = finite(row["window"], label + " taper")
        finite(row["rate_m_s2"], label + " rate")
        finite(row["adjoint_source"], label + " adjoint source")
        if win < 0 or win > 1 or (times and t <= times[-1]):
            raise ValueError(f"{label} trace has invalid time/window ordering")
        times.append(t); u.append(v); w.append(win)
    ids = [i for i, t in enumerate(times) if WINDOW[0] <= t <= WINDOW[1]]
    active = [i for i in ids if w[i] > 0]
    if len(ids) != 101 or len(active) != 99:
        raise ValueError(f"{label} does not have the frozen 101/99-point measurement window")
    if any(times[i] < WINDOW[0] or times[i] > WINDOW[1] for i in ids):
        raise ValueError(f"{label} window mismatch")
    return {"rows": trace, "times": times, "velocity": u, "window": w,
            "window_indices": ids, "active_window_indices": active}


def waveform_metric(reference: dict[str, Any], candidate: dict[str, Any], weighted: bool) -> dict[str, Any]:
    key = "window_indices" if weighted else "active_window_indices"
    ixs = reference[key]
    if ixs != candidate[key] or reference["times"] != candidate["times"] or reference["window"] != candidate["window"]:
        raise ValueError("measurement time/window grid differs")
    if weighted:
        a = [reference["velocity"][i] * reference["window"][i] for i in ixs]
        b = [candidate["velocity"][i] * candidate["window"][i] for i in ixs]
    else:
        a = [reference["velocity"][i] for i in ixs]
        b = [candidate["velocity"][i] for i in ixs]
    return rel_l2(a, b, "waveform")


def check_angles(data: Path, label: str) -> dict[str, Any]:
    out = {}
    node_count = len(rows(data / "nodes.csv"))
    for part, expected in (("along", 68), ("across", 10)):
        angle_rows = rows(data / f"{part}-angles.csv")
        if len(angle_rows) != expected:
            raise ValueError(f"{label} has wrong {part} angle count")
        values = (data / f"{part}.bin").read_bytes()
        if len(values) != node_count * expected * 8:
            raise ValueError(f"{label} {part} field dimensions differ from node/angle inventory")
        n = len(values) // 8
        unpacked = struct.unpack(f"{n}d", values)
        if any(not math.isfinite(x) for x in unpacked):
            raise ValueError(f"{label} {part} section fields contain nonfinite values")
        out[part] = angle_rows
    return out


def angular_workspace_bytes(degree: int, nodes: int, angles: int, cache: int) -> tuple[int,int]:
    # Exact existing visitor accounting, with its fixed tile and packing sizes.
    matrix=32*(degree+1+16)*16
    accumulation=16*nodes*2*6*16
    weight=25*(degree+1)*16*16
    build=25*(degree+1)*8
    recurrence=2*(degree+3)*8
    eigen=2*1024*1024
    generic=cache+matrix+accumulation+weight+build+recurrence+eigen
    full_fields=angles*nodes*2*6*16
    coefficients=4*nodes*(degree+1)*8
    real_tile=(degree+1+4*32)*64*8
    real_weights=(degree+1)*64*8
    real=cache+full_fields+coefficients+real_tile+real_weights+build+recurrence+eigen
    return max(generic,real), full_fields


def check_result(data: Path, job: dict[str, Any], helper, probes: list[dict[str, str]],
                 source_xyz: tuple[float, float, float], receiver_xyz: tuple[float, float, float]) -> dict[str, Any]:
    meta = read_json(data / "result.json")
    backend, degree, nq = job["backend"], job["degree"], job["nq"]
    if (meta.get("backend") != backend or meta.get("degree") != degree or meta.get("nq") != nq or
            meta.get("phase_d_angular_convergence_v1") is not True or
            meta.get("phase_d_study2_control") is not True or
            meta.get("phase_d_study2_sparse_sampling") is not True or
            meta.get("phase_d_study2_sampling_degree") != 100 or
            meta.get("phase_d_study2_absolute_h") is not (backend == "dspec") or
            meta.get("phase_d_batched_angular") is not (backend == "dspec") or
            meta.get("frequency_count") != 35 or meta.get("first_bin") != 7 or meta.get("last_bin") != 41 or
            meta.get("along_angles") != 68 or meta.get("across_angles") != 10 or
            meta.get("fixed_probe_count") != 168 or meta.get("fixed_probe_contact_terms_included") is not False or
            meta.get("kernel_adjoint_weight_source") != "this mesh's own waveform-derived adjoint source"):
        raise ValueError(f"{job['id']} metadata differs from exact angular-v1 workload")
    scope = meta.get("phase_d_angular_scope")
    expected_scope = {"case":"study2-transverse-velocity-traveltime", "filter_mHz":FILTER_MHZ,
                      "fft_bins":[7,41], "window_s":[1000,1200], "taper_s":40,
                      "component":"T", "quantity":"velocity", "weight_source":"own_waveform"}
    if scope != expected_scope:
        raise ValueError(f"{job['id']} filter/window scope metadata differs")
    node_count=int(meta.get("native_side_count",0)); observed_nodes=len(rows(data/'nodes.csv'))
    if not 0<node_count<=350 or node_count!=observed_nodes:
        raise ValueError(f"{job['id']} native-side metadata differs from the node inventory or reviewed bound")
    for key, upper in (("gravity_max_abs_nd_error",1e-10),):
        value=finite(meta.get(key),job['id']+"."+key)
        if value<0 or value>upper: raise ValueError(f"{job['id']} {key} fails")
    ident=finite(meta.get("time_shift_adjoint_identity"),job['id']+".time_shift_identity")
    if ident<=0 or abs(ident-1)>1e-10: raise ValueError(f"{job['id']} time-shift identity fails")
    if finite(meta.get("measurement_rate_norm"),job['id']+".measurement_rate_norm")<=0:
        raise ValueError(f"{job['id']} measurement norm is not positive")
    if meta.get("phase_d_selected_h") != 0.02035382429229197:
        raise ValueError(f"{job['id']} selected mesh-step parameter differs from fixed h")
    if meta.get("fixed_probe_regular_part_only") is not True:
        raise ValueError(f"{job['id']} fixed probes are not labeled regular-part samples")
    if backend=="dspec":
        if (meta.get("native_angular_route")!="real-legendre" or
            meta.get("native_radial_route")!="batched-direct-index" or
            meta.get("native_derivative_weights_mode")!="native-dref-over-jacobian" or
            meta.get("modal_operator_route")!="legacy" or meta.get("modal_operator_oracle") is not False or
            meta.get("phase_e_angular_oracle") is not False or meta.get("phase_c_zero_bin_skips") is not True or
            meta.get("frequency_workers")!=8 or meta.get("solve_workers")!=1 or meta.get("workers")!=1 or
            meta.get("phase_c_validate_direct_adjoint") is not False or
            meta.get("phase_c_validate_reconstruction") is not False or
            meta.get("phase_d_common_validation") is not False):
            raise ValueError(f"{job['id']} DSpec route/worker/oracle metadata differs")
        geom=finite(meta.get("phase_d_actual_cache_bytes"),job['id']+".geometry_bytes")
        peak=finite(meta.get("phase_d_actual_cache_plus_temporary_peak_bytes"),job['id']+".combined_bytes")
        product=finite(meta.get("phase_c2_product_storage_bytes"),job['id']+".product_bytes")
        angle_count=int(meta.get("along_angles",0))+int(meta.get("across_angles",0))
        expected_geometry=78*(5*11*(degree+1)*8+952)
        expected_peak,full_fields=angular_workspace_bytes(degree,node_count,angle_count,expected_geometry)
        expected_product=angle_count*node_count*35*16
        if meta.get("phase_d_angular_memory_policy_applied") is not True or \
           meta.get("phase_d_angular_geometry_cap_bytes")!=32*1024*1024 or \
           meta.get("phase_d_angular_cache_plus_temporary_cap_bytes")!=64*1024*1024 or \
           meta.get("phase_d_batched_cache_max_bytes")!=64*1024*1024 or \
           not 0<node_count<=350 or node_count!=observed_nodes or \
           geom!=expected_geometry or geom>32*1024*1024 or peak<=0 or peak>64*1024*1024 or \
           product!=expected_product or product<=0 or product>32*1024*1024 or peak!=expected_peak:
            raise ValueError(f"{job['id']} exceeds reviewed expanded angular memory bounds")
        if finite(meta.get("phase_d_wigner_recursion_table_scratch_bytes"),job['id']+".Wigner_scratch") != 2*(degree+3)*8:
            raise ValueError(f"{job['id']} Wigner local recurrence scratch report differs")
        if nq==10 and int(meta.get("native_side_count",0))>350:
            raise ValueError(f"{job['id']} radial-side count exceeds reviewed memory bound")
        weights=helper.read_weights(data)
        integrated=helper.fixed_table(data,probes,backend,source_xyz,receiver_xyz)
        helper.validate_dspec_bin_products(data,weights,integrated)
        helper.validate_fields(data,weights,probes)
        products=read_unit_products(data,weights,probes)
        cache={"policy":"geometry32MiB/cache+temporary64MiB", "geometry_bytes":geom,
               "expected_geometry_bytes":expected_geometry,"node_count":node_count,
               "observed_node_rows":observed_nodes,
               "cache_plus_temporary_peak_bytes":peak,"expected_formula_peak_bytes":expected_peak,
               "product_bytes":product,"expected_product_bytes":expected_product,
               "full_angle_fields_bytes":full_fields,
               "wigner_recursion_table_scratch_bytes":2*(degree+3)*8}
    else:
        if finite(meta.get("yspec_radial_solver_tolerance"),job['id']+".ODE_tol")!=1e-7:
            raise ValueError(f"{job['id']} YSpec ODE tolerance differs from frozen control")
        if meta.get("workers")!=8 or meta.get("frequency_workers")!=1 or meta.get("solve_workers")!=8:
            raise ValueError(f"{job['id']} YSpec worker metadata differs")
        if meta.get("phase_d_angular_memory_policy_applied") is not False:
            raise ValueError(f"{job['id']} must not claim the DSpec-only memory policy")
        if meta.get("phase_d_batched_cache_max_bytes")!=32*1024*1024:
            raise ValueError(f"{job['id']} YSpec legacy cache alias differs from the default combined cap")
        weights=helper.read_weights(data)
        integrated=helper.fixed_table(data,probes,backend,source_xyz,receiver_xyz)
        products=None; cache=None
    trace=check_trace(data,job['id'])
    angles=check_angles(data,job['id'])
    table=helper.fixed_table(data,probes,backend,source_xyz,receiver_xyz)
    boundary=[]
    for row in rows(data/'mesh-elements.csv'):
        boundary.append((int(row['layer']),finite(row['lower_radius_m'],'lower boundary'),
                         finite(row['upper_radius_m'],'upper boundary'),int(row['fluid'])))
    return {"metadata":meta,"trace":trace,"weights":weights,"integrated":table,
            "unit_products":products,"angles":angles,"boundaries":boundary,"memory":cache}


def compare_tables(ref: dict[str, Any], cand: dict[str, Any], probes: list[dict[str,str]]) -> dict[str, Any]:
    ids=sorted(ref.keys() & cand.keys())
    selected={"all_numeric":[],"far_solid_interior_100km":[],"source_shell":[],"interface":[],
              "positive_distance_near_100km":[]}
    byid={p['probe_id']:p for p in probes}
    for pid in ids:
        a,b=ref[pid],cand[pid]
        av,bv=a['numeric_kernel'],b['numeric_kernel']
        if not isinstance(av,float) or not isinstance(bv,float): continue
        selected['all_numeric'].append(pid)
        probe=byid[pid]; ds=float(a['source_distance_3d_m']); dr=float(a['receiver_distance_3d_m'])
        shell=probe['layer']=='9' and probe['radial_trace'] in ('inner','outer')
        if shell: selected['source_shell'].append(pid)
        if a['location_class']=='material_interface': selected['interface'].append(pid)
        if min(ds,dr)<100_000 and min(ds,dr)>0: selected['positive_distance_near_100km'].append(pid)
        if (probe['material_side']=='solid' and a['location_class']=='bulk' and
            b['location_class']=='bulk' and not shell and ds>=100_000 and dr>=100_000):
            selected['far_solid_interior_100km'].append(pid)
    out={}
    for name, pids in selected.items():
        out[name]=rel_l2([float(ref[p]['numeric_kernel']) for p in pids],
                         [float(cand[p]['numeric_kernel']) for p in pids],name)
        out[name]['count']=len(pids);out[name]['probe_ids']=pids
    if len(selected['far_solid_interior_100km'])!=89:
        raise ValueError(f"far-solid-interior cohort has {len(selected['far_solid_interior_100km'])}, expected89")
    depth_groups={}
    for depth,expected in DEPTH_STRATA.items():
        pids=[]
        for pid in selected['far_solid_interior_100km']:
            radius=float(byid[pid]['radius_m']); depth_km=round((6371000.0-radius)/1000.0)
            if depth_km==depth:pids.append(pid)
        if len(pids)!=expected: raise ValueError(f"depth stratum {depth}km has {len(pids)}, expected{expected}")
        depth_groups[str(depth)]=rel_l2([float(ref[p]['numeric_kernel']) for p in pids],
                                       [float(cand[p]['numeric_kernel']) for p in pids],f"depth {depth}")
        depth_groups[str(depth)]['probe_ids']=pids
    out['depth_strata_km']=depth_groups
    out['source_shell_raw_unmatched_trace_semantics']=sum(1 for p in probes if p['layer']=='9' and p['radial_trace'] in ('inner','outer'))
    return out


def relative_kernel_metric(metric: dict[str, Any]) -> dict[str, Any]:
    out=dict(metric)
    refmax=out.get('reference_max_abs',0.0)
    maxdiff=out.get('maximum_absolute_difference',0.0)
    if refmax==0:
        out['max_relative_to_reference_max']=0.0 if maxdiff==0 else None
    else: out['max_relative_to_reference_max']=maxdiff/refmax
    return out


def compare_weights(ref: dict[int,complex], cand: dict[int,complex], label: str) -> dict[str,Any]:
    if set(ref)!=set(cand): raise ValueError(f"{label} bin inventories differ")
    return rel_l2([ref[k] for k in sorted(ref)],[cand[k] for k in sorted(cand)],label)



def read_unit_products(data: Path, weights: dict[int, complex], expected_probes: list[dict[str, str]]) -> dict[tuple[str,int],complex]:
    table=rows(data/'fixed-probe-bin-kernel.csv')
    expected={p['probe_id'] for p in expected_probes}
    if len(table)!=168*35: raise ValueError(f"{data} fixed-probe/bin inventory count differs")
    found={}
    for row in table:
        pid=row['probe_id']; bin_id=int(row['bin']); key=(pid,bin_id)
        if pid not in expected or bin_id not in weights or key in found:
            raise ValueError(f"{data} has invalid or duplicate product key {key}")
        if row['unit_adjoint_computed']=='0':
            if bin_id not in (7,8,41) or weights[bin_id]!=0j:
                raise ValueError(f"{data} skipped a nonzero unit-adjoint bin")
            unit=0j
        elif row['unit_adjoint_computed']=='1':
            unit=complex(finite(row['unit_complex_product_re'],data.name+'.unit.real'),
                         finite(row['unit_complex_product_im'],data.name+'.unit.imag'))
        else: raise ValueError(f"{data} has unknown unit-adjoint marker")
        pval=complex(finite(row['P_re'],data.name+'.P.real'),finite(row['P_im'],data.name+'.P.imag'))
        if pval!=weights[bin_id] and abs(pval-weights[bin_id])>1e-12*max(abs(pval),abs(weights[bin_id])):
            raise ValueError(f"{data} per-bin P differs from its own waveform weights")
        found[key]=unit
    if len(found)!=168*35: raise ValueError(f"{data} product grid is incomplete")
    for bin_id in range(7,42):
        flags={row['unit_adjoint_computed'] for row in table if int(row['bin'])==bin_id}
        expected_flag='0' if bin_id in (7,8,41) else '1'
        if flags!={expected_flag}:
            raise ValueError(f"{data} unit-product compute/skip coverage differs at bin {bin_id}")
    return found


def read_spectrum(data: Path) -> dict[int,complex]:
    result={}
    for row in rows(data/'spectrum.csv'):
        k=int(row['bin'])
        if k in result: raise ValueError(f"{data} repeats spectrum bin {k}")
        result[k]=complex(finite(row['real_m'],data.name+'.spectrum.real'),
                          finite(row['imag_m'],data.name+'.spectrum.imag'))
    if set(result)!=set(range(7,42)): raise ValueError(f"{data} spectrum bins differ from 7..41")
    return result


def spectrum_tail(reference: dict[int,complex], candidate: dict[int,complex]) -> dict[str,Any]:
    per_bin={}
    for k in sorted(reference):
        a,b=reference[k],candidate[k];d=abs(a-b);den=abs(b)
        if den==0 and d!=0: rel=None
        else: rel=d/den if den else 0.0
        per_bin[str(k)]={"lower_cutoff":a.real,"lower_imag":a.imag,
                         "higher_cutoff":b.real,"higher_imag":b.imag,
                         "absolute_difference":d,"relative_to_higher_bin":rel,
                         "higher_bin_exact_zero":den==0}
    all_metric=rel_l2([candidate[k] for k in sorted(candidate)],
                      [reference[k] for k in sorted(reference)],"spectrum tail")
    return {"all_bins":all_metric,"per_bin":per_bin,
            "interpretation":"receiver-spectrum differences are additive cutoff increments over the added degree range"}


def product_tail(reference: dict[tuple[str,int],complex], candidate: dict[tuple[str,int],complex]) -> dict[str,Any]:
    if set(reference)!=set(candidate): raise ValueError("unit product keys differ across cutoffs")
    out={}
    for k in range(7,42):
        if k in (7,8,41):
            out[str(k)]={"status":"unit-adjoint skipped at exact-zero P; not a physical zero measurement"}
            continue
        ids=sorted(pid for pid,bin_id in reference if bin_id==k)
        out[str(k)]=rel_l2([candidate[(pid,k)] for pid in ids],
                           [reference[(pid,k)] for pid in ids],f"unit product bin {k}")
    out['active_unit_bins']=32
    out['skipped_unit_bins']=[7,8,41]
    out['interpretation']="active unit-product differences are complete cutoff differences and include cross terms; skipped exact-zero-P bins are placeholders, not physical zeros, and these are not isolated additive kernel degree blocks"
    return out


def section_difference(reference: Path, candidate: Path, part: str) -> dict[str,Any]:
    ref=(reference/f'{part}.bin').read_bytes(); cand=(candidate/f'{part}.bin').read_bytes()
    if len(ref)!=len(cand) or len(ref)%8: raise ValueError(f"section binary dimensions differ for {part}")
    n=len(ref)//8
    a=struct.unpack(f'{n}d',ref);b=struct.unpack(f'{n}d',cand)
    if any(not math.isfinite(v) for v in (*a,*b)): raise ValueError(f"nonfinite section value in {part}")
    return rel_l2(list(a),list(b),f"{part} section")

def compare_measurement(reference: dict[str,Any], candidate: dict[str,Any]) -> dict[str,Any]:
    if (reference['trace']['times']!=candidate['trace']['times'] or
        reference['trace']['window']!=candidate['trace']['window']):
        raise ValueError("full waveform time grids differ")
    return {"whole_trace_unweighted":rel_l2(reference['trace']['velocity'],candidate['trace']['velocity'],"whole trace"),
            "window_taper_weighted":waveform_metric(reference['trace'],candidate['trace'],True),
            "window_unweighted_99_diagnostic":waveform_metric(reference['trace'],candidate['trace'],False),
            "measurement_rate_norm":compare_scalar(candidate['metadata']['measurement_rate_norm'],reference['metadata']['measurement_rate_norm'],"measurement_rate_norm"),
            "P_weights":compare_weights(reference['weights'],candidate['weights'],"P weights")}


def validate_plan_job(job: dict[str,Any], job_id: str, stage: str, validation: str, build: str) -> None:
    expected={
      'angular-l400-dspec-nq6':('dspec',400,6), 'angular-l400-yspec-nq6':('yspec',400,6),
      'angular-l600-dspec-nq6':('dspec',600,6), 'angular-l600-yspec-nq6':('yspec',600,6),
      'angular-l800-dspec-nq6':('dspec',800,6), 'angular-l800-yspec-nq6':('yspec',800,6),
      'angular-l800-dspec-nq8':('dspec',800,8), 'angular-l800-dspec-nq10':('dspec',800,10),
    }
    backend,degree,nq=expected[job_id]
    meta=job.get('metadata',{})
    if (meta.get('backend'),meta.get('degree'),meta.get('nq'))!=(backend,degree,nq):
        raise ValueError(f"{job_id} metadata differs from its immutable queue slot")
    if job.get('memory_gib')!=16 or job.get('timeout_seconds')!=10800 or job.get('cwd')!=f'{stage}/candidate/dspec':
        raise ValueError(f"{job_id} resource/cwd declaration differs from the approved run")
    cmd=job.get('command',[])
    required=('DSFAST_RADIAL_ROUTE=batched','DSFAST_ANGULAR_ROUTE=real-legendre',
              'DSFAST_RADIAL_ORACLE=0','DSFAST_ANGULAR_ORACLE=0','DSFAST_NATIVE_DERIVATIVE=dref',
              'DSPEC_MODAL_OPERATOR_ROUTE=legacy','DSPEC_MODAL_OPERATOR_ORACLE=0',
              'DSPEC_PHASE_D_ANGULAR_CONVERGENCE_V1=1','DSPEC_PHASE_D_STUDY2_CONTROL=1',
              'DSPEC_PHASE_D_STUDY2_SPARSE_SECTIONS=1','DSPEC_PHASE_C_SKIP_ZERO_KERNEL_BINS=1',
              'DSPEC_PHASE_C_DIRECT_ADJOINT=1','DSPEC_PHASE_C_VALIDATE_DIRECT_ADJOINT=0',
              'DSPEC_PHASE_C_VALIDATE_RECONSTRUCTION=0','DSPEC_PHASE_C2_VALIDATE_WEIGHT_LINE=0',
              'DSPEC_PHASE_D_VALIDATE_BATCHED=0','DSPEC_PHASE_D_VALIDATE_COMMON=0',
              'DSPEC_PHASE_D_SCALAR_ONLY=0','OMP_NUM_THREADS=1','OMP_THREAD_LIMIT=8','OMP_DYNAMIC=FALSE')
    if any(item not in cmd for item in required): raise ValueError(f"{job_id} omits an exact route/validation environment control")
    if cmd.count('DSPEC_PHASE_D_ANGULAR_CONVERGENCE_V1=1')!=1: raise ValueError(f"{job_id} angular-v1 flag is ambiguous")
    if not {'DSPEC_PHASE_D_MESH_SCALE=0.9','DSPEC_PHASE_D_REQUESTED_RELATIVE_ACCURACY=0.001'}.issubset(set(cmd)):
        raise ValueError(f"{job_id} mesh selector/accuracy differs from the frozen case")
    if backend=='dspec':
        if 'DSPEC_FREQUENCY_WORKERS=8' not in cmd or 'DSPEC_PHASE_D_STUDY2_ABSOLUTE_H=0.02035382429229197' not in cmd:
            raise ValueError(f"{job_id} DSpec worker/fixed-h controls differ")
        if any(x.startswith('DSPEC_FREQUENCY_WORKERS=') for x in cmd if x!='DSPEC_FREQUENCY_WORKERS=8'):
            raise ValueError(f"{job_id} contains conflicting worker controls")
        case=f'{stage}/candidate/dspec/case/kernels-small.json'
        refinement='0.02035382429229197';workers='1'
        outputs_d=('data/parameters.txt','data/result.json','data/trace.csv','data/spectrum.csv','data/components.csv',
                   'data/nodes.csv','data/frequency-grid.csv','data/mesh-elements.csv','data/along.bin','data/along-angles.csv',
                   'data/across.bin','data/across-angles.csv','data/fixed-probe-kernel.csv','data/fixed-probe-frequency-weights.csv',
                   'data/fixed-probe-bin-kernel.csv','data/fixed-probe-unit-fields-bin9.csv','data/fixed-probe-unit-fields-bin25.csv',
                   'data/phase-c2-zero-bin-support.csv')
        if set(job.get('outputs',[]))!=set(outputs_d): raise ValueError(f"{job_id} DSpec output contract differs")
    else:
        if any(x.startswith('DSPEC_FREQUENCY_WORKERS=') for x in cmd) or cmd.count('-u')<2 or \
           'DSPEC_FREQUENCY_WORKERS' not in cmd or 'DSPEC_PHASE_D_STUDY2_ABSOLUTE_H' not in cmd:
            raise ValueError(f"{job_id} YSpec must explicitly unset DSpec frequency/fixed-h routes")
        case=f'{validation}/data/kernels-small-yspec-tight-remote.json'
        refinement='0.022615360324768855';workers='8'
        outputs_y=('data/parameters.txt','data/result.json','data/trace.csv','data/spectrum.csv','data/components.csv',
                   'data/nodes.csv','data/frequency-grid.csv','data/mesh-elements.csv','data/along.bin','data/along-angles.csv',
                   'data/across.bin','data/across-angles.csv','data/fixed-probe-kernel.csv','data/fixed-probe-frequency-weights.csv')
        if set(job.get('outputs',[]))!=set(outputs_y): raise ValueError(f"{job_id} YSpec output contract differs")
    i=cmd.index('/usr/bin/taskset')
    expected_argv=[backend,case,'{JOB_DIR}/data',str(degree),str(nq),refinement,workers,'sections','0','-','-',
                   f'{validation}/data/study1-fixed-probes-v2.csv']
    if cmd[i:i+4]!=['/usr/bin/taskset','-c','16-23',f'{build}/phase-d-driver'] or cmd[i+4:]!=expected_argv:
        raise ValueError(f"{job_id} executable argv differs from its exact degree/NQ/case controls")
    if any(('dspec-science-study1-v1' in x or 'dspec-science-study2-v1' in x) for x in job.get('inputs',[])):
        raise ValueError(f"{job_id} has stale candidate input paths")
    ys_dso='/space/adcm2/research/builds/phase-e-stage1-backend-profiling-v1/yspec-timer-disabled/libyspec-doc.so'
    needed={f'{build}/phase-d-driver',f'{build}/command.json',case,
            f'{stage}/candidate/dspec/case/prem.200',f'{validation}/data/study1-fixed-probes-v2.csv',
            '/usr/bin/env','/usr/bin/taskset',ys_dso,ys_dso.rsplit('/',1)[0]+'/command.json'}
    if set(job.get('inputs',[]))!=needed:
        raise ValueError(f"{job_id} input closure differs from the exact build/case/probe/runtime dependencies")


def main() -> None:
    ap=argparse.ArgumentParser()
    for name in ("run-root","run-plan","supervisor-state","probe-csv","probe-receipt","case","tight-case","model",
                 "reference-root","reference-manifest","helpers","out"):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args(); helper=load_helper(args.helpers)
    plan=read_json(args.run_plan); plan_jobs={j['id']:j for j in plan.get('jobs',[])}
    if set(plan_jobs)!=set(JOB_IDS) or len(plan_jobs)!=8: raise ValueError("run plan differs from exact 8-job queue")
    if plan.get('thread_cap')!=20: raise ValueError("run plan omits the reviewed 20-thread supervisor cap")
    base='/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya'
    stage=f'{base}/scientific-inputs/dspec-science-angular-v1'
    validation=f'{base}/scientific-inputs/dspec-science-angular-v1-validation'
    build='/space/adcm2/research/builds/dspec-science-angular-v1/dspec-timer-disabled'
    if plan.get('campaign_root')!=f'{base}/dspec-science-angular-v1-control':
        raise ValueError("campaign root is not the approved angular-v1 campaign service root")
    if plan.get('environment',{}).get('DSPEC_PHASE_D_ANGULAR_CONVERGENCE_V1')!='1':
        raise ValueError("plan does not pin angular-v1 mode")
    if plan.get('reference_manifest_sha256')!=sha(args.reference_manifest):
        raise ValueError("run plan does not bind the exact reference manifest")
    build_plan=f'{stage}/configs/angular-build-v1-plan.json'
    external=plan.get('external_campaigns',[])
    if len(external)!=1 or external[0].get('plan_path')!=build_plan or \
       external[0].get('plan_sha256')!='44131d560ce27aa42ff4c086a6df64e1ecd5d810056bd06c47bb59d827de258d':
        raise ValueError("producer plan is not bound to the reviewed angular build plan")
    for jid in JOB_IDS: validate_plan_job(plan_jobs[jid],jid,stage,validation,build)
    state_jobs=verify_state(args.run_root,args.run_plan,args.supervisor_state)
    ref_integrity=verify_manifest(args.reference_root,args.reference_manifest)
    probes=helper.rows(args.probe_csv); receipt=read_json(args.probe_receipt)
    if len(probes)!=168 or sha(args.probe_csv)!=receipt.get('probe_csv_sha256'):
        raise ValueError("fixed physical probe grid differs from its exact 168-point receipt")
    if sha(args.case)!=receipt.get('case_sha256') or sha(args.model)!="4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea":
        raise ValueError("case or PREM model hash differs from the frozen physical fixture")
    data_root=args.reference_root.parent/'data'
    if (sha(data_root/'generate_study1_fixed_probes.py')!=receipt.get('generator_sha256') or
        sha(data_root/'study1-phase-d-along-angles.csv')!=receipt.get('path_reference_sha256') or
        sha(data_root/'study1-phase-d-mesh-elements.csv')!=receipt.get('mesh_receipt_sha256')):
        raise ValueError("probe generator/path/mesh provenance differs from the fixed-probe receipt")
    case=read_json(args.case); tight=read_json(args.tight_case)
    expected_tight=json.loads(json.dumps(case)); expected_tight.pop('solver',None)
    tight_cmp=json.loads(json.dumps(tight));
    if tight_cmp.get('solver')!={"tolerance":1e-7}: raise ValueError("YSpec case does not pin ODE tolerance 1e-7")
    tight_cmp.pop('solver'); tight_cmp['model']=expected_tight['model']
    if tight_cmp!=expected_tight or tight.get('model')!=str(args.model):
        raise ValueError("YSpec case differs from the frozen case outside model path/ODE tolerance")
    source=receipt['source']; receiver=receipt['receiver']
    source_xyz=helper.cart(float(source['radius_m']),float(source['geocentric_theta_rad']),float(source['geocentric_phi_rad']))
    receiver_xyz=helper.cart(float(receiver['radius_m']),float(receiver['geocentric_theta_rad']),float(receiver['geocentric_phi_rad']))
    result={}; by_case={}
    for jid in JOB_IDS:
        cfg=plan_jobs[jid]['metadata']
        if cfg.get('purpose') not in ('angular degree ladder','L800 radial-order control'):
            raise ValueError(f"{jid} purpose metadata not approved")
        backend=cfg['backend']; degree=cfg['degree']; nq=cfg['nq']
        job={"id":jid,"backend":backend,"degree":degree,"nq":nq}
        checked=check_result(args.run_root/jid/'data',job,helper,probes,source_xyz,receiver_xyz)
        checked['data_path']=args.run_root/jid/'data'
        checked['completion']={"returncode":state_jobs[jid]['returncode'],"inputs_unchanged":state_jobs[jid]['inputs_unchanged'],
                              "cleanup_verified_gone":state_jobs[jid]['cleanup']['verified_gone']}
        result[jid]=checked
        if nq==6: by_case[(backend,degree)]=checked
    # Validate exact spectral/time/taper/section grids for every job, including NQ8/NQ10.
    reference_dir={"dspec":args.reference_root/'l300-dspec',"yspec":args.reference_root/'l300-yspec-tight'}
    baseline_traces={backend:check_trace(reference_dir[backend],backend+' L300 reference')
                     for backend in ('dspec','yspec')}
    d800=by_case[('dspec',800)]
    for backend in ('dspec','yspec'):
        ref=reference_dir[backend]
        base_grid=(ref/'frequency-grid.csv').read_bytes()
        ref_angles={part:rows(ref/f'{part}-angles.csv') for part in ('along','across')}
        ref_trace=baseline_traces[backend]
        for degree in (400,600,800):
            run=by_case[(backend,degree)]
            data=args.run_root/f'angular-l{degree}-{backend}-nq6'/'data'
            if (data/'frequency-grid.csv').read_bytes()!=base_grid:
                raise ValueError(f"{backend} L{degree} frequency grid differs from L300")
            if run['angles']!=ref_angles:
                raise ValueError(f"{backend} L{degree} section geometry differs from L300")
            if (data/'nodes.csv').read_bytes()!=(ref/'nodes.csv').read_bytes():
                raise ValueError(f"{backend} L{degree} NQ6 radial output nodes differ from L300")
            if (data/'mesh-elements.csv').read_bytes()!=(ref/'mesh-elements.csv').read_bytes():
                raise ValueError(f"{backend} L{degree} NQ6 element boundaries differ from L300")
            if run['trace']['times']!=ref_trace['times'] or run['trace']['window']!=ref_trace['window']:
                raise ValueError(f"{backend} L{degree} trace time/taper grid differs from L300")
    for nq in (8,10):
        jid=f'angular-l800-dspec-nq{nq}'; run=result[jid]; data=args.run_root/jid/'data'
        if ((data/'frequency-grid.csv').read_bytes()!=(reference_dir['dspec']/'frequency-grid.csv').read_bytes() or
            run['trace']['times']!=baseline_traces['dspec']['times'] or
            run['trace']['window']!=baseline_traces['dspec']['window']):
            raise ValueError(f"{jid} frequency/time/taper grid differs from the frozen fixture")
        if run['angles']!=d800['angles']:
            raise ValueError(f"{jid} section angles differ from L800 NQ6")
    # Verify high-L radial controls preserve the exact element boundaries and fixed angular grid.
    for nq in (8,10):
        j=f'angular-l800-dspec-nq{nq}'; r=result[j]
        if r['boundaries']!=d800['boundaries']:
            raise ValueError(f"L800 NQ{nq} changes element boundaries")
        if r['angles']!=d800['angles']:
            raise ValueError(f"L800 NQ{nq} changes section grid")
    # Load saved Study2 L200/L300 references and compare cutoff increments and candidates to L800.
    refs={}
    for backend in ('dspec','yspec'):
        refname='l200-dspec-sparse-bridge' if backend=='dspec' else 'l200-yspec-tight'
        refdata=args.reference_root/refname
        refs[(backend,200)]={"metadata":read_json(refdata/'result.json'),"trace":check_trace(refdata,backend+' L200'),
                            "weights":helper.read_weights(refdata),"integrated":helper.fixed_table(refdata,probes,backend,source_xyz,receiver_xyz),
                            "unit_products":read_unit_products(refdata,helper.read_weights(refdata),probes) if backend=='dspec' else None,
                            "data_path":refdata}
        refs[(backend,300)]=by_case[(backend,300)] if (backend,300) in by_case else None
        # This campaign does not rerun L300; use the hash-bound Study2 output.
        d300=args.reference_root/f'l300-{backend if backend=="dspec" else "yspec-tight"}'
        refs[(backend,300)]={"metadata":read_json(d300/'result.json'),"trace":check_trace(d300,backend+' L300'),
                            "weights":helper.read_weights(d300),"integrated":helper.fixed_table(d300,probes,backend,source_xyz,receiver_xyz),
                            "unit_products":read_unit_products(d300,helper.read_weights(d300),probes) if backend=='dspec' else None,
                            "data_path":d300}
    comparisons={}
    levels=(200,300,400,600,800)
    for backend in ('dspec','yspec'):
        comps={}
        for lo,hi in zip(levels,levels[1:]):
            a=refs[(backend,lo)] if lo<400 else by_case[(backend,lo)]
            b=refs[(backend,hi)] if hi<400 else by_case[(backend,hi)]
            key=f'{lo}_to_{hi}'
            # The higher degree is the fixed reference for each increment.
            kernel=compare_tables(b['integrated'],a['integrated'],probes)
            measurement=compare_measurement(b,a)
            comps[key]={"kernel":kernel,"measurement":measurement}
        # Every rung versus the highest available L800 reference.
        top=by_case[(backend,800)]
        for lo in (200,300,400,600):
            a=refs[(backend,lo)] if lo<400 else by_case[(backend,lo)]
            comps[f'{lo}_vs_800']={"kernel":compare_tables(top['integrated'],a['integrated'],probes),
                                  "measurement":compare_measurement(top,a)}
        comparisons[backend]=comps
    tails={}; section_diagnostics={}
    for backend in ('dspec','yspec'):
        paths={300:refs[(backend,300)]['data_path']}
        paths.update({level:by_case[(backend,level)]['data_path'] for level in (400,600,800)})
        tails[backend]={}; section_diagnostics[backend]={}
        for lo,hi in zip((300,400,600),(400,600,800)):
            lo_item=refs[(backend,lo)] if lo==300 else by_case[(backend,lo)]
            hi_item=by_case[(backend,hi)]
            key=f'{lo}_to_{hi}'
            tails[backend][key]={"receiver_spectrum":spectrum_tail(read_spectrum(paths[lo]),read_spectrum(paths[hi]))}
            if backend=='dspec':
                tails[backend][key]['fixed_probe_unit_products']=product_tail(lo_item['unit_products'],hi_item['unit_products'])
            section_diagnostics[backend][key]={part:section_difference(paths[lo],paths[hi],part) for part in ('along','across')}
        for lo in (300,400,600):
            lo_item=refs[(backend,lo)] if lo==300 else by_case[(backend,lo)]
            hi_item=by_case[(backend,800)]
            key=f'{lo}_vs_800'
            tails[backend][key]={"receiver_spectrum":spectrum_tail(read_spectrum(paths[lo]),read_spectrum(paths[800]))}
            if backend=='dspec':
                tails[backend][key]['fixed_probe_unit_products']=product_tail(lo_item['unit_products'],hi_item['unit_products'])
            section_diagnostics[backend][key]={part:section_difference(paths[lo],paths[800],part) for part in ('along','across')}
    # DSpec only: fixed P800 counterfactual from the saved unit complex products; keep own-P as primary.
    dspec_refs={}
    for level in (200,300): dspec_refs[level]=refs[('dspec',level)]
    dspec_refs.update({level:by_case[('dspec',level)] for level in (400,600,800)})
    p800=dspec_refs[800]['weights']
    for level in (200,300,400,600,800):
        item=dspec_refs[level]
        if item.get('unit_products') is None:
            data=(args.reference_root/('l200-dspec' if level==200 else 'l300-dspec')) if level<400 else args.run_root/f'angular-l{level}-dspec-nq6'/'data'
            w=helper.read_weights(data)
            units=read_unit_products(data,w,probes)
            item['unit_products']=units
        own={pid:sum(item['weights'][k]*item['unit_products'][(pid,k)] for k in range(7,42)).real for pid in {key[0] for key in item['unit_products']}}
        counter={pid:sum(p800[k]*item['unit_products'][(pid,k)] for k in range(7,42)).real for pid in own}
        dspec_refs[level]['own_P_counterfactual']={"P800_weighted_full_grid":counter,
                                                   "P_own_weighted_full_grid":own,
                                                   "label":"DSpec unit products combined with L800 P; counterfactual, not native kernel"}
    cross_method={}
    for degree in (400,600,800):
        d=by_case[('dspec',degree)]; y=by_case[('yspec',degree)]
        cross_method[str(degree)]={"reference":"DSpec","candidate":"YSpec",
             "fixed_probe_tables":compare_tables(d['integrated'],y['integrated'],probes),
             "measurement":compare_measurement(d,y),
             "weight_note":"each backend uses its own waveform-derived P; diagnostic only and does not isolate field differences"}
    # Required acceptance gates: 0.1% on each named measure/cohort, no floors.
    acceptance={}
    for backend,comps in comparisons.items():
        accepted={}
        for level in (400,600):
            gate_names=(f'{level}_vs_800',f'{level}_to_600',f'600_to_800') if level==400 else ('600_vs_800','600_to_800')
            checks=[]
            for name in gate_names:
                item=comps[name]
                measures=[item['measurement']['window_taper_weighted'],item['measurement']['P_weights']]
                norm=item['measurement']['measurement_rate_norm']['relative_difference']
                checks.extend(ratio_pass(m,TARGET,backend+' '+name) for m in measures)
                checks.append(norm is not None and norm<=TARGET)
                for cohort in ('far_solid_interior_100km',):
                    checks.append(ratio_pass(item['kernel'][cohort],TARGET,backend+' '+name+' '+cohort))
                    maxratio=item['kernel'][cohort]['max_relative_to_reference_max']
                    checks.append(maxratio is not None and maxratio<=TARGET)
                    for depth,m in item['kernel']['depth_strata_km'].items():
                        checks.append(ratio_pass(m,TARGET,backend+' '+name+' depth'+depth))
                        maxratio=m['max_relative_to_reference_max']
                        checks.append(maxratio is not None and maxratio<=TARGET)
            certified=(level==400 and all(checks))
            accepted[str(level)]={"passes":all(checks),"certified":certified,
                                  "successive_increment_count":2 if level==400 else 1,
                                  "required_comparison_sets":gate_names,
                                  "criterion":"every metric in each required comparison <=0.001, no denominator floors"}
        acceptance[backend]=accepted
    radial={}; radial_control_pass={}
    nq_results={nq:result[f'angular-l800-dspec-nq{nq}'] for nq in (6,8,10)}
    for low,high in ((6,8),(6,10),(8,10)):
        ref=nq_results[high]; cand=nq_results[low]
        label=f'nq{low}_vs_nq{high}'
        metric={"reference_nq":high,"candidate_nq":low,
                "kernel":compare_tables(ref['integrated'],cand['integrated'],probes),
                "measurement":compare_measurement(ref,cand)}
        radial[label]=metric
        checks=[ratio_pass(metric['measurement'][name],TARGET,label)
                for name in ('window_taper_weighted','P_weights')]
        norm=metric['measurement']['measurement_rate_norm']['relative_difference']
        checks.append(norm is not None and math.isfinite(norm) and 0<=norm<=TARGET)
        kernel=metric['kernel']['far_solid_interior_100km']
        checks.extend((ratio_pass(kernel,TARGET,label+' far-solid'),
                       kernel['max_relative_to_reference_max'] is not None and
                       kernel['max_relative_to_reference_max']<=TARGET))
        for depth,m in metric['kernel']['depth_strata_km'].items():
            checks.extend((ratio_pass(m,TARGET,label+' depth'+depth),
                           m['max_relative_to_reference_max'] is not None and
                           m['max_relative_to_reference_max']<=TARGET))
        radial_control_pass[label]=all(checks)
    radial_reference_validated=all(radial_control_pass.values())
    for backend in acceptance:
        for entry in acceptance[backend].values():
            entry['angular_metrics_pass']=entry['passes']
            entry['dspec_radial_controls_pass']=radial_reference_validated
            if backend=='dspec':
                entry['passes']=entry['passes'] and radial_reference_validated
                entry['certified']=entry['certified'] and radial_reference_validated
                entry['certification_note']=('DSpec certification also requires all L800 NQ6/8/10 radial comparisons to pass'
                                             if not radial_reference_validated else
                                             'all L800 NQ6/8/10 radial comparisons passed')
            else:
                entry['certification_note']='YSpec angular criteria are reported separately; DSpec radial controls do not alter this backend result'
    report={"schema":"angular-v1-result-v1","plan_sha256":sha(args.run_plan),"state_sha256":sha(args.supervisor_state),
            "reference_integrity":ref_integrity,"job_ids":list(JOB_IDS),"target_fraction":TARGET,
            "case_scope":{"filter_mHz":FILTER_MHZ,"bins":[7,41],"window_s":list(WINDOW),"taper_s":40,
                          "measurement":"transverse velocity traveltime; each degree uses own waveform P"},
            "comparisons":comparisons,"cutoff_tails":tails,"section_diagnostics":section_diagnostics,
            "cross_method_diagnostics":cross_method,
            "acceptance":acceptance,"radial_controls":radial,"radial_control_pass":radial_control_pass,
            "radial_reference_validated":radial_reference_validated,
            "jobs":{jid:{"metadata":result[jid]['metadata'],"memory":result[jid]['memory'],
                          "completion":result[jid]['completion']} for jid in JOB_IDS},
            "dspec_P800_counterfactual":{str(level):dspec_refs[level]['own_P_counterfactual'] for level in (200,300,400,600,800)},
            "unit_product_coverage":{"active_bins":32,"skipped_zero_adjoint_bins":[7,8,41],
                                     "skipped_values_are_placeholders_not_physical_zero_evidence":True},
            "limitations":["source-shell/contact values are retained diagnostics, not matched one-sided bounded-pointwise gates",
                           "the high-L independent full three-dimensional derivative action is unsupported by the current API and was not claimed",
                           "degree-block receiver spectrum increments are additive, but full unit-product/kernel differences include cross terms"]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({"out":str(args.out),"acceptance":acceptance,"radial_controls":{
        k:{"far_kernel_l2":v['kernel']['far_solid_interior_100km']['relative_l2'],
           "window_weighted":v['measurement']['window_taper_weighted']['relative_l2']} for k,v in radial.items()}},sort_keys=True))


if __name__ == "__main__":
    main()
