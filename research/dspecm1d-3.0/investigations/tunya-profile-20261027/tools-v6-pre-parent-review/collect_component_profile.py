#!/usr/bin/env python3
"""Validate and summarize a hash-verified profile-06 evidence snapshot."""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


EXPECTED_IDS = ["qualification-checks"] + [
    f"{prefix}-l{l}" for l in (128, 384, 1536)
    for prefix in ("dsv4-native", "dsv4-profiled", "dsv6-native", "dsv6-profiled",
                   "baseline", "ysv6-native-wave", "ysv6-native", "ysv6-profiled-wave", "ysv6-profiled")
]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def finite_tree(value, where="root"):
    if isinstance(value, dict):
        for k, v in value.items():
            finite_tree(v, f"{where}.{k}")
    elif isinstance(value, list):
        for i, v in enumerate(value):
            finite_tree(v, f"{where}[{i}]")
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        if not math.isfinite(value):
            raise ValueError(f"non-finite number at {where}")


def validate_snapshot(root: Path):
    state = read_json(root / "state.json")
    receipt = read_json(root / "snapshot-receipt.json")
    listed = receipt.get("sha256")
    if not isinstance(listed, dict):
        raise ValueError("snapshot receipt has no sha256 map")
    if "state.json" in listed and listed["state.json"] != sha(root / "state.json"):
        raise ValueError("snapshot state hash mismatch")
    for rel, expected in listed.items():
        p = root / rel
        if not p.is_file() or sha(p) != expected:
            raise ValueError(f"snapshot hash mismatch: {rel}")
    jobs = {j["id"]: j for j in state.get("jobs", [])}
    for jid, job in jobs.items():
        if job.get("status") != "completed":
            continue
        if job.get("returncode") != 0 or job.get("inputs_unchanged") is not True:
            raise ValueError(f"completed job lacks successful/invariant receipt: {jid}")
        if job.get("cleanup", {}).get("verified_gone") is not True:
            raise ValueError(f"completed job cleanup not verified: {jid}")
        logrel = f"{jid}/job.log"
        if logrel not in listed or listed[logrel] != job.get("log_sha256"):
            raise ValueError(f"snapshot missing or state log hash disagreement: {jid}")
        for rel, digest in job.get("outputs_sha256", {}).items():
            key = f"{jid}/{rel}"
            if key not in listed or listed[key] != digest:
                raise ValueError(f"snapshot missing or state output hash disagreement: {key}")
        if set(job.get("outputs_sha256", {})) != set(job.get("outputs", [])):
            raise ValueError(f"completed job output manifest is incomplete: {jid}")
        for rel in job.get("outputs", []):
            p = root / jid / rel
            if p.suffix == ".json":
                finite_tree(read_json(p), f"{jid}/{rel}")
            elif p.suffix == ".csv":
                raw = np.genfromtxt(p, delimiter=",", names=True, dtype=float, encoding=None)
                if raw.size == 0 or not np.isfinite(np.column_stack([raw[n] for n in raw.dtype.names])).all():
                    raise ValueError(f"empty or non-finite numeric CSV: {jid}/{rel}")
            elif p.suffix == ".bin":
                # These frozen probe artifacts are native float64 arrays; stream via mmap.
                if p.stat().st_size % 8:
                    raise ValueError(f"binary output is not a whole number of float64 values: {jid}/{rel}")
                data = np.memmap(p, dtype=np.float64, mode="r")
                if not np.isfinite(data).all():
                    raise ValueError(f"non-finite binary output: {jid}/{rel}")
                del data
    return state, jobs, listed


def file_pair(root, jobs, listed, a, b, paths):
    ja, jb = jobs.get(a), jobs.get(b)
    if not ja or not jb or ja.get("status") != "completed" or jb.get("status") != "completed":
        return {"status": "pending", "jobs": [a, b]}
    comparisons = []
    for rel in paths:
        ka, kb = f"{a}/{rel}", f"{b}/{rel}"
        if ka not in listed or kb not in listed:
            raise ValueError(f"required paired file absent from verified snapshot: {ka} or {kb}")
        comparisons.append({"file": rel, "byte_exact": listed[ka] == listed[kb],
                            "sha256_native": listed[ka], "sha256_profiled": listed[kb]})
    return {"status": "PASS" if all(x["byte_exact"] for x in comparisons) else "FAIL",
            "jobs": [a, b], "files": comparisons}


def profile_from_log(path: Path):
    found = []
    for line in path.read_text(errors="replace").splitlines():
        marker = "MATHPROFILE "
        if marker in line:
            found.append(json.loads(line.split(marker, 1)[1]))
    if len(found) != 1:
        raise ValueError(f"expected one MATHPROFILE record in {path}, found {len(found)}")
    p = found[0]
    stages = p.get("stages")
    if not isinstance(stages, dict) or not stages:
        raise ValueError(f"profile stages missing in {path}")
    finite_tree(p, str(path))
    roots = [v for k, v in stages.items() if k in ("physical_case_total", "physical_kernel_total", "physical_kernel_probe_total")]
    if not roots:
        roots = [v for k, v in stages.items() if k.endswith("_total")]
    if len(roots) != 1:
        raise ValueError(f"cannot identify unique profile root in {path}")
    exclusive = sum(float(v["exclusive_s"]) for v in stages.values())
    inclusive = float(roots[0]["inclusive_s"])
    return p, {"exclusive_sum_s": exclusive, "root_inclusive_s": inclusive,
               "difference_s": exclusive - inclusive,
               "closes_within_1e-6_s": abs(exclusive - inclusive) <= max(1e-6, 2e-12 * abs(inclusive))}


def trace_array(path):
    a = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding=None)
    if a.size == 0 or not a.dtype.names:
        raise ValueError(f"empty/malformed trace: {path}")
    arr = np.column_stack([a[n] for n in a.dtype.names])
    if not np.isfinite(arr).all():
        raise ValueError(f"non-finite trace: {path}")
    return a.dtype.names, arr


def waveform_metrics(dpath, ypath):
    _, d = trace_array(dpath)
    _, y = trace_array(ypath)
    n = min(len(d), len(y))
    if n == 0 or not np.array_equal(d[:n, 0], y[:n, 0]):
        raise ValueError("DSpec/YSpec waveform traces do not share an identical prefix grid")
    t, x, w = y[:n, 0], d[:n, 1], y[:n, 2]
    _, yv, yw = y[:n, 0], y[:n, 1], y[:n, 2]
    if not np.array_equal(w, yw):
        raise ValueError("window weights differ")
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    mask = (t >= 1000.0) & (t <= 1200.0)
    # Existing campaign energy definition: 0.5*dt*sum(Hann*vertical^2).
    ed = 0.5 * dt * float(np.sum(w[mask] * x[mask] ** 2))
    ey = 0.5 * dt * float(np.sum(w[mask] * yv[mask] ** 2))
    diff = x[mask] - yv[mask]
    weighted_l2 = math.sqrt(float(np.sum(w[mask] * diff ** 2)) * dt)
    relative_energy = abs(ed - ey) / max(abs(ed), np.finfo(float).tiny)
    return {"window_s": [1000, 1200], "energy_definition": "0.5*dt*sum(Hann1000-1200*vertical_m^2)",
            "dSpec_energy_m2_s": ed, "ySpec_energy_m2_s": ey,
            "relative_energy_difference": relative_energy, "weighted_trace_l2_m_sqrt_s": weighted_l2,
            "compared_samples": n, "dSpec_samples": len(d), "ySpec_samples": len(y),
            "interpretation": "descriptive cross-code waveform comparison on the shared prefix; existing physical energy scale is 1e-3; no new pass threshold applied"}


def cross_code_channels(root, jobs, listed, l):
    dsid, ysid = f"dsv6-native-l{l}", f"ysv6-native-l{l}"
    if jobs.get(dsid, {}).get("status") != "completed" or jobs.get(ysid, {}).get("status") != "completed":
        return {"status": "pending", "jobs": [dsid, ysid]}
    dmeta = read_json(root / dsid / "data/kernel_probe.json")
    ymeta = read_json(root / ysid / "data/result.json")
    expected_ds = "point-major; actual sampled angle_index outer, selected radial side inner, doubles[A,C,F,L,N,kappa] summed over frequency"
    if dmeta.get("binary_layout") != expected_ds:
        raise ValueError(f"unexpected DSpec binary layout for L{l}")
    if ymeta.get("kernel_binary_layout") != "float64 native endian; angle-major, then radial node, then A,C,F,L,N parameter; v6 layout":
        raise ValueError(f"unexpected YSpec binary layout for L{l}")
    points = int(dmeta.get("actual_sampled_point_count", -1))
    if points != int(ymeta.get("kernel_actual_points", -2)) or points <= 0:
        raise ValueError(f"cross-code sample count mismatch for L{l}")
    dp = root / dsid / "data/kernel_probe.bin"
    yp = root / ysid / "data/kernel_probe.bin"
    if dp.stat().st_size != points * 6 * 8 or yp.stat().st_size != points * 5 * 8:
        raise ValueError(f"cross-code binary shape mismatch for L{l}")
    d = np.memmap(dp, dtype=np.float64, mode="r").reshape(points, 6)
    y = np.memmap(yp, dtype=np.float64, mode="r").reshape(points, 5)
    names = ("A", "C", "F", "L", "N")
    out = {}
    for i, name in enumerate(names):
        delta = np.abs(d[:, i] - y[:, i])
        loc = int(np.argmax(delta))
        scale = max(float(np.max(np.abs(d[:, i]))), float(np.max(np.abs(y[:, i]))), np.finfo(float).tiny)
        out[name] = {"max_abs_difference": float(delta[loc]), "scale_max_abs": scale,
                     "scaled_max_abs_difference": float(delta[loc]) / scale,
                     "max_location_point_index": loc, "dSpec_value": float(d[loc, i]), "ySpec_value": float(y[loc, i])}
    del d, y
    return {"status": "descriptive_only", "L": l, "points": points,
        "channels": out,
        "limitation": "first-64-angle kernel diagnostic only; cross-code differences do not establish full-grid accuracy"}


def collect(root: Path):
    state, jobs, listed = validate_snapshot(root)
    pairs = {}
    for l in (128, 384, 1536):
        for ver in ("v4", "v6"):
            pairs[f"ds-{ver}-kernel-native-profile-L{l}"] = file_pair(root, jobs, listed,
                f"ds{ver}-native-l{l}", f"ds{ver}-profiled-l{l}",
                ["data/kernel_probe.bin", "data/trace.csv", "data/spectrum.csv"])
        pairs[f"ys-kernel-native-profile-L{l}"] = file_pair(root, jobs, listed,
            f"ysv6-native-l{l}", f"ysv6-profiled-l{l}", ["data/kernel_probe.bin"])
        pairs[f"ys-wave-native-profile-L{l}"] = file_pair(root, jobs, listed,
            f"ysv6-native-wave-l{l}", f"ysv6-profiled-wave-l{l}", ["data/trace.csv", "data/spectrum.csv"])
    profiles = {}
    for jid, job in jobs.items():
        if "profiled" not in jid or job.get("status") != "completed":
            continue
        p, closure = profile_from_log(root / jid / "job.log")
        if jid.startswith("ds"):
            expected_solver, expected_objective = "dspec", -1
        elif "-wave-" in jid:
            expected_solver, expected_objective = "YSpecpp.physical_waveform", 0
        else:
            expected_solver, expected_objective = "YSpecpp.physical_kernel-probe", 0
        if p.get("solver") != expected_solver or p.get("objective") != expected_objective:
            raise ValueError(f"unexpected profiler identity/objective for {jid}: {p.get('solver')}/{p.get('objective')}")
        profiles[jid] = {"solver": p.get("solver"), "objective": p.get("objective"),
                         "stages": p["stages"], "counters": p.get("counters", {}), "root_closure": closure}
    telemetry = {}
    for jid, job in jobs.items():
        if jid == "qualification-checks":
            continue
        telemetry[jid] = {"status": job.get("status", "missing"), "returncode": job.get("returncode"),
            "wall_elapsed_s": job.get("elapsed_seconds"), "child_user_cpu_s": job.get("user_cpu_seconds_children"),
            "child_system_cpu_s": job.get("system_cpu_seconds_children"),
            "peak_observed_rss_gib_userspace": job.get("peak_observed_rss_gib_userspace"),
            "peak_observed_threads": job.get("peak_observed_threads"),
            "cleanup_verified": job.get("cleanup", {}).get("verified")}
    waveforms = {}
    cross_code = {}
    for l in (128, 384, 1536):
        cross_code[str(l)] = cross_code_channels(root, jobs, listed, l)
        required = [f"dsv6-native-l{l}", f"ysv6-native-wave-l{l}"]
        if all(jobs.get(x, {}).get("status") == "completed" for x in required):
            waveforms[str(l)] = waveform_metrics(root / required[0] / "data/trace.csv",
                                                 root / required[1] / "data/trace.csv")
        else:
            waveforms[str(l)] = {"status": "pending", "jobs": required}
    pending = [jid for jid in EXPECTED_IDS if jobs.get(jid, {}).get("status") != "completed"]
    overhead = {}
    for l in (128, 384, 1536):
        for native, profiled, label in ((f"dsv6-native-l{l}", f"dsv6-profiled-l{l}", "dsv6-kernel"),
                                        (f"ysv6-native-l{l}", f"ysv6-profiled-l{l}", "ysv6-kernel"),
                                        (f"ysv6-native-wave-l{l}", f"ysv6-profiled-wave-l{l}", "ysv6-waveform")):
            a, b = jobs.get(native), jobs.get(profiled)
            key = f"{label}-L{l}"
            if not a or not b or a.get("status") != "completed" or b.get("status") != "completed":
                overhead[key] = {"status": "pending"}
            else:
                def delta(field):
                    av, bv = a.get(field), b.get(field)
                    return None if av is None or bv is None else bv - av
                overhead[key] = {"status": "apparent single-pass difference (profiled minus native)",
                    "wall_elapsed_s": delta("elapsed_seconds"),
                    "child_user_cpu_s": delta("user_cpu_seconds_children"),
                    "child_system_cpu_s": delta("system_cpu_seconds_children")}
    return {"schema": "component-profile-summary-v1", "snapshot_root": str(root.resolve()),
        "campaign_status_at_snapshot": state.get("status"), "task_id": state.get("task_id"),
        "thread_cap": state.get("thread_cap"), "claim_limits": [
            "Kernel probes use only the first 64 angular points for each L and are diagnostics, not full global kernels.",
            "L128/L384/L1536 cases are not a degree-convergence or representative global-throughput claim.",
            "The physical spatial grid is 64*179=11456 points; tensor_volume_point_count=595712 counts 52 frequency visits.",
            "Counters such as modal_lmax and solver_degree are accumulated counts; actual degree is in result metadata.",
            "Native/profiled wall-time differences are single-pass apparent measurements, not statistical estimates.",
            "MATHPROFILE inclusive scopes overlap; only exclusive stage sums are used for root-closure accounting."],
        "pending_jobs": pending, "telemetry": telemetry, "exact_parity": pairs,
        "profiled_stage_records": profiles, "native_profile_overhead_apparent": overhead,
        "waveform_cross_code_descriptive": waveforms, "kernel_cross_code_descriptive": cross_code,
        "snapshot_file_count": len(listed)}


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: collect_component_profile.py SNAPSHOT_DIR OUTPUT_JSON")
    root, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    if not root.is_dir() or out.exists():
        raise SystemExit("snapshot must exist and output must be fresh")
    result = collect(root)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x") as f:
        json.dump(result, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
