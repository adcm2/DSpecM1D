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
    f"{prefix}-l{l}" for l in (128, 384)
    for prefix in ("dsv4-native", "dsv4-profiled", "dsv6-native", "dsv6-profiled",
                   "baseline", "ysv6-native-wave", "ysv6-native", "ysv6-profiled-wave", "ysv6-profiled")
] + [f"{prefix}-l1536" for prefix in
     ("dsv6-native", "dsv6-profiled", "baseline", "ysv6-native-wave", "ysv6-native",
      "ysv6-profiled-wave", "ysv6-profiled")]


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
    if "state.json" not in listed or listed["state.json"] != sha(root / "state.json"):
        raise ValueError("snapshot state hash mismatch")
    for rel, expected in listed.items():
        p = root / rel
        if not p.is_file() or sha(p) != expected:
            raise ValueError(f"snapshot hash mismatch: {rel}")
    jobs = {j["id"]: j for j in state.get("jobs", [])}
    if len(jobs) != len(state.get("jobs", [])):
        raise ValueError("duplicate job ids in campaign state")
    unexpected = sorted(set(jobs) - set(EXPECTED_IDS))
    if unexpected:
        raise ValueError(f"unexpected job ids in campaign state: {unexpected}")
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
    if set(k for k in stages if k == "physical_case_total") != {"physical_case_total"}:
        raise ValueError(f"physical_case_total root missing or duplicated in {path}")
    root_stage = stages["physical_case_total"]
    for name, stage in stages.items():
        for field in ("inclusive_s", "exclusive_s", "calls"):
            if field not in stage or not math.isfinite(float(stage[field])) or float(stage[field]) < 0:
                raise ValueError(f"invalid nonnegative {field} in {name} at {path}")
        if int(stage["calls"]) != float(stage["calls"]):
            raise ValueError(f"noninteger calls in {name} at {path}")
    exclusive = sum(float(v["exclusive_s"]) for v in stages.values())
    inclusive = float(root_stage["inclusive_s"])
    tolerance = max(1e-8, 2e-12 * abs(inclusive))
    return p, {"exclusive_sum_s": exclusive, "root_inclusive_s": inclusive,
               "difference_s": exclusive - inclusive,
               "closure_tolerance_s": tolerance,
               "closes_within_existing_tolerance": abs(exclusive - inclusive) <= tolerance}


def trace_array(path):
    a = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding=None)
    if a.size == 0 or not a.dtype.names:
        raise ValueError(f"empty/malformed trace: {path}")
    arr = np.column_stack([a[n] for n in a.dtype.names])
    if not np.isfinite(arr).all():
        raise ValueError(f"non-finite trace: {path}")
    return a.dtype.names, arr


def waveform_metrics(dpath, ypath, dresult, yresult):
    _, d = trace_array(dpath)
    _, y = trace_array(ypath)
    n = min(len(d), len(y))
    if n == 0 or not np.array_equal(d[:n, 0], y[:n, 0]):
        raise ValueError("DSpec/YSpec waveform traces do not share an identical prefix grid")
    if not np.all(np.diff(d[:, 0]) == 5.0) or not np.all(np.diff(y[:, 0]) == 5.0):
        raise ValueError("waveform traces must use the exact 5-second sample spacing")
    if len(d) != 1024 or len(y) != 721:
        raise ValueError(f"unexpected waveform row counts: DSpec={len(d)} YSpec={len(y)}")
    t, x, dw = y[:, 0], d[:n, 1], d[:n, 2]
    _, yv, yw = y[:, 0], y[:, 1], y[:, 2]
    if float(np.max(np.abs(dw - yw))) > 1e-12:
        raise ValueError("DSpec/YSpec Hann weights differ by more than 1e-12")
    if np.any(d[721:, 2] != 0.0):
        raise ValueError("DSpec Hann weights are nonzero beyond the 721-sample YSpec prefix")
    expected = np.where((t >= 1000.0) & (t <= 1200.0),
                        0.5 * (1.0 - np.cos(2.0 * np.pi * (t - 1000.0) / 200.0)), 0.0)
    if float(np.max(np.abs(yw - expected))) > 1e-12:
        raise ValueError("saved waveform weights do not match the 1000-1200 s Hann window")
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    mask = (t >= 1000.0) & (t <= 1200.0)
    # Existing campaign definitions: energy=0.5*dt*sum(Hann*u^2); relative L2 uses YSpec as reference.
    ed = 0.5 * dt * float(np.sum(dw[mask] * x[mask] ** 2))
    ey = 0.5 * dt * float(np.sum(yw[mask] * yv[mask] ** 2))
    diff = x[mask] - yv[mask]
    y_norm = math.sqrt(float(np.sum(yw[mask] * yv[mask] ** 2)) * dt)
    weighted_l2 = math.sqrt(float(np.sum(yw[mask] * diff ** 2)) * dt)
    relative_l2 = weighted_l2 / max(y_norm, np.finfo(float).tiny)
    relative_energy = abs(ed - ey) / max(abs(ey), np.finfo(float).tiny)
    ds_saved, ys_saved = float(dresult["energy_m2_s"]), float(yresult["windowed_energy_m2_s"])
    ds_consistency = abs(ed - ds_saved) / max(abs(ds_saved), np.finfo(float).tiny)
    ys_consistency = abs(ey - ys_saved) / max(abs(ys_saved), np.finfo(float).tiny)
    gates = {"relative_weighted_l2_le_1e-3": relative_l2 <= 1e-3,
             "relative_energy_le_1e-3": relative_energy <= 1e-3,
             "dSpec_saved_energy_consistent_le_1e-3": ds_consistency <= 1e-3,
             "ySpec_saved_energy_consistent_le_1e-3": ys_consistency <= 1e-3}
    return {"window_s": [1000, 1200], "energy_definition": "0.5*dt*sum(Hann1000-1200*vertical_m^2)",
            "dSpec_energy_m2_s": ed, "ySpec_energy_m2_s": ey,
            "dSpec_saved_energy_m2_s": ds_saved, "ySpec_saved_energy_m2_s": ys_saved,
            "relative_energy_difference_using_ySpec_reference": relative_energy,
            "weighted_trace_l2_m_sqrt_s": weighted_l2, "ySpec_weighted_norm_m_sqrt_s": y_norm,
            "relative_weighted_l2_using_ySpec_reference": relative_l2,
            "saved_energy_relative_consistency": {"dSpec": ds_consistency, "ySpec": ys_consistency},
            "gates": gates, "status": "PASS" if all(gates.values()) else "FAIL",
            "compared_samples": n, "dSpec_samples": len(d), "ySpec_samples": len(y),
            "interpretation": "shared-prefix waveform comparison; 1e-3 limits are the existing campaign gates; relative L2 and energy use YSpec as reference"}


def cross_code_channels(root, jobs, listed, l):
    dsid, ysid = f"dsv6-native-l{l}", f"ysv6-native-l{l}"
    if jobs.get(dsid, {}).get("status") != "completed" or jobs.get(ysid, {}).get("status") != "completed":
        return {"status": "pending", "jobs": [dsid, ysid]}
    dmeta = read_json(root / dsid / "data/kernel_probe.json")
    ymeta = read_json(root / ysid / "data/result.json")
    ycmd = jobs[ysid].get("command", [])
    if "kernel-probe" not in ycmd or any("contact" in str(arg).lower() for arg in ycmd):
        raise ValueError(f"YSpec L{l} kernel command is not the expected contacts-off probe")
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


def radial_alignment(root, jobs, listed, l):
    yid, did = f"ysv6-native-l{l}", f"dsv6-native-l{l}"
    if jobs.get(yid, {}).get("status") != "completed" or jobs.get(did, {}).get("status") != "completed":
        return {"status": "pending", "jobs": [yid, did]}
    yjob, djob = jobs[yid], jobs[did]
    ymeta = read_json(root / yid / "data/result.json")
    remote = ymeta.get("kernel_sampling_input")
    expected_suffix = f"/profile-06/{did}/data/radial_nodes.csv"
    if not isinstance(remote, str) or not remote.endswith(expected_suffix):
        raise ValueError(f"YSpec radial input metadata does not map to the DSpec L{l} radial output")
    yd = yjob.get("inputs_sha256_before", {}).get(remote)
    dd = djob.get("outputs_sha256", {}).get("data/radial_nodes.csv")
    local = f"{did}/data/radial_nodes.csv"
    if not yd or yd != dd or local not in listed or listed[local] != dd:
        raise ValueError(f"state-bound exact radial-node identity failed for L{l}")
    if yjob.get("inputs_unchanged") is not True or yjob.get("inputs_sha256_after", {}).get(remote) != yd:
        raise ValueError(f"YSpec sampling input changed during L{l} job")
    return {"status": "PASS", "L": l, "snapshot_file": local, "sha256": dd,
            "mapping": "YSpec consumed this exact state-hash-bound DSpec radial-node file",
            "radial_interpolation": "none; no interpolation between radii or across interfaces"}


def dspec_qdot_checks(jobs, root):
    out = {}
    for jid, job in jobs.items():
        if not jid.startswith("dsv") or job.get("status") != "completed":
            continue
        meta = read_json(root / jid / "data/kernel_probe.json")
        if meta.get("source_contact_enabled") is not False:
            raise ValueError(f"expected contacts-off volume-only DSpec kernel for {jid}")
        value = float(meta["energy_adjoint_relative_error_vs_2E"])
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"invalid qdot/(2E) error for {jid}")
        out[jid] = {"energy_adjoint_relative_error_vs_2E": value,
                    "gate_le_1e-10": value <= 1e-10,
                    "actual_degree": read_json(root / jid / "data/result.json").get("degree")}
    return out


def collect(root: Path):
    state, jobs, listed = validate_snapshot(root)
    pairs = {}
    for l in (128, 384):
        for ver in ("v4", "v6"):
            pairs[f"ds-{ver}-kernel-native-profile-L{l}"] = file_pair(root, jobs, listed,
                f"ds{ver}-native-l{l}", f"ds{ver}-profiled-l{l}",
                ["data/kernel_probe.bin", "data/trace.csv", "data/spectrum.csv"])
    l = 1536
    pairs[f"ds-v6-kernel-native-profile-L{l}"] = file_pair(root, jobs, listed,
        f"dsv6-native-l{l}", f"dsv6-profiled-l{l}",
        ["data/kernel_probe.bin", "data/trace.csv", "data/spectrum.csv"])
    for l in (128, 384, 1536):
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
            "cleanup_verified": job.get("cleanup", {}).get("verified_gone")}
        if jid.startswith("dsv") and job.get("status") == "completed":
            telemetry[jid]["actual_degree_from_result"] = read_json(root / jid / "data/result.json").get("degree")
        elif jid.startswith("ysv6") and job.get("status") == "completed":
            telemetry[jid]["actual_degree_from_result"] = read_json(root / jid / "data/result.json").get("lmax")
    waveforms = {}
    cross_code = {}
    radial = {}
    for l in (128, 384, 1536):
        cross_code[str(l)] = cross_code_channels(root, jobs, listed, l)
        radial[str(l)] = radial_alignment(root, jobs, listed, l)
        required = [f"dsv6-native-l{l}", f"ysv6-native-wave-l{l}"]
        if all(jobs.get(x, {}).get("status") == "completed" for x in required):
            waveforms[str(l)] = waveform_metrics(root / required[0] / "data/trace.csv",
                root / required[1] / "data/trace.csv",
                read_json(root / required[0] / "data/result.json"),
                read_json(root / required[1] / "data/result.json"))
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
    parity_ok = all(x["status"] == "PASS" for x in pairs.values())
    qdot = dspec_qdot_checks(jobs, root)
    qdot_ok = all(v["gate_le_1e-10"] for v in qdot.values()) and len(qdot) == sum(
        1 for jid in EXPECTED_IDS if jid.startswith("dsv") and jobs.get(jid, {}).get("status") == "completed")
    radial_ok = all(v["status"] == "PASS" for v in radial.values())
    profile_ids = {jid for jid in EXPECTED_IDS if "profiled" in jid and jobs.get(jid, {}).get("status") == "completed"}
    profile_ok = set(profiles) == profile_ids and all(p["root_closure"]["closes_within_existing_tolerance"] for p in profiles.values())
    waveform_complete = all(v.get("status") == "PASS" for v in waveforms.values())
    complete = not pending and state.get("status") == "completed"
    summary_status = "INCOMPLETE" if not complete else (
        "PASS" if parity_ok and profile_ok and waveform_complete and qdot_ok and radial_ok else "FAIL")
    return {"schema": "component-profile-summary-v1", "status": summary_status,
        "snapshot_root": str(root.resolve()),
        "campaign_status_at_snapshot": state.get("status"), "task_id": state.get("task_id"),
        "thread_cap": state.get("thread_cap"), "claim_limits": [
            "Kernel probes use only the first 64 angular points for each L and are diagnostics, not full global kernels.",
            "L128/L384/L1536 cases are not a degree-convergence or representative global-throughput claim.",
            "The physical spatial grid is 64*179=11456 points; tensor_volume_point_count=595712 counts 52 frequency visits.",
            "All current kernel probes are contacts-off, volume-only diagnostics; separate source-shell contacts are omitted, so this is not a full physical global dJ certificate.",
            "YSpec waveform jobs are independent forward solves; YSpec kernel-probe energies use the shared DSpec baseline and are not treated as independent waveform energies.",
            "Radial comparisons use exact state-bound shared sampling files; no interpolation is performed between radii or across interfaces.",
            "Counters such as modal_lmax and solver_degree are accumulated counts; actual degree is in result metadata.",
            "Native/profiled wall-time differences are single-pass apparent measurements, not statistical estimates.",
            "MATHPROFILE inclusive scopes overlap; only exclusive stage sums are used for root-closure accounting."],
        "pending_jobs": pending, "telemetry": telemetry, "exact_parity": pairs,
        "profiled_stage_records": profiles, "native_profile_overhead_apparent": overhead,
        "waveform_cross_code_descriptive": waveforms, "kernel_cross_code_descriptive": cross_code,
        "radial_sampling_identity": radial, "dSpec_qdot_2E_checks": qdot,
        "acceptance_gates": {"all_26_expected_jobs_completed": complete,
            "all_native_profile_parity_passed": parity_ok, "profile_records_and_closure_passed": profile_ok,
            "waveform_1e-3_gates_passed": waveform_complete,
            "dSpec_qdot_2E_le_1e-10_passed": qdot_ok,
            "state_bound_radial_identity_passed": radial_ok},
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
