#!/usr/bin/env python3
"""Bounded same-grid Phase D timing and output-parity check."""
import argparse
import csv
import hashlib
import json
import math
import struct
from pathlib import Path


def read_json(path):
    return json.loads((Path(path) / "result.json").read_text())


def rows(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def rel_l2(path_a, path_b, columns):
    a, b = rows(path_a), rows(path_b)
    if len(a) != len(b):
        raise ValueError(f"row count differs: {path_a} vs {path_b}")
    name = Path(path_a).name
    grid_columns = {"trace.csv": ("time_s", "window"), "spectrum.csv": ("bin", "frequency_hz"), "components.csv": ("bin",)}[name]
    for ra, rb in zip(a, b):
        for key in grid_columns:
            if ra[key] != rb[key]:
                raise ValueError(f"grid/key differs in {key}: {ra[key]} vs {rb[key]}")
    num = den = 0.0
    for ra, rb in zip(a, b):
        for key in columns:
            x, y = float(ra[key]), float(rb[key])
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError(f"non-finite {key}")
            num += (x - y) ** 2
            den += y * y
    return math.sqrt(num / den) if den else (0.0 if num == 0 else math.inf)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def compare_csv_set(a, b):
    metrics = {}
    for key, name, cols in [
        ("trace_velocity", "trace.csv", ("velocity_m_s",)),
        ("trace_rate", "trace.csv", ("rate_m_s2",)),
        ("trace_adjoint_source", "trace.csv", ("adjoint_source",)),
        ("spectrum_complex", "spectrum.csv", ("real_m", "imag_m")),
        ("component_up", "components.csv", ("up_re", "up_im")),
        ("component_north", "components.csv", ("north_re", "north_im")),
        ("component_east", "components.csv", ("east_re", "east_im")),
    ]:
        metrics[key] = rel_l2(Path(a) / name, Path(b) / name, cols)
    return metrics


def rel_bin(path_a, path_b):
    ba, bb = Path(path_a).read_bytes(), Path(path_b).read_bytes()
    if len(ba) != len(bb) or len(ba) % 8:
        raise ValueError(f"binary output size differs or is not float64 aligned: {path_a}, {path_b}")
    aa = [x[0] for x in struct.iter_unpack("<d", ba)]
    ab = [x[0] for x in struct.iter_unpack("<d", bb)]
    if not all(math.isfinite(x) for x in aa + ab):
        raise ValueError("non-finite value in native angular output array")
    den = sum(x * x for x in ab)
    num = sum((x - y) ** 2 for x, y in zip(aa, ab))
    return math.sqrt(num / den) if den else (0.0 if num == 0 else math.inf)


def zero_support(path):
    data = rows(Path(path) / "phase-c2-zero-bin-support.csv")
    values = [int(r["bin"]) for r in data]
    if values != [7, 8, 41]:
        raise ValueError(f"exact-zero bins differ: {values}")
    for r in data:
        if r["taper_exact_zero"] != "1" or r["measured_P_exact_zero"] != "1" or r["forward_solve_retained"] != "1" or r["unit_adjoint_skipped"] != "1":
            raise ValueError(f"invalid exact-zero support row: {r}")
    return values


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--plan", required=True)
    p.add_argument("--batched", required=True)
    p.add_argument("--scalar", required=True)
    p.add_argument("--profile-batched", required=True)
    p.add_argument("--validated", required=True)
    p.add_argument("--yspec", required=True)
    p.add_argument("--accuracy-check", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    batch, scalar = read_json(args.batched), read_json(args.scalar)
    profiled, ys = read_json(args.profile_batched), read_json(args.yspec)
    accuracy = json.loads(Path(args.accuracy_check).read_text())
    plan = json.loads(Path(args.plan).read_text())
    if accuracy.get("status") != "pass":
        raise ValueError("selected .9h accuracy receipt is not PASS")
    for label, r, expected_batch, expected_scalar, freq, solve in [
        ("batched-native", batch, True, False, 8, 1),
        ("scalar-native", scalar, False, True, 8, 1),
        ("batched-profile", profiled, True, False, 8, 1),
        ("matched-yspec", ys, False, False, 1, 8),
    ]:
        if r.get("backend") != ("yspec" if label == "matched-yspec" else "dspec"):
            raise ValueError(f"{label}: backend mismatch")
        if bool(r.get("phase_d_batched_angular")) != expected_batch or bool(r.get("phase_d_scalar_only")) != expected_scalar:
            raise ValueError(f"{label}: route metadata mismatch")
        if r.get("degree") != 200 or r.get("nq") != 6 or r.get("frequency_count") != 35:
            raise ValueError(f"{label}: workload metadata mismatch")
        if r.get("frequency_workers") != freq or r.get("solve_workers") != solve:
            raise ValueError(f"{label}: worker metadata mismatch")
        if abs(r.get("phase_d_requested_relative_accuracy", -1) - 0.001) > 1e-15 or abs(r.get("relative_error", -1) - 1e-5) > 1e-15:
            raise ValueError(f"{label}: requested accuracy or legacy relative_error metadata mismatch")
        if abs(r.get("phase_d_selected_h", r.get("refinement", 0)) - 0.02035382429229197) > 1e-14:
            raise ValueError(f"{label}: selected h mismatch")
    commands = {j["id"]: j["command"] for j in plan["jobs"]}
    def driver_argv(job_id):
        command = commands[job_id]
        idx = next(i for i, value in enumerate(command) if value in ("/space/adcm2/research/builds/phase-d-nq6-batched-v2/phase-d-native/phase-d-native-driver", "/space/adcm2/research/builds/phase-d-nq6-batched-v2/phase-d-profile/phase-d-profile-driver"))
        return command, command[idx:]
    for job_id in ("phase-d-batched-native", "phase-d-scalar-native", "phase-d-batched-profile"):
        command, argv = driver_argv(job_id)
        if "DSPEC_FREQUENCY_WORKERS=8" not in command or argv[1] != "dspec" or argv[7] != "1":
            raise ValueError(f"{job_id}: expected DSpec frequency8/CLI worker1 command")
    ycmd, yargv = driver_argv("phase-d-matched-yspec")
    if "DSPEC_FREQUENCY_WORKERS=8" in ycmd or ycmd[1:3] != ["-u", "DSPEC_FREQUENCY_WORKERS"] or "OMP_NUM_THREADS=1" not in ycmd or "OMP_THREAD_LIMIT=8" not in ycmd or yargv[1] != "yspec" or yargv[7] != "8":
        raise ValueError("YSpec command does not explicitly unset DSpec frequency workers and use CLI workers8/library threads1")
    if batch.get("phase_c_zero_bin_skips") is not True or profiled.get("phase_c_zero_bin_skips") is not True:
        raise ValueError("DSpec exact-zero skips are not enabled")
    if batch.get("phase_d_common_validation") is not False or scalar.get("phase_d_common_validation") is not False:
        raise ValueError("timed jobs unexpectedly enabled common-field capture")
    if profiled.get("phase_d_actual_temporary_peak_bytes", 1 << 60) > 33554432:
        raise ValueError("profiled DSpec angular temporary bound exceeds 32 MiB")
    for name in ("nodes.csv", "along-angles.csv", "across-angles.csv", "frequency-grid.csv", "mesh-elements.csv"):
        hs = [digest(Path(d) / name) for d in (args.batched, args.scalar, args.profile_batched, args.yspec, args.validated)]
        if len(set(hs)) != 1:
            raise ValueError(f"requested sampling geometry differs for {name}")
    node_count = len(rows(Path(args.batched) / "nodes.csv"))
    along_count = len(rows(Path(args.batched) / "along-angles.csv"))
    across_count = len(rows(Path(args.batched) / "across-angles.csv"))
    yspec_binary_diagnostics = {}
    for name, angle_count in (("along.bin", along_count), ("across.bin", across_count)):
        path = Path(args.yspec) / name
        expected_bytes = node_count * angle_count * 8
        if path.stat().st_size != expected_bytes:
            raise ValueError(f"YSpec {name} shape differs from nodes x requested angles")
        yspec_binary_diagnostics[name] = {"expected_float64_values": node_count * angle_count, "relative_l2_to_dspec_diagnostic_only": rel_bin(path, Path(args.batched) / name)}
    bzero = zero_support(args.batched)
    szero = zero_support(args.scalar)
    pzero = zero_support(args.profile_batched)
    if not (bzero == szero == pzero):
        raise ValueError("exact-zero inventories differ")
    pair = compare_csv_set(args.scalar, args.batched)
    profile_parity = compare_csv_set(args.profile_batched, args.batched)
    validated_parity = compare_csv_set(args.batched, args.validated)
    field_pair = {name: rel_bin(Path(args.scalar) / name, Path(args.batched) / name) for name in ("along.bin", "across.bin")}
    field_profile = {name: rel_bin(Path(args.profile_batched) / name, Path(args.batched) / name) for name in ("along.bin", "across.bin")}
    field_validated = {name: rel_bin(Path(args.batched) / name, Path(args.validated) / name) for name in ("along.bin", "across.bin")}
    for label, metrics in (("scalar-v-batched", {**pair, **field_pair}), ("profile-v-native", {**profile_parity, **field_profile}), ("v2-batched-v-validated-v1", {**validated_parity, **field_validated})):
        if max(metrics.values()) > 1e-11:
            raise ValueError(f"{label} output parity exceeds 1e-11: {metrics}")
    prof = json.loads((Path(args.profile_batched) / "profile.json").read_text())
    counters = {(x["phase"], x["name"].removeprefix("frequency_worker.")): int(x["value"]) for x in prof.get("counters", [])}
    phase = "combined_forward_unit_adjoint_reconstruction_contraction"
    required = {
        "phase_c2_forward_spectrum_frequency_calls": 35,
        "phase_c2_unit_adjoint_frequency_calls": 32,
        "phase_c2_exact_zero_unit_adjoint_skips": 3,
    }
    observed = {k: counters.get((phase, k), 0) for k in required}
    if observed != required:
        raise ValueError(f"combined-pass work counts differ: {observed}")
    workspace_calls = sum(int(r["calls"]) for r in prof.get("rows", []) if r.get("phase") == phase and r.get("scope_name", "").removeprefix("frequency_worker.") == "frequency_workspace_prepare_and_modal_capture")
    if workspace_calls != 35:
        raise ValueError(f"combined-pass workspace preparations expected35 (all forward bins), got {workspace_calls}")
    route_scope = [r for r in prof.get("rows", []) if r.get("phase") == phase and r.get("scope_name", "").removeprefix("frequency_worker.") == "phase_d_tiled_eigen_angular_reconstruction"]
    worker_angular_seconds = sum(float(r["inclusive_s"]) for r in route_scope)
    if sum(int(r["calls"]) for r in route_scope) <= 0:
        raise ValueError("profile has no tiled angular work")
    yspec_diag = compare_csv_set(args.yspec, args.batched)
    out = {
        "status": "pass",
        "selected_mesh_accuracy_receipt": str(args.accuracy_check),
        "same_grid_scalar_vs_batched_relative_l2_by_channel": {**pair, **field_pair},
        "same_grid_profile_vs_native_relative_l2_by_channel": {**profile_parity, **field_profile},
        "same_grid_v2_batched_vs_validated_v1_relative_l2_by_channel": {**validated_parity, **field_validated},
        "exact_zero_bins": bzero,
        "combined_pass_work_counts": observed,
        "combined_pass_frequency_workspace_preparations_all_forward_bins": workspace_calls,
        "profiled_tiled_angular_scope_calls": sum(int(r["calls"]) for r in route_scope),
        "profiled_tiled_angular_worker_elapsed_sum_s_overlapping_not_master_wall": worker_angular_seconds,
        "timings": {"batched_native_driver_s": batch["elapsed_s"], "scalar_native_driver_s": scalar["elapsed_s"], "batched_profile_driver_s": profiled["elapsed_s"], "yspec_driver_s": ys["elapsed_s"]},
        "yspec_backend_output_relative_l2_diagnostics_only_by_channel": yspec_diag,
        "yspec_native_array_diagnostics_only": yspec_binary_diagnostics,
        "grid_sha256": {name: digest(Path(args.batched) / name) for name in ("nodes.csv", "along-angles.csv", "across-angles.csv", "frequency-grid.csv", "mesh-elements.csv")},
        "sampling_gate": "exact byte-equal nodes, angle lists, frequency grid and mesh elements across all four jobs",
        "worker_note": "DSpec uses frequency_workers=8/solve_workers=1; YSpec explicitly unsets DSPEC_FREQUENCY_WORKERS and uses CLI workers=8 with OMP_NUM_THREADS=1/OMP_THREAD_LIMIT=8; both plan command and result metadata are checked.",
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
