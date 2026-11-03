#!/usr/bin/env python3
"""Data-only checks for the six reviewed Tunya L=3 kernel jobs.

This program never launches scientific executables. It reads completed job
directories and emits a JSON report, including descriptive cross-code errors;
those cross-code differences are deliberately not a global accuracy gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

JOBS = (
    "dsv4-native-l3", "dsv4-profiled-l3", "dsv6-native-l3",
    "dsv6-profiled-l3", "ysv6-native-l3", "ysv6-profiled-l3",
)
DS = JOBS[:4]
YS = JOBS[4:]
RADIAL_COUNT = 179
PARITY_LIMIT = 2.0e-13
CHUNK_POINTS = 8192
KERNEL_CHANNELS = {"dspec": ("A", "C", "F", "L", "N", "kappa"),
                   "yspec": ("A", "C", "F", "L", "N")}
DS_PHYSICAL_FIELDS = ("mode", "degree", "sem_degree", "nq_gll_points_per_element",
                      "polynomial_degree", "radial_refinement_maxstep", "qex", "dt_s",
                      "duration_s", "fft_samples", "fft_period_s", "frequency_bins",
                      "epsilon_si_s_inverse", "df_si_hz", "source", "receiver",
                      "source_element", "source_element_fluid", "receiver_element",
                      "receiver_element_fluid", "dc_included", "nyquist_raw_solved",
                      "pole_only_nonradial_adjoint", "forward_only_baseline",
                      "source_contact_enabled")
ENERGY_RESULT_FIELDS = ("energy_m2_s", "energy_adjoint_qdotp",
                        "energy_adjoint_relative_error_vs_2E")
ENERGY_DOMAIN_FIELDS = ("energy_adjoint_qdotp", "energy_adjoint_relative_error_vs_2E")


class CheckError(Exception):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_file(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise CheckError(f"expected a JSON object: {path}")
    def check(value, where):
        if isinstance(value, float) and not math.isfinite(value):
            raise CheckError(f"non-finite JSON number at {where} in {path}")
        if isinstance(value, dict):
            for key, child in value.items():
                check(child, f"{where}.{key}")
        elif isinstance(value, list):
            for index, child in enumerate(value):
                check(child, f"{where}[{index}]")
    check(value, "root")
    return value


def finite_number(value, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise CheckError(f"{label} is not numeric") from exc
    if not math.isfinite(number):
        raise CheckError(f"{label} is non-finite")
    return number


def numeric_csv(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise CheckError(f"empty CSV: {path}") from exc
        rows = []
        for line_number, row in enumerate(reader, 2):
            if not row or all(not item.strip() for item in row):
                continue
            if len(row) != len(header):
                raise CheckError(f"{path}:{line_number}: expected {len(header)} columns, got {len(row)}")
            try:
                values = [float(item) for item in row]
            except ValueError as exc:
                raise CheckError(f"{path}:{line_number}: nonnumeric cell") from exc
            if not all(math.isfinite(value) for value in values):
                raise CheckError(f"{path}:{line_number}: non-finite cell")
            rows.append(values)
    if not rows:
        raise CheckError(f"CSV has no numeric rows: {path}")
    return header, np.asarray(rows, dtype=np.float64)


def profile_check(log: Path, expected: str, required: bool) -> dict:
    lines = [line.split("MATHPROFILE ", 1)[1].strip()
             for line in log.read_text(encoding="utf-8", errors="replace").splitlines()
             if "MATHPROFILE " in line]
    if not lines and required:
        raise CheckError(f"profiled job has no MATHPROFILE record: {log}")
    if len(lines) > 1:
        raise CheckError(f"expected at most one MATHPROFILE record in {log}, found {len(lines)}")
    if not lines:
        return {"present": False}
    try:
        record = json.loads(lines[0])
    except json.JSONDecodeError as exc:
        raise CheckError(f"invalid MATHPROFILE JSON in {log}: {exc}") from exc
    solver = record.get("solver")
    if expected == "dspec" and (solver != "dspec" or record.get("objective") != -1):
        raise CheckError(f"unexpected DSpec profile identity/objective in {log}")
    if expected == "yspec" and (solver != "YSpecpp.physical_kernel" or record.get("objective") != 0):
        raise CheckError(f"unexpected YSpec profile identity/objective in {log}")
    stages = record.get("stages")
    if not isinstance(stages, dict) or "physical_case_total" not in stages:
        raise CheckError(f"profile has no physical_case_total root in {log}")
    exclusive_sum = 0.0
    for name, stage in stages.items():
        if not isinstance(stage, dict):
            raise CheckError(f"malformed profile stage {name} in {log}")
        for key in ("inclusive_s", "exclusive_s"):
            value = finite_number(stage.get(key), f"{name}.{key}")
            if value < 0:
                raise CheckError(f"negative profile time {name}.{key} in {log}")
        exclusive_sum += float(stage["exclusive_s"])
        calls = finite_number(stage.get("calls"), f"{name}.calls")
        if calls < 0:
            raise CheckError(f"negative profile calls for {name} in {log}")
    root = float(stages["physical_case_total"]["inclusive_s"])
    delta = abs(exclusive_sum - root)
    tolerance = max(1.0e-8, root * 2.0e-12)
    if delta > tolerance:
        raise CheckError(f"profile exclusive sum {exclusive_sum:.17g}s != root inclusive {root:.17g}s in {log}")
    for name, counter in record.get("counters", {}).items():
        value = finite_number(counter, f"counter {name}")
        if value < 0:
            raise CheckError(f"negative profile counter {name} in {log}")
    return {"present": True, "solver": solver, "objective": record.get("objective"),
            "stage_count": len(stages), "exclusive_sum_s": exclusive_sum,
            "root_inclusive_s": root, "accounting_abs_error_s": delta}


def read_job(task: Path, job_id: str) -> dict:
    root = task / job_id
    data = root / "data"
    if not root.is_dir() or not data.is_dir():
        raise CheckError(f"missing expected job/data directory: {root}")
    result = json_file(data / "result.json")
    if job_id.startswith("dsv"):
        trace_header, trace = numeric_csv(data / "trace.csv")
        spectrum_header, spectrum = numeric_csv(data / "spectrum.csv")
        if trace_header != ["time_s", "vertical_m", "window_weight"]:
            raise CheckError(f"unexpected DSpec trace columns for {job_id}: {trace_header}")
        if spectrum_header != ["bin", "frequency_hz", "real", "imag"]:
            raise CheckError(f"unexpected DSpec spectrum columns for {job_id}: {spectrum_header}")
    else:
        # Kernel-mode YSpec writes neither table; the baseline path in result.json
        # is its shared input, read here only to reject non-finite input values.
        baseline = Path(str(result.get("kernel_baseline_trace", "")))
        if not baseline.is_file():
            raise CheckError(f"YSpec baseline trace is unavailable: {baseline}")
        trace_header, trace = numeric_csv(baseline)
        spectrum_header, spectrum = [], np.empty((0, 0), dtype=np.float64)
    if not (data / ("kernel.bin" if job_id.startswith("dsv") else "kernels.bin")).is_file():
        raise CheckError(f"missing kernel output for {job_id}")
    profiled = job_id.endswith("profiled-l3")
    profile = profile_check(root / "job.log", "dspec" if job_id.startswith("dsv") else "yspec", profiled)
    if finite_number(result.get("lmax", result.get("degree")), f"{job_id}.lmax") != 3:
        raise CheckError(f"{job_id} is not an L=3 result")
    if len(trace_header) != trace.shape[1] or len(spectrum_header) != spectrum.shape[1]:
        raise CheckError(f"CSV header/data shape mismatch for {job_id}")
    return {"root": root, "data": data, "result": result,
            "trace_header": trace_header, "trace": trace,
            "spectrum_header": spectrum_header, "spectrum": spectrum,
            "profile": profile}


def read_radial(path: Path, dspec: bool) -> tuple[list[str], np.ndarray]:
    header, data = numeric_csv(path)
    expected = ("radial_index", "element", "node", "radius_m", "element_width_m", "fluid") if dspec else (
        "element", "node", "radius_m", "element_width_m", "fluid")
    if tuple(header) != expected:
        raise CheckError(f"unexpected radial-node columns in {path}: {header}")
    if data.shape[0] != RADIAL_COUNT:
        raise CheckError(f"expected {RADIAL_COUNT} radial nodes in {path}, found {data.shape[0]}")
    return header, data


def load_kernel(job: dict, job_id: str) -> tuple[Path, int, int, int, np.ndarray]:
    data = job["data"]
    result = job["result"]
    ds = job_id.startswith("dsv")
    binary = data / ("kernel.bin" if ds else "kernels.bin")
    if ds:
        domain = json_file(data / "kernel_domain.json")
        layout = domain.get("binary_layout")
        if not isinstance(layout, str) or "point-major" not in layout or "doubles[A,C,F,L,N,kappa]" not in layout:
            raise CheckError(f"unexpected DSpec kernel binary type/order in {job_id}: {layout!r}")
        job["domain"] = domain
        if domain.get("diagnostic_angle_subset_only") is True or domain.get("sampled_angle_subset_is_complete") is False:
            raise CheckError(f"{job_id} is an angular subset, not a full L=3 grid")
        if domain.get("selected_radial_trace_count") != RADIAL_COUNT:
            raise CheckError(f"{job_id} does not include all 179 selected radial traces")
        if domain.get("available_radial_trace_sides") != RADIAL_COUNT:
            raise CheckError(f"{job_id} model does not expose 179 radial trace sides")
        if domain.get("full_target_point_count") != domain.get("actual_sampled_point_count"):
            raise CheckError(f"{job_id} sampled fewer points than its complete target")
        if domain.get("actual_kernel_binary_bytes") != binary.stat().st_size:
            raise CheckError(f"{job_id} kernel byte count disagrees with kernel_domain.json")
        nt, np_ = int(domain["n_theta"]), int(domain["n_phi"])
        points = int(domain["actual_sampled_point_count"])
        channels = 6
        radial_path = data / "kernel_radial_nodes.csv"
        read_radial(radial_path, True)
        if result.get("mode") != "kernel":
            raise CheckError(f"{job_id} mode is not kernel")
        energy = finite_number(result.get("energy_m2_s"), f"{job_id}.energy_m2_s")
        qdotp = finite_number(result.get("energy_adjoint_qdotp"), f"{job_id}.energy_adjoint_qdotp")
        if energy <= 0 or abs(qdotp - 2.0 * energy) / abs(2.0 * energy) > 1.0e-10:
            raise CheckError(f"{job_id} q dot P(B raw) / (2E) gate exceeds 1e-10")
    else:
        if result.get("mode") != "kernel" or result.get("kernel_probe_mode") is True or result.get("kernel_probe_subset") is True:
            raise CheckError(f"{job_id} is not a complete kernel run")
        if int(result.get("kernel_radial_blocks", -1)) != RADIAL_COUNT:
            raise CheckError(f"{job_id} does not include all 179 global radial nodes")
        nt, np_ = int(result["kernel_n_theta"]), int(result["kernel_n_phi"])
        points = int(result["kernel_point_count"])
        channels = 5
        layout = result.get("kernel_binary_layout")
        if (not isinstance(layout, str) or "float64 native endian" not in layout or
                "angle-major, then radial node, then A,C,F,L,N parameter" not in layout):
            raise CheckError(f"unexpected YSpec kernel binary type/order in {job_id}: {layout!r}")
        if result.get("kernel_output_file") != "kernels.bin":
            raise CheckError(f"{job_id} result points to a different kernel file")
        if not (data / "radial_nodes.csv").is_file():
            sampling = Path(str(result.get("kernel_sampling_input", "")))
            if not sampling.is_file():
                raise CheckError(f"{job_id} radial sampling input is unavailable")
            radial_path = sampling
        else:
            radial_path = data / "radial_nodes.csv"
        read_radial(radial_path, False)
    if nt != 7 or points != nt * np_ * RADIAL_COUNT:
        raise CheckError(f"{job_id} does not contain its complete L=3 angular grid x 179 nodes")
    expected_bytes = points * channels * 8
    if binary.stat().st_size != expected_bytes:
        raise CheckError(f"{job_id} kernel size {binary.stat().st_size} != expected {expected_bytes}")
    return binary, points, channels, np_, np.fromfile(binary, dtype=np.float64)


def compare_arrays(path_a: Path, path_b: Path, limit: float, exact: bool = False,
                   channel_names: tuple[str, ...] = ()) -> dict:
    size_a, size_b = path_a.stat().st_size, path_b.stat().st_size
    if size_a != size_b or size_a % 8:
        raise CheckError(f"comparison shape/byte mismatch: {path_a} vs {path_b}")
    width = len(channel_names)
    if width == 0 or (size_a // 8) % width:
        raise CheckError(f"kernel size is not divisible by {width} channels: {path_a}")
    exact_bytes = sha256(path_a) == sha256(path_b)
    points = size_a // (8 * width)
    max_abs = np.zeros(width, dtype=np.float64)
    max_reference = np.zeros(width, dtype=np.float64)
    squared = np.zeros(width, dtype=np.float64)
    if not exact_bytes:
        for point_offset in range(0, points, CHUNK_POINTS):
            npoints = min(CHUNK_POINTS, points - point_offset)
            nvalues = npoints * width
            a = np.fromfile(path_a, dtype=np.float64, count=nvalues, offset=point_offset * width * 8).reshape(npoints, width)
            b = np.fromfile(path_b, dtype=np.float64, count=nvalues, offset=point_offset * width * 8).reshape(npoints, width)
            if a.size != nvalues or b.size != nvalues or not np.isfinite(a).all() or not np.isfinite(b).all():
                raise CheckError(f"non-finite/truncated values in comparison: {path_a}, {path_b}")
            d = np.abs(a - b)
            max_abs = np.maximum(max_abs, d.max(axis=0))
            max_reference = np.maximum(max_reference, np.maximum(np.abs(a).max(axis=0), np.abs(b).max(axis=0)))
            squared += np.einsum("ij,ij->j", d, d)
    scales = np.maximum(max_reference, np.finfo(np.float64).tiny)
    scaled = np.zeros(width) if exact_bytes else max_abs / scales
    metrics = {name: {"max_abs_error": float(max_abs[i]),
                      "scale_max_abs_value": float(max_reference[i]),
                      "scaled_max_abs_error": float(scaled[i]),
                      "relative_l2_error": math.sqrt(float(squared[i])) / float(scales[i])}
               for i, name in enumerate(channel_names)}
    if exact and not exact_bytes:
        raise CheckError(f"required exact parity failed: {path_a} vs {path_b}")
    bad = {name: values["scaled_max_abs_error"] for name, values in metrics.items()
           if values["scaled_max_abs_error"] > limit}
    if not exact and bad:
        raise CheckError(f"per-channel scaled maximum error exceeds {limit:.6g}: {bad}; {path_a} vs {path_b}")
    return {"byte_exact": exact_bytes, "channels": metrics,
            "maximum_channel_scaled_error": float(scaled.max(initial=0.0))}


def compare_numeric(a: np.ndarray, b: np.ndarray, exact: bool, label: str,
                    column_names: list[str] | None = None) -> dict:
    if a.shape != b.shape:
        raise CheckError(f"numeric table shape mismatch for {label}: {a.shape} vs {b.shape}")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise CheckError(f"non-finite values in {label}")
    errors = {}
    for column in range(a.shape[1]):
        av, bv = a[:, column], b[:, column]
        d = np.abs(av - bv)
        max_abs = float(d.max(initial=0.0))
        scale = max(float(np.abs(av).max(initial=0.0)), float(np.abs(bv).max(initial=0.0)), np.finfo(float).tiny)
        scaled = max_abs / scale
        name = column_names[column] if column_names is not None else f"column_{column}"
        errors[name] = {"max_abs_error": max_abs, "scale_max_abs_value": scale,
                        "scaled_max_abs_error": scaled, "exact_values": bool(np.array_equal(av, bv))}
        if exact and not errors[name]["exact_values"]:
            raise CheckError(f"required exact numeric parity failed for {label} {name}")
        if not exact and scaled > PARITY_LIMIT:
            raise CheckError(f"per-column scaled error {scaled:.6g} exceeds {PARITY_LIMIT:.6g} for {label} {name}")
    return {"exact_values": all(value["exact_values"] for value in errors.values()),
            "maximum_column_scaled_error": max((v["scaled_max_abs_error"] for v in errors.values()), default=0.0),
            "columns": errors}


def compare_job_outputs(a: dict, b: dict, exact: bool, label: str) -> dict:
    channel_names = KERNEL_CHANNELS["yspec" if label.startswith("ysv6") else "dspec"]
    result = {"kernel": compare_arrays(a["kernel_path"], b["kernel_path"], PARITY_LIMIT,
                                       exact, channel_names)}
    if label.startswith("ysv6"):
        volatile = {"peak_rss_bytes", "timings_s", "profile_stages"}
        metadata_a = {k: v for k, v in a["result"].items() if k not in volatile}
        metadata_b = {k: v for k, v in b["result"].items() if k not in volatile}
        if metadata_a != metadata_b:
            raise CheckError(f"relevant YSpec result metadata differs in {label}")
        if a["result"].get("kernel_baseline_trace") != b["result"].get("kernel_baseline_trace"):
            raise CheckError(f"YSpec native/profile baseline paths differ in {label}")
        result["relevant_result_metadata_exact"] = True
        result["baseline_trace_finite_only"] = True
    else:
        for field in ("trace", "spectrum"):
            if a[f"{field}_header"] != b[f"{field}_header"]:
                raise CheckError(f"{field} headers differ in {label}")
            result[field] = compare_numeric(a[field], b[field], exact, f"{label} {field}",
                                            a[f"{field}_header"])
    return result


def cross_code_summary(dsjob: dict, ysjob: dict, ds_radial: np.ndarray, ys_radial: np.ndarray) -> dict:
    ds_points = dsjob["kernel_values"].reshape((-1, 6))
    ys_points = ysjob["kernel_values"].reshape((-1, 5))
    if ds_points.shape[0] != ys_points.shape[0]:
        raise CheckError("cross-code point counts differ")
    if ds_radial.shape[0] != ys_radial.shape[0]:
        raise CheckError("cross-code radial-node counts differ")
    for ds_idx, ys_idx, name in ((1, 0, "element"), (2, 1, "node")):
        if not np.array_equal(ds_radial[:, ds_idx], ys_radial[:, ys_idx]):
            raise CheckError(f"cross-code radial {name} indices differ")
    for ds_idx, ys_idx, name in ((3, 2, "radius_m"), (4, 3, "element_width_m")):
        if not np.allclose(ds_radial[:, ds_idx], ys_radial[:, ys_idx], rtol=0, atol=1e-8):
            raise CheckError(f"cross-code radial {name} differs")
    if not np.array_equal(ds_radial[:, 5], ys_radial[:, 4]):
        raise CheckError("cross-code radial fluid flags differ")
    per_parameter = {}
    for dcol, ycol, name in ((0, 0, "A"), (1, 1, "C"), (2, 2, "F"), (3, 3, "L"), (4, 4, "N")):
        a, b = ds_points[:, dcol], ys_points[:, ycol]
        d = np.abs(a - b)
        scale = max(float(np.abs(a).max(initial=0.0)), float(np.abs(b).max(initial=0.0)), np.finfo(float).tiny)
        per_parameter[name] = {"max_abs_error": float(d.max(initial=0.0)), "scaled_max_abs_error": float(d.max(initial=0.0)) / scale}
    # Point-major layout is angular point then radial node. Select the regular
    # outer-core shell (2200-3200 km), retaining every angle at each qualifying node.
    radial_mask = ((ds_radial[:, 3] >= 2.2e6) & (ds_radial[:, 3] <= 3.2e6) &
                   (ds_radial[:, 5] == 1) & (ys_radial[:, 4] == 1))
    if not radial_mask.any():
        raise CheckError("no fluid nodes in regular 2200-3200 km outer-core shell")
    radial_count = ds_radial.shape[0]
    angles = ds_points.shape[0] // radial_count
    if ds_points.shape[0] != angles * radial_count:
        raise CheckError("kernel point count is not divisible by radial-node count")
    mask = np.tile(radial_mask, angles)
    dsk = ds_points[mask, 5]
    ysk = ys_points[mask, 0] + ys_points[mask, 1] + ys_points[mask, 2]
    diff = np.abs(dsk - ysk)
    scale = max(float(np.abs(dsk).max(initial=0.0)), float(np.abs(ysk).max(initial=0.0)), np.finfo(float).tiny)
    per_parameter["kappa_vs_A_plus_C_plus_F_fluid_2200_3200km"] = {
        "comparison_points": int(mask.sum()), "max_abs_error": float(diff.max(initial=0.0)),
        "scaled_max_abs_error": float(diff.max(initial=0.0)) / scale,
        "gate": "descriptive_only_no_global_accuracy_claim"}
    return {"per_parameter": per_parameter,
            "interpretation": "descriptive cross-code differences only; no global scientific accuracy gate"}


def stable_metadata(record: dict, excluded: set[str]) -> dict:
    return {key: value for key, value in record.items()
            if key not in excluded and "path" not in key.lower() and key != "model"}


def require_close_energy(left: dict, right: dict, fields: tuple[str, ...], label: str) -> dict:
    metrics = {}
    for field in fields:
        a = finite_number(left.get(field), f"{label}.{field}")
        b = finite_number(right.get(field), f"{label}.{field}")
        difference = abs(a - b)
        if "relative_error" in field:
            tolerance = max(1.0e-14, PARITY_LIMIT * max(abs(a), abs(b)))
            scaled = difference / max(abs(a), abs(b), np.finfo(float).tiny)
        else:
            tolerance = PARITY_LIMIT * max(abs(a), abs(b), np.finfo(float).tiny)
            scaled = difference / max(abs(a), abs(b), np.finfo(float).tiny)
        if difference > tolerance:
            raise CheckError(f"cross-version energy field {label}.{field} differs: {a} vs {b} (abs {difference}, tol {tolerance})")
        metrics[field] = {"left": a, "right": b, "abs_difference": difference,
                          "scaled_difference": scaled, "absolute_tolerance": tolerance}
    return metrics


def dspec_metadata_gates(loaded: dict) -> dict:
    result_transient = {"elapsed_s"}
    domain_transient = {"configured_output_cap_bytes", "available_disk_bytes_at_preflight"}
    native_profile = {}
    for native, profiled in (("dsv4-native-l3", "dsv4-profiled-l3"),
                             ("dsv6-native-l3", "dsv6-profiled-l3")):
        left, right = loaded[native], loaded[profiled]
        if stable_metadata(left["result"], result_transient) != stable_metadata(right["result"], result_transient):
            raise CheckError(f"stable DSpec result metadata differs between {native} and {profiled}")
        if stable_metadata(left["domain"], domain_transient) != stable_metadata(right["domain"], domain_transient):
            raise CheckError(f"stable DSpec kernel-domain metadata differs between {native} and {profiled}")
        native_profile[f"{native}__{profiled}"] = True

    cross_version = {}
    domain_fields = ("selected_radial_trace_count", "available_radial_trace_sides", "n_theta",
                     "n_phi", "full_target_point_count", "actual_sampled_point_count",
                     "binary_layout", "kernel_is_full_earth_volume_integral", "source_contact_enabled")
    for native, profiled in (("dsv4-native-l3", "dsv6-native-l3"),
                             ("dsv4-profiled-l3", "dsv6-profiled-l3")):
        left, right = loaded[native], loaded[profiled]
        for field in DS_PHYSICAL_FIELDS:
            if field not in left["result"] or field not in right["result"]:
                raise CheckError(f"missing DSpec physical option {field} in cross-version metadata")
            if left["result"][field] != right["result"][field]:
                raise CheckError(f"cross-version DSpec physical option {field} differs: {left['result'][field]!r} vs {right['result'][field]!r}")
        for field in ("dc_included", "pole_only_nonradial_adjoint", "source_contact_enabled",
                      "forward_only_baseline"):
            if not isinstance(left["result"][field], bool) or not isinstance(right["result"][field], bool):
                raise CheckError(f"DSpec option {field} must be a JSON boolean")
        for field in domain_fields:
            if field not in left["domain"] or field not in right["domain"]:
                raise CheckError(f"missing DSpec kernel-domain field {field}")
            if left["domain"][field] != right["domain"][field]:
                raise CheckError(f"cross-version DSpec domain field {field} differs")
        energy_result = require_close_energy(left["result"], right["result"], ENERGY_RESULT_FIELDS,
                                             f"{native}__{profiled}.result")
        energy_domain = require_close_energy(left["domain"], right["domain"], ENERGY_DOMAIN_FIELDS,
                                             f"{native}__{profiled}.domain")
        cross_version[f"{native}__{profiled}"] = {"physical_options_exact": True,
                                                     "domain_options_exact": True,
                                                     "result_energy": energy_result,
                                                     "domain_energy": energy_domain}
    return {"native_profile_stable_result_and_domain_metadata_exact": native_profile,
            "v4_v6_physical_options_and_domain_consistent": cross_version}


def run(task: Path) -> dict:
    missing = [name for name in JOBS if not (task / name / "data").is_dir()]
    if missing:
        raise CheckError(f"task is missing completed job/data directories: {', '.join(missing)}")
    loaded = {name: read_job(task, name) for name in JOBS}
    for name, job in loaded.items():
        binary, points, channels, nphi, values = load_kernel(job, name)
        if not np.isfinite(values).all():
            raise CheckError(f"non-finite kernel value in {name}")
        expected = points * channels
        if values.size != expected:
            raise CheckError(f"kernel shape mismatch for {name}: {values.size} != {expected}")
        job["kernel_path"] = binary
        job["kernel_values"] = values
        job["points"] = points
        job["nphi"] = nphi
    dspec_metadata = dspec_metadata_gates(loaded)
    # Validate each metadata/CSV result against every matching source-code variant.
    comparisons = {}
    for left, right, exact in (
        ("dsv4-native-l3", "dsv4-profiled-l3", True),
        ("dsv6-native-l3", "dsv6-profiled-l3", True),
        ("dsv4-native-l3", "dsv6-native-l3", False),
        ("dsv4-profiled-l3", "dsv6-profiled-l3", False),
        ("ysv6-native-l3", "ysv6-profiled-l3", True),
    ):
        comparisons[f"{left}__{right}"] = compare_job_outputs(loaded[left], loaded[right], exact, f"{left} vs {right}")
    # DSpec radial identity across v4/v6 and profiles, and YSpec radial input parity.
    ds_radials = {}
    for name in DS:
        _, data = read_radial(loaded[name]["data"] / "kernel_radial_nodes.csv", True)
        ds_radials[name] = data
    for a, b in ((DS[0], DS[1]), (DS[0], DS[2]), (DS[2], DS[3])):
        compare_numeric(ds_radials[a], ds_radials[b], False, f"{a} vs {b} radial nodes")
    ys_radials = {}
    for name in YS:
        result = loaded[name]["result"]
        p = loaded[name]["data"] / "radial_nodes.csv"
        if not p.is_file():
            p = Path(str(result.get("kernel_sampling_input", "")))
        _, ys_radials[name] = read_radial(p, False)
    _, ds_radial_no_index = read_radial(loaded["dsv6-native-l3"]["data"] / "kernel_radial_nodes.csv", True)
    cross_code = {
        "native": cross_code_summary(loaded["dsv6-native-l3"], loaded["ysv6-native-l3"], ds_radials["dsv6-native-l3"], ys_radials["ysv6-native-l3"]),
        "profiled": cross_code_summary(loaded["dsv6-profiled-l3"], loaded["ysv6-profiled-l3"], ds_radials["dsv6-profiled-l3"], ys_radials["ysv6-profiled-l3"]),
    }
    compare_numeric(ys_radials["ysv6-native-l3"], ys_radials["ysv6-profiled-l3"], True, "YSpec native/profile radial nodes")
    return {"schema": 1, "status": "PASS", "task_dir": str(task),
            "jobs": {name: {"kernel_sha256": sha256(job["kernel_path"]),
                            "kernel_points": job["points"], "kernel_nphi": job["nphi"],
                            "kernel_channels": 6 if name.startswith("dsv") else 5,
                            "trace_rows": int(job["trace"].shape[0]),
                            "spectrum_rows": int(job["spectrum"].shape[0]),
                            "profile": job["profile"]}
                     for name, job in loaded.items()},
            "parity_gates": comparisons, "dspec_metadata_gates": dspec_metadata,
            "cross_code_descriptive": cross_code,
            "accuracy_claim": "No global cross-code accuracy claim is made; near-source differences remain unresolved.",
            "gates": {"all_full_l3_grids": True, "all_179_radial_nodes": True,
                      "native_profile_parity": True, "dsv4_dsv6_parity": True,
                      "yspec_native_profile_exact": True, "dspec_qdotp_over_2e_le_1e-10": True,
                      "finite_outputs": True}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("TASK_DIR", type=Path)
    parser.add_argument("OUTPUT_JSON", type=Path)
    args = parser.parse_args()
    task = args.TASK_DIR.resolve(strict=True)
    output = args.OUTPUT_JSON.resolve(strict=False)
    if output.exists():
        parser.error(f"refusing to overwrite existing output JSON: {output}")
    report = {"schema": 1, "status": "FAIL", "task_dir": str(task), "errors": []}
    try:
        report = run(task)
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(tmp_name, output)
    finally:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
    print(f"{report['status']}: {output}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
