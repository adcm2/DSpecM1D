#!/usr/bin/env python3
"""Validate Study 1 fixed probes, section diagnostics, and L3 directional checks."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import struct
from pathlib import Path
from typing import Any


JOB_IDS = (
    "l200-dspec-fixed-probe",
    "l200-yspec-fixed-probe",
    "l200-dspec-h08-fixed-probe",
    "l200-yspec-tol1e-7-fixed-probe",
    "l3-dspec-volume-audit",
    "l3-yspec-volume-audit",
    "l3-dspec-fd-e1e-4-plus",
    "l3-dspec-fd-e1e-4-minus",
    "l3-dspec-fd-e5e-5-plus",
    "l3-dspec-fd-e5e-5-minus",
)
BASE_D = JOB_IDS[0]
BASE_Y = JOB_IDS[1]
FINE_D = JOB_IDS[2]
TIGHT_Y = JOB_IDS[3]
AUDIT_D = JOB_IDS[4]
AUDIT_Y = JOB_IDS[5]
FD_PAIRS = ((JOB_IDS[6], JOB_IDS[7], 1e-4), (JOB_IDS[8], JOB_IDS[9], 5e-5))
ANGULAR_LIMIT_RADIAL_M = (25_000.0, 50_000.0, 100_000.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def finite(value: str | float, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} is non-finite")
    return number


def close_scaled(actual: float, expected: float, relative_tolerance: float, label: str,
                 scale: float | None = None) -> None:
    magnitude = max(abs(actual), abs(expected), abs(scale) if scale is not None else 0.0)
    if magnitude == 0.0:
        if actual != expected:
            raise ValueError(f"{label} differs in an exact-zero comparison")
        return
    if abs(actual - expected) > relative_tolerance * magnitude:
        raise ValueError(f"{label} differs beyond its operand-scaled relative tolerance")


def close_complex(actual: complex, expected: complex, relative_tolerance: float, label: str,
                  scale: float | None = None) -> None:
    magnitude = max(abs(actual), abs(expected), abs(scale) if scale is not None else 0.0)
    if magnitude == 0.0:
        if actual != expected:
            raise ValueError(f"{label} differs in an exact-zero comparison")
        return
    if abs(actual - expected) > relative_tolerance * magnitude:
        raise ValueError(f"{label} differs beyond its operand-scaled relative tolerance")


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} is not a JSON object")
    return value


def cart(radius: float, theta: float, phi: float) -> tuple[float, float, float]:
    return (radius * math.sin(theta) * math.cos(phi),
            radius * math.sin(theta) * math.sin(phi), radius * math.cos(theta))


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.hypot(math.hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2])


def complex_value(row: dict[str, str], re_key: str, im_key: str, label: str) -> complex:
    return complex(finite(row[re_key], f"{label}.{re_key}"), finite(row[im_key], f"{label}.{im_key}"))


def load_state(run_root: Path, plan_path: Path, state_path: Path) -> dict[str, Any]:
    state = read_json(state_path)
    if state.get("status") != "completed" or state.get("plan_sha256") != sha256(plan_path):
        raise ValueError("producer state is not a completed snapshot of the supplied exact run plan")
    by_id = {entry.get("id"): entry for entry in state.get("jobs", [])}
    if set(by_id) != set(JOB_IDS):
        raise ValueError("producer state job inventory differs from the exact ten-job Study 1 queue")
    for job_id in JOB_IDS:
        entry = by_id[job_id]
        if (entry.get("status") != "completed" or entry.get("returncode") != 0 or
                entry.get("failure_reasons") or entry.get("inputs_unchanged") is not True or
                not entry.get("cleanup", {}).get("verified_gone")):
            raise ValueError(f"producer job {job_id} lacks clean completion evidence")
        job_dir = run_root / job_id
        log = job_dir / entry["log"]
        if not log.is_file() or sha256(log) != entry.get("log_sha256"):
            raise ValueError(f"producer log hash mismatch for {job_id}")
        for relative in entry.get("outputs", []):
            path = job_dir / relative
            expected = entry.get("outputs_sha256", {}).get(relative)
            if not path.is_file() or not expected or sha256(path) != expected:
                raise ValueError(f"producer output hash mismatch for {job_id}/{relative}")
    return by_id


def job_data(run_root: Path, job_id: str) -> Path:
    return run_root / job_id / "data"


def validate_result(data: Path, backend: str, degree: int, nq: int, first: int, last: int,
                    selected_h: float, solver_tolerance: float | None = None,
                    audit: bool = False, epsilon: float | None = None,
                    frozen_weights: bool = False) -> dict[str, Any]:
    meta = read_json(data / "result.json")
    selected_h_actual = finite(meta.get("phase_d_selected_h", math.nan), f"{data}.phase_d_selected_h")
    time_shift_identity = finite(meta.get("time_shift_adjoint_identity", math.nan), f"{data}.time_shift_adjoint_identity")
    gravity_error = finite(meta.get("gravity_max_abs_nd_error", math.nan), f"{data}.gravity_max_abs_nd_error")
    measurement_norm = finite(meta.get("measurement_rate_norm", math.nan), f"{data}.measurement_rate_norm")
    if (meta.get("backend") != backend or meta.get("degree") != degree or meta.get("nq") != nq or
            meta.get("frequency_count") != last - first + 1 or meta.get("first_bin") != first or
            meta.get("last_bin") != last or selected_h_actual <= 0.0 or abs(selected_h_actual - selected_h) > 1e-12 or
            time_shift_identity <= 0.0 or abs(time_shift_identity - 1.0) > 1e-10 or
            gravity_error < 0.0 or gravity_error > 1e-10 or measurement_norm <= 0.0 or
            meta.get("fixed_probe_contact_terms_included") is not False):
        raise ValueError(f"{data} result metadata differs from the fixed Study 1 workload")
    if meta.get("fixed_probe_regular_part_only") != (meta.get("fixed_probe_count", 0) > 0):
        raise ValueError(f"{data} fixed-probe/contact metadata is inconsistent")
    if solver_tolerance is not None and abs(finite(meta.get("yspec_radial_solver_tolerance", math.nan), "yspec_radial_solver_tolerance") - solver_tolerance) > 1e-15:
        raise ValueError(f"{data} YSpec ODE tolerance differs from the planned control")
    if epsilon is not None and finite(meta.get("epsilon_perturbation", math.nan), "epsilon_perturbation") != epsilon:
        raise ValueError(f"{data} perturbation amplitude differs from the exact central pair")
    expected_weight_source = "frozen supplied baseline waveform adjoint source" if frozen_weights else "this mesh's own waveform-derived adjoint source"
    if meta.get("kernel_adjoint_weight_source") != expected_weight_source:
        raise ValueError(f"{data} adjoint-weight provenance differs from the planned functional")
    if audit:
        if meta.get("along_angles") != 91 or meta.get("across_angles") != 0:
            raise ValueError(f"{data} audit is not the frozen L3 91-angle full-volume geometry")
        finite(meta.get("integrated_profile_kernel", math.nan), f"{data}.integrated_profile_kernel")
        native_sides = meta.get("native_side_count")
        if not isinstance(native_sides, int) or native_sides <= 0 or len(rows(data / "nodes.csv")) != native_sides:
            raise ValueError(f"{data} audit side count differs from its explicit nodes table")
    else:
        if meta.get("along_angles") != 135 or meta.get("across_angles") != 19:
            raise ValueError(f"{data} section geometry differs from the matched L200 grid")
        native_sides = meta.get("native_side_count")
        if not isinstance(native_sides, int) or native_sides <= 0 or len(rows(data / "nodes.csv")) != native_sides:
            raise ValueError(f"{data} native side count differs from its explicit nodes table")
    if meta.get("fixed_probe_count", 0) not in (0, 168):
        raise ValueError(f"{data} fixed probe count is not the shared grid")
    return meta


def read_trace(path: Path) -> list[dict[str, str]]:
    result = rows(path)
    required = {"time_s", "velocity_m_s", "rate_m_s2", "window", "adjoint_source"}
    if not result or set(result[0]) != required:
        raise ValueError(f"{path} trace schema differs")
    for row in result:
        for key in required:
            finite(row[key], f"{path}.{key}")
    return result


def read_weights(data: Path) -> dict[int, complex]:
    result = rows(data / "fixed-probe-frequency-weights.csv")
    if len(result) != 35:
        raise ValueError(f"{data} own frequency-weight output must contain all 35 bins")
    found: dict[int, complex] = {}
    for row in result:
        bin_id = int(row["bin"])
        if bin_id in found:
            raise ValueError(f"{data} repeats frequency bin {bin_id}")
        weight = complex_value(row, "weight_re", "weight_im", "weight")
        force = complex_value(row, "force_re", "force_im", "force")
        p = complex_value(row, "P_re", "P_im", "P")
        close_complex(p, weight * force, 1e-13, f"{data} P at bin {bin_id}", abs(p) + abs(weight * force))
        found[bin_id] = p
    if set(found) != set(range(7, 42)):
        raise ValueError(f"{data} own P inventory is not exactly bins 7–41")
    return found


def fixed_table(data: Path, expected_probes: list[dict[str, str]], backend: str,
                source_xyz: tuple[float, float, float],
                receiver_xyz: tuple[float, float, float]) -> dict[str, dict[str, str | float]]:
    table = rows(data / "fixed-probe-kernel.csv")
    if len(table) != len(expected_probes):
        raise ValueError(f"{data} fixed integrated probe output count differs")
    expected = {row["probe_id"]: row for row in expected_probes}
    actual: dict[str, dict[str, str | float]] = {}
    for row in table:
        point_id = row["probe_id"]
        if point_id not in expected or point_id in actual:
            raise ValueError(f"{data} fixed probe ID is missing or duplicated: {point_id}")
        ref = expected[point_id]
        for key in ("layer", "material_side", "radial_trace"):
            if row[key] != ref[key]:
                raise ValueError(f"{data}/{point_id} changed requested {key}")
        for key in ("radius_m", "theta_rad", "phi_rad"):
            if abs(finite(row[key], key) - float(ref[key])) > 1e-9:
                raise ValueError(f"{data}/{point_id} changed physical {key}")
        point_xyz = cart(float(ref["radius_m"]), float(ref["theta_rad"]), float(ref["phi_rad"]))
        source_distance, receiver_distance = distance(point_xyz, source_xyz), distance(point_xyz, receiver_xyz)
        if (abs(finite(row["source_distance_3d_m"], "source_distance_3d_m") - source_distance) > 1e-5 or
                abs(finite(row["receiver_distance_3d_m"], "receiver_distance_3d_m") - receiver_distance) > 1e-5):
            raise ValueError(f"{data}/{point_id} stored source/receiver distance differs from Cartesian geometry")
        numeric_or_unexposed: str | float
        if row["kernel_vs_s_m-3"] == "not_exposed":
            numeric_or_unexposed = "not_exposed"
            if backend != "yspec" or row["radial_trace"] == "interior" or row["sample_semantics"] != "yspec_source_trace_not_exposed_by_kernel_sampling_api":
                raise ValueError(f"{data}/{point_id} has unexpected unavailable sample")
        else:
            numeric_or_unexposed = finite(row["kernel_vs_s_m-3"], f"{data}/{point_id}.kernel")
        if backend == "dspec":
            if row["phase_diagnostic_sem_element"] == "not_exposed" or row["resolved_element_trace"] in ("", "not_exposed"):
                raise ValueError(f"{data}/{point_id} DSpec element/trace resolution is absent")
            if row["radial_trace"] in ("inner", "outer") and row["resolved_element_trace"] in ("interior", "interior_resolved_to_lower_radius_element_endpoint"):
                raise ValueError(f"{data}/{point_id} one-sided trace was not resolved to its requested endpoint")
        elif row["phase_diagnostic_sem_element"] != "not_exposed" or row["resolved_element_trace"] != "not_exposed":
            raise ValueError(f"{data}/{point_id} YSpec exposes DSpec-only SEM metadata")
        if row["material_side"] == "fluid" and row["sample_semantics"] != "fluid_ln_vs_structural_zero_not_active_parameter":
            raise ValueError(f"{data}/{point_id} fluid lnVs row was not classified structurally")
        actual[point_id] = dict(row)
        actual[point_id]["numeric_kernel"] = numeric_or_unexposed
    if set(actual) != set(expected):
        raise ValueError(f"{data} fixed physical probe ID inventory differs")
    source_radius = math.hypot(math.hypot(source_xyz[0], source_xyz[1]), source_xyz[2])
    source_shell = [row for row in actual.values()
                    if row["layer"] == "9" and row["material_side"] == "solid" and
                    row["radial_trace"] in ("inner", "outer")]
    if len(source_shell) != 28:
        raise ValueError(f"{data} source-shell one-sided trace count differs from the frozen 28-row grid")
    source_coincident = [row for row in source_shell if row["location_class"] == "source_coincident"]
    off_source_shell = [row for row in source_shell if row["location_class"] == "source_shell_trace"]
    if len(source_coincident) != 2 or len(off_source_shell) != 26:
        raise ValueError(f"{data} does not distinguish point coincidence from off-angle source-shell traces")
    if any(abs(float(row["radius_m"]) - source_radius) > 1e-3 for row in source_shell):
        raise ValueError(f"{data} source-shell one-sided rows are not at the physical source radius")
    if backend == "yspec" and any(row["numeric_kernel"] != "not_exposed" for row in source_shell):
        raise ValueError(f"{data} must preserve all source-shell rows as unavailable one-sided traces")
    return actual


def validate_dspec_bin_products(data: Path, weights: dict[int, complex], integrated: dict[str, dict[str, str | float]]) -> dict[str, complex]:
    entries = rows(data / "fixed-probe-bin-kernel.csv")
    if len(entries) != 168 * 35:
        raise ValueError(f"{data} should retain one fixed probe record for every bin")
    sums = {probe_id: 0j for probe_id in integrated}
    sum_abs = {probe_id: 0.0 for probe_id in integrated}
    seen: set[tuple[str, int]] = set()
    for row in entries:
        probe_id, bin_id = row["probe_id"], int(row["bin"])
        if probe_id not in sums or bin_id not in weights or (probe_id, bin_id) in seen:
            raise ValueError(f"{data} fixed-probe per-bin key inventory is invalid")
        seen.add((probe_id, bin_id))
        p = complex(finite(row["P_re"], "P_re"), finite(row["P_im"], "P_im"))
        close_complex(p, weights[bin_id], 1e-13, f"{data} per-probe P", abs(p) + abs(weights[bin_id]))
        computed = row["unit_adjoint_computed"] == "1"
        if row["unit_adjoint_computed"] not in ("0", "1"):
            raise ValueError(f"{data} invalid zero-skip marker")
        if not computed:
            if bin_id not in (7, 8, 41) or p != 0j or row["unit_complex_product_re"] != "not_computed" or row["unit_complex_product_im"] != "not_computed":
                raise ValueError(f"{data} uncomputed unit product is not an exact taper-zero bin")
            unit_product = 0j
        else:
            if bin_id in (7, 8, 41):
                raise ValueError(f"{data} exact taper-zero bin was unexpectedly solved")
            unit_product = complex(finite(row["unit_complex_product_re"], "unit_re"),
                                   finite(row["unit_complex_product_im"], "unit_im"))
        weighted = complex(finite(row["weighted_complex_product_re"], "weighted_re"),
                           finite(row["weighted_complex_product_im"], "weighted_im"))
        product = p * unit_product
        close_complex(weighted, product, 1e-12, f"{data} weighted product", abs(weighted) + abs(product))
        close_scaled(finite(row["weighted_kernel_vs"], "weighted_kernel"), weighted.real,
                     1e-13, f"{data} real kernel projection")
        sums[probe_id] += weighted
        sum_abs[probe_id] += abs(weighted)
    for probe_id, value in sums.items():
        target = integrated[probe_id]["numeric_kernel"]
        if not isinstance(target, float):
            raise ValueError(f"{data}/{probe_id} integrated kernel is not numeric")
        if sum_abs[probe_id] == 0.0:
            close_scaled(value.real, target, 2e-12, f"{data}/{probe_id} exact-zero integrated kernel")
        elif abs(value.real - target) > 2e-12 * sum_abs[probe_id]:
            raise ValueError(f"{data}/{probe_id} integrated kernel does not equal its ordered bin sum")
    return sums


def validate_fields(data: Path, weights: dict[int, complex], expected_probes: list[dict[str, str]]) -> None:
    expected_ids = {row["probe_id"] for row in expected_probes}
    for bin_id in (9, 25):
        entries = rows(data / f"fixed-probe-unit-fields-bin{bin_id}.csv")
        if len(entries) != 168 or {row["probe_id"] for row in entries} != expected_ids:
            raise ValueError(f"{data} selected-bin field packet inventory differs")
        p = weights[bin_id]
        for row in entries:
            row_p = complex(finite(row["P_re"], "field.P_re"), finite(row["P_im"], "field.P_im"))
            close_complex(row_p, p, 1e-13, f"{data} bin {bin_id} fixed field packet P", abs(row_p) + abs(p))
            for c in range(6):
                unit = complex_value(row, f"unit_adjoint_dual_strain_{c}_re", f"unit_adjoint_dual_strain_{c}_im", "unit_adjoint")
                actual = complex_value(row, f"measurement_weighted_adjoint_dual_strain_{c}_re", f"measurement_weighted_adjoint_dual_strain_{c}_im", "weighted_adjoint")
                expected = p.conjugate() * unit
                close_complex(actual, expected, 1e-12,
                              f"{data} bin {bin_id} probe {row['probe_id']} adjoint strain weight",
                              abs(actual) + abs(expected))
                complex_value(row, f"forward_raw_strain_{c}_re", f"forward_raw_strain_{c}_im", "forward")


def _relative_l2(reference: list[float], candidate: list[float], label: str) -> float | None:
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(candidate, reference)))
    denominator = math.sqrt(sum(a * a for a in reference))
    if denominator == 0.0:
        if numerator == 0.0:
            return 0.0
        return None
    return finite(numerator / denominator, label)


def probe_control_summary(reference: dict[str, dict[str, str | float]],
                          candidate: dict[str, dict[str, str | float]],
                          probe_rows: list[dict[str, str]]) -> dict[str, Any]:
    ids = sorted(set(reference) & set(candidate))
    usable = [(pid, reference[pid]["numeric_kernel"], candidate[pid]["numeric_kernel"])
              for pid in ids if isinstance(reference[pid]["numeric_kernel"], float) and
              isinstance(candidate[pid]["numeric_kernel"], float)]
    ref_values = [float(item[1]) for item in usable]
    candidate_values = [float(item[2]) for item in usable]
    absolute = [(pid, abs(float(other) - float(ref))) for pid, ref, other in usable]
    max_item = max(absolute, key=lambda item: item[1]) if absolute else None
    expected = {row["probe_id"]: row for row in probe_rows}
    source_shell_ids = [pid for pid in ids if expected[pid]["layer"] == "9" and
                        expected[pid]["radial_trace"] in ("inner", "outer")]
    zero_pairs = sum(float(ref) == 0.0 and float(other) == 0.0 for _, ref, other in usable)
    maximum_relative = 0.0
    for _, ref, other in usable:
        magnitude = max(abs(float(ref)), abs(float(other)))
        relative = abs(float(other) - float(ref)) / magnitude if magnitude else 0.0
        maximum_relative = max(maximum_relative, relative)
    return {
        "common_numeric_probe_count": len(usable),
        "unmatched_source_shell_one_sided_count": len(source_shell_ids),
        "relative_l2_using_reference_kernel_norm": _relative_l2(ref_values, candidate_values, "fixed probe relative L2"),
        "reference_kernel_l2_norm": math.sqrt(sum(value * value for value in ref_values)),
        "candidate_kernel_l2_norm": math.sqrt(sum(value * value for value in candidate_values)),
        "max_absolute_difference": max_item[1] if max_item else None,
        "max_absolute_difference_probe_id": max_item[0] if max_item else None,
        "max_pointwise_relative_difference": maximum_relative,
        "exact_zero_pair_count": zero_pairs,
    }


def trace_control_summary(reference_data: Path, candidate_data: Path) -> dict[str, Any]:
    reference, candidate = read_trace(reference_data / "trace.csv"), read_trace(candidate_data / "trace.csv")
    if len(reference) != len(candidate):
        raise ValueError("control and reference traces differ in sample count")
    if any(a["time_s"] != b["time_s"] for a, b in zip(reference, candidate)):
        raise ValueError("control and reference traces differ in exact sample times")
    fields = ("velocity_m_s", "rate_m_s2", "adjoint_source")
    summary = {field: _relative_l2([float(row[field]) for row in reference],
                                   [float(row[field]) for row in candidate],
                                   f"trace {field} relative L2") for field in fields}
    summary["sample_count"] = len(reference)
    return summary


def native_section_difference(reference_data: Path, candidate_data: Path,
                              source_radius: float) -> dict[str, Any]:
    reference_nodes = rows(reference_data / "nodes.csv")
    candidate_nodes = rows(candidate_data / "nodes.csv")
    if reference_nodes != candidate_nodes:
        raise ValueError("native DSpec/YSpec comparison does not share identical element/node/radius coordinates")
    sections: dict[str, Any] = {}
    for section in ("along", "across"):
        ref_angles_path, cand_angles_path = reference_data / f"{section}-angles.csv", candidate_data / f"{section}-angles.csv"
        ref_values_path, cand_values_path = reference_data / f"{section}.bin", candidate_data / f"{section}.bin"
        if not ref_angles_path.exists() and not cand_angles_path.exists():
            continue
        if not ref_angles_path.is_file() or not cand_angles_path.is_file() or not ref_values_path.is_file() or not cand_values_path.is_file():
            raise ValueError(f"native DSpec/YSpec {section} section inventory differs")
        ref_angles, cand_angles = rows(ref_angles_path), rows(cand_angles_path)
        if ref_angles != cand_angles:
            raise ValueError(f"native DSpec/YSpec {section} angular coordinates differ")
        ref_raw, cand_raw = ref_values_path.read_bytes(), cand_values_path.read_bytes()
        if len(ref_raw) % 8 or len(cand_raw) % 8 or len(ref_raw) != len(cand_raw):
            raise ValueError(f"native DSpec/YSpec {section} arrays have different byte sizes")
        ref_values, cand_values = struct.unpack(f"{len(ref_raw)//8}d", ref_raw), struct.unpack(f"{len(cand_raw)//8}d", cand_raw)
        if len(ref_values) != len(reference_nodes) * len(ref_angles):
            raise ValueError(f"native DSpec/YSpec {section} array dimensions do not match the shared coordinates")
        source_indices, off_source_indices = [], []
        for n, node in enumerate(reference_nodes):
            target = source_indices if abs(float(node["radius_m"]) - source_radius) <= 1e-3 else off_source_indices
            target.extend(range(n * len(ref_angles), (n + 1) * len(ref_angles)))
        def stats(indices: list[int]) -> dict[str, Any]:
            finite_indices = [i for i in indices if math.isfinite(ref_values[i]) and math.isfinite(cand_values[i])]
            diffs = [cand_values[i] - ref_values[i] for i in finite_indices]
            denominator = math.sqrt(sum(ref_values[i] ** 2 for i in finite_indices))
            numerator = math.sqrt(sum(diff * diff for diff in diffs))
            return {"sample_count": len(indices), "finite_pair_count": len(finite_indices),
                    "nonfinite_pair_count": len(indices) - len(finite_indices), "reference_l2_norm": denominator,
                    "candidate_l2_norm": math.sqrt(sum(cand_values[i] ** 2 for i in finite_indices)),
                    "difference_l2_norm": numerator,
                    "relative_l2_using_reference_norm": numerator / denominator if denominator else (0.0 if numerator == 0.0 else None),
                    "max_absolute_difference": max(map(abs, diffs), default=None),
                    "max_abs_reference_value": max((abs(ref_values[i]) for i in finite_indices), default=None),
                    "max_abs_candidate_value": max((abs(cand_values[i]) for i in finite_indices), default=None)}
        sections[section] = {"raw_unmasked_all_samples": stats(list(range(len(ref_values)))),
                             "source_radius_unmatched_one_sided_cohort": stats(source_indices),
                             "off_source_radius_common_geometry": stats(off_source_indices)}
    return sections


def section_diagnostics(data: Path, source: tuple[float, float, float], receiver: tuple[float, float, float]) -> dict[str, Any]:
    nodes = rows(data / "nodes.csv")
    sections: dict[str, Any] = {}
    for section in ("along", "across"):
        angle_path = data / f"{section}-angles.csv"
        values_path = data / f"{section}.bin"
        if not angle_path.is_file() or not values_path.is_file():
            if section == "across" and not (data / "across-angles.csv").exists() and not (data / "across.bin").exists():
                continue
            raise ValueError(f"{data} missing one of {section} angle/value files")
        angle_rows = rows(angle_path)
        raw = values_path.read_bytes()
        values = struct.unpack(f"{len(raw)//8}d", raw)
        if len(raw) % 8 or len(values) != len(nodes) * len(angle_rows):
            raise ValueError(f"{data}/{section} sampled section dimensions differ")
        point_diagnostics = []
        for n, node in enumerate(nodes):
            radius = finite(node["radius_m"], "node radius")
            for a, angle in enumerate(angle_rows):
                value = values[n * len(angle_rows) + a]
                theta, phi = finite(angle["colatitude_rad"], "angle theta"), finite(angle["longitude_rad"], "angle phi")
                xyz = cart(radius, theta, phi)
                ds, dr = distance(xyz, source), distance(xyz, receiver)
                point_diagnostics.append((value, ds, dr, n, a, radius, theta, phi))
        finite_values = [item[0] for item in point_diagnostics if math.isfinite(item[0])]
        if len(finite_values) != len(values):
            nonfinite = len(values) - len(finite_values)
        else:
            nonfinite = 0
        summary: dict[str, Any] = {
            "sample_count": len(values), "nonfinite_count": nonfinite,
            "raw_min": min(finite_values) if finite_values else None,
            "raw_max": max(finite_values) if finite_values else None,
            "raw_max_abs": max(map(abs, finite_values)) if finite_values else None,
            "exclusion_regions_are_section_subsets_not_volume_norms": True,
            "exclusion_subsets": {},
        }
        def describe(item):
            node = nodes[item[3]]
            angle = angle_rows[item[4]]
            return {"value": item[0], "source_distance_m": item[1], "receiver_distance_m": item[2],
                    "node_row": item[3], "angle_row": item[4], "element": node["element"],
                    "node": node["node"], "layer": node["layer"], "radius_m": item[5],
                    "solid": node["solid"], "theta_rad": item[6], "phi_rad": item[7]}
        finite_points = [item for item in point_diagnostics if math.isfinite(item[0])]
        if finite_points:
            min_point = min(finite_points, key=lambda item: item[0])
            max_point = max(finite_points, key=lambda item: item[0])
            abs_point = max(finite_points, key=lambda item: abs(item[0]))
            summary["raw_extreme_locations"] = {"min": describe(min_point), "max": describe(max_point),
                                                 "max_abs": describe(abs_point)}
        source_radius = math.hypot(math.hypot(source[0], source[1]), source[2])
        source_shell_samples = [item for item in point_diagnostics if abs(item[5] - source_radius) <= 1e-3]
        finite_source_shell = [item[0] for item in source_shell_samples if math.isfinite(item[0])]
        summary["source_radius_native_samples"] = {
            "sample_count": len(source_shell_samples),
            "nonfinite_count": len(source_shell_samples) - len(finite_source_shell),
            "raw_min": min(finite_source_shell) if finite_source_shell else None,
            "raw_max": max(finite_source_shell) if finite_source_shell else None,
            "raw_max_abs": max(map(abs, finite_source_shell)) if finite_source_shell else None,
            "nearest_3d_sample_distance_m": min((item[1] for item in source_shell_samples), default=None),
            "interpretation": "raw native samples retained; YSpec does not expose the artificial source-shell inner/outer SEM trace, so these are not matched one-sided DSpec/YSpec comparisons",
        }
        for radius_m in ANGULAR_LIMIT_RADIAL_M:
            included = [item[0] for item in point_diagnostics
                        if math.isfinite(item[0]) and item[1] >= radius_m and item[2] >= radius_m]
            summary["exclusion_subsets"][str(int(radius_m / 1000)) + "km"] = {
                "point_count": len(included),
                "min": min(included) if included else None,
                "max": max(included) if included else None,
                "max_abs": max(map(abs, included)) if included else None,
            }
        def nearest_samples(which: int):
            if not finite_points:
                return None
            minimum = min(item[which] for item in finite_points)
            tied = [item for item in finite_points if abs(item[which] - minimum) <= 1e-6]
            return {"distance_m": minimum, "tie_tolerance_m": 1e-6,
                    "all_nearest_ties": [describe(item) for item in tied]}
        summary["nearest_source_samples"] = nearest_samples(1)
        summary["nearest_receiver_samples"] = nearest_samples(2)
        sections[section] = summary
    return sections


def validate_samples(run_root: Path, probe_rows: list[dict[str, str]], probe_receipt: dict[str, Any]) -> dict[str, Any]:
    manifests: dict[str, Any] = {}
    results: dict[str, dict[str, Any]] = {}
    source_xyz = cart(probe_receipt["source"]["radius_m"], probe_receipt["source"]["geocentric_theta_rad"], probe_receipt["source"]["geocentric_phi_rad"])
    receiver_xyz = cart(probe_receipt["receiver"]["radius_m"], probe_receipt["receiver"]["geocentric_theta_rad"], probe_receipt["receiver"]["geocentric_phi_rad"])
    for job_id, backend, degree, nq, h, tol, audit in (
        (BASE_D, "dspec", 200, 6, 0.02035382429229197, None, False),
        (BASE_Y, "yspec", 200, 6, 0.02035382429229197, 1e-5, False),
        (FINE_D, "dspec", 200, 6, 0.01809228825981428, None, False),
        (TIGHT_Y, "yspec", 200, 6, 0.02035382429229197, 1e-7, False),
        (AUDIT_D, "dspec", 3, 8, 0.05, None, True),
        (AUDIT_Y, "yspec", 3, 8, 0.05, 1e-5, True),
    ):
        data = job_data(run_root, job_id)
        meta = validate_result(data, backend, degree, nq, 7, 41, h,
                               solver_tolerance=tol, audit=audit)
        results[job_id] = meta
        if job_id in (BASE_D, BASE_Y, FINE_D, TIGHT_Y, AUDIT_D, AUDIT_Y):
            expected_count = 168
            if meta.get("fixed_probe_count") != expected_count:
                raise ValueError(f"{job_id} did not consume all 168 fixed probes")
            probes = fixed_table(data, probe_rows, backend, source_xyz, receiver_xyz)
            if backend == "dspec":
                storage_bytes = meta.get("fixed_probe_storage_bytes")
                if not isinstance(storage_bytes, int) or storage_bytes <= 0 or storage_bytes > (32 << 20):
                    raise ValueError(f"{job_id} DSpec fixed-probe storage violates its positive 32 MiB bound")
            weight_rows = read_weights(data)
            if backend == "dspec":
                validate_dspec_bin_products(data, weight_rows, probes)
                if degree == 200:
                    validate_fields(data, weight_rows, probe_rows)
            else:
                source_not_exposed = sum(v["numeric_kernel"] == "not_exposed" for v in probes.values())
                if source_not_exposed != 28:
                    raise ValueError(f"{job_id} must mark exactly the 28 YSpec source-shell traces unavailable")
            manifests[job_id] = {"fixed_probe_count": len(probes), "frequency_weight_count": len(weight_rows)}
    for baseline_job in (BASE_D, BASE_Y):
        if results[baseline_job].get("native_side_count") != 208:
            raise ValueError(f"{baseline_job} baseline must retain the 208 native one-sided radial samples")
    if probe_receipt.get("probe_count") != 168 or probe_receipt.get("probe_csv_sha256") != sha256(Path(probe_receipt["_csv_path"])):
        raise ValueError("fixed probe provenance receipt does not bind the staged CSV")

    ds = fixed_table(job_data(run_root, BASE_D), probe_rows, "dspec", source_xyz, receiver_xyz)
    ys = fixed_table(job_data(run_root, BASE_Y), probe_rows, "yspec", source_xyz, receiver_xyz)
    fine_d = fixed_table(job_data(run_root, FINE_D), probe_rows, "dspec", source_xyz, receiver_xyz)
    tight_y = fixed_table(job_data(run_root, TIGHT_Y), probe_rows, "yspec", source_xyz, receiver_xyz)
    p_ds = read_weights(job_data(run_root, BASE_D))
    p_ys = read_weights(job_data(run_root, BASE_Y))
    d_unit = {row["probe_id"]: {} for row in probe_rows}
    for row in rows(job_data(run_root, BASE_D) / "fixed-probe-bin-kernel.csv"):
        pid, k = row["probe_id"], int(row["bin"])
        if row["unit_complex_product_re"] != "not_computed":
            d_unit[pid][k] = complex(float(row["unit_complex_product_re"]), float(row["unit_complex_product_im"]))
    attribution = []
    abs_diffs = []
    crossweights = {}
    expected_by_id = {row["probe_id"]: row for row in probe_rows}
    for pid in sorted(ds):
        dv, yv = ds[pid]["numeric_kernel"], ys[pid]["numeric_kernel"]
        if not isinstance(dv, float) or not isinstance(yv, float):
            continue
        cross = sum((p_ys[k] * d_unit[pid][k]).real for k in range(7, 42) if k in d_unit[pid])
        crossweights[pid] = cross
        abs_diffs.append(abs(dv - yv))
        probe = expected_by_id[pid]
        xyz = cart(float(probe["radius_m"]), float(probe["theta_rad"]), float(probe["phi_rad"]))
        source_distance, receiver_distance = distance(xyz, source_xyz), distance(xyz, receiver_xyz)
        attribution.append({"probe_id": pid, "radius_m": float(probe["radius_m"]),
                            "source_distance_m": source_distance, "receiver_distance_m": receiver_distance,
                            "material_side": probe["material_side"], "radial_trace": probe["radial_trace"],
                            "dspec_own_weight_kernel": dv,
                            "yspec_own_weight_kernel": yv,
                            "dspec_unit_products_with_yspec_own_P_counterfactual": cross,
                            "dspec_minus_yspec_own_kernel": dv - yv})
    if not attribution:
        raise ValueError("no common exposed physical probes for DSpec/YSpec comparison")
    p_fine_d = read_weights(job_data(run_root, FINE_D))
    p_tight_y = read_weights(job_data(run_root, TIGHT_Y))
    def weight_change(reference: dict[int, complex], candidate: dict[int, complex]) -> dict[str, Any]:
        bins = sorted(set(reference) & set(candidate))
        ref_norm = math.sqrt(sum(abs(reference[k]) ** 2 for k in bins))
        delta_norm = math.sqrt(sum(abs(candidate[k] - reference[k]) ** 2 for k in bins))
        maximum = max(bins, key=lambda k: abs(candidate[k] - reference[k]))
        return {"bin_count": len(bins), "reference_P_l2_norm": ref_norm,
                "candidate_P_l2_norm": math.sqrt(sum(abs(candidate[k]) ** 2 for k in bins)),
                "P_difference_l2_norm": delta_norm,
                "relative_l2_using_reference_P_norm": delta_norm / ref_norm if ref_norm else (0.0 if delta_norm == 0.0 else None),
                "max_abs_P_difference": abs(candidate[maximum] - reference[maximum]),
                "max_difference_bin": maximum,
                "exact_zero_bins_in_reference": [k for k in bins if reference[k] == 0j],
                "exact_zero_bins_in_candidate": [k for k in bins if candidate[k] == 0j]}
    p_diff = [abs(p_ds[k] - p_ys[k]) for k in range(7, 42)]
    output = {
        "schema": 1,
        "study": "Study 1 discrepancy and point-probe diagnostics",
        "evidence_limits": {
            "section_exclusions_are_not_volume_norms": True,
            "fixed_probe_source_receiver_samples_are_diagnostics_not_contact_values": True,
            "native_source_radius_samples_are_unmatched_one_sided_traces": "Raw extrema are retained, but YSpec does not expose the artificial source-shell inner/outer SEM trace; no matched one-sided DSpec/YSpec conclusion is drawn at that radius.",
            "unavailable_fixed_source_shell_rows": 28,
            "common_fixed_point_comparison_uses_only_backend_exposed_numeric_values": True,
            "yspec_per_frequency_unit_kernel_fields_available": False,
            "crossweight_control": "DSpec unit products weighted by YSpec own P; not a YSpec unit-field comparison",
            "smooth_action": "true 3D integrated smooth radial lnVs profile compared against DSpec full-forward centered linearized-functional differences",
            "independent_delta_H_action": "not performed",
        },
        "probe_csv_sha256": sha256(Path(probe_receipt["_csv_path"])),
        "probe_receipt_sha256": probe_receipt["_receipt_sha256"],
        "probe_count": len(probe_rows),
        "angle_count": probe_receipt["unique_angle_count"],
        "result_metadata": results,
        "probe_validation": manifests,
        "pointwise_native_kernel_comparison": {
            "common_exposed_points": len(attribution),
            "max_abs_own_kernel_difference": max(abs_diffs),
            "mean_abs_own_kernel_difference": sum(abs_diffs) / len(abs_diffs),
            "max_abs_own_P_difference_by_bin": max(p_diff),
            "fixed_point_attribution": attribution,
            "unmatched_source_shell_trace_samples": [
                {"probe_id": pid, "radial_trace": expected_by_id[pid]["radial_trace"],
                 "DSpec_raw_kernel_sample": ds[pid]["numeric_kernel"],
                 "YSpec_value": ys[pid]["numeric_kernel"],
                 "location_class": ds[pid]["location_class"]}
                for pid in sorted(ds)
                if expected_by_id[pid]["layer"] == "9" and
                expected_by_id[pid]["radial_trace"] in ("inner", "outer")
            ],
        },
        "l200_control_effects": {
            "dspec_radial_mesh_refinement": {
                "fixed_probe_kernels": probe_control_summary(ds, fine_d, probe_rows),
                "waveform_traces": trace_control_summary(job_data(run_root, BASE_D), job_data(run_root, FINE_D)),
                "measurement_rate_norm_baseline": results[BASE_D]["measurement_rate_norm"],
                "measurement_rate_norm_refined": results[FINE_D]["measurement_rate_norm"],
                "own_frequency_P_change": weight_change(p_ds, p_fine_d),
                "section_extrema_are_separately_reported": True,
            },
            "yspec_ode_tolerance_1e-5_to_1e-7": {
                "fixed_probe_kernels": probe_control_summary(ys, tight_y, probe_rows),
                "waveform_traces": trace_control_summary(job_data(run_root, BASE_Y), job_data(run_root, TIGHT_Y)),
                "measurement_rate_norm_baseline": results[BASE_Y]["measurement_rate_norm"],
                "measurement_rate_norm_tight": results[TIGHT_Y]["measurement_rate_norm"],
                "own_frequency_P_change": weight_change(p_ys, p_tight_y),
                "section_extrema_are_separately_reported": True,
            },
        },
        "native_l200_dspec_yspec_array_comparison": native_section_difference(
            job_data(run_root, BASE_D), job_data(run_root, BASE_Y),
            math.hypot(math.hypot(source_xyz[0], source_xyz[1]), source_xyz[2])),
        "section_diagnostics": {},
        "full_3d_profile_actions": {
            "dspec_l3_integral": results[AUDIT_D]["integrated_profile_kernel"],
            "yspec_l3_integral": results[AUDIT_Y]["integrated_profile_kernel"],
            "fd_pairs": [],
        },
    }
    source_geo = probe_receipt["source"]
    receiver_geo = probe_receipt["receiver"]
    for job_id in (BASE_D, BASE_Y, FINE_D, TIGHT_Y):
        output["section_diagnostics"][job_id] = section_diagnostics(job_data(run_root, job_id), source_xyz, receiver_xyz)

    base_trace = read_trace(job_data(run_root, AUDIT_D) / "trace.csv")
    base_meta = results[AUDIT_D]
    dt_values = [float(base_trace[i + 1]["time_s"]) - float(base_trace[i]["time_s"]) for i in range(len(base_trace) - 1)]
    if not dt_values or max(dt_values) - min(dt_values) > 1e-12:
        raise ValueError("L3 baseline trace has a nonuniform sample interval")
    dt = dt_values[0]
    base_adjoint = [float(row["adjoint_source"]) for row in base_trace]
    base_velocity = [float(row["velocity_m_s"]) for row in base_trace]
    base_velocity_norm = math.sqrt(sum(value * value for value in base_velocity))
    base_adjoint_norm = math.sqrt(sum(value * value for value in base_adjoint))
    if base_velocity_norm <= 0.0 or base_adjoint_norm <= 0.0:
        raise ValueError("L3 baseline waveform and adjoint source must both have positive norms")
    support_min, support_max = 5.3e6, 5.8e6
    source_radius, receiver_radius = source_geo["radius_m"], receiver_geo["radius_m"]
    support_clearance = min(source_radius - support_max, receiver_radius - support_max)
    if support_clearance <= 0:
        raise ValueError("the smooth lnVs bump touches source or receiver radii")
    fd_failures: list[str] = []
    for plus_id, minus_id, epsilon in FD_PAIRS:
        plus_data, minus_data = job_data(run_root, plus_id), job_data(run_root, minus_id)
        plus_meta = validate_result(plus_data, "dspec", 3, 8, 7, 41, 0.05, audit=True,
                                    epsilon=epsilon, frozen_weights=False)
        minus_meta = validate_result(minus_data, "dspec", 3, 8, 7, 41, 0.05, audit=True,
                                     epsilon=-epsilon, frozen_weights=False)
        plus_trace, minus_trace = read_trace(plus_data / "trace.csv"), read_trace(minus_data / "trace.csv")
        if len(plus_trace) != len(base_trace) or len(minus_trace) != len(base_trace):
            raise ValueError("L3 baseline and +/- waveform traces have different lengths")
        if any(base["window"] != plus["window"] or base["window"] != minus["window"]
               for base, plus, minus in zip(base_trace, plus_trace, minus_trace)):
            raise ValueError("L3 baseline and +/- cases do not use the exact same measurement window")
        mesh_evidence = {}
        for relative in ("mesh-elements.csv", "nodes.csv", "frequency-grid.csv", "parameters.txt"):
            base_path = job_data(run_root, AUDIT_D) / relative
            plus_path, minus_path = plus_data / relative, minus_data / relative
            if not base_path.is_file() or not plus_path.is_file() or not minus_path.is_file():
                raise ValueError(f"L3 centered difference lacks required matched-space file {relative}")
            base_hash, plus_hash, minus_hash = sha256(base_path), sha256(plus_path), sha256(minus_path)
            if base_hash != plus_hash or base_hash != minus_hash:
                raise ValueError(f"L3 centered pair does not share exact mesh/grid/parameter file {relative}")
            mesh_evidence[relative] = base_hash
        delta = 0.0
        for i, (base, plus, minus) in enumerate(zip(base_trace, plus_trace, minus_trace)):
            t = i * dt
            if (float(base["time_s"]) != t or float(plus["time_s"]) != t or float(minus["time_s"]) != t):
                raise ValueError("L3 finite-difference traces do not share the exact uniform time grid")
            delta += dt * base_adjoint[i] * (float(plus["velocity_m_s"]) - float(minus["velocity_m_s"])) / (2.0 * epsilon)
        integral = float(base_meta["integrated_profile_kernel"])
        scale = max(abs(delta), abs(integral))
        relative = abs(delta - integral) / scale if scale else (0.0 if delta == integral == 0.0 else math.inf)
        relative = finite(relative, f"FD relative difference at epsilon {epsilon}")
        if scale == 0.0:
            fd_failures.append(f"FD and kernel action are both exactly zero at epsilon {epsilon}; nonzero action is required")
        if relative > 2e-4:
            fd_failures.append(f"FD/kernel action difference exceeds 2e-4 at epsilon {epsilon}")
        output["full_3d_profile_actions"]["fd_pairs"].append({
            "epsilon_ln_vs_peak": epsilon,
            "same_frozen_baseline_functional": True,
            "baseline_trace_job": AUDIT_D,
            "exact_shared_fe_and_frequency_space": True,
            "shared_space_file_sha256": mesh_evidence,
            "plus_job": plus_id,
            "minus_job": minus_id,
            "centered_full_forward_linearized_functional_derivative": delta,
            "dspec_true_3d_kernel_quadrature": integral,
            "absolute_difference": abs(delta - integral),
            "relative_difference_using_max_magnitude_scale": relative,
            "radial_bump_support_m": [support_min, support_max],
            "source_radius_m": source_radius,
            "receiver_radius_m": receiver_radius,
            "minimum_source_or_receiver_radius_clearance_m": support_clearance,
            "source_lift_or_receiver_projection_changes_expected": False,
            "baseline_waveform_l2_norm": base_velocity_norm,
            "baseline_adjoint_source_l2_norm": base_adjoint_norm,
            "receiver_functional_side": "unchanged solid side at the existing solid-fluid interface",
            "fd_dot_product_uses_baseline_adjoint_source": True,
            "plus_minus_runs_own_adjoint_sources_ignored_by_fd_dot_product": True,
        })
    fds = output["full_3d_profile_actions"]["fd_pairs"]
    output["full_3d_profile_actions"]["fd_step_pair_relative_change"] = (
        abs(fds[0]["centered_full_forward_linearized_functional_derivative"] -
            fds[1]["centered_full_forward_linearized_functional_derivative"]) /
        max(abs(fds[0]["centered_full_forward_linearized_functional_derivative"]),
            abs(fds[1]["centered_full_forward_linearized_functional_derivative"]))
        if max(abs(fds[0]["centered_full_forward_linearized_functional_derivative"]),
               abs(fds[1]["centered_full_forward_linearized_functional_derivative"])) else 0.0)
    output["full_3d_profile_actions"]["fd_step_pair_relative_change"] = finite(
        output["full_3d_profile_actions"]["fd_step_pair_relative_change"], "FD step-pair stability")
    if output["full_3d_profile_actions"]["fd_step_pair_relative_change"] > 2e-4:
        fd_failures.append("FD step-pair stability exceeds 2e-4")
    output["science_acceptance"] = {"status": "FAIL" if fd_failures else "PASS",
                                    "fd_relative_limit": 2e-4,
                                    "fd_failures": fd_failures}
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--run-plan", type=Path, required=True)
    parser.add_argument("--supervisor-state", type=Path, required=True)
    parser.add_argument("--probe-csv", type=Path, required=True)
    parser.add_argument("--probe-receipt", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--tight-case", type=Path, required=True)
    parser.add_argument("--case-model", type=Path, required=True)
    parser.add_argument("--mesh-receipt", type=Path, required=True)
    parser.add_argument("--path-reference", type=Path, required=True)
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    probe_rows = rows(args.probe_csv)
    receipt = read_json(args.probe_receipt)
    receipt["_csv_path"] = str(args.probe_csv)
    receipt["_receipt_sha256"] = sha256(args.probe_receipt)
    if sha256(args.probe_csv) != receipt.get("probe_csv_sha256"):
        raise ValueError("staged fixed-probe CSV differs from its hash-bound receipt")
    provenance = ((args.case, "case_sha256"), (args.mesh_receipt, "mesh_receipt_sha256"),
                  (args.path_reference, "path_reference_sha256"), (args.generator, "generator_sha256"))
    for path, key in provenance:
        if sha256(path) != receipt.get(key):
            raise ValueError(f"probe receipt provenance mismatch for {key}")
    base_case, tight_case = read_json(args.case), read_json(args.tight_case)
    expected_tight = copy.deepcopy(base_case)
    expected_tight.pop("solver", None)
    normalized_tight = copy.deepcopy(tight_case)
    if normalized_tight.get("solver") != {"tolerance": 1e-7}:
        raise ValueError("tight YSpec case does not set exactly the planned 1e-7 solver tolerance")
    normalized_tight.pop("solver")
    normalized_tight["model"] = expected_tight["model"]
    if normalized_tight != expected_tight:
        raise ValueError("tight YSpec case changed a scientific field besides solver tolerance/model path")
    expected_model_path = "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/scientific-inputs/dspec-science-study1-v1/candidate/dspec/case/prem.200"
    if str(args.case_model) != expected_model_path or tight_case.get("model") != expected_model_path:
        raise ValueError("tight case model reference does not target the exact staged frozen prem.200 deck")
    if sha256(args.case_model) != "4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea":
        raise ValueError("staged prem.200 deck differs from the hash-bound Study 1 build input")
    if len(probe_rows) != 168 or receipt.get("probe_count") != len(probe_rows):
        raise ValueError("staged common probe grid does not have the planned 168 rows")
    load_state(args.run_root, args.run_plan, args.supervisor_state)
    report = validate_samples(args.run_root, probe_rows, receipt)
    report["tight_case_sha256"] = sha256(args.tight_case)
    report["case_model_sha256"] = sha256(args.case_model)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    status = report["science_acceptance"]["status"]
    print(json.dumps({"status": status, "out": str(args.out),
                      "probes": report["probe_count"],
                      "common_exposed_points": report["pointwise_native_kernel_comparison"]["common_exposed_points"],
                      "fd_pairs": len(report["full_3d_profile_actions"]["fd_pairs"])}, sort_keys=True))
    if status != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
