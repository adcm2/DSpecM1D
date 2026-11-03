#!/usr/bin/env python3
"""Data-only six-job qualification plus a seventh analytic-reference control.

The frozen six-job checker is imported from tools/qualification_checks.py and
continues to own the original exact/native/profile, full-grid, radial, finite,
energy, and descriptive cross-code checks. Historical v4-v6 differences are
reported rather than hidden; analytic-control output must independently match
both v6 DSpec variants at the unchanged 2e-13 channel/column-scaled limit.
This program launches no scientific executable.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve()
CHECKER_PATH = HERE.parents[1] / "tools" / "qualification_checks.py"
SPEC = importlib.util.spec_from_file_location("frozen_qualification_checks", CHECKER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import frozen checker at {CHECKER_PATH}")
CHECKER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CHECKER
SPEC.loader.exec_module(CHECKER)

LIMIT = CHECKER.PARITY_LIMIT
CONTROL = "dsv4-analytic-control-l3"
ORIGINAL_PAIRS = (
    ("dsv4-native-l3", "dsv6-native-l3"),
    ("dsv4-profiled-l3", "dsv6-profiled-l3"),
)
DOMAIN_FIELDS = ("selected_radial_trace_count", "available_radial_trace_sides", "n_theta",
                 "n_phi", "full_target_point_count", "actual_sampled_point_count",
                 "binary_layout", "kernel_is_full_earth_volume_integral", "source_contact_enabled")


class QualificationError(Exception):
    pass


def _table_locations(a: np.ndarray, b: np.ndarray, columns: list[str]) -> dict:
    output = {}
    for index, name in enumerate(columns):
        row = int(np.argmax(np.abs(a[:, index] - b[:, index])))
        output[name] = {"row_index": row,
                        "left_value": float(a[row, index]),
                        "right_value": float(b[row, index])}
    return output


def _kernel_locations(left: Path, right: Path, names: tuple[str, ...]) -> dict:
    width = len(names)
    points = left.stat().st_size // (8 * width)
    max_abs = np.zeros(width, dtype=np.float64)
    location = np.zeros(width, dtype=np.int64)
    for offset in range(0, points, CHECKER.CHUNK_POINTS):
        count = min(CHECKER.CHUNK_POINTS, points - offset)
        values = count * width
        a = np.fromfile(left, dtype=np.float64, count=values, offset=offset * width * 8).reshape(count, width)
        b = np.fromfile(right, dtype=np.float64, count=values, offset=offset * width * 8).reshape(count, width)
        if a.size != values or b.size != values or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise QualificationError(f"truncated/non-finite kernel while locating historical error: {left}, {right}")
        errors = np.abs(a - b)
        chunk_max = errors.max(axis=0)
        for channel in range(width):
            if chunk_max[channel] > max_abs[channel]:
                max_abs[channel] = chunk_max[channel]
                location[channel] = offset + int(np.argmax(errors[:, channel]))
    return {name: {"point_index": int(location[i]), "max_abs_error": float(max_abs[i])}
            for i, name in enumerate(names)}


def _is_original_pair(label: str) -> bool:
    return any(label == f"{left} vs {right}" for left, right in ORIGINAL_PAIRS)


def _capture_historical_compare(a: dict, b: dict, exact: bool, label: str) -> dict:
    if exact or not _is_original_pair(label):
        return _ORIGINAL_COMPARE(a, b, exact, label)
    names = CHECKER.KERNEL_CHANNELS["dspec"]
    result = {"kernel": CHECKER.compare_arrays(a["kernel_path"], b["kernel_path"],
                                                math.inf, False, names)}
    result["kernel"]["max_error_locations"] = _kernel_locations(a["kernel_path"], b["kernel_path"], names)
    for field in ("trace", "spectrum"):
        if a[f"{field}_header"] != b[f"{field}_header"]:
            raise CHECKER.CheckError(f"{field} headers differ in historical comparison {label}")
        old_limit = CHECKER.PARITY_LIMIT
        try:
            CHECKER.PARITY_LIMIT = math.inf
            metrics = CHECKER.compare_numeric(a[field], b[field], False, f"{label} {field}",
                                                a[f"{field}_header"])
        finally:
            CHECKER.PARITY_LIMIT = old_limit
        metrics["max_error_locations"] = _table_locations(a[field], b[field], a[f"{field}_header"])
        result[field] = metrics
    return result


def _capture_historical_energy(left: dict, right: dict, fields: tuple[str, ...], label: str) -> dict:
    if not label.startswith(("dsv4-native-l3__dsv6-native-l3", "dsv4-profiled-l3__dsv6-profiled-l3")):
        return _ORIGINAL_ENERGY(left, right, fields, label)
    metrics = {}
    for field in fields:
        a = CHECKER.finite_number(left.get(field), f"{label}.{field}")
        b = CHECKER.finite_number(right.get(field), f"{label}.{field}")
        difference = abs(a - b)
        scale = max(abs(a), abs(b), np.finfo(float).tiny)
        tolerance = max(1.0e-14, LIMIT * scale) if "relative_error" in field else LIMIT * scale
        metrics[field] = {"left": a, "right": b, "abs_difference": difference,
                          "scaled_difference": difference / scale, "absolute_tolerance": tolerance,
                          "within_2e-13_gate": difference <= tolerance}
    return metrics


def _compare_control_outputs(control: dict, target: dict, label: str) -> dict:
    names = CHECKER.KERNEL_CHANNELS["dspec"]
    result = {"kernel": CHECKER.compare_arrays(control["kernel_path"], target["kernel_path"],
                                                math.inf, False, names)}
    for field in ("trace", "spectrum"):
        if control[f"{field}_header"] != target[f"{field}_header"]:
            raise CHECKER.CheckError(f"{field} headers differ for analytic-control comparison {label}")
        old_limit = CHECKER.PARITY_LIMIT
        try:
            CHECKER.PARITY_LIMIT = math.inf
            result[field] = CHECKER.compare_numeric(control[field], target[field], False,
                                                     f"analytic control {label} {field}",
                                                     control[f"{field}_header"])
        finally:
            CHECKER.PARITY_LIMIT = old_limit
    for field in CHECKER.DS_PHYSICAL_FIELDS:
        if field not in control["result"] or field not in target["result"]:
            raise CHECKER.CheckError(f"missing analytic-control physical option {field} ({label})")
        if control["result"][field] != target["result"][field]:
            raise CHECKER.CheckError(f"analytic-control physical option differs: {field} ({label})")
    for field in DOMAIN_FIELDS:
        if field not in control["domain"] or field not in target["domain"]:
            raise CHECKER.CheckError(f"missing analytic-control kernel domain field {field} ({label})")
        if control["domain"][field] != target["domain"][field]:
            raise CHECKER.CheckError(f"analytic-control kernel domain differs: {field} ({label})")
    result["physical_options_exact"] = True
    result["domain_options_exact"] = True
    return result


def _comparison_within(metrics: dict) -> bool:
    return all(float(item["scaled_max_abs_error"]) <= LIMIT
               for group in ("kernel", "trace", "spectrum")
               for item in (metrics[group].get("channels", {}).values()
                            if group == "kernel" else metrics[group].get("columns", {}).values()))


def _historical_status(report: dict) -> tuple[bool, dict]:
    details = {}
    parity_ok = True
    for left, right in ORIGINAL_PAIRS:
        key = f"{left}__{right}"
        # The frozen checker's key uses double underscores despite label text "vs".
        metrics = report["parity_gates"].get(key)
        if metrics is None:
            raise QualificationError(f"frozen checker omitted original comparison {key}")
        within = _comparison_within(metrics)
        details[key] = {"status": "PASS" if within else "FAIL",
                        "within_2e-13_gate": within,
                        "metrics": metrics,
                        "scope": "historical original v4-v6 result; not repaired by the analytic control"}
        parity_ok &= within
    energy_ok = True
    for pair, info in report["dspec_metadata_gates"]["v4_v6_physical_options_and_domain_consistent"].items():
        pair_energy_ok = True
        for section in ("result_energy", "domain_energy"):
            for field, metric in info[section].items():
                if not metric.get("within_2e-13_gate", True):
                    energy_ok = False
                    pair_energy_ok = False
        details.setdefault(pair, {})["energy_gate_within_2e-13"] = all(
            m.get("within_2e-13_gate", True)
            for section in ("result_energy", "domain_energy")
            for m in info[section].values())
        if not pair_energy_ok:
            details[pair]["status"] = "FAIL"
            details[pair]["within_2e-13_gate"] = False
            parity_ok = False
    parity_ok &= energy_ok
    report["original_dsv4_dsv6_parity"] = {"status": "PASS" if parity_ok else "FAIL",
                                           "within_2e-13_gate": parity_ok,
                                           "comparisons": details}
    report["gates"]["dsv4_dsv6_parity"] = parity_ok
    return parity_ok, details


def _require_receiver_predecessor(task: Path) -> dict:
    state_path = task / "state.json"
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise QualificationError(f"receiver control predecessor state unavailable or malformed: {state_path}") from exc
    jobs = state.get("jobs") if isinstance(state, dict) else None
    if not isinstance(jobs, list):
        raise QualificationError("campaign state has no job list")
    positions = {entry.get("id"): (index, entry) for index, entry in enumerate(jobs) if isinstance(entry, dict)}
    receiver = positions.get("receiver-gate")
    analytic_receiver = positions.get("analytic-control-receiver-gate")
    control = positions.get(CONTROL)
    if (receiver is None or analytic_receiver is None or control is None or
            receiver[0] >= control[0] or analytic_receiver[0] >= control[0]):
        raise QualificationError("both independent receiver gates must precede the analytic-control job")
    verified = {}
    for job_id, pair in (("receiver-gate", receiver),
                         ("analytic-control-receiver-gate", analytic_receiver)):
        entry = pair[1]
        if entry.get("status") != "completed" or entry.get("returncode") != 0 or entry.get("inputs_unchanged") is not True:
            raise QualificationError(f"{job_id} predecessor did not complete with unchanged inputs")
        log = task / job_id / "job.log"
        expected_hash = entry.get("log_sha256")
        if not isinstance(expected_hash, str) or not log.is_file() or CHECKER.sha256(log) != expected_hash:
            raise QualificationError(f"{job_id} predecessor log is missing or its state hash does not verify")
        verified[job_id] = {"completed_before_control": True,
                            "returncode": 0, "inputs_unchanged": True,
                            "log_sha256": expected_hash}
    return verified


def _old_control_diagnostic(original: dict, control: dict) -> dict:
    """Retain a three-way diagnostic without adding an old/control gate."""
    names = CHECKER.KERNEL_CHANNELS["dspec"]
    result = {"kernel": CHECKER.compare_arrays(original["kernel_path"], control["kernel_path"],
                                                math.inf, False, names)}
    result["kernel"]["max_error_locations"] = _kernel_locations(
        original["kernel_path"], control["kernel_path"], names)
    for field in ("trace", "spectrum"):
        if original[f"{field}_header"] != control[f"{field}_header"]:
            raise CHECKER.CheckError(f"{field} headers differ in old/control diagnostic")
        old_limit = CHECKER.PARITY_LIMIT
        try:
            CHECKER.PARITY_LIMIT = math.inf
            metrics = CHECKER.compare_numeric(original[field], control[field], False,
                                                f"original v4 vs analytic control {field}",
                                                original[f"{field}_header"])
        finally:
            CHECKER.PARITY_LIMIT = old_limit
        metrics["max_error_locations"] = _table_locations(
            original[field], control[field], original[f"{field}_header"])
        result[field] = metrics
    result["gate"] = "diagnostic_only"
    return result


def run(task: Path) -> dict:
    predecessor = _require_receiver_predecessor(task)
    global _ORIGINAL_COMPARE, _ORIGINAL_ENERGY
    _ORIGINAL_COMPARE = CHECKER.compare_job_outputs
    _ORIGINAL_ENERGY = CHECKER.require_close_energy
    CHECKER.compare_job_outputs = _capture_historical_compare
    CHECKER.require_close_energy = _capture_historical_energy
    try:
        report = CHECKER.run(task)
    finally:
        CHECKER.compare_job_outputs = _ORIGINAL_COMPARE
        CHECKER.require_close_energy = _ORIGINAL_ENERGY
    original_parity_ok, _ = _historical_status(report)
    report["original_parity_passed"] = original_parity_ok

    control = CHECKER.read_job(task, CONTROL)
    control_binary, points, channels, nphi, values = CHECKER.load_kernel(control, CONTROL)
    if not np.isfinite(values).all():
        raise CHECKER.CheckError("non-finite analytic-control kernel value")
    control.update({"kernel_path": control_binary, "kernel_values": values,
                    "points": points, "nphi": nphi})
    # Re-read the v6 targets using the frozen loader; the six original jobs were
    # already fully qualified by CHECKER.run above.
    targets = {}
    for name in ("dsv6-native-l3", "dsv6-profiled-l3"):
        job = CHECKER.read_job(task, name)
        binary, npoints, _, nphi, array = CHECKER.load_kernel(job, name)
        if not np.isfinite(array).all():
            raise CHECKER.CheckError(f"non-finite values in {name}")
        job.update({"kernel_path": binary, "kernel_values": array,
                    "points": npoints, "nphi": nphi})
        targets[name] = job
    control_comparisons = {}
    control_ok = True
    for name, target in targets.items():
        metrics = _compare_control_outputs(control, target, f"{CONTROL} vs {name}")
        within = _comparison_within(metrics)
        energy_metrics = {}
        for section, fields, left_key, right_key in (
                ("result", CHECKER.ENERGY_RESULT_FIELDS, "result", "result"),
                ("domain", CHECKER.ENERGY_DOMAIN_FIELDS, "domain", "domain")):
            try:
                energy_metrics[section] = CHECKER.require_close_energy(
                    control[left_key], target[right_key], fields, f"{CONTROL}__{name}.{section}")
            except CHECKER.CheckError as exc:
                within = False
                energy_metrics[section] = {"status": "FAIL", "error": str(exc)}
        metrics["energy"] = energy_metrics
        metrics["status"] = "PASS" if within else "FAIL"
        metrics["within_2e-13_gate"] = within
        control_comparisons[f"{CONTROL}__{name}"] = metrics
        control_ok &= within
    # Control results and radial geometry must still describe the same complete
    # physical sampling setup. The frozen reader validates all 179 nodes.
    control_radial_path = control["data"] / "kernel_radial_nodes.csv"
    _, control_radial = CHECKER.read_radial(control_radial_path, True)
    for name, target in targets.items():
        _, target_radial = CHECKER.read_radial(target["data"] / "kernel_radial_nodes.csv", True)
        old_limit = CHECKER.PARITY_LIMIT
        try:
            CHECKER.PARITY_LIMIT = math.inf
            radial = CHECKER.compare_numeric(control_radial, target_radial, False,
                                             f"{CONTROL} vs {name} radial nodes")
        finally:
            CHECKER.PARITY_LIMIT = old_limit
        control_comparisons[f"{CONTROL}__{name}"]["radial_nodes"] = radial
        if radial["maximum_column_scaled_error"] > LIMIT:
            control_ok = False
            control_comparisons[f"{CONTROL}__{name}"]["status"] = "FAIL"
            control_comparisons[f"{CONTROL}__{name}"]["within_2e-13_gate"] = False

    report["analytic_reference_control"] = {
        "job_id": CONTROL,
        "independent_receiver_predecessors": predecessor,
        "comparisons": control_comparisons,
        "status": "PASS" if control_ok else "FAIL",
        "within_2e-13_gate": control_ok,
        "scope": "analytic control compared separately with both original v6 native and profiled outputs"}
    original = CHECKER.read_job(task, "dsv4-native-l3")
    old_binary, old_points, _, old_nphi, old_values = CHECKER.load_kernel(original, "dsv4-native-l3")
    original.update({"kernel_path": old_binary, "kernel_values": old_values,
                     "points": old_points, "nphi": old_nphi})
    report["original_v4_vs_analytic_control_diagnostic"] = _old_control_diagnostic(original, control)
    if not control_ok:
        report["status"] = "FAIL"
        report.setdefault("errors", []).append("analytic reference control failed 2e-13 comparison against v6")
    elif not original_parity_ok:
        report["status"] = "PASS_WITH_ANALYTIC_REFERENCE_CONTROL"
    else:
        report["status"] = "PASS"
    report["schema"] = 2
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("TASK_DIR", type=Path)
    parser.add_argument("OUTPUT_JSON", type=Path)
    args = parser.parse_args()
    task = args.TASK_DIR.resolve(strict=True)
    output = args.OUTPUT_JSON.resolve(strict=False)
    if output.exists():
        parser.error(f"refusing to overwrite existing output JSON: {output}")
    report = {"schema": 2, "status": "FAIL", "task_dir": str(task), "errors": []}
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
    return 0 if report["status"] in ("PASS", "PASS_WITH_ANALYTIC_REFERENCE_CONTROL") else 1


if __name__ == "__main__":
    sys.exit(main())
