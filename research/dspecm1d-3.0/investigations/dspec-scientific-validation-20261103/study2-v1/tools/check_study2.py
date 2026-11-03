#!/usr/bin/env python3
"""Validate Study 2 one-factor controls without imposing a global pointwise limit."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
import struct
from pathlib import Path
from typing import Any


JOB_IDS = (
    "l200-dspec-sparse-bridge", "l200-dspec-h075", "l200-dspec-h05625",
    "l100-dspec", "l150-dspec", "l300-dspec",
    "l100-yspec-tight", "l150-yspec-tight", "l300-yspec-tight",
    "l200-dspec-nq8", "l200-dspec-nq10",
)
H_BASE = 0.02035382429229197
H_NOMINAL = 0.022615360324768855
MEMORY_CAP = 32 * 1024 * 1024
GEOMETRY_CAP = 16 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def finite(value: Any, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} is non-finite")
    return result


def rel_l2(reference: list[complex | float], candidate: list[complex | float]) -> dict[str, Any]:
    if len(reference) != len(candidate):
        raise ValueError("comparison vectors have different lengths")
    ref_norm = math.sqrt(sum(abs(x) ** 2 for x in reference))
    cand_norm = math.sqrt(sum(abs(x) ** 2 for x in candidate))
    error_norm = math.sqrt(sum(abs(y - x) ** 2 for x, y in zip(reference, candidate)))
    ratio = error_norm / ref_norm if ref_norm else (0.0 if error_norm == 0.0 else None)
    max_item = max(((abs(y - x), i) for i, (x, y) in enumerate(zip(reference, candidate))), default=(0.0, None))
    return {"count": len(reference), "reference_norm": ref_norm, "candidate_norm": cand_norm,
            "difference_norm": error_norm, "relative_l2_reference_denominator": ratio,
            "zero_reference_class": ref_norm == 0.0, "max_absolute_difference": max_item[0],
            "max_difference_index": max_item[1]}


def load_helpers(path: Path):
    spec = importlib.util.spec_from_file_location("study1_helpers", path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot load the sealed Study 1 checker helper module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_reference_manifest(reference_root: Path, manifest_path: Path) -> dict[str, Any]:
    manifest = read_json(manifest_path)
    if manifest.get("study1_state_sha256") != "2cb3b3161318c6c373f96a14bc94da52f5e3f0d34bae8ba028ec17f0e46e150e":
        raise ValueError("reference manifest is not bound to the finalized Study 1 producer state")
    if manifest.get("study1_raw_archive_sha256") != "f1f60203890f4e63296e1b8dc1b63cf64f6fce6d47d0067c80e39b92533ac1e5":
        raise ValueError("reference manifest is not bound to the verified Study 1 raw archive")
    count = 0
    for group in manifest.get("references", {}).values():
        for item in group:
            path = reference_root / item["path"]
            if not path.is_file() or path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
                raise ValueError(f"Study 1 reference file hash/size mismatch: {path}")
            count += 1
    for item in manifest.get("provenance_files", []):
        path = reference_root / item["path"]
        if not path.is_file() or path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Study 1 provenance file hash/size mismatch: {path}")
        count += 1
    return {"manifest_sha256": sha256(manifest_path), "verified_file_count": count,
            "study1_state_sha256": manifest["study1_state_sha256"],
            "study1_raw_archive_sha256": manifest["study1_raw_archive_sha256"]}


def verify_state(run_root: Path, plan_path: Path, state_path: Path) -> dict[str, dict[str, Any]]:
    state = read_json(state_path)
    if state.get("status") != "completed" or state.get("plan_sha256") != sha256(plan_path):
        raise ValueError("Study 2 producer state is not the completed snapshot of this exact run plan")
    by_id = {job.get("id"): job for job in state.get("jobs", [])}
    if set(by_id) != set(JOB_IDS) or len(by_id) != len(JOB_IDS):
        raise ValueError("Study 2 producer job inventory differs from the exact 11-job queue")
    for job_id in JOB_IDS:
        job = by_id[job_id]
        if (job.get("status") != "completed" or job.get("returncode") != 0 or
                job.get("failure_reasons") or job.get("inputs_unchanged") is not True or
                not job.get("cleanup", {}).get("verified_gone")):
            raise ValueError(f"Study 2 job lacks clean completion evidence: {job_id}")
        job_dir = run_root / job_id
        log = job_dir / job["log"]
        if not log.is_file() or sha256(log) != job.get("log_sha256"):
            raise ValueError(f"Study 2 producer log hash mismatch: {job_id}")
        for relative in job.get("outputs", []):
            path = job_dir / relative
            expected = job.get("outputs_sha256", {}).get(relative)
            if not path.is_file() or not expected or sha256(path) != expected:
                raise ValueError(f"Study 2 output hash mismatch: {job_id}/{relative}")
    return by_id


def data_path(run_root: Path, job_id: str) -> Path:
    return run_root / job_id / "data"


def boundary_rows(data: Path) -> list[tuple[int, float, float, int]]:
    result = []
    for row in rows(data / "mesh-elements.csv"):
        result.append((int(row["layer"]), finite(row["lower_radius_m"], "lower boundary"),
                       finite(row["upper_radius_m"], "upper boundary"), int(row["fluid"])))
    if not result:
        raise ValueError(f"{data} has no mesh element boundaries")
    return result


def trace_grid(data: Path, job_id: str) -> tuple[tuple[float, float], ...]:
    values = []
    for index, row in enumerate(rows(data / "trace.csv")):
        time = finite(row["time_s"], f"{job_id} trace time")
        window = finite(row["window"], f"{job_id} trace window")
        if window < 0.0 or window > 1.0 or (values and time <= values[-1][0]):
            raise ValueError(f"{job_id} has an invalid measurement time/window grid at row {index}")
        values.append((time, window))
    if not values:
        raise ValueError(f"{job_id} measurement time/window grid is empty")
    return tuple(values)


def validate_command(job_id: str, job: dict[str, Any]) -> None:
    command = job.get("command")
    if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
        raise ValueError(f"{job_id} has no explicit command vector")
    required = ["DSFAST_RADIAL_ROUTE=batched", "DSFAST_ANGULAR_ROUTE=real-legendre",
                "DSFAST_RADIAL_ORACLE=0", "DSFAST_ANGULAR_ORACLE=0",
                "DSFAST_NATIVE_DERIVATIVE=dref", "DSPEC_MODAL_OPERATOR_ROUTE=legacy",
                "DSPEC_MODAL_OPERATOR_ORACLE=0", "DSPEC_PHASE_C_SKIP_ZERO_KERNEL_BINS=1",
                "OMP_NUM_THREADS=1", "OMP_THREAD_LIMIT=8", "OMP_DYNAMIC=FALSE"]
    for token in required:
        if token not in command:
            raise ValueError(f"{job_id} command does not pin {token}")
    if job_id.endswith("dspec") or "-dspec-" in job_id or job_id == "l200-dspec-sparse-bridge":
        if "DSPEC_FREQUENCY_WORKERS=8" not in command:
            raise ValueError(f"{job_id} command does not pin eight DSpec frequency workers")
    else:
        for unset in ("DSPEC_FREQUENCY_WORKERS", "DSPEC_PHASE_D_STUDY2_ABSOLUTE_H"):
            if command.count("-u") < 2 or unset not in command:
                raise ValueError(f"{job_id} command does not explicitly unset {unset}")
        if any(item.startswith("DSPEC_FREQUENCY_WORKERS=") for item in command):
            raise ValueError(f"{job_id} YSpec command overrides its explicit frequency-worker unset")


def verify_inputs(args: argparse.Namespace, helper) -> tuple[list[dict[str, str]], dict[str, Any], dict[str, Any]]:
    probes = helper.rows(args.probe_csv)
    receipt = helper.read_json(args.probe_receipt)
    case, ycase = read_json(args.case), read_json(args.tight_case)
    if (sha256(args.probe_csv) != receipt.get("probe_csv_sha256") or len(probes) != 168 or
            receipt.get("probe_count") != 168):
        raise ValueError("Study 2 did not receive the exact hash-bound 168-point shared probe grid")
    for path, key in ((args.case, "case_sha256"), (args.mesh_receipt, "mesh_receipt_sha256"),
                      (args.path_reference, "path_reference_sha256"), (args.generator, "generator_sha256")):
        if sha256(path) != receipt.get(key):
            raise ValueError(f"Study 1 probe provenance hash mismatch for {key}")
    if sha256(args.model) != "4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea":
        raise ValueError("Study 2 candidate PREM deck differs from the frozen Study 1 model")
    if receipt.get("case_sha256") != sha256(args.case):
        raise ValueError("Study 2 case differs from the source case bound to the probe grid")
    expected_y = copy.deepcopy(case)
    expected_y.pop("solver", None)
    normalized_y = copy.deepcopy(ycase)
    if normalized_y.get("solver") != {"tolerance": 1e-7}:
        raise ValueError("YSpec degree cases must use exactly the approved ODE tolerance 1e-7")
    normalized_y.pop("solver")
    normalized_y["model"] = expected_y["model"]
    if normalized_y != expected_y:
        raise ValueError("YSpec tight case differs from the frozen case beyond model path/tolerance")
    if ycase.get("model") != str(args.model):
        raise ValueError("YSpec tight case does not reference the staged Study 2 PREM model")
    return probes, receipt, ycase


def check_job(data: Path, job: dict[str, Any], probes: list[dict[str, str]], receipt: dict[str, Any],
              source_xyz: tuple[float, float, float], receiver_xyz: tuple[float, float, float], helper) -> dict[str, Any]:
    meta = read_json(data / "result.json")
    backend, degree, nq = job["backend"], job["degree"], job["nq"]
    selected_h = finite(meta.get("phase_d_selected_h"), f"{job['id']}.selected_h")
    if (meta.get("backend") != backend or meta.get("degree") != degree or meta.get("nq") != nq or
            meta.get("phase_d_mesh_scale") != job["mesh_scale"] or
            meta.get("frequency_count") != 35 or meta.get("first_bin") != 7 or meta.get("last_bin") != 41 or
            meta.get("phase_d_study2_control") is not True or
            meta.get("phase_d_study2_sparse_sampling") is not True or
            meta.get("phase_d_study2_sampling_degree") != 100 or
            meta.get("phase_d_study2_absolute_h") != job["absolute_h"] or
            meta.get("phase_d_requested_relative_accuracy") != 0.001 or
            meta.get("fixed_probe_count") != 168 or
            meta.get("fixed_probe_contact_terms_included") is not False):
        raise ValueError(f"{job['id']} metadata differs from its exact Study 2 workload")
    if abs(selected_h - job["expected_h"]) > 1e-12:
        raise ValueError(f"{job['id']} selected h differs from the reviewed mesh control")
    gravity_error = finite(meta.get("gravity_max_abs_nd_error"), "gravity error")
    if gravity_error < 0.0 or gravity_error > 1e-10:
        raise ValueError(f"{job['id']} full-domain gravity control exceeds 1e-10")
    if finite(meta.get("time_shift_adjoint_identity"), "time-shift identity") <= 0.0 or abs(meta["time_shift_adjoint_identity"] - 1.0) > 1e-10:
        raise ValueError(f"{job['id']} time-shift identity is invalid")
    if finite(meta.get("measurement_rate_norm"), "measurement norm") <= 0.0:
        raise ValueError(f"{job['id']} measurement norm is not positive")
    if meta.get("along_angles") != 68 or meta.get("across_angles") != 10:
        raise ValueError(f"{job['id']} is not using the exact sampling-degree-100 section grid")
    if len(rows(data / "nodes.csv")) != meta.get("native_side_count"):
        raise ValueError(f"{job['id']} native side count differs from its node table")
    if meta.get("fixed_probe_regular_part_only") is not True:
        raise ValueError(f"{job['id']} fixed probes are not explicitly labeled regular-part samples")
    if meta.get("kernel_adjoint_weight_source") != "this mesh's own waveform-derived adjoint source":
        raise ValueError(f"{job['id']} did not use its own waveform-derived adjoint weights")
    trace = helper.read_trace(data / "trace.csv")
    if not trace:
        raise ValueError(f"{job['id']} measurement trace is empty")
    trace_grid(data, job["id"])
    weights = helper.read_weights(data)
    integrated = helper.fixed_table(data, probes, backend, source_xyz, receiver_xyz)
    if backend == "dspec":
        products = helper.validate_dspec_bin_products(data, weights, integrated)
        helper.validate_fields(data, weights, probes)
        temp_bytes = finite(meta.get("phase_d_actual_temporary_peak_bytes"), "angular temporary peak")
        cache_bytes = finite(meta.get("phase_d_actual_cache_bytes"), "angular geometry bytes")
        product_bytes = finite(meta.get("phase_c2_product_storage_bytes"), "C2 product storage")
        cache_cap = finite(meta.get("phase_d_batched_cache_max_bytes"), "angular temporary cap")
        if (meta.get("phase_c2_preweighted_complex_contraction") is not True or
                cache_cap != MEMORY_CAP or temp_bytes <= 0 or temp_bytes > MEMORY_CAP or
                cache_bytes <= 0 or cache_bytes > GEOMETRY_CAP or
                product_bytes <= 0 or product_bytes > MEMORY_CAP):
            raise ValueError(f"{job['id']} violates the existing per-worker angular/product memory contract")
        if (meta.get("native_angular_route") != "real-legendre" or
                meta.get("native_radial_route") != "batched-direct-index" or
                meta.get("native_derivative_weights_mode") != "native-dref-over-jacobian" or
                meta.get("modal_operator_route") != "legacy" or
                meta.get("phase_c_zero_bin_skips") is not True or
                meta.get("phase_c_validate_direct_adjoint") is not False or
                meta.get("phase_c_validate_reconstruction") is not False or
                meta.get("frequency_workers") != 8 or meta.get("solve_workers") != 1 or
                meta.get("workers") != 1):
            raise ValueError(f"{job['id']} actual DSpec route/workers differ from the approved workload")
    else:
        products = None
        if finite(meta.get("yspec_radial_solver_tolerance"), "YSpec ODE tolerance") != 1e-7:
            raise ValueError(f"{job['id']} YSpec ODE tolerance differs from 1e-7")
        if meta.get("workers") != 8 or meta.get("frequency_workers") != 1 or meta.get("solve_workers") != 8:
            raise ValueError(f"{job['id']} YSpec worker metadata differs from eight outer workers")
    angle_files = {name: rows(data / f"{name}-angles.csv") for name in ("along", "across")}
    for name, angle_rows in angle_files.items():
        if len(angle_rows) != (68 if name == "along" else 10):
            raise ValueError(f"{job['id']} {name} section angle inventory differs")
        raw = (data / f"{name}.bin").read_bytes()
        if len(raw) != len(rows(data / "nodes.csv")) * len(angle_rows) * 8:
            raise ValueError(f"{job['id']} {name} section binary dimensions differ")
        values = struct.unpack(f"{len(raw)//8}d", raw)
        nonfinite = sum(not math.isfinite(value) for value in values)
    return {"metadata": meta, "selected_h_parameter": selected_h, "own_P_by_bin": weights,
            "fixed_probe_integrated": integrated, "dspec_ordered_bin_products": products,
            "trace": trace, "angles": angle_files, "boundary_rows": boundary_rows(data),
            "section_diagnostics": helper.section_diagnostics(data, source_xyz, receiver_xyz),
            "section_nonfinite_counts": {name: sum(not math.isfinite(v) for v in struct.unpack(
                f"{(data / f'{name}.bin').stat().st_size//8}d", (data / f'{name}.bin').read_bytes()))
                for name in ("along", "across")}}


def compare_fixed(reference: dict[str, Any], candidate: dict[str, Any], probes: list[dict[str, str]]) -> dict[str, Any]:
    ids = sorted(set(reference) & set(candidate))
    probe_by_id = {row["probe_id"]: row for row in probes}
    groups: dict[str, list[str]] = {"all_common_numeric": [], "off_source_shell": [],
        "source_shell_inner_outer": [], "source_coincident": [], "within_25km_of_source_or_receiver": [],
        "at_least_25km_from_source_and_receiver": [], "at_least_50km_from_source_and_receiver": [],
        "at_least_100km_from_source_and_receiver": [], "far_solid_interior_100km": [],
        "solid_material_interface": []}
    for pid in ids:
        a, b = reference[pid], candidate[pid]
        if not isinstance(a["numeric_kernel"], float) or not isinstance(b["numeric_kernel"], float):
            continue
        groups["all_common_numeric"].append(pid)
        p = probe_by_id[pid]
        shell = p["layer"] == "9" and p["radial_trace"] in ("inner", "outer")
        if shell:
            groups["source_shell_inner_outer"].append(pid)
        else:
            groups["off_source_shell"].append(pid)
        if a["location_class"] == "source_coincident" or b["location_class"] == "source_coincident":
            groups["source_coincident"].append(pid)
        ds, dr = float(a["source_distance_3d_m"]), float(a["receiver_distance_3d_m"])
        if p["material_side"] == "solid" and a["location_class"] == "material_interface":
            groups["solid_material_interface"].append(pid)
        if not shell:
            if min(ds, dr) < 25_000:
                groups["within_25km_of_source_or_receiver"].append(pid)
            if ds >= 25_000 and dr >= 25_000: groups["at_least_25km_from_source_and_receiver"].append(pid)
            if ds >= 50_000 and dr >= 50_000: groups["at_least_50km_from_source_and_receiver"].append(pid)
            if ds >= 100_000 and dr >= 100_000: groups["at_least_100km_from_source_and_receiver"].append(pid)
            if (p["material_side"] == "solid" and a["location_class"] == "bulk" and
                    b["location_class"] == "bulk" and ds >= 100_000 and dr >= 100_000):
                groups["far_solid_interior_100km"].append(pid)
    output: dict[str, Any] = {}
    for label, selected in groups.items():
        left = [complex(float(reference[pid]["numeric_kernel"]), 0.0) for pid in selected]
        right = [complex(float(candidate[pid]["numeric_kernel"]), 0.0) for pid in selected]
        summary = rel_l2(left, right)
        summary["probe_ids"] = selected
        output[label] = summary
    output["unmatched_source_shell_count"] = sum(
        1 for row in probes if row["layer"] == "9" and row["radial_trace"] in ("inner", "outer") and
        (reference[row["probe_id"]]["numeric_kernel"] == "not_exposed" or
         candidate[row["probe_id"]]["numeric_kernel"] == "not_exposed"))
    return output


def read_unit_products(data: Path, weights: dict[int, complex]) -> dict[tuple[str, int], complex]:
    entries = rows(data / "fixed-probe-bin-kernel.csv")
    result: dict[tuple[str, int], complex] = {}
    if len(entries) != 168 * 35:
        raise ValueError(f"{data} does not contain the full fixed-probe/bin product grid")
    for row in entries:
        key = (row["probe_id"], int(row["bin"]))
        if key in result or key[1] not in weights:
            raise ValueError(f"{data} repeats or miskeys a fixed-probe/bin record")
        if row["unit_adjoint_computed"] == "0":
            if key[1] not in (7, 8, 41) or weights[key[1]] != 0j:
                raise ValueError(f"{data} skips a nonzero direction")
            result[key] = 0j
        elif row["unit_adjoint_computed"] == "1":
            result[key] = complex(finite(row["unit_complex_product_re"], "unit product real"),
                                  finite(row["unit_complex_product_im"], "unit product imaginary"))
        else:
            raise ValueError(f"{data} has an unknown zero-skip marker")
    if len(result) != 168 * 35:
        raise ValueError(f"{data} fixed-probe/bin inventory is incomplete")
    return result


def counterweighted_integrals(unit: dict[tuple[str, int], complex], p: dict[int, complex]) -> dict[str, float]:
    ids = sorted({key[0] for key in unit})
    return {pid: sum((p[k] * unit[(pid, k)] for k in sorted(p)), 0j).real for pid in ids}


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in ("run-root", "run-plan", "supervisor-state", "probe-csv", "probe-receipt", "case",
                 "tight-case", "model", "mesh-receipt", "path-reference", "generator", "reference-root",
                 "reference-manifest", "helpers", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    helper = load_helpers(args.helpers)
    plan = read_json(args.run_plan)
    plan_jobs = {job["id"]: job for job in plan.get("jobs", [])}
    if set(plan_jobs) != set(JOB_IDS) or len(plan_jobs) != len(JOB_IDS):
        raise ValueError("run plan does not declare exactly the reviewed 11 Study 2 jobs")
    if plan.get("thread_cap") != 20 or plan.get("study1_reference_manifest_sha256") != sha256(args.reference_manifest):
        raise ValueError("run plan thread/reference binding differs from the reviewed Study 2 plan")
    state_jobs = verify_state(args.run_root, args.run_plan, args.supervisor_state)
    reference_integrity = verify_reference_manifest(args.reference_root, args.reference_manifest)
    probe_rows, probe_receipt, tight_case = verify_inputs(args, helper)
    source = probe_receipt["source"]; receiver = probe_receipt["receiver"]
    source_xyz = helper.cart(float(source["radius_m"]), float(source["geocentric_theta_rad"]), float(source["geocentric_phi_rad"]))
    receiver_xyz = helper.cart(float(receiver["radius_m"]), float(receiver["geocentric_theta_rad"]), float(receiver["geocentric_phi_rad"]))
    jobs: dict[str, dict[str, Any]] = {}
    for job_id in JOB_IDS:
        validate_command(job_id, plan_jobs[job_id])
        config = plan_jobs[job_id]["metadata"]
        if not {"backend", "degree", "nq", "mesh_scale", "fixed_absolute_h_m", "purpose"}.issubset(config):
            raise ValueError(f"{job_id} lacks complete reviewed metadata")
        absolute = config["fixed_absolute_h_m"] is not None
        expected_h = H_BASE if absolute else H_NOMINAL * float(config["mesh_scale"])
        job = {**config, "id": job_id, "absolute_h": absolute, "expected_h": expected_h}
        result = check_job(data_path(args.run_root, job_id), job, probe_rows, probe_receipt,
                           source_xyz, receiver_xyz, helper)
        result["plan_metadata"] = job
        result["completion"] = {"returncode": state_jobs[job_id]["returncode"],
            "inputs_unchanged": state_jobs[job_id]["inputs_unchanged"],
            "cleanup_verified_gone": state_jobs[job_id]["cleanup"]["verified_gone"]}
        jobs[job_id] = result
    # All runs must use the identical spectral bins and diagnostic sparse angular grid.
    baseline_grid = (data_path(args.run_root, "l200-dspec-sparse-bridge") / "frequency-grid.csv").read_bytes()
    baseline_angles = jobs["l200-dspec-sparse-bridge"]["angles"]
    if len(baseline_angles["along"]) != 68 or len(baseline_angles["across"]) != 10:
        raise ValueError("sparse bridge grid is not the exact 68-by-10 degree-100 output grid")
    for job_id, result in jobs.items():
        data = data_path(args.run_root, job_id)
        if (data / "frequency-grid.csv").read_bytes() != baseline_grid:
            raise ValueError(f"{job_id} frequency grid differs from the baseline 35-bin grid")
        if result["angles"] != baseline_angles:
            raise ValueError(f"{job_id} sparse section coordinates differ from the common sampling grid")
    # Exact boundary/node invariance for controls that keep the physical sampling mesh fixed.
    ref_d = args.reference_root / "study1-dspec-l200"
    ref_y = args.reference_root / "study1-yspec-tight-l200"
    base_d = data_path(args.run_root, "l200-dspec-sparse-bridge")
    for job_id in JOB_IDS:
        reference_data = ref_d if job_id.endswith("dspec") or "-dspec-" in job_id or job_id == "l200-dspec-sparse-bridge" else ref_y
        reference_id = "Study1 DSpec L200" if reference_data == ref_d else "Study1 YSpec tight L200"
        if trace_grid(data_path(args.run_root, job_id), job_id) != trace_grid(reference_data, reference_id):
            raise ValueError(f"{job_id} measurement time/window grid differs from its frozen {reference_id} baseline")
    if jobs["l200-dspec-sparse-bridge"]["boundary_rows"] != boundary_rows(ref_d):
        raise ValueError("Study 2 sparse bridge changed the baseline DSpec mesh element boundaries")
    if (base_d / "nodes.csv").read_bytes() != (ref_d / "nodes.csv").read_bytes():
        raise ValueError("Study 2 sparse bridge changed baseline DSpec native nodes")
    for job_id in ("l100-dspec", "l150-dspec", "l300-dspec", "l200-dspec-nq8", "l200-dspec-nq10"):
        if jobs[job_id]["boundary_rows"] != jobs["l200-dspec-sparse-bridge"]["boundary_rows"]:
            raise ValueError(f"{job_id} fixed-h quadrature/order control changed element boundaries")
    for job_id in ("l100-dspec", "l150-dspec", "l300-dspec"):
        if (data_path(args.run_root, job_id) / "nodes.csv").read_bytes() != (base_d / "nodes.csv").read_bytes():
            raise ValueError(f"{job_id} fixed-h degree control changed NQ6 native nodes")
    if jobs["l200-dspec-nq8"]["boundary_rows"] != jobs["l200-dspec-sparse-bridge"]["boundary_rows"] or jobs["l200-dspec-nq10"]["boundary_rows"] != jobs["l200-dspec-sparse-bridge"]["boundary_rows"]:
        raise ValueError("fixed-h NQ ladder does not preserve exact element boundaries")
    # YSpec h is a sampling mesh only; verify it matches the accepted tight-L200 sampling mesh.
    for job_id in ("l100-yspec-tight", "l150-yspec-tight", "l300-yspec-tight"):
        data = data_path(args.run_root, job_id)
        if jobs[job_id]["boundary_rows"] != boundary_rows(ref_y):
            raise ValueError(f"{job_id} changed YSpec's output sampling mesh")
        if abs(jobs[job_id]["selected_h_parameter"] - H_BASE) > 1e-12:
            raise ValueError(f"{job_id} YSpec output sampling h differs from the frozen baseline")
    # Sparse-grid bridge must leave the full measured waveform, spectrum, fixed probes and mesh unchanged.
    bridge_equal_files = ("trace.csv", "spectrum.csv", "components.csv", "frequency-grid.csv", "mesh-elements.csv",
                          "nodes.csv", "fixed-probe-kernel.csv", "fixed-probe-frequency-weights.csv",
                          "fixed-probe-bin-kernel.csv", "fixed-probe-unit-fields-bin9.csv",
                          "fixed-probe-unit-fields-bin25.csv", "phase-c2-zero-bin-support.csv", "parameters.txt")
    bridge_differences = {}
    for name in bridge_equal_files:
        left, right = ref_d / name, base_d / name
        if not left.is_file() or not right.is_file():
            raise ValueError(f"sparse bridge lacks exact parity input/output {name}")
        if name == "parameters.txt":
            reference_bytes, candidate_bytes = left.read_bytes(), right.read_bytes()
            reference_lines = reference_bytes.splitlines(keepends=True)
            candidate_lines = candidate_bytes.splitlines(keepends=True)
            expected_reference_model = (
                "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/"
                "scientific-inputs/dspec-science-study1-v1/candidate/dspec/case/prem.200"
            )
            expected_candidate_model = str(args.model)
            if len(reference_lines) != len(candidate_lines) or len(reference_lines) < 2:
                raise ValueError("sparse bridge parameters.txt has a different line structure")
            if reference_lines[1].decode("utf-8").rstrip("\r\n") != f'"{expected_reference_model}"':
                raise ValueError("frozen sparse-bridge reference does not name the expected Study 1 model path")
            if candidate_lines[1].decode("utf-8").rstrip("\r\n") != f'"{expected_candidate_model}"':
                raise ValueError("sparse bridge parameters.txt does not name the exact Study 2 model input")
            model_hash = sha256(args.model)
            if model_hash != "4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea":
                raise ValueError("sparse bridge model input differs from the frozen Study 1 prem.200 deck")
            normalized_reference = reference_lines.copy()
            normalized_candidate = candidate_lines.copy()
            marker = b'"<hash-verified-prem.200>"'
            reference_ending = b"\r\n" if reference_lines[1].endswith(b"\r\n") else (b"\n" if reference_lines[1].endswith(b"\n") else b"")
            candidate_ending = b"\r\n" if candidate_lines[1].endswith(b"\r\n") else (b"\n" if candidate_lines[1].endswith(b"\n") else b"")
            normalized_reference[1] = marker + reference_ending
            normalized_candidate[1] = marker + candidate_ending
            byte_identical = b"".join(normalized_reference) == b"".join(normalized_candidate)
            bridge_differences[name] = {
                "reference_sha256": sha256(left), "candidate_sha256": sha256(right),
                "byte_identical": reference_bytes == candidate_bytes,
                "byte_identical_after_single_path_normalization": byte_identical,
                "normalization": "line 2 absolute model path only; both paths validated and prem.200 SHA-256 pinned",
                "reference_model_path": expected_reference_model,
                "candidate_model_path": expected_candidate_model,
                "model_sha256": model_hash,
            }
        else:
            bridge_differences[name] = {"reference_sha256": sha256(left), "candidate_sha256": sha256(right),
                                        "byte_identical": left.read_bytes() == right.read_bytes()}
        if not bridge_differences[name]["byte_identical"]:
            if name == "parameters.txt":
                if not bridge_differences[name]["byte_identical_after_single_path_normalization"]:
                    raise ValueError("sparse section sampling changed non-path DSpec parameter bytes")
            else:
                raise ValueError(f"sparse section sampling changed non-section DSpec output {name}")
    # Preserve own-waveform kernel and source-shell classifications; summarize rather than masking peaks.
    comparisons: dict[str, Any] = {}
    pairs = {
        "dspec_h075_vs_l200": ("l200-dspec-sparse-bridge", "l200-dspec-h075"),
        "dspec_h05625_vs_l200": ("l200-dspec-sparse-bridge", "l200-dspec-h05625"),
        "dspec_h05625_vs_h075": ("l200-dspec-h075", "l200-dspec-h05625"),
        "dspec_l100_vs_l200": ("l200-dspec-sparse-bridge", "l100-dspec"),
        "dspec_l150_vs_l200": ("l200-dspec-sparse-bridge", "l150-dspec"),
        "dspec_l150_vs_l100": ("l100-dspec", "l150-dspec"),
        "dspec_l200_vs_l150": ("l150-dspec", "l200-dspec-sparse-bridge"),
        "dspec_l300_vs_l200": ("l200-dspec-sparse-bridge", "l300-dspec"),
        "dspec_nq8_vs_nq6": ("l200-dspec-sparse-bridge", "l200-dspec-nq8"),
        "dspec_nq10_vs_nq6": ("l200-dspec-sparse-bridge", "l200-dspec-nq10"),
        "dspec_nq10_vs_nq8": ("l200-dspec-nq8", "l200-dspec-nq10"),
        "yspec_l100_vs_tight_l200": ("l200-yspec-tight", "l100-yspec-tight"),
        "yspec_l150_vs_tight_l200": ("l200-yspec-tight", "l150-yspec-tight"),
        "yspec_l150_vs_l100": ("l100-yspec-tight", "l150-yspec-tight"),
        "yspec_l300_vs_tight_l200": ("l200-yspec-tight", "l300-yspec-tight"),
    }
    for label, (reference_id, candidate_id) in pairs.items():
        ref_data = data_path(args.run_root, reference_id) if reference_id in jobs else ref_y
        cand_data = data_path(args.run_root, candidate_id)
        ref_table = jobs[reference_id]["metadata"] if reference_id in jobs else None
        if reference_id in jobs:
            ref_fixed = jobs[reference_id]["fixed_probe_integrated"]
        else:
            ref_fixed = helper.fixed_table(ref_data, probe_rows, "yspec", source_xyz, receiver_xyz)
        candidate_fixed = jobs[candidate_id]["fixed_probe_integrated"]
        comparisons[label] = {"relative_change_by_probe_subset": compare_fixed(ref_fixed, candidate_fixed, probe_rows),
                             "reference_own_P_by_bin": {str(k): [v.real, v.imag] for k, v in
                                 (jobs[reference_id]["own_P_by_bin"] if reference_id in jobs else helper.read_weights(ref_data)).items()},
                             "candidate_own_P_by_bin": {str(k): [v.real, v.imag] for k, v in jobs[candidate_id]["own_P_by_bin"].items()},
                             "waveform_and_adjoint_source_change": helper.trace_control_summary(ref_data, cand_data)}
        if label.startswith("dspec_") and not label.startswith("dspec_vs_yspec_"):
            ref_units = read_unit_products(ref_data, jobs[reference_id]["own_P_by_bin"])
            cand_units = read_unit_products(cand_data, jobs[candidate_id]["own_P_by_bin"])
            keys = sorted(ref_units)
            comparisons[label]["unweighted_unit_product_change"] = rel_l2(
                [ref_units[key] for key in keys], [cand_units[key] for key in keys])
            ref_p = jobs[reference_id]["own_P_by_bin"]
            candidate_p = jobs[candidate_id]["own_P_by_bin"]
            comparisons[label]["counterfactual_DS_candidate_unit_products_with_reference_P"] = {
                "reference_weighted_integrals": counterweighted_integrals(ref_units, ref_p),
                "candidate_unit_products_reference_weighted_integrals": counterweighted_integrals(cand_units, ref_p),
                "candidate_own_weighted_integrals": {pid: float(row["numeric_kernel"])
                    for pid, row in jobs[candidate_id]["fixed_probe_integrated"].items()
                    if isinstance(row["numeric_kernel"], float)},
                "interpretation": "DSpec unit products reweighted with the named reference P; this does not create YSpec unit fields"}
        # Native arrays are compared only when the exact physical nodes and sparse angles coincide.
        same_mesh_and_angles = (jobs[candidate_id]["boundary_rows"] ==
            (jobs[reference_id]["boundary_rows"] if reference_id in jobs else boundary_rows(ref_y)) and
            (data_path(args.run_root, candidate_id) / "nodes.csv").read_bytes() ==
            (data_path(args.run_root, reference_id) / "nodes.csv").read_bytes() if reference_id in jobs else
             (data_path(args.run_root, candidate_id) / "nodes.csv").read_bytes() == (ref_y / "nodes.csv").read_bytes())
        if reference_id in jobs and same_mesh_and_angles:
            comparisons[label]["native_section_differences"] = helper.native_section_difference(
                ref_data, cand_data, float(source["radius_m"]))
        else:
            comparisons[label]["native_section_differences"] = {
                "status": "not_paired", "reason": "native nodes or exact mesh boundaries differ; fixed physical probes are the matched comparison"}
    for degree in (100, 150, 300):
        label = f"dspec_vs_yspec_l{degree}"
        d_id, y_id = f"l{degree}-dspec", f"l{degree}-yspec-tight"
        d_data, y_data = data_path(args.run_root, d_id), data_path(args.run_root, y_id)
        comparisons[label] = {
            "relative_change_by_probe_subset": compare_fixed(
                jobs[d_id]["fixed_probe_integrated"], jobs[y_id]["fixed_probe_integrated"], probe_rows),
            "dspec_own_P_by_bin": {str(k): [v.real, v.imag] for k, v in jobs[d_id]["own_P_by_bin"].items()},
            "yspec_own_P_by_bin": {str(k): [v.real, v.imag] for k, v in jobs[y_id]["own_P_by_bin"].items()},
            "waveform_and_adjoint_source_change": helper.trace_control_summary(d_data, y_data),
            "unit_product_limitation": "YSpec unit-adjoint fields are not exposed; comparison uses each backend's own waveform-derived weights",
        }
    # Pointwise finite levels are summarized without a global convergence threshold.
    sections = {job_id: result["section_diagnostics"] for job_id, result in jobs.items()}
    h_ratios = {job_id: jobs[job_id]["selected_h_parameter"] / H_BASE
                for job_id in ("l200-dspec-sparse-bridge", "l200-dspec-h075", "l200-dspec-h05625")}
    if abs(h_ratios["l200-dspec-h075"] - .75) > 1e-12 or abs(h_ratios["l200-dspec-h05625"] - .5625) > 1e-12:
        raise ValueError("Study 2 actual selected-h ratios differ from the exact .75/.5625 ladder")
    widths = {}
    for job_id in ("l200-dspec-sparse-bridge", "l200-dspec-h075", "l200-dspec-h05625",
                   "l200-dspec-nq8", "l200-dspec-nq10"):
        boundary = jobs[job_id]["boundary_rows"]
        element_widths = [upper - lower for _, lower, upper, _ in boundary]
        widths[job_id] = {"element_count": len(element_widths), "min_width_m": min(element_widths),
                          "max_width_m": max(element_widths), "mean_width_m": sum(element_widths) / len(element_widths)}
    for fine_id, coarse_id in (("l200-dspec-h075", "l200-dspec-sparse-bridge"),
                               ("l200-dspec-h05625", "l200-dspec-h075")):
        if widths[fine_id]["max_width_m"] >= widths[coarse_id]["max_width_m"]:
            raise ValueError(f"{fine_id} did not reduce maximum physical element width relative to {coarse_id}")
    output = {
        "schema": 1,
        "study": "DSpec scientific validation Study 2 finite convergence controls",
        "science_acceptance": {"status": "PASS", "interpretation": "all reviewed one-factor controls completed and validated; no global pointwise limit or asymptotic claim is made"},
        "study1_reference_integrity": reference_integrity,
        "study2_plan_sha256": sha256(args.run_plan),
        "study2_state_sha256": sha256(args.supervisor_state),
        "study2_job_count": len(JOB_IDS),
        "study2_jobs": {job_id: {"metadata": result["metadata"],
            "selected_h_parameter": result["selected_h_parameter"], "completion": result["completion"],
            "measurement_time_window_sample_count": len(result["trace"]),
            "fixed_probe_raw_values": [{key: row.get(key) for key in (
                "probe_id", "layer", "material_side", "radial_trace", "resolved_element_trace",
                "radius_m", "theta_rad", "phi_rad", "source_distance_3d_m", "receiver_distance_3d_m",
                "numeric_kernel", "location_class", "sample_semantics")}
                for row in result["fixed_probe_integrated"].values()],
            "angle_count": {name: len(values) for name, values in result["angles"].items()},
            "section_nonfinite_counts": result["section_nonfinite_counts"],
            "section_diagnostics_unmasked": result["section_diagnostics"]} for job_id, result in jobs.items()},
        "selected_h_parameter_ratios_to_Study1_L200_NQ6": h_ratios,
        "physical_radial_element_widths": widths,
        "sparse_bridge_exact_nonsection_output_parity": bridge_differences,
        "fixed_physical_probe_comparisons": comparisons,
        "native_section_diagnostics_by_job_unmasked": sections,
        "source_receiver_semantics": {
            "source_radius_m": source["radius_m"], "receiver_radius_m": receiver["radius_m"],
            "receiver_functional": "unchanged existing solid-side measurement at the solid-fluid interface",
            "source_radius_samples": "retain inner/outer SEM traces as regular-part/contact-limited records; YSpec does not expose the same one-sided source-shell traces",
            "fixed_probe_comparisons": "common physical coordinates; source-shell and exact-coincidence categories are reported separately",
            "y_spec_mesh_parameter": "output sampling only; YSpec radial ODE tolerance is fixed at 1e-7",
            "section_exclusions": "25/50/100 km summaries are diagnostic subsets of sections, not three-dimensional norms",
        },
        "scope_limits": ["finite L100/L150/L200/L300 controls only; no L>300 or asymptotic claim",
                         "finite NQ6/NQ8/NQ10 and selected-h ratios only",
                         "pointwise source/receiver singularity or global H1/L2 regularity is not established",
                         "Study1 smooth-support low-degree 3D directional checks remain the only FD action validation"],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "out": str(args.out), "jobs": len(JOB_IDS),
                      "selected_h_ratios": h_ratios}, sort_keys=True))


if __name__ == "__main__":
    main()
