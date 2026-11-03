#!/usr/bin/env python3
"""Audit the bounded same-workload Phase E profiling queue."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import struct
from pathlib import Path


CSV_CHANNELS = {
    "trace_velocity": ("trace.csv", ("velocity_m_s",)),
    "trace_rate": ("trace.csv", ("rate_m_s2",)),
    "trace_adjoint_source": ("trace.csv", ("adjoint_source",)),
    "spectrum_complex": ("spectrum.csv", ("real_m", "imag_m")),
    "component_up": ("components.csv", ("up_re", "up_im")),
    "component_north": ("components.csv", ("north_re", "north_im")),
    "component_east": ("components.csv", ("east_re", "east_im")),
}
GEOMETRY_FILES = ("nodes.csv", "along-angles.csv", "across-angles.csv",
                  "frequency-grid.csv", "mesh-elements.csv")


def json_file(path):
    return json.loads((Path(path) / "result.json").read_text())


def csv_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_l2(a_path, b_path, name, columns):
    a, b = csv_rows(Path(a_path) / name), csv_rows(Path(b_path) / name)
    if len(a) != len(b):
        raise ValueError(f"{name}: row count differs")
    grids = {
        "trace.csv": ("time_s", "window"),
        "spectrum.csv": ("bin", "frequency_hz"),
        "components.csv": ("bin",),
    }[name]
    numerator = denominator = 0.0
    for ra, rb in zip(a, b):
        for field in grids:
            if ra[field] != rb[field]:
                raise ValueError(f"{name}: exact grid differs in {field}")
        for field in columns:
            x, y = float(ra[field]), float(rb[field])
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError(f"{name}: non-finite {field}")
            numerator += (x - y) ** 2
            denominator += y * y
    return math.sqrt(numerator / denominator) if denominator else (0.0 if numerator == 0 else math.inf)


def compare_physical(a, b):
    return {key: relative_l2(a, b, file, columns)
            for key, (file, columns) in CSV_CHANNELS.items()}


def relative_binary(a_path, b_path):
    a_bytes, b_bytes = Path(a_path).read_bytes(), Path(b_path).read_bytes()
    if len(a_bytes) != len(b_bytes) or len(a_bytes) % 8:
        raise ValueError("binary output length differs or is not float64 aligned")
    a = [value[0] for value in struct.iter_unpack("<d", a_bytes)]
    b = [value[0] for value in struct.iter_unpack("<d", b_bytes)]
    if not all(math.isfinite(value) for value in a + b):
        raise ValueError("non-finite angular output value")
    numerator = sum((x - y) ** 2 for x, y in zip(a, b))
    denominator = sum(y * y for y in b)
    return math.sqrt(numerator / denominator) if denominator else (0.0 if numerator == 0 else math.inf)


def all_output_errors(a, b):
    return {
        **compare_physical(a, b),
        "along.bin": relative_binary(Path(a) / "along.bin", Path(b) / "along.bin"),
        "across.bin": relative_binary(Path(a) / "across.bin", Path(b) / "across.bin"),
    }


def validate_result(path, backend):
    result = json_file(path)
    if result.get("backend") != backend:
        raise ValueError(f"{path}: backend mismatch")
    if (result.get("degree"), result.get("nq"), result.get("frequency_count"),
            result.get("native_side_count"), result.get("along_angles"),
            result.get("across_angles")) != (200, 6, 35, 208, 135, 19):
        raise ValueError(f"{path}: selected workload mismatch")
    if abs(float(result.get("phase_d_selected_h", -1)) - 0.02035382429229197) > 1e-14:
        raise ValueError(f"{path}: selected mesh step differs")
    if abs(float(result.get("relative_error", -1)) - 1e-5) > 1e-15:
        raise ValueError(f"{path}: parsed legacy solver relative_error differs")
    for key in ("elapsed_s", "measurement_rate_norm", "time_shift_adjoint_identity"):
        if not math.isfinite(float(result.get(key, math.nan))):
            raise ValueError(f"{path}: non-finite {key}")
    return result


def check_geometry(paths):
    hashes = {}
    for name in GEOMETRY_FILES:
        values = {sha256(Path(path) / name) for path in paths}
        if len(values) != 1:
            raise ValueError(f"requested physical geometry differs for {name}")
        hashes[name] = next(iter(values))
    return hashes


def check_yspec_profile(path):
    metadata = json.loads((Path(path) / "yspec-profile.json").read_text())
    profile = metadata["profile"]
    if profile.get("enabled") is not True:
        raise ValueError("YSpec timer-enabled library did not export an enabled recorder")
    expected_bins = list(range(7, 42))
    work = profile["work_items"]
    if work.get("forward_count") != 35 or work.get("forward_completed_bins") != expected_bins:
        raise ValueError("YSpec forward profile coverage is not all 35 bins")
    if work.get("kernel_count") != 35 or work.get("kernel_completed_bins") != expected_bins:
        raise ValueError("YSpec kernel profile coverage is not all 35 work bins")
    skipped = work.get("kernel_exact_zero_skipped_bins")
    if skipped != [7, 8, 41]:
        raise ValueError(f"YSpec exact-zero skip coverage differs: {skipped}")
    thread_rows = profile.get("thread_coverage", [])
    if not 1 <= len(thread_rows) <= 8:
        raise ValueError(f"YSpec OpenMP thread coverage count is outside [1,8]: {len(thread_rows)}")
    forward_thread_bins = []
    kernel_thread_bins = []
    kernel_thread_skips = []
    for row in thread_rows:
        if not 0 <= int(row["thread_id"]) < 8:
            raise ValueError(f"YSpec worker thread id is outside configured range: {row}")
        if int(row["forward_work_items"]) != len(row["forward_completed_bins"]):
            raise ValueError(f"YSpec per-thread forward bin/item count mismatch: {row}")
        if int(row["kernel_work_items"]) != len(row["kernel_completed_bins"]):
            raise ValueError(f"YSpec per-thread kernel bin/item count mismatch: {row}")
        forward_thread_bins.extend(row["forward_completed_bins"])
        kernel_thread_bins.extend(row["kernel_completed_bins"])
        kernel_thread_skips.extend(row["kernel_exact_zero_skipped_bins"])
    if sorted(forward_thread_bins) != expected_bins or sorted(kernel_thread_bins) != expected_bins:
        raise ValueError("YSpec per-thread bin partition does not exactly cover all 35 forward/kernel bins")
    if sorted(kernel_thread_skips) != [7, 8, 41]:
        raise ValueError("YSpec per-thread exact-zero skip bins differ from [7,8,41]")
    names = metadata.get("returned_parameter_names", [])
    if "vs" not in names or metadata.get("requested_parameter_set") != "Velocities":
        raise ValueError("YSpec actual parameter-set/output metadata is incomplete")
    worker = profile["worker_scope_sums"]
    caller = profile["caller_thread_scope_sums"]
    positive_worker = ("forward_frequency", "forward_degree_synthesis", "kernel_frequency",
                       "radial_solve_R", "radial_solve_T", "radial_solve_S",
                       "strain_table", "angular_wigner_table", "angular_table_packing",
                       "angular_evaluation", "material_contraction", "reduction")
    for name in positive_worker:
        row = worker.get(name, {})
        if int(row.get("calls", 0)) <= 0 or not float(row.get("inclusive_s", 0)) > 0:
            raise ValueError(f"YSpec expected worker scope is absent: {name}")
    positive_caller = ("outer_api", "model_and_depth_setup", "driver_setup",
                       "measurement", "time_transform", "output_io")
    for name in positive_caller:
        row = caller.get(name, {})
        if int(row.get("calls", 0)) <= 0 or not float(row.get("inclusive_s", 0)) > 0:
            raise ValueError(f"YSpec expected caller scope is absent: {name}")
    for domain in (worker, caller):
        for name, row in domain.items():
            inclusive = float(row["inclusive_s"])
            exclusive = float(row["exclusive_s"])
            if (not math.isfinite(inclusive) or not math.isfinite(exclusive)
                    or inclusive < 0 or exclusive < 0 or exclusive > inclusive + 1e-9):
                raise ValueError(f"YSpec profile has invalid {name} inclusive/exclusive times")
    bridge = profile.get("master_thread_work_item_root_scope_bridge", {})
    if not 0 <= int(bridge.get("work_item_count", -1)) <= 70:
        raise ValueError("YSpec master-thread root-scope bridge count is invalid")
    if not math.isfinite(float(bridge.get("elapsed_s", math.nan))) or float(bridge["elapsed_s"]) < 0:
        raise ValueError("YSpec master-thread root-scope bridge time is invalid")
    solver = metadata.get("radial_solver_options", {})
    if solver != {"relative_tolerance": 1e-5, "start_decay_efolds": 11.0, "step_factor": 0.25}:
        raise ValueError(f"YSpec unchanged solver defaults differ: {solver}")
    return metadata


def check_dspec_profile(path):
    profile = json.loads((Path(path) / "profile.json").read_text())
    phase = "combined_forward_unit_adjoint_reconstruction_contraction"
    rows = [row for row in profile.get("rows", []) if row.get("phase") == phase]
    by_name = {}
    for row in rows:
        name = row.get("scope_name", "").removeprefix("frequency_worker.")
        by_name[name] = row
    expected = ("combined_frequency_workspace_forward_unit_adjoint_and_reconstruction",
                "frequency_workspace_prepare_and_modal_capture",
                "phase_d_stream_projected_radial_modes_to_strain_tables",
                "phase_d_tiled_eigen_angular_reconstruction",
                "phase_e_angular_plan_and_weight_preparation",
                "phase_e_angular_mode_lookup_materialization",
                "phase_e_angular_tile_weight_packing",
                "phase_e_angular_coefficient_row_lookup_and_packing",
                "phase_e_angular_eigen_gemm",
                "phase_e_angular_tile_scatter",
                "phase_e_angular_material_evaluation_and_contraction")
    for name in expected:
        row = by_name.get(name)
        if not row or int(row.get("calls", 0)) <= 0 or not float(row.get("inclusive_s", 0)) > 0:
            raise ValueError(f"DSpec expected profile scope is absent: {name}")
        inclusive, exclusive = float(row["inclusive_s"]), float(row["exclusive_s"])
        if (not math.isfinite(inclusive) or not math.isfinite(exclusive)
                or exclusive < 0 or inclusive < exclusive):
            raise ValueError(f"DSpec scope has invalid inclusive/exclusive time: {name}")
    root = by_name["combined_frequency_workspace_forward_unit_adjoint_and_reconstruction"]
    if int(root["calls"]) != 35:
        raise ValueError(f"DSpec combined forward work must cover all 35 bins: {root['calls']}")
    workspace = by_name["frequency_workspace_prepare_and_modal_capture"]
    if int(workspace["calls"]) != 35:
        raise ValueError(f"DSpec forward workspace preparation must cover all 35 bins: {workspace['calls']}")
    counters = {(row.get("phase"), row.get("name")): int(row.get("value", 0))
                for row in profile.get("counters", [])}
    if counters.get((phase, "phase_e_angular_weight_tiles"), 0) <= 0:
        raise ValueError("DSpec angular weight-tile counter is not positive")
    if counters.get((phase, "phase_c2_exact_zero_unit_adjoint_skips"), 0) != 3:
        raise ValueError("DSpec unit-adjoint exact-zero skip count must be exactly 3")
    if counters.get((phase, "phase_c2_forward_spectrum_frequency_calls"), 0) != 35:
        raise ValueError("DSpec forward spectrum must retain all 35 frequency calls")
    if counters.get((phase, "phase_c2_unit_adjoint_frequency_calls"), 0) != 32:
        raise ValueError("DSpec active unit-adjoint work must cover the other 32 bins")
    if not profile.get("frequency_route"):
        raise ValueError("DSpec profile is not the combined frequency route")
    wall = profile.get("wall_partition", {})
    for key in ("profile_pipeline_wall_s", "master_exclusive_partition_s", "residual_s", "tolerance_s"):
        if not math.isfinite(float(wall.get(key, math.nan))):
            raise ValueError(f"DSpec profile wall partition is missing finite {key}")
    if abs(float(wall["residual_s"])) > float(wall["tolerance_s"]):
        raise ValueError("DSpec master-exclusive wall partition does not close")
    worker_exclusive = sum(float(row["exclusive_s"]) for row in rows
                           if row.get("execution_domain") == "worker_elapsed_sum")
    worker_roots = float(root["inclusive_s"])
    worker_closure = worker_exclusive - worker_roots
    if abs(worker_closure) > max(1e-6, 1e-3 * max(worker_roots, 1e-9)):
        raise ValueError(f"DSpec worker scope exclusive-time sum does not close: {worker_closure}")
    return profile, {name: by_name[name] for name in expected}


def main():
    parser = argparse.ArgumentParser()
    for option in ("d-native-1", "d-native-2", "d-native-3", "d-profile",
                   "y-native-1", "y-native-2", "y-native-3", "y-profile",
                   "d-reference", "y-reference"):
        parser.add_argument("--" + option, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    d_native = [getattr(args, f"d_native_{i}") for i in (1, 2, 3)]
    y_native = [getattr(args, f"y_native_{i}") for i in (1, 2, 3)]
    d_reference, y_reference = args.d_reference, args.y_reference
    d_paths = d_native + [args.d_profile, d_reference]
    y_paths = y_native + [args.y_profile, y_reference]
    d_results = [validate_result(path, "dspec") for path in d_paths]
    y_results = [validate_result(path, "yspec") for path in y_paths]
    geometry = check_geometry(d_paths + y_paths)

    parity = {}
    for label, paths, ref in (("dspec", d_paths[:-1], d_reference),
                              ("yspec", y_paths[:-1], y_reference)):
        all_errors = {}
        for index, path in enumerate(paths):
            errors = all_output_errors(path, ref)
            all_errors[f"run_{index + 1}"] = errors
            if max(errors.values()) > 1e-11:
                raise ValueError(f"{label} output parity to immutable Phase D reference exceeds 1e-11: {errors}")
        parity[label] = all_errors

    d_native_errors = [all_output_errors(path, d_native[0]) for path in d_native[1:]]
    y_native_errors = [all_output_errors(path, y_native[0]) for path in y_native[1:]]
    if max(v for result in d_native_errors + y_native_errors for v in result.values()) > 1e-11:
        raise ValueError("same-backend native repeat outputs exceed the 1e-11 parity gate")

    d_profile, d_scopes = check_dspec_profile(args.d_profile)
    y_profile = check_yspec_profile(args.y_profile)
    d_times = [float(validate_result(path, "dspec")["elapsed_s"]) for path in d_native]
    y_times = [float(validate_result(path, "yspec")["elapsed_s"]) for path in y_native]
    def distribution(times):
        median = statistics.median(times)
        return {"seconds": times, "median_s": median, "min_s": min(times), "max_s": max(times),
                "range_s": max(times) - min(times),
                "relative_range_to_median": (max(times) - min(times)) / median}

    y_record = y_profile["profile"]
    y_worker = y_record["worker_scope_sums"]
    y_caller = y_record["caller_thread_scope_sums"]
    y_worker_exclusive = sum(float(row["exclusive_s"]) for row in y_worker.values())
    y_worker_root_inclusive = sum(float(y_worker[name]["inclusive_s"])
                                  for name in ("forward_frequency", "kernel_frequency"))
    y_worker_closure_residual = y_worker_exclusive - y_worker_root_inclusive
    if abs(y_worker_closure_residual) > max(1e-6, 1e-3 * max(y_worker_root_inclusive, 1e-9)):
        raise ValueError(f"YSpec worker scope exclusive-time sum does not close: {y_worker_closure_residual}")
    y_caller_exclusive = sum(float(row["exclusive_s"]) for row in y_caller.values())
    y_bridge = float(y_record["master_thread_work_item_root_scope_bridge"]["elapsed_s"])
    y_driver_elapsed = float(y_results[3]["elapsed_s"])
    y_api_elapsed = float(y_caller.get("outer_api", {}).get("inclusive_s", 0.0))
    y_unscoped_driver_residual = y_driver_elapsed - y_caller_exclusive - y_bridge
    d_native_median = statistics.median(d_times)
    y_native_median = statistics.median(y_times)

    report = {
        "status": "pass",
        "scope": "Stage 1 matched instrumentation-only profiling; no optimization or Stage 2",
        "selected_case": {"degree": 200, "nq": 6, "mesh_step_h": 0.02035382429229197,
                          "frequency_bins": list(range(7, 42)), "physical_sides": 208,
                          "angles": {"along": 135, "across": 19}, "max_output_relative_l2": 1e-11},
        "requested_same_case_is_not_identical_backend_work": {
            "dspec_parameter": "relative lnVs contracted directly",
            "yspec_requested_set": y_profile["requested_parameter_set"],
            "yspec_actual_returned_parameter_names": y_profile["returned_parameter_names"],
            "yspec_selected_parameter": y_profile["selected_output_parameter"],
            "yspec_forward_channels": y_profile["forward_output_channels"],
            "yspec_radial_solver_settings": y_profile["radial_solver_options"],
            "yspec_internal_radial_mesh_note": y_profile["yspec_internal_radial_mesh"],
        },
        "native_timing_repeat_distribution_driver_elapsed_s": {
            "dspec": distribution(d_times), "yspec": distribution(y_times)},
        "profile_run_driver_elapsed_s_single_observation": {
            "dspec": float(d_results[3]["elapsed_s"]), "yspec": float(y_results[3]["elapsed_s"])},
        "instrumentation_overhead_vs_fresh_native_median": {
            "dspec_driver_elapsed_s": float(d_results[3]["elapsed_s"]) - d_native_median,
            "dspec_relative": float(d_results[3]["elapsed_s"]) / d_native_median - 1.0,
            "yspec_driver_elapsed_s": y_driver_elapsed - y_native_median,
            "yspec_relative": y_driver_elapsed / y_native_median - 1.0,
            "context": "single enabled observation versus three timer-disabled same-source/compiler controls; native min/max/range are reported separately"},
        "scope_time_accounting": {
            "yspec_worker_scope_exclusive_sum_s": y_worker_exclusive,
            "yspec_forward_and_kernel_root_inclusive_sum_s": y_worker_root_inclusive,
            "yspec_worker_closure_residual_s": y_worker_closure_residual,
            "yspec_caller_scope_exclusive_sum_s": y_caller_exclusive,
            "yspec_master_frequency_root_bridge_s": y_bridge,
            "yspec_outer_api_inclusive_s": y_api_elapsed,
            "yspec_profile_driver_elapsed_s": y_driver_elapsed,
            "yspec_unscoped_driver_residual_after_caller_and_bridge_s": y_unscoped_driver_residual,
            "dspec_profile_wall_partition": d_profile.get("wall_partition"),
            "dspec_combined_worker_scope_exclusive_sum_s": sum(
                float(row["exclusive_s"]) for row in d_profile.get("rows", [])
                if row.get("phase") == "combined_forward_unit_adjoint_reconstruction_contraction"
                and row.get("execution_domain") == "worker_elapsed_sum"),
            "dspec_combined_frequency_root_inclusive_s": float(d_scopes[
                "combined_frequency_workspace_forward_unit_adjoint_and_reconstruction"]["inclusive_s"]),
            "interpretation": "Worker sums are elapsed across frequency workers and overlap the caller frequency-pass interval; caller API and driver residuals are separate domains, not serial critical-path additions."},
        "same_backend_output_parity_to_immutable_phase_d_reference_by_run": parity,
        "same_backend_native_repeat_relative_l2": {"dspec": d_native_errors, "yspec": y_native_errors},
        "requested_grid_sha256": geometry,
        "cross_backend_output_relative_l2_diagnostic_only": all_output_errors(d_native[0], y_native[0]),
        "dspec_profile": {"scope_rows": d_profile.get("rows", []),
                          "counters": d_profile.get("counters", []),
                          "required_angular_scopes": d_scopes},
        "yspec_profile": y_profile["profile"],
        "interpretation_limits": [
            "Supervisor wall, driver elapsed, caller-thread scope totals, and summed worker scope times are distinct measurements.",
            "Worker elapsed totals sum across frequency tasks and must not be subtracted from master or supervisor wall.",
            "YSpec computes ParameterSet::Velocities and the adapter selects vs; DSpec contracts lnVs directly.",
            "YSpec radial integration uses its own DepthSolves controls; the DSpec SEM mesh is only the shared sampled-output grid.",
            "Cross-backend raw angular arrays are diagnostic only; parity is gated within each backend against its own Phase D reference.",
        ],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
