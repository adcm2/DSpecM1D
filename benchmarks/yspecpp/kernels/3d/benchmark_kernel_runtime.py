#!/usr/bin/env python3
"""Run a matched DSpecM1D/YSpecpp 3D kernel timing or replay saved measurements."""
from __future__ import annotations

import argparse
import csv
import hashlib
from html import escape
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MODULI = "ACFLN"
OBJECTIVES = ("J_Z", "J_N")
TIMEOUT_PER_CALL_S = 300
HARD_RUNTIME_S = 900


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def strict_json_loads(text: str):
    def reject_constant(value):
        raise ValueError(f"nonstandard JSON number {value}")
    return json.loads(text, parse_constant=reject_constant)


def ensure_finite_tree(value, label: str) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            ensure_finite_tree(item, f"{label}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            ensure_finite_tree(item, f"{label}[{index}]")
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        if not np.isfinite(value):
            raise ValueError(f"nonfinite numeric value in {label}: {value!r}")


def validate_samples_csv(path: Path, expected_count: int) -> dict[str, int]:
    required = ("objective", "modulus", "radius_m", "theta_rad", "phi_rad",
                "element", "fluid", "kernel_objective_per_Pa_m3", "deltaX_Pa", "dr_m", "dOmega")
    counts = {modulus: 0 for modulus in MODULI}
    for row_number, row in enumerate(read_csv(path), start=2):
        if row.get("objective") != "J_Z" or row.get("modulus") not in counts:
            raise ValueError(f"unexpected objective/modulus in shared samples row {row_number}")
        for key in required[2:]:
            finite_float(row, key)
        radius = float(row["radius_m"])
        theta = float(row["theta_rad"])
        phi = float(row["phi_rad"])
        dr, domega = float(row["dr_m"]), float(row["dOmega"])
        if int(row["element"]) < 0 or int(row["fluid"]) not in (0, 1):
            raise ValueError(f"invalid element/fluid metadata in shared samples row {row_number}")
        dvolume = radius * radius * dr * domega
        if (radius <= 0 or not 0.0 <= theta <= np.pi or not -2 * np.pi <= phi <= 2 * np.pi or
                dr <= 0 or domega <= 0 or not np.isfinite(dvolume) or dvolume <= 0):
            raise ValueError(f"invalid radius, angle, or volume weight in shared samples row {row_number}")
        counts[row["modulus"]] += 1
    if not counts or any(count != expected_count for count in counts.values()):
        raise ValueError(f"shared sample row counts differ from tutorial metadata: {counts}")
    return counts


def finite_float(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not np.isfinite(value):
        raise ValueError(f"nonfinite {key} in {row}")
    return value


def kernel_vectors(path: Path) -> dict[tuple[str, str], tuple[list[tuple[float, float, float]], np.ndarray]]:
    grouped: dict[tuple[str, str], tuple[list[tuple[float, float, float]], list[float]]] = {}
    for row in read_csv(path):
        objective, modulus = row["objective"], row["modulus"]
        if objective not in OBJECTIVES or modulus not in MODULI:
            raise ValueError(f"unexpected objective or modulus in {path}: {row}")
        key = (objective, modulus)
        coords, values = grouped.setdefault(key, ([], []))
        coords.append(tuple(finite_float(row, name) for name in
                            ("radius_m", "theta_rad", "phi_rad")))
        values.append(finite_float(row, "kernel_per_Pa_m3"))
    if not grouped:
        raise ValueError(f"empty kernel output: {path}")
    return {key: (coords, np.asarray(values)) for key, (coords, values) in grouped.items()}


def accuracy_gate(dspec_csv: Path, yspec_csv: Path, samples_csv: Path,
                  tolerance: float, measurement_count: int,
                  dspec_result: dict, tutorial_metadata: dict | None = None) -> dict:
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("accuracy tolerance must be finite and positive")
    dspec = kernel_vectors(dspec_csv)
    yspec = kernel_vectors(yspec_csv)
    objectives = ("J_Z", "J_N")[:measurement_count]
    expected = {(objective, modulus) for objective in objectives for modulus in MODULI}
    if set(dspec) != expected or set(yspec) != expected:
        raise ValueError("timed solvers did not return every requested objective and modulus")
    reference: dict[str, list[dict]] = {modulus: [] for modulus in MODULI}
    for row in read_csv(samples_csv):
        if row["objective"] != "J_Z":
            continue
        modulus = row["modulus"]
        if modulus not in MODULI:
            raise ValueError(f"unexpected modulus in tutorial sample grid: {modulus}")
        coord = tuple(finite_float(row, name) for name in
                      ("radius_m", "theta_rad", "phi_rad"))
        radius = coord[0]
        dvolume = radius * radius * finite_float(row, "dr_m") * finite_float(row, "dOmega")
        reference[modulus].append({
            "coord": coord,
            "dvolume": dvolume,
            "deltaX_Pa": finite_float(row, "deltaX_Pa"),
            "tutorial_kernel": finite_float(row, "kernel_objective_per_Pa_m3"),
        })
    if any(not rows for rows in reference.values()):
        raise ValueError("tutorial sample grid is missing one or more moduli")
    scales = dspec_result.get("objective_scale_S", {})
    if any(not np.isfinite(float(scales.get(obj, 0.0))) or float(scales.get(obj, 0.0)) <= 0
           for obj in objectives):
        raise ValueError("DSpec timed output lacks a finite positive objective scale S")

    baseline_moduli = {"A": 352.96e9, "C": 352.96e9, "F": 129.6025e9,
                       "L": 111.67875e9, "N": 111.67875e9}
    profile_rows = []
    integral_rows = []
    tutorial_rows = []
    tutorial_integrals = []
    profile_pass = True
    integral_pass = True
    tutorial_pass = True

    def compare_profiles(first: np.ndarray, second: np.ndarray, weights: np.ndarray,
                         scale: float) -> dict:
        total_weight = float(np.sum(weights))
        if not total_weight > 0:
            raise ValueError("profile quadrature has zero physical weight")
        first_rms = float(np.sqrt(np.dot(weights, first * first) / total_weight))
        second_rms = float(np.sqrt(np.dot(weights, second * second) / total_weight))
        difference_rms = float(np.sqrt(np.dot(weights, (first - second) ** 2) / total_weight))
        informative = min(first_rms, second_rms) > 1e-6 * scale
        return {"weighted_rms_first": first_rms, "weighted_rms_second": second_rms,
                "weighted_rms_difference": difference_rms, "informative": informative,
                "relative_difference": difference_rms / max(first_rms, second_rms)
                    if informative else None,
                "scaled_absolute_difference": difference_rms / scale
                    if not informative else None,
                "gate_pass": bool(informative and
                    difference_rms / max(first_rms, second_rms) <= tolerance)}

    for objective, modulus in sorted(expected):
        dcoords, dvalues = dspec[(objective, modulus)]
        ycoords, yvalues = yspec[(objective, modulus)]
        expected_coords = [row["coord"] for row in reference[modulus]]
        if dcoords != expected_coords or ycoords != expected_coords:
            raise ValueError(f"sample count, ordering, or coordinates differ for {objective}/{modulus}")
        dvolumes = np.asarray([row["dvolume"] for row in reference[modulus]])
        deltas = np.asarray([row["deltaX_Pa"] for row in reference[modulus]])
        weights = dvolumes * deltas * deltas
        support_volume = float(np.sum(dvolumes))
        signal_scale = float(scales[objective])
        kernel_scale = signal_scale / (baseline_moduli[modulus] * support_volume)
        profile = compare_profiles(dvalues, yvalues, weights, kernel_scale)
        profile_rows.append({"objective": objective, "modulus": modulus,
                             "sample_count": len(dvalues), "signal_scale_S": signal_scale,
                             "Kscale_per_Pa_m3": kernel_scale, **profile})
        profile_pass &= profile["gate_pass"]

        integral_d = float(np.dot(dvolumes * deltas, dvalues))
        integral_y = float(np.dot(dvolumes * deltas, yvalues))
        informative_integral = min(abs(integral_d), abs(integral_y)) > 1e-6 * signal_scale
        relative_integral = abs(integral_d - integral_y) / max(
            abs(integral_d), abs(integral_y), 1e-14 * signal_scale)
        scaled_integral = abs(integral_d - integral_y) / signal_scale
        passed_integral = bool(informative_integral and relative_integral <= tolerance)
        integral_pass &= passed_integral
        integral_rows.append({"objective": objective, "modulus": modulus,
                              "I_D": integral_d, "I_Y": integral_y,
                              "informative": informative_integral,
                              "relative_difference": relative_integral if informative_integral else None,
                              "scaled_absolute_difference": scaled_integral
                                  if not informative_integral else None,
                              "gate_pass": passed_integral})

        if objective == "J_Z":
            tutorial = np.asarray([row["tutorial_kernel"] for row in reference[modulus]])
            tutorial_profile = compare_profiles(dvalues, tutorial, weights, kernel_scale)
            tutorial_rows.append({"modulus": modulus, **tutorial_profile})
            tutorial_pass &= tutorial_profile["gate_pass"]
            tutorial_integral = float(np.dot(dvolumes * deltas, tutorial))
            tutorial_integral_informative = min(abs(integral_d), abs(tutorial_integral)) > 1e-6 * signal_scale
            tutorial_integral_relative = abs(integral_d - tutorial_integral) / max(
                abs(integral_d), abs(tutorial_integral), 1e-14 * signal_scale)
            tutorial_integral_pass = bool(tutorial_integral_informative and
                                          tutorial_integral_relative <= tolerance)
            tutorial_pass &= tutorial_integral_pass
            tutorial_integrals.append({"modulus": modulus, "I_timed_D": integral_d,
                                       "I_tutorial": tutorial_integral,
                                       "informative": tutorial_integral_informative,
                                       "relative_difference": tutorial_integral_relative
                                           if tutorial_integral_informative else None,
                                       "scaled_absolute_difference": abs(integral_d - tutorial_integral) /
                                           signal_scale if not tutorial_integral_informative else None,
                                       "gate_pass": tutorial_integral_pass})

    if not profile_pass or not integral_pass or not tutorial_pass:
        raise ValueError("actual timed-output accuracy gate failed or was uninformative "
                         f"(profiles={profile_pass}, integrals={integral_pass}, tutorial={tutorial_pass})")
    return {"passed": True, "profile_tolerance": tolerance,
            "integral_tolerance": tolerance, "profile_metric": "dV*deltaX^2 weighted RMS",
            "weak_signal_floor": "1e-6 of Kscale or S; weak cases have no relative PASS",
            "profiles": profile_rows, "integrals": integral_rows,
            "dspec_tutorial_profile_comparison": tutorial_rows,
            "dspec_tutorial_integral_comparison": tutorial_integrals}



def git_info(source: Path) -> dict:
    def run(args: list[str]) -> str:
        return subprocess.run(args, cwd=source, check=True, text=True,
                              capture_output=True).stdout.strip()
    head = run(["git", "rev-parse", "HEAD"])
    diff = subprocess.run(["git", "diff", "--binary", "HEAD"], cwd=source,
                          check=True, capture_output=True).stdout
    untracked = run(["git", "ls-files", "--others", "--exclude-standard"]).splitlines()
    untracked_hashes = {name: sha256(source / name) for name in untracked
                        if (source / name).is_file()}
    return {"path": str(source.resolve()), "head": head,
            "dirty_paths": run(["git", "status", "--short"]),
            "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
            "untracked_file_sha256": untracked_hashes}


def build_info(source: Path, build: Path, executable: Path, target: str,
               wrapper_source_root: Path | None = None) -> dict:
    cache = build / "CMakeCache.txt"
    if not cache.is_file():
        raise ValueError(f"missing CMakeCache.txt: {cache}")
    cache_text = cache.read_text(encoding="utf-8", errors="replace")
    values = {}
    for line in cache_text.splitlines():
        match = re.match(r"([^:#]+):[^=]+=(.*)$", line)
        if match:
            values[match.group(1)] = match.group(2)
    build_type = values.get("CMAKE_BUILD_TYPE", "")
    if build_type.lower() != "release":
        raise ValueError(f"{build} is not configured as Release (got {build_type!r})")
    compiler = values.get("CMAKE_CXX_COMPILER", "")
    cache_source = Path(values.get("CMAKE_HOME_DIRECTORY", "")).resolve()
    flags = values.get("CMAKE_CXX_FLAGS_RELEASE", "")
    flag_files = list(build.rglob(f"CMakeFiles/{target}.dir/flags.make"))
    if len(flag_files) != 1:
        raise ValueError(f"expected one target-specific flags.make for {target}")
    target_flags = flag_files[0].read_text(errors="replace")
    target_compile_flags = next((line.split("=", 1)[1].strip()
                                 for line in target_flags.splitlines()
                                 if line.startswith("CXX_FLAGS =")), "")
    if "-O3" not in target_compile_flags:
        raise ValueError(f"timing target {target} does not show explicit -O3 flags")
    build_make_files = list(build.rglob(f"CMakeFiles/{target}.dir/build.make"))
    if len(build_make_files) != 1:
        raise ValueError(f"expected one target-specific build.make for {target}")
    target_build = build_make_files[0].read_text(errors="replace")
    expected_source = source.resolve()
    wrapper_root = (wrapper_source_root or source).resolve()
    target_source = (wrapper_root / "tutorials/sensitivity_kernels_3d/kernel_timing.cpp"
                     if target == "dspec_3d_kernel_timing" else
                     wrapper_root / "benchmarks/yspecpp/kernels/3d/yspec_kernel_timing.cpp")
    expected_cmake_source = target_source.parent
    if cache_source != expected_cmake_source:
        raise ValueError(f"build cache source differs from timing target source: {cache_source}")
    if str(target_source) not in target_build:
        raise ValueError(f"target build metadata does not reference the supplied source checkout: {target}")
    if target == "yspec_kernel_timing":
        yspec_source = Path(values.get("YSPEC_SOURCE_DIR", "")).resolve()
        if yspec_source != expected_source:
            raise ValueError(f"YSPEC_SOURCE_DIR differs from supplied YSpec checkout: {yspec_source}")
        library_builds = list(build.rglob("CMakeFiles/YSpec.dir/build.make"))
        if len(library_builds) != 1 or str(expected_source / "src/") not in library_builds[0].read_text(
                encoding="utf-8", errors="replace"):
            raise ValueError("YSpec library build metadata does not reference the supplied checkout")
    package_config = None
    if target == "dspec_3d_kernel_timing":
        package_config = Path(values.get("DSpecM1D_DIR", "")).resolve()
        package_file = package_config / "DSpecM1DConfig.cmake"
        if not package_file.is_file():
            raise ValueError(f"installed DSpecM1D package config is missing: {package_file}")
        if str(package_config.parent.parent.parent / "include") not in target_flags:
            raise ValueError("DSpec timing target flags do not include the installed DSpecM1D headers")
    if not executable.is_file():
        raise ValueError(f"timing executable does not exist: {executable}")
    compiler_version = subprocess.run([compiler, "--version"], check=True, text=True,
                                       capture_output=True).stdout.splitlines()[0]
    dependencies = {}
    for name, value in values.items():
        if name.startswith("FETCHCONTENT_SOURCE_DIR_") and value:
            dep = Path(value)
            if dep.is_dir():
                try:
                    dependencies[name] = git_info(dep)
                except subprocess.CalledProcessError:
                    dependencies[name] = {"path": str(dep.resolve()), "git": "unavailable"}
    for dep in build.rglob("*-src"):
        if dep.is_dir() and dep.name not in dependencies:
            try:
                dependencies[dep.name] = git_info(dep)
            except subprocess.CalledProcessError:
                dependencies[dep.name] = {"path": str(dep.resolve()), "git": "unavailable"}
    pin_prefixes = ("DSPECM1D_", "GSHTRANS_", "PLANETMODEL_") if target == "dspec_3d_kernel_timing" else (
        "YSPEC_", "PLANETMODEL_", "GSHTRANS_")
    selected_dependency_pins = {
        key: value for key, value in values.items()
        if key.endswith("_TAG") and key.startswith(pin_prefixes)
    }
    return {"source": git_info(source), "wrapper_source": git_info(wrapper_root),
            "build_dir": str(build.resolve()),
            "cmake_build_type": build_type, "compiler": compiler,
            "compiler_version": compiler_version, "release_flags": flags,
            "target_compile_flags": target_compile_flags,
            "cache_source_directory": str(cache_source),
            "installed_package_config": str(package_config) if package_config else None,
            "installed_package_config_sha256": (sha256(package_config / "DSpecM1DConfig.cmake")
                                                 if package_config else None),
            "target_flags_sha256": hashlib.sha256(target_flags.encode()).hexdigest(),
            "selected_dependency_pins": selected_dependency_pins,
            "dependencies": dependencies, "executable": str(executable.resolve()),
            "executable_sha256": sha256(executable),
            "cmake_cache_sha256": sha256(cache)}


def thread_environment() -> dict[str, str]:
    env = os.environ.copy()
    limits = {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
              "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
              "VECLIB_MAXIMUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
    env.update(limits)
    return {"env": env, "limits": limits}


def validate_tutorial_metadata(path: Path, lmax: int, model_hash: str,
                               samples_hash: str) -> dict:
    if not path.is_file():
        raise ValueError(f"missing tutorial run metadata: {path}")
    metadata = strict_json_loads(path.read_text(encoding="utf-8"))
    ensure_finite_tree(metadata, "tutorial metadata")
    expected = {"solver": "DSpecM1D", "objective": "J_Z=Re(sum(conj(q)*P u))",
                "lmax": lmax, "frequency_bins_zero_based": [1, 2, 3],
                "nt": 128, "gravity": "full", "attenuation": False}
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"tutorial metadata mismatch for {key}: {metadata.get(key)!r}")
    expected_angles = max(8, 2 * lmax + 2) * max(24, 4 * lmax + 4)
    if metadata.get("angular_nodes") != expected_angles:
        raise ValueError("tutorial angular node count differs from the requested lmax")
    if not isinstance(metadata.get("sample_count"), int) or metadata["sample_count"] <= 0:
        raise ValueError("tutorial metadata has invalid sample_count")

    def close(path_name: str, actual, expected_value: float) -> None:
        try:
            value = float(actual)
        except (TypeError, ValueError) as error:
            raise ValueError(f"tutorial metadata is missing numeric {path_name}") from error
        if not np.isfinite(value) or not np.isclose(value, expected_value, rtol=1e-12, atol=0.0):
            raise ValueError(f"tutorial metadata mismatch for {path_name}: {actual!r}")

    for key, value in {"dt_s": 20.0, "df_hz": 1.0 / 2560.0,
                       "epsilon_s_inv": 1.0 / 360.0,
                       "physical_G_SI": 6.6723e-11}.items():
        close(key, metadata.get(key), value)
    for name, location_expected in {
        "source": {"latitude_deg": 10.0, "longitude_deg": 20.0, "depth_km": 50.0},
        "receiver": {"latitude_deg": 45.0, "longitude_deg": 90.0, "depth_km": 3.0},
    }.items():
        location = metadata.get(name, {})
        for key, value in location_expected.items():
            close(f"{name}.{key}", location.get(key), value)
    moment = np.asarray(metadata.get("moment_tensor_Nm", []), dtype=float)
    if moment.shape != (6,) or not np.isfinite(moment).all() or not np.allclose(
            moment, [1e19, 2e18, -3e18, -5e18, 4e18, -5e18], rtol=1e-12, atol=0):
        raise ValueError("tutorial metadata mismatch for moment_tensor_Nm")
    if metadata.get("q_nonzero") != {"Z,1": [0.7, 0.2], "Z,2": [-0.25, 0.45]}:
        raise ValueError(f"tutorial metadata mismatch for q_nonzero: {metadata.get('q_nonzero')!r}")
    window = metadata.get("time_window", {})
    expected_window = {"type": "Hann", "sample_start_inclusive": 0,
                       "sample_end_inclusive": 90, "sample_count": 128}
    if any(window.get(key) != value for key, value in expected_window.items()):
        raise ValueError(f"tutorial metadata mismatch for time_window: {window!r}")
    provenance = metadata.get("provenance", {})
    if provenance.get("model_sha256") != model_hash:
        raise ValueError("tutorial metadata mismatch for provenance.model_sha256")
    sample_hash = provenance.get("samples_sha256")
    if sample_hash is None:
        sample_hash = metadata.get("render_provenance", {}).get(
            "numeric_inputs_sha256", {}).get("samples.csv")
    if sample_hash != samples_hash:
        raise ValueError("tutorial metadata mismatch for samples.csv SHA-256")
    return metadata


def assert_configs_match(first: dict, second: dict) -> None:
    ensure_finite_tree(first, "DSpec resolved config")
    ensure_finite_tree(second, "YSpec resolved config")
    def compare(left, right, path="config"):
        if isinstance(left, dict) and isinstance(right, dict):
            if set(left) != set(right):
                raise ValueError(f"resolved configuration keys differ at {path}")
            for key in left:
                compare(left[key], right[key], f"{path}.{key}")
        elif isinstance(left, list) and isinstance(right, list):
            if len(left) != len(right):
                raise ValueError(f"resolved configuration lengths differ at {path}")
            for index, (a, b) in enumerate(zip(left, right)):
                compare(a, b, f"{path}[{index}]")
        elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
            if not np.isclose(left, right, rtol=1e-12, atol=0.0):
                raise ValueError(f"resolved configuration mismatch at {path}: {left} != {right}")
        elif left != right:
            raise ValueError(f"resolved configuration mismatch at {path}: {left!r} != {right!r}")
    compare(first, second)


def assert_config_matches_tutorial(config: dict, metadata: dict, lmax: int,
                                   measurement_count: int, sample_count: int) -> None:
    ensure_finite_tree(config, "resolved timing config")
    grid = config["frequency_grid"]
    expected_grid = {"nt": metadata["nt"], "dt_s": metadata["dt_s"],
                     "df_hz": metadata["df_hz"],
                     "epsilon_s_inv": metadata["epsilon_s_inv"],
                     "solver_bins_zero_based": [1, 2, 3],
                     "degrees_inclusive": [0, lmax]}
    for key, expected in expected_grid.items():
        actual = grid.get(key)
        if isinstance(expected, list):
            if actual != expected:
                raise ValueError(f"timing configuration mismatch with tutorial at frequency_grid.{key}")
        elif not np.isfinite(float(actual)) or not np.isclose(
                float(actual), float(expected), rtol=1e-12, atol=0.0):
            raise ValueError(f"timing configuration mismatch with tutorial at frequency_grid.{key}")
    physics = config["physics"]
    if physics.get("gravity") != metadata["gravity"] or physics.get("attenuation") != metadata["attenuation"]:
        raise ValueError("timing physics configuration differs from frozen tutorial metadata")
    if not np.isclose(float(physics.get("physical_G")), float(metadata["physical_G_SI"]),
                      rtol=1e-12, atol=0.0):
        raise ValueError("timing physical G differs from frozen tutorial metadata")
    for config_name, metadata_name in (("source", "source"), ("receiver", "receiver")):
        actual = config[config_name]
        expected = metadata[metadata_name]
        values = {"geocentric_latitude_deg": expected["latitude_deg"],
                  "longitude_deg": expected["longitude_deg"],
                  "depth_m": float(expected["depth_km"]) * 1000.0}
        for key, value in values.items():
            if not np.isclose(float(actual.get(key)), float(value), rtol=1e-12, atol=0.0):
                raise ValueError(f"timing {config_name} differs from tutorial at {key}")
    if not np.allclose(np.asarray(config["source"]["moment_tensor_Nm"], dtype=float),
                       np.asarray(metadata["moment_tensor_Nm"], dtype=float), rtol=1e-12, atol=0.0):
        raise ValueError("timing moment tensor differs from frozen tutorial metadata")
    window = config["window"]
    expected_window = {"type": metadata["time_window"]["type"],
                       "first_sample": metadata["time_window"]["sample_start_inclusive"],
                       "last_sample": metadata["time_window"]["sample_end_inclusive"],
                       "nt": metadata["time_window"]["sample_count"]}
    if any(window.get(key) != value for key, value in expected_window.items()):
        raise ValueError("timing window differs from frozen tutorial metadata")
    q = {"bin_1": metadata["q_nonzero"]["Z,1"],
         "bin_2": metadata["q_nonzero"]["Z,2"]}
    names = config["objective"]["names"]
    expected_names = ["J_Z", "J_N"][:measurement_count]
    if [item.get("name") for item in names] != expected_names:
        raise ValueError("timing objective list differs from requested measurements")
    if names[0].get("component") != "Z" or names[0].get("q_nonzero") != q:
        raise ValueError("timed J_Z weights differ from frozen tutorial metadata")
    if measurement_count == 2 and (names[1].get("component") != "N" or
                                   names[1].get("q_nonzero") != q):
        raise ValueError("timed J_N weights differ from frozen tutorial metadata")
    expected_angles = int(metadata["angular_nodes"])
    if config.get("angular_node_count") != expected_angles:
        raise ValueError("timed angular node count differs from frozen tutorial metadata")
    if sample_count != int(metadata["sample_count"]):
        raise ValueError("timed sample count differs from frozen tutorial metadata")
    if config.get("radial_node_count") * expected_angles != sample_count:
        raise ValueError("timed radial/angular counts do not match tutorial sample count")


def run_one(executable: Path, solver: str, args: argparse.Namespace,
            run_dir: Path, env: dict[str, str], started: float,
            model_hash: str, samples_hash: str) -> dict:
    run_dir.mkdir(parents=True, exist_ok=False)
    command = [str(executable.resolve()), str(args.model.resolve()),
               str(args.samples.resolve()), str(run_dir.resolve()), str(args.lmax),
               str(args.measurements), "1"]
    remaining = HARD_RUNTIME_S - (time.monotonic() - started)
    if remaining <= 0:
        raise TimeoutError("hard 900-second campaign limit reached")
    completed = subprocess.run(command, check=False, text=True, capture_output=True,
                               env=env, timeout=min(TIMEOUT_PER_CALL_S, remaining))
    (run_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, command,
                                            output=completed.stdout, stderr=completed.stderr)
    result = strict_json_loads(completed.stdout)
    ensure_finite_tree(result, f"{solver} timing output")
    expected_api_count = args.measurements if solver == "DSpecM1D" else 1
    if (result.get("solver") != solver or result.get("lmax") != args.lmax or
            result.get("measurements") != args.measurements or result.get("repeats") != 1 or
            result.get("warmup_batches") != 1 or
            result.get("warmup_api_calls") != (args.measurements if solver == "DSpecM1D" else 1) or
            len(result.get("api_s", [])) != expected_api_count):
        raise ValueError(f"timing executable metadata mismatch for {solver}: {result}")
    timing_keys = ("setup_s", "warmup_batch_s", "first_workload_api_s",
                   "first_model_to_arrays_s", "calculation_s", "end_to_end_s", "peak_rss_bytes")
    for key in timing_keys:
        value = float(result[key])
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"invalid {key} from {solver}: {value}")
    if any(not np.isfinite(float(value)) or float(value) < 0 for value in result["api_s"]):
        raise ValueError(f"invalid API timing returned by {solver}")
    if result["first_model_to_arrays_s"] < result["setup_s"]:
        raise ValueError(f"first model-to-arrays time is less than setup time for {solver}")
    kernels = run_dir / "kernels.csv"
    if not kernels.is_file() or kernels.stat().st_size == 0:
        raise ValueError(f"{solver} produced no timed kernel output")
    config = result.get("config")
    if not isinstance(config, dict):
        raise ValueError(f"{solver} did not emit resolved configuration metadata")
    result["peak_rss_label"] = (
        "diagnostic getrusage high-water RSS; may inherit launcher pre-exec high-water; "
        "not a solver-memory comparison")
    result["profiles_sha256"] = sha256(kernels)
    result["config"].setdefault("angular_node_count", max(8, 2 * args.lmax + 2) *
                                max(24, 4 * args.lmax + 4))
    result["config"].setdefault("radial_node_count", result.get("sample_count", result.get("point_count")) //
                                result["config"]["angular_node_count"])
    metadata = {"schema": "matched-kernel-timing-run-v1", "solver": solver,
                "lmax": args.lmax, "config": result["config"],
                "provenance": {"model_sha256": model_hash,
                               "samples_sha256": samples_hash,
                               "kernels_sha256": result["profiles_sha256"]}}
    (run_dir / "meta.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return result


def save_plot(raw: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    colors = {"DSpecM1D": "#386cb0", "YSpecpp": "#f0027f"}
    summary = {}
    for index, solver in enumerate(("DSpecM1D", "YSpecpp")):
        values = np.asarray([row["api_s"] for row in raw if row["solver"] == solver])
        median = float(np.median(values))
        low, high = float(values.min()), float(values.max())
        summary[solver] = {"api_s": values.tolist(), "median_s": median,
                           "min_s": low, "max_s": high,
                           "spread_s": high - low}
        ax.errorbar(index, median, yerr=[[median - low], [high - median]],
                    fmt="o", capsize=6, color=colors[solver])
        ax.scatter(np.full(values.size, index), values, color=colors[solver], alpha=.6, s=26)
    per_call = {}
    for solver in ("DSpecM1D", "YSpecpp"):
        rows = [row for row in raw if row["solver"] == solver]
        max_calls = max((len(row.get("api_s_by_objective", [])) for row in rows), default=0)
        per_call[solver] = {
            f"call_{index + 1}_s": [row["api_s_by_objective"][index]
                                    for row in rows if len(row.get("api_s_by_objective", [])) > index]
            for index in range(max_calls)}
    ax.set_xticks([0, 1], ["DSpecM1D", "YSpecpp"])
    ax.set_ylabel("Kernel API compute time (s, steady clock)")
    ax.set_title("Matched sampled-kernel API time; median and measured spread")
    ax.grid(axis="y", alpha=.25)
    for extension in ("png", "svg", "pdf"):
        fig.savefig(output / f"kernel_runtime.{extension}", dpi=300)
    plt.close(fig)
    with (output / "timings.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["replicate", "solver", "api_s", "api_s_by_objective", "first_workload_api_s",
                  "first_model_to_arrays_s", "calculation_s", "end_to_end_s", "setup_s",
                  "warmup_batch_s", "peak_rss_bytes", "peak_rss_label", "measurements", "point_count"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in raw:
            serialized = {key: row.get(key, "") for key in fields}
            serialized["api_s_by_objective"] = json.dumps(row.get("api_s_by_objective", []))
            writer.writerow(serialized)
    summary["api_s_by_objective_or_batch"] = per_call
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def write_report(raw: list[dict], output: Path, metadata: dict,
                 accuracy: dict | None) -> None:
    def summary(key: str) -> str:
        values = {solver: [float(row[key]) for row in raw if row["solver"] == solver]
                  for solver in ("DSpecM1D", "YSpecpp")}
        return "".join(
            f"<tr><th>{solver}</th><td>{np.median(items):.6g}</td>"
            f"<td>{min(items):.6g}–{max(items):.6g}</td></tr>"
            for solver, items in values.items())

    config = metadata.get("resolved_config", {})
    grid = config.get("frequency_grid", {})
    physics = config.get("physics", {})
    source = config.get("source", {})
    receiver = config.get("receiver", {})
    moment = source.get("moment_tensor_Nm", [])
    window = config.get("window", {})
    objectives = config.get("objective", {}).get("names", [])
    bins = grid.get("solver_bins_zero_based", [])
    df_hz = grid.get("df_hz")
    frequencies = ([f"{int(index) * float(df_hz):.7g} Hz "
                    f"({int(index) * float(df_hz) * 1000:.7g} mHz)"
                    for index in bins] if df_hz is not None else [])
    degrees = grid.get("degrees_inclusive", [])
    lmax = degrees[-1] if degrees else config.get("lmax", "not recorded")
    input_info = metadata.get("inputs", {})
    model_path = input_info.get("model") or "not recorded"
    model_name = Path(model_path).name
    model_identity = (f"{model_name} · SHA-256 {input_info.get('model_sha256', 'not recorded')}")
    point_counts = input_info.get("sample_row_counts_by_modulus", {})
    point_count_text = ", ".join(f"{key}: {value}" for key, value in point_counts.items())
    if not point_count_text:
        point_count_text = str(config.get("sample_count", "not recorded"))

    def location_text(value: dict) -> str:
        fields = (("geocentric_latitude_deg", "lat"), ("longitude_deg", "lon"),
                  ("depth_m", "depth m"))
        return ", ".join(f"{label} {value[key]}" for key, label in fields if key in value)

    measurement_text = ", ".join(item.get("name", "?") for item in objectives) or "not recorded"
    api_protocol = ("DSpecM1D: one scalar call per objective; YSpecpp: one batched call"
                    if len(objectives) == 2 else
                    "one J_Z workload; DSpecM1D scalar and YSpecpp batched")
    experiment_rows = [
        ("Degree limit (lmax)", lmax),
        ("Points per modulus", point_count_text),
        ("Objectives and API calls", f"{measurement_text}; {api_protocol}"),
        ("Frequency bins", ", ".join(f"{int(index)}: {frequency}"
                                      for index, frequency in zip(bins, frequencies)) or "not recorded"),
        ("Time grid and window", f"nt={grid.get('nt', 'n/a')}, dt={grid.get('dt_s', 'n/a')} s; "
         f"{window.get('type', 'n/a')} samples {window.get('first_sample', 'n/a')}–"
         f"{window.get('last_sample', 'n/a')} of {window.get('nt', 'n/a')}"),
        ("Source / receiver", f"{location_text(source) or 'not recorded'} / "
         f"{location_text(receiver) or 'not recorded'}"),
        ("Moment tensor", ", ".join(str(value) for value in moment) or "not recorded"),
        ("Gravity physics", f"{physics.get('gravity', 'not recorded')}; "
         f"G={physics.get('physical_G', 'not recorded')} m³ kg⁻¹ s⁻²"),
        ("Model identity", model_identity),
    ]
    experiment_html = "".join(
        f"<tr><th>{escape(str(label))}</th><td>{escape(str(value))}</td></tr>"
        for label, value in experiment_rows)

    if accuracy is None:
        gate_html = "<p>Accuracy gate results are unavailable in this saved report.</p>"
    else:
        rows = "".join(
            f"<tr><td>{item['objective']}</td><td>{item['modulus']}</td>"
            f"<td>{item['weighted_rms_difference']:.5g}</td>"
            f"<td>{item['relative_difference'] if item['relative_difference'] is not None else 'weak signal'}</td>"
            f"<td>PASS</td></tr>" for item in accuracy["profiles"])
        gate_html = ("<p>Actual timed arrays passed the DSpec–YSpec weighted profile and volume-integral gates. "
                     "DSpec timed J_Z values also passed comparison with the saved tutorial values.</p>"
                     "<table><thead><tr><th>Objective</th><th>Modulus</th><th>Weighted RMS difference</th>"
                     f"<th>Relative difference</th><th>Gate</th></tr></thead><tbody>{rows}</tbody></table>")
        integral_rows = "".join(
            f"<tr><td>{item['objective']}</td><td>{item['modulus']}</td>"
            f"<td>{item['I_D']:.5g}</td><td>{item['I_Y']:.5g}</td>"
            f"<td>{item['relative_difference'] if item['relative_difference'] is not None else 'weak signal'}</td>"
            f"<td>PASS</td></tr>" for item in accuracy["integrals"])
        gate_html += ("<h3>Volume integral checks</h3><table><thead><tr><th>Objective</th><th>Modulus</th>"
                      "<th>DSpec integral</th><th>YSpec integral</th><th>Relative difference</th><th>Gate</th>"
                      f"</tr></thead><tbody>{integral_rows}</tbody></table>")
    protocol = metadata.get("protocol", "")
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Matched kernel timing</title>
<style>body{{font:15px system-ui,sans-serif;max-width:980px;margin:2rem auto;padding:0 1rem;color:#222}}
table{{border-collapse:collapse;margin:1rem 0}}th,td{{border:1px solid #bbb;padding:.4rem .65rem;text-align:left}}
img{{max-width:100%}}code{{background:#eee;padding:.1rem .25rem}}</style></head><body>
<h1>Matched 3D kernel timing</h1>
<p>{protocol}</p>
<h2>Experiment configuration</h2>
<table><tbody>{experiment_html}</tbody></table>
<p>Each timing is specific to the recorded builds, model and exact sample grid. The process-order alternated by replicate.
API time excludes model loading, setup, warmup, and CSV export. Calculation time is setup plus the measured API calls;
first model-to-arrays is the observed setup plus the first cold workload call. The raw getrusage high-water RSS field
is retained as a diagnostic; it may inherit launcher pre-exec high-water and is not a solver-memory comparison.</p>
<table><thead><tr><th>Solver</th><th>Median (s)</th><th>Observed min–max (s)</th></tr></thead><tbody>
<tr><th colspan="3">Steady API workload</th></tr>{summary('api_s')}
<tr><th colspan="3">First model-to-arrays workload</th></tr>{summary('first_model_to_arrays_s')}
<tr><th colspan="3">Setup</th></tr>{summary('setup_s')}
<tr><th colspan="3">Setup + measured API calls</th></tr>{summary('calculation_s')}
<tr><th colspan="3">End-to-end process</th></tr>{summary('end_to_end_s')}
</tbody></table>
<h2>Numerical gates</h2>{gate_html}
<p>Profile metric: volume-weighted RMS with weights <code>dV × ΔX²</code>, using each modulus' own ΔX;
weak-signal cases are not accepted as relative passes. Integral gates use the matching <code>dV × ΔX</code> objective integral.</p>
<h2>Files</h2><ul><li><a href="timing.json">timing metadata and provenance</a></li>
<li><a href="build_info.json">compiler, flags, dependency and revision evidence</a></li>
<li><a href="timings.csv">raw timings</a></li><li><a href="accuracy_gate.json">gate details</a></li>
<li><a href="kernel_runtime.png">PNG plot</a>, <a href="kernel_runtime.svg">SVG</a>, <a href="kernel_runtime.pdf">PDF</a></li></ul>
<p>Environment and exact input hashes are in <a href="timing.json">timing.json</a>.</p>
<p>Dependency pins are recorded separately for each library build. Library-specific differences, including PlanetModel,
are retained as built; the timing workflow does not upgrade dependencies.</p>
</body></html>"""
    (output / "report.html").write_text(html, encoding="utf-8")


def replay(args: argparse.Namespace) -> None:
    source = args.replay.resolve()
    timings = read_csv(source / "timings.csv")
    raw = []
    for row in timings:
        raw.append({"replicate": int(row["replicate"]) if row.get("replicate") else "",
                    "solver": row["solver"], "api_s": float(row["api_s"]),
                    "api_s_by_objective": json.loads(row["api_s_by_objective"])
                        if row.get("api_s_by_objective") else [float(row["api_s"])],
                    "first_workload_api_s": float(row.get("first_workload_api_s", row["api_s"])),
                    "first_model_to_arrays_s": float(row["first_model_to_arrays_s"]
                                                       if row.get("first_model_to_arrays_s")
                                                       else row.get("calculation_s", row["end_to_end_s"])),
                    "calculation_s": float(row.get("calculation_s", row["end_to_end_s"])),
                    "end_to_end_s": float(row["end_to_end_s"]),
                    "setup_s": float(row["setup_s"]),
                    "warmup_batch_s": float(row.get("warmup_batch_s", 0)),
                    "peak_rss_bytes": int(row["peak_rss_bytes"]),
                    "peak_rss_label": ("diagnostic getrusage high-water RSS; may inherit launcher "
                                       "pre-exec high-water; not a solver-memory comparison"),
                    "measurements": int(row["measurements"]),
                    "point_count": int(row["point_count"])})
    if not raw:
        raise ValueError("saved timing data is empty")
    ensure_finite_tree(raw, "saved timings")
    args.output.mkdir(parents=True, exist_ok=True)
    save_plot(raw, args.output)
    saved_metadata = source / "timing.json"
    metadata = strict_json_loads(saved_metadata.read_text(encoding="utf-8")) if saved_metadata.is_file() else {}
    ensure_finite_tree(metadata, "saved timing metadata")
    if "resolved_config" not in metadata:
        first_run = source / "run_01_dspecm1d" / "meta.json"
        if first_run.is_file():
            sidecar = strict_json_loads(first_run.read_text(encoding="utf-8"))
            ensure_finite_tree(sidecar, "saved DSpec run metadata")
            metadata["resolved_config"] = sidecar.get("config", {})
    metadata["rendering"] = {"mode": "replay", "source": str(source),
                             "source_timing_sha256": sha256(saved_metadata) if saved_metadata.is_file() else None}
    (args.output / "timing.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    for name in ("build_info.json", "accuracy_gate.json"):
        original = source / name
        if original.is_file():
            shutil.copy2(original, args.output / name)
    accuracy = metadata.get("accuracy_gate")
    write_report(raw, args.output, metadata, accuracy)
    print(f"Replayed plot from saved timing data into {args.output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dspec", type=Path)
    parser.add_argument("--yspec", type=Path)
    parser.add_argument("--dspec-source", type=Path)
    parser.add_argument("--yspec-source", type=Path)
    parser.add_argument("--dspec-build", type=Path)
    parser.add_argument("--yspec-build", type=Path)
    parser.add_argument("--model", type=Path)
    parser.add_argument("--samples", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lmax", type=int, choices=(3, 12), default=3)
    parser.add_argument("--measurements", type=int, choices=(1, 2), default=1)
    parser.add_argument("--repeats", type=int, choices=(1, 3, 5), default=5)
    parser.add_argument("--accuracy-tolerance", type=float, default=1e-3)
    parser.add_argument("--replay", type=Path, help="replot saved timings.csv without solver calls")
    args = parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        raise SystemExit(f"refusing nonempty output directory: {args.output}")
    if args.replay:
        replay(args)
        return
    if not np.isfinite(args.accuracy_tolerance) or args.accuracy_tolerance <= 0:
        raise SystemExit("accuracy tolerance must be finite and positive")
    required = (args.dspec, args.yspec, args.dspec_source, args.yspec_source,
                args.dspec_build, args.yspec_build, args.model, args.samples)
    if any(value is None for value in required):
        parser.error("timing mode requires both executables, source/build dirs, model, and samples")
    for path in (args.dspec, args.yspec, args.dspec_source, args.yspec_source,
                 args.dspec_build, args.yspec_build, args.model, args.samples):
        if not path.exists():
            raise SystemExit(f"required input does not exist: {path}")
    model_hash, samples_hash = sha256(args.model), sha256(args.samples)
    tutorial_json = args.samples.parent / "run.json"
    tutorial_metadata = validate_tutorial_metadata(tutorial_json, args.lmax,
                                                   model_hash, samples_hash)
    sample_counts = validate_samples_csv(args.samples, int(tutorial_metadata["sample_count"]))
    tutorial_metadata_hash = sha256(tutorial_json)
    # Inspect both toolchains before creating any report output.
    dspec_build_info = build_info(args.dspec_source, args.dspec_build, args.dspec, "dspec_3d_kernel_timing")
    yspec_build_info = build_info(args.yspec_source, args.yspec_build, args.yspec,
                                  "yspec_kernel_timing", wrapper_source_root=args.dspec_source)
    if (Path(dspec_build_info["compiler"]).resolve() != Path(yspec_build_info["compiler"]).resolve() or
            dspec_build_info["compiler_version"] != yspec_build_info["compiler_version"]):
        raise SystemExit("DSpecM1D and YSpecpp timing builds do not use the same C++ compiler identity/version")
    env_info = thread_environment()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "build_info.json").write_text(json.dumps(
        {"DSpecM1D": dspec_build_info, "YSpecpp": yspec_build_info}, indent=2) + "\n",
        encoding="utf-8")
    raw = []
    accuracy = None
    started = time.monotonic()
    for replicate in range(args.repeats):
        order = (("DSpecM1D", args.dspec), ("YSpecpp", args.yspec)) if replicate % 2 == 0 else (
            ("YSpecpp", args.yspec), ("DSpecM1D", args.dspec))
        current = {}
        for solver, executable in order:
            run_dir = args.output / f"run_{replicate + 1:02d}_{solver.lower()}"
            current[solver] = run_one(executable, solver, args, run_dir,
                                      env_info["env"], started, model_hash, samples_hash)
            actual_count = int(current[solver].get("sample_count", current[solver].get("point_count", -1)))
            assert_config_matches_tutorial(current[solver]["config"], tutorial_metadata,
                                           args.lmax, args.measurements, actual_count)
        assert_configs_match(current["DSpecM1D"]["config"], current["YSpecpp"]["config"])
        if replicate == 0:
            accuracy = accuracy_gate(args.output / "run_01_dspecm1d" / "kernels.csv",
                                     args.output / "run_01_yspecpp" / "kernels.csv",
                                     args.samples, args.accuracy_tolerance,
                                     args.measurements, current["DSpecM1D"], tutorial_metadata)
            (args.output / "accuracy_gate.json").write_text(
                json.dumps(accuracy, indent=2) + "\n", encoding="utf-8")
        for solver in ("DSpecM1D", "YSpecpp"):
            result = current[solver]
            # DSpec evaluates each scalar q separately; sum those API calls to
            # compare the same two-objective workload as YSpec's batched call.
            api_time = float(sum(result["api_s"]))
            raw.append({"replicate": replicate + 1, "solver": solver,
                        "api_s": api_time, "calculation_s": result["calculation_s"],
                        "api_s_by_objective": [float(value) for value in result["api_s"]],
                        "first_workload_api_s": result["first_workload_api_s"],
                        "first_model_to_arrays_s": result["first_model_to_arrays_s"],
                        "warmup_batch_s": result["warmup_batch_s"],
                        "end_to_end_s": result["end_to_end_s"],
                        "setup_s": result["setup_s"],
                        "peak_rss_bytes": result["peak_rss_bytes"],
                        "peak_rss_label": result["peak_rss_label"],
                        "measurements": args.measurements,
                        "point_count": result.get("sample_count", result.get("point_count"))})
    save_plot(raw, args.output)
    metadata = {"protocol": "shared tutorial sample grid; homogeneous model; full gravity; bins 1-3; one or two real windowed SI objectives; no source contacts; volume API only; each process performs one warmup batch. DSpecM1D makes one timed scalar API call per objective; YSpecpp makes one timed batched call.",
               "accuracy_gate": accuracy,
               "timing_summary": json.loads((args.output / "summary.json").read_text()),
               "timing_order": "alternating executable process order",
               "resource_limits": {"per_process_timeout_s": TIMEOUT_PER_CALL_S,
                                   "hard_campaign_timeout_s": HARD_RUNTIME_S,
                                   "thread_environment": env_info["limits"]},
               "environment": {"host": platform.node(), "platform": platform.platform(),
                               "processor": platform.processor(),
                               "machine": platform.machine(), "cpu_count": os.cpu_count(),
                               "python": platform.python_version(), "numpy": np.__version__,
                               "matplotlib": matplotlib.__version__},
               "inputs": {"model": str(args.model.resolve()), "model_sha256": model_hash,
                          "samples": str(args.samples.resolve()), "samples_sha256": samples_hash,
                          "tutorial_metadata": str(tutorial_json.resolve()),
                          "tutorial_metadata_sha256": tutorial_metadata_hash,
                          "sample_row_counts_by_modulus": sample_counts},
               "resolved_config": current["DSpecM1D"]["config"],
               "timing_definitions": {"api_s": "sum of measured API calls; excludes setup, warmup, CSV export",
                                      "first_model_to_arrays_s": "actual setup plus the first cold workload API call",
                                      "calculation_s": "setup plus all measured API calls; excludes warmup and export",
                                      "end_to_end_s": "whole process interval including warmup and output export"},
               "build_info": {"DSpecM1D": dspec_build_info, "YSpecpp": yspec_build_info},
               "library_specific_dependency_pins": {
                   "DSpecM1D": dspec_build_info["selected_dependency_pins"],
                   "YSpecpp": yspec_build_info["selected_dependency_pins"]},
               "runtime_warning": "Results describe these builds, this model, and this shared sampling grid; they do not imply a universal speed ranking."}
    (args.output / "timing.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                             encoding="utf-8")
    write_report(raw, args.output, metadata, accuracy)
    print(f"Timed-output accuracy gate passed; results saved to {args.output}")


if __name__ == "__main__":
    main()
