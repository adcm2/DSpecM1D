#!/usr/bin/env python3
"""Run the paper-matched high-frequency Fig. 5 DSpecM1D record section."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

import run_paper_suite as suite

BASE = Path(__file__).resolve().parent
DSPEC_DEFAULT = BASE.parents[1]
CASE_FILE = BASE / "cases/paper-suite-cases.json"
CRITERIA_DEFAULT = BASE / "cases/paper-suite-protocol.json"
OUTPUT_DEFAULT = DSPEC_DEFAULT / "build/phase3/paper-suite-20261020"
UNIT_BY_QUANTITY = {"displacement": "m", "velocity": "m/s", "acceleration": "m/s^2"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_all_receiver_outputs(series_path: Path, spectra_path: Path,
                                 meta_path: Path, receivers: list[list[float]],
                                 settings: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    """Validate the complete all-receiver SI record and spectral grid."""
    meta = json.loads(meta_path.read_text())
    expected_series_header = ["time_s"] + [
        f"r{receiver}_{component}_SI"
        for receiver in range(len(receivers)) for component in ("Z", "N", "E")
    ]
    with series_path.open(newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, None)
        if header != expected_series_header:
            raise SystemExit(f"all-receiver time header mismatch: {header}")
        records = list(reader)
    dt = float(settings["time_step_s"])
    expected_samples = int(round(settings["duration_min"] * 60.0 / dt)) + 1
    if len(records) != expected_samples or expected_samples != int(meta["output_samples"]):
        raise SystemExit("all-receiver sample count does not match the inclusive duration grid")
    try:
        series_values = np.asarray([[float(value) for value in row] for row in records], dtype=float)
    except (TypeError, ValueError) as error:
        raise SystemExit("all-receiver time series contains a malformed value") from error
    if series_values.shape != (expected_samples, len(expected_series_header)) or not np.isfinite(series_values).all():
        raise SystemExit("all-receiver time series has invalid shape or non-finite SI values")
    time_s = series_values[:, 0]
    expected_time = np.arange(expected_samples, dtype=float) * dt
    if (not np.all(np.diff(time_s) > 0) or time_s[0] != 0.0 or
            not np.allclose(time_s, expected_time, rtol=0.0, atol=1e-10) or
            not math.isclose(time_s[-1], settings["duration_min"] * 60.0,
                             rel_tol=0.0, abs_tol=1e-10)):
        raise SystemExit("all-receiver time grid is not monotone from zero through the inclusive duration")

    expected_spectra_header = ["frequency_bin", "frequency_hz", "omega_real_rad_s",
                               "omega_imag_rad_s", "receiver_index", "component",
                               "real", "imag"]
    first_bin = int(meta["frequency_first_index"])
    last_bin = int(meta["frequency_last_index_exclusive"])
    expected_keys = {(index, receiver, component)
                     for index in range(first_bin, last_bin)
                     for receiver in range(len(receivers)) for component in range(3)}
    seen = set()
    frequency_by_bin = {}
    df = float(meta["frequency_spacing_hz"])
    damping = float(meta["damping_epsilon_rad_s"])
    with spectra_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_spectra_header:
            raise SystemExit(f"all-receiver spectrum header mismatch: {reader.fieldnames}")
        for row in reader:
            try:
                index = int(row["frequency_bin"])
                receiver = int(row["receiver_index"])
                component = int(row["component"])
                frequency = float(row["frequency_hz"])
                omega_real = float(row["omega_real_rad_s"])
                omega_imag = float(row["omega_imag_rad_s"])
                real = float(row["real"])
                imag = float(row["imag"])
            except (TypeError, ValueError, KeyError) as error:
                raise SystemExit("all-receiver spectrum contains a malformed value") from error
            key = (index, receiver, component)
            if key not in expected_keys or key in seen:
                raise SystemExit(f"unexpected or duplicate receiver/component/bin tuple: {key}")
            if not all(math.isfinite(value) for value in
                       (frequency, omega_real, omega_imag, real, imag)):
                raise SystemExit(f"non-finite all-receiver spectrum value at {key}")
            if not math.isclose(frequency, index * df, rel_tol=1e-12, abs_tol=1e-15):
                raise SystemExit(f"frequency grid mismatch at bin {index}")
            if not math.isclose(omega_real, 2.0 * math.pi * frequency,
                                rel_tol=1e-10, abs_tol=1e-14):
                raise SystemExit(f"real angular-frequency grid mismatch at bin {index}")
            if not math.isclose(omega_imag, -damping,
                                rel_tol=1e-12, abs_tol=1e-15):
                raise SystemExit(f"imaginary angular-frequency grid mismatch at bin {index}")
            values = (frequency, omega_real, omega_imag)
            previous = frequency_by_bin.setdefault(index, values)
            if previous != values:
                raise SystemExit(f"inconsistent spectral grid across receivers at bin {index}")
            seen.add(key)
    if seen != expected_keys:
        raise SystemExit(f"all-receiver spectra are incomplete: {len(seen)} of {len(expected_keys)} tuples")
    return time_s, series_values[:, 1:], {
        "time_header": expected_series_header,
        "time_samples": expected_samples,
        "time_step_s": dt,
        "time_start_s": float(time_s[0]),
        "time_end_s": float(time_s[-1]),
        "spectral_header": expected_spectra_header,
        "frequency_bins": last_bin - first_bin,
        "receivers": len(receivers),
        "components": 3,
        "spectral_tuples": len(seen),
        "grid_validated": True,
        "all_si_values_finite": True,
    }


def plot_record_section(path: Path, time_s: np.ndarray, traces: np.ndarray,
                        receiver_lons: list[float], component: str,
                        unit: str) -> list[float]:
    import matplotlib.pyplot as plt

    rows = 0 if component == "Z" else 1
    data = traces[:, rows, :]
    peaks = np.max(np.abs(data), axis=1)
    safe = np.where(peaks > 0, peaks, 1.0)
    normalized = data / safe[:, None]
    fig, axis = plt.subplots(figsize=(12, 6))
    mesh = axis.pcolormesh(time_s, receiver_lons, normalized, shading="auto",
                           cmap="seismic", vmin=-1, vmax=1, rasterized=True)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Equatorial angular distance (degrees)")
    axis.set_title(f"Figure 5 adopted acceleration {component} record section")
    fig.colorbar(mesh, ax=axis, label=f"Normalized by each receiver peak ({unit})")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return [float(value) for value in peaks]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true",
                        help="use the top-frequency sentinel and three receivers")
    parser.add_argument("--diagnostic-max-degree", type=int,
                        help="optional explicitly non-paper pilot degree cap")
    parser.add_argument("--dspec-source", type=Path, default=DSPEC_DEFAULT)
    parser.add_argument("--dependency-root", type=Path,
                        default=suite.DEPENDENCY_DEFAULT)
    parser.add_argument("--criteria", type=Path, default=CRITERIA_DEFAULT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args()
    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"

    case = json.loads(CASE_FILE.read_text())["cases"]["fig5"]
    criteria = json.loads(args.criteria.read_text())
    if "ex4/Fig5" not in {entry["case"] for entry in criteria["case_inventory"]}:
        raise SystemExit("ex4/Fig5 is absent from the frozen criteria")
    input_path = args.dspec_source / case["source_parameter_file"]
    expected_hash = criteria["source_basis"]["source_parameter_hashes"]["ex4"]
    if sha256(input_path) != expected_hash:
        raise SystemExit(f"source parameter hash differs from reviewed criteria: {input_path}")

    settings = dict(case)
    if args.diagnostic_max_degree is not None:
        if not args.pilot or not (1 <= args.diagnostic_max_degree < case["degree_max"]):
            raise SystemExit("--diagnostic-max-degree requires --pilot and a cap below 7500")
        settings["degree_max"] = args.diagnostic_max_degree
        settings["diagnostic_degree_cap"] = True
    if args.pilot:
        settings["frequency_min_mHz"], settings["frequency_max_mHz"] = case["pilot_frequency_band_mHz"]
        settings["receivers_geocentric"] = [case["receivers_geocentric"][i] for i in (0, 45, 90)]
        settings["run_scope"] = ("non-paper diagnostic capped below required degree; not paper result"
                                  if args.diagnostic_max_degree is not None else
                                  "bounded top-frequency, three-receiver pilot; not paper result")
    else:
        settings["run_scope"] = "full published-style Figure 5 record section"

    output_root = args.output_root.resolve()
    case_root = output_root / "fig5" / ("pilot" if args.pilot else "full")
    if case_root.exists() and any(case_root.iterdir()):
        raise SystemExit(f"refusing to overwrite existing run artifacts: {case_root}")
    build = output_root / "build-dspec-fig5"
    case_root.mkdir(parents=True, exist_ok=True)
    model_source = args.dspec_source / case["model_file"]
    model_copy = case_root / "inputs/data/models" / model_source.name
    model_copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(model_source, model_copy)
    result_dir = case_root / "results"
    result_dir.mkdir(parents=True, exist_ok=True)
    params_path = case_root / "inputs/data/params/paper-suite-fig5.txt"
    params_path.parent.mkdir(parents=True, exist_ok=True)
    prefix = result_dir / "fig5"
    params_path.write_text(suite.make_dspec_parameters(model_copy.name, prefix, settings))

    suite.run(["cmake", "-S", str(BASE), "-B", str(build),
               "-DCMAKE_BUILD_TYPE=Release", "-DBENCHMARK_SOLVER=dspec",
               f"-DDSPEC_SOURCE_DIR={args.dspec_source.resolve()}",
               f"-DVERIFIED_DEPENDENCY_ROOT={args.dependency_root.resolve()}"])
    suite.run(["cmake", "--build", str(build), "--target", "dspec_adapter", "--parallel", "2"])
    binary = build / "dspec_adapter"
    snapshot = case_root / "source_snapshot"
    snapshot.mkdir(exist_ok=True)
    source_paths = (BASE / "run_fig5_record.py", BASE / "run_paper_suite.py",
                    BASE / "run_paper2.py", BASE / "dspec_adapter.cpp",
                    BASE / "CMakeLists.txt", CASE_FILE)
    source_hashes = {}
    for source in source_paths:
        shutil.copyfile(source, snapshot / source.name)
        source_hashes[source.name] = sha256(snapshot / source.name)
    shutil.copyfile(binary, snapshot / binary.name)

    adapter_command = [str(binary), str(params_path), str(prefix),
                       str(case["dspec"]["nq"]), "10", "0.05",
                       "paper_suite", "dynamic"]
    try:
        wall, rss_kib = suite.run_adapter(
            adapter_command, result_dir / "dspec.stdout.log",
            result_dir / "dspec.resource.txt")
    except subprocess.CalledProcessError as error:
        adapter_output = (result_dir / "dspec.stdout.log").read_text()
        if "lMax =" in adapter_output and "refused" in adapter_output:
            status = "solver_rejected_by_degree_safety_guard"
        else:
            status = "solver_process_failed; calculation_status_unknown"
        failure = {
            "case": "fig5",
            "run_scope": settings["run_scope"],
            "status": status,
            "return_code": error.returncode,
            "command": adapter_command,
            "adapter_output": adapter_output,
            "settings": settings,
            "source_parameter_sha256": sha256(input_path),
            "model_sha256": sha256(model_source),
            "adapter_binary_sha256": sha256(snapshot / binary.name),
            "harness_source_sha256": source_hashes,
            "resource_log": str(result_dir / "dspec.resource.txt"),
        }
        (case_root / "failure.json").write_text(json.dumps(failure, indent=2) + "\n")
        if status == "solver_rejected_by_degree_safety_guard":
            message = "Fig. 5 stopped at the degree safety guard; see failure.json"
        else:
            message = "Fig. 5 adapter process failed; calculation stage is unknown; see failure.json"
        raise SystemExit(message) from error
    all_series_path = prefix.with_name(prefix.name + ".all_receivers.timeseries.csv")
    all_spectra_path = prefix.with_name(prefix.name + ".all_receivers.spectra.csv")
    if not all_series_path.is_file() or not all_spectra_path.is_file():
        raise SystemExit("adapter did not produce the expected all-receiver output")
    time_s, values, output_validation = validate_all_receiver_outputs(
        all_series_path, all_spectra_path, prefix.with_suffix(".meta.json"),
        settings["receivers_geocentric"], settings)
    traces = values.reshape(len(time_s), len(settings["receivers_geocentric"]), 3).transpose(1, 2, 0)
    lons = [location[1] for location in settings["receivers_geocentric"]]
    unit = UNIT_BY_QUANTITY[settings["quantity"]]
    peaks_z = plot_record_section(case_root / "record_Z.png", time_s, traces, lons, "Z", unit)
    peaks_n = plot_record_section(case_root / "record_N.png", time_s, traces, lons, "N", unit)

    lines = (result_dir / "dspec.stdout.log").read_text().splitlines()
    maxsteps = [float(line.split(":", 1)[1]) for line in lines if line.startswith("Max step:")]
    provenance = {
        "case": "fig5",
        "paper_figure": 5,
        "run_scope": settings["run_scope"],
        "criteria_sha256": sha256(args.criteria),
        "case_table_sha256": sha256(CASE_FILE),
        "source_parameter_file": str(input_path),
        "source_parameter_sha256": sha256(input_path),
        "model_file": str(model_source),
        "model_sha256": sha256(model_source),
        "DSpecM1D_HEAD_at_run": subprocess.check_output(
            ["git", "-C", str(args.dspec_source), "rev-parse", "HEAD"], text=True).strip(),
        "harness_source_sha256": source_hashes,
        "harness_source_snapshot": "source_snapshot/",
        "adapter_binary_sha256": sha256(snapshot / binary.name),
        "settings": settings,
        "quantity_and_units": {"quantity": settings["quantity"], "unit": unit,
                                "paper_quantity_explicit": False},
        "source_time": "Heaviside step moment; no extra source filter",
        "coordinate_convention": "equatorial geocentric locations; 0 to 180 degrees in 2-degree increments",
        "mesh_strategy": "dynamic frequency-chunked multi-SEM",
        "actual_maxstep_per_frequency_chunk": maxsteps,
        "runtime_environment": {"OMP_NUM_THREADS": "2", "OMP_DYNAMIC": "FALSE"},
        "runtime_s": wall,
        "peak_child_rss_kib": rss_kib,
        "output_validation": output_validation,
        "arrival_markers": "Not overlaid; no independent TauP marker dataset is present in the frozen case inputs.",
    }
    (case_root / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    summary = {
        "case": "fig5",
        "paper_figure": 5,
        "run_scope": settings["run_scope"],
        "quantity": settings["quantity"],
        "unit": unit,
        "record_geometry": {"receiver_count": len(lons), "longitude_deg": lons,
                            "latitude_geocentric_deg": 0.0},
        "grid": {"time_samples": len(time_s), "time_step_s": settings["time_step_s"],
                 "duration_s": float(time_s[-1] - time_s[0]),
                 "all_finite": True, "all_receiver_spectra_csv": all_spectra_path.name,
                 "validation": output_validation},
        "peak_amplitude_by_receiver_SI": {"Z": peaks_z, "N": peaks_n},
        "runtime": {"seconds": wall, "peak_rss_kib": rss_kib},
        "record_plots": ["record_Z.png", "record_N.png"],
        "provenance": "provenance.json",
        "scientific_scope_note": "No cross-solver error score is reported; acceleration is adopted from ex4, as Fig. 5 does not name the observable.",
    }
    (case_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"case_output": str(case_root), "runtime_s": wall,
                      "peak_rss_kib": rss_kib, "receivers": len(lons),
                      "samples": len(time_s)}, indent=2))


if __name__ == "__main__":
    main()
