#!/usr/bin/env python3
"""Run one paper-matched Fig. 1, 3, or 4 comparison case."""

from __future__ import annotations

import argparse
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

import run_paper2 as common

BASE = Path(__file__).resolve().parent
DSPEC_DEFAULT = BASE.parents[1]
YSPEC_DEFAULT = common.YSPEC_DEFAULT
DEPENDENCY_DEFAULT = common.DEPENDENCY_DEFAULT
CASE_FILE = BASE / "cases/paper-suite-cases.json"
CRITERIA_DEFAULT = BASE / "cases/paper-suite-protocol.json"
UNIT_BY_QUANTITY = {
    "displacement": "m",
    "velocity": "m/s",
    "acceleration": "m/s^2",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], log: Path | None = None) -> float:
    started = time.perf_counter()
    if log is None:
        subprocess.run(command, check=True)
    else:
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(completed.stdout)
        if completed.returncode:
            print(completed.stdout, end="", flush=True)
            raise subprocess.CalledProcessError(completed.returncode, command)
    return time.perf_counter() - started


def run_adapter(command: list[str], log: Path, resource_log: Path) -> tuple[float, int | None]:
    timer = Path("/usr/bin/time")
    if not timer.is_file():
        return run(command, log), None
    elapsed = run([str(timer), "-f", "MAX_RSS_KIB=%M", "-o",
                   str(resource_log), *command], log)
    line = resource_log.read_text().strip()
    if not line.startswith("MAX_RSS_KIB="):
        raise SystemExit(f"could not read peak RSS from {resource_log}")
    return elapsed, int(line.split("=", 1)[1])


def settings_for(case: dict, pilot: bool) -> dict:
    settings = dict(case)
    if pilot:
        lo, hi = case["pilot_frequency_band_mHz"]
        settings["frequency_min_mHz"] = lo
        settings["frequency_max_mHz"] = hi
        settings["run_scope"] = "bounded top-frequency pilot; not paper result"
    else:
        settings["run_scope"] = "full published-style band"
    return settings


def make_dspec_parameters(model_name: str, prefix: Path, settings: dict) -> str:
    m = settings["moment_tensor_Nm"]
    gravity = {"none": 0, "cowling": 1, "full": 2}[settings["gravity"]]
    output = {"displacement": 0, "velocity": 1, "acceleration": 2}[
        settings["quantity"]]
    receivers = settings["receivers_geocentric"]
    values = [
        str(prefix), f"models/{model_name}", "4",
        "1" if settings["attenuation"] else "0", str(gravity), str(output),
        "0", f"{settings['dspec']['relative_error']:.17g}",
        str(settings["degree_min"]), str(settings["degree_max"]),
        f"{settings['frequency_min_mHz']:.17g}",
        f"{settings['frequency_max_mHz']:.17g}",
        f"{settings['duration_min']:.17g}", f"{settings['time_step_s']:.17g}",
        *(f"{corner:.17g}" for corner in settings["filter_mHz"]),
        f"{settings['source_depth_km']:.17g}",
        f"{settings['source_geocentric_latitude_deg']:.17g}",
        f"{settings['source_longitude_deg']:.17g}",
        *(f"{m[key]:.17g}" for key in ("rr", "rt", "rp", "tt", "tp", "pp")),
        f"{settings['receiver_depth_km']:.17g}", str(len(receivers)),
        *(f"{lat:.17g} {lon:.17g}" for lat, lon in receivers),
    ]
    return "\n".join(values) + "\n"


def make_yspec_input(model: Path, settings: dict) -> dict:
    m = settings["moment_tensor_Nm"]
    receivers = settings["receivers_geocentric"]
    return {
        "model": str(model.resolve()),
        "attenuation": settings["attenuation"],
        "gravity": settings["gravity"],
        "output": {"prefix": "unused", "quantity": settings["quantity"],
                   "potential_and_tilt_corrections": False},
        "degrees": {"min": settings["degree_min"], "max": settings["degree_max"]},
        "frequencies_mHz": {"min": settings["frequency_min_mHz"],
                             "max": settings["frequency_max_mHz"]},
        "time_series": {"length_min": settings["duration_min"],
                        "time_step_s": settings["time_step_s"]},
        "filter_mHz": settings["filter_mHz"],
        "source": {
            "depth_km": settings["source_depth_km"],
            "latitude": common.geographic_latitude(
                settings["source_geocentric_latitude_deg"]),
            "longitude": settings["source_longitude_deg"],
            "moment_tensor_Nm": m,
        },
        "receivers": {
            "depth_km": settings["receiver_depth_km"],
            "locations": [[common.geographic_latitude(lat), lon]
                          for lat, lon in receivers],
        },
        "solver": {"tolerance": settings["yspec"]["tolerance"],
                   "start_decay": settings["yspec"]["start_decay"]},
    }


def compare_series(reference: np.ndarray, candidate: np.ndarray,
                   unit: str) -> dict:
    diff = candidate - reference
    peak = float(np.max(np.abs(reference)))
    norm = float(np.linalg.norm(reference))
    maxdiff = float(np.max(np.abs(diff)))
    return {
        "reference_peak": peak,
        "unit": unit,
        "mean_absolute_difference": float(np.mean(np.abs(diff))),
        "mean_difference_percent_of_reference_peak":
            100 * float(np.mean(np.abs(diff))) / peak if peak else None,
        "maximum_absolute_difference": maxdiff,
        "maximum_difference_percent_of_reference_peak":
            100 * maxdiff / peak if peak else None,
        "relative_l2": float(np.linalg.norm(diff) / norm) if norm else None,
    }


def normalize_inverse_check_units(checks: list[dict], unit: str) -> list[dict]:
    """Replace the inherited displacement-specific diagnostic suffixes."""
    converted = []
    for check in checks:
        row = {key.replace("_m", "_SI"): value for key, value in check.items()}
        row["unit"] = unit
        converted.append(row)
    return converted


def compare_spectrum_amplitude(reference: dict, candidate: dict,
                               points: dict, band: tuple[float, float],
                               component: int) -> dict:
    keys = [key for key in sorted(reference)
            if key[1] == component and
            band[0] / 1000 <= points[key]["frequency_hz"] <= band[1] / 1000]
    if not keys:
        raise SystemExit("no spectral samples in requested comparison band")
    ref = np.asarray([abs(reference[key]) for key in keys])
    cand = np.asarray([abs(candidate[key]) for key in keys])
    delta = cand - ref
    peak = float(np.max(ref))
    complex_ref = np.asarray([reference[key] for key in keys])
    complex_cand = np.asarray([candidate[key] for key in keys])
    return {
        "quantity": "raw complex spectrum amplitude",
        "unit": "SI output unit / Hz",
        "frequency_band_mHz": list(band),
        "frequency_bins": len(keys),
        "mean_absolute_amplitude_difference": float(np.mean(np.abs(delta))),
        "mean_amplitude_difference_percent_of_reference_peak":
            100 * float(np.mean(np.abs(delta))) / peak if peak else None,
        "maximum_absolute_amplitude_difference": float(np.max(np.abs(delta))),
        "maximum_amplitude_difference_percent_of_reference_peak":
            100 * float(np.max(np.abs(delta))) / peak if peak else None,
        "complex_relative_l2_diagnostic":
            float(np.linalg.norm(complex_cand - complex_ref) /
                  np.linalg.norm(complex_ref)) if np.linalg.norm(complex_ref) else None,
    }


def legacy_processed_amplitude_spectrum(spectra: dict, meta: dict,
                                       filter_corners: list[float],
                                       window_fraction: float) -> tuple[np.ndarray, np.ndarray]:
    """Apply ex1's adopted pre-window/filter/final-window chain to Z traces."""
    nfft = int(meta["samples"])
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    output_samples = int(meta["output_samples"])
    frequency = np.arange(nfft // 2 + 1, dtype=float) * df
    bins = np.zeros(nfft // 2 + 1, dtype=np.complex128)
    for (index, component), value in spectra.items():
        if component == 0:
            bins[index] = value

    time = np.arange(nfft, dtype=float) * dt
    duration = max((output_samples - 1) * dt, dt)

    def hann_ends(fraction: float) -> np.ndarray:
        edge = fraction * duration
        window = np.ones(nfft)
        if edge <= 0:
            return window
        left = time < edge
        right = time > duration - edge
        window[left] = 0.5 * (1.0 - np.cos(np.pi * time[left] / edge))
        window[right] = 0.5 * (1.0 - np.cos(np.pi * (duration - time[right]) / edge))
        window[time >= duration] = 0.0
        return window

    raw_time = np.fft.irfft(bins, n=nfft) * nfft * df * np.exp(epsilon * time)
    pre_spectrum = np.fft.rfft(raw_time * hann_ends(0.05), n=nfft)
    filtered_time = np.fft.irfft(
        pre_spectrum * common.half_cosine(frequency, filter_corners), n=nfft
    ) * nfft * df
    final_spectrum = np.fft.rfft(filtered_time * hann_ends(window_fraction), n=nfft)
    return frequency, np.abs(final_spectrum)


def plot_timeseries(path: Path, time: np.ndarray, d: np.ndarray, y: np.ndarray,
                    quantity: str, window: tuple[float, float] | None,
                    metrics: dict) -> None:
    import matplotlib.pyplot as plt

    selector = np.ones(time.shape, dtype=bool)
    if window:
        selector = (time >= window[0]) & (time <= window[1])
    labels = ("Z", "N", "E")
    unit = UNIT_BY_QUANTITY[quantity]
    fig, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    for component, label in enumerate(labels):
        peak = metrics[label]["reference_peak"] or 1.0
        axes[component, 0].plot(time[selector], d[component, selector] / peak,
                                label="Current DSpecM1D", linewidth=0.9)
        axes[component, 0].plot(time[selector], y[component, selector] / peak,
                                label="Pinned YSpecpp", linestyle="--", linewidth=0.9)
        axes[component, 0].set_ylabel(f"{label} / YSpecpp peak")
        axes[component, 1].plot(time[selector],
                                100 * (d[component, selector] - y[component, selector]) / peak,
                                color="black", linewidth=0.8)
        axes[component, 1].set_ylabel("difference (% of peak)")
        axes[component, 0].grid(alpha=0.25)
        axes[component, 1].grid(alpha=0.25)
    axes[0, 0].set_title(f"Filtered {quantity}")
    axes[0, 1].set_title("DSpecM1D − YSpecpp")
    axes[2, 0].set_xlabel("Time (s)")
    axes[2, 1].set_xlabel("Time (s)")
    axes[0, 0].legend(loc="upper right")
    fig.supylabel(f"{quantity} / YSpecpp peak ({unit})")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def plot_spectrum(path: Path, points: dict, dspec: dict, yspec: dict,
                  band: tuple[float, float], unit: str) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for component, label in ((0, "Z"),):
        keys = sorted(key for key in dspec if key[1] == component and
                      band[0] / 1000 <= points[key]["frequency_hz"] <= band[1] / 1000)
        f = np.asarray([points[key]["frequency_hz"] * 1000 for key in keys])
        d = np.asarray([abs(dspec[key]) for key in keys])
        y = np.asarray([abs(yspec[key]) for key in keys])
        axes[0].plot(f, d, label="Current DSpecM1D")
        axes[0].plot(f, y, "--", label="Pinned YSpecpp")
        peak = np.max(y)
        axes[1].plot(f, 100 * (d - y) / peak if peak else d-y, color="black")
        axes[0].set_ylabel(f"|Z| ({unit})")
        axes[1].set_ylabel("amplitude difference (% of YSpecpp band peak)")
    axes[0].set_yscale("log")
    axes[0].legend()
    axes[1].set_xlabel("Frequency (mHz)")
    axes[0].grid(alpha=0.25)
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("fig1", "fig3", "fig4"), required=True)
    parser.add_argument("--pilot", action="store_true",
                        help="run the top-frequency sentinel band only")
    parser.add_argument("--dspec-source", type=Path, default=DSPEC_DEFAULT)
    parser.add_argument("--yspec-source", type=Path, default=YSPEC_DEFAULT)
    parser.add_argument("--dependency-root", type=Path, default=DEPENDENCY_DEFAULT)
    parser.add_argument("--criteria", type=Path, default=CRITERIA_DEFAULT)
    parser.add_argument("--output-root", type=Path,
                        default=DSPEC_DEFAULT / "build/phase3/paper-suite-20261020")
    parser.add_argument("--build-root", type=Path,
                        help="optional isolated build directory root (for concurrent cases)")
    args = parser.parse_args()
    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"

    case_table = json.loads(CASE_FILE.read_text())["cases"]
    case = case_table[args.case]
    criteria = json.loads(args.criteria.read_text())
    reviewed = {entry["case"] for entry in criteria["case_inventory"]}
    expected_case = {"fig1": "ex1/Fig1", "fig3": "ex7/Fig3", "fig4": "ex3/Fig4"}[args.case]
    if expected_case not in reviewed:
        raise SystemExit(f"{expected_case} is absent from the frozen criteria")

    input_path = args.dspec_source / case["source_parameter_file"]
    expected_hash = criteria["source_basis"]["source_parameter_hashes"][
        Path(case["source_parameter_file"]).stem]
    if sha256(input_path) != expected_hash:
        raise SystemExit(f"source parameter hash differs from reviewed criteria: {input_path}")
    if not args.yspec_source.is_dir():
        raise SystemExit(f"missing pinned YSpecpp source: {args.yspec_source}")
    common.verify_dependencies(args.dependency_root)
    yspec_tree, yspec_count = common.filetree_sha256(args.yspec_source)
    if yspec_tree != common.YSPEC_FILETREE_SHA256 or yspec_count != 76:
        raise SystemExit("YSpecpp source tree differs from its pinned snapshot")

    settings = settings_for(case, args.pilot)
    suite_root = args.output_root.resolve()
    case_root = suite_root / args.case / ("pilot" if args.pilot else "full")
    if case_root.exists() and any(case_root.iterdir()):
        raise SystemExit(f"refusing to overwrite existing run artifacts: {case_root}")
    build_root = (args.build_root.resolve() if args.build_root else suite_root)
    build_d, build_y = build_root / "build-dspec", build_root / "build-yspecpp"
    model_source = args.dspec_source / case["model_file"]
    model_hash = sha256(model_source)
    model_copy = case_root / "inputs/data/models" / model_source.name
    model_copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(model_source, model_copy)
    result_dir = case_root / "results"
    result_dir.mkdir(parents=True, exist_ok=True)
    params_path = case_root / "inputs/data/params/paper-suite.txt"
    params_path.parent.mkdir(parents=True, exist_ok=True)
    dspec_prefix = result_dir / "dspec"
    params_path.write_text(make_dspec_parameters(model_copy.name,
                                                  dspec_prefix, settings))
    yinput = case_root / "inputs/yspec.json"
    yinput.write_text(json.dumps(make_yspec_input(model_copy, settings), indent=2) + "\n")

    run(["cmake", "-S", str(BASE), "-B", str(build_d),
         "-DCMAKE_BUILD_TYPE=Release", "-DBENCHMARK_SOLVER=dspec",
         f"-DDSPEC_SOURCE_DIR={args.dspec_source.resolve()}",
         f"-DVERIFIED_DEPENDENCY_ROOT={args.dependency_root.resolve()}"])
    run(["cmake", "--build", str(build_d), "--target", "dspec_adapter", "--parallel", "2"])
    run(["cmake", "-S", str(BASE), "-B", str(build_y),
         "-DCMAKE_BUILD_TYPE=Release", "-DBENCHMARK_SOLVER=yspecpp",
         f"-DYSPEC_SOURCE_DIR={args.yspec_source.resolve()}",
         f"-DVERIFIED_DEPENDENCY_ROOT={args.dependency_root.resolve()}"])
    run(["cmake", "--build", str(build_y), "--target", "yspecpp_adapter", "--parallel", "2"])
    planetmodel = common.verify_planetmodel(build_y)

    # Keep the exact harness and executables beside each result. Later suite
    # additions must not make an earlier run's provenance ambiguous.
    snapshot_dir = case_root / "source_snapshot"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    source_paths = (BASE / "run_paper_suite.py", BASE / "run_paper2.py",
                    BASE / "dspec_adapter.cpp", BASE / "yspecpp_adapter.cpp",
                    BASE / "CMakeLists.txt", CASE_FILE)
    source_hashes = {}
    for source_path in source_paths:
        snapshot_path = snapshot_dir / source_path.name
        shutil.copyfile(source_path, snapshot_path)
        source_hashes[source_path.name] = sha256(snapshot_path)

    suite_mesh = settings["dspec"]["mesh"]
    if suite_mesh == "single":
        maxstep = settings["dspec"]["maxstep_earth_radii"]
    else:
        maxstep = 0.05
    ns = settings["dspec"]["nq"]
    dspec_log, yspec_log = result_dir / "dspec.stdout.log", result_dir / "yspecpp.stdout.log"
    dspec_wall, dspec_rss_kib = run_adapter(
        [str(build_d / "dspec_adapter"), str(params_path), str(dspec_prefix),
         str(ns), "10", f"{maxstep:.17g}", "paper_suite", suite_mesh],
        dspec_log, result_dir / "dspec.resource.txt")
    yprefix = result_dir / "yspecpp"
    shutil.copyfile(build_d / "dspec_adapter", snapshot_dir / "dspec_adapter")
    shutil.copyfile(build_y / "yspecpp_adapter", snapshot_dir / "yspecpp_adapter")
    binary_hashes = {name: sha256(snapshot_dir / name)
                     for name in ("dspec_adapter", "yspecpp_adapter")}
    yspec_wall, yspec_rss_kib = run_adapter(
        [str(build_y / "yspecpp_adapter"), str(yinput), str(yprefix)],
        yspec_log, result_dir / "yspecpp.resource.txt")

    spectra_d, points_d = common.read_spectra(dspec_prefix.with_suffix(".spectra.csv"))
    spectra_y, points_y = common.read_spectra(yprefix.with_suffix(".spectra.csv"))
    meta_d = common.read_meta(dspec_prefix.with_suffix(".meta.json"))
    meta_y = common.read_meta(yprefix.with_suffix(".meta.json"))
    grid_checks = common.assert_same_grid(meta_d, meta_y, points_d, points_y)
    output_samples = int(meta_d["output_samples"])
    time_d, native_d = common.read_timeseries(dspec_prefix.with_suffix(".timeseries.csv"))
    time_y, native_y = common.read_timeseries(yprefix.with_suffix(".timeseries.csv"))
    time_d, time_y = np.asarray(time_d), np.asarray(time_y)
    native_d, native_y = np.asarray(native_d), np.asarray(native_y)
    shared_time_d, shared_d = common.shared_inverse_transform(
        spectra_d, meta_d, settings["filter_mHz"], output_samples)
    shared_time_y, shared_y = common.shared_inverse_transform(
        spectra_y, meta_y, settings["filter_mHz"], output_samples)
    direct_d = common.direct_inverse_checks(spectra_d, meta_d,
                                             settings["filter_mHz"], shared_d)
    direct_y = common.direct_inverse_checks(spectra_y, meta_y,
                                             settings["filter_mHz"], shared_y)
    unit = UNIT_BY_QUANTITY[settings["quantity"]]
    direct_d = normalize_inverse_check_units(direct_d, unit)
    direct_y = normalize_inverse_check_units(direct_y, unit)
    for label, native, shared in (("DSpecM1D", native_d, shared_d),
                                  ("YSpecpp", native_y, shared_y)):
        if native.shape != shared.shape or not np.isfinite(native).all():
            raise SystemExit(f"{label}: invalid native time series")
        peak = float(np.max(np.abs(shared)))
        error = float(np.max(np.abs(native - shared)))
        if error / peak > 2.0e-10 if peak else error > 1e-20:
            raise SystemExit(f"{label}: native/shared transform mismatch")

    channel_names = ("Z", "N", "E")
    metric_windows = {"full_output": None}
    if settings.get("display_time_window_s"):
        metric_windows["published_display_window"] = tuple(settings["display_time_window_s"])
    time_metrics = {}
    for window_name, window in metric_windows.items():
        mask = np.ones(shared_time_d.shape, dtype=bool)
        if window:
            mask = (shared_time_d >= window[0]) & (shared_time_d <= window[1])
        time_metrics[window_name] = {
            label: compare_series(shared_y[i, mask], shared_d[i, mask], unit)
            for i, label in enumerate(channel_names)
        }
    spectrum_band = tuple(settings.get("spectrum_plot_band_mHz",
                                      [settings["frequency_min_mHz"],
                                       settings["frequency_max_mHz"]]))
    raw_spectrum_metrics = {
        "Z": compare_spectrum_amplitude(spectra_y, spectra_d, points_y,
                                         spectrum_band, 0)
    }
    processed_frequency_d, processed_d = legacy_processed_amplitude_spectrum(
        spectra_d, meta_d, settings["filter_mHz"], 0.5)
    processed_frequency_y, processed_y = legacy_processed_amplitude_spectrum(
        spectra_y, meta_y, settings["filter_mHz"], 0.5)
    if not np.array_equal(processed_frequency_d, processed_frequency_y):
        raise SystemExit("processed-spectrum frequency grids differ")
    processed_mask = ((processed_frequency_d >= spectrum_band[0] / 1000) &
                      (processed_frequency_d <= spectrum_band[1] / 1000))
    processed_delta = np.abs(processed_d[processed_mask] - processed_y[processed_mask])
    processed_peak = float(np.max(processed_y[processed_mask]))
    processed_spectrum_metrics = {
        "quantity": "legacy processed Z spectrum amplitude",
        "unit": "adopted common pipeline arbitrary scale; SI proportional",
        "filter_mHz": settings["filter_mHz"],
        "pre_window_fraction": 0.05,
        "final_window_fraction": 0.5,
        "mean_absolute_amplitude_difference": float(np.mean(processed_delta)),
        "mean_amplitude_difference_percent_of_YSpecpp_peak":
            100 * float(np.mean(processed_delta)) / processed_peak if processed_peak else None,
        "maximum_absolute_amplitude_difference": float(np.max(processed_delta)),
        "maximum_amplitude_difference_percent_of_YSpecpp_peak":
            100 * float(np.max(processed_delta)) / processed_peak if processed_peak else None,
        "frequency_bins": int(np.count_nonzero(processed_mask)),
    }
    import matplotlib.pyplot as plt
    processed_plot, processed_axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    display_frequency = processed_frequency_d[processed_mask] * 1000
    display_d, display_y = processed_d[processed_mask], processed_y[processed_mask]
    processed_axes[0].plot(display_frequency, display_d, label="Current DSpecM1D")
    processed_axes[0].plot(display_frequency, display_y, "--", label="Pinned YSpecpp")
    processed_axes[0].set_ylabel("Processed |Z| (common scale)")
    processed_axes[0].legend()
    processed_axes[1].plot(display_frequency,
                           100 * (display_d - display_y) / processed_peak if processed_peak else display_d-display_y,
                           color="black")
    processed_axes[1].set_ylabel("Difference (% of YSpecpp peak)")
    processed_axes[1].set_xlabel("Frequency (mHz)")
    for axis in processed_axes:
        axis.grid(alpha=0.25)
    processed_plot.tight_layout()
    processed_plot.savefig(case_root / "processed_spectrum.png", dpi=170)
    plt.close(processed_plot)
    plot_spectrum(case_root / "raw_spectrum.png", points_d, spectra_d,
                  spectra_y, spectrum_band, f"{unit}/Hz")
    plot_timeseries(case_root / "comparison_full.png", shared_time_d,
                    shared_d, shared_y, settings["quantity"], None,
                    time_metrics["full_output"])
    if "published_display_window" in time_metrics:
        plot_timeseries(case_root / "comparison_display.png", shared_time_d,
                        shared_d, shared_y, settings["quantity"],
                        metric_windows["published_display_window"],
                        time_metrics["published_display_window"])

    source_head = subprocess.check_output(
        ["git", "-C", str(args.dspec_source), "rev-parse", "HEAD"], text=True).strip()
    maxsteps = [float(line.split(":", 1)[1]) for line in
                dspec_log.read_text().splitlines() if line.startswith("Max step:")]
    provenance = {
        "case": args.case,
        "paper_figure": settings["paper_figure"],
        "run_scope": settings["run_scope"],
        "criteria_file": str(args.criteria.resolve()),
        "criteria_sha256": sha256(args.criteria),
        "case_table": str(CASE_FILE),
        "case_table_sha256": sha256(CASE_FILE),
        "source_parameter_file": str(input_path),
        "source_parameter_sha256": sha256(input_path),
        "model_file": str(model_source),
        "model_sha256": model_hash,
        "DSpecM1D_HEAD_at_run": source_head,
        "harness_source_sha256": {
            **source_hashes,
        },
        "harness_source_snapshot": "source_snapshot/",
        "adapter_binary_sha256": binary_hashes,
        "yspecpp_commit": common.YSPEC_COMMIT,
        "yspecpp_tree_sha256": yspec_tree,
        "planetmodel": planetmodel,
        "dependency_shas": common.PINNED_DEPS,
        "settings": settings,
        "quantity_and_units": {"quantity": settings["quantity"], "unit": unit,
                                "spectrum_unit": f"{unit}/Hz"},
        "source_time": "Heaviside step moment; no extra source filter",
        "coordinate_convention": "case latitudes geocentric; converted to YSpecpp geographic inputs at adapter boundary",
        "model_interpolation_note": "DSpecM1D interpolates velocities before moduli; YSpecpp splines moduli.",
        "dspec_dynamic_mesh_steps": maxsteps,
        "dspec_mesh_strategy": suite_mesh,
        "dspec_mesh_cap_earth_radii": settings["dspec"].get("maxstep_earth_radii"),
        "runtime_environment": {"OMP_NUM_THREADS": "2", "OMP_DYNAMIC": "FALSE"},
        "runtime_s": {"dspec": dspec_wall, "yspecpp": yspec_wall},
        "peak_child_rss_kib": {"dspec": dspec_rss_kib, "yspecpp": yspec_rss_kib},
        "interpolation_caveat": "Continuous model representations differ; comparison is not an exact-model identity test.",
    }
    (case_root / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    df = float(meta_d["frequency_spacing_hz"])
    full_frequency_bins = int(round(
        (case["frequency_max_mHz"] - case["frequency_min_mHz"]) / 1000 / df))
    pilot_frequency_bins = max(1, len(spectra_d) // 3)
    pilot_scale = full_frequency_bins / pilot_frequency_bins if args.pilot else 1.0
    summary = {
        "case": args.case,
        "paper_figure": settings["paper_figure"],
        "run_scope": settings["run_scope"],
        "quantity": settings["quantity"],
        "unit": unit,
        "spectrum_unit": f"{unit}/Hz",
        "grid_checks": grid_checks,
        "native_inverse_checks": {"unit": unit,
                                  "dspec": direct_d, "yspecpp": direct_y},
        "time_metrics": time_metrics,
        "raw_spectrum_metrics": raw_spectrum_metrics,
        "legacy_processed_spectrum_metrics": processed_spectrum_metrics,
        "runtime": {"seconds": provenance["runtime_s"],
                    "peak_child_rss_kib": provenance["peak_child_rss_kib"],
                    "pilot_linear_frequency_extrapolation_seconds":
                        {k: v * pilot_scale for k, v in provenance["runtime_s"].items()}},
        "dspec_mesh": {"strategy": suite_mesh, "maxstep_per_chunk": maxsteps,
                       "cap_earth_radii": provenance["dspec_mesh_cap_earth_radii"]},
        "paper_context": criteria["case_inventory"][
            [c["case"] for c in criteria["case_inventory"]].index(expected_case)
        ]["paper_result"],
        "artifact_paths": {
            "raw_spectra": ["results/dspec.spectra.csv", "results/yspecpp.spectra.csv"],
            "filtered_traces": ["results/dspec.timeseries.csv", "results/yspecpp.timeseries.csv"],
            "raw_spectrum_plot": "raw_spectrum.png",
            "legacy_processed_spectrum_plot": "processed_spectrum.png",
            "full_time_plot": "comparison_full.png",
            "display_time_plot": "comparison_display.png" if "published_display_window" in time_metrics else None,
            "adapter_stdout": ["results/dspec.stdout.log", "results/yspecpp.stdout.log"],
            "provenance": "provenance.json",
        },
    }
    (case_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"case_output": str(case_root), "run_scope": settings["run_scope"],
                      "runtime": summary["runtime"],
                      "max_channel_mean_percent": max(
                          row["mean_difference_percent_of_reference_peak"]
                          for case_metrics in time_metrics.values()
                          for row in case_metrics.values()
                          if row["mean_difference_percent_of_reference_peak"] is not None)}, indent=2))


if __name__ == "__main__":
    main()
