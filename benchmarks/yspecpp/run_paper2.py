#!/usr/bin/env python3
"""Run the paper Figure 2 DSpecM1D/YSpecpp displacement comparison."""

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

import numpy as np

from run_benchmark import (
    PINNED_DEPS,
    YSPEC_ARCHIVE_SHA256,
    YSPEC_COMMIT,
    YSPEC_FILETREE_SHA256,
    filetree_sha256,
    geographic_latitude,
    read_meta,
    read_spectra,
    read_timeseries,
    verify_dependencies,
)


BASE = Path(__file__).resolve().parent
DSPEC_DEFAULT = BASE.parents[1]
YSPEC_DEFAULT = Path("/tmp/yspecpp-scope/da380-YSpecpp-84be5a6")
DEPENDENCY_DEFAULT = DSPEC_DEFAULT / "build/phase12-fix/fetched-pinned/_deps"
PLANETMODEL_COMMIT = "fa6e230fde25152d40f27d9839de44c4f3182a24"
PUBLISHED = {
    "source_depth_km": 647.1,
    "source_geocentric_latitude_deg": -13.82,
    "source_longitude_deg": -67.25,
    "receiver_geocentric_latitude_deg": 80.0,
    "receiver_longitude_deg": 0.0,
    "receiver_depth_km": 3.0,
    "duration_min": 300.0,
    "time_step_s": 1.0,
    "degree_min": 0,
    "degree_max": 300,
    "frequency_min_mHz": 1.5,
    "frequency_max_mHz": 21.0,
    "filter_mHz": [1.9, 2.0, 20.0, 20.1],
    "attenuation": False,
    "gravity": "full",
    "quantity": "displacement",
    "dspec_sem_points": 5,
    "dspec_relative_error": 1.0e-5,
    "source_time_convention": "step moment; no extra source filter",
}
MOMENT_NM = {
    "rr": -7.59e20,
    "rt": -2.50e21,
    "rp": 4.20e19,
    "tt": 7.75e20,
    "tp": -2.48e20,
    "pp": -1.60e19,
}


def run(command: list[str], *, cwd: Path | None = None,
        env: dict[str, str] | None = None,
        log_path: Path | None = None) -> None:
    print("+", " ".join(map(str, command)), flush=True)
    if log_path is None:
        subprocess.run(command, cwd=cwd, env=env, check=True)
        return
    completed = subprocess.run(command, cwd=cwd, env=env, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    log_path.write_text(completed.stdout)
    print(completed.stdout, end="", flush=True)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, command)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_planetmodel(build_dir: Path) -> dict[str, str]:
    checkout = build_dir / "_deps/planetmodel-src"
    if not checkout.is_dir() or not (checkout / ".git").exists():
        raise SystemExit("YSpecpp did not create its pinned planetmodel checkout")
    actual = subprocess.check_output(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    if actual != PLANETMODEL_COMMIT:
        raise SystemExit(f"planetmodel is {actual}, expected {PLANETMODEL_COMMIT}")
    dirty = subprocess.check_output(
        ["git", "-C", str(checkout), "status", "--porcelain",
         "--untracked-files=normal"], text=True)
    if dirty:
        raise SystemExit("planetmodel has uncommitted or untracked files")
    return {"commit": actual, "clean": "true"}


def dspec_parameters(model_name: str, prefix: Path,
                     settings: dict) -> str:
    m = MOMENT_NM
    values = [
        str(prefix), f"models/{model_name}", "4", "0", "2", "0", "0",
        f"{settings['dspec_relative_error']:.17g}",
        str(settings["degree_min"]), str(settings["degree_max"]),
        f"{settings['frequency_min_mHz']:.17g}",
        f"{settings['frequency_max_mHz']:.17g}",
        f"{settings['duration_min']:.17g}",
        f"{settings['time_step_s']:.17g}",
        *(f"{corner:.17g}" for corner in settings["filter_mHz"]),
        f"{settings['source_depth_km']:.17g}",
        f"{settings['source_geocentric_latitude_deg']:.17g}",
        f"{settings['source_longitude_deg']:.17g}",
        *(f"{MOMENT_NM[key]:.17g}" for key in ("rr", "rt", "rp", "tt", "tp", "pp")),
        f"{settings['receiver_depth_km']:.17g}", "1",
        f"{settings['receiver_geocentric_latitude_deg']:.17g} "
        f"{settings['receiver_longitude_deg']:.17g}",
    ]
    return "\n".join(values) + "\n"


def yspec_input(model: Path, settings: dict) -> dict:
    return {
        "model": "../../inputs/data/models/" + model.name,
        "attenuation": False,
        "gravity": "full",
        "output": {
            "prefix": "unused",
            "quantity": "displacement",
            "potential_and_tilt_corrections": False,
        },
        "degrees": {"min": settings["degree_min"],
                    "max": settings["degree_max"]},
        "frequencies_mHz": {
            "min": settings["frequency_min_mHz"],
            "max": settings["frequency_max_mHz"],
        },
        "time_series": {"length_min": settings["duration_min"],
                        "time_step_s": settings["time_step_s"]},
        "filter_mHz": settings["filter_mHz"],
        "source": {
            "depth_km": settings["source_depth_km"],
            "latitude": geographic_latitude(
                settings["source_geocentric_latitude_deg"]),
            "longitude": settings["source_longitude_deg"],
            "moment_tensor_Nm": MOMENT_NM,
        },
        "receivers": {
            "depth_km": settings["receiver_depth_km"],
            "locations": [[geographic_latitude(
                settings["receiver_geocentric_latitude_deg"]),
                settings["receiver_longitude_deg"]]],
        },
        "solver": {"tolerance": 1.0e-5, "start_decay": 11.0},
    }


def half_cosine(frequency_hz: np.ndarray,
                corners_mhz: list[float]) -> np.ndarray:
    f1, f2, f3, f4 = [value / 1000.0 for value in corners_mhz]
    result = np.zeros_like(frequency_hz)
    rising = (frequency_hz >= f1) & (frequency_hz < f2)
    flat = (frequency_hz >= f2) & (frequency_hz < f3)
    falling = (frequency_hz >= f3) & (frequency_hz < f4)
    result[rising] = 0.5 * (1.0 - np.cos(
        np.pi * (frequency_hz[rising] - f1) / (f2 - f1)))
    result[flat] = 1.0
    result[falling] = 0.5 * (1.0 - np.cos(
        np.pi * (f4 - frequency_hz[falling]) / (f4 - f3)))
    return result


def shared_inverse_transform(spectra: dict, meta: dict,
                             filter_corners: list[float],
                             output_samples: int) -> tuple[np.ndarray, np.ndarray]:
    nfft = int(meta["samples"])
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    frequency = np.arange(nfft // 2 + 1, dtype=float) * df
    filt = half_cosine(frequency, filter_corners)
    bins = np.zeros((3, nfft // 2 + 1), dtype=np.complex128)
    for (index, component), value in spectra.items():
        bins[component, index] = value
    time = np.arange(output_samples, dtype=float) * dt
    # FFTW's c2r transform is unnormalized; NumPy irfft divides by nfft.
    traces = np.fft.irfft(bins * filt[None, :], n=nfft, axis=1)
    traces = traces[:, :output_samples] * (
        nfft * df * np.exp(epsilon * time)[None, :])
    if not np.isfinite(traces).all():
        raise SystemExit("shared inverse transform produced non-finite values")
    return time, traces


def compare_traces(reference: np.ndarray, candidate: np.ndarray) -> dict:
    difference = candidate - reference
    peak = float(np.max(np.abs(reference)))
    max_abs = float(np.max(np.abs(difference)))
    return {
        "reference_peak_m": peak,
        "mean_absolute_difference_m": float(np.mean(np.abs(difference))),
        "mean_absolute_difference_percent_of_reference_peak":
            100.0 * float(np.mean(np.abs(difference))) / peak if peak else None,
        "max_absolute_difference_m": max_abs,
        "peak_difference_percent_of_reference_peak":
            100.0 * max_abs / peak if peak else None,
        "relative_l2": float(np.linalg.norm(difference) /
                              np.linalg.norm(reference))
            if np.linalg.norm(reference) else None,
    }


def assert_same_grid(meta_d: dict, meta_y: dict,
                     points_d: dict, points_y: dict) -> dict:
    checks = {}
    for key in ("samples", "output_samples"):
        if meta_d[key] != meta_y[key]:
            raise SystemExit(f"grid mismatch for {key}: "
                             f"{meta_d[key]} vs {meta_y[key]}")
        checks[key] = meta_d[key]
    for key in ("sample_interval_s", "frequency_spacing_hz",
                "damping_epsilon_rad_s"):
        left, right = float(meta_d[key]), float(meta_y[key])
        if not (math.isfinite(left) and math.isfinite(right)):
            raise SystemExit(f"non-finite grid metadata: {key}")
        delta = abs(left - right)
        if delta > 1e-12 * max(abs(left), abs(right), 1e-30):
            raise SystemExit(f"grid mismatch for {key}: {left} vs {right}")
        checks[key] = {"dspec": left, "yspecpp": right,
                       "absolute_difference": delta}
    if set(points_d) != set(points_y):
        raise SystemExit("frequency/component samples do not match")
    for key in sorted(points_d):
        for field in ("frequency_hz", "omega_real_rad_s",
                      "omega_imag_rad_s"):
            left, right = points_d[key][field], points_y[key][field]
            tolerance = 1e-10 if field == "omega_real_rad_s" else 1e-12
            if abs(left - right) > tolerance * max(abs(left), abs(right), 1e-30):
                raise SystemExit(f"bin {key} mismatch for {field}: "
                                 f"{left} vs {right}")
    checks["frequency_component_samples"] = len(points_d)
    return checks


def direct_inverse_checks(spectra: dict, meta: dict,
                          corners: list[float],
                          traces: np.ndarray) -> list[dict]:
    """Check selected samples using a direct positive-frequency DFT sum."""
    nfft = int(meta["samples"])
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    indices = sorted({0, 1, min(137, traces.shape[1] - 1),
                      traces.shape[1] // 2, traces.shape[1] - 1})
    checks = []
    for component in range(3):
        bins = sorted(index for index, channel in spectra
                      if channel == component)
        for sample in indices:
            total = 0.0
            for frequency_bin in bins:
                frequency = frequency_bin * df
                weight = half_cosine(np.asarray([frequency]), corners)[0]
                value = spectra[(frequency_bin, component)] * weight
                phase = 2.0 * math.pi * frequency_bin * sample / nfft
                if frequency_bin == 0 or (
                        nfft % 2 == 0 and frequency_bin == nfft // 2):
                    total += value.real * (-1.0 if frequency_bin == nfft // 2
                                            and sample % 2 else 1.0)
                else:
                    total += 2.0 * (value * complex(math.cos(phase),
                                                     math.sin(phase))).real
            expected = df * math.exp(epsilon * sample * dt) * total
            actual = float(traces[component, sample])
            peak = float(np.max(np.abs(traces[component])))
            error = abs(expected - actual)
            scaled = error / peak if peak else error
            if not math.isfinite(expected) or not math.isfinite(actual) or scaled > 2e-10:
                raise SystemExit("direct inverse-transform oracle failed at "
                                 f"component {component}, sample {sample}: "
                                 f"peak-scaled error={scaled:.3g}")
            checks.append({"component": component, "sample": sample,
                           "expected_m": expected, "actual_m": actual,
                           "absolute_error_m": error,
                           "error_over_channel_peak": scaled})
    return checks


def write_comparison_csv(path: Path, time: np.ndarray,
                         dspec: np.ndarray, yspec: np.ndarray) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "Z_dspec_m", "Z_yspecpp_m",
                         "N_dspec_m", "N_yspecpp_m", "E_dspec_m",
                         "E_yspecpp_m"])
        for index, t in enumerate(time):
            writer.writerow([f"{t:.17g}",
                             *(f"{value:.17g}" for value in
                               (dspec[0, index], yspec[0, index],
                                dspec[1, index], yspec[1, index],
                                dspec[2, index], yspec[2, index]))])


def make_plot(path: Path, time: np.ndarray,
              dspec: np.ndarray, yspec: np.ndarray,
              metrics: dict[str, dict], pilot: bool) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    for component, label in enumerate(("Z", "N", "E")):
        peak = metrics[label]["reference_peak_m"]
        peak = peak if peak else 1.0
        axes[component, 0].plot(time / 3600, dspec[component] / peak,
                               label="Current DSpecM1D", linewidth=1.0)
        axes[component, 0].plot(time / 3600, yspec[component] / peak,
                               label="Pinned YSpecpp", linewidth=1.0,
                               linestyle="--")
        axes[component, 0].set_ylabel(f"{label} / peak")
        axes[component, 0].text(
            0.02, 0.96, f"YSpecpp peak = {peak:.3g} m",
            transform=axes[component, 0].transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"})
        axes[component, 0].grid(alpha=0.25)
        axes[component, 1].plot(
            time / 3600,
            100 * (dspec[component] - yspec[component]) / peak,
            color="black", linewidth=0.8)
        axes[component, 1].set_ylabel("Difference (% of YSpecpp peak)")
        axes[component, 1].grid(alpha=0.25)
    axes[0, 0].set_title("Filtered displacement")
    axes[0, 1].set_title("DSpecM1D − YSpecpp")
    fig.supylabel("Displacement / YSpecpp peak (dimensionless)")
    axes[2, 0].set_xlabel("Time (hours)")
    axes[2, 1].set_xlabel("Time (hours)")
    axes[0, 0].legend(loc="upper right")
    title = ("Diagnostic narrow-band pilot, 11.15–11.35 mHz "
             "(not the full Figure 2 result)" if pilot else
             "Figure 2 full-band comparison: Bolivia source, PREM-TI ocean, "
             "degrees 0–300")
    fig.suptitle(title)
    fig.tight_layout(rect=(0.06, 0.0, 1.0, 0.97))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dspec-source", type=Path, default=DSPEC_DEFAULT)
    parser.add_argument("--yspec-source", type=Path, default=YSPEC_DEFAULT)
    parser.add_argument("--dependency-root", type=Path,
                        default=DEPENDENCY_DEFAULT)
    parser.add_argument("--protocol", type=Path,
                        default=BASE / "cases/paper2-protocol.json",
                        help="protocol JSON (defaults to the tracked Figure 2 protocol)")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--full", action="store_true",
                        help="run the published 1.5–21 mHz band")
    args = parser.parse_args()
    default_root = args.dspec_source / "build/phase3/paper2-20261020"
    outdir = (args.output or default_root / ("full" if args.full else "pilot")).resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"output directory must be empty: {outdir}")
    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"
    if not args.dspec_source.is_dir():
        parser.error(f"missing DSpecM1D checkout: {args.dspec_source}")
    if not args.yspec_source.is_dir():
        parser.error(f"missing pinned YSpecpp snapshot: {args.yspec_source}")
    verify_dependencies(args.dependency_root)
    yspec_tree_hash, yspec_file_count = filetree_sha256(args.yspec_source)
    if yspec_tree_hash != YSPEC_FILETREE_SHA256 or yspec_file_count != 76:
        raise SystemExit("YSpecpp source tree differs from its pinned snapshot")

    settings = dict(PUBLISHED)
    if not args.full:
        # Four to six frequencies straddle the multi-SEM split near 11.25 mHz.
        # The degree cap, 5 h FFT, nq=5, and 1e-5 accuracy remain at paper scale.
        settings["frequency_min_mHz"] = 11.15
        settings["frequency_max_mHz"] = 11.35
    protocol_path = args.protocol.resolve()
    protocol = json.loads(protocol_path.read_text())
    source_head = subprocess.check_output(
        ["git", "-C", str(args.dspec_source), "rev-parse", "HEAD"],
        text=True).strip()
    if protocol.get("status") not in (None, "TRACKED_FIGURE2_PROTOCOL",
                                       "FROZEN_BEFORE_PAPER2_RUN"):
        raise SystemExit("unsupported paper2 protocol status")
    criteria = protocol.get("criteria", {})
    if criteria.get("paper_mean_misfit_percent_upper_exclusive") != 1.0:
        raise SystemExit("paper2 mean-misfit gate differs from the frozen 1% limit")
    if criteria.get("normalization") != (
            "max absolute YSpecpp trace per component on inclusive 0..18000 s"):
        raise SystemExit("paper2 comparison normalization differs from protocol")
    if criteria.get("published_pointwise_max_percent_context_only") != 0.4:
        raise SystemExit("paper2 pointwise context differs from protocol")
    expected_source_head = protocol.get("source_commit")
    if expected_source_head and source_head != expected_source_head:
        raise SystemExit(f"source HEAD {source_head} differs from frozen protocol")
    if sha256(args.dspec_source / "data/params/ex2.txt") != protocol[
            "preserved_inputs_sha256"]["data/params/ex2.txt"]:
        raise SystemExit("ex2 source input differs from frozen protocol")
    if PUBLISHED["frequency_min_mHz"] != protocol["solve_band_mHz"][0] or \
            PUBLISHED["frequency_max_mHz"] != protocol["solve_band_mHz"][1]:
        raise SystemExit("runner frequency band differs from frozen protocol")
    frozen_checks = {
        "duration_s": int(PUBLISHED["duration_min"] * 60),
        "dt_s": PUBLISHED["time_step_s"],
        "degree_range": [PUBLISHED["degree_min"], PUBLISHED["degree_max"]],
        "filter_corners_mHz": PUBLISHED["filter_mHz"],
        "nq": PUBLISHED["dspec_sem_points"],
        "relative_error": PUBLISHED["dspec_relative_error"],
        "gravity": PUBLISHED["gravity"],
        "attenuation": PUBLISHED["attenuation"],
        "receiver_depth_km": PUBLISHED["receiver_depth_km"],
    }
    for key, expected in frozen_checks.items():
        if protocol.get(key) != expected:
            raise SystemExit(f"runner setting {key} differs from frozen protocol")
    inputs = outdir / "inputs"
    models = inputs / "data/models"
    params_dir = inputs / "data/params"
    yjson_dir = inputs / "yspec"
    result_dir = outdir / "results"
    for path in (models, params_dir, yjson_dir, result_dir):
        path.mkdir(parents=True, exist_ok=True)
    model_source = args.dspec_source / "data/models/prem.200.noatten.txt"
    if sha256(model_source) != protocol["preserved_inputs_sha256"][
            "data/models/prem.200.noatten.txt"]:
        raise SystemExit("PREM-TI model differs from frozen protocol")
    model_target = models / model_source.name
    shutil.copyfile(model_source, model_target)

    build_d = outdir / "build-dspec"
    build_y = outdir / "build-yspecpp"
    dspec_build_args = [
        "cmake", "-S", str(BASE), "-B", str(build_d),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBENCHMARK_SOLVER=dspec",
        f"-DDSPEC_SOURCE_DIR={args.dspec_source.resolve()}",
        f"-DVERIFIED_DEPENDENCY_ROOT={args.dependency_root.resolve()}",
    ]
    yspec_build_args = [
        "cmake", "-S", str(BASE), "-B", str(build_y),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBENCHMARK_SOLVER=yspecpp",
        f"-DYSPEC_SOURCE_DIR={args.yspec_source.resolve()}",
        f"-DVERIFIED_DEPENDENCY_ROOT={args.dependency_root.resolve()}",
    ]
    run(dspec_build_args)
    run(["cmake", "--build", str(build_d), "--target", "dspec_adapter",
         "--parallel", "2"])
    run(yspec_build_args)
    run(["cmake", "--build", str(build_y), "--target", "yspecpp_adapter",
         "--parallel", "2"])
    planetmodel = verify_planetmodel(build_y)

    params_path = params_dir / "paper2.txt"
    dspec_prefix = result_dir / "dspec"
    params_path.write_text(dspec_parameters(model_target.name, dspec_prefix,
                                            settings))
    yinput_path = yjson_dir / "paper2.json"
    yinput_path.write_text(json.dumps(yspec_input(model_target, settings),
                                      indent=2) + "\n")
    provenance = {
        "run_scope": "full" if args.full else "bounded pilot",
        "protocol": "published-style Figure 2 comparison; new output, no fit",
        "protocol_file": str(protocol_path),
        "protocol_sha256": sha256(protocol_path),
        "settings": settings,
        "moment_tensor_Nm_from": "data/params/ex2.txt",
        "moment_tensor_Nm": MOMENT_NM,
        "moment_tensor_input_sha256": sha256(
            args.dspec_source / "data/params/ex2.txt"),
        "model_source": str(model_source),
        "model_sha256": sha256(model_source),
        "model_note": ("DSpecM1D interpolates velocities before constructing "
                       "moduli; YSpecpp splines elastic moduli. The two "
                       "continuous interpolants are not identical."),
        "dspec_source_head": subprocess.check_output(
            ["git", "-C", str(args.dspec_source), "rev-parse", "HEAD"],
            text=True).strip(),
        "dspec_tracked_source_diff_sha256": hashlib.sha256(
            subprocess.check_output(["git", "-C", str(args.dspec_source),
                                     "diff", "--binary", "HEAD", "--",
                                     "DSpecM1D"])).hexdigest(),
        "yspecpp_commit": YSPEC_COMMIT,
        "yspecpp_tree_sha256": yspec_tree_hash,
        "yspecpp_file_count": yspec_file_count,
        "yspecpp_archive_sha256": YSPEC_ARCHIVE_SHA256,
        "planetmodel": planetmodel,
        "dependency_shas": PINNED_DEPS,
        "harness_sources": {
            path.name: sha256(path)
            for path in (BASE / "run_paper2.py", BASE / "dspec_adapter.cpp",
                         BASE / "yspecpp_adapter.cpp", BASE / "CMakeLists.txt",
                         BASE / "README.md", BASE / "cases/paper2-protocol.json")
        },
        "runtime_environment": {key: os.environ[key] for key in (
            "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_MAX_ACTIVE_LEVELS")},
        "latitude_convention": {
            "dspec_input": "geocentric degrees",
            "yspecpp_input": "geographic degrees converted with the formula "
                             "documented in run_benchmark.py",
            "paper_ambiguity": "paper gives site coordinates only; this uses "
                               "DSpecM1D's declared geocentric coordinates",
        },
        "frequency_grid": {
            "expected_samples": 32768,
            "expected_frequency_spacing_hz": 1.0 / 32768.0,
            "damping_e_folds": 5,
        },
        "filter": {
            "corners_mHz": settings["filter_mHz"],
            "application": "one direct half-cosine taper to raw displacement "
                           "spectra before inverse transform",
            "spectra_saved": "raw, unfiltered physical displacement per Hz",
            "native_timeseries_saved": "filtered physical displacement in m",
            "shared_transform": "NumPy irfft with physical df and damping undo",
        },
        "git_status_before_run": subprocess.check_output(
            ["git", "-C", str(args.dspec_source), "status", "--short",
             "--branch"], text=True).splitlines(),
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    env = dict(os.environ)
    dspec_meta_prefix = result_dir / "dspec"
    yspec_meta_prefix = result_dir / "yspecpp"
    dspec_log = result_dir / "dspec.stdout.log"
    yspec_log = result_dir / "yspecpp.stdout.log"
    run([str(build_d / "dspec_adapter"), str(params_path),
         str(dspec_meta_prefix), str(settings["dspec_sem_points"]), "10",
         "0.05", "paper2"], env=env, log_path=dspec_log)
    run([str(build_y / "yspecpp_adapter"), str(yinput_path),
         str(yspec_meta_prefix)], env=env, log_path=yspec_log)

    spectra_d, points_d = read_spectra(dspec_meta_prefix.with_suffix(".spectra.csv"))
    spectra_y, points_y = read_spectra(yspec_meta_prefix.with_suffix(".spectra.csv"))
    meta_d = read_meta(dspec_meta_prefix.with_suffix(".meta.json"))
    meta_y = read_meta(yspec_meta_prefix.with_suffix(".meta.json"))
    grid_checks = assert_same_grid(meta_d, meta_y, points_d, points_y)
    mesh_steps = [float(line.split(":", 1)[1])
                  for line in dspec_log.read_text().splitlines()
                  if line.startswith("Max step:")]
    if not mesh_steps or not all(math.isfinite(value) and value > 0
                                 for value in mesh_steps):
        raise SystemExit("DSpecM1D dynamic mesh steps missing from run log")
    output_samples = int(meta_d["output_samples"])
    time_d, native_d = read_timeseries(
        dspec_meta_prefix.with_suffix(".timeseries.csv"))
    time_y, native_y = read_timeseries(
        yspec_meta_prefix.with_suffix(".timeseries.csv"))
    time_d = np.asarray(time_d, dtype=float)
    time_y = np.asarray(time_y, dtype=float)
    native_d = np.asarray(native_d, dtype=float)
    native_y = np.asarray(native_y, dtype=float)
    if len(time_d) != output_samples or len(time_y) != output_samples:
        raise SystemExit("native time-series length differs from metadata")
    time_shared_d, common_d = shared_inverse_transform(
        spectra_d, meta_d, settings["filter_mHz"], output_samples)
    time_shared_y, common_y = shared_inverse_transform(
        spectra_y, meta_y, settings["filter_mHz"], output_samples)
    direct_checks = {
        "dspec": direct_inverse_checks(spectra_d, meta_d,
                                        settings["filter_mHz"], common_d),
        "yspecpp": direct_inverse_checks(spectra_y, meta_y,
                                          settings["filter_mHz"], common_y),
    }
    if not np.allclose(time_shared_d, time_shared_y, rtol=0.0, atol=1e-12):
        raise SystemExit("shared physical time axes differ")
    if not np.allclose(time_d, time_y, rtol=0.0, atol=1e-12):
        raise SystemExit("native time axes differ")
    native_checks = {}
    for label, native, common in (("dspec", native_d, common_d),
                                  ("yspecpp", native_y, common_y)):
        if native.shape != common.shape:
            raise SystemExit(f"{label}: native/shared time-series shape differs")
        peak = float(np.max(np.abs(common)))
        error = float(np.max(np.abs(native - common)))
        scaled = error / peak if peak else error
        if not math.isfinite(scaled):
            raise SystemExit(f"{label}: native/shared comparison is non-finite")
        if scaled > 2.0e-10:
            raise SystemExit(f"{label}: native and shared filtered traces differ "
                             f"by {scaled:.3g} of peak")
        native_checks[label] = {
            "max_abs_difference_m": error,
            "max_difference_over_shared_peak": scaled,
        }
    metrics = {
        label: compare_traces(common_y[index], common_d[index])
        for index, label in enumerate(("Z", "N", "E"))
    }
    for label, row in metrics.items():
        if row["reference_peak_m"] <= 0.0 or not all(
                math.isfinite(value) for value in row.values()
                if isinstance(value, float)):
            raise SystemExit(f"{label}: invalid or non-finite comparison metrics")
    max_mean_percent = max(
        row["mean_absolute_difference_percent_of_reference_peak"]
        for row in metrics.values() if row["mean_absolute_difference_percent_of_reference_peak"] is not None)
    max_pointwise_percent = max(
        row["peak_difference_percent_of_reference_peak"]
        for row in metrics.values() if row["peak_difference_percent_of_reference_peak"] is not None)
    write_comparison_csv(outdir / "filtered_traces.csv", time_shared_d,
                         common_d, common_y)
    make_plot(outdir / "comparison.png", time_shared_d,
              common_d, common_y, metrics, pilot=not args.full)
    full_gate_passed = max_mean_percent < 1.0 if args.full else None
    summary = {
        "provenance": "provenance.json",
        "settings": settings,
        "observable": "three-component displacement (ZNE), metres",
        "grid_checks": grid_checks,
        "native_vs_shared_transform_checks": native_checks,
        "direct_inverse_transform_checks": direct_checks,
        "run_metadata": {"dspec": meta_d, "yspecpp": meta_y},
        "dspec_dynamic_mesh": {
            "strategy": meta_d["mesh_strategy"],
            "derived_maxstep_per_frequency_chunk": mesh_steps,
            "derived_nskip_from_frequency_grid": max(
                1, (int(meta_d["frequency_last_index_exclusive"]) -
                    int(meta_d["frequency_first_index"])) // 20),
        },
        "metrics": metrics,
        "maximum_channel_mean_percent": max_mean_percent,
        "maximum_pointwise_difference_percent_of_each_YSpecpp_peak":
            max_pointwise_percent,
        "full_band_mean_misfit_gate": {
            "applicable": bool(args.full),
            "limit_percent_exclusive": 1.0,
            "maximum_channel_mean_percent": max_mean_percent,
            "passed": full_gate_passed,
        },
        "published_comparison_statement": {
            "reported_mean_misfit_bound_percent": 1.0,
            "reported_three_channel_pointwise_max_percent": 0.4,
            "current_result_is_a_new_pinned_YSpecpp_comparison": True,
            "accepted_full_figure2_result": bool(args.full),
            "pilot_is_diagnostic_only": not args.full,
        },
        "artifacts": {
            "figure": "comparison.png",
            "filtered_traces": "filtered_traces.csv",
            "raw_spectra": ["results/dspec.spectra.csv",
                             "results/yspecpp.spectra.csv"],
            "adapter_stdout": ["results/dspec.stdout.log",
                               "results/yspecpp.stdout.log"],
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (outdir / "summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["component", "reference_peak_m",
                         "mean_absolute_difference_m",
                         "mean_difference_percent_of_reference_peak",
                         "max_absolute_difference_m",
                         "peak_difference_percent_of_reference_peak",
                         "relative_l2"])
        for label, row in metrics.items():
            writer.writerow([label, row["reference_peak_m"],
                             row["mean_absolute_difference_m"],
                             row["mean_absolute_difference_percent_of_reference_peak"],
                             row["max_absolute_difference_m"],
                             row["peak_difference_percent_of_reference_peak"],
                             row["relative_l2"]])
    if args.full and not full_gate_passed:
        raise SystemExit(f"full-band mean misfit gate failed: "
                         f"{max_mean_percent:.6g}% is not below 1%")
    print(json.dumps({"output": str(outdir), "settings": settings,
                      "metrics": metrics,
                      "maximum_channel_mean_percent": max_mean_percent,
                      "maximum_pointwise_percent": max_pointwise_percent},
                     indent=2))


if __name__ == "__main__":
    main()
