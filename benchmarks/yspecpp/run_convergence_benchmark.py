#!/usr/bin/env python3
"""Run the paper Figure 6 internal SEM spatial-convergence benchmark.

The finest mesh is a declared numerical reference, not an exact solution. Run
--pilot first to measure high-frequency coarse/fine cost; --full starts the
sequential 4/5/6-point sweep and should be launched only after that pilot has
been reviewed.
"""
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
from datetime import datetime, timezone

import numpy as np

BASE = Path(__file__).resolve().parent
DSPEC_DEFAULT = BASE.parents[1]
CRITERIA_DEFAULT = BASE / "cases/paper-suite-protocol.json"
MODEL_REL = Path("data/models/prem.200.no.noatten.txt")
PARAM_REL = Path("data/params/ex5.txt")
MOMENT_NM = {
    "rr": -7.59e20, "rt": -2.50e21, "rp": 4.20e19,
    "tt": 7.75e20, "tp": -2.48e20, "pp": -1.60e19,
}
PUBLISHED = {
    "source_depth_km": 647.1,
    "source_geocentric_latitude_deg": -13.82,
    "source_longitude_deg": -67.25,
    "receiver_depth_km": 0.0,
    "receiver_geocentric_latitude_deg": 51.68,
    "receiver_longitude_deg": 103.64,
    "degree_min": 0, "degree_max": 750,
    "frequency_min_mHz": 0.1, "frequency_max_mHz": 50.0,
    "duration_min": 500.0, "time_step_s": 1.0,
    "filter_mHz": [0.1, 0.15, 45.0, 50.0],
    "quantity": "displacement", "gravity": "full", "attenuation": False,
    "relative_error": 1e-5, "nskip": 1,
}
EPW_VALUES = (0.5, 1.0, 2.0, 3.0, 7.0, 14.0)
REFERENCE = (6, 14.0)
PAPER_TARGETS = {(4, 3.0), (5, 2.0), (6, 1.0)}
CHANNELS = ("Z", "N", "E")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def filetree_sha256(root: Path) -> tuple[str, int]:
    entries = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        entries.append((path.relative_to(root).as_posix(), sha256(path)))
    return hashlib.sha256(json.dumps(entries, separators=(",", ":")).encode()).hexdigest(), len(entries)


def parse_model_min_s_velocity(path: Path) -> tuple[float, float]:
    """Return minimum positive Vsv/Vsh (m/s) and outer radius (m)."""
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 8:
            continue
        try:
            rows.append([float(value) for value in fields[:8]])
        except ValueError:
            continue
    if not rows:
        raise SystemExit(f"could not parse radial model rows from {path}")
    # DSpecM1D model columns are radius, density, Vpv, Vsv, Qkappa, Qmu,
    # Vph, Vsh, eta. Zero shear speeds mark fluid regions and are excluded.
    solid_speeds = [row[index] for row in rows for index in (3, 7)
                    if row[index] > 0 and math.isfinite(row[index])]
    if not solid_speeds:
        raise SystemExit("model contains no positive solid shear speeds")
    radius = max(row[0] for row in rows)
    return min(solid_speeds), radius


def make_parameters(model_name: str, prefix: Path, settings: dict) -> str:
    values = [
        str(prefix), f"models/{model_name}", "4", "0", "2", "0", "0",
        f"{settings['relative_error']:.17g}",
        str(settings["degree_min"]), str(settings["degree_max"]),
        f"{settings['frequency_min_mHz']:.17g}",
        f"{settings['frequency_max_mHz']:.17g}",
        f"{settings['duration_min']:.17g}", f"{settings['time_step_s']:.17g}",
        *(f"{corner:.17g}" for corner in settings["filter_mHz"]),
        f"{settings['source_depth_km']:.17g}",
        f"{settings['source_geocentric_latitude_deg']:.17g}",
        f"{settings['source_longitude_deg']:.17g}",
        *(f"{MOMENT_NM[key]:.17g}" for key in ("rr", "rt", "rp", "tt", "tp", "pp")),
        f"{settings['receiver_depth_km']:.17g}", "1",
        f"{settings['receiver_geocentric_latitude_deg']:.17g} {settings['receiver_longitude_deg']:.17g}",
    ]
    return "\n".join(values) + "\n"


def make_settings(args, nq: int, epw: float, pilot: bool,
                  minimum_wavelength_m: float, earth_radius_m: float) -> dict:
    settings = dict(PUBLISHED)
    settings["filter_mHz"] = list(PUBLISHED["filter_mHz"])
    settings["mesh_order_points"] = nq
    settings["elements_per_minimum_wavelength"] = epw
    settings["mesh_maxstep_earth_radii"] = minimum_wavelength_m / (
        epw * earth_radius_m)
    if pilot:
        # A very narrow top-frequency sentinel measures fine-mesh cost at full
        # degree and time-grid resolution. It is diagnostic only, not a Fig. 6 result.
        settings["frequency_min_mHz"] = 49.9
        settings["frequency_max_mHz"] = 50.0
        settings["run_scope"] = "49.9-50.0 mHz cost pilot; not a paper result"
    else:
        settings["run_scope"] = "full published-style 0.1-50 mHz solve band"
    return settings


def run_command(command: list[str], log: Path, resource: Path,
                cwd: Path, env: dict[str, str]) -> tuple[float, int | None]:
    timer = Path("/usr/bin/time")
    timed = timer.is_file()
    actual = [str(timer), "-f", "MAX_RSS_KIB=%M", "-o", str(resource), *command] if timed else command
    started = time.perf_counter()
    completed = subprocess.run(actual, cwd=cwd, env=env, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    log.write_text(completed.stdout)
    if completed.returncode:
        print(completed.stdout, end="", flush=True)
        raise subprocess.CalledProcessError(completed.returncode, command)
    rss = None
    if timed and resource.is_file():
        line = resource.read_text().strip()
        if line.startswith("MAX_RSS_KIB="):
            rss = int(line.split("=", 1)[1])
    return elapsed, rss


def read_csv_matrix(path: Path, names: list[str]) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"empty adapter output: {path}")
    time_s = np.asarray([float(row["time_s"]) for row in rows], dtype=float)
    values = np.asarray([[float(row[name]) for name in names] for row in rows],
                        dtype=float).T
    if not np.isfinite(time_s).all() or not np.isfinite(values).all():
        raise SystemExit(f"non-finite values in {path}")
    return time_s, values


def read_spectra(path: Path) -> tuple[dict, dict]:
    spectra, points = {}, {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["frequency_bin"]), int(row["component"]))
            spectra[key] = complex(float(row["real"]), float(row["imag"]))
            points[key] = {
                "frequency_hz": float(row["frequency_hz"]),
                "omega_real_rad_s": float(row["omega_real_rad_s"]),
                "omega_imag_rad_s": float(row["omega_imag_rad_s"]),
            }
    if not spectra:
        raise SystemExit(f"empty raw spectrum: {path}")
    return spectra, points


def read_meta(path: Path) -> dict:
    meta = json.loads(path.read_text())
    for key in ("samples", "output_samples", "sample_interval_s",
                "frequency_spacing_hz", "damping_epsilon_rad_s"):
        if key not in meta:
            raise SystemExit(f"adapter metadata missing {key}: {path}")
    return meta


def half_cosine(frequency_hz: np.ndarray, corners_mhz: list[float]) -> np.ndarray:
    f1, f2, f3, f4 = [v / 1000.0 for v in corners_mhz]
    result = np.zeros_like(frequency_hz)
    rising = (frequency_hz >= f1) & (frequency_hz < f2)
    flat = (frequency_hz >= f2) & (frequency_hz < f3)
    falling = (frequency_hz >= f3) & (frequency_hz < f4)
    result[rising] = 0.5 * (1.0 - np.cos(np.pi * (frequency_hz[rising]-f1)/(f2-f1)))
    result[flat] = 1.0
    result[falling] = 0.5 * (1.0 - np.cos(np.pi * (f4-frequency_hz[falling])/(f4-f3)))
    return result


def inverse_transform(spectra: dict, meta: dict, corners: list[float]) -> tuple[np.ndarray, np.ndarray]:
    nfft = int(meta["samples"])
    count = int(meta["output_samples"])
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    eps = float(meta["damping_epsilon_rad_s"])
    bins = np.zeros((3, nfft // 2 + 1), dtype=np.complex128)
    for (index, component), value in spectra.items():
        bins[component, index] = value
    frequencies = np.arange(nfft // 2 + 1) * df
    bins *= half_cosine(frequencies, corners)[None, :]
    # irfft is the normalized inverse of FFTW's unnormalized c2r transform.
    traces = np.fft.irfft(bins, n=nfft, axis=1)[:, :count]
    time_s = np.arange(count, dtype=float) * dt
    traces *= nfft * df * np.exp(eps * time_s)[None, :]
    if not np.isfinite(traces).all():
        raise SystemExit("independent inverse transform produced non-finite values")
    return time_s, traces


def direct_inverse_checks(spectra: dict, meta: dict, corners: list[float],
                          traces: np.ndarray) -> list[dict]:
    nfft = int(meta["samples"])
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    eps = float(meta["damping_epsilon_rad_s"])
    sample_ids = sorted({0, 1, min(137, traces.shape[1] - 1),
                         traces.shape[1] // 2, traces.shape[1] - 1})
    checks = []
    for component in range(3):
        for sample in sample_ids:
            total = 0.0
            for (index, channel), value in spectra.items():
                if channel != component:
                    continue
                weight = half_cosine(np.asarray([index * df]), corners)[0]
                weighted = value * weight
                phase = 2 * math.pi * index * sample / nfft
                if index == 0 or (nfft % 2 == 0 and index == nfft // 2):
                    total += weighted.real * (-1.0 if index == nfft // 2 and sample % 2 else 1.0)
                else:
                    total += 2 * (weighted * complex(math.cos(phase), math.sin(phase))).real
            expected = df * math.exp(eps * sample * dt) * total
            actual = float(traces[component, sample])
            peak = float(np.max(np.abs(traces[component])))
            error = abs(actual - expected)
            scaled = error / peak if peak else error
            if not math.isfinite(expected) or scaled > 2e-10:
                raise SystemExit(f"direct inverse-transform check failed: component={component} sample={sample} scaled_error={scaled:.3g}")
            checks.append({"component": CHANNELS[component], "sample": sample,
                           "absolute_error_si": error, "error_over_channel_peak": scaled})
    return checks


def eq4_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict:
    diff = candidate - reference
    peak = float(np.max(np.abs(reference)))
    mean_abs = float(np.mean(np.abs(diff)))
    max_abs = float(np.max(np.abs(diff)))
    norm = float(np.linalg.norm(reference))
    return {
        "reference_peak_SI": peak,
        "mean_absolute_difference_SI": mean_abs,
        "mean_error_percent_of_reference_peak": 100 * mean_abs / peak if peak else None,
        "maximum_absolute_difference_SI": max_abs,
        "maximum_error_percent_of_reference_peak": 100 * max_abs / peak if peak else None,
        "relative_L2": float(np.linalg.norm(diff) / norm) if norm else None,
        "normalization": "max absolute value of NQ6 EPW14 filtered channel",
    }


def plot_convergence(path: Path, metrics: dict) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)
    colors = dict(zip(CHANNELS, ("#1f77b4", "#ff7f0e", "#2ca02c")))
    for nq, ax in zip((4, 5, 6), axes):
        for channel in CHANNELS:
            points = sorted((epw, metrics[f"nq{nq}_epw{epw:g}"]["channels"][channel]
                             ["mean_error_percent_of_reference_peak"])
                            for epw in EPW_VALUES if f"nq{nq}_epw{epw:g}" in metrics)
            # The reference compares exactly with itself and has zero error,
            # which is undefined on the logarithmic y-axis.
            if nq == REFERENCE[0]:
                points = [point for point in points if point[0] != REFERENCE[1]]
            if points:
                ax.plot([p[0] for p in points], [p[1] for p in points],
                        marker="o", color=colors[channel], label=channel)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{nq}-point SEM")
        ax.set_xlabel("Elements per minimum wavelength")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("Mean Eq. (4) error (%)")
    axes[-1].legend(title="Component")
    fig.suptitle("Figure 6 spatial convergence vs NQ6, 14 EPW numerical reference")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_traces(path: Path, time_s: np.ndarray, reference: np.ndarray,
                selected: dict[str, np.ndarray]) -> None:
    import matplotlib.pyplot as plt

    selector = (time_s <= 3600.0)
    fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
    for component, label in enumerate(CHANNELS):
        axes[component].plot(time_s[selector], reference[component, selector],
                             color="black", linewidth=1.2, label="NQ6 EPW14 reference")
        for name, values in selected.items():
            axes[component].plot(time_s[selector], values[component, selector],
                                 linewidth=0.8, alpha=0.85, label=name)
        axes[component].set_ylabel(f"{label} displacement (m)")
        axes[component].grid(alpha=0.25)
    axes[-1].set_xlabel("Time (s)")
    axes[0].legend(ncol=2, fontsize=8)
    fig.suptitle("Filtered SI traces, first hour")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--pilot", action="store_true")
    mode.add_argument("--full", action="store_true")
    parser.add_argument("--adapter", type=Path, required=True,
                        help="existing dspec_adapter binary; no build is started")
    parser.add_argument("--dspec-source", type=Path, default=DSPEC_DEFAULT)
    parser.add_argument("--criteria", type=Path, default=CRITERIA_DEFAULT)
    parser.add_argument("--output-root", type=Path,
                        default=DSPEC_DEFAULT / "build/phase3/paper-suite-20261020/fig6-convergence")
    parser.add_argument("--include-fine-control", action="store_true",
                        help="add NQ6 at 28 elements per wavelength after the pilot is reviewed")
    args = parser.parse_args()

    dspec = args.dspec_source.resolve()
    adapter = args.adapter.resolve()
    model_source = dspec / MODEL_REL
    params_source = dspec / PARAM_REL
    if not adapter.is_file():
        parser.error(f"adapter binary does not exist: {adapter}")
    if not model_source.is_file() or not params_source.is_file():
        parser.error("missing required ex5 parameter/model inputs")
    criteria = json.loads(args.criteria.read_text())
    reviewed = {item.get("case") for item in criteria.get("case_inventory", [])}
    if "ex5/Fig6" not in reviewed:
        raise SystemExit("ex5/Fig6 absent from frozen criteria")
    expected_hash = criteria["source_basis"]["source_parameter_hashes"]["ex5"]
    if sha256(params_source) != expected_hash:
        raise SystemExit("data/params/ex5.txt differs from the hash in frozen criteria")
    criteria_hash = sha256(args.criteria)
    min_vs, earth_radius = parse_model_min_s_velocity(model_source)
    minimum_wavelength = min_vs / (PUBLISHED["frequency_max_mHz"] / 1000.0)
    root = args.output_root.resolve()
    run_name = "pilot" if args.pilot else "full"
    run_root = root / run_name
    if run_root.exists():
        raise SystemExit(f"refusing to overwrite existing run directory: {run_root}")
    run_root.mkdir(parents=True)
    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OMP_MAX_ACTIVE_LEVELS"] = "1"
    # Keep common BLAS libraries at one thread, independently of OpenMP.
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = "1"

    cases = [(4, 3.0), REFERENCE] if args.pilot else [
        (nq, epw) for nq in (4, 5, 6) for epw in EPW_VALUES
    ]
    if args.full and args.include_fine_control:
        cases.append((6, 28.0))
    input_dir = run_root / "inputs"
    model_copy = input_dir / "data/models" / model_source.name
    model_copy.parent.mkdir(parents=True)
    shutil.copyfile(model_source, model_copy)
    output_dir = run_root / "results"
    output_dir.mkdir()
    snapshot = run_root / "source_snapshot"
    snapshot.mkdir()
    adapter_copy = snapshot / "dspec_adapter"
    shutil.copyfile(adapter, adapter_copy)
    adapter_hash = sha256(adapter_copy)
    adapter_source = BASE / "dspec_adapter.cpp"
    if not adapter_source.is_file():
        raise SystemExit(f"adapter source is missing: {adapter_source}")
    shutil.copyfile(adapter_source, snapshot / adapter_source.name)
    source_hash = sha256(snapshot / adapter_source.name)
    shutil.copyfile(args.criteria, snapshot / "criteria.json")
    shutil.copyfile(Path(__file__), snapshot / Path(__file__).name)
    input_copy = input_dir / PARAM_REL
    input_copy.parent.mkdir(parents=True)
    shutil.copyfile(params_source, input_copy)

    min_vs_m_s = min_vs
    source_head = subprocess.check_output(
        ["git", "-C", str(dspec), "rev-parse", "HEAD"], text=True).strip()
    env_note = {
        "OMP_NUM_THREADS": env["OMP_NUM_THREADS"],
        "OPENBLAS_NUM_THREADS": env["OPENBLAS_NUM_THREADS"],
        "MKL_NUM_THREADS": env["MKL_NUM_THREADS"],
        "jobs_used_by_this_runner": 1,
    }
    provenance = {
        "case": "ex5/Fig6", "paper_figure": 6,
        "criteria_file": str(args.criteria.resolve()),
        "criteria_sha256": criteria_hash,
        "source_parameter_file": str(params_source),
        "source_parameter_sha256": sha256(params_source),
        "source_parameter_copy": "inputs/data/params/ex5.txt",
        "model_file": str(model_source), "model_sha256": sha256(model_source),
        "model_copy": f"inputs/data/models/{model_source.name}",
        "DSpecM1D_HEAD_at_run": source_head,
        "adapter_binary_sha256": adapter_hash,
        "adapter_binary_snapshot": "source_snapshot/dspec_adapter",
        "adapter_source_sha256": source_hash,
        "adapter_source_snapshot": "source_snapshot/dspec_adapter.cpp",
        "runner_sha256": sha256(Path(__file__)),
        "YSpecpp": None,
        "configuration": {
            **PUBLISHED,
            "source_time_convention": "step moment; no extra source-time filter",
            "coordinates": "source and receiver input latitudes are geocentric",
            "receiver_site": "TLY; rounded paper coordinates (51.68, 103.64)",
            "published_text_caption_filter_difference":
                "Section 3.3 says (0.2,50) mHz; Fig. 6 caption corners are (0.1,0.15)/(45,50) mHz. Runner uses the complete captioned half-cosine taper.",
            "reference_definition": "NQ6 at 14 elements per minimum wavelength; numerical reference, not exact solution",
            "published_convergence_context":
                "about 0.1% mean Eq. (4) error near 3, 2, and 1 EPW for NQ4, NQ5, and NQ6 respectively",
        },
        "mesh_scale": {
            "minimum_positive_shear_speed_m_s": min_vs_m_s,
            "earth_radius_m": earth_radius,
            "maximum_frequency_hz": PUBLISHED["frequency_max_mHz"] / 1000.0,
            "minimum_wavelength_m": minimum_wavelength,
            "derivation": "minimum positive Vsv/Vsh in ocean-less PREM divided by 50 mHz",
        },
        "runtime_environment": env_note,
        "run_scope": "49.9-50.0 mHz cost pilot; not a paper result" if args.pilot else "full sweep",
        "runs": [],
        "status": "RUNNING",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (run_root / "summary.json").write_text(json.dumps(provenance, indent=2) + "\n")

    # DSpec parameter paths are interpreted relative to the copied data tree.
    cwd = input_dir
    run_results = {}
    for nq, epw in cases:
        settings = make_settings(args, nq, epw, args.pilot,
                                 minimum_wavelength, earth_radius)
        label = f"nq{nq}_epw{epw:g}"
        case_out = output_dir / label
        case_out.mkdir()
        prefix = case_out / "dspec"
        params_file = input_dir / "data/params" / f"{label}.txt"
        params_file.write_text(make_parameters(model_copy.name, prefix, settings))
        log = case_out / "dspec.stdout.log"
        resource = case_out / "dspec.resource.txt"
        command = [str(adapter), str(params_file), str(prefix), str(nq),
                   str(PUBLISHED["nskip"]),
                   f"{settings['mesh_maxstep_earth_radii']:.17g}",
                   "paper_suite", "single"]
        wall_s, rss_kib = run_command(command, log, resource, cwd, env)
        raw, points = read_spectra(prefix.with_suffix(".spectra.csv"))
        meta = read_meta(prefix.with_suffix(".meta.json"))
        native_time, native = read_csv_matrix(
            prefix.with_suffix(".timeseries.csv"), ["Z_SI", "N_SI", "E_SI"])
        time_grid, filtered = inverse_transform(raw, meta, settings["filter_mHz"])
        if not np.allclose(native_time, time_grid, rtol=0, atol=1e-9):
            raise SystemExit(f"{label}: native and independent time grids differ")
        if native.shape != filtered.shape:
            raise SystemExit(f"{label}: native and independent trace shapes differ")
        native_peak = float(np.max(np.abs(filtered)))
        native_diff = float(np.max(np.abs(native - filtered)))
        native_ratio = native_diff / native_peak if native_peak else native_diff
        if native_ratio > 2e-10:
            raise SystemExit(f"{label}: native/IFFT check failed ({native_ratio:.3g})")
        direct_checks = direct_inverse_checks(raw, meta, settings["filter_mHz"], filtered)
        if int(meta["degree_max"]) != PUBLISHED["degree_max"]:
            raise SystemExit(f"{label}: unexpected degree cap in adapter metadata")
        run_results[label] = {
            "nq": nq, "elements_per_minimum_wavelength": epw,
            "maxstep_earth_radii": settings["mesh_maxstep_earth_radii"],
            "frequency_band_mHz": [settings["frequency_min_mHz"], settings["frequency_max_mHz"]],
            "raw_spectrum_file": str(prefix.with_suffix(".spectra.csv").relative_to(run_root)),
            "native_filtered_timeseries_file": str(prefix.with_suffix(".timeseries.csv").relative_to(run_root)),
            "metadata_file": str(prefix.with_suffix(".meta.json").relative_to(run_root)),
            "parameter_file": str(params_file.relative_to(run_root)),
            "frequency_component_samples": len(raw),
            "fft_samples": int(meta["samples"]),
            "output_samples": int(meta["output_samples"]),
            "dt_s": float(meta["sample_interval_s"]),
            "df_hz": float(meta["frequency_spacing_hz"]),
            "wall_time_s": wall_s,
            "peak_rss_kib": rss_kib,
            "raw_spectra_finite": all(math.isfinite(v.real) and math.isfinite(v.imag) for v in raw.values()),
            "filtered_native_finite": bool(np.isfinite(native).all()),
            "independent_ifft_max_error_over_trace_peak": native_ratio,
            "direct_ifft_checks": direct_checks,
            "_time": time_grid,
            "_filtered": filtered,
            "_native": native,
        }
        print(f"{label}: {wall_s:.1f}s, RSS={rss_kib if rss_kib is not None else 'n/a'} KiB",
              flush=True)
        # Refresh a compact progress file after every expensive mesh.
        progress = dict(provenance)
        progress["runs"] = [{k: v for k, v in item.items()
                             if not k.startswith("_")} for item in run_results.values()]
        progress["status"] = "RUNNING"
        (run_root / "summary.json").write_text(json.dumps(progress, indent=2) + "\n")

    clean_runs = [{k: v for k, v in item.items() if not k.startswith("_")}
                  for item in run_results.values()]
    metrics = {}
    if args.full:
        ref_label = f"nq{REFERENCE[0]}_epw{REFERENCE[1]:g}"
        ref_time = run_results[ref_label]["_time"]
        ref_trace = run_results[ref_label]["_filtered"]
        for label, item in run_results.items():
            if not np.array_equal(item["_time"], ref_time):
                raise SystemExit(f"{label}: time samples do not match reference")
            channels = {}
            for index, channel in enumerate(CHANNELS):
                channels[channel] = eq4_metrics(ref_trace[index], item["_filtered"][index])
            metrics[label] = {
                "nq": item["nq"],
                "elements_per_minimum_wavelength": item["elements_per_minimum_wavelength"],
                "reference": ref_label,
                "channels": channels,
                "paper_target_case": (item["nq"], item["elements_per_minimum_wavelength"]) in PAPER_TARGETS,
            }
        plot_convergence(run_root / "fig6_convergence.png", metrics)
        selected = {
            "NQ4 EPW3": run_results["nq4_epw3"]["_filtered"],
            "NQ5 EPW2": run_results["nq5_epw2"]["_filtered"],
            "NQ6 EPW1": run_results["nq6_epw1"]["_filtered"],
        }
        plot_traces(run_root / "fig6_reference_traces.png", ref_time, ref_trace, selected)
    provenance["runs"] = clean_runs
    provenance["metrics_vs_reference"] = metrics
    provenance["plots"] = (["fig6_convergence.png", "fig6_reference_traces.png"]
                           if args.full else [])
    provenance["status"] = "COMPLETE"
    provenance["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (run_root / "summary.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Saved {run_root / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
