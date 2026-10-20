#!/usr/bin/env python3
"""Run or reprocess the four paired paper benchmark cases (Figures 1–4)."""
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

BASE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(BASE))
from run_benchmark import (  # noqa: E402
    PINNED_DEPS, YSPEC_ARCHIVE_SHA256, YSPEC_COMMIT, YSPEC_FILETREE_SHA256,
    filetree_sha256, geographic_latitude, inverse_samples, compare_vectors, read_meta, read_spectra,
    read_timeseries, verify_dependencies,
)

CASES = json.loads((BASE / "cases/paper-cases.json").read_text())["cases"]
UNITS = {"displacement": "m", "velocity": "m/s", "acceleration": "m/s^2"}
COMPONENTS = ("Z", "N", "E")
MOMENT_KEYS = ("rr", "rt", "rp", "tt", "tp", "pp")


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def run(
    command: list[str],
    *,
    log: Path | None = None,
    env: dict[str, str] | None = None,
) -> float:
    started = time.perf_counter()
    if log is None:
        subprocess.run(command, check=True, env=env)
    else:
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            env=env,
        )
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(result.stdout)
        if result.returncode:
            print(result.stdout, end="", flush=True)
            raise subprocess.CalledProcessError(result.returncode, command)
    return time.perf_counter() - started


def verify_planetmodel(build_dir: Path) -> dict:
    checkout = build_dir / "_deps" / "planetmodel-src"
    expected = "fa6e230fde25152d40f27d9839de44c4f3182a24"
    if not (checkout / ".git").exists():
        raise SystemExit("YSpecpp build lacks its pinned planetmodel checkout")
    commit = subprocess.check_output(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        [
            "git", "-C", str(checkout), "status", "--porcelain",
            "--untracked-files=normal",
        ],
        text=True,
    )
    if commit != expected or dirty:
        raise SystemExit("planetmodel checkout differs from the pinned clean commit")
    return {"commit": commit, "clean": True}


def parameter_text(case: dict, model: str, prefix: Path, settings: dict) -> str:
    gravity_code = {"none": 0, "cowling": 1, "full": 2}[case["gravity"]]
    output_code = {
        "displacement": 0,
        "velocity": 1,
        "acceleration": 2,
    }[case["quantity"]]
    moment = case["moment_tensor_Nm"]
    receivers = case["receivers_geocentric"]
    values = [
        str(prefix),
        f"models/{model}",
        "4",
        "1" if case["attenuation"] else "0",
        str(gravity_code),
        str(output_code),
        "0",
        f"{case['dspec']['relative_error']:.17g}",
        str(case["degree_min"]),
        str(case["degree_max"]),
        f"{settings['frequency_min_mHz']:.17g}",
        f"{settings['frequency_max_mHz']:.17g}",
        f"{case['duration_min']:.17g}",
        f"{case['time_step_s']:.17g}",
        *(f"{corner:.17g}" for corner in case["filter_mHz"]),
        f"{case['source_depth_km']:.17g}",
        f"{case['source_geocentric_latitude_deg']:.17g}",
        f"{case['source_longitude_deg']:.17g}",
        *(f"{moment[key]:.17g}" for key in MOMENT_KEYS),
        f"{case['receiver_depth_km']:.17g}",
        str(len(receivers)),
        *(f"{latitude:.17g} {longitude:.17g}"
          for latitude, longitude in receivers),
    ]
    return "\n".join(values) + "\n"


def yspec_input(case: dict, model: Path, settings: dict) -> dict:
    source_latitude = geographic_latitude(
        case["source_geocentric_latitude_deg"]
    )
    receiver_locations = [
        [geographic_latitude(latitude), longitude]
        for latitude, longitude in case["receivers_geocentric"]
    ]
    return {
        "model": str(model.resolve()),
        "attenuation": case["attenuation"],
        "gravity": case["gravity"],
        "output": {
            "prefix": "unused",
            "quantity": case["quantity"],
            "potential_and_tilt_corrections": False,
        },
        "degrees": {
            "min": case["degree_min"],
            "max": case["degree_max"],
        },
        "frequencies_mHz": {
            "min": settings["frequency_min_mHz"],
            "max": settings["frequency_max_mHz"],
        },
        "time_series": {
            "length_min": case["duration_min"],
            "time_step_s": case["time_step_s"],
        },
        "filter_mHz": case["filter_mHz"],
        "source": {
            "depth_km": case["source_depth_km"],
            "latitude": source_latitude,
            "longitude": case["source_longitude_deg"],
            "moment_tensor_Nm": case["moment_tensor_Nm"],
        },
        "receivers": {
            "depth_km": case["receiver_depth_km"],
            "locations": receiver_locations,
        },
        "solver": case["yspec"],
    }


def half_cosine(freq: np.ndarray, corners_mhz: list[float]) -> np.ndarray:
    f1, f2, f3, f4 = np.asarray(corners_mhz, dtype=float) / 1000.0
    weights = np.zeros(freq.shape, dtype=float)
    rising = (freq >= f1) & (freq < f2)
    flat = (freq >= f2) & (freq < f3)
    falling = (freq >= f3) & (freq < f4)
    weights[rising] = 0.5 * (
        1 - np.cos(np.pi * (freq[rising] - f1) / (f2 - f1))
    )
    weights[flat] = 1.0
    weights[falling] = 0.5 * (
        1 - np.cos(np.pi * (f4 - freq[falling]) / (f4 - f3))
    )
    return weights


def physical_inverse(
    spectra: dict,
    meta: dict,
    corners: list[float],
    samples: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the filter once, then undo spectral damping once."""
    nfft = int(meta["samples"])
    output_samples = int(
        meta["output_samples"] if samples is None else samples
    )
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    frequency = np.arange(nfft // 2 + 1) * df
    bins = np.zeros((3, nfft // 2 + 1), dtype=np.complex128)
    for (index, component), value in spectra.items():
        bins[component, index] = value

    filtered = bins * half_cosine(frequency, corners)[None, :]
    traces = np.fft.irfft(filtered, n=nfft, axis=1)
    time = np.arange(output_samples, dtype=float) * dt
    traces = traces[:, :output_samples] * (
        nfft * df * np.exp(epsilon * time)[None, :]
    )
    if not np.isfinite(traces).all():
        raise SystemExit("physical inverse transform produced non-finite values")
    return time, traces


def hann_ends(time: np.ndarray, duration: float, fraction: float) -> np.ndarray:
    edge = fraction * duration
    window = np.ones(time.size)
    if edge <= 0:
        return window
    left = time < edge
    right = time > duration - edge
    window[left] = 0.5 * (
        1 - np.cos(np.pi * time[left] / edge)
    )
    window[right] = 0.5 * (
        1 - np.cos(np.pi * (duration - time[right]) / edge)
    )
    window[time >= duration] = 0.0
    return window


def processed_spectrum(
    spectra: dict,
    meta: dict,
    corners: list[float],
    final_fraction: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Figure 1 chain: physical trace, 5% Hann, filter, 50% Hann."""
    nfft = int(meta["samples"])
    output_samples = int(meta["output_samples"])
    dt = float(meta["sample_interval_s"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    bins = np.zeros(nfft // 2 + 1, dtype=np.complex128)
    for (index, component), value in spectra.items():
        if component == 0:
            bins[index] = value

    frequency = np.arange(nfft // 2 + 1) * df
    time = np.arange(nfft) * dt
    duration = max((output_samples - 1) * dt, dt)
    physical = np.fft.irfft(bins, n=nfft) * nfft * df
    physical *= np.exp(epsilon * time)
    pre_window_spectrum = np.fft.rfft(
        physical * hann_ends(time, duration, 0.05), n=nfft
    )
    filtered_time = np.fft.irfft(
        pre_window_spectrum * half_cosine(frequency, corners), n=nfft
    ) * nfft * df
    processed = np.fft.rfft(
        filtered_time * hann_ends(time, duration, final_fraction), n=nfft
    )
    return frequency, np.abs(processed)


def metric(reference: np.ndarray, candidate: np.ndarray, unit: str) -> dict:
    difference = candidate - reference
    reference_peak = float(np.max(np.abs(reference)))
    comparison = compare_vectors(
        [complex(value) for value in reference],
        [complex(value) for value in candidate],
    )
    mean_difference = float(np.mean(np.abs(difference)))
    maximum_difference = float(np.max(np.abs(difference)))
    return {
        "reference_peak": reference_peak,
        "unit": unit,
        "mean_absolute_difference": mean_difference,
        "mean_difference_percent_of_reference_peak": (
            100 * mean_difference / reference_peak if reference_peak else None
        ),
        "maximum_absolute_difference": maximum_difference,
        "maximum_difference_percent_of_reference_peak": (
            100 * maximum_difference / reference_peak if reference_peak else None
        ),
        "relative_l2": comparison["relative_l2"],
    }


def spectrum_metric(
    reference: dict,
    candidate: dict,
    points: dict,
    band_mhz: list[float],
) -> dict:
    lower_hz, upper_hz = (value / 1000 for value in band_mhz)
    keys = [
        key for key in sorted(reference)
        if key[1] == 0
        and lower_hz <= points[key]["frequency_hz"] <= upper_hz
    ]
    candidate_keys = {
        key for key in candidate
        if key[1] == 0
        and lower_hz <= points[key]["frequency_hz"] <= upper_hz
    }
    if not keys or set(keys) != candidate_keys:
        raise SystemExit("raw spectrum bins differ or empty in comparison band")

    reference_amplitudes = np.asarray([abs(reference[key]) for key in keys])
    candidate_amplitudes = np.asarray([abs(candidate[key]) for key in keys])
    difference = np.abs(candidate_amplitudes - reference_amplitudes)
    reference_peak = float(np.max(reference_amplitudes))
    mean_difference = float(np.mean(difference))
    maximum_difference = float(np.max(difference))
    return {
        "frequency_band_mHz": band_mhz,
        "frequency_bins": len(keys),
        "unit": "SI observable / Hz",
        "mean_peak_normalized_percent": (
            100 * mean_difference / reference_peak if reference_peak else None
        ),
        "maximum_peak_normalized_percent": (
            100 * maximum_difference / reference_peak if reference_peak else None
        ),
    }


def plot_timeseries(
    path: Path,
    time: np.ndarray,
    dspec: np.ndarray,
    yspec: np.ndarray,
    quantity: str,
    window: list[float] | None,
    metrics: dict,
) -> None:
    import matplotlib.pyplot as plt

    if window is None:
        selected = np.ones(time.size, dtype=bool)
    else:
        selected = (time >= window[0]) & (time <= window[1])

    figure, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
    for component, label in enumerate(COMPONENTS):
        peak = metrics[label]["reference_peak"] or 1.0
        axes[component, 0].plot(
            time[selected], dspec[component, selected] / peak, label="DSpecM1D"
        )
        axes[component, 0].plot(
            time[selected], yspec[component, selected] / peak,
            "--", label="Pinned YSpecpp",
        )
        axes[component, 0].set_ylabel(f"{label} / YSpecpp peak")
        axes[component, 1].plot(
            time[selected],
            100 * (dspec[component, selected] - yspec[component, selected]) / peak,
            color="black",
        )
        axes[component, 1].set_ylabel("difference (% peak)")
        axes[component, 0].grid(alpha=0.25)
        axes[component, 1].grid(alpha=0.25)

    figure.suptitle(f"Filtered {quantity}")
    axes[0, 0].legend()
    axes[2, 0].set_xlabel("Time (s)")
    axes[2, 1].set_xlabel("Time (s)")
    figure.supylabel(UNITS[quantity])
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_raw(
    path: Path,
    points: dict,
    dspec: dict,
    yspec: dict,
    band_mhz: list[float],
    unit: str,
) -> None:
    import matplotlib.pyplot as plt

    lower_hz, upper_hz = (value / 1000 for value in band_mhz)
    keys = [
        key for key in sorted(dspec)
        if key[1] == 0
        and lower_hz <= points[key]["frequency_hz"] <= upper_hz
    ]
    frequency_mhz = np.asarray(
        [points[key]["frequency_hz"] * 1000 for key in keys]
    )
    dspec_amplitude = np.asarray([abs(dspec[key]) for key in keys])
    yspec_amplitude = np.asarray([abs(yspec[key]) for key in keys])
    yspec_peak = float(np.max(yspec_amplitude))

    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(frequency_mhz, dspec_amplitude, label="DSpecM1D")
    axes[0].plot(frequency_mhz, yspec_amplitude, "--", label="Pinned YSpecpp")
    axes[0].set_yscale("log")
    axes[0].set_ylabel(f"|Z| ({unit}/Hz)")
    axes[0].legend()

    if yspec_peak:
        difference = 100 * (dspec_amplitude - yspec_amplitude) / yspec_peak
    else:
        difference = dspec_amplitude - yspec_amplitude
    axes[1].plot(frequency_mhz, difference, color="black")
    axes[1].set_ylabel("difference (% YSpecpp peak)")
    axes[1].set_xlabel("Frequency (mHz)")
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_processed(
    path: Path,
    frequency_hz: np.ndarray,
    dspec: np.ndarray,
    yspec: np.ndarray,
    band_mhz: list[float],
) -> None:
    import matplotlib.pyplot as plt

    lower_hz, upper_hz = (value / 1000 for value in band_mhz)
    selected = (frequency_hz >= lower_hz) & (frequency_hz <= upper_hz)
    yspec_peak = float(np.max(yspec[selected]))

    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    frequency_mhz = frequency_hz[selected] * 1000
    axes[0].plot(frequency_mhz, dspec[selected], label="DSpecM1D")
    axes[0].plot(frequency_mhz, yspec[selected], "--", label="Pinned YSpecpp")
    axes[0].set_ylabel("Processed |Z| (common scale)")
    axes[0].legend()

    if yspec_peak:
        difference = 100 * (dspec[selected] - yspec[selected]) / yspec_peak
    else:
        difference = dspec[selected] - yspec[selected]
    axes[1].plot(frequency_mhz, difference, color="black")
    axes[1].set_ylabel("difference (% peak)")
    axes[1].set_xlabel("Frequency (mHz)")
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def outputs(result_dir: Path, solver: str) -> tuple[dict, dict, dict, np.ndarray, np.ndarray]:
    prefix = result_dir / solver
    spectra, points = read_spectra(prefix.with_suffix(".spectra.csv"))
    meta = read_meta(prefix.with_suffix(".meta.json"))
    times, channels = read_timeseries(prefix.with_suffix(".timeseries.csv"))
    return spectra, points, meta, np.asarray(times), np.asarray(channels)


def selected_sample_oracle(
    spectra: dict,
    meta: dict,
    corners: list[float],
    time: np.ndarray,
    native: np.ndarray,
) -> list[dict]:
    """Use the shared sparse-bin DFT oracle at the original five samples."""
    sample_indices = sorted({
        0,
        1,
        min(137, len(time) - 1),
        len(time) // 2,
        len(time) - 1,
    })
    frequency_spacing = float(meta["frequency_spacing_hz"])
    filtered_spectra = {
        (index, component): value * half_cosine(
            np.asarray([index * frequency_spacing]), corners
        )[0]
        for (index, component), value in spectra.items()
    }
    sample_times = [float(time[index]) for index in sample_indices]
    expected = inverse_samples(
        filtered_spectra, sample_times, meta, "paper selected samples"
    )

    checks = []
    for component in range(3):
        channel_peak = float(np.max(np.abs(native[component])))
        for column, sample_index in enumerate(sample_indices):
            expected_value = float(expected[component][column])
            actual_value = float(native[component, sample_index])
            error = abs(expected_value - actual_value)
            scaled_error = error / channel_peak if channel_peak else error
            limit = 2e-10 if channel_peak else 1e-20
            if (not math.isfinite(expected_value)
                    or not math.isfinite(actual_value)
                    or scaled_error > limit):
                raise SystemExit(
                    "selected inverse oracle failed at component "
                    f"{component}, sample {sample_index}: {scaled_error:.3g}"
                )
            checks.append({
                "component": component,
                "sample": sample_index,
                "expected": expected_value,
                "actual": actual_value,
                "error_over_channel_peak": scaled_error,
            })
    return checks


def process(
    case_id: str,
    case: dict,
    result_dir: Path,
    output: Path,
    parent_provenance: dict | None = None,
) -> dict:
    dspec, dspec_points, dspec_meta, dspec_time, dspec_native = outputs(
        result_dir, "dspec"
    )
    yspec, yspec_points, yspec_meta, yspec_time, yspec_native = outputs(
        result_dir, "yspecpp"
    )
    expected_quantity = case["quantity"]
    if (dspec_meta["observable"].split()[0] != expected_quantity
            or yspec_meta["observable"].split()[0] != expected_quantity):
        raise SystemExit("saved observable differs from case configuration")
    if (not np.array_equal(dspec_time, yspec_time)
            or len(dspec_time) != int(dspec_meta["output_samples"])
            or len(yspec_time) != int(yspec_meta["output_samples"])):
        raise SystemExit("native time grids differ or metadata sample counts are wrong")

    for key in ("samples", "output_samples"):
        if int(dspec_meta[key]) != int(yspec_meta[key]):
            raise SystemExit(f"FFT grid mismatch for {key}")
    for key in (
        "sample_interval_s",
        "frequency_spacing_hz",
        "damping_epsilon_rad_s",
    ):
        if not math.isclose(
            float(dspec_meta[key]),
            float(yspec_meta[key]),
            rel_tol=1e-12,
            abs_tol=1e-14,
        ):
            raise SystemExit(f"FFT grid mismatch for {key}")

    settings_match = (
        dspec_meta["gravity_selector"] == case["gravity"]
        and yspec_meta["gravity_selector"] == case["gravity"]
        and bool(dspec_meta["attenuation"]) == case["attenuation"]
        and bool(yspec_meta["attenuation"]) == case["attenuation"]
        and [dspec_meta["degree_min"], dspec_meta["degree_max"]]
        == [case["degree_min"], case["degree_max"]]
        and [yspec_meta["degree_min"], yspec_meta["degree_max"]]
        == [case["degree_min"], case["degree_max"]]
        and int(dspec_meta["nq"]) == int(case["dspec"]["nq"])
        and math.isclose(
            float(dspec_meta["relative_error"]),
            float(case["dspec"]["relative_error"]),
            rel_tol=1e-12,
        )
        and math.isclose(
            float(yspec_meta["solver_tolerance"]),
            float(case["yspec"]["tolerance"]),
            rel_tol=1e-12,
        )
    )
    if not settings_match:
        raise SystemExit("saved solver settings differ from paper case configuration")

    for metadata in (dspec_meta, yspec_meta):
        if metadata["spectra_filtered"] or not metadata["native_time_series_filtered"]:
            raise SystemExit("expected raw spectra and native filtered traces")
        if not np.allclose(
            metadata["filter_mHz"], case["filter_mHz"], rtol=1e-12, atol=1e-14
        ):
            raise SystemExit("saved native filter differs from paper case")

    if set(dspec_points) != set(yspec_points):
        raise SystemExit("frequency/component bins differ")
    for key in dspec_points:
        for field in ("frequency_hz", "omega_real_rad_s", "omega_imag_rad_s"):
            if not math.isclose(
                dspec_points[key][field],
                yspec_points[key][field],
                rel_tol=1e-10,
                abs_tol=1e-14,
            ):
                raise SystemExit(f"frequency grid mismatch at {key} ({field})")

    spectrum_band = case.get(
        "spectrum_plot_band_mHz",
        [case["frequency_min_mHz"], case["frequency_max_mHz"]],
    )
    unit = UNITS[expected_quantity]
    dspec_time, dspec_shared = physical_inverse(
        dspec, dspec_meta, case["filter_mHz"]
    )
    yspec_time, yspec_shared = physical_inverse(
        yspec, yspec_meta, case["filter_mHz"]
    )
    if not np.array_equal(dspec_time, yspec_time):
        raise SystemExit("physical inverse time axes differ")

    inverse_checks = {}
    for label, spectra, metadata, shared, native in (
        ("dspec", dspec, dspec_meta, dspec_shared, dspec_native),
        ("yspecpp", yspec, yspec_meta, yspec_shared, yspec_native),
    ):
        if native.shape != shared.shape:
            raise SystemExit(f"{label}: native/shared shape mismatch")
        channel_peak = float(np.max(np.abs(shared)))
        maximum_error = float(np.max(np.abs(native - shared)))
        scaled_error = maximum_error / channel_peak if channel_peak else maximum_error
        limit = 2e-10 if channel_peak else 1e-20
        if not math.isfinite(scaled_error) or scaled_error > limit:
            raise SystemExit(
                f"{label}: native/physical inverse mismatch {scaled_error:.3g}"
            )
        # Preserve the paper runners' separate selected-sample and native/shared
        # 2e-10 peak-scaled gates; the core's 1e-12 gate is unchanged.
        inverse_checks[label] = selected_sample_oracle(
            spectra, metadata, case["filter_mHz"], dspec_time, shared
        )

    windows = {"full_output": None}
    display_window = case.get("display_time_window_s")
    if display_window:
        windows["display_window"] = display_window
    time_metrics = {}
    for window_name, window in windows.items():
        if window is None:
            selected = np.ones(dspec_time.size, dtype=bool)
        else:
            selected = (dspec_time >= window[0]) & (dspec_time <= window[1])
        time_metrics[window_name] = {
            label: metric(
                yspec_shared[index, selected],
                dspec_shared[index, selected],
                unit,
            )
            for index, label in enumerate(COMPONENTS)
        }

    raw_spectrum = spectrum_metric(yspec, dspec, yspec_points, spectrum_band)
    processed_frequency_d, processed_dspec = processed_spectrum(
        dspec, dspec_meta, case["filter_mHz"]
    )
    processed_frequency_y, processed_yspec = processed_spectrum(
        yspec, yspec_meta, case["filter_mHz"]
    )
    if not np.array_equal(processed_frequency_d, processed_frequency_y):
        raise SystemExit("processed spectrum grids differ")

    lower_hz, upper_hz = (value / 1000 for value in spectrum_band)
    spectrum_selected = (
        (processed_frequency_d >= lower_hz)
        & (processed_frequency_d <= upper_hz)
    )
    processed_difference = np.abs(
        processed_dspec[spectrum_selected] - processed_yspec[spectrum_selected]
    )
    processed_peak = float(np.max(processed_yspec[spectrum_selected]))
    mean_processed_difference = float(np.mean(processed_difference))
    maximum_processed_difference = float(np.max(processed_difference))
    processed_spectrum_metrics = {
        "quantity": "processed Z amplitude",
        "frequency_band_mHz": spectrum_band,
        "filter_mHz": case["filter_mHz"],
        "pre_window_fraction": 0.05,
        "final_window_fraction": 0.5,
        "mean_peak_normalized_percent": (
            100 * mean_processed_difference / processed_peak if processed_peak else None
        ),
        "maximum_peak_normalized_percent": (
            100 * maximum_processed_difference / processed_peak
            if processed_peak else None
        ),
        "frequency_bins": int(np.count_nonzero(spectrum_selected)),
    }

    output.mkdir(parents=True, exist_ok=True)
    plot_raw(
        output / "raw_spectrum.png",
        yspec_points,
        dspec,
        yspec,
        spectrum_band,
        unit,
    )
    plot_processed(
        output / "processed_spectrum.png",
        processed_frequency_d,
        processed_dspec,
        processed_yspec,
        spectrum_band,
    )
    plot_timeseries(
        output / "comparison_full.png",
        dspec_time,
        dspec_shared,
        yspec_shared,
        expected_quantity,
        None,
        time_metrics["full_output"],
    )
    if "display_window" in time_metrics:
        plot_timeseries(
            output / "comparison_display.png",
            dspec_time,
            dspec_shared,
            yspec_shared,
            expected_quantity,
            windows["display_window"],
            time_metrics["display_window"],
        )

    summary = {
        "case": case_id,
        "paper_figure": case["paper_figure"],
        "run_scope": case.get("run_scope", "full published band"),
        "quantity": expected_quantity,
        "unit": unit,
        "spectrum_unit": unit + "/Hz",
        "filter_mHz": case["filter_mHz"],
        "time_filter_period_window_s": case.get("time_filter_period_window_s"),
        "grid_samples": int(dspec_meta["samples"]),
        "grid_checks": {
            "frequency_component_bins": len(dspec_points),
            "native_time_samples": len(dspec_time),
        },
        "native_inverse_checks": inverse_checks,
        "time_metrics": time_metrics,
        "raw_spectrum_metrics": {"Z": raw_spectrum},
        "processed_spectrum_metrics": processed_spectrum_metrics,
        "parent_provenance_link": parent_provenance,
        "artifacts": {
            "raw_spectrum": "raw_spectrum.png",
            "processed_spectrum": "processed_spectrum.png",
            "full_time": "comparison_full.png",
            "display_time": (
                "comparison_display.png" if "display_window" in time_metrics else None
            ),
        },
    }
    if case_id == "fig2":
        maximum_mean = max(
            row["mean_difference_percent_of_reference_peak"]
            for row in time_metrics["full_output"].values()
        )
        full_run = case.get("run_scope", "full published band") not in (
            "bounded pilot", "pilot", "diagnostic"
        )
        summary["fig2_mean_gate"] = {
            "applicable": full_run,
            "limit_percent_exclusive": 1.0,
            "maximum_channel_mean_percent": maximum_mean,
            "passed": maximum_mean < 1.0 if full_run else None,
        }
        if full_run and maximum_mean >= 1.0:
            raise SystemExit(
                f"Fig. 2 per-channel mean gate failed: {maximum_mean:g}%"
            )

    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def validate_saved_identity(case_id: str, provenance: dict) -> None:
    """Reject saved evidence that identifies a different configured source/case."""
    case = CASES[case_id]
    source_hash = provenance.get("source_parameter_sha256")
    source_hash = source_hash or provenance.get("moment_tensor_input_sha256")
    if source_hash != case["source_parameter_sha256"]:
        raise SystemExit(f"saved provenance source hash does not match {case_id}")
    declared_case = provenance.get("case")
    if declared_case is not None and declared_case != case_id:
        raise SystemExit(f"saved provenance case {declared_case} does not match {case_id}")
    declared_figure = provenance.get("paper_figure")
    if declared_figure is not None and declared_figure != case["paper_figure"]:
        raise SystemExit(
            f"saved provenance figure {declared_figure} does not match {case_id}"
        )


def find_saved(root: Path, case_id: str) -> tuple[Path, Path, Path]:
    candidates = [
        root,
        root / "full",
        root / "pilot",
        root / case_id / "full",
        root / case_id / "pilot",
    ]
    for candidate in candidates:
        results = candidate / "results"
        if ((results / "dspec.spectra.csv").is_file()
                and (results / "yspecpp.spectra.csv").is_file()):
            provenance_path = candidate / "provenance.json"
            if not provenance_path.is_file():
                raise SystemExit(f"saved run lacks provenance: {candidate}")
            validate_saved_identity(case_id, json.loads(provenance_path.read_text()))
            return candidate, results, provenance_path

    matches = []
    for spectrum_file in root.rglob("dspec.spectra.csv"):
        results = spectrum_file.parent
        if (results.name != "results"
                or not (results / "yspecpp.spectra.csv").is_file()):
            continue

        candidate = results.parent
        provenance_path = candidate / "provenance.json"
        provenance = (
            json.loads(provenance_path.read_text())
            if provenance_path.is_file() else {}
        )
        try:
            validate_saved_identity(case_id, provenance)
        except SystemExit:
            continue
        if provenance:
            matches.append((candidate, results, provenance_path))

    preferred = [item for item in matches if item[0].name == "full"]
    choices = preferred or matches
    if len(choices) == 1:
        return choices[0]
    if not choices:
        raise SystemExit(f"saved run lacks paired raw spectra for {case_id}: {root}")
    rendered = "; ".join(str(item[0]) for item in choices)
    raise SystemExit(
        f"multiple saved runs match {case_id}; specify one exact run directory: {rendered}"
    )


def write_new_run_provenance(
    output: Path,
    source: Path,
    source_parameter: Path,
    source_parameter_copy: Path,
    model: Path,
    settings: dict,
    case: dict,
    yspec_tree: str,
    yspec_file_count: int,
    planetmodel: dict,
    environment: dict,
    adapter_mode: str,
) -> dict:
    dspec_head = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    tracked_changes = subprocess.check_output(
        ["git", "-C", str(source), "diff", "--binary", "HEAD", "--", "DSpecM1D"]
    )
    return {
        "mode": "new paper run",
        "case": case["case"],
        "run_scope": case["run_scope"],
        "settings": settings,
        "case_configuration": case,
        "source_parameter_sha256": digest(source_parameter),
        "source_parameter_snapshot": str(source_parameter_copy.relative_to(output)),
        "model_sha256": digest(model),
        "DSpecM1D_HEAD": dspec_head,
        "DSpecM1D_tracked_code_diff_sha256": hashlib.sha256(
            tracked_changes
        ).hexdigest(),
        "YSpecpp_commit": YSPEC_COMMIT,
        "YSpecpp_tree_sha256": yspec_tree,
        "YSpecpp_tree_file_count": yspec_file_count,
        "YSpecpp_archive_sha256": YSPEC_ARCHIVE_SHA256,
        "dependency_shas": PINNED_DEPS,
        "planetmodel": planetmodel,
        "geocentric_to_geographic": (
            "YSpecpp latitude is converted at the input boundary using "
            "run_benchmark.geographic_latitude"
        ),
        "selected_adapter_mode": adapter_mode,
        "runtime_environment": {
            key: environment[key]
            for key in ("OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_MAX_ACTIVE_LEVELS")
        },
    }


def snapshot_harness(
    output: Path,
    dspec_build: Path,
    yspec_build: Path,
    provenance: dict,
) -> None:
    snapshot = output / "source_snapshot"
    snapshot.mkdir()
    source_names = (
        "run_paper.py",
        "run_benchmark.py",
        "dspec_adapter.cpp",
        "yspecpp_adapter.cpp",
        "CMakeLists.txt",
    )
    for name in source_names:
        source_path = BASE / name
        if source_path.is_file():
            shutil.copyfile(source_path, snapshot / name)

    for name in ("dspec_adapter", "yspecpp_adapter"):
        build = dspec_build if name == "dspec_adapter" else yspec_build
        binary = build / name
        if binary.is_file():
            shutil.copyfile(binary, snapshot / name)

    case_file = BASE / "cases" / "paper-cases.json"
    if case_file.is_file():
        shutil.copyfile(case_file, snapshot / case_file.name)
    readme = BASE / "README.md"
    if readme.is_file():
        shutil.copyfile(readme, snapshot / readme.name)

    provenance["harness_source_sha256"] = {
        path.name: digest(path)
        for path in snapshot.iterdir()
        if path.is_file() and path.suffix in {".py", ".cpp", ".txt", ".json", ".md"}
    }
    provenance["adapter_binary_sha256"] = {
        name: digest(snapshot / name)
        for name in ("dspec_adapter", "yspecpp_adapter")
        if (snapshot / name).is_file()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=tuple(CASES),
        help="omit with --saved-run to reprocess all four cases",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--saved-run",
        type=Path,
        help="reprocess immutable existing spectra and time-series data",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="use and label the configured bounded pilot band",
    )
    parser.add_argument("--dspec-source", type=Path)
    parser.add_argument("--yspec-source", type=Path)
    parser.add_argument("--dependency-root", type=Path)
    args = parser.parse_args()

    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise SystemExit(f"output directory must be empty: {output}")

    if args.saved_run:
        if args.pilot:
            parser.error("--pilot applies only to new runs")
        case_ids = [args.case] if args.case else list(CASES)
        saved_runs = {
            case_id: find_saved(args.saved_run.resolve(), case_id)
            for case_id in case_ids
        }
        for case_id, (parent, results, provenance_path) in saved_runs.items():
            parent_provenance = (
                json.loads(provenance_path.read_text())
                if provenance_path.is_file() else {}
            )
            case = dict(CASES[case_id])
            case["run_scope"] = parent_provenance.get(
                "run_scope", "full published band"
            )
            saved_settings = parent_provenance.get("settings", {})
            for key in ("frequency_min_mHz", "frequency_max_mHz"):
                case[key] = saved_settings.get(key, case[key])

            case_output = output if args.case else output / case_id
            case_output.mkdir(parents=True, exist_ok=True)
            link = {
                "path": str(parent),
                "provenance_sha256": (
                    digest(provenance_path) if provenance_path.is_file() else None
                ),
                "raw_spectra_sha256": {
                    name: digest(results / f"{name}.spectra.csv")
                    for name in ("dspec", "yspecpp")
                },
            }
            process(case_id, case, results, case_output, link)
            record = {
                "mode": "saved-run reprocessing",
                "parent": link,
                "parent_provenance": parent_provenance,
                "case_configuration": case,
            }
            (case_output / "provenance.json").write_text(
                json.dumps(record, indent=2) + "\n"
            )

        print(json.dumps({
            "output": str(output),
            "mode": "saved-run reprocessing",
            "cases": case_ids,
        }, indent=2))
        return

    if args.case is None:
        parser.error("--case is required for a new run")
    source_arguments = (args.dspec_source, args.yspec_source, args.dependency_root)
    if any(path is None for path in source_arguments):
        parser.error(
            "new runs require explicit --dspec-source, --yspec-source, "
            "and --dependency-root"
        )

    case = dict(CASES[args.case])
    case["dspec"] = dict(case["dspec"])
    case["yspec"] = dict(case["yspec"])
    case["case"] = args.case
    case["run_scope"] = "bounded pilot" if args.pilot else "full published band"
    if args.pilot:
        pilot_band = case.get("pilot_frequency_band_mHz")
        if not pilot_band:
            raise SystemExit(f"{args.case} has no predeclared pilot band")
        settings = {
            "frequency_min_mHz": pilot_band[0],
            "frequency_max_mHz": pilot_band[1],
        }
    else:
        settings = {
            "frequency_min_mHz": case["frequency_min_mHz"],
            "frequency_max_mHz": case["frequency_max_mHz"],
        }

    source = args.dspec_source.resolve()
    yspec_source = args.yspec_source.resolve()
    dependency_root = args.dependency_root.resolve()
    for path, label in ((source, "DSpecM1D"), (yspec_source, "YSpecpp")):
        if not path.is_dir():
            parser.error(f"missing {label} source: {path}")

    source_parameter = source / case["source_parameter_file"]
    if digest(source_parameter) != case["source_parameter_sha256"]:
        raise SystemExit("source parameter SHA differs from paper case configuration")
    verify_dependencies(dependency_root)
    yspec_tree, yspec_file_count = filetree_sha256(yspec_source)
    if yspec_tree != YSPEC_FILETREE_SHA256 or yspec_file_count != 76:
        raise SystemExit("YSpecpp source differs from its pinned snapshot")

    output.mkdir(parents=True, exist_ok=True)
    model = source / case["model_file"]
    input_root = output / "inputs"
    model_copy = input_root / "data" / "models" / model.name
    model_copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(model, model_copy)

    results = output / "results"
    results.mkdir()
    source_parameter_copy = input_root / "source_parameters" / source_parameter.name
    source_parameter_copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_parameter, source_parameter_copy)

    dspec_prefix = results / "dspec"
    parameter_path = input_root / "data" / "params" / "paper.txt"
    parameter_path.parent.mkdir(parents=True, exist_ok=True)
    parameter_path.write_text(
        parameter_text(case, model.name, dspec_prefix, settings)
    )
    yspec_input_path = input_root / "yspec.json"
    yspec_input_path.write_text(
        json.dumps(yspec_input(case, model_copy, settings), indent=2) + "\n"
    )

    dspec_build = output / "build-dspec"
    yspec_build = output / "build-yspecpp"
    common_cmake = [
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DVERIFIED_DEPENDENCY_ROOT={dependency_root}",
    ]
    run([
        "cmake", "-S", str(BASE), "-B", str(dspec_build),
        *common_cmake,
        "-DBENCHMARK_SOLVER=dspec",
        f"-DDSPEC_SOURCE_DIR={source}",
    ])
    run(["cmake", "--build", str(dspec_build), "--target", "dspec_adapter", "--parallel", "2"])
    run([
        "cmake", "-S", str(BASE), "-B", str(yspec_build),
        *common_cmake,
        "-DBENCHMARK_SOLVER=yspecpp",
        f"-DYSPEC_SOURCE_DIR={yspec_source}",
    ])
    run(["cmake", "--build", str(yspec_build), "--target", "yspecpp_adapter", "--parallel", "2"])

    planetmodel = verify_planetmodel(yspec_build)
    mesh = case["dspec"]["mesh"]
    maxstep = case["dspec"].get("maxstep_earth_radii", 0.05)
    adapter_mode = "paper-dynamic" if mesh == "dynamic" else "paper-single"
    environment = dict(
        os.environ,
        OMP_NUM_THREADS="2",
        OMP_DYNAMIC="FALSE",
        OMP_MAX_ACTIVE_LEVELS="1",
    )
    run([
        str(dspec_build / "dspec_adapter"),
        str(parameter_path),
        str(dspec_prefix),
        str(case["dspec"]["nq"]),
        "10",
        f"{maxstep:.17g}",
        adapter_mode,
    ], log=results / "dspec.stdout.log", env=environment)
    run([
        str(yspec_build / "yspecpp_adapter"),
        str(yspec_input_path),
        str(results / "yspecpp"),
    ], log=results / "yspecpp.stdout.log", env=environment)

    provenance = write_new_run_provenance(
        output,
        source,
        source_parameter,
        source_parameter_copy,
        model,
        settings,
        case,
        yspec_tree,
        yspec_file_count,
        planetmodel,
        environment,
        adapter_mode,
    )
    snapshot_harness(output, dspec_build, yspec_build, provenance)
    (output / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    summary = process(
        args.case,
        case,
        results,
        output,
        {"provenance": "provenance.json"},
    )
    print(json.dumps({
        "output": str(output),
        "case": args.case,
        "metrics": summary["time_metrics"],
    }, indent=2))


if __name__ == "__main__":
    main()
