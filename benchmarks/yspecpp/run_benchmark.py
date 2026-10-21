#!/usr/bin/env python3
"""Build isolated solver adapters, run shared cases, and compare outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from html import escape
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

BASE = Path(__file__).resolve().parent
PINNED_DEPS = {
    "numericconcepts-src": "e03a2889d43a03660d0879f46dda1406d90c9517",
    "interpolation-src": "d5a1711ca6a89fc8ad8148fe13212398fe569bc0",
    "fftwpp-src": "cb152c0e85f3822493c2d9e85a56307ac7b4f1c0",
    "gaussquad-src": "7d4786d192f68077cd99a34dcdf1379ed1eeac9a",
    "gshtrans-src": "8dd116df8df4222fe9b10a40c7213d34f2a274fb",
}
YSPEC_COMMIT = "84be5a6cbe442324c164e5193a66b0911fca5860"
YSPEC_FILETREE_SHA256 = "5e77d93a8b2dbe469813ff3d190a664752158d7c66008a5b4c61a20ced38979e"
YSPEC_ARCHIVE_SHA256 = "6d79288770cd6b73f8679e532d8810ddc6bbc79e8fe2a06cba04e8d0ae4183c3"


def run(command: list[str], *, cwd: Path | None = None) -> None:
    print("+", " ".join(map(str, command)), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def sha(path: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def filetree_sha256(root: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    files = sorted(path for path in root.rglob("*")
                   if path.is_file() and ".git" not in path.parts)
    for path in files:
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest(), len(files)


def verify_dependencies(root: Path) -> None:
    for repo, expected in PINNED_DEPS.items():
        checkout = root / repo
        actual = sha(checkout)
        if actual != expected:
            raise SystemExit(f"dependency {repo} is {actual}, expected {expected}")
        dirty = subprocess.check_output(
            ["git", "-C", str(checkout), "status", "--porcelain",
             "--untracked-files=normal"], text=True)
        if dirty:
            raise SystemExit(f"dependency {repo} has uncommitted or untracked files")


def dspec_parameters(case: dict, defaults: dict, model_name: str, prefix: Path,
                     degree_min: int, degree_max: int,
                     mode_type: int = 4) -> str:
    source = defaults["source"]
    m = source["moment_tensor_Nm"]
    gravity = 1 if case["gravity"] == "cowling" else 2
    values = [
        str(prefix), f"models/{model_name}", str(mode_type),
        "1" if case["attenuation"] else "0", str(gravity), "1", "0",
        str(defaults["solver"]["dspec"].get("relative_error", 1e-3)),
        str(degree_min), str(degree_max),
        str(defaults["frequencies_mHz"][0]),
        str(defaults["frequencies_mHz"][1]),
        str(defaults["duration_min"]), str(defaults["time_step_s"]),
        "0", "0", "0", "0",
        str(source["depth_km"]), str(source["latitude"]),
        str(source["longitude"]), str(m["rr"]), str(m["rt"]),
        str(m["rp"]), str(m["tt"]), str(m["tp"]), str(m["pp"]),
        str(case["receiver_depth_km"]), "1",
        f"{case['receiver'][0]} {case['receiver'][1]}",
    ]
    return "\n".join(values) + "\n"


def geographic_latitude(geocentric_degrees: float) -> float:
    """Convert the shared physical latitude to YSpecpp's geographic input."""
    geocentric = math.radians(geocentric_degrees)
    return math.degrees(math.atan(math.tan(geocentric) / 0.993277))


def yspec_input(case: dict, defaults: dict, model: Path,
                degree_min: int, degree_max: int) -> dict:
    source = defaults["source"]
    return {
        "model": str(model),
        "attenuation": case["attenuation"],
        "gravity": case["gravity"],
        "output": {"prefix": "unused", "quantity": "velocity",
                   "potential_and_tilt_corrections": False},
        "degrees": {"min": degree_min, "max": degree_max},
        "frequencies_mHz": {"min": defaults["frequencies_mHz"][0],
                             "max": defaults["frequencies_mHz"][1]},
        "time_series": {"length_min": defaults["duration_min"],
                         "time_step_s": defaults["time_step_s"]},
        "source": {
            "depth_km": source["depth_km"],
            "latitude": geographic_latitude(source["latitude"]),
            "longitude": source["longitude"],
            "moment_tensor_Nm": source["moment_tensor_Nm"],
        },
        "receivers": {
            "depth_km": case["receiver_depth_km"],
            "locations": [[geographic_latitude(case["receiver"][0]),
                            case["receiver"][1]]],
        },
        "solver": defaults["solver"]["yspec"],
    }


def read_spectra(path: Path) -> tuple[dict[tuple[int, int], complex], dict]:
    values: dict[tuple[int, int], complex] = {}
    points: dict[tuple[int, int], dict] = {}
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle)
        if rows.fieldnames != ["frequency_bin", "frequency_hz",
                               "omega_real_rad_s", "omega_imag_rad_s",
                               "component", "real", "imag"]:
            raise SystemExit(f"{path}: unexpected spectrum columns")
        count = 0
        for row in rows:
            key = (int(row["frequency_bin"]), int(row["component"]))
            if key in values:
                raise SystemExit(f"{path}: duplicate frequency/component key {key}")
            scalars = [float(row[name]) for name in
                       ("frequency_hz", "omega_real_rad_s", "omega_imag_rad_s",
                        "real", "imag")]
            if not all(math.isfinite(value) for value in scalars):
                raise SystemExit(f"{path}: non-finite value at {key}")
            values[key] = complex(scalars[3], scalars[4])
            points[key] = {name: float(row[name]) for name in
                           ("frequency_hz", "omega_real_rad_s", "omega_imag_rad_s")}
            count += 1
        if count == 0:
            raise SystemExit(f"{path}: spectrum contains no samples")
    return values, points


def assert_close(name: str, left: float, right: float,
                 relative_tolerance: float = 1e-12) -> None:
    scale = max(abs(left), abs(right), 1e-30)
    if abs(left - right) > relative_tolerance * scale:
        raise SystemExit(f"{name} mismatch: {left:.17g} vs {right:.17g}")


def read_meta(path: Path) -> dict:
    return json.loads(path.read_text())


def read_timeseries(path: Path) -> tuple[list[float], list[list[float]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or reader.fieldnames[0] != "time_s" or len(reader.fieldnames) != 4:
            raise SystemExit(f"{path}: expected time_s and three components")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"{path}: time series contains no samples")
    time = [float(row["time_s"]) for row in rows]
    keys = [key for key in rows[0] if key != "time_s"]
    channels = [[float(row[key]) for row in rows] for key in keys]
    if not all(math.isfinite(value) for value in time + [v for c in channels for v in c]):
        raise SystemExit(f"{path}: time series contains a non-finite value")
    return time, channels


def compare_vectors(reference: list[complex], candidate: list[complex]) -> dict:
    if not reference or len(reference) != len(candidate):
        raise SystemExit("cannot compare empty or differently sized vectors")
    if not all(math.isfinite(abs(value)) for value in reference + candidate):
        raise SystemExit("cannot compare vectors containing non-finite values")
    difference = [a - b for a, b in zip(candidate, reference)]
    ref_norm = math.sqrt(sum(abs(value) ** 2 for value in reference))
    peak = max((abs(value) for value in reference), default=0.0)
    return {
        "relative_l2": math.sqrt(sum(abs(value) ** 2 for value in difference))
        / ref_norm if ref_norm else None,
        "max_abs_difference": max((abs(value) for value in difference),
                                    default=0.0),
        "max_difference_over_reference_peak": max(
            (abs(value) for value in difference), default=0.0
        ) / peak if peak else None,
        "reference_peak": peak,
    }


def inverse_samples(spectra: dict[tuple[int, int], complex],
                    time: list[float], meta: dict,
                    label: str = "inverse transform") -> list[list[float]]:
    """Reconstruct ZNE samples directly from the saved physical spectrum."""
    nfft = int(meta["samples"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    dt = float(meta["sample_interval_s"])
    if nfft < 2:
        raise SystemExit(f"{label}: invalid time-transform dimensions")
    channels = []
    for component in range(3):
        frequency = {
            bin_index: value for (bin_index, channel), value in spectra.items()
            if channel == component and value != 0j
        }
        positive_bins = sorted(k for k in frequency if 1 <= k <= nfft // 2)
        expected = []
        for sample, t in enumerate(time):
            n = round(t / dt)
            if abs(t - n * dt) > 1e-12 * max(abs(t), dt, 1.0):
                raise SystemExit(f"{label}: time sample is off the DFT grid")
            value = frequency.get(0, 0j).real
            last_positive = nfft // 2
            if nfft % 2 == 0:
                value += frequency.get(last_positive, 0j).real * ((-1) ** n)
                last_positive -= 1
            for bin_index in positive_bins:
                if bin_index > last_positive:
                    break
                angle = 2.0 * math.pi * bin_index * n / nfft
                value += 2.0 * (frequency[bin_index] *
                                complex(math.cos(angle), math.sin(angle))).real
            expected.append(df * math.exp(epsilon * t) * value)
        channels.append(expected)
    return channels


def validate_inverse_transform(spectra: dict[tuple[int, int], complex],
                               time: list[float], channels: list[list[float]],
                               meta: dict, label: str) -> list[dict]:
    """Check each native time trace against a direct inverse-DFT oracle."""
    if len(channels) != 3 or any(len(values) != len(time) for values in channels):
        raise SystemExit(f"{label}: invalid time-transform dimensions")
    expected_channels = inverse_samples(spectra, time, meta, label)
    checks = []
    for component, (expected, actual_values) in enumerate(
            zip(expected_channels, channels)):
        metric = compare_vectors([complex(value) for value in expected],
                                 [complex(value) for value in actual_values])
        if metric["reference_peak"] == 0.0:
            if metric["max_abs_difference"] != 0.0:
                raise SystemExit(f"{label} component {component}: zero spectrum "
                                 "produced a nonzero time series")
        elif (metric["relative_l2"] is None or
              metric["relative_l2"] > 1e-12 or
              metric["max_difference_over_reference_peak"] > 1e-12):
            raise SystemExit(f"{label} component {component}: native inverse "
                             "transform disagrees with direct DFT oracle: "
                             f"L2={metric['relative_l2']}, "
                             f"peak-scaled max="
                             f"{metric['max_difference_over_reference_peak']}")
        checks.append({"component": component, **metric})
    return checks


def compare_saved_run(case_id: str, result_dir: Path, label: str,
                      degree: int | str) -> dict:
    """Compare one sum or degree with the same spectrum, time, and DFT gates."""
    suffix = "all" if degree == "all" else label
    context = case_id if degree == "all" else f"{case_id} l={degree}"
    dspec, dpoints = read_spectra(result_dir / f"dspec-{suffix}.spectra.csv")
    yspec, ypoints = read_spectra(result_dir / f"yspecpp-{suffix}.spectra.csv")
    dmeta = read_meta(result_dir / f"dspec-{suffix}.meta.json")
    ymeta = read_meta(result_dir / f"yspecpp-{suffix}.meta.json")
    if dmeta["samples"] != ymeta["samples"]:
        raise SystemExit(f"{context}: FFT length mismatch: "
                         f"{dmeta['samples']} vs {ymeta['samples']}")
    if degree == "all" and dmeta["output_samples"] != ymeta["output_samples"]:
        raise SystemExit(f"{case_id}: output sample count mismatch")
    assert_close(f"{context}: frequency spacing",
                 dmeta["frequency_spacing_hz"], ymeta["frequency_spacing_hz"])
    if degree == "all":
        assert_close(f"{case_id}: sample_interval_s",
                     dmeta["sample_interval_s"], ymeta["sample_interval_s"])
    assert_close(f"{context}: damping epsilon",
                 dmeta["damping_epsilon_rad_s"],
                 ymeta["damping_epsilon_rad_s"])

    keys = sorted(dspec.keys() & yspec.keys())
    if len(keys) != len(dspec) or len(keys) != len(yspec):
        raise SystemExit(f"{context}: frequency/component grids do not match exactly")
    grid = None
    if degree == "all":
        grid = {
            "case": case_id,
            "nfft_dspec": dmeta["samples"],
            "nfft_yspecpp": ymeta["samples"],
            "frequency_spacing_abs_difference_hz": abs(
                dmeta["frequency_spacing_hz"] - ymeta["frequency_spacing_hz"]),
            "damping_epsilon_abs_difference_rad_s": abs(
                dmeta["damping_epsilon_rad_s"] - ymeta["damping_epsilon_rad_s"]),
            "max_frequency_abs_difference_hz": 0.0,
            "max_omega_real_abs_difference_rad_s": 0.0,
            "max_omega_imag_abs_difference_rad_s": 0.0,
        }
    for key in keys:
        for field in ("frequency_hz", "omega_real_rad_s", "omega_imag_rad_s"):
            difference = abs(dpoints[key][field] - ypoints[key][field])
            if grid is not None:
                metric_name = {
                    "frequency_hz": "max_frequency_abs_difference_hz",
                    "omega_real_rad_s": "max_omega_real_abs_difference_rad_s",
                    "omega_imag_rad_s": "max_omega_imag_abs_difference_rad_s",
                }[field]
                grid[metric_name] = max(grid[metric_name], difference)
            assert_close(f"{context} bin {key} {field}", dpoints[key][field],
                         ypoints[key][field],
                         1e-10 if field == "omega_real_rad_s" else 1e-12)

    spectra_records = []
    for component in range(3):
        component_keys = [key for key in keys if key[1] == component]
        metric = compare_vectors([yspec[key] for key in component_keys],
                                 [dspec[key] for key in component_keys])
        record = {"case": case_id}
        if degree != "all":
            record["degree"] = degree
        record.update({
            "component": component,
            "frequency_hz": [dpoints[key]["frequency_hz"]
                             for key in component_keys],
            "relative_l2": metric["relative_l2"],
            "max_abs_difference": metric["max_abs_difference"],
            "difference_over_peak": metric["max_difference_over_reference_peak"],
            "yspec_peak": metric["reference_peak"],
        })
        spectra_records.append(record)
        if degree == "all":
            for key in component_keys:
                difference = abs(dspec[key] - yspec[key])
                spectra_records.append({
                    "case": case_id, "component": component,
                    "frequency_hz": dpoints[key]["frequency_hz"],
                    "frequency_bin": key[0],
                    "point_abs_difference": difference,
                    "difference_over_peak": difference / metric["reference_peak"]
                    if metric["reference_peak"] else None,
                })

    time_d, series_d = read_timeseries(
        result_dir / f"dspec-{suffix}.timeseries.csv")
    time_y, series_y = read_timeseries(
        result_dir / f"yspecpp-{suffix}.timeseries.csv")
    inverse = [
        {"case": case_id, "degree": degree, "solver": solver,
         "components": validate_inverse_transform(
             spectra, times, series, meta, f"{case_id} {solver} {label}")}
        for solver, spectra, times, series, meta in (
            ("dspec", dspec, time_d, series_d, dmeta),
            ("yspecpp", yspec, time_y, series_y, ymeta),
        )
    ]
    if len(time_d) != len(time_y):
        raise SystemExit(f"{context}: time grid sample counts do not match")
    for index, (left, right) in enumerate(zip(time_d, time_y)):
        assert_close(f"{context} time {index}", left, right)
    series_records = []
    for component, name in enumerate(("Z", "N", "E")):
        metric = compare_vectors(
            [complex(value) for value in series_y[component]],
            [complex(value) for value in series_d[component]])
        record = {"case": case_id}
        if degree != "all":
            record["degree"] = degree
        record.update({"component": name, **metric})
        series_records.append(record)
    return {"spectra": spectra_records, "time_series": series_records,
            "inverse_transform_checks": inverse, "grid_check": grid}


def make_plots(case_id: str, outdir: Path, spectrum_metrics: list[dict]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; CSV and JSON results were still written")
        return

    time_d, dspec = read_timeseries(outdir / case_id / "dspec-all.timeseries.csv")
    time_y, yspec = read_timeseries(outdir / case_id / "yspecpp-all.timeseries.csv")
    trace_file = outdir / f"{case_id}-traces.png"
    fig, axes = plt.subplots(3, 2, figsize=(10, 7), sharex="col")
    labels = ("Z", "N", "E")
    for channel, label in enumerate(labels):
        common = min(len(time_d), len(time_y))
        residual = [dspec[channel][i] - yspec[channel][i]
                    for i in range(common)]
        axes[channel, 0].plot(time_y[:common], yspec[channel][:common],
                              label="YSpecpp", lw=1.1)
        axes[channel, 0].plot(time_d[:common], dspec[channel][:common],
                              label="DSpecM1D", lw=1.1, alpha=0.8)
        axes[channel, 1].plot(time_y[:common], residual, lw=1.0, color="black")
        axes[channel, 0].set_ylabel(f"{label} (m/s)")
        axes[channel, 1].set_ylabel("residual (m/s)")
        axes[channel, 0].grid(alpha=0.25)
        axes[channel, 1].grid(alpha=0.25)
    axes[0, 0].legend(frameon=False)
    axes[2, 0].set_xlabel("time (s)")
    axes[2, 1].set_xlabel("time (s)")
    fig.suptitle(f"{case_id}: matched velocity traces and residuals")
    fig.tight_layout()
    fig.savefig(trace_file, dpi=180)
    plt.close(fig)

    rows = [row for row in spectrum_metrics
            if row["case"] == case_id and "point_abs_difference" in row]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for component, label in enumerate(labels):
        selected = [row for row in rows if row["component"] == component]
        ax.plot([row["frequency_hz"] * 1000 for row in selected],
                [row["difference_over_peak"] for row in selected],
                marker="o", label=label)
    ax.set_xlabel("frequency (mHz)")
    ax.set_ylabel("|DSpecM1D − YSpecpp| / YSpecpp component peak")
    ax.set_title(f"{case_id}: complex spectrum error")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / f"{case_id}-spectrum-error.png", dpi=180)
    plt.close(fig)


def make_attenuation_plot(outdir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; attenuation comparison remains in CSV")
        return
    pairs = (("solid-full-off", "solid-full-on", "all-solid"),
             ("ocean-seafloor-full-off", "ocean-seafloor-full-on", "ocean"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (off_id, on_id, label) in zip(axes, pairs):
        for solver, style in (("dspec", "-"), ("yspecpp", "--")):
            off, off_points = read_spectra(outdir / off_id / f"{solver}-all.spectra.csv")
            on, on_points = read_spectra(outdir / on_id / f"{solver}-all.spectra.csv")
            for component, name in enumerate(("Z", "N", "E")):
                bins = sorted(key[0] for key in off if key[1] == component)
                ratios = []
                for bin_index in bins:
                    key = (bin_index, component)
                    before = abs(off[key])
                    ratios.append(abs(on[key]) / before if before else float("nan"))
                ax.plot([off_points[(b, component)]["frequency_hz"] * 1000 for b in bins], ratios,
                        linestyle=style, marker="o" if solver == "dspec" else "s",
                        label=f"{solver}, {name}")
        ax.axhline(1.0, color="black", lw=0.8, alpha=0.6)
        ax.set_title(label)
        ax.set_xlabel("frequency (mHz)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("|attenuation on| / |attenuation off|")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Attenuation response at matched frequencies")
    fig.tight_layout()
    fig.savefig(outdir / "attenuation-before-after.png", dpi=180)
    plt.close(fig)


def render_report(summary: dict, command: str) -> str:
    """Render a standalone report from the exact resolved benchmark output."""
    esc = lambda value: escape(str(value), quote=True)
    report_only = command.startswith("REPORT-ONLY:")
    display_command = command.removeprefix("REPORT-ONLY:").strip()

    def relative_metric(value: object, peak: float) -> str:
        if peak == 0.0:
            return "undefined (zero reference)"
        return "not applicable" if value is None else str(value)
    provenance = summary["provenance"]
    defaults = summary["defaults"]
    source = defaults["source"]
    moment = ", ".join(f"{key}={value:g} N m"
                        for key, value in source["moment_tensor_Nm"].items())
    case_ids = [case["id"] for case in summary["cases"]]
    lines = [
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>DSpecM1D and YSpecpp benchmark</title>",
        "<style>body{font:16px system-ui;max-width:1150px;margin:2rem auto;padding:0 1rem;line-height:1.5}table{border-collapse:collapse}th,td{border:1px solid #bbb;padding:.35rem .55rem;text-align:left}img{max-width:100%;height:auto}code{overflow-wrap:anywhere}</style></head><body>",
        "<h1>DSpecM1D and YSpecpp compact benchmark</h1>",
        ("<p><strong>Report-only replay:</strong> this HTML was rendered from archived JSON and plots; no solver ran to create this report. Provenance below identifies the original benchmark.</p>" if report_only else ""),
        f"<p>Run profile: <strong>{esc(provenance['refinement'])}</strong>; concurrent jobs {esc(provenance['jobs'])}; OpenMP threads {esc(provenance['runtime_environment']['OMP_NUM_THREADS'])}. All values below come from this run's resolved <code>summary.json</code>.</p>",
        f"<h2>Physical setup and solver settings</h2><p>This report covers {len(case_ids)} selected homogeneous fixture(s): {esc(case_ids)}. The available set exercises solid/ocean models, full/Cowling gravity, attenuation off/on, and ocean seafloor receivers. Fixtures intentionally avoid the known difference between YSpecpp's modulus spline and DSpecM1D's velocity spline. Source is the configured six-component moment tensor at the shared source depth/location; each receiver is specified in geocentric coordinates and converted to geographic latitude at the YSpecpp input boundary. Output is ZNE velocity in m/s; spectra are velocity per Hz.</p>",
        f"<p>Source: depth {esc(source['depth_km'])} km; geocentric latitude {esc(source['latitude'])}°, longitude {esc(source['longitude'])}°; tensor (N m): {esc(moment)}. Degrees {esc(defaults['degrees'])}; frequency band {esc(defaults['frequencies_mHz'])} mHz; duration {esc(defaults['duration_min'])} min; time step {esc(defaults['time_step_s'])} s. DSpecM1D uses nq={esc(defaults['solver']['dspec']['nq'])}, max step={esc(defaults['solver']['dspec']['maxstep'])}, relative error={esc(defaults['solver']['dspec'].get('relative_error'))}, nskip={esc(defaults['solver']['dspec']['nskip'])}; YSpecpp uses tolerance={esc(defaults['solver']['yspec']['tolerance'])}, start_decay={esc(defaults['solver']['yspec']['start_decay'])}.</p>",
        f"<p>Each fixture is run at degrees {esc(defaults['degrees'])} individually and as their sum. The comparison uses matched raw complex spectra and native unfiltered time traces. Frequency bins, complex angular frequencies, damping, and time grids are checked separately. Every native time series is also checked against an independent direct inverse-DFT calculation using the saved spectrum, physical damping correction, and sample interval.</p>",
        "<h2>Measured comparisons</h2><p>Relative L2 is the norm of (DSpecM1D − YSpecpp) divided by the YSpecpp norm. Peak-scaled maximum is the maximum absolute difference divided by the YSpecpp component peak. Limits for the fine measured profile are 5e−6 for solid full-gravity fixtures, and 1e−5 for ocean full-gravity and Cowling fixtures; the same limits cover attenuation-off cases and their attenuation-on full-gravity regressions. Coarse and medium profiles are diagnostic and are not acceptance measurements. When the reference peak is zero, normalized metrics are shown as undefined and the absolute error remains visible with its physical units. Such a channel passes only when its maximum absolute error is exactly zero; even a subnormal nonzero difference fails. This report displays measurements and limits; cross-code agreement is not a scientific accuracy claim.</p>",
        "<p>Trace plots overlay native Z/N/E velocity traces and show DSpecM1D-minus-YSpecpp residuals in m/s. Spectrum plots show pointwise complex-spectrum magnitude difference divided by the YSpecpp component peak. When present, the attenuation plot shows on/off amplitude ratios by solver and component.</p>",
        "<table><thead><tr><th>Case/model</th><th>Receiver (geocentric lat, lon; depth km)</th><th>Gravity</th><th>Attenuation</th><th>Component/degree</th><th>FFT/output samples</th><th>df (Hz)</th><th>Damping epsilon (rad/s)</th><th>Spectrum relative L2</th><th>Spectrum max / peak</th><th>Spectrum max absolute (m/s/Hz)</th><th>Time relative L2</th><th>Time max / peak</th><th>Time max absolute (m/s)</th><th>Fine gate</th></tr></thead><tbody>",
    ]
    for case in summary["cases"]:
        case_id = case["id"]
        spectra = {(row["component"], row.get("degree", "all")): row
                   for row in summary["degree_spectra"] + summary["spectra"]
                   if row["case"] == case_id and "relative_l2" in row}
        traces = {(row["component"], row.get("degree", "all")): row
                  for row in summary["degree_time_series"] + summary["time_series"]
                  if row["case"] == case_id}
        grid = summary["run_metadata"][case_id]["dspec"]
        for key in sorted(spectra, key=lambda item: (str(item[1]), item[0])):
            component = ("Z", "N", "E")[key[0]]
            spec, trace = spectra[key], traces[(component, key[1])]
            degree = "sum" if key[1] == "all" else f"l={key[1]}"
            limit = (1e-5 if case["ocean"] or case["gravity"] == "cowling"
                     else 5e-6)
            component_checks = (
                (spec["yspec_peak"], spec["max_abs_difference"],
                 spec["relative_l2"], spec["difference_over_peak"]),
                (trace["reference_peak"], trace["max_abs_difference"],
                 trace["relative_l2"], trace["max_difference_over_reference_peak"]),
            )
            passed = all(
                maximum == 0.0 if peak == 0.0 else
                l2 is not None and peak_scaled is not None
                and l2 <= limit and peak_scaled <= limit
                for peak, maximum, l2, peak_scaled in component_checks
            )
            gate = ("not measured" if provenance["refinement"] != "fine" else
                    "PASS" if passed else "FAIL")
            lines.append("<tr><td>" + "</td><td>".join(map(esc, (
                f"{case_id} / {case['model']}",
                f"{case['receiver'][0]}°, {case['receiver'][1]}°; {case['receiver_depth_km']}",
                case["gravity"],
                "on" if case["attenuation"] else "off", f"{degree}, {component}",
                f"{grid['samples']}/{grid['output_samples']}",
                grid["frequency_spacing_hz"], grid["damping_epsilon_rad_s"],
                relative_metric(spec["relative_l2"], spec["yspec_peak"]),
                relative_metric(spec["difference_over_peak"], spec["yspec_peak"]),
                spec["max_abs_difference"],
                relative_metric(trace["relative_l2"], trace["reference_peak"]),
                relative_metric(trace["max_difference_over_reference_peak"],
                                trace["reference_peak"]),
                trace["max_abs_difference"],
                f"{gate} (limit {limit:g})",
            ))) + "</td></tr>")
    lines += ["</tbody></table><h2>Reproduction and provenance</h2>",
        f"<p>Command: <code>{esc(display_command)}</code>. YSpecpp commit {esc(provenance['yspecpp_commit'])}; file tree SHA-256 {esc(provenance['yspecpp_filetree_sha256'])}; archive SHA-256 {esc(provenance['yspecpp_source_archive_sha256'])}. DSpecM1D HEAD {esc(provenance['dspec_source_head'])}; its tracked-code diff SHA-256 is {esc(provenance['dspec_tracked_source_diff_sha256'])}. The case specification SHA-256 is {esc(provenance['case_spec_sha256'])}; fixture identities are in <a href=\"provenance.json\">provenance.json</a>.</p>",
        f"<p>Dependency revisions: <code>{esc(provenance['dependency_shas'])}</code>. YSpecpp's separate planetmodel pin: <code>{esc(provenance['planetmodel_sha'])}</code>. The runner checks these identities before solving.</p>",
        "<p>Limits: compact homogeneous fixtures only; no paper-scale runs, fluid-source support, convergence campaign, or scientific accuracy claim. The solver projects and spline choices remain those of their independently pinned projects.</p><h2>Plots</h2><ul>"]
    for case in summary["cases"]:
        case_id = case["id"]
        for suffix, label in (("traces", "time traces and residuals"),
                              ("spectrum-error", "pointwise spectrum error")):
            name = f"{case_id}-{suffix}.png"
            lines.append(f"<li><a href=\"{esc(name)}\">{esc(case_id)} {label}</a></li><figure><img src=\"{esc(name)}\" alt=\"{esc(case_id)} {label}\"></figure>")
    available = {case["id"] for case in summary["cases"]}
    if {"solid-full-off", "solid-full-on", "ocean-seafloor-full-off",
            "ocean-seafloor-full-on"}.issubset(available):
        lines.append("<li><a href=\"attenuation-before-after.png\">attenuation response</a></li><figure><img src=\"attenuation-before-after.png\" alt=\"Attenuation response ratios\"></figure>")
    lines += ["</ul><p>Machine-readable results: <a href=\"summary.json\">summary.json</a>, <a href=\"summary.csv\">summary.csv</a>, and <a href=\"provenance.json\">provenance.json</a>.</p></body></html>\n"]
    return "".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dspec-source", type=Path, required=True)
    parser.add_argument("--yspec-source", type=Path, required=True)
    parser.add_argument("--dependency-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--case", action="append", dest="case_ids")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--refinement", choices=("coarse", "medium", "fine"),
                        default="coarse")
    parser.add_argument("--configure-only", action="store_true")
    args = parser.parse_args()
    outdir = args.output or args.dspec_source / "build/phase3/benchmark-results"
    outdir = outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"output directory must be empty: {outdir}")
    if not 1 <= args.jobs <= 2:
        parser.error("--jobs must be 1 or 2")
    os.environ["OMP_NUM_THREADS"] = "2"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    if not args.yspec_source.is_dir():
        parser.error("--yspec-source must point to the immutable YSpecpp snapshot")
    if not args.dspec_source.is_dir():
        parser.error("--dspec-source must point to the DSpecM1D checkout")
    verify_dependencies(args.dependency_root)
    yspec_tree_hash, yspec_file_count = filetree_sha256(args.yspec_source)
    if yspec_tree_hash != YSPEC_FILETREE_SHA256 or yspec_file_count != 76:
        raise SystemExit("YSpecpp snapshot does not match the pinned source archive")

    inputs = outdir / "inputs"
    models = inputs / "data/models"
    params_dir = inputs / "data/params"
    yjson_dir = inputs / "yspec"
    build_d = outdir / "build-dspec"
    build_y = outdir / "build-yspecpp"
    for path in (models, params_dir, yjson_dir, outdir):
        path.mkdir(parents=True, exist_ok=True)

    for solver, build_dir, extra in (
        ("dspec", build_d, [f"-DDSPEC_SOURCE_DIR={args.dspec_source.resolve()}"]),
        ("yspecpp", build_y, [f"-DYSPEC_SOURCE_DIR={args.yspec_source.resolve()}"]),
    ):
        run(["cmake", "-S", str(BASE), "-B", str(build_dir),
             f"-DBENCHMARK_SOLVER={solver}",
             f"-DVERIFIED_DEPENDENCY_ROOT={args.dependency_root.resolve()}",
             *extra])
    if args.configure_only:
        return

    expected_ys = args.yspec_source / "YSpec"
    if not expected_ys.is_dir():
        raise SystemExit(f"not a YSpecpp source tree: {args.yspec_source}")
    for build_dir, target in ((build_d, "dspec_adapter"),
                              (build_y, "yspecpp_adapter")):
        run(["cmake", "--build", str(build_dir), "--target", target,
             "--parallel", str(args.jobs)])
    fetched_planetmodel = build_y / "_deps/planetmodel-src"
    if not fetched_planetmodel.is_dir() or not (fetched_planetmodel / ".git").exists():
        raise SystemExit("YSpecpp did not create its pinned planetmodel checkout")
    actual = sha(fetched_planetmodel)
    expected = "fa6e230fde25152d40f27d9839de44c4f3182a24"
    if actual != expected:
        raise SystemExit(f"planetmodel is {actual}, expected {expected}")
    if subprocess.check_output(
            ["git", "-C", str(fetched_planetmodel), "status", "--porcelain",
             "--untracked-files=normal"], text=True):
        raise SystemExit("planetmodel has uncommitted or untracked files")

    with (BASE / "cases/cases.json").open() as handle:
        specification = json.load(handle)
    defaults = specification["defaults"]
    defaults = json.loads(json.dumps(defaults))
    if args.refinement == "medium":
        defaults["solver"]["dspec"].update({"nq": 5, "maxstep": 0.0125})
        defaults["solver"]["yspec"]["tolerance"] = 1e-7
    elif args.refinement == "fine":
        defaults["solver"]["dspec"].update({"nq": 7, "maxstep": 0.00625})
        defaults["solver"]["yspec"]["tolerance"] = 1e-9
    cases = specification["cases"]
    if args.case_ids:
        wanted = set(args.case_ids)
        cases = [case for case in cases if case["id"] in wanted]
        missing = wanted - {case["id"] for case in cases}
        if missing:
            raise SystemExit(f"unknown case id(s): {', '.join(sorted(missing))}")

    provenance = {
        "case_spec": "benchmarks/yspecpp/cases/cases.json",
        "case_spec_sha256": hashlib.sha256(
            (BASE / "cases/cases.json").read_bytes()).hexdigest(),
        "fixture_sha256": {
            name: hashlib.sha256((BASE / "cases" / name).read_bytes()).hexdigest()
            for name in sorted({case["model"] for case in cases})
        },
        "dspec_source_head": subprocess.check_output(
            ["git", "-C", str(args.dspec_source), "rev-parse", "HEAD"], text=True
        ).strip(),
        "command": shlex.join(
            ["python3", str(Path(__file__).resolve()), *sys.argv[1:]]
        ),
        "yspecpp_commit": YSPEC_COMMIT,
        "yspecpp_filetree_sha256": yspec_tree_hash,
        "yspecpp_file_count": yspec_file_count,
        "yspecpp_source_archive_sha256": YSPEC_ARCHIVE_SHA256,
        "dependency_shas": PINNED_DEPS,
        "planetmodel_sha": "fa6e230fde25152d40f27d9839de44c4f3182a24",
        "jobs": args.jobs,
        "refinement": args.refinement,
        "runtime_environment": {
            "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"],
            "OMP_DYNAMIC": os.environ["OMP_DYNAMIC"],
        },
        "dspec_tracked_source_files_changed": subprocess.check_output(
            ["git", "-C", str(args.dspec_source), "diff", "--name-only",
             "HEAD", "--", "DSpecM1D"], text=True).splitlines(),
        "dspec_tracked_source_diff_sha256": hashlib.sha256(subprocess.check_output(
            ["git", "-C", str(args.dspec_source), "diff", "--binary",
             "HEAD", "--", "DSpecM1D"])).hexdigest(),
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    spectrum_metrics: list[dict] = []
    series_metrics: list[dict] = []
    degree_spectrum_metrics: list[dict] = []
    degree_series_metrics: list[dict] = []
    grid_checks: list[dict] = []
    inverse_transform_checks: list[dict] = []
    for case in cases:
        case_id = case["id"]
        result_dir = outdir / case_id
        result_dir.mkdir(parents=True, exist_ok=True)
        model_source = BASE / "cases" / case["model"]
        model_target = models / case["model"]
        shutil.copyfile(model_source, model_target)
        dspec_exe = build_d / "dspec_adapter"
        yspec_exe = build_y / "yspecpp_adapter"

        degree_runs = [(0, defaults["degrees"][-1], 0,
                        defaults["degrees"][-1], "all", 4)]
        degree_runs.extend(
            (0, 1, 0, 0, "l0", 1) if degree == 0 else
            (degree, degree, degree, degree, f"l{degree}", 4)
            for degree in defaults["degrees"])
        for d_min, d_max, y_min, y_max, label, mode_type in degree_runs:
            prefix_d = result_dir / f"dspec-{label}"
            dspec_settings = defaults["solver"]["dspec"]
            params_path = params_dir / f"{case_id}-{label}.txt"
            params_path.write_text(dspec_parameters(
                case, defaults, model_target.name, prefix_d,
                d_min, d_max, mode_type))
            yinput_path = yjson_dir / f"{case_id}-{label}.json"
            yinput_path.write_text(json.dumps(yspec_input(
                case, defaults, model_target, y_min, y_max), indent=2) + "\n")
            prefix_y = result_dir / f"yspecpp-{label}"
            run([str(dspec_exe), str(params_path), str(prefix_d),
                 str(dspec_settings["nq"]), str(dspec_settings["nskip"]),
                 str(dspec_settings["maxstep"])])
            run([str(yspec_exe), str(yinput_path), str(prefix_y)])

        for label, degree in [("all", "all"),
                              *((f"l{value}", value)
                                for value in defaults["degrees"])]:
            compared = compare_saved_run(case_id, result_dir, label, degree)
            if degree == "all":
                spectrum_metrics.extend(compared["spectra"])
                series_metrics.extend(compared["time_series"])
                grid_checks.append(compared["grid_check"])
            else:
                degree_spectrum_metrics.extend(compared["spectra"])
                degree_series_metrics.extend(compared["time_series"])
            inverse_transform_checks.extend(
                compared["inverse_transform_checks"])

    run_metadata = {}
    for case in cases:
        run_metadata[case["id"]] = {
            solver: read_meta(outdir / case["id"] / f"{solver}-all.meta.json")
            for solver in ("dspec", "yspecpp")
        }
    (outdir / "summary.json").write_text(json.dumps({
        "provenance": provenance,
        "defaults": defaults,
        "cases": cases,
        "model_decks": {
            name: (BASE / "cases" / name).read_text()
            for name in sorted({case["model"] for case in cases})
        },
        "run_metadata": run_metadata,
        "spectra": spectrum_metrics,
        "time_series": series_metrics,
        "degree_spectra": degree_spectrum_metrics,
        "degree_time_series": degree_series_metrics,
        "grid_checks": grid_checks,
        "inverse_transform_checks": inverse_transform_checks,
    }, indent=2) + "\n")
    with (outdir / "summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["kind", "case", "component", "frequency_hz",
                         "relative_l2", "max_abs_difference",
                         "max_difference_over_yspec_peak"])
        for item in spectrum_metrics:
            if "relative_l2" in item:
                writer.writerow(["spectrum", item["case"], item["component"], "",
                                 item["relative_l2"], item["max_abs_difference"],
                                 item["difference_over_peak"]])
            else:
                writer.writerow(["spectrum_point", item["case"],
                                 item["component"], item["frequency_hz"], "",
                                 item["point_abs_difference"],
                                 item["difference_over_peak"]])
        for item in series_metrics:
            writer.writerow(["time_series", item["case"], item["component"], "",
                             item["relative_l2"], item["max_abs_difference"],
                             item["max_difference_over_reference_peak"]])
        for item in degree_spectrum_metrics:
            writer.writerow(["degree_spectrum", item["case"],
                             f"l{item['degree']}:component{item['component']}", "",
                             item["relative_l2"], item["max_abs_difference"],
                             item["difference_over_peak"]])
        for item in degree_series_metrics:
            writer.writerow(["degree_time_series", item["case"],
                             f"l{item['degree']}:{item['component']}", "",
                             item["relative_l2"], item["max_abs_difference"],
                             item["max_difference_over_reference_peak"]])

    for case in cases:
        make_plots(case["id"], outdir, spectrum_metrics)
    available = {case["id"] for case in cases}
    if {"solid-full-off", "solid-full-on", "ocean-seafloor-full-off",
            "ocean-seafloor-full-on"}.issubset(available):
        make_attenuation_plot(outdir)
    command = provenance["command"]
    (outdir / "report.html").write_text(
        render_report(json.loads((outdir / "summary.json").read_text()), command)
    )
    print(f"Results: {outdir}")


if __name__ == "__main__":
    main()
