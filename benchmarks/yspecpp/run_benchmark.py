#!/usr/bin/env python3
"""Build isolated solver adapters, run shared cases, and compare outputs."""

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


def validate_inverse_transform(spectra: dict[tuple[int, int], complex],
                               time: list[float], channels: list[list[float]],
                               meta: dict, label: str) -> list[dict]:
    """Check each native time trace against a direct inverse-DFT oracle."""
    nfft = int(meta["samples"])
    df = float(meta["frequency_spacing_hz"])
    epsilon = float(meta["damping_epsilon_rad_s"])
    dt = float(meta["sample_interval_s"])
    if nfft < 2 or len(time) != len(channels[0]):
        raise SystemExit(f"{label}: invalid time-transform dimensions")
    checks = []
    for component, actual_values in enumerate(channels):
        frequency = [spectra.get((bin_index, component), 0j)
                     for bin_index in range(nfft // 2 + 1)]
        expected = []
        for sample, t in enumerate(time):
            n = round(t / dt)
            if abs(t - n * dt) > 1e-12 * max(abs(t), dt, 1.0):
                raise SystemExit(f"{label}: time sample is off the DFT grid")
            value = frequency[0].real
            last_positive = nfft // 2
            if nfft % 2 == 0:
                value += frequency[last_positive].real * ((-1) ** n)
                last_positive -= 1
            for bin_index in range(1, last_positive + 1):
                angle = 2.0 * math.pi * bin_index * n / nfft
                value += 2.0 * (frequency[bin_index] *
                                complex(math.cos(angle), math.sin(angle))).real
            expected.append(df * math.exp(epsilon * t) * value)
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

        dspec_full, dspec_points = read_spectra(result_dir / "dspec-all.spectra.csv")
        yspec_full, yspec_points = read_spectra(result_dir / "yspecpp-all.spectra.csv")
        meta_d = read_meta(result_dir / "dspec-all.meta.json")
        meta_y = read_meta(result_dir / "yspecpp-all.meta.json")
        for field in ("samples", "frequency_spacing_hz"):
            if field == "samples":
                if meta_d[field] != meta_y[field]:
                    raise SystemExit(f"{case_id}: FFT length mismatch: {meta_d[field]} vs {meta_y[field]}")
            else:
                assert_close(f"{case_id}: {field}", meta_d[field], meta_y[field])
        for field in ("output_samples", "sample_interval_s"):
            if field == "output_samples":
                if meta_d[field] != meta_y[field]:
                    raise SystemExit(f"{case_id}: output sample count mismatch")
            else:
                assert_close(f"{case_id}: {field}", meta_d[field], meta_y[field])
        assert_close(f"{case_id}: damping epsilon",
                     meta_d["damping_epsilon_rad_s"],
                     meta_y["damping_epsilon_rad_s"])
        common_keys = sorted(dspec_full.keys() & yspec_full.keys())
        if len(common_keys) != len(dspec_full) or len(common_keys) != len(yspec_full):
            raise SystemExit(f"{case_id}: frequency/component grids do not match exactly")
        grid_check = {
            "case": case_id,
            "nfft_dspec": meta_d["samples"],
            "nfft_yspecpp": meta_y["samples"],
            "frequency_spacing_abs_difference_hz": abs(
                meta_d["frequency_spacing_hz"] - meta_y["frequency_spacing_hz"]),
            "damping_epsilon_abs_difference_rad_s": abs(
                meta_d["damping_epsilon_rad_s"] - meta_y["damping_epsilon_rad_s"]),
            "max_frequency_abs_difference_hz": 0.0,
            "max_omega_real_abs_difference_rad_s": 0.0,
            "max_omega_imag_abs_difference_rad_s": 0.0,
        }
        for key in common_keys:
            for field in ("frequency_hz", "omega_real_rad_s", "omega_imag_rad_s"):
                difference = abs(dspec_points[key][field] - yspec_points[key][field])
                metric_name = {
                    "frequency_hz": "max_frequency_abs_difference_hz",
                    "omega_real_rad_s": "max_omega_real_abs_difference_rad_s",
                    "omega_imag_rad_s": "max_omega_imag_abs_difference_rad_s",
                }[field]
                grid_check[metric_name] = max(grid_check[metric_name], difference)
                assert_close(f"{case_id} bin {key} {field}",
                             dspec_points[key][field], yspec_points[key][field],
                             1e-10 if field == "omega_real_rad_s" else 1e-12)
        grid_checks.append(grid_check)
        for component in range(3):
            keys = [key for key in common_keys if key[1] == component]
            metric = compare_vectors([yspec_full[key] for key in keys],
                                     [dspec_full[key] for key in keys])
            spectrum_metrics.append({
                "case": case_id, "component": component,
                "frequency_hz": [dspec_points[key]["frequency_hz"] for key in keys],
                "relative_l2": metric["relative_l2"],
                "max_abs_difference": metric["max_abs_difference"],
                "difference_over_peak": metric["max_difference_over_reference_peak"],
                "yspec_peak": metric["reference_peak"],
            })
            for key in keys:
                difference = abs(dspec_full[key] - yspec_full[key])
                spectrum_metrics.append({
                    "case": case_id, "component": component,
                    "frequency_hz": dspec_points[key]["frequency_hz"],
                    "frequency_bin": key[0], "point_abs_difference": difference,
                    "difference_over_peak": difference / metric["reference_peak"]
                    if metric["reference_peak"] else None,
                })

        time_d, dseries = read_timeseries(result_dir / "dspec-all.timeseries.csv")
        time_y, yseries = read_timeseries(result_dir / "yspecpp-all.timeseries.csv")
        inverse_transform_checks.append({
            "case": case_id, "degree": "all", "solver": "dspec",
            "components": validate_inverse_transform(
                dspec_full, time_d, dseries, meta_d, f"{case_id} DSpec all")})
        inverse_transform_checks.append({
            "case": case_id, "degree": "all", "solver": "yspecpp",
            "components": validate_inverse_transform(
                yspec_full, time_y, yseries, meta_y, f"{case_id} YSpecpp all")})
        if len(time_d) != len(time_y):
            raise SystemExit(f"{case_id}: time grid sample counts do not match")
        for index, (time_left, time_right) in enumerate(zip(time_d, time_y)):
            assert_close(f"{case_id} time sample {index}", time_left, time_right)
        for component, label in enumerate(("Z", "N", "E")):
            metric = compare_vectors(
                [complex(value) for value in yseries[component]],
                [complex(value) for value in dseries[component]],
            )
            series_metrics.append({"case": case_id, "component": label,
                                   **metric})

        for degree in defaults["degrees"]:
            label = f"l{degree}"
            dspec_degree, dspec_degree_points = read_spectra(
                result_dir / f"dspec-{label}.spectra.csv")
            yspec_degree, yspec_degree_points = read_spectra(
                result_dir / f"yspecpp-{label}.spectra.csv")
            meta_d_degree = read_meta(result_dir / f"dspec-{label}.meta.json")
            meta_y_degree = read_meta(result_dir / f"yspecpp-{label}.meta.json")
            if meta_d_degree["samples"] != meta_y_degree["samples"]:
                raise SystemExit(f"{case_id} l={degree}: FFT length mismatch")
            assert_close(f"{case_id} l={degree}: frequency spacing",
                         meta_d_degree["frequency_spacing_hz"],
                         meta_y_degree["frequency_spacing_hz"])
            assert_close(f"{case_id} l={degree}: damping epsilon",
                         meta_d_degree["damping_epsilon_rad_s"],
                         meta_y_degree["damping_epsilon_rad_s"])
            degree_keys = sorted(dspec_degree.keys() & yspec_degree.keys())
            if len(degree_keys) != len(dspec_degree) or len(degree_keys) != len(yspec_degree):
                raise SystemExit(f"{case_id} l={degree}: Fourier bins differ")
            for key in degree_keys:
                for field in ("frequency_hz", "omega_real_rad_s", "omega_imag_rad_s"):
                    assert_close(f"{case_id} l={degree} bin {key} {field}",
                                 dspec_degree_points[key][field],
                                 yspec_degree_points[key][field],
                                 1e-10 if field == "omega_real_rad_s" else 1e-12)
            for component in range(3):
                keys = [key for key in degree_keys if key[1] == component]
                metric = compare_vectors(
                    [yspec_degree[key] for key in keys],
                    [dspec_degree[key] for key in keys])
                degree_spectrum_metrics.append({
                    "case": case_id, "degree": degree, "component": component,
                    "frequency_hz": [dspec_degree_points[key]["frequency_hz"]
                                     for key in keys],
                    "relative_l2": metric["relative_l2"],
                    "max_abs_difference": metric["max_abs_difference"],
                    "difference_over_peak":
                        metric["max_difference_over_reference_peak"],
                    "yspec_peak": metric["reference_peak"],
                })
            time_d_degree, dseries_degree = read_timeseries(
                result_dir / f"dspec-{label}.timeseries.csv")
            time_y_degree, yseries_degree = read_timeseries(
                result_dir / f"yspecpp-{label}.timeseries.csv")
            inverse_transform_checks.append({
                "case": case_id, "degree": degree, "solver": "dspec",
                "components": validate_inverse_transform(
                    dspec_degree, time_d_degree, dseries_degree, meta_d_degree,
                    f"{case_id} DSpec l={degree}")})
            inverse_transform_checks.append({
                "case": case_id, "degree": degree, "solver": "yspecpp",
                "components": validate_inverse_transform(
                    yspec_degree, time_y_degree, yseries_degree, meta_y_degree,
                    f"{case_id} YSpecpp l={degree}")})
            if len(time_d_degree) != len(time_y_degree):
                raise SystemExit(f"{case_id} l={degree}: time sample counts differ")
            for index, (left, right) in enumerate(zip(time_d_degree, time_y_degree)):
                assert_close(f"{case_id} l={degree} time {index}", left, right)
            for component, component_name in enumerate(("Z", "N", "E")):
                metric = compare_vectors(
                    [complex(value) for value in yseries_degree[component]],
                    [complex(value) for value in dseries_degree[component]])
                degree_series_metrics.append({
                    "case": case_id, "degree": degree,
                    "component": component_name, **metric,
                })

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
    print(f"Results: {outdir}")


if __name__ == "__main__":
    main()
