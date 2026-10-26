#!/usr/bin/env python3
"""Run and report the bounded DSpecM1D/YSpecpp elastic-kernel comparison."""
from __future__ import annotations

import argparse
import base64
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

HERE = Path(__file__).resolve().parent
DSPEC_DEFAULT = Path(__file__).resolve().parents[3]
YSPEC_DEFAULT = DSPEC_DEFAULT.parent / "YSpecpp"
DEPS_DSPEC = DSPEC_DEFAULT / "build/phase4-20261021/_deps"
DEPS_YSPEC = DSPEC_DEFAULT / "build/phase12-fix/fetched-pinned/_deps"
PLANETMODEL = YSPEC_DEFAULT / "build/library/_deps/planetmodel-src"
PLANETMODEL_DSPEC = DEPS_DSPEC / "planetmodel-src"
YSPEC_COMMIT = "a294fb321cbcf9b6736710eb7ebdcec2c694dc5c"
PLANETMODEL_COMMITS = {
    "dspec": "35be4ee54b3d84df1483e647878c5446a886c90b",
    "yspecpp": "fa6e230fde25152d40f27d9839de44c4f3182a24",
}
DEP_PINS = {
    "numericconcepts-src": "e03a2889d43a03660d0879f46dda1406d90c9517",
    "interpolation-src": "d5a1711ca6a89fc8ad8148fe13212398fe569bc0",
    "fftwpp-src": "cb152c0e85f3822493c2d9e85a56307ac7b4f1c0",
    "gaussquad-src": "7d4786d192f68077cd99a34dcdf1379ed1eeac9a",
    "gshtrans-src": "8dd116df8df4222fe9b10a40c7213d34f2a274fb",
}
GATE = {"metadata_rel": 1e-12, "omega_rel": 5e-11, "spectra_l2": 1e-5,
        "spectra_peak_max": 1e-5, "objective_rel": 1e-5,
        "derivative_rel": 1e-3, "kernel_fd_rel": 1e-4,
        "source_local_kernel_fd_rel": 1e-6,
        "informative_floor_fraction": 1e-6}


def command(args: list[str], env: dict[str, str] | None = None) -> None:
    print("+", shlex.join(args), flush=True)
    subprocess.run(args, check=True, env=env)


def git(root: Path, *args: str) -> str:
    value = subprocess.check_output(["git", "-C", str(root), *args], text=True)
    return value.rstrip("\n") if args and args[0] == "status" else value.strip()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_identity(root: Path) -> dict:
    status = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain",
                                      "--untracked-files=all"], text=True).rstrip("\n")
    dirty_files = {}
    for line in status.splitlines():
        path = root / line[3:]
        if path.is_file():
            dirty_files[line[3:]] = file_hash(path)
    return {"path": str(root), "commit": git(root, "rev-parse", "HEAD"),
            "status_porcelain": status, "dirty_file_sha256": dirty_files}


def clean_dependency(path: Path, expected_commit: str) -> dict:
    commit = git(path, "rev-parse", "HEAD")
    status = git(path, "status", "--porcelain", "--untracked-files=all")
    if commit != expected_commit or status:
        raise SystemExit(f"dependency is not the expected clean checkout: {path} ({commit}, {status!r})")
    return {"path": str(path), "commit": commit, "status_porcelain": status}


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as f:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def rel(a: float, b: float) -> float:
    return abs(a-b) / max(abs(a), abs(b), 1e-300)


def spectrum(path: Path) -> dict[tuple[str, int], complex]:
    rows = read_csv(path)
    result = {}
    for r in rows:
        key = (r["component"], int(r["bin"]))
        if key in result:
            raise ValueError(f"duplicate spectrum row: {key}")
        value = complex(float(r["real_m_s"]), float(r["imag_m_s"]))
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValueError(f"non-finite spectrum row: {key}")
        result[key] = value
    expected = {(c,b) for c in "ZNE" for b in (1,2,3)}
    if set(result) != expected:
        raise ValueError(f"unexpected spectrum grid in {path}: {set(result)}")
    return result


def derivatives(path: Path) -> dict[tuple[str, str], dict]:
    result = {}
    for r in read_csv(path):
        key=(r["observable"],r["modulus"])
        if key in result:
            raise ValueError(f"duplicate derivative row: {key}")
        for field in ("objective","kernel","volume_kernel","fd_h_0_01","fd_h_0_005","fd_h_0_0025",
                      "richardson","relative_kernel_fd","source_kernel","source_delta",
                      "source_contribution","total_kernel"):
            if not math.isfinite(float(r[field])):
                raise ValueError(f"non-finite derivative value {key}/{field}")
        result[key]=r
    if not result:
        raise ValueError(f"no derivative rows in {path}")
    return result


def metrics(left: dict, right: dict) -> tuple[dict, list[dict]]:
    if set(left) != set(right):
        raise ValueError("spectrum keys differ")
    by_component = {}
    for component in "ZNE":
        keys = sorted(k for k in left if k[0] == component)
        numerator = sum(abs(left[k]-right[k])**2 for k in keys)
        denom = sum(abs(right[k])**2 for k in keys)
        peak = max((abs(left[k]-right[k]) for k in keys), default=0)
        scale = max((abs(right[k]) for k in keys), default=0)
        by_component[component] = {
            "l2": math.sqrt(numerator / max(denom, 1e-300)),
            "peak_scaled_max": peak / max(scale, 1e-300),
            "pass": math.sqrt(numerator / max(denom, 1e-300)) <= GATE["spectra_l2"]
                    and peak / max(scale, 1e-300) <= GATE["spectra_peak_max"]}
    detail = [{"component": k[0], "bin": k[1], "left_real": left[k].real,
               "left_imag": left[k].imag, "right_real": right[k].real,
               "right_imag": right[k].imag,
               "complex_relative_error": abs(left[k]-right[k]) /
                   max(abs(left[k]), abs(right[k]), 1e-300)}
              for k in sorted(left)]
    return by_component, detail


def check_metadata(left: dict, right: dict) -> tuple[list[dict], float]:
    exact_fields = ("nt", "dt_s", "df_hz", "damping_s", "gravity", "attenuation",
                    "degrees", "q_bins", "q_values", "window_samples", "window",
                    "source_geocentric", "receiver_geocentric", "moment_tensor_Nm",
                    "shell_A_C_F_L_N_Pa", "physical_G", "radius_m",
                    "density_kg_m3", "vp_m_s", "vs_m_s", "source_radius_m",
                    "support_layer", "support_bounds_m")
    checks = []
    for key in exact_fields:
        x, y = left[key], right[key]
        if isinstance(x, list):
            if len(x) != len(y):
                error = math.inf
            else:
                error = max((rel(float(a), float(b)) if isinstance(a, (int, float)) else
                             (0.0 if a == b else math.inf) for a, b in zip(x, y)), default=0.0)
        elif isinstance(x, (int, float)):
            error = rel(float(x), float(y))
        else:
            error = 0.0 if x == y else math.inf
        ok = error <= GATE["metadata_rel"]
        checks.append({"quantity": key, "left": x, "right": y,
                       "relative_or_exact_error": error, "pass": ok})
    omega = rel(float(left["omega_bin1"]), float(right["omega_bin1"]))
    checks.append({"quantity": "omega_bin1_rad_s", "left": left["omega_bin1"],
                   "right": right["omega_bin1"], "relative_or_exact_error": omega,
                   "pass": omega <= GATE["omega_rel"]})
    return checks, omega


def compare_pair(label: str, left_dir: Path, right_dir: Path,
                 kernel_fd_limit: float | None = None,
                 reference_scales: dict[str, tuple[float, float]] | None = None) -> tuple[dict, list[dict]]:
    left_meta = json.loads((left_dir / "meta.json").read_text())
    right_meta = json.loads((right_dir / "meta.json").read_text())
    metadata, omega = check_metadata(left_meta, right_meta)
    pair_pass = all(item["pass"] for item in metadata)
    spectra_metrics, spectrum_rows = metrics(
        spectrum(left_dir / "spectra.csv"), spectrum(right_dir / "spectra.csv"))
    spectra_pass = all(item["pass"] for item in spectra_metrics.values())
    pair_pass &= spectra_pass
    left_derivatives = derivatives(left_dir / "derivatives.csv")
    right_derivatives = derivatives(right_dir / "derivatives.csv")
    if set(left_derivatives) != set(right_derivatives):
        raise ValueError(f"derivative keys differ for {label}")
    merged, rows = [], []
    for key in sorted(left_derivatives):
        x, y = left_derivatives[key], right_derivatives[key]
        left_kernel, right_kernel = float(x["kernel"]), float(y["kernel"])
        observable = key[0]
        scale = max(
            max(abs(float(v["kernel"])) for v in left_derivatives.values()
                if v["observable"] == observable),
            max(abs(float(v["kernel"])) for v in right_derivatives.values()
                if v["observable"] == observable), 1e-300)
        if reference_scales is not None:
            scale = reference_scales[observable][0]
        informative = max(abs(left_kernel), abs(right_kernel)) >= GATE["informative_floor_fraction"] * scale
        derivative_relative = rel(left_kernel, right_kernel)
        derivative_scaled = abs(left_kernel - right_kernel) / scale
        derivative_pass = (derivative_relative <= GATE["derivative_rel"] if informative
                           else derivative_scaled <= GATE["derivative_rel"])

        left_objective, right_objective = float(x["objective"]), float(y["objective"])
        objective_scale = max(
            max(abs(float(v["objective"])) for v in left_derivatives.values()
                if v["observable"] == observable),
            max(abs(float(v["objective"])) for v in right_derivatives.values()
                if v["observable"] == observable), 1e-300)
        if reference_scales is not None:
            objective_scale = reference_scales[observable][1]
        objective_relative = rel(left_objective, right_objective)
        objective_informative = max(abs(left_objective), abs(right_objective)) >= (
            GATE["informative_floor_fraction"] * objective_scale)
        objective_scaled = abs(left_objective - right_objective) / objective_scale
        objective_pass = (objective_relative <= GATE["objective_rel"] if objective_informative
                          else objective_scaled <= GATE["objective_rel"])

        left_richardson, right_richardson = float(x["richardson"]), float(y["richardson"])
        left_fd = abs(left_kernel - left_richardson)
        right_fd = abs(right_kernel - right_richardson)
        left_fd_metric = (left_fd / max(abs(left_kernel), abs(left_richardson), 1e-300)
                          if informative else left_fd / scale)
        right_fd_metric = (right_fd / max(abs(right_kernel), abs(right_richardson), 1e-300)
                           if informative else right_fd / scale)
        fd_limit = GATE["kernel_fd_rel"] if kernel_fd_limit is None else kernel_fd_limit
        fd_pass = (left_fd_metric <= fd_limit and right_fd_metric <= fd_limit)
        row_pass = derivative_pass and objective_pass and fd_pass
        pair_pass &= row_pass
        row = {
            "pair": label, "observable": observable, "modulus": key[1],
            "informative": informative,
            "dspec_objective": x["objective"], "yspec_objective": y["objective"],
            "objective_relative": objective_relative,
            "objective_scaled_absolute": objective_scaled,
            "objective_informative": objective_informative,
            "dspec_kernel": left_kernel, "yspec_kernel": right_kernel,
            "derivative_relative": derivative_relative,
            "derivative_scaled_absolute": derivative_scaled,
            "dspec_fd_richardson": x["richardson"],
            "yspec_fd_richardson": y["richardson"],
            "dspec_kernel_fd_relative": x["relative_kernel_fd"],
            "yspec_kernel_fd_relative": y["relative_kernel_fd"],
            "dspec_fd_metric": left_fd_metric, "yspec_fd_metric": right_fd_metric,
            "dspec_source_delta": x["source_delta"],
            "dspec_volume_kernel": x["volume_kernel"],
            "dspec_positive_volume_kernel": x["positive_volume_kernel"],
            "dspec_negative_volume_kernel": x["negative_volume_kernel"],
            "dspec_source_kernel": x["source_kernel"],
            "dspec_source_contribution": x["source_contribution"],
            "yspec_source_kernel": y["source_kernel"],
            "yspec_source_delta": y["source_delta"],
            "yspec_source_contribution": y["source_contribution"],
            "yspec_volume_kernel": y["volume_kernel"],
            "yspec_positive_volume_kernel": y["positive_volume_kernel"],
            "yspec_negative_volume_kernel": y["negative_volume_kernel"],
            "dspec_total_kernel": x["total_kernel"],
            "yspec_total_kernel": y["total_kernel"],
            "objective_pass": objective_pass, "kernel_fd_pass": fd_pass,
            "derivative_pass": derivative_pass, "pass": row_pass,
        }
        for field in ("fd_h_0_08", "fd_h_0_04", "fd_h_0_02", "fd_h_0_01",
                      "fd_h_0_005", "fd_h_0_0025", "r4_0_08_0_04",
                      "r4_0_04_0_02", "r4_0_02_0_01", "r4_0_01_0_005",
                      "r4_0_005_0_0025", "r6_0_08_0_04_0_02",
                      "r6_0_04_0_02_0_01", "r6_0_02_0_01_0_005",
                      "r6_0_01_0_005_0_0025"):
            if x.get(field, "") and y.get(field, ""):
                row["dspec_"+field] = x[field]
                row["yspec_"+field] = y[field]
        merged.append(row)
        rows.append(row)

    return ({"pass": pair_pass, "metadata": metadata,
             "omega_bin1_relative": omega,
             "spectra_by_component": spectra_metrics,
             "spectra_pass": spectra_pass, "spectra_detail": spectrum_rows,
             "directions": merged}, rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dspec-source", type=Path, default=DSPEC_DEFAULT)
    ap.add_argument("--yspec-source", type=Path, default=YSPEC_DEFAULT)
    ap.add_argument("--dspec-deps", type=Path, default=DEPS_DSPEC)
    ap.add_argument("--yspec-deps", type=Path, default=DEPS_YSPEC)
    ap.add_argument("--planetmodel", type=Path, default=PLANETMODEL)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--build-root", type=Path, default=Path("/tmp/matched-elastic-kernel-build"))
    ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args()
    out = a.output.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SystemExit(f"output must be empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    if git(a.yspec_source, "rev-parse", "HEAD") != YSPEC_COMMIT:
        raise SystemExit("YSpecpp source is not the requested immutable commit")
    if git(a.yspec_source, "status", "--porcelain", "--untracked-files=all"):
        raise SystemExit("YSpecpp source checkout must be clean")
    dependency_roots = {
        "dspec": (a.dspec_deps, a.dspec_deps / "planetmodel-src"),
        "yspecpp": (a.yspec_deps, a.planetmodel),
    }
    dependencies = {}
    for solver, (dependency_root, planetmodel_root) in dependency_roots.items():
        shared = {name: clean_dependency(dependency_root / name, pin)
                  for name, pin in DEP_PINS.items()}
        planetmodel = clean_dependency(planetmodel_root, PLANETMODEL_COMMITS[solver])
        dependencies[solver] = {"dependency_root": str(dependency_root),
                                "shared": shared, "planetmodel": planetmodel}
    env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
               MKL_NUM_THREADS="2")
    model = HERE / "homogeneous_shell.deck"
    same_element_model = HERE / "homogeneous_same_element.deck"
    source_node_model = HERE / "homogeneous_source_node.deck"
    source_local_model = HERE / "homogeneous_source_local.deck"
    source_local_za_model = HERE / "homogeneous_source_local_za.deck"
    inputs = out / "inputs"
    inputs.mkdir()
    saved_model = inputs / model.name
    saved_same_element_model = inputs / same_element_model.name
    saved_source_node_model = inputs / source_node_model.name
    saved_source_local_model = inputs / source_local_model.name
    saved_source_local_za_model = inputs / source_local_za_model.name
    for source_fixture, saved_fixture in ((model, saved_model),
                                           (same_element_model, saved_same_element_model),
                                           (source_node_model, saved_source_node_model),
                                           (source_local_model, saved_source_local_model),
                                           (source_local_za_model, saved_source_local_za_model)):
        shutil.copy2(source_fixture, saved_fixture)
        if file_hash(saved_fixture) != file_hash(source_fixture):
            raise RuntimeError(f"saved input fixture does not match source fixture: {source_fixture.name}")
    historical_root = DSPEC_DEFAULT / "build/sensitivity-kernels-yspecpp-20261023/accepted"
    historical = inputs / "historical_legacy_unaligned"
    historical.mkdir()
    historical_files = {}
    for name in ("summary.json", "provenance.json"):
        source_file = historical_root / name
        saved_file = historical / name
        shutil.copy2(source_file, saved_file)
        source_hash, saved_hash = file_hash(source_file), file_hash(saved_file)
        if source_hash != saved_hash:
            raise RuntimeError(f"historical reference was not copied byte-for-byte: {name}")
        historical_files[name] = {"sha256": saved_hash, "source_sha256": source_hash,
                                  "saved_sha256": saved_hash, "source": str(source_file),
                                  "saved": str(saved_file.relative_to(out))}
    build_commands, runs = [], {}
    configs = [("dspec", "coarse", {"maxstep": .25, "nq": 4, "support_layer": 1}),
               ("dspec", "fine", {"maxstep": .05, "nq": 8, "support_layer": 1}),
               ("yspecpp", "baseline", {"tolerance": 1e-4, "support_layer": 1}),
               ("yspecpp", "tight", {"tolerance": 1e-9, "support_layer": 1})]
    for solver in ("dspec", "yspecpp"):
        build = a.build_root.resolve() / solver
        cfg = ["cmake", "-S", str(HERE), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release", "-DSOLVER="+solver,
               "-DDSPEC_SOURCE_DIR="+str(a.dspec_source),
               "-DYSPEC_SOURCE_DIR="+str(a.yspec_source)]
        dependency = a.dspec_deps if solver == "dspec" else a.yspec_deps
        cfg += ["-DDEPENDENCY_ROOT="+str(dependency)]
        planetmodel_root = dependency_roots[solver][1]
        cfg += ["-DPLANETMODEL_SOURCE_DIR="+str(planetmodel_root)]
        build_commands.append(shlex.join(cfg))
        command(cfg, env)
        build_cmd = ["cmake", "--build", str(build), "--target", "kernel_adapter",
                     "-j", str(a.jobs)]
        build_commands.append(shlex.join(build_cmd)); command(build_cmd, env)
        exe = build / "kernel_adapter"
        for _, label, values in [x for x in configs if x[0] == solver]:
            dest = out / "adapters" / solver / label
            dest.mkdir(parents=True)
            if solver == "dspec":
                cmd = [str(exe), str(saved_model), str(dest), str(values["maxstep"]), str(values["nq"]), str(values["support_layer"])]
            else:
                cmd = [str(exe), str(saved_model), str(dest), str(values["tolerance"]), label, str(values["support_layer"])]
            build_commands.append(shlex.join(cmd)); command(cmd, env)
            runs[f"{solver}_{label}"] = {"path": str(dest.relative_to(out)), **values}
        diagnostic_solver_dir = out / "diagnostics/same-element/adapters" / solver
        diagnostic_solver_dir.mkdir(parents=True)
        diagnostic_label = "fine" if solver == "dspec" else "tight"
        diagnostic_dir = diagnostic_solver_dir / diagnostic_label
        if solver == "dspec":
            cmd = [str(exe), str(saved_same_element_model), str(diagnostic_dir), "0.05", "8", "1"]
        else:
            cmd = [str(exe), str(saved_same_element_model), str(diagnostic_dir), "1e-9", "tight", "1"]
        build_commands.append(shlex.join(cmd))
        command(cmd, env)
        runs[f"{solver}_same_element_{diagnostic_label}"] = {
            "path": str(diagnostic_dir.relative_to(out)),
            **({"maxstep": 0.05, "nq": 8} if solver == "dspec" else {"tolerance": 1e-9})}
        for case, fixture, layer, selector in (("source-node", saved_source_node_model, 1, None),
                                               ("source-local-other", saved_source_local_model, 3, "source-local-other"),
                                               ("source-local-za", saved_source_local_za_model, 4, "source-local-za")):
            case_solver_dir = out / f"diagnostics/{case}/adapters" / solver
            case_solver_dir.mkdir(parents=True)
            case_label = "fine" if solver == "dspec" else "tight"
            case_dir = case_solver_dir / case_label
            if solver == "dspec":
                cmd = [str(exe), str(fixture), str(case_dir), "0.05", "8", str(layer)]
            else:
                cmd = [str(exe), str(fixture), str(case_dir), "1e-9", "tight", str(layer)]
            if selector:
                cmd.append(selector)
            build_commands.append(shlex.join(cmd)); command(cmd, env)
            runs[f"{solver}_{case.replace('-', '_')}_{case_label}"] = {
                "path": str(case_dir.relative_to(out)), "support_layer": layer,
                **({"maxstep": 0.05, "nq": 8} if solver == "dspec" else {"tolerance": 1e-9})}
    baseline_d = out/"adapters/dspec/coarse"
    fine_d = out/"adapters/dspec/fine"
    baseline_y = out/"adapters/yspecpp/baseline"
    tight_y = out/"adapters/yspecpp/tight"
    pairs = [("coarse_vs_baseline", baseline_d, baseline_y),
             ("fine_vs_tight", fine_d, tight_y)]
    same_element_d = out / "diagnostics/same-element/adapters/dspec/fine"
    same_element_y = out / "diagnostics/same-element/adapters/yspecpp/tight"
    diagnostic_pair = ("same_element", same_element_d, same_element_y)
    source_node_d = out / "diagnostics/source-node/adapters/dspec/fine"
    source_node_y = out / "diagnostics/source-node/adapters/yspecpp/tight"
    source_local_other_d = out / "diagnostics/source-local-other/adapters/dspec/fine"
    source_local_other_y = out / "diagnostics/source-local-other/adapters/yspecpp/tight"
    source_local_za_d = out / "diagnostics/source-local-za/adapters/dspec/fine"
    source_local_za_y = out / "diagnostics/source-local-za/adapters/yspecpp/tight"
    comparisons, rows = {}, []
    primary_pass = True
    source_local_maps = [derivatives(path / "derivatives.csv") for path in
                         (source_local_other_d, source_local_other_y, source_local_za_d, source_local_za_y)]
    source_local_scales = {
        obs: (max(abs(float(v["kernel"])) for values in source_local_maps for v in values.values()
                  if v["observable"] == obs),
              max(abs(float(v["objective"])) for values in source_local_maps for v in values.values()
                  if v["observable"] == obs))
        for obs in ("Z", "N")
    }
    comparison_specs = [(label, left_dir, right_dir, None) for label, left_dir, right_dir in pairs]
    comparison_specs += [("same_element", same_element_d, same_element_y, None),
                         ("source_node", source_node_d, source_node_y, None),
                         ("source_local_other", source_local_other_d, source_local_other_y,
                          GATE["source_local_kernel_fd_rel"]),
                         ("source_local_za", source_local_za_d, source_local_za_y,
                          GATE["source_local_kernel_fd_rel"])]
    for label, left_dir, right_dir, fd_limit in comparison_specs:
        scales = source_local_scales if label.startswith("source_local_") else None
        comparison, pair_rows = compare_pair(label, left_dir, right_dir, fd_limit, scales)
        comparisons[label] = comparison
        rows.extend(pair_rows)
        if label in ("fine_vs_tight", "same_element", "source_local_other", "source_local_za"):
            primary_pass &= comparison["pass"]
    source_node_diag = comparisons["source_node"]
    source_node_objectives = {obs: max(row["objective_relative"] for row in source_node_diag["directions"]
                                        if row["observable"] == obs) for obs in ("Z", "N")}
    source_node_diag["waveform_objectives_by_observable"] = source_node_objectives
    source_node_diag["waveform_metadata_pass"] = all(item["pass"] for item in source_node_diag["metadata"])
    source_node_diag["waveform_frequency_pass"] = source_node_diag["omega_bin1_relative"] <= GATE["omega_rel"]
    source_node_diag["waveform_pass"] = (source_node_diag["waveform_metadata_pass"] and
                                          source_node_diag["waveform_frequency_pass"] and
                                          source_node_diag["spectra_pass"] and
                                          all(v <= GATE["objective_rel"] for v in source_node_objectives.values()))
    primary_pass &= source_node_diag["waveform_pass"]
    for case_dir in (fine_d, same_element_d, source_node_d, source_local_other_d, source_local_za_d):
        case_meta = json.loads((case_dir / "meta.json").read_text())
        if not case_meta.get("source_aligned", False):
            raise RuntimeError(f"DSpec source-aware mesh was not reported aligned: {case_dir}")
        bounds = case_meta["source_element_bounds_m"]
        if not (bounds[0] < bounds[1] and
                min(abs(case_meta["source_radius_m"]-bounds[0]),
                    abs(case_meta["source_radius_m"]-bounds[1])) <= 1e-6):
            raise RuntimeError(f"DSpec source is not a boundary of a nonzero-width element: {case_dir}")
    write_csv(out/"comparison.csv",rows)
    fd_rows=[]
    fd_fields=("fd_h_0_08", "fd_h_0_04", "fd_h_0_02", "fd_h_0_01", "fd_h_0_005", "fd_h_0_0025",
               "r4_0_08_0_04", "r4_0_04_0_02", "r4_0_02_0_01", "r4_0_01_0_005", "r4_0_005_0_0025",
               "r6_0_08_0_04_0_02", "r6_0_04_0_02_0_01", "r6_0_02_0_01_0_005", "r6_0_01_0_005_0_0025")
    source_rows = comparisons["source_local_other"]["directions"] + comparisons["source_local_za"]["directions"]
    if {(row["observable"], row["modulus"]) for row in source_rows} != {(o, m) for o in ("Z", "N") for m in "ACFLN"}:
        raise RuntimeError("split source-local cases do not provide exactly ten unique observable/modulus directions")
    comparisons["source_local"] = {"pass": comparisons["source_local_other"]["pass"] and comparisons["source_local_za"]["pass"],
                                   "directions": source_rows,
                                   "setups": {"other_nine": comparisons["source_local_other"],
                                              "revised_za": comparisons["source_local_za"]}}
    for row in source_rows:
        for solver in ("dspec", "yspec"):
            for field in fd_fields:
                fd_rows.append({"observable":row["observable"],"modulus":row["modulus"],
                                "solver":solver,"sequence":field,"value":row.get(solver+"_"+field,"")})
    write_csv(out/"source-local-fd-sequence.csv",fd_rows)
    fine_meta = json.loads((fine_d / "meta.json").read_text())
    tight_meta = json.loads((tight_y / "meta.json").read_text())
    source_node_meta = json.loads((source_node_d / "meta.json").read_text())
    source_local_other_meta = json.loads((source_local_other_d / "meta.json").read_text())
    source_local_other_meta_y = json.loads((source_local_other_y / "meta.json").read_text())
    source_local_meta = json.loads((source_local_za_d / "meta.json").read_text())
    source_local_meta_y = json.loads((source_local_za_y / "meta.json").read_text())
    expected_local_bounds = [6200000.0, 6371000.0]
    if any(abs(float(x)-y) > 1e-6 for x,y in zip(source_local_other_meta["support_bounds_m"], expected_local_bounds)):
        raise RuntimeError(f"original source-local support bounds differ from fixture: {source_local_other_meta['support_bounds_m']}")
    if any(abs(float(x)-y) > 1e-6 for x,y in zip(source_local_meta["support_bounds_m"], [6250000.0, 6371000.0])):
        raise RuntimeError(f"revised ZA support bounds differ from fixture: {source_local_meta['support_bounds_m']}")
    if not source_local_meta["support_bounds_m"][0] < source_local_meta["source_radius_m"] < source_local_meta["support_bounds_m"][1]:
        raise RuntimeError("source-local support does not contain the source strictly inside its boundaries")
    if source_local_meta["support_bounds_m"] != source_local_meta_y["support_bounds_m"]:
        raise RuntimeError("source-local physical support differs between solvers")
    if not source_local_meta_y["source_local_support_case"] or source_local_meta_y["kernel_mesh_max_element_m"] != 1000.0:
        raise RuntimeError("YSpec source-local receiver-break workaround was not applied")
    if abs((source_local_meta_y["receiver_radius_m"]-source_local_meta_y["source_radius_m"])/1000.0-47.0) > 1e-8:
        raise RuntimeError("source-local YSpec one-kilometre mesh no longer spans 47 source-to-receiver intervals")
    if abs(source_node_meta["source_radius_m"]-6321000.0) > 1e-6:
        raise RuntimeError("source-node fixture no longer places its duplicated boundary at the source")
    # Report physical setup from the exact adapter metadata and saved deck used
    # above; the Python harness does not maintain a second copy of the model.
    deck_rows = [line.split() for line in saved_model.read_text().splitlines()[3:]
                 if len(line.split()) >= 6]
    model_boundaries = sorted({float(row[0]) for row in deck_rows})
    configuration = {
        "radius_m": fine_meta["radius_m"],
        "model_boundaries_m": model_boundaries,
        "shell_m": [model_boundaries[1], model_boundaries[2]],
        "density_kg_m3": fine_meta["density_kg_m3"],
        "vp_m_s": fine_meta["vp_m_s"], "vs_m_s": fine_meta["vs_m_s"],
        "moduli_Pa": dict(zip("ACFLN", fine_meta["shell_A_C_F_L_N_Pa"])),
        "source_geocentric": fine_meta["source_geocentric"],
        "receiver_geocentric": fine_meta["receiver_geocentric"],
        "source_geographic": dict(zip(("latitude_deg", "longitude_deg", "depth_m"),
                                       tight_meta["source_geographic"])),
        "receiver_geographic": dict(zip(("latitude_deg", "longitude_deg", "depth_m"),
                                         tight_meta["receiver_geographic"])),
        "moment_tensor_Nm": fine_meta["moment_tensor_Nm"],
        "dt_s": fine_meta["dt_s"], "fft_samples": fine_meta["nt"],
        "frequency_bins": sorted({int(row["bin"]) for row in
                                   read_csv(fine_d / "spectra.csv")}),
        "frequency_bin_hz": [float(row["frequency_hz"])
                             for row in {int(item["bin"]): item for item in
                                         read_csv(fine_d / "spectra.csv")}.values()],
        "requested_band_hz": tight_meta["requested_band_hz"],
        "q_bins": tight_meta["q_bins"], "q_values": tight_meta["q_values"],
        "window_samples": fine_meta["window_samples"], "window": fine_meta["window"],
        "damping_s": fine_meta["damping_s"], "physical_G_m3_kg_s2": fine_meta["physical_G"],
        "gravity": fine_meta["gravity"], "attenuation": fine_meta["attenuation"],
        "degrees": fine_meta["degrees"],
        "source_local_support_m": source_local_meta["support_bounds_m"],
        "source_local_source_radius_m": source_local_meta["source_radius_m"],
        "source_local_original_nine_support_m": source_local_other_meta["support_bounds_m"],
        "source_local_za_support_m": source_local_meta["support_bounds_m"],
        "source_node_radius_m": source_node_meta["source_radius_m"],
        "finite_difference_fractional_steps": [0.08, 0.04, 0.02, 0.01, 0.005, 0.0025],
        "finite_difference_extrapolants": {
            "R4(h,h/2)": "(4 D(h/2)-D(h))/3 for every adjacent pair",
            "R6(h,h/2,h/4)": "(16 R4(h/2,h/4)-R4(h,h/2))/15 for all four triplets",
            "kernel_gate_estimate": "retained R4(0.005,0.0025), computed from h=0.01,0.005,0.0025; larger-step R4/R6 are diagnostics"},
        "yspec_source_local_mesh_workaround": {
            "commit": YSPEC_COMMIT, "max_element_m": source_local_meta_y["kernel_mesh_max_element_m"],
            "order": source_local_meta_y["kernel_mesh_order"], "production_source_modified": False,
            "reason": "Pinned YSpec Mesh break collection compares later candidates against a growing edges.back(); after appending an interior source, the receiver above it is omitted. The benchmark-only 1 km maximum element resolves the 47 km source-to-receiver interval without editing YSpecpp."},
    }
    (out/"provenance.json").write_text(json.dumps({
        "experiment":"matched 1D full-gravity elastic modulus kernels",
        "created_utc":subprocess.check_output(["date","-u","+%Y-%m-%dT%H:%M:%SZ"],text=True).strip(),
        "dspec":source_identity(a.dspec_source),"yspecpp":source_identity(a.yspec_source),
        "verified_dependency_roots":dependencies,
        "commands":build_commands,"runs":runs,
        "fixtures":{"primary":{"sha256":file_hash(saved_model),"source":str(model),"saved":str(saved_model)},
                    "same_element":{"sha256":file_hash(saved_same_element_model),
                                    "source":str(same_element_model),"saved":str(saved_same_element_model)},
                    "source_node":{"sha256":file_hash(saved_source_node_model),
                                   "source":str(source_node_model),"saved":str(saved_source_node_model)},
                    "source_local":{"sha256":file_hash(saved_source_local_model),
                                    "source":str(source_local_model),"saved":str(saved_source_local_model)},
                    "source_local_za":{"sha256":file_hash(saved_source_local_za_model),
                                        "source":str(source_local_za_model),"saved":str(saved_source_local_za_model)},
                    "historical_legacy_unaligned": historical_files},
        "gates":GATE,
        "configuration": configuration}, indent=2)+"\n")
    # Plots are descriptive; the tables above carry all gate values.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for pair,left_dir,right_dir in pairs:
        dl,dr=derivatives(left_dir/"derivatives.csv"),derivatives(right_dir/"derivatives.csv")
        fig,ax=plt.subplots(figsize=(9,4.5)); labels=[]; x=[]; y=[]
        for key in sorted(dl):
            labels.append(key[0]+" "+key[1]); x.append(float(dl[key]["kernel"])); y.append(float(dr[key]["kernel"]))
        idx=list(range(len(labels))); ax.plot(idx,x,"o",label="DSpecM1D integrated kernel"); ax.plot(idx,y,"x",label="YSpecpp integrated kernel")
        ax.set_xticks(idx,labels,rotation=35,ha="right")
        ax.set_ylabel("dJ/dh [m s]")
        ax.grid(True,alpha=.25); ax.legend(); fig.tight_layout()
        fig.savefig(out/(pair+"-derivatives.png"),dpi=150); plt.close(fig)
    ddata=read_csv(fine_d/"samples.csv"); ydata=read_csv(tight_y/"samples.csv")
    fig,axes=plt.subplots(2,5,figsize=(18,6))
    for ax,(obs,mod) in zip(axes.flat,[(o,m) for o in ("Z","N") for m in "ACFLN"]):
        for data,label in ((ddata,"DSpecM1D fine"),(ydata,"YSpecpp tight")):
            group=[r for r in data if r["observable"]==obs and r["modulus"]==mod]
            ax.plot([float(r["radius_m"])/1000 for r in group],
                    [float(r["kernel_dr_per_Pa_m"]) for r in group],linestyle="None",marker=".",markersize=2,label=label)
        ax.set_title(obs+" "+mod); ax.set_xlabel("radius (km)"); ax.grid(True,alpha=.2)
        ax.set_ylabel("K_X [s Pa$^{-1}$]")
    axes[0,0].legend(fontsize=8); fig.tight_layout(); fig.savefig(out/"fine-tight-sampled-kernels.png",dpi=150); plt.close(fig)
    local_ddata=read_csv(source_local_other_d/"samples.csv") + read_csv(source_local_za_d/"samples.csv")
    local_ydata=read_csv(source_local_other_y/"samples.csv") + read_csv(source_local_za_y/"samples.csv")
    fig,axes=plt.subplots(2,5,figsize=(18,6))
    for ax,(obs,mod) in zip(axes.flat,[(o,m) for o in ("Z","N") for m in "ACFLN"]):
        for data,label in ((local_ddata,"DSpecM1D fine"),(local_ydata,"YSpecpp tight")):
            group=[r for r in data if r["observable"]==obs and r["modulus"]==mod]
            ax.plot([float(r["radius_m"])/1000 for r in group],
                    [float(r["kernel_dr_per_Pa_m"]) for r in group],linestyle="None",marker=".",markersize=2,label=label)
        if obs == "Z" and mod == "A":
            ax.axvspan(6250,6371,color="orange",alpha=.12,label="revised Z/A support")
        else:
            ax.axvspan(6200,6371,color="orange",alpha=.12,label="original support")
        ax.axvline(6321,color="black",linestyle=":",linewidth=.8,label="source")
        ax.set_title(obs+" "+mod); ax.set_xlabel("radius (km)"); ax.grid(True,alpha=.2)
        ax.set_ylabel("K_X [s Pa$^{-1}$]")
    axes[0,0].legend(fontsize=7); fig.tight_layout(); fig.savefig(out/"source-local-sampled-kernels.png",dpi=150); plt.close(fig)
    diagnostic_d = json.loads((same_element_d / "meta.json").read_text())
    diagnostic_compare = comparisons["same_element"]
    diagnostic_objectives = {}
    for obs in ("Z", "N"):
        diagnostic_objectives[obs] = max(
            row["objective_relative"] for row in diagnostic_compare["directions"]
            if row["observable"] == obs)
    same_element_diagnostic = {
        "fixture_sha256": file_hash(saved_same_element_model),
        "source_element": diagnostic_d["source_element"],
        "source_element_bounds_m": diagnostic_d["source_element_bounds_m"],
        "receiver_element": diagnostic_d["receiver_element"],
        "receiver_element_bounds_m": diagnostic_d["receiver_element_bounds_m"],
        "spectra_by_component": diagnostic_compare["spectra_by_component"],
        "objective_relative_max_by_observable": diagnostic_objectives,
        "max_integrated_derivative_relative": max(
            row["derivative_relative"] for row in diagnostic_compare["directions"]),
        "source_aligned": diagnostic_d["source_aligned"],
        "source_radius_m": diagnostic_d["source_radius_m"],
        "historical_summary_sha256": historical_files["summary.json"]["sha256"],
        "interpretation": "This fresh run applies the repaired source-aware mesh to the original homogeneous fixture where source and receiver formerly shared one radial element. It compares against frozen legacy unaligned evidence and makes no lower-level contact-mechanism claim.",
    }
    summary={"gate_status":"PASS" if primary_pass else "FAIL",
             "primary_pairs":["fine_vs_tight","same_element","source_local_other","source_local_za","source_node_waveform"],
             "same_element_source_aware":same_element_diagnostic,
             "source_node_diagnostic":{"source_radius_m":source_node_meta["source_radius_m"],
                "source_element_bounds_m":source_node_meta["source_element_bounds_m"],
                "receiver_element_bounds_m":source_node_meta["receiver_element_bounds_m"],
                "identical_material_partition":True,"waveform_pass":source_node_diag["waveform_pass"]},
             "source_local_original_nine_support_m":[6200000.0,6371000.0],
             "source_local_za_support_m":[6250000.0,6371000.0],
             "source_local_source_radius_m":source_local_meta["source_radius_m"],
             "same_element_diagnostic":same_element_diagnostic,"comparisons":comparisons}
    za_row=next(row for row in source_rows if row["observable"]=="Z" and row["modulus"]=="A")
    def cancellation(positive: float, negative: float, net: float) -> dict:
        return {"positive_volume_integral_ms":positive,"negative_volume_integral_ms":negative,
                "net_volume_integral_ms":net,
                "absolute_contribution_over_net":(positive+abs(negative))/max(abs(net),1e-300)}
    summary["source_local_za_fd_replacement"]={"supersedes":"noise-dominated original-support Z/A own-FD row",
        "partition_added_m":6250000.0,"support_bounds_m":[6250000.0,6371000.0],
        "source_radius_m":source_local_meta["source_radius_m"],"perturbed_modulus":"A",
        "gate_unchanged":True}
    summary["source_local_za_signal_conditioning"]={
        "dspec":cancellation(float(za_row["dspec_positive_volume_kernel"]),
                              float(za_row["dspec_negative_volume_kernel"]),
                              float(za_row["dspec_volume_kernel"])),
        "yspec":cancellation(float(za_row["yspec_positive_volume_kernel"]),
                              float(za_row["yspec_negative_volume_kernel"]),
                              float(za_row["yspec_volume_kernel"]))}
    za_conditioning = summary["source_local_za_signal_conditioning"]
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    status="PASS" if primary_pass else "FAIL"

    def sig(value: float) -> str:
        return f"{value:.4g}"

    def table_rows(rows: list[dict]) -> str:
        return "".join(
            "<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in values) + "</tr>"
            for values in ([r["pair"], r["observable"], r["modulus"],
                            sig(r["dspec_kernel"]), sig(r["yspec_kernel"]),
                            sig(r["derivative_relative"]), sig(r["dspec_fd_metric"]),
                            sig(r["yspec_fd_metric"]), "informative" if r["informative"] else "weak"]
                           for r in rows))

    source_rows = comparisons["source_local"]["directions"]
    def source_point_rows(rows: list[dict]) -> str:
        return "".join("<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in values) + "</tr>"
                       for r in rows for values in [[r["pair"], r["observable"], r["modulus"],
                           sig(float(r["dspec_source_kernel"])), sig(float(r["dspec_source_delta"])),
                           sig(float(r["dspec_source_contribution"])),
                           sig(float(r["yspec_source_kernel"])), sig(float(r["yspec_source_delta"])),
                           sig(float(r["yspec_source_contribution"]))]])
    def source_total_rows(rows: list[dict]) -> str:
        return "".join("<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in values) + "</tr>"
                       for r in rows for values in [[r["pair"], r["observable"], r["modulus"],
                           sig(float(r["dspec_volume_kernel"])), sig(float(r["dspec_total_kernel"])),
                           sig(float(r["dspec_fd_richardson"])),
                           sig(float(r["yspec_volume_kernel"])), sig(float(r["yspec_total_kernel"])),
                           sig(float(r["yspec_fd_richardson"])),
                           sig(float(r["derivative_relative"]))]])
    fd_steps=("fd_h_0_08","fd_h_0_04","fd_h_0_02","fd_h_0_01","fd_h_0_005","fd_h_0_0025")
    fd_r4=("r4_0_08_0_04","r4_0_04_0_02","r4_0_02_0_01","r4_0_01_0_005","r4_0_005_0_0025")
    fd_r6=("r6_0_08_0_04_0_02","r6_0_04_0_02_0_01","r6_0_02_0_01_0_005","r6_0_01_0_005_0_0025")
    fd_sequence_html="".join("<tr>"+"".join(f"<td>{escape(str(v))}</td>" for v in
        [r["pair"],r["observable"],r["modulus"],solver,*[sig(float(r[solver+"_"+f])) for f in fd_steps+fd_r4+fd_r6]])+"</tr>"
        for r in source_rows for solver in ("dspec","yspec"))

    primary_rows = comparisons["fine_vs_tight"]["directions"]
    coarse_rows = comparisons["coarse_vs_baseline"]["directions"]
    coarse_meta = json.loads((baseline_d / "meta.json").read_text())
    baseline_meta = json.loads((baseline_y / "meta.json").read_text())
    dspec_meta = fine_meta
    yspec_meta = tight_meta
    bins = configuration["frequency_bin_hz"]
    moduli = configuration["moduli_Pa"]
    geo_src = configuration["source_geographic"]
    geo_rec = configuration["receiver_geographic"]
    diagnostic_spectra = same_element_diagnostic["spectra_by_component"]
    primary_spectra = comparisons["fine_vs_tight"]["spectra_by_component"]
    source_node_spectra = comparisons["source_node"]["spectra_by_component"]
    source_node_obj = comparisons["source_node"]["waveform_objectives_by_observable"]
    source_local_meta_d = json.loads((source_local_za_d/"meta.json").read_text())
    source_local_meta_y = json.loads((source_local_za_y/"meta.json").read_text())
    historical_summary = json.loads((historical/"summary.json").read_text())
    historical_element = historical_summary["same_element_diagnostic"]
    def embedded_image(name: str, alt: str) -> str:
        payload = base64.b64encode((out/name).read_bytes()).decode("ascii")
        return f'<img alt="{escape(alt)}" src="data:image/png;base64,{payload}">'
    primary_plot = embedded_image("fine_vs_tight-derivatives.png", "Remote-shell derivative comparison")
    fig,ax=plt.subplots(figsize=(9,4.5)); labels=[]; x=[]; y=[]
    for row in source_rows:
        labels.append(row["observable"]+" "+row["modulus"])
        x.append(float(row["dspec_kernel"])); y.append(float(row["yspec_kernel"]))
    idx=list(range(len(labels))); ax.plot(idx,x,"o",label="DSpecM1D integrated kernel"); ax.plot(idx,y,"x",label="YSpecpp integrated kernel")
    ax.set_xticks(idx,labels,rotation=35,ha="right"); ax.set_ylabel("dJ/dh [m s]")
    ax.grid(True,alpha=.25); ax.legend(); fig.tight_layout()
    fig.savefig(out/"source_local-derivatives.png",dpi=150); plt.close(fig)
    local_plot = embedded_image("source_local-derivatives.png", "Source-local total derivative comparison")
    primary_curve_plot = embedded_image("fine-tight-sampled-kernels.png", "Remote-shell full radial kernel curves")
    local_curve_plot = embedded_image("source-local-sampled-kernels.png", "Source-local full radial kernel curves")
    report = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Matched elastic kernels</title><style>
body{{font:16px system-ui;max-width:1280px;margin:2rem auto;padding:0 1rem;color:#222}}
.table-wrap{{max-width:100%;overflow-x:auto}}table{{border-collapse:collapse;font-size:14px;min-width:1050px}}
td,th{{border:1px solid #bbb;padding:.35rem .5rem;white-space:nowrap}}th{{background:#eee}}
img{{max-width:100%}}code{{background:#eee;padding:.1rem .2rem}}
</style></head><body>
<h1>Matched 1D source-treatment and elastic-kernel comparison</h1>
<p><strong>Campaign result: {status}</strong>. The gate covers the source-aware fine/tight remote-shell pair, the repaired original same-element fixture, nine original source-local directions, the revised source-local Z/A direction, and the source-node waveform diagnostic. The original noise-dominated Z/A finite-difference case is superseded by the user-requested setup change; its saved failure evidence remains archived and was not rerun.</p>
<p>Numerical settings from the adapter run metadata: DSpecM1D coarse uses SEM maximum-step fraction {sig(coarse_meta['maxstep_fraction'])} with nq={coarse_meta['nq']}, and fine uses {sig(dspec_meta['maxstep_fraction'])} with nq={dspec_meta['nq']}. YSpecpp baseline uses ODE tolerance {sig(baseline_meta['tolerance'])}; tight uses {sig(yspec_meta['tolerance'])}. The primary gate compares the latter fine/tight pair.</p>
<p>The saved homogeneous solid deck has radius {sig(configuration['radius_m'])} m and remote-shell boundaries {configuration['model_boundaries_m']} m. For that shell case, a fractional change h is applied independently to A/C/F/L/N only in {configuration['shell_m']} m; density and other moduli are fixed. Actual properties from DSpec metadata: density {sig(configuration['density_kg_m3'])} kg m⁻³, Vp {sig(configuration['vp_m_s'])} m s⁻¹, Vs {sig(configuration['vs_m_s'])} m s⁻¹, moduli A/C/F/L/N {', '.join(f'{k}={sig(v)}' for k,v in moduli.items())} Pa. Source and receiver lie outside this remote shell, so its source modulus change is zero.</p>
<p>Source geocentric [lat,lon,depth] {configuration['source_geocentric']}; YSpec geographic [lat,lon,depth] [{sig(geo_src['latitude_deg'])}, {sig(geo_src['longitude_deg'])}, {sig(geo_src['depth_m'])}]. Receiver geocentric {configuration['receiver_geocentric']}; geographic [{sig(geo_rec['latitude_deg'])}, {sig(geo_rec['longitude_deg'])}, {sig(geo_rec['depth_m'])}]. Tensor [rr,rt,rp,tt,tp,pp] {configuration['moment_tensor_Nm']} N m. Gravity {configuration['gravity']}, attenuation {configuration['attenuation']}, degrees {configuration['degrees']}, physical G {sig(configuration['physical_G_m3_kg_s2'])} m³ kg⁻¹ s⁻².</p>
<p>Grid: {configuration['fft_samples']} samples at {sig(configuration['dt_s'])} s; raw bins {configuration['frequency_bins']} at {', '.join(sig(v*1e3) for v in bins)} mHz. Requested YSpec band {', '.join(sig(v*1e3) for v in configuration['requested_band_hz'])} mHz. q occupies bins {configuration['q_bins']} with real/imaginary pairs {configuration['q_values']}; the {configuration['window_samples']}-sample {configuration['window']} measurement mixes raw bin 3. Damping is {sig(configuration['damping_s'])} s⁻¹. Physical omega relative difference is {sig(comparisons['fine_vs_tight']['omega_bin1_relative'])} (limit 5e-11, allowing the legacy DSpec π literal).</p>
<p>The mapped YSpec adjoint is <code>j[n]=W[n] Σ(qᵣ cosθ − qᵢ sinθ)</code>, so <code>dt Σ j[n]u[n]</code> matches the DSpec real spectral scalar. The radial integral uses native YSpec K with its r²dr weights; sampled curves show r²K in the DSpec dr convention. Source-local finite differences use centered fractional steps 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025. All five adjacent R4 and four triplet R6 extrapolants are reported; the existing kernel gate remains tied to R4(0.005,0.0025), with no threshold change.</p>
<h2>Fine / tight: primary rows</h2><div class="table-wrap"><table><thead><tr><th>pair</th><th>obs</th><th>modulus</th><th>DSpec dJ/dh [m s]</th><th>YSpec dJ/dh [m s]</th><th>cross relative</th><th>DSpec kernel/FD</th><th>YSpec kernel/FD</th><th>signal</th></tr></thead><tbody>{table_rows(primary_rows)}</tbody></table></div>
<p>Primary raw spectrum errors (limits apply separately to L2 and peak-scaled maximum for each component): {'; '.join(f'{c}: L2 {sig(v["l2"])}, peak {sig(v["peak_scaled_max"])}' for c,v in primary_spectra.items())}. Objective relative errors by observable: {', '.join(f'{o} {sig(max(r["objective_relative"] for r in primary_rows if r["observable"]==o))}' for o in ('Z','N'))}.</p>
<h2>Gates</h2><p>Shared metadata relative/exact error ≤1e-12; physical angular frequency ≤5e-11; each Z/N/E spectrum L2 and peak-scaled maximum ≤1e-5; objective ≤1e-5; informative cross derivatives ≤1e-3; remote-shell kernel/Richardson ≤1e-4; source-local kernel/Richardson ≤1e-6. Informative means ≥1e-6 of the largest absolute derivative in that observable across the ten reported directions. Weak rows use scaled absolute errors against that common maximum. All ten directions remain visible; weak directions have no automatic pass. No limit changed.</p>
<h2>Coarse / baseline diagnostic</h2><div class="table-wrap"><table><thead><tr><th>pair</th><th>obs</th><th>modulus</th><th>DSpec dJ/dh</th><th>YSpec dJ/dh</th><th>cross relative</th><th>DSpec FD metric</th><th>YSpec FD metric</th><th>signal</th></tr></thead><tbody>{table_rows(coarse_rows)}</tbody></table></div>
<h2>Source in its original radial element</h2><p>The fresh repaired run uses the unpartitioned homogeneous fixture. The source-aware mesh reports source element {same_element_diagnostic['source_element']} {same_element_diagnostic['source_element_bounds_m']} m and receiver element {same_element_diagnostic['receiver_element']} {same_element_diagnostic['receiver_element_bounds_m']} m; both source and receiver align to nonzero-width cells. Each raw spectrum L2 error is {'; '.join(f'{c} {sig(v["l2"])}' for c,v in diagnostic_spectra.items())}; maximum Z/N objective errors are {', '.join(f'{o} {sig(v)}' for o,v in diagnostic_objectives.items())}. The frozen legacy unaligned evidence recorded Z/N/E raw spectrum L2 errors {', '.join(f'{c} {sig(v["l2"])}' for c,v in historical_element['spectra_by_component'].items())} and objective errors {', '.join(f'{o} {sig(v)}' for o,v in historical_element['objective_relative_max_by_observable'].items())}. Its unchanged summary SHA-256 is {historical_files['summary.json']['sha256']}; summary and provenance are copied byte-for-byte under inputs/historical_legacy_unaligned/.</p>
<h2>Source already on a mesh node</h2><p>An identical-material duplicated partition lies exactly at the source radius {sig(source_node_meta['source_radius_m'])} m, without a physical material contrast. Source and receiver cells are {source_node_meta['source_element_bounds_m']} m and {source_node_meta['receiver_element_bounds_m']} m. Spectrum L2 errors are {'; '.join(f'{c} {sig(v["l2"])}' for c,v in source_node_spectra.items())}; objective errors are {', '.join(f'{o} {sig(v)}' for o,v in source_node_obj.items())}. This checks the existing-node path and is not an Earth-interface test.</p>
<h2>Source-local modulus perturbations</h2><p>The nine non-Z/A checks retain their original 6200–6371 km support. The revised Z/A check uses an identical-material partition at 6250 km and A-only support 6250–6371 km. The source at {sig(source_local_meta_d['source_radius_m'])} m remains strictly inside both supports; the receiver is unchanged at {sig(source_local_meta_y['receiver_radius_m'])} m. The material remains homogeneous, and density and other moduli stay fixed. Source-local finite differences and derivative comparisons below are labelled by fixture so the revised A row is not attributed to the original setup.</p>
<p>For the revised Z/A direction, the signed positive and negative kernel-volume integrals are retained in the adapter CSVs and the signal-conditioning summary. This reports cancellation directly; it was not used to choose the fixture. YSpec uses order {source_local_meta_y['kernel_mesh_order']} and maximum element {sig(source_local_meta_y['kernel_mesh_max_element_m'])} m. Its existing 1 km benchmark-only setting resolves the 47 km source-to-receiver span into 47 intervals; the YSpecpp source remains unchanged.</p>
<h3>Source-point coefficients and contributions</h3><div class="table-wrap"><table><thead><tr><th>fixture</th><th>obs</th><th>modulus</th><th>DSpec k(source) [m s Pa⁻¹]</th><th>DSpec δX [Pa]</th><th>DSpec point [m s]</th><th>YSpec k(source) [m s Pa⁻¹]</th><th>YSpec δX [Pa]</th><th>YSpec point [m s]</th></tr></thead><tbody>{source_point_rows(source_rows)}</tbody></table></div>
<h3>Volume plus source-point totals and Richardson checks</h3><div class="table-wrap"><table><thead><tr><th>fixture</th><th>obs</th><th>modulus</th><th>DSpec volume [m s]</th><th>DSpec total [m s]</th><th>DSpec Richardson [m s]</th><th>YSpec volume [m s]</th><th>YSpec total [m s]</th><th>YSpec Richardson [m s]</th><th>cross relative</th></tr></thead><tbody>{source_total_rows(source_rows)}</tbody></table></div>
<p>For revised Z/A, signed volume-kernel contributions are DSpec positive {sig(za_conditioning['dspec']['positive_volume_integral_ms'])}, negative {sig(za_conditioning['dspec']['negative_volume_integral_ms'])}, net {sig(za_conditioning['dspec']['net_volume_integral_ms'])} m s; YSpec positive {sig(za_conditioning['yspec']['positive_volume_integral_ms'])}, negative {sig(za_conditioning['yspec']['negative_volume_integral_ms'])}, net {sig(za_conditioning['yspec']['net_volume_integral_ms'])} m s. Absolute-contribution/net ratios are {sig(za_conditioning['dspec']['absolute_contribution_over_net'])} and {sig(za_conditioning['yspec']['absolute_contribution_over_net'])}, respectively.</p>
<h3>Centered FD, R4, and R6 sequence</h3><p>The table combines nine original-support directions and the revised 6250–6371 km Z/A direction. The former 6200–6371 km Z/A finite difference is excluded from this run; its frozen failure remains in the prior evidence. For required rows, all six centered differences, five adjacent R4 estimates, and four triplet R6 estimates are reported.</p><div class="table-wrap"><table><thead><tr><th>fixture</th><th>obs</th><th>modulus</th><th>solver</th><th>D .08</th><th>D .04</th><th>D .02</th><th>D .01</th><th>D .005</th><th>D .0025</th><th>R4 .08/.04</th><th>R4 .04/.02</th><th>R4 .02/.01</th><th>R4 .01/.005</th><th>R4 .005/.0025</th><th>R6 .08/.04/.02</th><th>R6 .04/.02/.01</th><th>R6 .02/.01/.005</th><th>R6 .01/.005/.0025</th></tr></thead><tbody>{fd_sequence_html}</tbody></table></div>
<h2>Curves and data</h2><p>{primary_plot}{local_plot}</p><p>{primary_curve_plot}</p><p>{local_curve_plot}</p><p>Full precision data are in <a href="comparison.csv">comparison.csv</a>, <a href="source-local-fd-sequence.csv">source-local-fd-sequence.csv</a>, adapter CSVs under <code>adapters/</code> and <code>diagnostics/</code>, <a href="summary.json">summary.json</a>, and <a href="provenance.json">provenance.json</a>. Plots are embedded in this report. No broad convergence claim is made.</p>
</body></html>"""
    (out/"report.html").write_text(report)
    print(f"Benchmark artifacts: {out}\nGate result: {status}")
    return 0 if primary_pass else 1

if __name__ == "__main__":
    raise SystemExit(main())
