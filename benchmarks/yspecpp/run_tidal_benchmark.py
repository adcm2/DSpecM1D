#!/usr/bin/env python3
"""Run the isolated, degree-2 Figure 7 tidal-resolution diagnostic."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CRITERIA = Path(__file__).resolve().parent / "cases/paper-suite-protocol.json"
DEFAULT_MODEL = ROOT / "data/models/prem.200.no.noatten.txt"
DEFAULT_OUTPUT = ROOT / "build/phase3/paper-suite-20261020/fig7-tide"
PROFILES = ((0.01, "coarse"), (0.001, "medium"), (0.0001, "fine"))
PERIODS_H = (12.0, 48.0, 168.0, 1680.0)
PAPER_ROUNDED_BV_RATIOS = (1.3, 5.2, 18.2, 182.2)
NQ = 5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_criteria(path: Path) -> dict:
    """Validate the tracked paper-suite protocol and its Fig. 7 input anchor."""
    try:
        criteria = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read criteria file {path}: {exc}") from exc
    if not isinstance(criteria, dict):
        raise ValueError("criteria must be a JSON object")
    cases = criteria.get("case_inventory")
    if not isinstance(cases, list) or not any(
            isinstance(case, dict) and case.get("case") == "ex6/Fig7" for case in cases):
        raise ValueError("criteria case_inventory must include ex6/Fig7")
    source_hashes = criteria.get("source_basis", {}).get("source_parameter_hashes", {})
    expected_hash = source_hashes.get("ex6") if isinstance(source_hashes, dict) else None
    input_path = ROOT / "data/params/ex6.txt"
    if not isinstance(expected_hash, str) or not input_path.is_file():
        raise ValueError("criteria must provide the ex6 parameter hash and data/params/ex6.txt must exist")
    actual_hash = sha256(input_path)
    if actual_hash != expected_hash:
        raise ValueError(f"ex6 parameter hash does not match criteria: {actual_hash}")
    return criteria


def estimate_elements(model_path: Path, maxstep: float) -> dict:
    lines = model_path.read_text().splitlines()
    n_nodes = int(lines[2].split()[0])
    radii = []
    for line in lines[3:3 + n_nodes]:
        fields = line.split()
        if fields:
            radii.append(float(fields[0]))
    layers: list[list[float]] = []
    for radius in radii:
        if layers and radius == layers[-1][-1]:
            layers.append([radius])
        else:
            if not layers:
                layers.append([])
            layers[-1].append(radius)
    outer = max(radii)
    widths = [layer[-1] - layer[0] for layer in layers]
    counts = [math.ceil(width / (maxstep * outer)) for width in widths]
    ne = sum(counts)
    return {
        "model_nodes": n_nodes,
        "model_layers": len(layers),
        "layer_element_counts": counts,
        "elements": ne,
        "local_gll_nodes": ne * NQ,
        "unique_radial_dofs_approx": ne * (NQ - 1) + 1,
        "full_gravity_spheroidal_dofs_approx": 3 * (ne * (NQ - 1) + 1),
        "element_width_cap_m": maxstep * outer,
        "outer_radius_m": outer,
    }


def available_memory() -> dict:
    result: dict[str, int | None] = {"mem_available_kib": None, "cgroup_current_bytes": None}
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                result["mem_available_kib"] = int(line.split()[1])
                break
    except OSError:
        pass
    for path in (Path("/sys/fs/cgroup/memory.current"),
                 Path("/sys/fs/cgroup/memory/memory.usage_in_bytes")):
        try:
            result["cgroup_current_bytes"] = int(path.read_text().strip())
            break
        except (OSError, ValueError):
            pass
    return result


def write_synthetic_deck(path: Path) -> None:
    rho = 2.0 * 5515.0
    rows = [
        "Manufactured constant-density sphere for Fig. 7 tidal-load quadrature",
        "1 -1.0 1 2",
        "2 0 0",
        f"0 {rho:.17g} 8000 4500 0 0 8000 4500 1",
        f"2000 {rho:.17g} 8000 4500 0 0 8000 4500 1",
    ]
    path.write_text("\n".join(rows) + "\n")


def gll_nodes_weights(n: int) -> tuple[np.ndarray, np.ndarray]:
    degree = n - 1
    interior = np.polynomial.legendre.Legendre.basis(degree).deriv().roots()
    nodes = np.concatenate(([-1.0], interior, [1.0]))
    weights = 2.0 / (degree * (degree + 1) *
                     np.polynomial.legendre.Legendre.basis(degree)(nodes) ** 2)
    return nodes, weights


def manufactured_expected(a: float, elements: int) -> list[dict]:
    xi, weights = gll_nodes_weights(NQ)
    rows = []
    for element in range(elements):
        left = a * element / elements
        right = a * (element + 1) / elements
        jacobian = (right - left) / 2.0
        radii = (left + right) / 2.0 + jacobian * xi
        for node in range(NQ):
            # Independent polynomial-antiderivative oracle for int ell_n(r) r^3 dr.
            basis = np.polynomial.Polynomial([1.0])
            for other, xother in enumerate(xi):
                if other == node:
                    continue
                basis *= np.polynomial.Polynomial([-xother, 1.0]) / (xi[node] - xother)
            # Integrate in the reference coordinate, retaining the Jacobian explicitly.
            rpoly = np.polynomial.Polynomial([(left + right) / 2.0, jacobian])
            primitive = (basis * rpoly ** 3).integ()
            integral = jacobian * (primitive(1.0) - primitive(-1.0))
            rho = 2.0
            a_factor = a ** 2
            rows.append({
                "element": element,
                "node": node,
                "radius": radii[node],
                "u_local": 2.0 * rho / a_factor * integral,
                "v_local": 6.0 * rho / a_factor * integral,
            })
    return rows


def run_timed(adapter: Path, args: list[str], prefix: Path, memory_limit_gib: float) -> dict:
    time_file = prefix.with_suffix(".time.txt")
    log_file = prefix.with_suffix(".stdout.log")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["OMP_DYNAMIC"] = "FALSE"
    env["OPENBLAS_NUM_THREADS"] = "1"
    limit_bytes = int(memory_limit_gib * 1024 ** 3)

    def impose_memory_limit() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))

    command = ["/usr/bin/time", "-v", "-o", str(time_file), str(adapter), *args]
    started = time.perf_counter()
    with log_file.open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=env, preexec_fn=impose_memory_limit,
                                   check=False, text=True)
    wall_s = time.perf_counter() - started
    time_text = time_file.read_text() if time_file.exists() else ""
    rss_kib = None
    for line in time_text.splitlines():
        if "Maximum resident set size" in line:
            rss_kib = int(line.split(":", 1)[1].strip())
    if completed.returncode != 0:
        raise RuntimeError(f"adapter failed ({completed.returncode}); see {log_file}")
    meta_path = prefix.with_suffix(".meta.json")
    metadata = json.loads(meta_path.read_text())
    if args[0] == "run":
        required = ("linear_residual_relative_l2", "normwise_backward_error_inf",
                    "normwise_backward_error_limit", "linear_solver", "global_dofs",
                    "matrix_nonzeros")
        missing = [key for key in required if key not in metadata]
        if missing:
            raise RuntimeError(f"adapter metadata missing required solve checks: {missing}")
        if metadata.get("status") != "success":
            raise RuntimeError(f"adapter solve status is not success: {metadata.get('status')}")
        for key in ("linear_residual_relative_l2", "normwise_backward_error_inf"):
            if not math.isfinite(float(metadata[key])):
                raise RuntimeError(f"non-finite solve diagnostic {key} in {meta_path}")
        if metadata["normwise_backward_error_inf"] > metadata["normwise_backward_error_limit"]:
            raise RuntimeError(f"normwise backward error exceeds frozen limit in {meta_path}")
    metadata["wall_time_s_measured"] = wall_s
    metadata["peak_rss_kib_measured"] = rss_kib
    metadata["address_space_limit_gib"] = memory_limit_gib
    meta_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def verify_manufactured(adapter: Path, run_dir: Path, memory_limit_gib: float) -> list[dict]:
    fixture_dir = run_dir / "manufactured"
    fixture_dir.mkdir()
    results = []
    for elements, maxstep in ((1, 1.0), (2, 0.5)):
        deck = fixture_dir / f"constant-density-{elements}-element.deck"
        write_synthetic_deck(deck)
        prefix = fixture_dir / f"nq5-{elements}-element"
        meta = run_timed(adapter, ["manufactured", str(deck), str(maxstep), str(prefix)],
                         prefix, memory_limit_gib)
        local_path = prefix.with_suffix(".local-loads.csv")
        with local_path.open(newline="") as stream:
            actual = list(csv.DictReader(stream))
        expected = manufactured_expected(2.0, elements)
        if len(actual) != len(expected):
            raise RuntimeError(f"manufactured local row count mismatch for {elements} elements")
        local_errors = []
        expected_u: dict[int, float] = {}
        expected_v: dict[int, float] = {}
        for row, oracle in zip(actual, expected):
            if int(row["element"]) != oracle["element"] or int(row["node"]) != oracle["node"]:
                raise RuntimeError("manufactured local row order/index mismatch")
            uidx, vidx = int(row["u_global"]), int(row["v_global"])
            expected_u[uidx] = expected_u.get(uidx, 0.0) + oracle["u_local"]
            expected_v[vidx] = expected_v.get(vidx, 0.0) + oracle["v_local"]
            local_errors.append(max(
                abs(float(row["u_local"]) - oracle["u_local"]),
                abs(float(row["v_local"]) - oracle["v_local"]),
            ))
        global_errors = []
        for row in actual:
            uidx, vidx = int(row["u_global"]), int(row["v_global"])
            global_errors.extend((
                abs(float(row["u_global_rhs"]) - expected_u[uidx]),
                abs(float(row["v_global_rhs"]) - expected_v[vidx]),
            ))
        sum_u = sum(expected_u.values())
        sum_v = sum(expected_v.values())
        if max(local_errors, default=0.0) > 2e-12 or max(global_errors, default=0.0) > 2e-12:
            raise RuntimeError("manufactured local/global RHS differs from analytic polynomial oracle")
        if not math.isclose(sum_u, 4.0, rel_tol=2e-13, abs_tol=2e-13):
            raise RuntimeError(f"manufactured U sum is {sum_u}, expected 4")
        if not math.isclose(sum_v, 12.0, rel_tol=2e-13, abs_tol=2e-13):
            raise RuntimeError(f"manufactured V sum is {sum_v}, expected 12")
        if not math.isclose(sum_v / sum_u, 3.0, rel_tol=2e-13, abs_tol=2e-13):
            raise RuntimeError("manufactured V/U ratio is not 3")
        if elements == 2:
            left_end = next(row for row in actual if int(row["element"]) == 0 and int(row["node"]) == NQ - 1)
            right_start = next(row for row in actual if int(row["element"]) == 1 and int(row["node"]) == 0)
            if int(left_end["u_global"]) != int(right_start["u_global"]):
                raise RuntimeError("two-element manufactured fixture did not share the U interface DOF")
            if int(left_end["v_global"]) != int(right_start["v_global"]):
                raise RuntimeError("two-element manufactured fixture did not share the V interface DOF")
        result = {
            "elements": elements,
            "nq": NQ,
            "analytic_u_sum": sum_u,
            "analytic_v_sum": sum_v,
            "v_over_u": sum_v / sum_u,
            "max_local_absolute_error": max(local_errors, default=0.0),
            "max_global_assembly_absolute_error": max(global_errors, default=0.0),
            "adapter_meta": meta,
            "deck_sha256": sha256(deck),
        }
        (fixture_dir / f"nq5-{elements}-element.check.json").write_text(json.dumps(result, indent=2) + "\n")
        results.append(result)
    return results


def read_radial(path: Path) -> dict:
    rows = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            values = [float(row[key]) for key in ("radius_m", "U_scaled_m_real", "U_scaled_m_imag", "V_scaled_m_real", "V_scaled_m_imag")]
            if not all(math.isfinite(value) for value in values):
                raise RuntimeError(f"non-finite radial output in {path}")
            rows.append({
                "r_over_a": values[0] / 6371000.0,
                "U": complex(values[1], values[2]),
                "V": complex(values[3], values[4]),
            })
    return {"r": np.array([row["r_over_a"] for row in rows]),
            "U": np.array([row["U"] for row in rows]),
            "V": np.array([row["V"] for row in rows])}


def generate_figures(run_dir: Path, summaries: list[dict], n2_paths: dict[str, Path]) -> None:
    profiles = [name for _, name in PROFILES]
    radial = {}
    max_amp_by_period = {period: 0.0 for period in PERIODS_H}
    for record in summaries:
        if record["mode"] != "physical":
            continue
        series = read_radial(Path(record["radial_csv"]))
        key = (record["resolution"], record["period_hours"])
        radial[key] = series
        period = record["period_hours"]
        max_amp_by_period[period] = max(max_amp_by_period[period],
                                        float(np.max(np.abs(series["U"]))),
                                        float(np.max(np.abs(series["V"]))))
    if not all(max_amp_by_period.values()):
        raise RuntimeError("no radial values available for plotting")

    fig, axes = plt.subplots(3, 4, figsize=(13, 8), sharey=True)
    colors = {12.0: "#0072B2", 48.0: "#D55E00", 168.0: "#009E73", 1680.0: "#CC79A7"}
    for row, resolution in enumerate(profiles):
        for col, period in enumerate(PERIODS_H):
            ax = axes[row, col]
            values = radial[(resolution, period)]
            scale = max_amp_by_period[period]
            ax.plot(values["U"].real / scale, values["r"], color=colors[period], lw=1.0, label="U")
            ax.plot(values["V"].real / scale, values["r"], color=colors[period], lw=1.0, ls="--", label="V")
            meta = next(item for item in summaries if item.get("mode") == "physical" and item["resolution"] == resolution and item["period_hours"] == period)
            ax.axhline(meta["icb_radius_m"] / 6371000.0, color="0.3", lw=0.7, ls=":")
            ax.axhline(meta["cmb_radius_m"] / 6371000.0, color="0.3", lw=0.7, ls=":")
            ax.axvline(0.0, color="0.5", lw=0.6)
            ax.set_title(f"{period:g} h")
            if col == 0:
                ax.set_ylabel(f"{resolution}\nr/a")
            if row == 2:
                ax.set_xlabel("real(U,V) / period peak")
            ax.grid(True, alpha=0.2)
            ax.set_ylim(0, 1)
            ax.set_xlim(-1.05, 1.05)
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("Figure 7 degree-2 tidal response shapes (each period scaled to its maximum across meshes)\nUnit c20; arbitrary external amplitude/sign")
    fig.tight_layout()
    fig.savefig(run_dir / "fig7-radial-response.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.8, 6.2))
    path = n2_paths["fine"]
    data = np.genfromtxt(path, delimiter=",", names=True)
    ax.plot(data["n_squared_s_minus2"], data["radius_m"] / 6371000.0, color="#009E73", lw=1.3)
    fine_meta = next(item for item in summaries if item.get("mode") == "physical" and item["resolution"] == "fine")
    ax.axhline(fine_meta["icb_radius_m"] / 6371000.0, color="0.3", lw=0.8, ls=":", label="ICB")
    ax.axhline(fine_meta["cmb_radius_m"] / 6371000.0, color="0.3", lw=0.8, ls="--", label="CMB")
    ax.set_xlabel(r"$N^2$ (s$^{-2}$)")
    ax.set_ylabel("r/a")
    ax.set_title("PREM outer-core buoyancy profile")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(run_dir / "fig7-n2-profile.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adapter", type=Path, required=True)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA,
                        help="tracked paper-suite criteria JSON (default: cases/paper-suite-protocol.json)")
    parser.add_argument("--paper-text", type=Path, default=None,
                        help="optional local paper text to hash for provenance; not needed to run")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-id", default=time.strftime("run-%Y%m%d-%H%M%S"))
    parser.add_argument("--pilot-only", action="store_true", help="run manufactured gates plus one selected-resolution pilot")
    parser.add_argument("--pilot-resolution", choices=tuple(name for _, name in PROFILES),
                        default="coarse", help="resolution for --pilot-only (default: coarse)")
    parser.add_argument("--pilot-period-hours", type=float, default=168.0)
    parser.add_argument("--memory-limit-gib", type=float, default=6.0,
                        help="per-process address-space cap; defaults to 6 GiB")
    args = parser.parse_args()
    adapter = args.adapter.resolve()
    model = args.model.resolve()
    if not adapter.is_file() or not os.access(adapter, os.X_OK):
        parser.error(f"adapter is not an executable file: {adapter}")
    if not model.is_file():
        parser.error(f"model not found: {model}")
    criteria_path = args.criteria.resolve()
    try:
        validate_criteria(criteria_path)
    except ValueError as exc:
        parser.error(str(exc))
    paper_text_path = args.paper_text.resolve() if args.paper_text is not None else None
    if paper_text_path is not None and not paper_text_path.is_file():
        parser.error(f"paper text not found: {paper_text_path}")
    if args.memory_limit_gib <= 0 or args.memory_limit_gib > 6.0:
        parser.error("memory limit must be in (0, 6] GiB to preserve the aggregate budget")
    run_dir = args.output_root.resolve() / args.run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    preflight = {
        "date": time.strftime("%Y-%m-%d"),
        "cpu_count_logical": os.cpu_count(),
        "available_memory": available_memory(),
        "aggregate_compile_jobs": 1,
        "runtime_omp_threads": 2,
        "runtime_jobs": 1,
        "memory_limit_gib_per_adapter": args.memory_limit_gib,
        "model_path": str(model),
        "model_sha256": sha256(model),
        "adapter_path": str(adapter),
        "adapter_sha256": sha256(adapter),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "criteria_path": str(criteria_path),
        "criteria_sha256": sha256(criteria_path),
        "paper_text_path": str(paper_text_path) if paper_text_path is not None else None,
        "paper_text_sha256": sha256(paper_text_path) if paper_text_path is not None else None,
        "paper_reference_doi": "10.1093/gji/ggag211",
        "dspec_source_hashes": {},
        "profiles": {name: estimate_elements(model, step) for step, name in PROFILES},
        "periods_h": list(PERIODS_H),
        "pilot_resolution": args.pilot_resolution if args.pilot_only else None,
        "paper_rounded_bv_ratios_context_only": list(PAPER_ROUNDED_BV_RATIOS),
        "physical_settings": {
            "degree": 2, "nq": NQ, "gravity": "full self-gravity",
            "attenuation": False, "ocean": False,
            "external_potential": "cbar_20*(r/a)^2*Y_20",
            "cbar_20_dimensionless": 1.0, "force_sign": 1,
            "amplitude_claim": "normalized response shapes only; paper does not specify c20 amplitude or sign",
        },
    }
    for relative in (
        "DSpecM1D/src/SEM/SEMConstructor.h",
        "DSpecM1D/src/SEM/SEMMatrices.h",
        "DSpecM1D/src/SEM/SEMLG.h",
        "DSpecM1D/src/model_info/ModelInput.h",
        "DSpecM1D/src/model_info/MeshModel.h",
        "DSpecM1D/src/model_info/earthmesh/RadialMeshDefinition.h",
    ):
        path = ROOT / relative
        preflight["dspec_source_hashes"][relative] = sha256(path)
    (run_dir / "preflight.json").write_text(json.dumps(preflight, indent=2) + "\n")

    manufactured = verify_manufactured(adapter, run_dir, args.memory_limit_gib)
    cases = []
    n2_paths: dict[str, Path] = {}
    run_plan = [(step, name, args.pilot_period_hours)
                for step, name in PROFILES if name == args.pilot_resolution] if args.pilot_only else [
        (step, resolution, period)
        for step, resolution in PROFILES
        for period in PERIODS_H
    ]
    for step, resolution, period in run_plan:
        case_dir = run_dir / resolution / f"period-{period:g}h"
        case_dir.mkdir(parents=True, exist_ok=False)
        prefix = case_dir / "fig7-tide"
        started = time.perf_counter()
        meta = run_timed(adapter, ["run", str(model), str(step), str(period), str(prefix)],
                         prefix, args.memory_limit_gib)
        elapsed = time.perf_counter() - started
        radial_path = prefix.with_suffix(".radial.csv")
        n2_path = prefix.with_suffix(".n2.csv")
        if not n2_paths.get(resolution):
            n2_paths[resolution] = n2_path
        n2 = np.genfromtxt(n2_path, delimiter=",", names=True)
        min_n2 = float(np.min(n2["n_squared_s_minus2"]))
        # Match the paper's/ex6 convention: the maximum BV period comes from
        # the most negative N^2, via sqrt(-min(N^2)), not max positive N^2.
        bv_period_h = 2.0 * math.pi / math.sqrt(-min_n2) / 3600.0 if min_n2 < 0 else None
        item = {
            "mode": "physical",
            "resolution": resolution,
            "maxstep_fraction_of_a": step,
            "period_hours": period,
            "period_over_discrete_max_bv_period": period / bv_period_h if bv_period_h else None,
            "discrete_max_bv_period_hours": bv_period_h,
            "minimum_n_squared_s_minus2": min_n2,
            "mesh_elements": meta["mesh_elements"],
            "local_gll_nodes": meta["local_gll_nodes"],
            "global_dofs": meta["global_dofs"],
            "matrix_nonzeros": meta["matrix_nonzeros"],
            "linear_residual_relative_l2": meta["linear_residual_relative_l2"],
            "normwise_backward_error_inf": meta["normwise_backward_error_inf"],
            "normwise_backward_error_limit": meta["normwise_backward_error_limit"],
            "icb_radius_m": meta["icb_radius_m"],
            "cmb_radius_m": meta["cmb_radius_m"],
            "peak_rss_kib": meta["peak_rss_kib_measured"],
            "wall_time_s": elapsed,
            "adapter_meta_path": str(prefix.with_suffix(".meta.json")),
            "radial_csv": str(radial_path),
            "n2_csv": str(n2_path),
            "icb_radius_m": meta["icb_radius_m"],
            "cmb_radius_m": meta["cmb_radius_m"],
            "max_abs_U_scaled_m": meta["max_abs_U_scaled_m"],
            "max_abs_V_scaled_m": meta["max_abs_V_scaled_m"],
        }
        cases.append(item)
        print(f"{resolution} {period:g}h: NE={meta['mesh_elements']} dofs={meta['global_dofs']} residual={meta['linear_residual_relative_l2']:.3e} RSS={meta['peak_rss_kib_measured']} KiB wall={elapsed:.1f}s", flush=True)

    summaries = manufactured + cases
    (run_dir / "summary.json").write_text(json.dumps({
        "schema": 1,
        "run_dir": str(run_dir),
        "preflight": preflight,
        "manufactured_checks": manufactured,
        "physical_cases": cases,
        "interpretation": "Fig. 7 resolution/stability diagnostic only; no YSpecpp accuracy claim and no scalar accuracy gate. The maximum BV period follows the paper/ex6 convention 2*pi/sqrt(-min(N^2)) over the outer core. Unit cbar_20 and sign=+1 produce normalized shapes because the paper does not give the Sun/Moon amplitude or sign.",
    }, indent=2) + "\n")
    with (run_dir / "summary.csv").open("w", newline="") as stream:
        fields = ["resolution", "maxstep_fraction_of_a", "period_hours", "period_over_discrete_max_bv_period", "discrete_max_bv_period_hours", "minimum_n_squared_s_minus2", "mesh_elements", "global_dofs", "matrix_nonzeros", "linear_residual_relative_l2", "normwise_backward_error_inf", "normwise_backward_error_limit", "peak_rss_kib", "wall_time_s", "max_abs_U_scaled_m", "max_abs_V_scaled_m"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: item.get(field) for field in fields} for item in cases)
    if not args.pilot_only:
        generate_figures(run_dir, cases, n2_paths)
    print(run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
