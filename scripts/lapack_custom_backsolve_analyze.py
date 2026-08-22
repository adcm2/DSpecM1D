#!/usr/bin/env python3
"""Validate and summarize the focused post-zgbtrf solve experiment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value: float) -> bool:
    return math.isfinite(value)


def read_serial_profile(path: Path, label: str) -> dict:
    fields = path.read_text().strip().split("\t", 6)
    if len(fields) != 7 or fields[0] != "PROFILE":
        raise RuntimeError(f"invalid serial profile: {path}")
    profile = json.loads(fields[6])
    timings = profile["timings_seconds"]["all"]
    result = {
        "backend": label,
        "wall_seconds": float(profile["wall_seconds"]),
        "factorization_worker_seconds": float(timings["factorization"]),
        "solve_worker_seconds": float(timings["solve"]),
        "output_norm": float(fields[5]),
        "counts": profile["counts"],
    }
    if any(not finite(result[key]) or result[key] < 0.0 for key in (
            "wall_seconds", "factorization_worker_seconds",
            "solve_worker_seconds", "output_norm")):
        raise RuntimeError(f"non-finite serial profile: {path}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--check", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--serial-zgbtrs", type=Path, required=True)
    parser.add_argument("--serial-custom", type=Path, required=True)
    parser.add_argument("--serial-eigen", type=Path, required=True)
    args = parser.parse_args()
    raw_lines = [line.rstrip("\n") for line in args.raw.read_text().splitlines()]
    values: dict[str, str] = {}
    timings: list[dict] = []
    benchmark_differences: list[dict] = []
    for line in raw_lines:
        fields = line.split("\t")
        if fields[0] in {"validation_cases", "validation_passed",
                         "pivot_exercising_cases", "max_normalized_difference",
                         "max_normalized_residual"}:
            values[fields[0]] = fields[1]
        elif fields[0].isdigit() and len(fields) == 9:
            timings.append({
                "threads": int(fields[0]), "case": fields[1],
                "backend": fields[2], "seconds": float(fields[3]),
                "n": int(fields[4]), "kl": int(fields[5]),
                "ku": int(fields[6]), "nrhs": int(fields[7]),
                "systems": int(fields[8])})
        elif fields[0] == "benchmark_max_normalized_difference":
            benchmark_differences.append({"case": fields[1],
                                          "threads": int(fields[2]),
                                          "value": float(fields[3])})
    expected = int(values["validation_cases"])
    passed = int(values["validation_passed"])
    if passed != expected or expected < 1:
        raise RuntimeError("validation cases did not all pass")
    if int(values["pivot_exercising_cases"]) < 1:
        raise RuntimeError("no validation case exercised a row pivot")
    max_difference = float(values["max_normalized_difference"])
    max_residual = float(values["max_normalized_residual"])
    if not finite(max_difference) or not finite(max_residual) or max_difference >= 1e-12 or max_residual >= 1e-12:
        raise RuntimeError("validation discrepancy or residual is too large")
    if len(timings) != 24 or len(benchmark_differences) != 12:
        raise RuntimeError("unexpected timing row count")
    if any(not finite(row["seconds"]) or row["seconds"] <= 0.0
           for row in timings):
        raise RuntimeError("non-finite timing")
    by_key = {(row["case"], row["threads"], row["backend"]): row
              for row in timings}
    speedups = []
    scaling = []
    for case in sorted({row["case"] for row in timings}):
        for threads in (1, 10, 32, 50):
            lapack = by_key[(case, threads, "lapacke_zgbtrs")]["seconds"]
            custom = by_key[(case, threads, "eigen_custom")]["seconds"]
            speedups.append({"case": case, "threads": threads,
                             "lapack_over_custom": lapack / custom})
        lapack_scale = (by_key[(case, 1, "lapacke_zgbtrs")]["seconds"] /
                        by_key[(case, 32, "lapacke_zgbtrs")]["seconds"])
        custom_scale = (by_key[(case, 1, "eigen_custom")]["seconds"] /
                        by_key[(case, 32, "eigen_custom")]["seconds"])
        scaling.append({"case": case, "lapack_1_to_32": lapack_scale,
                        "custom_1_to_32": custom_scale})
    pilot_records = []
    for line in args.pilot.read_text().splitlines():
        if line.startswith("fmax_mHz"):
            continue
        fields = line.split("\t", 6)
        profile = json.loads(fields[6])
        pilot_records.append({"backend": fields[1], "wall_seconds": profile["wall_seconds"],
                              "solve_worker_seconds": profile["timings_seconds"]["all"]["solve"],
                              "norm": float(fields[5])})
    if len(pilot_records) != 2:
        raise RuntimeError("expected one pilot record per backend")
    check = args.check.read_text().strip().splitlines()
    if not any(line.startswith("CHECK\t") for line in check):
        raise RuntimeError("missing production numerical comparison")
    check_fields = next(line.split("\t") for line in check
                        if line.startswith("CHECK\t"))
    full_difference = float(check_fields[5])
    if not finite(full_difference) or full_difference > 1e-9:
        raise RuntimeError("production discrepancy exceeds tolerance")

    serial = {
        "zgbtrs": read_serial_profile(args.serial_zgbtrs, "zgbtrs"),
        "custom": read_serial_profile(args.serial_custom, "custom"),
        "eigen": read_serial_profile(args.serial_eigen, "eigen"),
    }
    count_keys = ("frequency_systems", "solves", "rhs", "dimension_min",
                  "dimension_max", "kl_min", "kl_max", "ku_min", "ku_max")
    for key in count_keys:
        if serial["zgbtrs"]["counts"][key] != serial["custom"]["counts"][key]:
            raise RuntimeError(f"serial LAPACK workload mismatch: {key}")
    zgbtrs_norm = serial["zgbtrs"]["output_norm"]
    custom_norm = serial["custom"]["output_norm"]
    serial["custom_over_zgbtrs_wall"] = (
        serial["custom"]["wall_seconds"] / serial["zgbtrs"]["wall_seconds"])
    serial["custom_over_zgbtrs_solve"] = (
        serial["custom"]["solve_worker_seconds"] /
        serial["zgbtrs"]["solve_worker_seconds"])
    serial["custom_over_zgbtrs_factorization"] = (
        serial["custom"]["factorization_worker_seconds"] /
        serial["zgbtrs"]["factorization_worker_seconds"])
    serial["output_norm_absolute_difference"] = abs(custom_norm - zgbtrs_norm)
    serial["output_norm_normalized_difference"] = (
        abs(custom_norm - zgbtrs_norm) /
        (1.0 + abs(custom_norm) + abs(zgbtrs_norm)))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "lapack_custom_backsolve_summary.json"
    summary = {
        "validation": {"cases": expected, "passed": passed,
                        "pivot_exercising_cases": int(values["pivot_exercising_cases"]),
                        "max_normalized_difference": max_difference,
                        "max_normalized_residual": max_residual},
        "timings": timings, "speedups": speedups, "scaling_1_to_32": scaling,
        "benchmark_max_normalized_differences": benchmark_differences,
        "production_pilot": {"records": pilot_records,
                              "full_output_max_abs_difference": full_difference},
        "serial_f20_gate": serial,
        "policy": "one warm-up and one measured run per timing point; no threshold",
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    metadata_path = args.output_dir / "lapack_custom_backsolve_metadata.json"
    profiling_executable = Path(
        "/tmp/dspecm1d-custom-backsolve-gcc15/bin/stage_profile_solver")
    final_check_executable = Path(
        "/tmp/dspecm1d-custom-backsolve-gcc15/bin/final_solver_benchmark")
    pilot_metadata_path = args.output_dir / "custom_f35_t32_pilot_metadata.json"
    solver_path = Path("DSpecM1D/src/LapackBandSolver.h")
    analyzer_path = Path(__file__)
    benchmark_path = Path("benchmarks/lapack_custom_backsolve_probe.cpp")
    loaded_runtimes = {
        "lapacke": Path("/usr/lib/x86_64-linux-gnu/liblapacke.so.3.11.0"),
        "lapack": Path(
            "/usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3"),
        "openblas": Path(
            "/usr/lib/x86_64-linux-gnu/openblas-pthread/"
            "libopenblasp-r0.3.21.so"),
        "blas": Path(
            "/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3"),
    }
    metadata = {
        "campaign": "focused post-zgbtrf Eigen/C++ backsolve experiment",
        "algorithm": "LAPACK 3.11.0 ZGBTRS TRANS=N replay of ZSWAP/ZGERU/ ZTBSV",
        "compiler": "/opt/gcc-15.2.0/bin/g++",
        "threads": [1, 10, 32, 50], "environment": {
            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1"},
        "raw_artifacts": {str(args.raw): sha256(args.raw),
                          str(args.pilot): sha256(args.pilot),
                          str(args.check): sha256(args.check),
                          str(args.serial_zgbtrs): sha256(args.serial_zgbtrs),
                          str(args.serial_custom): sha256(args.serial_custom),
                          str(args.serial_eigen): sha256(args.serial_eigen)},
        "retained_baseline": {
            "path": "results/final_benchmark/final_production_timings.tsv",
            "sha256": sha256(Path("results/final_benchmark/final_production_timings.tsv")),
            "comparison": "primary_f35_t32 LAPACK median 22.496827324 s",
        },
        "executable": {"path": str(args.executable.resolve()),
                       "sha256": sha256(args.executable)},
        "source": {"path": "benchmarks/lapack_custom_backsolve_probe.cpp",
                   "sha256": sha256(benchmark_path)},
        "executions": {
            "serial_f20_gate": {
                "policy": "one untimed warm-up and one measured run per backend",
                "parameter_sha256": "dad298c571ce0f3a14f6c130b61b4dbf0665cefea718bc229e6d9c7fe616a569",
                "environment": {
                    "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
                    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                    "BLIS_NUM_THREADS": "1",
                },
                "zgbtrs": {
                    "revision": "8f7a5d3c03b14425e5644a5ffd2faaeb996a3129",
                    "executable": "/tmp/dspecm1d-zgbtrs-serial-gcc15/bin/stage_profile_solver",
                    "executable_sha256": "17ba92135bec221c5b43a7ca858c0c41b4802aab1e74475978748f094a6c0c43",
                    "command": "env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ./bin/stage_profile_solver ./data/params/custom_backsolve_serial_f20.txt lapack 1 1",
                },
                "custom": {
                    "executable": "/tmp/dspecm1d-custom-backsolve-gcc15/bin/stage_profile_solver",
                    "executable_sha256": "bceb4f3607bab9a9231930fb8423343cae830299e0748da74bfa72770922fb3b",
                    "command": "env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ./bin/stage_profile_solver ./data/params/custom_backsolve_serial_f20.txt lapack 1 1",
                },
                "eigen_control": {
                    "executable": "/tmp/dspecm1d-custom-backsolve-gcc15/bin/stage_profile_solver",
                    "executable_sha256": "bceb4f3607bab9a9231930fb8423343cae830299e0748da74bfa72770922fb3b",
                    "command": "env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ./bin/stage_profile_solver ./data/params/custom_backsolve_serial_f20.txt eigen 1 1",
                },
            },
            "profiling_pilot": {
                "command": "env OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 python3 scripts/stage_profile.py --eigen-build /tmp/dspecm1d-stage-g-off-verify --lapack-build /tmp/dspecm1d-custom-backsolve-gcc15 --output-dir results/lapack_custom_backsolve --fmax 35 --lmax 750 --threads 32 --repetitions 1 --tag custom_f35_t32_pilot",
                "executable": str(profiling_executable),
                "executable_sha256": sha256(profiling_executable),
                "raw_path": str(args.pilot.resolve()),
                "raw_sha256": sha256(args.pilot),
            },
            "full_output_check": {
                "command": "env OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 /tmp/dspecm1d-custom-backsolve-gcc15/bin/final_solver_benchmark /tmp/dspecm1d-custom-backsolve-gcc15/data/params/custom_f35_t32_pilot_35.txt lapack 1 1e-9 > results/lapack_custom_backsolve/custom_f35_t32_final_check.txt 2>&1",
                "executable": str(final_check_executable),
                "executable_sha256": sha256(final_check_executable),
                "output_path": str(args.check.resolve()),
                "output_sha256": sha256(args.check),
            },
        },
        "retained_hashes": {
            "pilot_metadata": sha256(pilot_metadata_path),
            "production_lapack_band_solver": sha256(solver_path),
            "analyzer": sha256(analyzer_path),
            "benchmark_source": sha256(benchmark_path),
            "loaded_runtime_libraries": {
                name: {"path": str(path), "sha256": sha256(path)}
                for name, path in loaded_runtimes.items()
            },
            "zgbtrs_source": {
                "version": "3.11.0",
                "path": "Reference-LAPACK/lapack v3.11.0/SRC/zgbtrs.f",
                "sha256": "96e35c7c9ceb47aa7cc51bd21f16365a5f713b44ceb635802fa5eb1ac58ae550",
            },
        },
        "summary": {"path": str(summary_path.resolve()),
                    "sha256": sha256(summary_path)},
        "validation_gate": "passed; universal experimental production pilot run",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"summary": str(summary_path),
                      "metadata": str(metadata_path),
                      "max_full_output_difference": full_difference}, indent=2))


if __name__ == "__main__":
    main()
