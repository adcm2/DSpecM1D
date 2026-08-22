#!/usr/bin/env python3
"""Normalize the focused LAPACK scaling diagnostic artifacts.

The expensive runs are intentionally performed outside this parser.  It
records the exact linked libraries, OpenBLAS introspection, source hashes,
category profiles, and the two focused band-call probes in machine-readable
metadata and summary files.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(command: list[str]) -> str:
    return subprocess.check_output(command, text=True).strip()


def profile(path: Path, threads: int) -> dict[str, Any]:
    lines = [line for line in path.read_text().splitlines()
             if line.startswith("PROFILE\t")]
    if len(lines) != 1:
        raise RuntimeError(f"expected one PROFILE line in {path}")
    fields = lines[0].split("\t", 6)
    if len(fields) != 7 or fields[1] != "lapack":
        raise RuntimeError(f"invalid profile record in {path}")
    data = json.loads(fields[6])
    counts = data["counts"]
    if counts["lapack_factorize"] != counts["solves"]:
        raise RuntimeError(f"factor/solve count mismatch in {path}")
    if not (data["wall_seconds"] > 0.0 and all(
            value >= 0.0 for value in data["timings_seconds"]["all"].values())):
        raise RuntimeError(f"non-finite or non-positive timing in {path}")
    all_times = data["timings_seconds"]["all"]
    worker = data["accounting"]["workers"]["all"]
    return {
        "threads": threads,
        "wall_seconds": data["wall_seconds"],
        "worker_total_seconds": worker["total"],
        "worker_categorized_seconds": worker["categorized"],
        "worker_unclassified_seconds": worker["unclassified"],
        "base_operator_worker_seconds": all_times["base_operator_preparation"],
        "dynamic_matrix_worker_seconds": all_times["frequency_matrix_construction"],
        "factorization_worker_seconds": all_times["factorization"],
        "solve_worker_seconds": all_times["solve"],
        "projection_worker_seconds": all_times["source_receiver_projection"],
        "unclassified_worker_seconds": all_times["unclassified"],
        "spheroidal_factorization_worker_seconds":
            data["timings_seconds"]["spheroidal"]["factorization"],
        "spheroidal_solve_worker_seconds":
            data["timings_seconds"]["spheroidal"]["solve"],
        "counts": counts,
        "source": str(path),
    }


def openblas_info() -> dict[str, Any]:
    path = "/lib/x86_64-linux-gnu/libopenblas.so.0"
    library = ctypes.CDLL(path)
    library.openblas_get_parallel.restype = ctypes.c_int
    library.openblas_get_num_threads.restype = ctypes.c_int
    library.openblas_get_config.restype = ctypes.c_char_p
    parallel = library.openblas_get_parallel()
    return {
        "path": path,
        "parallel_code": parallel,
        "parallel_kind": {0: "sequential", 1: "pthread", 2: "openmp"}.get(
            parallel, "unknown"),
        "runtime_threads": library.openblas_get_num_threads(),
        "config": library.openblas_get_config().decode(),
        "environment": {key: os.environ.get(key, "") for key in (
            "OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
    }


def resolved_library_targets() -> dict[str, str]:
    names = {
        "openblas": "/lib/x86_64-linux-gnu/libopenblas.so.0",
        "lapacke": "/lib/x86_64-linux-gnu/liblapacke.so.3",
        "lapack": "/lib/x86_64-linux-gnu/liblapack.so.3",
        "blas": "/lib/x86_64-linux-gnu/libblas.so.3",
    }
    return {key: run(["readlink", "-f", value])
            for key, value in names.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--profile-t10", type=Path, required=True)
    parser.add_argument("--profile-t32", type=Path, required=True)
    parser.add_argument("--thread-probe", type=Path, required=True)
    parser.add_argument("--band-probe", type=Path, required=True)
    parser.add_argument("--allocation-probe", type=Path, required=True)
    parser.add_argument("--band-executable", type=Path, required=True)
    parser.add_argument("--allocation-executable", type=Path, required=True)
    parser.add_argument("--final-executable", type=Path, required=True)
    parser.add_argument("--profile-executable", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    thread_lines = args.thread_probe.read_text().splitlines()
    thread_header = ["sample_seconds", "process_nlwp", "OMP_NUM_THREADS",
                     "OPENBLAS_NUM_THREADS", "OMP_DYNAMIC", "backend"]
    if not thread_lines or thread_lines[0].split("\t") != thread_header:
        raise RuntimeError("unexpected effective-thread probe header")
    thread_samples = []
    for line in thread_lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(thread_header):
            raise RuntimeError("malformed effective-thread probe row")
        sample = dict(zip(thread_header, fields))
        if (int(sample["process_nlwp"]) < 1 or
                sample["OMP_NUM_THREADS"] != "32" or
                sample["OPENBLAS_NUM_THREADS"] != "1" or
                sample["OMP_DYNAMIC"] != "FALSE" or
                sample["backend"] != "lapack"):
            raise RuntimeError("unexpected effective-thread probe value")
        thread_samples.append(sample)
    observed_nlwp = sorted({int(sample["process_nlwp"])
                            for sample in thread_samples})
    if observed_nlwp != [1, 32]:
        raise RuntimeError(f"unexpected process thread counts: {observed_nlwp}")

    profiles = [profile(args.profile_t10, 10), profile(args.profile_t32, 32)]
    profile_header = [
        "threads", "wall_seconds", "worker_total_seconds",
        "worker_categorized_seconds", "worker_unclassified_seconds",
        "base_operator_worker_seconds", "dynamic_matrix_worker_seconds",
        "factorization_worker_seconds", "solve_worker_seconds",
        "projection_worker_seconds",
        "unclassified_worker_seconds", "spheroidal_factorization_worker_seconds",
        "spheroidal_solve_worker_seconds",
    ]
    profiles_path = args.output_dir / "lapack_scaling_profiles.tsv"
    with profiles_path.open("w") as target:
        target.write("\t".join(profile_header) + "\n")
        for row in profiles:
            target.write("\t".join(str(row[key]) for key in profile_header) + "\n")

    expected_headers = {
        args.band_probe: ["threads", "metric", "seconds", "n", "kl", "ku",
                          "nrhs", "systems"],
        args.allocation_probe: ["threads", "metric", "seconds", "n", "kl",
                                "ku", "nrhs", "systems", "repeats"],
    }
    expected_metrics = {
        args.band_probe: {"zgbtrf_only", "zgbtrs_only_allocate",
                          "zgbtrs_only_reused", "zgbtrf_plus_zgbtrs"},
        args.allocation_probe: {"zgbtrs_allocate", "zgbtrs_reused"},
    }
    micro_header: list[str] = ["probe"]
    micro_rows: list[dict[str, str]] = []
    for path in (args.band_probe, args.allocation_probe):
        lines = path.read_text().splitlines()
        if not lines:
            raise RuntimeError(f"empty probe output: {path}")
        header = lines[0].split("\t")
        if header != expected_headers[path]:
            raise RuntimeError(f"unexpected probe header in {path}: {header}")
        for name in header:
            if name not in micro_header:
                micro_header.append(name)
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) != len(header):
                raise RuntimeError(f"malformed probe row: {path}")
            if fields[1] not in expected_metrics[path]:
                raise RuntimeError(f"unexpected probe metric in {path}: {fields[1]}")
            try:
                numbers = [float(fields[index]) for index in (0, 2, 3, 4, 5, 6, 7)]
                if any(not math.isfinite(value) or value <= 0.0
                       for value in numbers):
                    raise ValueError
                if len(fields) == 9:
                    repeats = float(fields[8])
                    if not math.isfinite(repeats) or repeats <= 0.0:
                        raise ValueError
            except ValueError as error:
                raise RuntimeError(f"non-positive/non-finite probe row: {path}") from error
            micro_rows.append(dict(zip(["probe", *header], [path.stem, *fields])))
    microbenchmarks_path = args.output_dir / "lapack_scaling_microbenchmarks.tsv"
    with microbenchmarks_path.open("w") as target:
        target.write("\t".join(micro_header) + "\n")
        for row in micro_rows:
            target.write("\t".join(row.get(key, "") for key in micro_header) + "\n")

    source_paths = [
        "benchmarks/lapack_band_concurrency_probe.cpp",
        "benchmarks/lapack_solve_allocation_probe.cpp",
        "benchmarks/stage_profile_solver.cpp",
        "scripts/lapack_scaling_diagnostic.py",
        "DSpecM1D/src/LapackBandSolver.h",
    ]
    probe_sources = {
        "lapack_band_concurrency_probe.cpp": args.repo /
        "benchmarks/lapack_band_concurrency_probe.cpp",
        "lapack_solve_allocation_probe.cpp": args.repo /
        "benchmarks/lapack_solve_allocation_probe.cpp",
    }
    # The command-line probe outputs are accepted only when the corresponding
    # executable was rebuilt after its source.  This catches stale raw TSVs.
    # The absolute build paths are recorded below; existence is checked there.
    build_probe_paths = {
        "lapack_band_concurrency_probe.cpp": args.band_executable,
        "lapack_solve_allocation_probe.cpp": args.allocation_executable,
    }
    for name, source in probe_sources.items():
        executable = build_probe_paths[name]
        if not executable.exists() or executable.stat().st_mtime_ns < source.stat().st_mtime_ns:
            raise RuntimeError(f"probe executable is stale or missing: {executable}")
    links = run(["ldd", str(args.final_executable)])
    profile_links = run(["ldd", str(args.profile_executable)])
    t10, t32 = profiles
    summary = {
        "diagnosis": "zgbtrs/solve granularity and concurrent tiny-call overhead",
        "factorization_scaling_ratio_t32_over_t10":
            t32["factorization_worker_seconds"] / t10["factorization_worker_seconds"],
        "solve_scaling_ratio_t32_over_t10":
            t32["solve_worker_seconds"] / t10["solve_worker_seconds"],
        "wall_speedup_t10_over_t32": t10["wall_seconds"] / t32["wall_seconds"],
        "allocation_effect": "negligible in the focused minimum-size solve probe",
        "next_optimization": "reduce tiny zgbtrs call count or batch/reuse solves before affinity tuning",
        "profiles": profiles,
    }
    summary_path = args.output_dir / "lapack_scaling_diagnostic_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    executable_paths = {
        "final": args.final_executable.resolve(),
        "profiling": args.profile_executable.resolve(),
        "band_probe": args.band_executable.resolve(),
        "allocation_probe": args.allocation_executable.resolve(),
    }
    raw_artifact_paths = {
        "profile_t10": args.profile_t10.resolve(),
        "profile_t32": args.profile_t32.resolve(),
        "effective_thread_probe": args.thread_probe.resolve(),
        "band_probe": args.band_probe.resolve(),
        "allocation_probe": args.allocation_probe.resolve(),
    }
    derived_artifact_paths = {
        "profiles": profiles_path.resolve(),
        "microbenchmarks": microbenchmarks_path.resolve(),
        "summary": summary_path.resolve(),
    }
    metadata = {
        "revision": run(["git", "rev-parse", "HEAD"]),
        "git_status_short": run(["git", "status", "--short"]),
        "machine": {
            "hostname": platform.node(),
            "lscpu": run(["lscpu"]),
        },
        "controls": {
            key: os.environ.get(key, "") for key in (
                "OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "BLIS_NUM_THREADS")
        },
        "openblas": openblas_info(),
        "linked_libraries": {
            "final_executable": links,
            "profiling_executable": profile_links,
        },
        "resolved_library_targets": resolved_library_targets(),
        "executables": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in executable_paths.items()
        },
        "raw_artifacts": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in raw_artifact_paths.items()
        },
        "derived_artifacts": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in derived_artifact_paths.items()
        },
        "source_sha256": {path: sha256(args.repo / path) for path in source_paths},
        "probe_source_sha256": {
            "lapack_band_concurrency_probe.cpp": sha256(
                args.repo / "benchmarks/lapack_band_concurrency_probe.cpp"),
            "lapack_solve_allocation_probe.cpp": sha256(
                args.repo / "benchmarks/lapack_solve_allocation_probe.cpp"),
        },
        "probe_build_coherence": {
            name: {
                "source": str(probe_sources[name]),
                "source_mtime_ns": probe_sources[name].stat().st_mtime_ns,
                "executable": str(build_probe_paths[name]),
                "executable_mtime_ns": build_probe_paths[name].stat().st_mtime_ns,
            }
            for name in probe_sources
        },
        "effective_thread_probe": {
            "nominal_omp_threads": 32,
            "openblas_num_threads": 1,
            "observed_process_nlwp": observed_nlwp,
            "samples": thread_samples,
            "interpretation": "32 outer OpenMP workers, no hidden BLAS worker threads",
        },
        "profile_records": profiles,
        "microbenchmark_rows": micro_rows,
    }
    (args.output_dir / "lapack_scaling_diagnostic_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
