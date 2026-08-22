#!/usr/bin/env python3
"""Validate and summarize the final production solver benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        return list(csv.DictReader(source, delimiter="\t"))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def timing_summary(records: list[dict[str, str]], series: str,
                   fmax: float, threads: int) -> dict[str, Any]:
    selected = [record for record in records
                if record["series"] == series
                and float(record["fmax_mHz"]) == fmax
                and int(record["threads"]) == threads]
    grouped = {backend: [float(record["seconds"]) for record in selected
                         if record["backend"] == backend]
               for backend in ("eigen", "lapack")}
    if not grouped["eigen"] or not grouped["lapack"]:
        raise RuntimeError(f"missing timing group {series}, {fmax}, {threads}")
    eigen_median = statistics.median(grouped["eigen"])
    lapack_median = statistics.median(grouped["lapack"])
    return {
        "series": series, "fmax_mHz": fmax, "threads": threads,
        "eigen_repetitions": len(grouped["eigen"]),
        "lapack_repetitions": len(grouped["lapack"]),
        "eigen_min_s": min(grouped["eigen"]),
        "eigen_median_s": eigen_median,
        "eigen_max_s": max(grouped["eigen"]),
        "lapack_min_s": min(grouped["lapack"]),
        "lapack_median_s": lapack_median,
        "lapack_max_s": max(grouped["lapack"]),
        "speedup_eigen_over_lapack": eigen_median / lapack_median,
        "lapack_wall_time_reduction_percent":
            100.0 * (1.0 - lapack_median / eigen_median),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timings", type=Path, required=True)
    parser.add_argument("--comparisons", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tag", default="final_production")
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--build", type=Path, required=True)
    args = parser.parse_args()

    timings = load_tsv(args.timings)
    comparisons = load_tsv(args.comparisons)
    metadata = json.loads(args.metadata.read_text())
    if metadata.get("state") != "complete":
        raise RuntimeError("campaign metadata is not complete")
    if sha256(args.timings) != metadata.get("raw_sha256"):
        raise RuntimeError("raw timing SHA-256 does not match metadata")
    if sha256(args.comparisons) != metadata.get("comparisons_sha256"):
        raise RuntimeError("comparison SHA-256 does not match metadata")

    retained_hashes = {
        "template_sha256": sha256(args.repo / "data/params/ex3.txt"),
        "harness_source_sha256": sha256(
            args.repo / "benchmarks/final_solver_benchmark.cpp"),
        "driver_sha256": sha256(args.repo / "scripts/final_solver_benchmark.py"),
        "executable_sha256": sha256(
            args.build / "bin/final_solver_benchmark"),
    }
    for key, value in retained_hashes.items():
        if value != metadata["request"].get(key):
            raise RuntimeError(f"{key} does not match campaign metadata")
    compile_commands = subprocess.check_output(
        ["ninja", "-C", str(args.build), "-t", "commands",
         "final_solver_benchmark"], text=True)
    if "DSPECM1D_ENABLE_LAPACK_BAND_SOLVER" not in compile_commands:
        raise RuntimeError("retained executable lacks LAPACK capability")
    if "DSPECM1D_ENABLE_PROFILING" in compile_commands:
        raise RuntimeError("retained executable unexpectedly enables profiling")

    requested = metadata["request"]["points"]
    requested_ids = {point["point_id"] for point in requested}
    if set(metadata["completed_points"]) != requested_ids:
        raise RuntimeError("completed-point set differs from request")
    if {record["point_id"] for record in comparisons} != requested_ids:
        raise RuntimeError("comparison-point set differs from request")
    if {record["point_id"] for record in timings} != requested_ids:
        raise RuntimeError("timing-point set differs from request")

    comparisons_by_id = {record["point_id"]: record for record in comparisons}
    for point in requested:
        point_records = [record for record in timings
                         if record["point_id"] == point["point_id"]]
        expected = 2 * int(point["repetitions"])
        if len(point_records) != expected:
            raise RuntimeError(
                f"{point['point_id']} has {len(point_records)}/{expected} timings")
        if point_records[0]["backend"] != point["first_backend"]:
            raise RuntimeError(f"launch order mismatch in {point['point_id']}")
        for backend in ("eigen", "lapack"):
            backend_records = [record for record in point_records
                               if record["backend"] == backend]
            if sorted(int(record["repetition"]) for record in backend_records) != \
                    list(range(1, int(point["repetitions"]) + 1)):
                raise RuntimeError(f"bad repetitions for {point['point_id']} {backend}")
        shapes = {(int(record["rows"]), int(record["cols"]))
                  for record in point_records}
        if len(shapes) != 1:
            raise RuntimeError(f"shape mismatch within {point['point_id']}")
        for record in point_records:
            values = [float(record["seconds"]), float(record["norm"])]
            if not all(math.isfinite(value) and value > 0.0 for value in values):
                raise RuntimeError(f"invalid timing/norm in {point['point_id']}")
            if record["first_backend"] != point["first_backend"]:
                raise RuntimeError(f"backend-order mismatch in {point['point_id']}")
        comparison = comparisons_by_id[point["point_id"]]
        difference = float(comparison["max_abs_difference"])
        tolerance = float(comparison["tolerance"])
        if (not math.isfinite(difference) or difference < 0.0 or
                difference > tolerance or
                tolerance != float(metadata["request"]["tolerance"])):
            raise RuntimeError(f"invalid comparison for {point['point_id']}")
        if (int(comparison["rows"]), int(comparison["cols"])) != next(iter(shapes)):
            raise RuntimeError(f"comparison shape mismatch in {point['point_id']}")
        for backend in ("eigen", "lapack"):
            final_record = [record for record in point_records
                            if record["backend"] == backend][-1]
            if not math.isclose(float(final_record["norm"]),
                                float(comparison[f"{backend}_norm"]),
                                rel_tol=1e-14, abs_tol=0.0):
                raise RuntimeError(
                    f"comparison norm mismatch in {point['point_id']} {backend}")

    primary_fmax = [5.0, 10.0, 20.0, 35.0, 55.0]
    primary = [timing_summary(timings, "primary", fmax, 32)
               for fmax in primary_fmax]
    single = [timing_summary(timings, "single_thread", fmax, 1)
              for fmax in (5.0, 10.0, 20.0)]
    scaling = [timing_summary(
        timings, "primary" if threads == 32 else "scaling", 35.0, threads)
        for threads in (10, 20, 32, 50)]
    numerical = [{
        "point_id": record["point_id"], "series": record["series"],
        "fmax_mHz": float(record["fmax_mHz"]),
        "threads": int(record["threads"]),
        "rows": int(record["rows"]), "cols": int(record["cols"]),
        "max_abs_difference": float(record["max_abs_difference"]),
        "tolerance": float(record["tolerance"]),
    } for record in comparisons]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / f"{args.tag}_summary.json"
    primary_csv = args.output_dir / f"{args.tag}_primary.csv"
    single_csv = args.output_dir / f"{args.tag}_single_thread.csv"
    scaling_csv = args.output_dir / f"{args.tag}_scaling.csv"
    numerical_csv = args.output_dir / f"{args.tag}_numerical.csv"
    summary_path.write_text(json.dumps({
        "validation": {
            "timing_records": len(timings),
            "paired_comparisons": len(comparisons),
            "raw_sha256": sha256(args.timings),
            "comparisons_sha256": sha256(args.comparisons),
            "retained_hashes": retained_hashes,
            "lapack_capability_enabled": True,
            "profiling_enabled": False,
            "maximum_discrepancy": max(
                row["max_abs_difference"] for row in numerical),
        },
        "primary_32_threads": primary,
        "single_thread": single,
        "scaling_35_mHz": scaling,
        "numerical": numerical,
    }, indent=2, sort_keys=True) + "\n")
    write_csv(primary_csv, primary)
    write_csv(single_csv, single)
    write_csv(scaling_csv, scaling)
    write_csv(numerical_csv, numerical)

    import matplotlib.pyplot as plt

    primary_figure = args.output_dir / f"{args.tag}_wall_time.png"
    plt.figure(figsize=(7.0, 4.5))
    plt.plot(primary_fmax, [row["eigen_median_s"] for row in primary], "o-",
             label="Eigen SparseLU")
    plt.plot(primary_fmax, [row["lapack_median_s"] for row in primary], "o-",
             label="LAPACK direct-band LU")
    plt.xlabel(r"$f_{\max}$ [mHz]")
    plt.ylabel("median complete spectra wall time [s]")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(primary_figure, dpi=180)
    plt.close()

    scaling_figure = args.output_dir / f"{args.tag}_thread_scaling.png"
    threads = [row["threads"] for row in scaling]
    plt.figure(figsize=(7.0, 4.5))
    plt.plot(threads, [row["eigen_median_s"] for row in scaling], "o-",
             label="Eigen SparseLU")
    plt.plot(threads, [row["lapack_median_s"] for row in scaling], "o-",
             label="LAPACK direct-band LU")
    plt.xlabel("OpenMP threads")
    plt.ylabel("median complete spectra wall time [s]")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(scaling_figure, dpi=180)
    plt.close()

    print(json.dumps({
        "summary": str(summary_path), "primary_figure": str(primary_figure),
        "scaling_figure": str(scaling_figure), "timing_records": len(timings),
        "paired_comparisons": len(comparisons),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
