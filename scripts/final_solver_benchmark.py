#!/usr/bin/env python3
"""Run or resume the final same-binary Eigen/LAPACK benchmark campaign."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
from typing import Any


PRIMARY_FMAX = [5.0, 10.0, 20.0, 35.0, 55.0]
SINGLE_FMAX = [5.0, 10.0, 20.0]
SCALING_THREADS = [10, 20, 32, 50]
RAW_HEADER = ["point_id", "series", "fmax_mHz", "threads", "backend",
              "repetition", "seconds", "rows", "cols", "norm",
              "first_backend"]
CHECK_HEADER = ["point_id", "series", "fmax_mHz", "threads", "rows", "cols",
                "eigen_norm", "lapack_norm", "max_abs_difference", "tolerance"]


def run_text(command: list[str], *, cwd: Path | None = None) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def replace_parameter(template: Path, destination: Path, fmax: float,
                      lmax: int) -> None:
    lines = template.read_text().splitlines(keepends=True)

    def replace_after(marker: str, value: str) -> None:
        for index, line in enumerate(lines):
            if marker in line:
                for target in range(index + 1, len(lines)):
                    if (lines[target].strip() and
                            not lines[target].lstrip().startswith("#")):
                        ending = "\n" if lines[target].endswith("\n") else ""
                        lines[target] = f"  {value}{ending}"
                        return
        raise RuntimeError(f"parameter marker not found: {marker}")

    replace_after("# fmax", f"{fmax:.12g}")
    replace_after("# f21 filter", f"{max(0.1, 0.9 * fmax):.12g}")
    replace_after("# f22 filter", f"{fmax:.12g}")
    replace_after("# lmax", str(lmax))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(lines))


def cpu_metadata() -> dict[str, str]:
    wanted = {"Model name", "CPU(s)", "Thread(s) per core", "Core(s) per socket",
              "Socket(s)"}
    result: dict[str, str] = {}
    for line in run_text(["lscpu"]).splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            if key.strip() in wanted:
                result[key.strip()] = value.strip()
    return result


def build_metadata(build: Path) -> dict[str, str]:
    wanted = {"CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
              "CMAKE_CXX_FLAGS_RELEASE", "DSPECM1D_ENABLE_LAPACK_BAND_SOLVER"}
    result: dict[str, str] = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(("#", "//")) or ":" not in line or "=" not in line:
            continue
        key_with_type, value = line.split("=", 1)
        key = key_with_type.split(":", 1)[0]
        if key in wanted:
            result[key] = value
    return result


def point_id(series: str, fmax: float, threads: int) -> str:
    return f"{series}_f{fmax:g}_t{threads}"


def requested_points(args: argparse.Namespace) -> list[dict[str, Any]]:
    points: list[dict[str, Any]] = []
    for fmax in args.primary_fmax:
        points.append({"series": "primary", "fmax": fmax,
                       "threads": args.primary_threads,
                       "repetitions": args.primary_repetitions})
    for fmax in args.single_fmax:
        points.append({"series": "single_thread", "fmax": fmax,
                       "threads": 1, "repetitions": args.single_repetitions})
    for threads in args.scaling_threads:
        if threads == args.primary_threads and args.scaling_fmax in args.primary_fmax:
            continue
        points.append({"series": "scaling", "fmax": args.scaling_fmax,
                       "threads": threads,
                       "repetitions": args.scaling_repetitions})
    for index, point in enumerate(points):
        point["point_id"] = point_id(point["series"], point["fmax"],
                                     point["threads"])
        point["first_backend"] = "eigen" if index % 2 == 0 else "lapack"
    return points


def load_tsv(path: Path, header: list[str]) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames != header:
            raise RuntimeError(f"unexpected header in {path}")
        return list(reader)


def write_tsv(path: Path, header: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=header, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        target.flush()
        os.fsync(target.fileno())


def append_rows(path: Path, header: list[str], rows: list[dict[str, Any]]) -> None:
    new_file = not path.exists()
    with path.open("a", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=header, delimiter="\t",
                                lineterminator="\n")
        if new_file:
            writer.writeheader()
        writer.writerows(rows)
        target.flush()
        os.fsync(target.fileno())


def parse_output(output: str, repetitions: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    measurements: list[dict[str, Any]] = []
    check: dict[str, Any] | None = None
    for line in output.splitlines():
        fields = line.split("\t")
        if fields[0] == "RESULT" and len(fields) == 7:
            measurements.append({
                "backend": fields[1], "repetition": int(fields[2]),
                "seconds": float(fields[3]), "rows": int(fields[4]),
                "cols": int(fields[5]), "norm": float(fields[6])})
        elif fields[0] == "CHECK" and len(fields) == 7:
            check = {"rows": int(fields[1]), "cols": int(fields[2]),
                     "eigen_norm": float(fields[3]),
                     "lapack_norm": float(fields[4]),
                     "max_abs_difference": float(fields[5]),
                     "tolerance": float(fields[6])}
    if check is None or len(measurements) != 2 * repetitions:
        raise RuntimeError(f"incomplete benchmark output:\n{output}")
    for value in [*(row["seconds"] for row in measurements),
                  *(row["norm"] for row in measurements),
                  check["eigen_norm"], check["lapack_norm"],
                  check["max_abs_difference"]]:
        if not math.isfinite(value):
            raise RuntimeError("non-finite benchmark result")
    for backend in ("eigen", "lapack"):
        rows = [row for row in measurements if row["backend"] == backend]
        if len(rows) != repetitions:
            raise RuntimeError(f"expected {repetitions} {backend} measurements")
    return measurements, check


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tag", default="final_production")
    parser.add_argument("--lmax", type=int, default=750)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--primary-fmax", type=float, nargs="+",
                        default=PRIMARY_FMAX)
    parser.add_argument("--primary-threads", type=int, default=32)
    parser.add_argument("--primary-repetitions", type=int, default=3)
    parser.add_argument("--single-fmax", type=float, nargs="*",
                        default=SINGLE_FMAX)
    parser.add_argument("--single-repetitions", type=int, default=1)
    parser.add_argument("--scaling-fmax", type=float, default=35.0)
    parser.add_argument("--scaling-threads", type=int, nargs="*",
                        default=SCALING_THREADS)
    parser.add_argument("--scaling-repetitions", type=int, default=2)
    args = parser.parse_args()

    positive = [args.lmax, args.primary_threads, args.primary_repetitions,
                args.single_repetitions, args.scaling_repetitions,
                *args.scaling_threads]
    if any(value < 1 for value in positive) or args.tolerance <= 0.0:
        parser.error("dimensions, thread counts, repetitions, and tolerance must be positive")

    args.repo = args.repo.resolve()
    args.build = args.build.resolve()
    args.output_dir = args.output_dir.resolve()
    executable = args.build / "bin/final_solver_benchmark"
    template = args.repo / "data/params/ex3.txt"
    source = args.repo / "benchmarks/final_solver_benchmark.cpp"
    if not executable.is_file():
        raise RuntimeError(f"benchmark executable not found: {executable}")

    compile_commands = run_text(
        ["ninja", "-C", str(args.build), "-t", "commands",
         "final_solver_benchmark"])
    if "DSPECM1D_ENABLE_LAPACK_BAND_SOLVER" not in compile_commands:
        raise RuntimeError("benchmark build does not enable LAPACK capability")
    if "DSPECM1D_ENABLE_PROFILING" in compile_commands:
        raise RuntimeError("primary benchmark executable unexpectedly enables profiling")

    points = requested_points(args)
    request = {
        "lmax": args.lmax, "tolerance": args.tolerance, "points": points,
        "template_sha256": sha256(template),
        "executable_sha256": sha256(executable),
        "harness_source_sha256": sha256(source),
        "driver_sha256": sha256(Path(__file__)),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / f"{args.tag}_timings.tsv"
    checks_path = args.output_dir / f"{args.tag}_comparisons.tsv"
    metadata_path = args.output_dir / f"{args.tag}_metadata.json"

    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("request") != request:
            raise RuntimeError("existing tag has different campaign settings or hashes")
    else:
        metadata = {
            "campaign": "final production same-binary solver benchmark",
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "state": "running", "request": request,
            "measurement": "complete preferred adaptive multi-SEM spectra() call",
            "physics": "ex3/PREM all modes, attenuation enabled, identical per backend",
            "warmup": "one untimed spectra() call per backend and point",
            "backend_order": "alternating first backend across ordered points",
            "scaling_32_threads": "reuse primary fmax=35 mHz, 32-thread measurements",
            "blas_threads": 1, "omp_dynamic": "FALSE",
            "build": build_metadata(args.build),
            "machine": {"hostname": platform.node(), **cpu_metadata()},
            "compiler_version": run_text(
                ["/opt/gcc-15.2.0/bin/g++", "--version"]).splitlines()[0],
            "git_revision": run_text(["git", "rev-parse", "HEAD"], cwd=args.repo),
            "git_status_at_start": run_text(
                ["git", "status", "--short"], cwd=args.repo),
            "completed_points": [],
        }
        atomic_json(metadata_path, metadata)

    raw_rows = load_tsv(raw_path, RAW_HEADER)
    check_rows = load_tsv(checks_path, CHECK_HEADER)
    completed = {row["point_id"] for row in check_rows}
    point_by_id = {point["point_id"]: point for point in points}
    for completed_id in completed:
        if completed_id not in point_by_id:
            raise RuntimeError(f"unexpected completed point {completed_id}")
        expected = 2 * point_by_id[completed_id]["repetitions"]
        actual = sum(row["point_id"] == completed_id for row in raw_rows)
        if actual != expected:
            raise RuntimeError(f"completed point {completed_id} has {actual}/{expected} rows")

    incomplete_ids = {row["point_id"] for row in raw_rows} - completed
    if incomplete_ids:
        raw_rows = [row for row in raw_rows if row["point_id"] not in incomplete_ids]
        write_tsv(raw_path, RAW_HEADER, raw_rows)
        print(f"discarded incomplete rows for: {sorted(incomplete_ids)}",
              file=sys.stderr, flush=True)

    for point in points:
        if point["point_id"] in completed:
            print(f"already complete: {point['point_id']}", file=sys.stderr,
                  flush=True)
            continue
        parameter = (args.build / "data/params" /
                     f"{args.tag}_{point['point_id']}.txt")
        replace_parameter(template, parameter, point["fmax"], args.lmax)
        environment = os.environ.copy()
        environment.update({
            "OMP_NUM_THREADS": str(point["threads"]),
            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
        })
        command = [str(executable), str(parameter), point["first_backend"],
                   str(point["repetitions"]), str(args.tolerance)]
        print(f"starting {point['point_id']} first={point['first_backend']}",
              file=sys.stderr, flush=True)
        output = subprocess.check_output(command, text=True, env=environment)
        measurements, check = parse_output(output, point["repetitions"])

        decorated = [{"point_id": point["point_id"],
                      "series": point["series"],
                      "fmax_mHz": f"{point['fmax']:.12g}",
                      "threads": point["threads"], **measurement,
                      "first_backend": point["first_backend"]}
                     for measurement in measurements]
        append_rows(raw_path, RAW_HEADER, decorated)
        append_rows(checks_path, CHECK_HEADER, [{
            "point_id": point["point_id"], "series": point["series"],
            "fmax_mHz": f"{point['fmax']:.12g}",
            "threads": point["threads"], **check}])
        completed.add(point["point_id"])
        metadata["completed_points"] = sorted(completed)
        metadata["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        atomic_json(metadata_path, metadata)
        medians = {backend: statistics.median(
            row["seconds"] for row in measurements if row["backend"] == backend)
            for backend in ("eigen", "lapack")}
        print(f"completed {point['point_id']}: medians={medians}, "
              f"max_difference={check['max_abs_difference']:.6g}",
              file=sys.stderr, flush=True)

    metadata["state"] = "complete"
    metadata["completed_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    metadata["raw_sha256"] = sha256(raw_path)
    metadata["comparisons_sha256"] = sha256(checks_path)
    atomic_json(metadata_path, metadata)
    print(json.dumps({"raw": str(raw_path), "comparisons": str(checks_path),
                      "metadata": str(metadata_path)}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
