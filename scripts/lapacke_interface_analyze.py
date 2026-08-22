#!/usr/bin/env python3
"""Validate and summarize the focused LAPACKE high-level/work probe."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

THREADS = (1, 10, 32, 50)
CONFIGURATIONS = ((17, 4, 4, 512), (256, 8, 8, 32),
                  (1554, 15, 15, 1))
OPERATIONS = ("zgbtrs", "zgbtrf")
INTERFACES = ("high_level", "work")
HEADER = ["threads", "operation", "interface", "seconds", "n", "kl",
          "ku", "nrhs", "systems", "repeats", "max_abs_difference"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != HEADER:
            raise RuntimeError(f"unexpected header: {reader.fieldnames}")
        rows = list(reader)

    expected = {
        (threads, operation, interface, n)
        for n, _, _, _ in CONFIGURATIONS
        for threads in THREADS
        for operation in OPERATIONS
        for interface in INTERFACES
    }
    observed = set()
    for row in rows:
        threads = int(row["threads"])
        n = int(row["n"])
        key = (threads, row["operation"], row["interface"], n)
        if key in observed:
            raise RuntimeError(f"duplicate row: {key}")
        observed.add(key)
        config = next((item for item in CONFIGURATIONS if item[0] == n), None)
        if config is None:
            raise RuntimeError(f"unexpected dimension: {n}")
        _, kl, ku, solve_repeats = config
        expected_repeats = solve_repeats if row["operation"] == "zgbtrs" else 1
        if (threads not in THREADS or int(row["kl"]) != kl
                or int(row["ku"]) != ku or int(row["nrhs"]) != 3
                or int(row["systems"]) != 128
                or int(row["repeats"]) != expected_repeats):
            raise RuntimeError(f"invalid configuration row: {row}")
        seconds = float(row["seconds"])
        difference = float(row["max_abs_difference"])
        if not math.isfinite(seconds) or seconds <= 0.0:
            raise RuntimeError(f"invalid timing: {row}")
        if not math.isfinite(difference) or difference > 1.0e-13:
            raise RuntimeError(f"interface output mismatch: {row}")
    if observed != expected:
        raise RuntimeError(f"missing/unexpected rows: {expected ^ observed}")
    return rows


def summarize(rows: list[dict]) -> dict:
    timings = {
        (int(row["n"]), int(row["threads"]), row["operation"],
         row["interface"]): float(row["seconds"])
        for row in rows
    }
    comparisons = []
    for n, kl, ku, repeats in CONFIGURATIONS:
        for threads in THREADS:
            for operation in OPERATIONS:
                high = timings[(n, threads, operation, "high_level")]
                work = timings[(n, threads, operation, "work")]
                comparisons.append({
                    "n": n,
                    "kl": kl,
                    "ku": ku,
                    "nrhs": 3,
                    "systems": 128,
                    "repeats": repeats if operation == "zgbtrs" else 1,
                    "threads": threads,
                    "operation": operation,
                    "high_level_seconds": high,
                    "work_seconds": work,
                    "work_speedup": high / work,
                })
    tiny_solve = [item for item in comparisons
                  if item["n"] == 17 and item["operation"] == "zgbtrs"]
    return {
        "validated_rows": len(rows),
        "maximum_absolute_difference": max(
            float(row["max_abs_difference"]) for row in rows),
        "comparisons": comparisons,
        "tiny_solve_work_speedups": {
            str(item["threads"]): item["work_speedup"] for item in tiny_solve
        },
        "production_pilot_run": False,
        "decision": (
            "LAPACKE_zgbtrs_work does not materially improve concurrent "
            "tiny-system solves; do not change the production LAPACKE call "
            "on this evidence, and proceed only in a later approved stage "
            "to the proposed custom small-band-LU backsolve experiment."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--lapacke", type=Path, required=True)
    args = parser.parse_args()

    rows = read_rows(args.raw)
    summary = summarize(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "lapacke_interface_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    metadata = {
        "revision": "711eb5d36c15475093d4b18793f3995d89d2265b",
        "timing_policy": "one untimed warm-up and one measured run per point",
        "thread_environment": {
            "OMP_DYNAMIC": "FALSE",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "BLIS_NUM_THREADS": "1",
        },
        "source_package": {
            "name": "Debian lapack 3.11.0-2",
            "upstream_tag": "v3.11.0",
            "archive_sha256": (
                "4b9ba79bfd4921ca820e83979db76ab3363155709444a787979e81c22285ffa9"
            ),
            "archive_url": (
                "http://deb.debian.org/debian/pool/main/l/lapack/"
                "lapack_3.11.0.orig.tar.gz"
            ),
            "wrapper_source_sha256": {
                "lapacke_zgbtrs.c": (
                    "4f5f757ba97a9dfa35064a50d4ca652d8abe9c68474297108fdf29c3762531b3"
                ),
                "lapacke_zgbtrs_work.c": (
                    "7685fe064c3605159c8ea8d693ffbb55ddbe5301997ea40b117981e8b38406cc"
                ),
                "lapacke_zgbtrf.c": (
                    "466ae668e6f25247e801783dd1242c5a54b06d3b0956db57b04658f36a6f15f9"
                ),
                "lapacke_zgbtrf_work.c": (
                    "4450fd4b785d7d6a556d0186d7b3a6b1db77d97da624ffb010a54c93c7621167"
                ),
            },
        },
        "executable": {
            "path": str(args.executable.resolve()),
            "sha256": sha256(args.executable),
        },
        "lapacke": {
            "path": str(args.lapacke.resolve()),
            "sha256": sha256(args.lapacke),
        },
        "raw": {"path": str(args.raw.resolve()), "sha256": sha256(args.raw)},
        "summary": {
            "path": str(summary_path.resolve()),
            "sha256": sha256(summary_path),
        },
    }
    metadata_path = args.output_dir / "lapacke_interface_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
