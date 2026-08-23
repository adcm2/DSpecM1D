#!/usr/bin/env python3
"""Validate and summarize the structure-aware SEM factorization probe."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


EXPECTED = {
    "radial_0S0",
    "toroidal_0T12",
    "spheroidal_0S12",
    "fluid_solid_0S40",
    "large_0S40",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key} for {row['name']}")
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("raw", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("table", type=Path)
    args = parser.parse_args()

    with args.raw.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 5 or {row["name"] for row in rows} != EXPECTED:
        raise ValueError("expected exactly the five documented real SEM systems")

    timing_keys = [
        "one_time_base_extraction_s", "local_dynamic_prep_s",
        "local_factor_s", "schur_s", "reduced_assembly_s",
        "reduced_band_prep_s", "interface_factor_s", "rhs_condense_s", "interface_solve_s",
        "recovery_s", "condensed_setup_intrinsic_s",
        "condensed_setup_with_scaffolding_s", "condensed_solve_recovery_s",
        "eigen_factor_s", "eigen_solve_s", "lapack_factor_s",
        "lapack_solve_s",
    ]
    for row in rows:
        for key in timing_keys:
            if number(row, key) <= 0:
                raise ValueError(f"non-positive {key} for {row['name']}")
        for key in ("condensed_info", "eigen_info", "lapack_info"):
            if int(row[key]) != 0:
                raise ValueError(f"failed solver status for {row['name']}")
        if int(row["factors_finite"]) != 1:
            raise ValueError(f"non-finite factor storage for {row['name']}")
        if int(row["finite"]) != 1:
            raise ValueError(f"non-finite solution for {row['name']}")
        for key in ("condensed_residual", "eigen_residual", "lapack_residual"):
            if number(row, key) >= 1e-12:
                raise ValueError(f"backward residual failed for {row['name']}")
        for key in ("condensed_vs_eigen", "condensed_vs_lapack",
                    "eigen_vs_lapack"):
            if number(row, key) >= 1e-7:
                raise ValueError(f"solution comparison failed for {row['name']}")
        if int(row["total_interior"]) + int(row["interface_dofs"]) != int(row["n"]):
            raise ValueError(f"partition does not cover {row['name']}")
        if int(row["reduced_n"]) != int(row["interface_dofs"]):
            raise ValueError(f"reduced dimension mismatch for {row['name']}")
        if int(row["reduced_kl"]) > 2 * int(row["components"]) or \
           int(row["reduced_ku"]) > 2 * int(row["components"]):
            raise ValueError(f"unexpectedly wide reduced system for {row['name']}")

    def total(keys):
        return sum(sum(number(row, key) for key in keys) for row in rows)

    lapack_factor = total(["lapack_factor_s"])
    lapack_total = total(["lapack_factor_s", "lapack_solve_s"])
    intrinsic_setup = total(["condensed_setup_intrinsic_s"])
    intrinsic_total = total([
        "condensed_setup_intrinsic_s", "condensed_solve_recovery_s"
    ])
    scaffold_total = total([
        "condensed_setup_with_scaffolding_s", "condensed_solve_recovery_s"
    ])
    summary = {
        "decision": "C_NO_WIN",
        "systems": len(rows),
        "warmups_per_solver_system": 1,
        "measured_repeats_per_solver_system": int(rows[0]["repeats"]),
        "aggregate": {
            "lapack_factor_s": lapack_factor,
            "lapack_factor_plus_solve_s": lapack_total,
            "condensed_intrinsic_setup_s": intrinsic_setup,
            "condensed_intrinsic_setup_plus_full_solve_s": intrinsic_total,
            "condensed_with_benchmark_scaffolding_plus_full_solve_s": scaffold_total,
            "setup_speedup_vs_lapack_factor": lapack_factor / intrinsic_setup,
            "intrinsic_total_speedup_vs_lapack": lapack_total / intrinsic_total,
            "scaffold_total_speedup_vs_lapack": lapack_total / scaffold_total,
        },
        "ranges": {
            "n": [min(int(row["n"]) for row in rows),
                  max(int(row["n"]) for row in rows)],
            "reduced_n": [min(int(row["reduced_n"]) for row in rows),
                          max(int(row["reduced_n"]) for row in rows)],
            "intrinsic_total_speedup_vs_lapack": [
                min(number(row, "total_speedup_vs_lapack") for row in rows),
                max(number(row, "total_speedup_vs_lapack") for row in rows),
            ],
        },
        "validation": {
            "max_condensed_residual": max(number(row, "condensed_residual")
                                            for row in rows),
            "max_condensed_vs_eigen": max(number(row, "condensed_vs_eigen")
                                            for row in rows),
            "max_condensed_vs_lapack": max(number(row, "condensed_vs_lapack")
                                             for row in rows),
        },
        "raw_sha256": sha256(args.raw),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")

    columns = [
        "name", "n", "nnz", "elements", "total_interior",
        "interface_dofs", "reduced_kl", "reduced_ku",
        "condensed_setup_intrinsic_s", "condensed_solve_recovery_s",
        "lapack_factor_s", "lapack_solve_s", "total_speedup_vs_lapack",
        "condensed_residual", "condensed_vs_lapack",
    ]
    with args.table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
