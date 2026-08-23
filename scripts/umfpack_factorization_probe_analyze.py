#!/usr/bin/env python3
"""Validate and summarize the focused UMFPACK factorization probe."""

import argparse
import csv
import json
import math
from pathlib import Path


CONFIGS = [
    ("robust", "natural"),
    ("robust", "default_amd"),
    ("default", "natural"),
    ("default", "default_amd"),
]
SYSTEMS = {
    "radial_small",
    "toroidal_small",
    "spheroidal_small",
    "spheroidal_medium",
    "spheroidal_large",
    "spheroidal_fluid_solid",
    "spheroidal_max",
}
FLOAT_FIELDS = {
    "eigen_setup_s", "eigen_factor_s", "eigen_solve_s",
    "umf_symbolic_s", "umf_numeric_s", "umf_raw_solve_s",
    "umf_refined_solve_s", "umf_strategy_used", "umf_ordering_used",
    "umf_was_scaled", "umf_noff_diag", "umf_lnz", "umf_unz",
    "umf_flops", "umf_numeric_size", "umf_peak_memory", "umf_rcond",
    "umf_raw_ir_steps", "umf_refined_ir_steps", "eigen_residual",
    "umf_raw_residual", "umf_refined_residual", "lapack_residual",
    "umf_raw_vs_eigen", "umf_refined_vs_eigen", "eigen_vs_lapack",
    "umf_raw_vs_lapack", "umf_refined_vs_lapack",
    "transfer_raw_vs_eigen", "transfer_refined_vs_eigen",
}


def location(row):
    return {
        "system": row["system"],
        "frequency_index": int(row["freq_index"]),
        "frequency_mhz": float(row["freq_mhz"]),
        "policy": row["policy"],
        "ordering": row["ordering"],
    }


def maximum(rows, field):
    row = max(rows, key=lambda item: float(item[field]))
    return {"value": float(row[field]), **location(row)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("raw", type=Path)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("summary_tsv", type=Path)
    args = parser.parse_args()

    with args.raw.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 56:
        raise SystemExit(f"expected 56 rows, found {len(rows)}")
    if {row["system"] for row in rows} != SYSTEMS:
        raise SystemExit("physical system set differs from the accepted seven-system set")
    for row in rows:
        for field in FLOAT_FIELDS:
            if not math.isfinite(float(row[field])):
                raise SystemExit(f"non-finite {field}: {location(row)}")
        for field in (
            "eigen_factor_info", "eigen_solve_info", "umf_symbolic_status",
            "umf_numeric_status", "umf_raw_status", "umf_refined_status",
            "lapack_factor_info", "lapack_solve_info",
        ):
            if int(row[field]) != 0:
                raise SystemExit(f"non-success {field}: {location(row)}")
        if int(row["finite"]) != 1:
            raise SystemExit(f"failed finite flag: {location(row)}")
        if int(row["repeats"]) != 16:
            raise SystemExit("unexpected measured batch size")
        if float(row["umf_was_scaled"]) != 1:
            raise SystemExit("UMFPACK did not use documented row-sum scaling")

    keys = {(r["system"], r["freq_index"], r["policy"], r["ordering"])
            for r in rows}
    expected = {(s, str(f), p, o) for s in SYSTEMS for f in (0, 1)
                for p, o in CONFIGS}
    if keys != expected:
        raise SystemExit("missing or duplicate system/frequency/configuration rows")

    aggregates = []
    for policy, ordering in CONFIGS:
        selected = [r for r in rows
                    if r["policy"] == policy and r["ordering"] == ordering]
        eigen_factor = sum(float(r["eigen_factor_s"]) for r in selected)
        umf_numeric = sum(float(r["umf_numeric_s"]) for r in selected)
        eigen_solve = sum(float(r["eigen_solve_s"]) for r in selected)
        raw_solve = sum(float(r["umf_raw_solve_s"]) for r in selected)
        refined_solve = sum(float(r["umf_refined_solve_s"]) for r in selected)
        aggregate = {
            "policy": policy,
            "ordering": ordering,
            "matrix_count": len(selected),
            # One Symbolic object serves both frequency rows for each system.
            "symbolic_seconds_sum": sum(float(r["umf_symbolic_s"])
                                          for r in selected
                                          if int(r["freq_index"]) == 0),
            "eigen_factor_seconds_sum": eigen_factor,
            "umf_numeric_seconds_sum": umf_numeric,
            "eigen_over_umf_numeric": eigen_factor / umf_numeric,
            "umf_numeric_percent_slower": 100.0 * (umf_numeric / eigen_factor - 1.0),
            "eigen_solve_seconds_sum": eigen_solve,
            "umf_raw_solve_seconds_sum": raw_solve,
            "eigen_over_umf_raw_solve": eigen_solve / raw_solve,
            "umf_refined_solve_seconds_sum": refined_solve,
            "strategy_used": sorted({int(float(r["umf_strategy_used"]))
                                      for r in selected}),
            "ordering_used": sorted({int(float(r["umf_ordering_used"]))
                                      for r in selected}),
            "total_off_diagonal_pivots": sum(int(float(r["umf_noff_diag"]))
                                               for r in selected),
        }
        aggregates.append(aggregate)

    summary = {
        "decision": "C",
        "decision_text": "robust UMFPACK is slower than Eigen",
        "records": len(rows),
        "systems": len(SYSTEMS),
        "dimensions": {
            "minimum": min(int(r["n"]) for r in rows),
            "maximum": max(int(r["n"]) for r in rows),
        },
        "aggregate": aggregates,
        "maxima": {field: maximum(rows, field) for field in (
            "eigen_residual", "umf_raw_residual", "umf_refined_residual",
            "lapack_residual", "umf_raw_vs_eigen", "umf_refined_vs_eigen",
            "eigen_vs_lapack", "umf_raw_vs_lapack",
            "umf_refined_vs_lapack", "transfer_raw_vs_eigen",
            "transfer_refined_vs_eigen",
        )},
        "validation": {
            "all_statuses_success": True,
            "all_outputs_finite": True,
            "maximum_residual_below_1e-12": max(
                float(r[f]) for r in rows for f in (
                    "eigen_residual", "umf_raw_residual",
                    "umf_refined_residual", "lapack_residual")) < 1e-12,
            "maximum_solution_discrepancy_below_1e-6": max(
                float(r[f]) for r in rows for f in (
                    "umf_raw_vs_eigen", "umf_refined_vs_eigen",
                    "eigen_vs_lapack", "umf_raw_vs_lapack",
                    "umf_refined_vs_lapack")) < 1e-6,
        },
        "concurrency_probe": "not run: robust serial UMFPACK failed the performance gate",
        "production_spectra": "not run",
    }
    if not all(summary["validation"].values()):
        raise SystemExit("numerical validation failed")

    args.summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    with args.summary_tsv.open("w", newline="") as stream:
        fields = [
            "policy", "ordering", "matrix_count", "symbolic_seconds_sum",
            "eigen_factor_seconds_sum", "umf_numeric_seconds_sum",
            "eigen_over_umf_numeric", "umf_numeric_percent_slower",
            "eigen_solve_seconds_sum", "umf_raw_solve_seconds_sum",
            "eigen_over_umf_raw_solve", "umf_refined_solve_seconds_sum",
            "strategy_used", "ordering_used", "total_off_diagonal_pivots",
        ]
        writer = csv.DictWriter(stream, fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for aggregate in aggregates:
            writer.writerow(aggregate)


if __name__ == "__main__":
    main()
