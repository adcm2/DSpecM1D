#!/usr/bin/env python3
"""Validate and summarize the fixed-pattern E2 kernel rows."""

import csv
import json
import math
import pathlib
import sys


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: eigen_sparse_matrix_e2_analyze.py RAW.tsv SUMMARY.json")
    raw_path, summary_path = map(pathlib.Path, sys.argv[1:])
    with raw_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 14:
        raise ValueError(f"expected 14 rows, got {len(rows)}")
    for row in rows:
        numeric = [float(row[key]) for key in (
            "current_s", "direct_s", "speedup", "max_coeff_diff",
            "max_relative_coeff_diff", "solve_discrepancy", "residual")]
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError(f"non-finite row: {row['system']} attenuation={row['attenuation']}")
        if row["structure_ok"] != "1":
            raise ValueError(f"structure mismatch: {row['system']}")
        if (row["current_analyze_completed"] != "1" or
                row["direct_analyze_completed"] != "1"):
            raise ValueError(
                f"SparseLU analysis failed: {row['system']} attenuation={row['attenuation']}")
        info_fields = (
            "current_factorize_info", "current_solve_info",
            "direct_factorize_info", "direct_solve_info")
        if any(int(row[key]) != 0 for key in info_fields):
            raise ValueError(
                f"SparseLU failure: {row['system']} attenuation={row['attenuation']}")
        if float(row["max_coeff_diff"]) > 1e-12 or float(row["solve_discrepancy"]) > 1e-12:
            raise ValueError(f"numerical mismatch: {row['system']} attenuation={row['attenuation']}")
    summary = {
        "rows": len(rows),
        "systems": sorted({row["system"] for row in rows}),
        "attenuation_rows": {key: sum(row["attenuation"] == key for row in rows)
                              for key in ("0", "1")},
        "min_speedup": min(float(row["speedup"]) for row in rows),
        "max_speedup": max(float(row["speedup"]) for row in rows),
        "max_residual": max(float(row["residual"]) for row in rows),
        "max_solution_discrepancy": max(float(row["solve_discrepancy"]) for row in rows),
        "max_coefficient_discrepancy": max(float(row["max_coeff_diff"]) for row in rows),
        "all_analysis_completed": True,
        "all_factorize_solve_info_zero": True,
        "all_union_counts_equal": all(
            row["nnz_union_off"] == row["nnz_union_on"] for row in rows),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
