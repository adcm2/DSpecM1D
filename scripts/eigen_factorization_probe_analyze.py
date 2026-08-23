#!/usr/bin/env python3
"""Validate and summarize the benchmark-only SparseLU factorization probe."""
from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def finite(row: dict[str, str]) -> bool:
    keys = ("factor_s", "solve_s", "pivot_fraction", "discrepancy", "residual")
    return row["finite"] == "1" and all(math.isfinite(float(row[k])) for k in keys)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: eigen_factorization_probe_analyze.py RAW SUMMARY", file=sys.stderr)
        return 2
    raw, out = map(Path, sys.argv[1:])
    with raw.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    required = {
        "factor_s", "solve_s", "nnz_l", "nnz_u", "factor_info", "solve_info",
        "discrepancy", "residual", "finite", "repeats",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("raw file has missing columns")
    if any(not finite(row) for row in rows):
        raise SystemExit("non-finite benchmark row")
    if any(row["factor_info"] != "0" or row["solve_info"] != "0" for row in rows):
        raise SystemExit("factorization or solve failure")
    if any(int(row["repeats"]) != 4 for row in rows):
        raise SystemExit("unexpected measured batch repeat count")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["config"]].append(row)
    expected_configs = {
        "pivot_1", "pivot_0p5", "pivot_0p1", "pivot_0p01", "pivot_0",
        "natural_default", "panel_1", "panel_4", "panel_8", "panel_32",
        "relax_4", "maxsuper_64", "maxsuper_256", "panel_4_relax_4",
        "panel_4_maxsuper_64",
    }
    expected_systems = {
        "radial_small", "toroidal_small", "spheroidal_small",
        "spheroidal_medium", "spheroidal_large", "spheroidal_fluid_solid",
        "spheroidal_max",
    }
    if set(grouped) != expected_configs:
        raise SystemExit(f"unexpected config set: {sorted(grouped)}")
    if {r["system"] for r in rows} != expected_systems:
        raise SystemExit("unexpected system set")
    baseline_keys = {(r["system"], r["freq"]) for r in grouped["pivot_1"]}
    if len(rows) != len(baseline_keys) * len(expected_configs):
        raise SystemExit("unexpected measured row count")
    if any({(r["system"], r["freq"]) for r in group} != baseline_keys
           for group in grouped.values()):
        raise SystemExit("configurations do not have identical coverage")

    def stats(group: list[dict[str, str]]) -> dict[str, float | int]:
        factors = [float(r["factor_s"]) for r in group]
        solves = [float(r["solve_s"]) for r in group]
        fills = [(int(r["nnz_l"]) + int(r["nnz_u"])) / int(r["nnz_a"]) for r in group]
        return {
            "rows": len(group),
            "factor_mean_s": statistics.mean(factors),
            "factor_median_s": statistics.median(factors),
            "solve_mean_s": statistics.mean(solves),
            "solve_median_s": statistics.median(solves),
            "factor_sum_s": sum(factors),
            "solve_sum_s": sum(solves),
            "fill_ratio_mean": statistics.mean(fills),
            "nnz_l_min": min(int(r["nnz_l"]) for r in group),
            "nnz_l_max": max(int(r["nnz_l"]) for r in group),
            "nnz_u_min": min(int(r["nnz_u"]) for r in group),
            "nnz_u_max": max(int(r["nnz_u"]) for r in group),
            "pivot_fraction_mean": statistics.mean(float(r["pivot_fraction"]) for r in group),
            "max_pivot_displacement": max(int(r["max_pivot_displacement"]) for r in group),
            "max_discrepancy": max(float(r["discrepancy"]) for r in group),
            "max_backward_residual": max(float(r["residual"]) for r in group),
        }

    max_discrepancy = max(float(r["discrepancy"]) for r in rows)
    max_residual = max(float(r["residual"]) for r in rows)
    if max_discrepancy > 1.e-9 or max_residual > 1.e-12:
        raise SystemExit("numerical validation gate failed")
    summary = {
        "experiment": "eigen_factorization_probe",
        "rows": len(rows),
        "systems": sorted({r["system"] for r in rows}),
        "frequencies_by_system": {
            system: sorted({float(r["freq"]) for r in rows if r["system"] == system})
            for system in sorted({r["system"] for r in rows})
        },
        "configs": {name: stats(grouped[name]) for name in sorted(grouped)},
        "baseline": "pivot_1",
        "config_speedups_vs_pivot_1": {
            name: stats(grouped["pivot_1"])["factor_sum_s"] / stats(group)["factor_sum_s"]
            for name, group in sorted(grouped.items())
        },
        "validation": {
            "all_finite": True,
            "all_factor_info_zero": True,
            "all_solve_info_zero": True,
            "max_discrepancy": max_discrepancy,
            "max_backward_residual": max_residual,
            "discrepancy_gate": 1.e-9,
            "backward_residual_gate": 1.e-12,
        },
    }
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
