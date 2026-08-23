#!/usr/bin/env python3
"""Validate and summarize the focused Eigen pivot-threshold experiment."""
import csv
import json
import math
import sys
from collections import Counter


def f(row, key):
    return float(row[key])


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: eigen_pivot0_validation_analyze.py RAW SUMMARY")
    with open(sys.argv[1], newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or len(rows) % 2:
        raise SystemExit("expected paired threshold-0/1 rows")
    required = {
        "family", "name", "freq_mhz", "epsilon", "attenuation", "threshold",
        "factor_info", "solve_info", "finite", "max_abs_a", "min_abs_diag_u", "max_abs_u",
        "max_abs_l", "growth_u_over_a", "residual", "solution_vs_threshold1",
        "transfer_vs_threshold1", "transfer_relative_vs_threshold1", "lapack_factor_info",
        "lapack_solve_info", "lapack_finite", "lapack_residual", "lapack_solution_vs_threshold1",
        "lapack_transfer_vs_threshold1", "lapack_transfer_relative_vs_threshold1",
        "solution_vs_lapack", "transfer_vs_lapack", "transfer_relative_vs_lapack",
    }
    if not required.issubset(rows[0]):
        raise SystemExit("raw TSV is missing required diagnostics")
    expected_names = {
        "radial_0S0", "toroidal_0T2", "toroidal_0T12",
        "spheroidal_0S2", "spheroidal_0S12", "fluid_solid_0S40",
    }
    if len(rows) != 216 or len({r["name"] for r in rows}) != 6:
        raise SystemExit("expected exactly 216 rows and six physical systems")
    if {r["name"] for r in rows} != expected_names:
        raise SystemExit("unexpected physical system set")
    if {r["threshold"] for r in rows} != {"0", "1"}:
        raise SystemExit("expected both thresholds")
    if {r["attenuation"] for r in rows} != {"0", "1"}:
        raise SystemExit("expected both attenuation states")
    if len({r["epsilon"] for r in rows}) != 3 or any(
        len({r["freq_mhz"] for r in rows if r["name"] == name}) != 3
        for name in expected_names
    ):
        raise SystemExit("expected three epsilons and three grid points per system")
    for row in rows:
        for key in ("factor_info", "solve_info", "finite", "lapack_factor_info", "lapack_solve_info", "lapack_finite"):
            if int(row[key]) != 0 and key in {"factor_info", "solve_info", "lapack_factor_info", "lapack_solve_info"}:
                raise SystemExit(f"nonzero {key}: {row}")
            if key == "finite" and int(row[key]) != 1:
                raise SystemExit(f"nonfinite result: {row}")
            if key == "lapack_finite" and int(row[key]) != 1:
                raise SystemExit(f"nonfinite LAPACK result: {row}")
        for key in required - {"family", "name", "freq_mhz", "epsilon", "attenuation", "threshold"}:
            if key in {"factor_info", "solve_info", "finite", "lapack_factor_info", "lapack_solve_info"}:
                continue
            if not math.isfinite(f(row, key)):
                raise SystemExit(f"nonfinite {key}: {row}")
    zero = [r for r in rows if r["threshold"] == "0"]
    one = [r for r in rows if r["threshold"] == "1"]
    key_fields = ("name", "freq_mhz", "epsilon", "attenuation")
    if len(zero) != len(one) or {tuple(r[k] for k in key_fields) for r in zero} != {tuple(r[k] for k in key_fields) for r in one}:
        raise SystemExit("threshold rows are not paired")
    max_relative_discrepancy = max(
        max(f(r, "solution_vs_threshold1") for r in zero),
        max(f(r, "transfer_relative_vs_threshold1") for r in zero),
        max(f(r, "solution_vs_lapack") for r in zero),
        max(f(r, "transfer_relative_vs_lapack") for r in zero),
    )
    def identity(row):
        return {
            "name": row["name"], "freq_mhz": float(row["freq_mhz"]),
            "epsilon": float(row["epsilon"]), "attenuation": int(row["attenuation"]),
            "threshold": int(row["threshold"]),
        }

    max_l_zero = max(f(r, "max_abs_l") for r in zero)
    max_l_one = max(f(r, "max_abs_l") for r in one)
    l_multiplier_ratio = max_l_zero / max_l_one
    growth_rejected = max_l_zero > 100.0 or l_multiplier_ratio > 100.0

    summary = {
        "rows": len(rows),
        "expected_rows": 216,
        "threshold_zero_rows": len(zero),
        "threshold_one_rows": len(one),
        "grid_points_per_system": 3,
        "epsilon_count": 3,
        "attenuation_count": 2,
        "families": dict(Counter(r["family"] for r in zero)),
        "systems": sorted(set(r["name"] for r in zero)),
        "epsilons": sorted(set(f(r, "epsilon") for r in zero)),
        "attenuation_values": sorted(set(int(r["attenuation"]) for r in zero)),
        "max_threshold_zero": {
            "residual": max(f(r, "residual") for r in zero),
            "solution_vs_threshold1": max(f(r, "solution_vs_threshold1") for r in zero),
            "transfer_vs_threshold1": max(f(r, "transfer_vs_threshold1") for r in zero),
            "transfer_relative_vs_threshold1": max(f(r, "transfer_relative_vs_threshold1") for r in zero),
            "solution_vs_lapack": max(f(r, "solution_vs_lapack") for r in zero),
            "transfer_vs_lapack": max(f(r, "transfer_vs_lapack") for r in zero),
            "transfer_relative_vs_lapack": max(f(r, "transfer_relative_vs_lapack") for r in zero),
            "growth_u_over_a": max(f(r, "growth_u_over_a") for r in zero),
            "max_abs_u": max(f(r, "max_abs_u") for r in zero),
            "max_abs_l": max(f(r, "max_abs_l") for r in zero),
            "max_abs_a": max(f(r, "max_abs_a") for r in zero),
        },
        "max_threshold_one": {
            "max_abs_l": max(f(r, "max_abs_l") for r in one),
            "max_abs_u": max(f(r, "max_abs_u") for r in one),
            "growth_u_over_a": max(f(r, "growth_u_over_a") for r in one),
            "min_abs_diag_u": min(f(r, "min_abs_diag_u") for r in one),
            "max_abs_a": max(f(r, "max_abs_a") for r in one),
        },
        "max_lapack_reference": {
            "residual": max(f(r, "lapack_residual") for r in zero),
            "solution_vs_threshold1": max(f(r, "lapack_solution_vs_threshold1") for r in zero),
            "transfer_vs_threshold1": max(f(r, "lapack_transfer_vs_threshold1") for r in zero),
            "transfer_relative_vs_threshold1": max(f(r, "lapack_transfer_relative_vs_threshold1") for r in zero),
        },
        "threshold_zero_nonidentity_pivots": max(int(r["nonidentity_pivots"]) for r in zero),
        "l_multiplier_ratio_zero_to_one": l_multiplier_ratio,
        "growth_ratio_rejection_threshold": 100.0,
        "growth_absolute_rejection_threshold": 100.0,
        "growth_rejection_observed": growth_rejected,
        "worst_cases": {
            "threshold_zero_max_l": identity(max(zero, key=lambda r: f(r, "max_abs_l"))),
            "threshold_zero_max_solution_difference": identity(max(zero, key=lambda r: f(r, "solution_vs_threshold1"))),
            "threshold_zero_max_transfer_difference": identity(max(zero, key=lambda r: f(r, "transfer_relative_vs_threshold1"))),
        },
        "decision": "PASS" if all(int(r["finite"]) == 1 and int(r["factor_info"]) == 0 and
                                    int(r["solve_info"]) == 0 and int(r["lapack_factor_info"]) == 0 and
                                    int(r["lapack_solve_info"]) == 0 and int(r["lapack_finite"]) == 1
                                    for r in rows) and not growth_rejected
                         else "FAIL",
    }
    with open(sys.argv[2], "w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
