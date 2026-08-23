#!/usr/bin/env python3
"""Validate and summarize the focused Eigen pivot-threshold 0.01 experiment."""
import csv
import json
import math
import sys
from collections import Counter


def number(row, key):
    return float(row[key])


def percentile(values, fraction):
    values = sorted(values)
    position = fraction * (len(values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1.0 - weight) + values[upper] * weight


def distribution(rows, field):
    values = [number(row, field) for row in rows]
    return {
        "minimum": min(values),
        "median": percentile(values, 0.5),
        "p95": percentile(values, 0.95),
        "maximum": max(values),
    }


def identity(row):
    return {
        "name": row["name"],
        "freq_mhz": float(row["freq_mhz"]),
        "epsilon": float(row["epsilon"]),
        "attenuation": int(row["attenuation"]),
        "threshold": float(row["threshold"]),
    }


def main():
    if len(sys.argv) != 4 or sys.argv[3] not in {"PASS", "FAIL"}:
        raise SystemExit(
            "usage: eigen_pivot001_validation_analyze.py RAW SUMMARY PASS|FAIL"
        )
    with open(sys.argv[1], newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    required = {
        "family", "name", "freq_mhz", "epsilon", "attenuation", "threshold", "n",
        "factor_info", "solve_info", "finite", "max_abs_a", "min_abs_diag_u",
        "max_abs_diag_u", "max_abs_u", "max_abs_l", "growth_u_over_a",
        "nonidentity_pivots", "pivot_fraction", "max_pivot_displacement",
        "residual", "solution_vs_threshold1", "transfer_vs_threshold1",
        "transfer_relative_vs_threshold1", "lapack_factor_info",
        "lapack_solve_info", "lapack_finite", "lapack_residual",
        "lapack_solution_vs_threshold1", "lapack_transfer_vs_threshold1",
        "lapack_transfer_relative_vs_threshold1", "solution_vs_lapack",
        "transfer_vs_lapack", "transfer_relative_vs_lapack",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("raw TSV is empty or missing required diagnostics")

    expected_names = {
        "radial_0S0", "toroidal_0T2", "toroidal_0T12",
        "spheroidal_0S2", "spheroidal_0S12", "fluid_solid_0S40",
    }
    if len(rows) != 216 or {row["name"] for row in rows} != expected_names:
        raise SystemExit("expected exactly 216 rows and the six accepted physical systems")
    if {row["threshold"] for row in rows} != {"0.01", "1"}:
        raise SystemExit("expected exactly thresholds 0.01 and 1")
    if {row["attenuation"] for row in rows} != {"0", "1"}:
        raise SystemExit("expected both attenuation states")
    if len({row["epsilon"] for row in rows}) != 3 or any(
        len({row["freq_mhz"] for row in rows if row["name"] == name}) != 3
        for name in expected_names
    ):
        raise SystemExit("expected three epsilons and three grid points per system")

    integer_status = {
        "factor_info": 0, "solve_info": 0, "finite": 1,
        "lapack_factor_info": 0, "lapack_solve_info": 0, "lapack_finite": 1,
    }
    numeric_fields = required - {
        "family", "name", "factor_info", "solve_info", "finite",
        "lapack_factor_info", "lapack_solve_info", "lapack_finite",
    }
    for row in rows:
        for field, expected in integer_status.items():
            if int(row[field]) != expected:
                raise SystemExit(f"unexpected {field}: {row}")
        for field in numeric_fields:
            if not math.isfinite(number(row, field)):
                raise SystemExit(f"nonfinite {field}: {row}")

    candidate = [row for row in rows if row["threshold"] == "0.01"]
    baseline = [row for row in rows if row["threshold"] == "1"]
    key_fields = ("name", "freq_mhz", "epsilon", "attenuation")
    candidate_by_key = {tuple(row[key] for key in key_fields): row for row in candidate}
    baseline_by_key = {tuple(row[key] for key in key_fields): row for row in baseline}
    if len(candidate) != 108 or candidate_by_key.keys() != baseline_by_key.keys():
        raise SystemExit("threshold rows are not exactly paired")

    same_case_l_ratios = []
    same_case_u_growth_ratios = []
    for key, row in candidate_by_key.items():
        reference = baseline_by_key[key]
        same_case_l_ratios.append(
            number(row, "max_abs_l") / max(number(reference, "max_abs_l"), 1.e-300)
        )
        same_case_u_growth_ratios.append(
            number(row, "growth_u_over_a") /
            max(number(reference, "growth_u_over_a"), 1.e-300)
        )

    max_candidate = {
        "residual": max(number(row, "residual") for row in candidate),
        "solution_vs_threshold1": max(number(row, "solution_vs_threshold1") for row in candidate),
        "transfer_vs_threshold1": max(number(row, "transfer_vs_threshold1") for row in candidate),
        "transfer_relative_vs_threshold1": max(
            number(row, "transfer_relative_vs_threshold1") for row in candidate
        ),
        "solution_vs_lapack": max(number(row, "solution_vs_lapack") for row in candidate),
        "transfer_vs_lapack": max(number(row, "transfer_vs_lapack") for row in candidate),
        "transfer_relative_vs_lapack": max(
            number(row, "transfer_relative_vs_lapack") for row in candidate
        ),
        "min_abs_diag_u": min(number(row, "min_abs_diag_u") for row in candidate),
        "max_abs_diag_u": max(number(row, "max_abs_diag_u") for row in candidate),
        "max_abs_u": max(number(row, "max_abs_u") for row in candidate),
        "max_abs_l": max(number(row, "max_abs_l") for row in candidate),
        "max_abs_a": max(number(row, "max_abs_a") for row in candidate),
        "growth_u_over_a": max(number(row, "growth_u_over_a") for row in candidate),
        "nonidentity_pivots": max(int(row["nonidentity_pivots"]) for row in candidate),
        "pivot_fraction": max(number(row, "pivot_fraction") for row in candidate),
        "max_pivot_displacement": max(
            int(row["max_pivot_displacement"]) for row in candidate
        ),
    }
    max_baseline = {
        "residual": max(number(row, "residual") for row in baseline),
        "min_abs_diag_u": min(number(row, "min_abs_diag_u") for row in baseline),
        "max_abs_diag_u": max(number(row, "max_abs_diag_u") for row in baseline),
        "max_abs_u": max(number(row, "max_abs_u") for row in baseline),
        "max_abs_l": max(number(row, "max_abs_l") for row in baseline),
        "max_abs_a": max(number(row, "max_abs_a") for row in baseline),
        "growth_u_over_a": max(number(row, "growth_u_over_a") for row in baseline),
        "nonidentity_pivots": max(int(row["nonidentity_pivots"]) for row in baseline),
        "pivot_fraction": max(number(row, "pivot_fraction") for row in baseline),
        "max_pivot_displacement": max(
            int(row["max_pivot_displacement"]) for row in baseline
        ),
    }
    per_system_growth = {}
    for name in sorted(expected_names):
        candidate_system = [row for row in candidate if row["name"] == name]
        baseline_system = [row for row in baseline if row["name"] == name]
        candidate_max_l = max(number(row, "max_abs_l") for row in candidate_system)
        baseline_max_l = max(number(row, "max_abs_l") for row in baseline_system)
        per_system_growth[name] = {
            "n": max(int(row["n"]) for row in candidate_system),
            "threshold_001_max_abs_l": candidate_max_l,
            "threshold_one_max_abs_l": baseline_max_l,
            "max_l_ratio_001_to_one": candidate_max_l / baseline_max_l,
            "threshold_001_max_u_over_a": max(
                number(row, "growth_u_over_a") for row in candidate_system
            ),
            "threshold_one_max_u_over_a": max(
                number(row, "growth_u_over_a") for row in baseline_system
            ),
            "threshold_001_max_pivot_fraction": max(
                number(row, "pivot_fraction") for row in candidate_system
            ),
        }

    summary = {
        "decision": sys.argv[3],
        "decision_semantics": (
            "Explicit scientific interpretation supplied after inspecting the full "
            "factor-growth distributions; it is not an arbitrary automated cutoff."
        ),
        "rows": len(rows),
        "candidate_rows": len(candidate),
        "baseline_rows": len(baseline),
        "grid_points_per_system": 3,
        "epsilon_count": 3,
        "attenuation_count": 2,
        "families": dict(Counter(row["family"] for row in candidate)),
        "systems": sorted({row["name"] for row in candidate}),
        "epsilons": sorted({number(row, "epsilon") for row in candidate}),
        "max_threshold_001": max_candidate,
        "max_threshold_one": max_baseline,
        "factor_growth_distributions": {
            "threshold_001_max_abs_l": distribution(candidate, "max_abs_l"),
            "threshold_one_max_abs_l": distribution(baseline, "max_abs_l"),
            "threshold_001_u_over_a": distribution(candidate, "growth_u_over_a"),
            "threshold_one_u_over_a": distribution(baseline, "growth_u_over_a"),
            "same_case_l_ratio_001_to_one": {
                "median": percentile(same_case_l_ratios, 0.5),
                "p95": percentile(same_case_l_ratios, 0.95),
                "maximum": max(same_case_l_ratios),
            },
            "same_case_u_growth_ratio_001_to_one": {
                "median": percentile(same_case_u_growth_ratios, 0.5),
                "p95": percentile(same_case_u_growth_ratios, 0.95),
                "maximum": max(same_case_u_growth_ratios),
            },
        },
        "per_system_growth": per_system_growth,
        "rejected_threshold_zero_comparison": {
            "max_abs_l": 5013.549795341349,
            "threshold_one_max_abs_l": 0.999999998013822,
            "ratio_zero_to_one": 5013.549805299151,
        },
        "max_lapack_reference": {
            "residual": max(number(row, "lapack_residual") for row in candidate),
            "solution_vs_threshold1": max(
                number(row, "lapack_solution_vs_threshold1") for row in candidate
            ),
            "transfer_vs_threshold1": max(
                number(row, "lapack_transfer_vs_threshold1") for row in candidate
            ),
            "transfer_relative_vs_threshold1": max(
                number(row, "lapack_transfer_relative_vs_threshold1") for row in candidate
            ),
        },
        "worst_cases": {
            "threshold_001_max_l": identity(
                max(candidate, key=lambda row: number(row, "max_abs_l"))
            ),
            "threshold_001_max_u_growth": identity(
                max(candidate, key=lambda row: number(row, "growth_u_over_a"))
            ),
            "threshold_001_max_solution_difference": identity(
                max(candidate, key=lambda row: number(row, "solution_vs_threshold1"))
            ),
            "threshold_001_max_transfer_difference": identity(
                max(candidate, key=lambda row: number(row, "transfer_relative_vs_threshold1"))
            ),
        },
    }
    with open(sys.argv[2], "w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
