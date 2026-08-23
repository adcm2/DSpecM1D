#!/usr/bin/env python3
"""Analyze the focused three-threshold Eigen numerical-error characterization."""
import csv
import json
import math
import sys
from collections import Counter


NAMES = {"radial_0S0", "toroidal_0T2", "toroidal_0T12",
         "spheroidal_0S2", "spheroidal_0S12", "fluid_solid_0S40"}
THRESHOLDS = (1.0, 0.1, 0.01)
REQUIRED = {
    "family", "name", "freq_mhz", "grid_offset", "stencil_member", "response_class",
    "epsilon", "attenuation", "threshold", "n", "factor_info", "solve_info", "finite",
    "min_abs_diag_u", "max_abs_diag_u", "max_abs_u", "max_abs_l", "growth_u_over_a",
    "nonidentity_pivots", "pivot_fraction", "max_pivot_displacement", "residual",
    "solution_vs_threshold1", "solution_abs_vs_threshold1", "transfer_norm",
    "transfer_abs_error_vs_threshold1", "transfer_relative_error_vs_threshold1",
    "local_response_scale", "local_response_normalized_error", "transfer_max_row",
    "transfer_max_col", "transfer_re", "transfer_im", "threshold1_transfer_re",
    "threshold1_transfer_im", "lapack_transfer_re", "lapack_transfer_im",
    "transfer_component_abs_error_vs_threshold1", "lapack_factor_info", "lapack_solve_info",
    "lapack_finite", "lapack_residual", "lapack_solution_vs_threshold1",
    "lapack_solution_abs_vs_threshold1", "lapack_transfer_abs_error_vs_threshold1",
    "lapack_transfer_relative_error_vs_threshold1", "solution_vs_lapack",
    "transfer_abs_error_vs_lapack", "transfer_relative_error_vs_lapack",
}


def f(row, key):
    return float(row[key])


def finite_rows(rows):
    skip = {"family", "name", "response_class"}
    for row in rows:
        for key, value in row.items():
            if key in skip:
                continue
            try:
                if not math.isfinite(float(value)):
                    return False, (key, row)
            except ValueError:
                return False, (key, row)
    return True, None


def dist(rows, key):
    values = sorted(f(row, key) for row in rows)
    return {"minimum": values[0], "median": values[len(values)//2],
            "maximum": values[-1]}


def ident(row):
    return {k: row[k] for k in ("name", "freq_mhz", "grid_offset", "epsilon",
                                "attenuation", "threshold", "response_class")}


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: eigen_pivot_error_characterization_analyze.py RAW SUMMARY")
    with open(sys.argv[1], newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or not REQUIRED.issubset(rows[0]):
        raise SystemExit("raw TSV is empty or missing diagnostics")
    if len(rows) != 396 or {r["name"] for r in rows} != NAMES:
        raise SystemExit("expected 396 rows and six accepted physical systems")
    observed_thresholds = {float(r["threshold"]) for r in rows}
    if observed_thresholds != set(THRESHOLDS):
        raise SystemExit("expected thresholds 1, 0.1, and 0.01")
    if {r["attenuation"] for r in rows} != {"0", "1"}:
        raise SystemExit("expected attenuation on and off")
    if {r["response_class"] for r in rows} - {"resonance", "flank", "antiresonance", "elsewhere"}:
        raise SystemExit("unexpected response class")
    ok, detail = finite_rows(rows)
    if not ok:
        raise SystemExit(f"nonfinite or nonnumeric raw field: {detail}")
    for row in rows:
        for key, expected in (("factor_info", 0), ("solve_info", 0), ("finite", 1),
                              ("lapack_factor_info", 0), ("lapack_solve_info", 0),
                              ("lapack_finite", 1)):
            if int(row[key]) != expected:
                raise SystemExit(f"unexpected {key}: {row}")
    key = ("name", "freq_mhz", "epsilon", "attenuation", "threshold")
    if len({tuple(r[k] for k in key) for r in rows}) != len(rows):
        raise SystemExit("duplicate case/threshold rows")
    def at_threshold(group, threshold):
        return [r for r in group if math.isclose(float(r["threshold"]), threshold,
                                                  rel_tol=0.0, abs_tol=1.e-15)]

    candidate = at_threshold(rows, 0.01)
    middle = at_threshold(rows, 0.1)
    baseline = at_threshold(rows, 1.0)
    pair_key = ("name", "freq_mhz", "epsilon", "attenuation")
    bmap = {tuple(r[k] for k in pair_key): r for r in baseline}
    if len(candidate) != len(middle) or len(candidate) != len(baseline):
        raise SystemExit("threshold sets are not balanced")

    def maxima(group):
        return {k: max(f(r, k) for r in group) for k in (
            "residual", "solution_vs_threshold1", "transfer_abs_error_vs_threshold1",
            "transfer_relative_error_vs_threshold1", "local_response_normalized_error",
            "solution_vs_lapack", "transfer_abs_error_vs_lapack",
            "transfer_relative_error_vs_lapack", "max_abs_u", "max_abs_l",
            "growth_u_over_a")}

    def growth(group):
        result = {k: max(f(r, k) for r in group) for k in (
            "max_abs_diag_u", "max_abs_u", "max_abs_l", "growth_u_over_a",
            "pivot_fraction", "max_pivot_displacement")}
        result["min_abs_diag_u"] = min(f(r, "min_abs_diag_u") for r in group)
        return result

    prior = max(candidate, key=lambda r: f(r, "transfer_relative_error_vs_threshold1"))
    sharp_names = {"radial_0S0", "spheroidal_0S2"}
    peak_shape = {}
    for name in sorted(sharp_names):
        peak_shape[name] = {}
        for ep in sorted({r["epsilon"] for r in rows if r["name"] == name}):
            for att in ("0", "1"):
                case = [r for r in rows if r["name"] == name and r["epsilon"] == ep
                        and r["attenuation"] == att]
                by_threshold = {}
                for threshold in ("1", "0.10000000000000001", "0.01"):
                    selected = at_threshold(case, float(threshold))
                    best = max(selected, key=lambda r: f(r, "transfer_norm"))
                    by_threshold[threshold] = {
                        "peak_offset": int(best["grid_offset"]),
                        "peak_freq_mhz": f(best, "freq_mhz"),
                        "peak_transfer_norm": f(best, "transfer_norm"),
                        "offset_transfer_norm": {r["grid_offset"]: f(r, "transfer_norm")
                                                   for r in selected},
                    }
                peak_shape[name][f"epsilon={ep},attenuation={att}"] = by_threshold
    summary = {
        "decision": "characterization_only_no_production_adoption",
        "rows": len(rows), "systems": sorted(NAMES),
        "threshold_rows": {"1": len(baseline), "0.1": len(middle), "0.01": len(candidate)},
        "stencil_rows": {name: len({r["grid_offset"] for r in rows if r["name"] == name})
                          for name in sorted(NAMES)},
        "response_classes_threshold001": dict(Counter(r["response_class"] for r in candidate)),
        "maxima_threshold_001": maxima(candidate),
        "maxima_threshold_01": maxima(middle),
        "maxima_threshold_1": maxima(baseline),
        "growth_threshold_001": growth(candidate),
        "growth_threshold_01": growth(middle),
        "growth_threshold_1": growth(baseline),
        "worst_threshold001_transfer_relative": ident(prior),
        "worst_threshold001_transfer_values": {
            "transfer_norm": f(prior, "transfer_norm"),
            "local_response_scale": f(prior, "local_response_scale"),
            "absolute_matrix_error": f(prior, "transfer_abs_error_vs_threshold1"),
            "local_scale_normalized_error": f(prior, "local_response_normalized_error"),
            "component_abs_error": f(prior, "transfer_component_abs_error_vs_threshold1"),
            "candidate_complex": [f(prior, "transfer_re"), f(prior, "transfer_im")],
            "threshold1_complex": [f(prior, "threshold1_transfer_re"), f(prior, "threshold1_transfer_im")],
            "lapack_complex": [f(prior, "lapack_transfer_re"), f(prior, "lapack_transfer_im")],
            "component": [int(prior["transfer_max_row"]), int(prior["transfer_max_col"])],
        },
        "peak_shape_and_location_sensitivity": peak_shape,
        "threshold001_vs_threshold1": {
            "max_relative_solution": max(f(r, "solution_vs_threshold1") for r in candidate),
            "max_absolute_solution": max(f(r, "solution_abs_vs_threshold1") for r in candidate),
            "max_absolute_transfer": max(f(r, "transfer_abs_error_vs_threshold1") for r in candidate),
            "max_relative_transfer": max(f(r, "transfer_relative_error_vs_threshold1") for r in candidate),
            "max_local_scale_normalized_transfer": max(f(r, "local_response_normalized_error") for r in candidate),
        },
        "threshold01_vs_threshold1": {
            "max_relative_solution": max(f(r, "solution_vs_threshold1") for r in middle),
            "max_relative_transfer": max(f(r, "transfer_relative_error_vs_threshold1") for r in middle),
            "max_local_scale_normalized_transfer": max(f(r, "local_response_normalized_error") for r in middle),
        },
        "interpretation": {
            "local_response_scale": "maximum threshold-1 receiver/output infinity norm over the same case's production-grid +/-2 stencil; it prevents a near-zero point from making an absolute error look artificially large",
            "complex_transfer": "transfer_re/im are the actual complex receiver/output entry with largest threshold-1 output magnitude in that row; indices identify its receiver and RHS column",
            "peak_shape": "the accepted +/-4 threshold-1 scan selects the local physical peak; +/-2 rows are retained only for radial 0S0 and spheroidal 0S2",
        },
    }
    with open(sys.argv[2], "w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
