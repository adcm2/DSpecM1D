#!/usr/bin/env python3
"""Validate and summarize the focused sequential SuperLU factorization probe."""

import csv
import json
import math
import pathlib
import sys

EXPECTED = {
    "radial_small", "toroidal_small", "spheroidal_small",
    "spheroidal_medium", "spheroidal_large",
    "spheroidal_fluid_solid", "spheroidal_max",
}


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise SystemExit(f"non-finite {key} in {row['system']}/{row['freq_index']}")
    return value


def main():
    if len(sys.argv) != 4:
        raise SystemExit("usage: superlu_factorization_probe_analyze.py RAW SUMMARY_JSON SUMMARY_TSV")
    raw, summary_json, summary_tsv = map(pathlib.Path, sys.argv[1:])
    with raw.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 14 or {r["system"] for r in rows} != EXPECTED:
        raise SystemExit("expected two rows for each of seven SEM systems")
    if any(sorted(int(r["freq_index"]) for r in rows if r["system"] == name) != [0, 1]
           for name in EXPECTED):
        raise SystemExit("frequency-index coverage is incomplete")
    integer_zero = ["eigen_factor_info", "eigen_solve_info", "superlu_dofact_info",
                    "superlu_factor_info", "superlu_solve_info",
                    "lapack_factor_info", "lapack_solve_info"]
    numeric = [
        "freq_mhz", "eigen_setup_s", "eigen_factor_s", "eigen_solve_s",
        "superlu_conversion_s", "superlu_setup_s", "superlu_dofact_s",
        "superlu_samepattern_s", "superlu_solve_s", "eigen_residual",
        "superlu_residual", "lapack_residual", "solution_discrepancy",
        "eigen_vs_lapack", "superlu_vs_lapack",
    ]
    for row in rows:
        if int(row["repeats"]) != 16 or int(row["finite"]) != 1:
            raise SystemExit("unexpected repeat count or non-finite result")
        if any(int(row[key]) != 0 for key in integer_zero):
            raise SystemExit(f"solver failure in {row['system']}/{row['freq_index']}")
        for key in numeric:
            number(row, key)
        if max(number(row, key) for key in
               ("eigen_residual", "superlu_residual", "lapack_residual")) > 1e-14:
            raise SystemExit("backward residual exceeds validation tolerance")
        if number(row, "solution_discrepancy") > 1e-6:
            raise SystemExit("solution discrepancy exceeds validation tolerance")

    per_system = []
    for name in sorted(EXPECTED):
        selected = [r for r in rows if r["system"] == name]
        eigen = sum(number(r, "eigen_factor_s") for r in selected)
        superlu = sum(number(r, "superlu_samepattern_s") for r in selected)
        per_system.append({
            "system": name,
            "n": int(selected[0]["n"]),
            "eigen_factor_seconds_sum": eigen,
            "superlu_samepattern_seconds_sum": superlu,
            "eigen_over_superlu_factor_speedup": eigen / superlu,
            "eigen_over_superlu_solve_speedup":
                sum(number(r, "eigen_solve_s") for r in selected) /
                sum(number(r, "superlu_solve_s") for r in selected),
        })
    eigen_total = sum(number(r, "eigen_factor_s") for r in rows)
    superlu_total = sum(number(r, "superlu_samepattern_s") for r in rows)
    payload = {
        "experiment": "superlu_factorization_probe",
        "decision": "B",
        "decision_text": "SuperLU SamePattern is not materially faster overall; the added dependency is not justified.",
        "rows": len(rows),
        "systems": len(EXPECTED),
        "frequencies_per_system": 2,
        "repeats": 16,
        "aggregate": {
            "eigen_factor_seconds_sum": eigen_total,
            "superlu_samepattern_seconds_sum": superlu_total,
            "eigen_over_superlu_factor_speedup": eigen_total / superlu_total,
            "max_eigen_residual": max(number(r, "eigen_residual") for r in rows),
            "max_superlu_residual": max(number(r, "superlu_residual") for r in rows),
            "max_lapack_residual": max(number(r, "lapack_residual") for r in rows),
            "max_solution_discrepancy": max(number(r, "solution_discrepancy") for r in rows),
            "max_eigen_vs_lapack": max(number(r, "eigen_vs_lapack") for r in rows),
            "max_superlu_vs_lapack": max(number(r, "superlu_vs_lapack") for r in rows),
        },
        "per_system": per_system,
        "concurrency_probe": "not run because the serial SamePattern gate did not show a material factorization gain",
        "production_spectra": "not run",
    }
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with summary_tsv.open("w", newline="") as stream:
        fields = list(per_system[0])
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader(); writer.writerows(per_system)
    print(json.dumps(payload["aggregate"], sort_keys=True))


if __name__ == "__main__":
    main()
