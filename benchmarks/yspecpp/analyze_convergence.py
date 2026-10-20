#!/usr/bin/env python3
"""Compare adjacent coarse/medium/fine runs without fitting or resampling."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from run_benchmark import assert_close, compare_vectors, read_meta, read_spectra, read_timeseries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--medium", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    roots = {"coarse": args.coarse, "medium": args.medium, "fine": args.fine}
    args.output.mkdir(parents=True, exist_ok=True)
    metrics = []
    settings = {}
    physical_contract = None
    source_contract = None

    for refinement, root in roots.items():
        summary = json.loads((root / "summary.json").read_text())
        try:
            case = next(case for case in summary["cases"]
                        if case["id"] == args.case)
        except StopIteration:
            raise SystemExit(f"{root}: case {args.case} is absent")
        defaults = json.loads(json.dumps(summary["defaults"]))
        defaults.pop("solver", None)
        provenance = summary["provenance"]
        physical = {
            "case": case,
            "defaults": defaults,
            "fixture_sha256": provenance["fixture_sha256"],
        }
        sources = {key: provenance[key] for key in (
            "dspec_source_head", "dspec_tracked_source_files_changed",
            "dspec_tracked_source_diff_sha256", "yspecpp_commit",
            "yspecpp_filetree_sha256", "yspecpp_source_archive_sha256",
            "dependency_shas", "planetmodel_sha")}
        if physical_contract is None:
            physical_contract = physical
            source_contract = sources
        elif physical != physical_contract:
            raise SystemExit(f"{root}: physical case/model/source/output settings differ")
        elif sources != source_contract:
            raise SystemExit(f"{root}: solver source or dependency provenance differs")
        if provenance.get("refinement") != refinement:
            raise SystemExit(f"{root}: expected {refinement} refinement, found "
                             f"{provenance.get('refinement')}")
        settings[refinement] = {
            "provenance": provenance,
            "defaults": summary["defaults"],
            "case": case,
        }

    for coarse_name, fine_name in (("coarse", "medium"), ("medium", "fine")):
        for solver in ("dspec", "yspecpp"):
            for mode in ("all", "l0", "l1", "l2", "l3", "l4"):
                left_dir = roots[coarse_name] / args.case
                right_dir = roots[fine_name] / args.case
                left_prefix = left_dir / f"{solver}-{mode}"
                right_prefix = right_dir / f"{solver}-{mode}"
                left_meta = read_meta(left_prefix.with_suffix(".meta.json"))
                right_meta = read_meta(right_prefix.with_suffix(".meta.json"))
                for field in ("samples", "output_samples"):
                    if left_meta[field] != right_meta[field]:
                        raise SystemExit(f"{solver} {mode}: {field} differs across refinements")
                for field in ("frequency_spacing_hz", "damping_epsilon_rad_s",
                              "sample_interval_s"):
                    assert_close(f"{solver} {mode}: {field}",
                                 left_meta[field], right_meta[field])

                left, left_points = read_spectra(left_prefix.with_suffix(".spectra.csv"))
                right, right_points = read_spectra(right_prefix.with_suffix(".spectra.csv"))
                keys = sorted(left.keys() & right.keys())
                if len(keys) != len(left) or len(keys) != len(right):
                    raise SystemExit(f"{solver} {mode}: spectral bins differ")
                for key in keys:
                    for field in ("frequency_hz", "omega_real_rad_s",
                                  "omega_imag_rad_s"):
                        assert_close(f"{solver} {mode} bin {key} {field}",
                                     left_points[key][field],
                                     right_points[key][field],
                                     1e-10 if field == "omega_real_rad_s" else 1e-12)
                for component in range(3):
                    selected = [key for key in keys if key[1] == component]
                    result = compare_vectors([right[key] for key in selected],
                                             [left[key] for key in selected])
                    metrics.append({
                        "family": "spectrum", "solver": solver,
                        "case": args.case, "mode": mode,
                        "component": component,
                        "refinement_pair": f"{coarse_name}_to_{fine_name}",
                        "relative_l2": result["relative_l2"],
                        "max_abs_difference": result["max_abs_difference"],
                        "max_difference_over_finer_peak":
                            result["max_difference_over_reference_peak"],
                    })

                left_time, left_series = read_timeseries(left_prefix.with_suffix(".timeseries.csv"))
                right_time, right_series = read_timeseries(right_prefix.with_suffix(".timeseries.csv"))
                if len(left_time) != len(right_time):
                    raise SystemExit(f"{solver} {mode}: time sample counts differ")
                for index, (a, b) in enumerate(zip(left_time, right_time)):
                    assert_close(f"{solver} {mode} time sample {index}", a, b)
                for component, name in enumerate(("Z", "N", "E")):
                    result = compare_vectors(
                        [complex(value) for value in right_series[component]],
                        [complex(value) for value in left_series[component]])
                    metrics.append({
                        "family": "time_series", "solver": solver,
                        "case": args.case, "mode": mode,
                        "component": name,
                        "refinement_pair": f"{coarse_name}_to_{fine_name}",
                        "relative_l2": result["relative_l2"],
                        "max_abs_difference": result["max_abs_difference"],
                        "max_difference_over_finer_peak":
                            result["max_difference_over_reference_peak"],
                    })

    (args.output / "convergence.json").write_text(json.dumps({
        "settings": settings,
        "metrics": metrics,
    }, indent=2) + "\n")
    with (args.output / "convergence.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    print(f"Convergence results: {args.output.resolve()}")


if __name__ == "__main__":
    main()
