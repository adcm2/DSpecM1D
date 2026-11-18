#!/usr/bin/env python3
"""Validate finite, nonzero multi-receiver pilot waveforms."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def check(path: Path, expected_receivers: int = 3) -> dict[str, object]:
    rows: list[list[float]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter=";"), 1):
            if len(row) != 1 + 3 * expected_receivers:
                raise ValueError(f"{path}:{line_number}: expected 10 columns, got {len(row)}")
            values = [float(value) for value in row]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path}:{line_number}: non-finite value")
            rows.append(values)
    if len(rows) < 2:
        raise ValueError(f"{path}: fewer than two samples")
    times = [row[0] for row in rows]
    if not math.isclose(times[0], 0.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(f"{path}: first time is not zero")
    intervals = [b - a for a, b in zip(times, times[1:])]
    if any(not math.isclose(value, 2.0, rel_tol=0.0, abs_tol=1e-9) for value in intervals):
        raise ValueError(f"{path}: sample interval is not 2 seconds")
    if times[-1] > 602.0 + 1e-9:
        raise ValueError(f"{path}: output extends beyond the requested 600 seconds")
    peaks = [max(abs(row[column]) for row in rows) for column in range(1, len(rows[0]))]
    family = "toroidal" if "toroidal" in path.name else "spheroidal"
    for receiver in range(expected_receivers):
        vertical = peaks[3 * receiver]
        horizontal = peaks[3 * receiver + 1:3 * receiver + 3]
        if max(horizontal) == 0.0:
            raise ValueError(f"{path}: receiver {receiver} has no horizontal excitation")
        if family == "toroidal" and vertical > 1e-8 * max(horizontal):
            raise ValueError(f"{path}: receiver {receiver} toroidal vertical is not negligible")
        if family == "spheroidal" and vertical == 0.0:
            raise ValueError(f"{path}: receiver {receiver} spheroidal vertical is identically zero")
    return {"rows": len(rows), "columns": len(rows[0]), "time_end_s": times[-1], "peaks": peaks}


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit("usage: check_outputs.py OUTPUT [OUTPUT ...]")
    for argument in sys.argv[1:]:
        path = Path(argument)
        summary = check(path)
        print(f"PASS {path}: {summary['rows']} rows, {summary['columns']} columns, end={summary['time_end_s']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
