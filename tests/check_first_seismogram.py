#!/usr/bin/env python3
"""Check the public single-receiver time-series output without extra packages."""
import math
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_first_seismogram.py OUTPUT")
    path = Path(sys.argv[1])
    rows = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.rstrip("\n").split(";")
            if len(fields) != 4:
                raise SystemExit(f"{path}:{line_number}: expected four columns, found {len(fields)}")
            try:
                row = [float(field) for field in fields]
            except ValueError as error:
                raise SystemExit(f"{path}:{line_number}: invalid numeric value: {error}")
            if not all(math.isfinite(value) for value in row):
                raise SystemExit(f"{path}:{line_number}: output contains non-finite values")
            rows.append(row)
    if len(rows) < 2:
        raise SystemExit(f"{path}: expected at least two data rows")
    times = [row[0] for row in rows]
    if len(rows) != 152 or times[0] != 0.0 or times[-1] != 302.0:
        raise SystemExit(f"{path}: expected 152 samples from 0 through 302 seconds")
    intervals = [later - earlier for earlier, later in zip(times, times[1:])]
    if any(not math.isclose(interval, 2.0, rel_tol=0.0, abs_tol=1e-9) for interval in intervals):
        raise SystemExit(f"{path}: expected a constant 2-second sample interval")
    if any(not any(row[column] != 0.0 for row in rows) for column in (1, 2, 3)):
        raise SystemExit(f"{path}: one or more receiver components are identically zero")
    print(f"Checked {len(rows)} finite rows with three nonzero receiver components: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
