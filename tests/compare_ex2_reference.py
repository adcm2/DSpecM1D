#!/usr/bin/env python3

import sys
import numpy as np


def summarize(expected_col: np.ndarray, actual_col: np.ndarray) -> tuple[float, float]:
    norm = float(np.max(np.abs(expected_col))) + 1e-12
    rel = np.abs(actual_col - expected_col) / norm * 100.0
    return float(np.mean(rel)), float(np.max(rel))


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: compare_ex2_reference.py EXPECTED ACTUAL", file=sys.stderr)
        return 2

    expected = np.loadtxt(sys.argv[1], delimiter=";")
    actual = np.loadtxt(sys.argv[2], delimiter=";")

    if expected.shape != actual.shape:
        print(f"shape mismatch: expected {expected.shape}, actual {actual.shape}", file=sys.stderr)
        return 1

    # Columns: time, DSpecM1D Z/N/E, YSpec Z/N/E, MINEOS Z/N/E, specnm Z/N/E
    component_names = ["Z", "N", "E"]
    component_cols = [1, 2, 3]

    avg_limit = 0.50
    max_limit = 5.00

    failed = False
    for name, col in zip(component_names, component_cols):
        avg_diff, max_diff = summarize(expected[:, col], actual[:, col])
        print(f"{name}: average={avg_diff:.4f}% max={max_diff:.4f}%")
        if avg_diff > avg_limit or max_diff > max_limit:
            failed = True

    if failed:
        print(
            f"ex2 migration check exceeded tolerance "
            f"(avg>{avg_limit:.2f}% or max>{max_limit:.2f}%)",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
