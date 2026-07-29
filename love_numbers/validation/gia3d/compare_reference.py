#!/usr/bin/env python3

import math
import subprocess
import sys


COMPONENTS = ("h_u", "k_u", "h_phi", "k_phi", "h_t", "k_t")
OUTPUT_COLUMNS = (1, 2, 3, 4, 7, 8)


def read_reference(path):
    reference = {}
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.split()
            reference[int(fields[0])] = tuple(map(float, fields[1:]))
    return reference


def relative_error(calculated, reference):
    if reference != 0.0:
        return abs((calculated - reference) / reference)
    return 0.0 if calculated == reference else math.inf


def main():
    executable, reference_path, *degrees = sys.argv[1:]
    completed = subprocess.run(
        [executable, *degrees],
        check=False,
        capture_output=True,
        text=True,
    )
    sys.stderr.write(completed.stderr)
    if completed.returncode != 0:
        return completed.returncode

    diagnostics = {}
    for line in completed.stderr.splitlines():
        fields = line.split()
        if not fields or fields[0] != "diagnostic":
            continue
        degree = int(fields[1])
        values = tuple(map(float, fields[2:]))
        if len(values) != 8 or not all(map(math.isfinite, values)):
            raise RuntimeError(f"invalid diagnostic row: {line}")
        diagnostics[degree] = values

    calculated = {}
    for line in completed.stdout.splitlines():
        fields = line.split()
        degree = int(fields[0])
        values = tuple(map(float, fields[1:]))
        if len(values) != 8 or not all(map(math.isfinite, values)):
            raise RuntimeError(f"invalid driver row: {line}")
        calculated[degree] = values
        print(line)

    reference = read_reference(reference_path)
    for degree_text in degrees:
        degree = int(degree_text)
        if degree not in diagnostics:
            raise RuntimeError(f"missing diagnostic row for degree {degree}")
        values = calculated[degree]
        expected = reference[degree]
        for name, column, published in zip(
            COMPONENTS, OUTPUT_COLUMNS, expected
        ):
            value = values[column - 1]
            difference = value - published
            error = relative_error(value, published)
            print(
                f"published l={degree} component={name} "
                f"calculated={value:.17e} reference={published:.17e} "
                f"signed_difference={difference:.17e} "
                f"relative_error={error:.17e}"
            )

    degree_one = calculated.get(1)
    if degree_one is not None:
        h_phi = degree_one[2]
        k_u = degree_one[1]
        k_phi = degree_one[3]
        k_t = degree_one[7]
        if (h_phi, k_u, k_phi, k_t) != (0.0, 0.0, 0.0, 0.0):
            raise RuntimeError("degree-one constrained responses are not zero")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
