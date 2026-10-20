#!/usr/bin/env python3
"""Compare two public semicolon-delimited seismograms with a small tolerance."""
from __future__ import annotations

import sys
import numpy as np


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_seismogram_outputs.py REFERENCE CANDIDATE")
    reference = np.loadtxt(sys.argv[1], delimiter=";")
    candidate = np.loadtxt(sys.argv[2], delimiter=";")
    if reference.shape != candidate.shape:
        raise SystemExit(f"shape mismatch: {reference.shape} != {candidate.shape}")
    if not np.isfinite(reference).all() or not np.isfinite(candidate).all():
        raise SystemExit("comparison contains non-finite values")
    if not np.allclose(reference, candidate, rtol=1e-12, atol=1e-15):
        difference = np.max(np.abs(reference - candidate))
        raise SystemExit(f"outputs differ beyond rtol=1e-12, atol=1e-15; max abs diff {difference:g}")
    print(f"Outputs agree within rtol=1e-12, atol=1e-15: {reference.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
