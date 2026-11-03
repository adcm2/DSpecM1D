# Global GSH physical coefficient probe

This isolated compile-only prototype is based on frozen DSpec v4 headers. It extracts one frequency's source and receiver-adjoint pole coefficients at one selected native SEM side, rotates each of the six complex strain-spin channels into the global frame, and synthesizes them over the complete native `ForBand(L,2,2.0)` angular point set with a band-`L` `SphericalGrid`, generated Wigner values, and the loop transform.

The custom grid uses the same `2L+1` Gauss-Legendre colatitudes and `FastFFTSize(4L+1)` longitudes as the complete-2L native rule while declaring only the field band `L`. `result.json` records the low-degree exact point/weight parity check, rotation convention, resource estimates, field layout, and per-stage wall times. `gsh_fields.bin` contains twelve channel-major complex128 arrays, forward channels 0–5 followed by receiver-adjoint channels 0–5. With `--reference on`, the program evaluates the accepted complex tensor sampler on the complete radial-by-angle tensor grid and writes the same layout to `reference.bin`; the diagnostic also records channel and pointwise modulus-kernel residuals.

For a resource-bounded high-degree pass, omit the reference route. Low-degree equivalence runs should set `--reference on`, use L=3, 6, and 12, and select multiple representative native sides in separate fresh output directories. No run has been made by the implementation agent; the parent owns validation and physical execution.

Example interface:

```text
physical_global_gsh MODEL EMPTY_OUTPUT --lmax 12 --nq 6 --radial-side-index 0 --radial-refinement 0.05 --qex 1 --frequency-bin 2 --reference on --max-memory-gib 8 --max-output-gib 8
```

All outputs are refused when `EMPTY_OUTPUT` is nonempty. Preflight reports an estimate of the driver's known arrays and available system RAM/disk; this is not a process peak-RSS guarantee and excludes solver and transform internals.
