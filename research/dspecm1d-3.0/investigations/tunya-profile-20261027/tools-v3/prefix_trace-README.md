# Trace-prefix adapter

`prefix_trace.py INPUT_FULL_CSV OUTPUT_PREFIX_CSV RECEIPT_JSON` validates the
fixed DSpec trace contract (header `time_s,vertical_m,window_weight`, 1024 rows,
5 s grid, finite values, and Hann weights for 1000–1200 s). It requires every
removed row 721–1023 to carry an exact numeric zero weight, verifies identical
positive `math.fsum` Hann-weighted energies for full and prefix inputs using
`0.5*dt*sum(weight*vertical_m^2)` (units m² s), and writes the original header
plus rows 0–720 byte-for-byte. It refuses either
existing output, never edits the input, and records before/after input hashes,
output hash, counts, grid/window range, and energy in a fresh receipt.
If an error occurs after output creation, the partial output and any receipt
artifact are preserved as evidence; retries must use fresh paths.

`test_prefix_trace.py` uses only synthetic temporary CSVs; it checks exact-byte
prefix preservation, removed-row zero enforcement, header/grid/finiteness
rejection, and refusal to replace an output. It does not invoke a solver or
modify campaign data.
