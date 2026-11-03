# Forward v1 checkpoint provenance

`physical_driver-profiled-forward-v1` remains preserved with binary SHA-256 `edc7270cab7da718a46d8328a973cdbcdce78413525a9d4adc531488f3fd1bf1`. The contemporaneous `freeze-forward-v1.sha256` records the source/header hashes at compilation time. The exact v1 source tree was not copied into a versioned snapshot before v2 edits, so those sources are no longer retained and the v1 executable is not reproducible from the current working tree. Treat v1 as an exploratory binary checkpoint only; final comparisons should use a frozen-source build.
