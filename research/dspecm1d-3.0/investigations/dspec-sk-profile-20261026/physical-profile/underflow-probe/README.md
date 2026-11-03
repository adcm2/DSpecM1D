# Optional FTZ/DAZ profiling preload

This isolated `LD_PRELOAD` diagnostic is intended only to test whether gradual underflow materially affects profiling speed. It does not edit DSpec or either numerical prototype. The constructor is inert unless `DSPEC_PROFILE_FTZ=1` is present. In that mode, on x86 with SSE/SSE2 it enables MXCSR FTZ (bit 15) and DAZ (bit 6), then logs the before/after values to standard error. On unsupported architectures or CPUs it exits with status 78 and a clear diagnostic. Any run using it changes handling of subnormal values and must be compared against the unmodified same binary for both timing and numerical differences before interpretation.

Example invocation for the parent-owned diagnostic only:

```sh
DSPEC_PROFILE_FTZ=1 LD_PRELOAD=/absolute/path/libdspec_profile_ftz.so <same frozen benchmark command>
```

No executable or preload test is run by the implementation agent. Build: `python3 investigations/dspec-sk-profile-20261026/algorithm-audit/run_compile.py -- /usr/bin/cc -std=c11 -O2 -Wall -Wextra -Werror -fPIC -shared ...`.
