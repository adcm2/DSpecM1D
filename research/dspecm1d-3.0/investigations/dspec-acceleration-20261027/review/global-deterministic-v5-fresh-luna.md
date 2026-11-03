# Fresh static review: combined global v5 deterministic FFT candidate

**Result: PASS for bounded source/config review; eligible for the separately required Sol review and build phase.** This review did not run a build, connect by SSH, or run science.

The v4→v5 implementation delta is limited to changing `FFTWpp::Measure` to `FFTWpp::Estimate` when constructing `CompleteBandGrid`, plus recording that planner policy in the output manifest. The backend diagnostic remains explicit: `GSHTrans::TransformKernel::Loop` with FFT longitude synthesis; Eigen GEMM is not used for angular synthesis. Local pinned FFTW++ `Options.h` declares both `Estimate` and `Measure` as valid `Flag` values (`FFTW_ESTIMATE` and `FFTW_MEASURE` respectively), confirming the selected symbol and constructor argument type.

The deterministic build configuration selects one `global-profiled` variant from the staged `combined-global-v5` source into the new `kernel-acceleration-20261027-global-deterministic-01` build root. Its build plan has one build job and an 8-thread cap. The benchmark v2 plan has 23 jobs and a 12-thread cap. Compared with the prior benchmark plan, all 23 job IDs and the full normalized plan are identical after substituting only the fresh staging directory and the new deterministic binary path. Thus the benchmark inputs, job arguments, resource limits, and acceptance checks are unchanged; staging and executable paths are fresh.

This changes only FFT planner policy and its diagnostic. It does not relax the strict byte comparison gate or alter cross-route acceptance criteria. The motivation recorded for the new candidate is the prior v4 L128 parallel byte mismatch (5.42e-16) with exact waveform-energy agreement; that observation is contextual and was not independently rerun here. No speedup or scientific result is established by this static review.

## SHA-256

- `candidates/combined-global-v5/global/global_gsh_tile.hpp`: `064b2783a3f8840dc8309b40a413bb23c18a8d884fd862a973e895377a1b2286`
- `candidates/combined-global-v5/global/global_kernel_consumer.hpp`: `1478b3629f6e8868c9cf4735c0bfb8d07ef88093b52071d8f22277611874af4f`
- `candidates/combined-global-v5/implementation_status.md`: `2837f5ec2376539e82514b961269eb535068f599d77f869a637f5af1d285ac2c`
- `configs/global-deterministic-build.json`: `0e08f4ec62ce1746d0c1425221ebc5bb70a7eb517c64b660b3692fd5f7b1fa67`
- `configs/global-deterministic-build-plan.json`: `2931e58fbc9b9decf14174c7ab123b10a37eb74f717f8d4b8c42b59ab232197a`
- `configs/global-benchmark-v2-plan.json`: `70e90f5bfef004233144054eb2668739ac392b1d79190eea1d477613b393296f`
