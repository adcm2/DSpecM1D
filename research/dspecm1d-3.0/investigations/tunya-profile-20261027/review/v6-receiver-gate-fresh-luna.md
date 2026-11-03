# Fresh Luna review: v6 receiver gate

**Verdict: PASS by static inspection.** The active gate is suitable for a reviewed Tunya build/run. This verdict does not claim successful compilation or numerical execution; neither was performed.

- The generated source has 17 syntactically complete `#include` lines. The v4 and v6 header SHA-256 values match their frozen manifests and the active source bindings (`ea276794…d8ff`, `4d054875…b6abca9`). The four extracted bodies—v4 `harmonic`, v4/v6 `receiverAtPole`, and v6 `angularWignerBytesPerAngle`—are verbatim substrings of the corresponding hash-verified headers.
- The oracle independently encodes all six S/T entries from `Y^L_nm(0,0)=rho delta(n,m)`: S has `(+rho)` at radial `(row 0,m=0)`, `(-rho,+rho)` in row 1 at `m=(-1,+1)`, and `(-i rho,-i rho)` in row 2; T has `(+i rho,+i rho)` in row 1 and `(-rho,+rho)` in row 2. This matches the source channel signs.
- For each tested degree, the old harmonic helper is checked against every degree-valid pair with `n,m` in `[-2,2]`. Nonzero components use an eight-ulp gate. V4 forbidden entries and zero phase components use the documented `5e-15 rho` bound; v6 forbidden entries and zero phase components require exact zero. The additional scaled absolute-error and old/new gates are present.
- The payload check independently sums `GSHIndices<All>(L,mMax,n).Size()` for capped and full `mMax`, then compares exact byte totals with the v6 estimator. It requires the full estimate to exceed the capped estimate when `L>2`; the test allocates no Wigner table.
- The source is plausible C++23 by inspection: headers cover the used standard/Eigen/GSHTrans types and algorithms; `std::bit_cast` is enabled by the campaign's `-std=gnu++23`; no malformed preprocessor lines or unresolved template tokens remain.

No blocker found for the standalone gate. Compilation and execution remain required on Tunya before interpreting any numerical or timing result.
