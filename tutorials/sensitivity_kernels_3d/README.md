# Windowed 3D material-kernel consumer

This small consumer uses the installed `DSpecM1D::DSpecM1D` package to evaluate the all-five-modulus `(A,C,F,L,N)` spatial kernel for one real scalar objective, `J_Z = Re(sum(conj(q) * P u))`. The single measurement has `q(Z,1)=0.7+0.2i` and `q(Z,2)=-0.25+0.45i`; the raw-bin solver interval is `[1,4)`, `nt=128`, `dt=20 s`, and the time window is Hann on samples 0–90. `run.json` records the resolved grid. Coordinates use geocentric colatitude `theta` and east longitude `phi`, both in radians. Volume kernels are in objective/(Pa m³) and integrate as `K_X deltaX r² dr dOmega`. Source-shell contacts are saved separately in objective/Pa per steradian and integrate over `dOmega`; do not multiply them by `r²`.

The `small` preset uses retained degrees 0–3, solver bins 1–3, and a shell-partitioned homogeneous sphere matching the accepted Phase D geometry. Its objective weights are nonzero at bins 1 and 2; bin 3 remains in the solver interval because window multiplication couples frequencies. The model deck is byte-identical to `build/3d-kernels/phase-d/full-02/common_homogeneous.deck` (SHA-256 `63bedc3e12811707ecf087c38f13a81b5de27643026e63f354239efd24340501`). The `pilot` preset uses degrees 0–12 with the same three-bin band and model. It is a bounded illustration and makes no convergence claim. Metadata records the full objective, source, receiver, Hann window, q weights and their output-derived scale, model units and reference moduli, angular nodes, and physical frequency grid.

The YSpecpp sampled-volume comparison covers matched volume fields away from the source; its sampled-volume API excludes source-shell contact. DSpecM1D also saves a separate source-local split for C, F, and L: a volume perturbation `deltaX = X0 cos²(pi (r-rs)/20000 m) chi(theta,phi)` over `|r-rs| <= 10000 m`, plus the source-shell contact for `deltaX(rs)=X0 chi`. The volume part integrates with `r² dr dOmega`, the contact part with `dOmega`, and the plot shows both signed terms and their sum in objective per unit fractional perturbation. A and N are omitted from this split rather than implying zero source-local volume contributions.

Build and run from an installed package (CMake 3.24+, Python 3, NumPy, and Matplotlib). For timing, point `CXX` at the compiler used for the YSpecpp timing build and build both installed-consumer targets as Release with explicit `-O3` and BLAS disabled:

```sh
cmake -S tutorials/sensitivity_kernels_3d -B /tmp/dspec-3d-consumer \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG' -DGSHTRANS_WITH_BLAS=OFF \
  -DCMAKE_PREFIX_PATH=/path/to/dspecm1d/install
cmake --build /tmp/dspec-3d-consumer \
  --target dspec_3d_kernel_consumer dspec_3d_kernel_timing --parallel 1
python3 tutorials/sensitivity_kernels_3d/run.py \
  --consumer /tmp/dspec-3d-consumer/dspec_3d_kernel_consumer \
  --model tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck \
  --lmax 3 --output /tmp/dspec-3d-small
```

Use `--lmax 12` and a fresh output directory for the bounded pilot. The numerical consumer has a 300-second subprocess limit; if it reaches the limit, the partial output is preserved. Plotting and replay are outside that limit, so saved output can be rendered separately. The runner refuses nonempty destinations before invoking the consumer. It writes full-precision `samples.csv`, `sections.csv`, source contact values and `source_split.csv`, run metadata and SHA-256 provenance, plus `integrals.json` with signed direction-integrated kernels and relative integral errors when YSpecpp data are supplied. It also writes `report.html` and 300 dpi PNG, SVG and PDF plots: per-modulus signed depth maps, great-circle and cross-plane sections, radial-layer globe cutaways, radial profiles, and the source-local volume/contact split. Globe titles identify the sampled-shell range; regions outside it are unsampled. When YSpecpp data are supplied, each modulus depth-map comparison shares a zero-centred scale between solvers and uses a separate labelled signed range for the difference. Sections are evaluated on points on the exact planes. A comparison requires both `--yspec-profiles PATH` and its adjacent timed-run `--yspec-meta PATH`; the runner checks the model, sampled-input and kernel CSV hashes, resolved model/physics/frequency/source/receiver/objective/window configuration, sample multiplicity, and coordinates before plotting. It preserves duplicate interface rows and does not compare source-shell contacts as volume samples. Use `python3 tutorials/sensitivity_kernels_3d/run.py --replay /path/to/previous/run --output /path/to/fresh/plots` to recreate plots without a solver run.

To render a matched comparison from saved data without rerunning either solver:

```sh
python3 tutorials/sensitivity_kernels_3d/run.py \
  --replay /path/to/dspec-run \
  --yspec-profiles /path/to/run_01_yspecpp/kernels.csv \
  --yspec-meta /path/to/run_01_yspecpp/meta.json \
  --output /path/to/fresh-comparison
```

## Matched Earth slices

`earth_slices.py` queries the existing DSpecM1D and YSpecpp timing adapters on
one shared source–receiver great-circle plane. The angular grid covers the full
circle and includes the exact receiver direction; radial samples run from
100 km to the surface and include the exact source and receiver radii. It writes
matched circular cross-sections and oblique 3D cut-plane views for A/C/F/L/N,
with shared DSpecM1D/YSpecpp limits and a separate signed difference scale.
The source and receiver markers use the radii recorded in the source run.

The fixture is a homogeneous sphere, so its cut-plane has no internal interface
contrast. The innermost 100 km is left blank because it is not sampled. These
plots show the computed low-frequency `J_Z` objective kernel; they do not imply
a ray path or a banana–doughnut sensitivity pattern. Previous matched-volume
checks establish DSpecM1D/YSpecpp agreement only over 4,200–6,200 km, where the
validated perturbation vanishes at both ends. Values queried elsewhere are
shown with their raw differences and are explicitly unvalidated; they are not
source-shell contact evidence. The run metadata reports differences by radius
and source distance. The runner checks the saved model hash and the adapters'
resolved source, receiver, moment tensor, physics, frequency grid, objective,
and time window before it accepts their arrays.

With the DSpec and YSpec timing executables built as described above and in the
[3D benchmark README](../../benchmarks/yspecpp/kernels/3d/README.md), generate
and replay the matched slices with:

```sh
python3 tutorials/sensitivity_kernels_3d/earth_slices.py \
  --source-run build/kernel-examples-phase12/comparison-illustration-12-02/source_run.json \
  --dspec-exe /path/to/dspec_3d_kernel_timing \
  --yspec-exe /path/to/yspec_kernel_timing \
  --output build/kernels-plot-redesign/run-01

python3 tutorials/sensitivity_kernels_3d/earth_slices.py \
  --replay build/kernels-plot-redesign/run-01 \
  --output build/kernels-plot-redesign/replay-01
```

The replay command uses saved point coordinates and kernel arrays; it does not
run either solver. `report.html` is the primary review page; it includes linked
PNG, SVG, and PDF figures, provenance, and the all-five-modulus comparison
summary. Numeric inputs, adapter metadata, and SHA-256 hashes are retained in
each fresh output directory. The generated
`report.html` identifies the source/model/window and links every figure export,
input array, kernel array, and full-plane/shell/near-receiver comparison table.
