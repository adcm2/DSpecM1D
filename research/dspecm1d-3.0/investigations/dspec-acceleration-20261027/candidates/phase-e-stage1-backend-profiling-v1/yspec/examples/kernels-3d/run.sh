#!/bin/sh
# A three-dimensional kernel for the named transverse cross-correlation
# window at 60 degrees (1000–1200 s). Direct S, sS and ScS all fall in that
# window, so this is a mixed packet rather than an isolated S sensitivity.
# The relative isotropic speed and inverse-Q fields at fixed density are
# drawn on an along-path section and on maps at four depths. Kernel and grid
# convergence remain unverified. The kernels are written as the planetmodel
# file s-wave.nc; runtime depends on the requested mesh and build.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels s-wave.json
plot kernels s-wave.json -o s-wave.png
plot volume s-wave.json --parameter vs -o s-wave-volume.png
