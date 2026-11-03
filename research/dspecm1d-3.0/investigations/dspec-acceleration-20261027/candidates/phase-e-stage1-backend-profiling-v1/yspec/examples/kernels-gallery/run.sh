#!/bin/sh
# A gallery of three-dimensional kernels for cross-correlation windows
# targeted on PREM TauP ray times (make_input.py wrote the inputs): at 60
# degrees, P and PP on Z and S, ScS and SS on T; at 110 degrees, Pdiff on Z
# and Sdiff on T. Window labels identify the targeted ray arrivals; depth
# phases and other branches can share a finite-frequency packet, so no
# single-phase isolation is claimed. The broad 8–40 mHz and longer-period
# 5–25 mHz bands serve different packet demonstrations; neither replaces the
# separately specified Bolivia comparison setup. Each run
# computes its kernels from one forward field and an adjoint field for each
# component; each measurement has its own kernels, drawn on sections along
# the path and across it at its middle, below its seismogram and window.
# Runtime depends on the build, degree and adjoint fields requested.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels gallery.json
plot gallery gallery.json -o gallery.png
plot banana gallery.json --measurement S -o banana-S.png
run yspec-kernels diffracted.json
plot gallery diffracted.json -o diffracted.png
