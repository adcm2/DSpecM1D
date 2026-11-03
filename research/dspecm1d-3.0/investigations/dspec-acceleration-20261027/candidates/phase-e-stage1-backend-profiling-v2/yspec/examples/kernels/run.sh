#!/bin/sh
# Radial kernels for two named cross-correlation windows at 60 degrees. The
# transverse 1000–1200 s window contains direct S (about 1094 s), sS (about
# 1109 s) and ScS (about 1191 s): it is a mixed packet, not an isolated S
# phase. The long-period vertical window is a packet target, not a validated
# Rayleigh measurement. The windows and resulting kernels remain
# scientifically unresolved. The two share the receiver and are computed
# together.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels traveltimes.json
plot kernels traveltimes.json --measurement S -o s-wave.png
plot kernels traveltimes.json --measurement long-period-Z-packet --deepest 700 \
  -o rayleigh.png
