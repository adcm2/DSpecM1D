#!/bin/sh
# One place, 40 degrees from the source, at depths from the sea floor through
# the mantle to the outer and inner cores. In the fluid outer core there is no
# toroidal motion, but the horizontal motion still has a transverse part
# wherever the field varies with azimuth.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec input.json
plot traces input.json --components ZRT -o depths.png
