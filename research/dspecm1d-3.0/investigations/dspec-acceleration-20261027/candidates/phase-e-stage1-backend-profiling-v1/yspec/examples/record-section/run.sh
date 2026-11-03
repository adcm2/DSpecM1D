#!/bin/sh
# A record section: velocity at the sea floor every 5 degrees along a great
# circle from the source, to periods of 40 s. The vertical and transverse
# components, each trace scaled to its peak and amplified twice, clipped
# where it would reach the next, to bring out the body waves.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec input.json
plot section input.json --component Z --gain 2 -o vertical.png
plot section input.json --component T --gain 2 -o transverse.png
