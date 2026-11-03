#!/bin/sh
# A point force, of the kind a landslide exerts: 1e11 N to the north, 5 km
# deep, recorded along a great circle 30 degrees east of north. Its radiation
# is that of a single force: Rayleigh waves on the vertical and radial
# components, Love waves on the transverse.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec input.json
plot section input.json --component Z -o vertical.png
plot section input.json --component T -o transverse.png
