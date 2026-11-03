#!/bin/sh
# One place, recorded twice: the ground displacement, and the output in counts
# of a broadband seismometer, which responds to velocity and whose response
# falls away at periods beyond its 120 s.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec input.json
plot traces input.json --receivers 1 2 --components ZRT -o station.png
