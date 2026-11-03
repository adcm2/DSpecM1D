#!/bin/sh
# The wavefield on a great circle through the source and the whole earth,
# from the sea floor to near the centre, at six times: body waves through the
# mantle and core, and surface waves along the top. The computed band includes
# the full 0.1–0.5 mHz lower filter ramp. Each time panel has its own disclosed
# percentile colour scale unless --common is requested.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-snapshot snapshot.json
plot snapshot snapshot.json --component up -o up.png
plot snapshot snapshot.json --component away -o away.png
