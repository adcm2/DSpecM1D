#!/bin/sh
# The same long-period seismogram without gravity, in the Cowling
# approximation, and with full self-gravitation, which matters most at the
# lowest frequencies.
set -e
cd "$(dirname "$0")"
. ../common.sh
for gravity in none cowling full; do
  run yspec "$gravity.json"
done
plot compare none.json cowling.json full.json --labels none Cowling full \
  --component Z -o gravity.png
