#!/bin/sh
# Runs every example in turn, each writing its output and figures in its own
# directory.
set -e
cd "$(dirname "$0")"
for example in record-section instrument depths gravity force cross-section \
  map kernels kernels-3d kernels-gallery; do
  echo "== $example"
  sh "$example/run.sh"
done
