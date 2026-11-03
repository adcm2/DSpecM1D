# Sourced by the scripts that run the examples, from their own directories.
#
# The programs are found in YSPEC_BIN if it is set, else in the build tree
# beside the examples (they are copied to build/examples), else on the PATH.
# The calculations use four threads unless OMP_NUM_THREADS says otherwise.

# In the build tree, where the sources are.
if [ -f ../source.sh ]; then
  . ../source.sh
fi

: "${OMP_NUM_THREADS:=4}"
export OMP_NUM_THREADS

program() {
  for directory in "${YSPEC_BIN:-}" ../../bin ../../build/bin; do
    if [ -n "$directory" ] && [ -x "$directory/$1" ]; then
      echo "$directory/$1"
      return
    fi
  done
  if command -v "$1" > /dev/null; then
    command -v "$1"
    return
  fi
  echo "cannot find $1: build YSpec, or set YSPEC_BIN to where it is" >&2
  exit 1
}

# Runs a program on an input file.
run() {
  name=$1
  shift
  executable=$(program "$name")
  echo "$name $*"
  "$executable" "$@"
}

# Plots with yspec-plot, which the Python package installs with its plot
# extra: from the PATH, in an activated environment, say.
plot() {
  if command -v yspec-plot > /dev/null; then
    yspec-plot "$@"
  else
    echo "not plotting: install the Python package with its plotting," \
      "pip install -e '${YSPEC_SOURCE_DIR:-../..}[plot]', in an environment" \
      "on the PATH" >&2
  fi
}
