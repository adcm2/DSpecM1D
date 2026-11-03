#!/usr/bin/env bash
# Regenerates the reference seismograms in tests/regression/*/reference from
# the Fortran YSpec, for the four cases of the Fortran's own test suite.
#
#   tests/fortran/make_regression_reference.sh /path/to/YSpec [commit]
#
# The Fortran is exported at the given commit (default a3ae638, the cleaned
# code) and changed in three ways, so that what is compared is the physics and
# not the Fortran's known faults:
#   - its step density is raised sixteen-fold (nw = 4 -> 64), so that it is
#     accurate to about 1e-6 of the peak instead of 1e-3;
#   - rad_derivs_ng no longer keeps the gravity term of rad_amat, which with
#     gravity switched off left l = 0 with gravity;
#   - rad_solver_ng starts its downward integration at the surface in both
#     branches; with the receiver below the source it started at the sea
#     floor, leaving the ocean out at l = 0.
# The inputs, tests/regression/*/input.json, are the JSON equivalents of the
# Fortran's tests/*/input.txt and are kept by hand.
set -euo pipefail
yspec=$(realpath "$1")
commit=${2:-a3ae638}
here=$(cd "$(dirname "$0")" && pwd)
root="$here/../.."
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

git -C "$yspec" archive "$commit" src tests examples | tar -x -C "$work"
cd "$work/src"
sed -i 's/parameter, private :: nw = 4 /parameter, private :: nw = 64 /' module_int.f90
python3 - <<'PY'
s = open('module_int.f90').read()
start = s.index('    subroutine rad_derivs_ng(rr,y,dydx)')
end = s.index('end subroutine rad_derivs_ng')
s = s[:start] + s[start:end].replace(
    'call rad_amat(rr,t,k,s)',
    'call rad_amat(rr,t,k,s)\n      s = s+4.0_dp*spline_value(r,rho,rho_cs,rr,il)'
    '*spline_value(r,grav,grav_cs,rr,il)/rr') + s[end:]
s = s.replace('    subroutine rad_derivs_ng(rr,y,dydx)\n',
              '    subroutine rad_derivs_ng(rr,y,dydx)\n      use module_model\n', 1)
open('module_int.f90', 'w').write(s)
s = open('module_solver.f90').read()
old = '''         if(rs /= r(nsl)) then
            call yint(nsl-1,is,r(nsl),rs,y2,rad_derivs_ng,scale2)
         end if'''
new = '''         if(rs /= r(nknot)) then
            call yint(nknot-1,is,r(nknot),rs,y2,rad_derivs_ng,scale2)
         end if'''
assert s.count(old) == 1
s = s.replace(old, new)
open('module_solver.f90', 'w').write(s)
PY
gfortran -O2 -fopenmp -I/usr/include -o yspec module_kinds.f90 module_util.f90 \
  module_spline.f90 module_model.f90 module_function.f90 module_fourier.f90 \
  module_start.f90 module_int.f90 module_solver.f90 yspec.f90 -lfftw3

for dir in "$work"/tests/*/; do
  name=$(basename "$dir")
  [[ -f "$dir/input.txt" ]] || continue
  case="$root/tests/regression/$name"
  mkdir -p "$case/reference"
  rm -f "$case"/reference/*
  model=$(python3 -I - "$dir/input.txt" <<'PY'
import sys
lines = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.strip().startswith('#')]
print(lines[1].strip('"\''))
PY
)
  run="$work/run-$name"
  mkdir -p "$run"
  cp "$dir/input.txt" "$work/examples/$model" "$run/"
  (cd "$run" && ./../src/yspec input.txt > /dev/null)
  cp "$run"/out.* "$case/reference/"
  echo "$name: model $model, $(ls "$case/reference" | wc -l) receivers"
done
