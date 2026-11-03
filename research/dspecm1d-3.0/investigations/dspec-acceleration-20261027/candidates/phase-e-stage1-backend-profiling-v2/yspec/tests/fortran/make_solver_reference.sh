#!/usr/bin/env bash
# Regenerates tests/data/solver/*.txt: the responses of the Fortran YSpec's
# radial solvers to unit jumps, for a list of degrees and frequencies, against
# which tests/TestSolver.cpp checks the C++ solvers.
#
#   tests/fortran/make_solver_reference.sh /path/to/YSpec [commit]
#
# The Fortran is exported at the given commit (default a3ae638, the cleaned
# code), its step density raised sixteen-fold (nw = 4 -> 64) so that it is
# accurate to about 1e-8. Two of its faults affect l = 0 without gravity:
# rad_derivs_ng keeps the gravity term of rad_amat, which is patched out here;
# and rad_solver_ng starts its downward integration at the sea floor, not the
# surface, when the receiver is below the source, which is not patched: those
# l = 0 lines are left out instead.
set -euo pipefail
yspec=$(realpath "$1")
commit=${2:-a3ae638}
here=$(cd "$(dirname "$0")" && pwd)
out="$here/../data/solver"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

git -C "$yspec" archive "$commit" src | tar -x -C "$work"
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
PY
gfortran -O2 -o driver module_kinds.f90 module_util.f90 module_spline.f90 \
  module_model.f90 module_function.f90 module_start.f90 module_int.f90 \
  module_solver.f90 "$here/driver.f90"

degrees="12
0 0.5
0 6.0
1 2.0
1 9.0
2 0.4
2 5.0
10 3.0
30 10.0
80 12.0
150 25.0
300 30.0
600 40.0"

# name model attenuation gravity source_depth_km receiver_depth_km
while read -r name model ats grav sd rd; do
  file="$out/$name.txt"
  {
    echo "# model $model"
    echo "# attenuation $ats"
    echo "# gravity $grav"
    echo "# source_depth_km $sd"
    echo "# receiver_depth_km $rd"
    echo "# time_series_min 60"
    printf '%s\n%s %s 0\n%s %s 60\n%s\n' "$yspec/examples/$model" "$ats" \
      "$grav" "$sd" "$rd" "$degrees" | ./driver
  } > "$file"
  # The receiver-below-source, no-gravity, l = 0 response is wrong in the
  # Fortran when there is an ocean.
  if [[ $grav == 0 && $(echo "$rd > $sd" | bc) == 1 && $model != prem.200.no ]]; then
    sed -i '/^R /d' "$file"
  fi
done <<'CASES'
full prem.200 1 2 29.5 3
full_below prem.200 1 2 10 50
cowling prem.200 1 1 150 3
cowling_below prem.200 1 1 10 50
none prem.200 1 0 29.5 3
none_below prem.200 0 0 10 50
fic_full prem.200.fic 1 2 29.5 3
fic_cowling prem.200.fic 1 1 600 3
fic_none prem.200.fic 0 0 10 50
no_ocean_full prem.200.no 1 2 29.5 0
no_ocean_none prem.200.no 1 0 29.5 0
CASES
echo "wrote $out"
