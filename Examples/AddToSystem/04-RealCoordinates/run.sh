#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Demonstrates the two coordinate spaces of the -cx/-cy/-cz constraints
#
# By default, the constraints are fractions of the box; --real switches them to
# the coordinate space of the input file. For an orthogonal box the two are
# interchangeable. A tilted box has no Cartesian sub-box, so only the fractional
# form is available there.
#
# Files:
#   ortho.lammpstrj ... system in orthogonal box with negative box coordinates
#   tilted.lammpstrj ... system in tilted box with negative box coordinates
#   W.FIELD ... W beads to add
################################################################################
bin="AddToSystem"
in_field="${ROOT}/W.FIELD"
ltrj_ortho="${ROOT}/ortho.lammpstrj"
ltrj_tilt="${ROOT}/tilted.lammpstrj"

# 1) place beads in a slab in orthogonal box: -2 to 2 in real units is 0.4 to
#    0.6 in fractions, so the two commands are equivalent and generate identical
#    coordinates
${bin} ${ltrj_ortho} ${in_field} 01.vtf -s 1 --add --real -cx -2 2
${bin} ${ltrj_ortho} ${in_field} 02.vtf -s 1 --add -cx 0.4 0.6

# 2) in the tilted box, the fractional constraint carves out a smaller cell of
#    the same shape instead of a Cartesian slab
${bin} ${ltrj_tilt} ${in_field} 03.vtf -s 2 --add -cx 0.4 0.6
#    ... using --real -cx -2 2 would error-out, as no Cartesian sub-box of a
#    tilted cell is periodic
