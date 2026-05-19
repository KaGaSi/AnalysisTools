#!/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Demonstrates exchange mode (-xb), bonded-bead distance constraints
# (--bonded), placement anchored to the tail bead (--tail), and reproducible
# placement (-s).
#
# Files:
#   init.FIELD - initial system with two-component solvent and one molecule type
#   C.FIELD    - C particles to be exchanged
#   D10.FIELD  - 10-bead molecules to be added
################################################################################
bin="AddToSystem"
in_field1="${ROOT}/init.FIELD"
in_field2="${ROOT}/C.FIELD"
in_field3="${ROOT}/D10.FIELD"

# 1) Build initial system from scratch
${bin} - ${in_field1} 1.vtf -b 10 10 10

# 2) exchange the most numerous bead type (S1)
${bin} 1.vtf ${in_field2} 2.vtf

# 3) exchange specified bead type
${bin} 1.vtf ${in_field2} 3.vtf -xb S2

# 4) --bonded
# Use all bonded beads in the system as the reference set for -ld/-hd, instead
# of naming specific bead types with -bt.  Here 100 E beads are placed at
# least 1.5 units away from any A or B bead in an AB molecule.
# 4) place
${bin} 1.vtf ${in_field2} 4.vtf --add --bonded -ld 1.5

# 5) --tail
# By default, -cx/-cy/-cz constrain the molecule's geometric centre.  --head
# uses the first bead; --tail uses the last bead.  Here AB molecules are placed
# with their tail (B bead) in the right half of the box (x in [5, 10]).
${bin} 1.vtf ${in_field2} 5.vtf --add --tail -cx 0.5 1
