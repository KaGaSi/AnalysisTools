#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Demonstrates exchange mode (-xb), bonded-bead distance constraints
# (--bonded), placement anchored to the tail bead (--tail), and reproducible
# placement (-s).
#
# Note that while added unbonded beads are always wrapped inside the box, added
# molecules may stick out.
#
# Files:
#   init.FIELD ... initial system (two-component solvent and one molecule type)
#   C.FIELD ... C particles to be exchanged
#   D1E9.FIELD ... 10-bead molecules to be added
################################################################################
bin="AddToSystem"
in_field1="${ROOT}/init.FIELD"
in_field2="${ROOT}/C.FIELD"
in_field3="${ROOT}/D1E9.FIELD"
vtf="01.vtf"

# 1) Build initial system from scratch
${bin} - ${in_field1} ${vtf} -s 1 -b 20 20 20

# 2) exchange either the most numerous bead type (S1) or the -xb-specified one
${bin} ${vtf} ${in_field2} 02.vtf -s 2
${bin} ${vtf} ${in_field2} 03.vtf -s 3 -xb S2

# 3) place beads specific distance from any bonded beads
${bin} ${vtf} ${in_field2} 04.vtf -s 4 --add --bonded -ld 0.1 -hd 0.2

# 4) use molecule's com/first/last bead for -cx constraint (don't rotate the
#    added molecule, so the effect is clearly visible)
${bin} ${vtf} ${in_field3} 05.vtf -s 5 --add -cx 0 0.01 --no-rotate
${bin} ${vtf} ${in_field3} 06.vtf -s 6 --add -cx 0 0.01 --no-rotate --head
${bin} ${vtf} ${in_field3} 07.vtf -s 7 --add -cx 0 0.01 --no-rotate --tail

# 5) use molecule's com/first/last bead for -hd constraint
${bin} ${vtf} ${in_field3} 08.vtf -s  8 --add --bonded -hd 0.1
${bin} ${vtf} ${in_field3} 09.vtf -s  9 --add --bonded -hd 0.1 --head
${bin} ${vtf} ${in_field3} 10.vtf -s 10 --add --bonded -hd 0.1 --tail
