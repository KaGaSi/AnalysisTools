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
#   D1E9.FIELD - 10-bead molecules to be added
################################################################################
bin="AddToSystem"
in_field1="${ROOT}/init.FIELD"
in_field2="${ROOT}/C.FIELD"
in_field3="${ROOT}/D1E9.FIELD"

# 1) Build initial system from scratch
${bin} - ${in_field1} 1.vtf -b 10 10 10

# 2) exchange the most numerous bead type (S1)
${bin} 1.vtf ${in_field2} 2.vtf

# 3) exchange specified bead type
${bin} 1.vtf ${in_field2} 3.vtf -xb S2

# 4) place beads specific distance from any bonded beads
${bin} 1.vtf ${in_field2} 4.vtf --add --bonded -ld 1.5 -hd 3

# 5) use molecule's first/last bead for -cx/-cy/-cz constraint instead of com
${bin} 1.vtf ${in_field2} 5.vtf --add --tail -cx 0.5 1
${bin} 1.vtf ${in_field2} 6.vtf --add --head -cx 0.5 1
