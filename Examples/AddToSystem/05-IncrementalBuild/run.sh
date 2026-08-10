#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Simple example of sequential additions to create a new system (reproducible
# through using -s option).
#
# Files:
#   CD.FIELD ... two-particle molecule to add
#   A.FIELD and B.FIELD ... unbonded beads to add
################################################################################

# path to the executable - default assumes AddToSystem is in the $PATH
bin=AddToSystem

# 1) randomly place CD 2-particle molecules into the whole box
${bin} - ${ROOT}/CD.FIELD 1.vtf -s 1
# 2) add A beads into two-thirds of the box
${bin} 1.vtf ${ROOT}/A.FIELD 2.vtf -s 2 --add -cx 0 0.66
# 3) add B beads into the remaining one-third of the box
${bin} 2.vtf ${ROOT}/B.FIELD 3.vtf -s 3 --add -cx 0.67 1
