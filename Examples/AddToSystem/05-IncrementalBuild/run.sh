#!/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Simple example of sequential additions to create a new system.
################################################################################

# path to the executable - default assumes AddToSystem is in the $PATH
bin=AddToSystem

# 1) randomly place CD 2-particle molecules into the whole box
${bin} - ${ROOT}/CD.FIELD 1.vtf
# 2) add A fluid into two-thirds of the box
${bin} 1.vtf ${ROOT}/A.FIELD 2.vtf --add -cx 0 0.66
# 3) add W fluid into the remaining one-third of the box
${bin} 2.vtf ${ROOT}/B.FIELD 3.vtf --add -cx 0.67 1
