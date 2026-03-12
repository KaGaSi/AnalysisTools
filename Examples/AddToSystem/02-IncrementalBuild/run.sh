#!/bin/env bash

################################################################################
# Simple example of sequential additions to create a new system.
################################################################################

# path to the executable - default assumes AddToSystem is in the $PATH
bin=AddToSystem

# 1) randomly place CD 2-particle molecules into the whole box
${bin} - CD.FIELD 1.vtf
# 2) add A fluid into two-thirds of the box
${bin} 1.vtf A.FIELD 2.vtf --add -cx 0 0.66
# 3) add W fluid into the remaining one-third of the box
${bin} 2.vtf B.FIELD 3.vtf --add -cx 0.67 1
