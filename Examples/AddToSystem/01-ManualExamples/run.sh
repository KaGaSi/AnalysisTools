#!/bin/env bash

################################################################################
# Commands to generate examples from the figure 3.1 in the manual
#
# Files:
#   in.lammpstrj ... input coordinate file (LAMMPS custom format)
#   FIELD ... DL_MESO-inspired file containing beeds to add
################################################################################


# Figure 3.1a
AddToSystem in.lammpstrj FIELD fig3_1a.vtf --add -ld 3 -bt A
# Figure 3.1b
AddToSystem in.lammpstrj FIELD fig3_1b.vtf --add -hd 4 -bt A
# Figure 3.1c
AddToSystem in.lammpstrj FIELD fig3_1c.vtf --add -ld 3 -hd 4 -bt A
# Figure 3.1d
AddToSystem in.lammpstrj FIELD fig3_1d.vtf --add -ld 3 -hd 4 -bt A -cx 0.5 1
# Figure 3.1e
AddToSystem in.lammpstrj FIELD fig3_1e.vtf --add -ld 3 -hd 4 -bt A -cx 0.5 1 -b 30 20 25
# Figure 3.1f
AddToSystem in.lammpstrj FIELD fig3_1f.vtf --add -ld 3 -hd 4 -bt A -cx 0.5 1 -b 30 20 25 -off -0.2 0.2 0
