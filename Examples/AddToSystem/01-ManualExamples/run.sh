#!/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Commands to generate examples from the figure 3.1 in the manual
#
# Files:
#   in.lammpstrj ... input coordinate file (LAMMPS custom format)
#   FIELD ... DL_MESO-inspired file containing beeds to add
################################################################################
in_ltrj="${ROOT}/in.lammpstrj"
in_field="${ROOT}/FIELD"

# Figure 3.1a
AddToSystem ${in_ltrj} ${in_field} fig3_1a.vtf -s 1 --add -ld 3 -bt A
# Figure 3.1b
AddToSystem ${in_ltrj} ${in_field} fig3_1b.vtf -s 2 --add -hd 4 -bt A
# Figure 3.1c
AddToSystem ${in_ltrj} ${in_field} fig3_1c.vtf -s 3 --add -ld 3 -hd 4 -bt A
# Figure 3.1d
AddToSystem ${in_ltrj} ${in_field} fig3_1d.vtf -s 4 --add -ld 3 -hd 4 -bt A -cx 0.5 1
# Figure 3.1e
AddToSystem ${in_ltrj} ${in_field} fig3_1e.vtf -s 5 --add -ld 3 -hd 4 -bt A -cx 0.5 1 -b 30 20 25
# Figure 3.1f
AddToSystem ${in_ltrj} ${in_field} fig3_1f.vtf -s 5 --add -ld 3 -hd 4 -bt A -cx 0.5 1 -b 30 20 25 -off -0.2 0.2 0
