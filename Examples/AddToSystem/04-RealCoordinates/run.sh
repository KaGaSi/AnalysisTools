#!/bin/env bash

################################################################################
# Demonstrates --real with negative box coordinates (LAMMPS-style centered box)
# and -ebt to reserve extra bead type slots in a LAMMPS data file output.
#
# LAMMPS simulations often use boxes centered at the origin, e.g., -10 to 10
# in each dimension.  Using --real lets you specify -cx/-cy/-cz constraints in
# the same coordinate space as the input file, including negative values.
#
# Files:
#   W.FIELD — 100 W solvent beads to add
################################################################################

bin=AddToSystem

# Generate a starting system: 200 W beads in a 20x20x20 box centered at the
# origin (coordinates from -10 to 10 in each dimension).
awk 'BEGIN {
  srand(1)
  print "ITEM: TIMESTEP"
  print "0"
  print "ITEM: NUMBER OF ATOMS"
  print "200"
  print "ITEM: BOX BOUNDS pp pp pp"
  print "-10.0 10.0"
  print "-10.0 10.0"
  print "-10.0 10.0"
  print "ITEM: ATOMS id element x y z"
  for (i = 1; i <= 200; i++) {
    printf "%d W %.4f %.4f %.4f\n", i, rand()*20-10, rand()*20-10, rand()*20-10
  }
}' > system.lammpstrj

# Add 100 W beads to the left half of the box (x from -10 to 0).
# Without --real, fractional coordinates 0-1 would be needed; with --real the
# values are interpreted directly in the file's coordinate space.
${bin} system.lammpstrj W.FIELD 1.vtf --add --real -cx -10 0

# Add 100 W beads to a central slab (-1 < z < 1).
${bin} system.lammpstrj W.FIELD 2.vtf --add --real -cz -1 1

# -ebt: reserve extra bead type entries in a LAMMPS data file.  Here the system
# has 1 type (W); -ebt 2 adds 2 placeholder entries, giving 3 types total in
# the Masses section.  Useful when additional types will be assigned later via
# LAMMPS pair_coeff or when running a series of builds with a fixed type count.
${bin} system.lammpstrj W.FIELD out.data --add --real -cx -10 0 -ebt 2
