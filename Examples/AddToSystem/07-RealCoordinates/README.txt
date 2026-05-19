This example shows --real with negative box coordinates and the -ebt option for
LAMMPS data file output.

The starting system is generated inline by the script: 200 W beads in a
20x20x20 box centered at the origin, with coordinates from -10 to 10 in each
dimension.  This is typical of LAMMPS simulations that use a symmetric box.

--real with negative coordinates:
  Without --real, -cx/-cy/-cz values are fractional (0 to 1, relative to the
  output box length).  --real switches them to the same coordinate space as
  the input file, including negative values:

    AddToSystem system.lammpstrj W.FIELD 1.vtf --add --real -cx -10 0

  Places 100 new W beads in the left half of the box (x in [-10, 0]).

    AddToSystem system.lammpstrj W.FIELD 2.vtf --add --real -cz -1 1

  Places 100 W beads in the central z-slab (z in [-1, 1]).

  Internally, bead coordinates are stored relative to the box origin (0 to
  box length).  --real values are automatically shifted by box.Low on input,
  so you do not need to account for the offset manually.

Extra bead types for LAMMPS output (-ebt):
  When the output is a LAMMPS data file (.data), -ebt <n> adds n placeholder
  bead type entries to the Masses section.  This is useful when the type count
  must match a fixed template, or when additional types will be defined later
  via LAMMPS pair_coeff commands:

    AddToSystem system.lammpstrj W.FIELD out.data --add --real -cx -10 0 -ebt 2

  The output file will have 3 entries in the Masses section (1 real + 2
  placeholders named "extra"), even though only W beads are present.
