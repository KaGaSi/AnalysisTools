In this example, a wire-like aggregate stretching the z-axis length of the
simulation box is created. The bash script Wire.sh creates the structure by
adding an A5B5 molecule (defined in A5B5.FIELD file) one by one; in the end, the
script also adds solvent particles outside the structure to obtain a fully
usable lammps data file (overall bead number density is 3 in keeping with
classical dissipative particle dynamics). See the Wire.sh script for details.

A simple input lammps script (lmp.in) is also included, so a lammps simulation
can be run out of the box should anyone be interested. The picture Wire.jpg
shows the resulting simulation box.
