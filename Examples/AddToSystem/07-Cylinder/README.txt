This example constructs a wire-like aggregate stretching the length of the
simulation box. The bash script Wire.sh adds molecules one by one to create the
aggregates and then fills the rest of the box with solvent particles.

The picture Wire.jpg shows the resulting simulation box; it can be generated
anew via the vmd program using the supplied vmd.tcl script.

The created lammps data file can be used with the supplied lammps input script
(lmp.in) in a fully workable dissipative particle dynamics simulation using
lammps.
