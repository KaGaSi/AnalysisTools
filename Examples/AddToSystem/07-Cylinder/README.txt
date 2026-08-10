This example constructs a cylindrical aggregate stretching the length of the
simulation box: 20 layers of 6 A4B6 diblock copolymers (defined in A4B6.FIELD),
each layer a 6-point star lying in the xy plane, stacked along the z-axis. The
rest of the box is filled with solvent (overall number density of 3 was used in
keeping with classical dissipative particle dynamics). A snapshot of the final
system is shown in system.jpg.

The 10x10x16 box comes from the first line of A4B6.FIELD, its z length chosen so
that the 20 layers end up 0.8 units apart. The molecules are added one at a time
by a pair of nested loops, each call reading the file written by the previous
one, so the aggregate is built molecule by molecule. Unlike in example 06, the
constraints are left in the default fractional units (0 to 1) and all molecules
are anchored by their first bead (--head), which is the free end of the A4
block, placing the A beads on the axis and the B ones outside.

1) the aggregate: 120 A4B6 molecules, one per call, each confined to a narrow
   column on the box axis and to the z-slab of its own layer; the first one
   creates the system from scratch, the rest are rotated around the z-axis, by
   60 degrees between the molecules of one layer and by another 25 degrees
   between successive layers, which twists the stack into a screw

   relevant options: -cx 0.49 0.51 -cy 0.49 0.51 ... the column along the axis
                     -cz <lo> <lo+0.05>           ... the layer, 0.8 units thick
                     --head --no-rotate           ... the very first molecule
                     -a 0 0 <angle>               ... all the others

2) solvent: 3600 W beads outside the aggregate, keeping a minimum distance from
   it instead of a constraint on the coordinates

   relevant options: -ld 0.5 -bt A B ... at least 0.5 from any A or B bead
                     -o system.vtf   ... a vtf alongside the data file output

The result is an A4 core reaching 1.5 from the axis and a B6 corona out to 3.

Run all steps:
  ./run.sh

The new system's composition should be investigated using the Info utility, and
the positions of the added species through vmd system.vtf -e vmd.tcl command.

The created lammps data file can be used with the supplied lammps input script
(lmp.in) in a fully workable dissipative particle dynamics simulation.
