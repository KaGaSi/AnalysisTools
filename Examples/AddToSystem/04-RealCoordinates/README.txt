This example demonstrates the two coordinate spaces of the -cx/-cy/-cz
constraints and the one case where they are not interchangeable.

Starting systems:
  Two ready-made lammps trajectories, each holding 200 A beads in a box with
  negative box coordinates to show AddToSystem deals correctly with those.

  ortho.lammpstrj  ... orthogonal box
  tilted.lammpstrj ... the same box but tilted (gamma = 60 degrees)

1) Orthogonal box (result shown in ortho.jpg):
  --real option and the default (fractional units 0-1) are equivalent:

  relevant options: --real -cx -2 2 ... 01.vtf, central slab, four units thick
                    -cx 0.4 0.6     ... 02.vtf, the same slab, as a fraction

2) Tilted box (result shown in tilted.jpg):
  as there is no periodic Cartesian slab in a tilted/triclinic box, --real
  cannot be used. The provided fractions cut along the cell vectors, carving out
  a smaller cell of the same shape:

  relevant options: -cx 0.4 0.6     ... 03.vtf, the same fractions, slanted
                    --real -cx -2 2 ... would error out

Run all steps:
  ./run.sh

The new systems' composition should be investigated using the Info utility, and
the positions of added species through vmd <file>.vtf -e vmd.tcl command.
