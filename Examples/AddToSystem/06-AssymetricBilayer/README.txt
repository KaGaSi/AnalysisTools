In this example, a bilayer is created composed primarily of one molecule type
(defined in A5B1.FIELD) with a 'notch' in one of the layers composed of a
different molecule type (defined in E5D1.FIELD). The rest of the simulation box
is filled with pure water on one side and water with salt on the other (overall
number density of 3 was used in keeping with classical dissipative particle
dynamics). A snapshot of the final system is shown in system.jpg.

The 20x10x10 box comes from the first line of the FIELD files and the bilayer is
built perpendicular to the x-axis, i.e., in the yz plane. Every step adds to the
result of the previous one, writing <count>.vtf, with the molecule and bead
counts substituted into <count>.FIELD by sed. All constraints are in real units
(--real) and anchored to the molecule's first bead (--head), which is the far
end of the A5/E5 tail; every step is seeded (-s) to keep the script
reproducible.

1) first layer: 200 A5B1 molecules, unrotated, so that they all point towards
   higher x, with their tail ends in a thin slab

   relevant options: -cx 8.9 9.0 --no-rotate --head

2) second layer: 190 A5B1 molecules turned around, so that the tails of the two
   layers meet in the middle of the box and the B beads end up facing the water.
   It is added as two rectangles that leave a notch at one corner:

   relevant options: -a 0 0 180                 ... turn the molecules around
                     -cx 8 8.1 -cy 2 10         ... the bigger rectangle
                     -cx 8 8.1 -cy 0 2 -cz 2 10 ... the smaller one

3) the notch: 10 E5D1 molecules into the corner left free by the two rectangles

   relevant options: -cx 8 8.1 -cy 0 2 -cz 0 2

4) water on one side of the bilayer: 2100 W beads

   relevant options: -cx 12 20

5) water and salt on the other side: 1300 W beads plus 100 P and 100 N ions

   relevant options: -cx 0 5

Run all steps:
  ./run.sh

The new systems' composition should be investigated using the Info utility, and
the positions of added species through vmd <file>.vtf -e vmd.tcl command.
