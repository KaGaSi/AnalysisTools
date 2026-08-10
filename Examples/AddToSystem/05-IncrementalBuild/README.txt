This example uses a simple sequential addition of species to create a box filled
two separated solvent phases and one molecule type spread throughout both
phases. A snapshot of the final system is shown in system.jpg.

1) Create a system from scratch by adding molecules to the whole box

2) Add first solvent to two-thirds of the box (increase the number of beads)
  relevant options: --add -cx 0 0.66

3) Add second solvent to one-third of the box (increase the number of beads)
  relevant options: --add -cx 0.67 1

Run all steps:
  ./run.sh

The new systems' composition should be investigated using the Info utility, and
the positions of added species through vmd <file>.vtf -e vmd.tcl command.
