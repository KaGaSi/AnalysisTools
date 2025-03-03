Examples of using coordinate manipulation options -cx/y/z, -m, and -sc. The
snapshots in the figure within the manual were obtained from data in
../AngleMolecules. To get all the resulting coordinate files, just run the
run.sh script that contains the command below.

To prepare a coordinate file traj.vtf in this directory, run:

Selected ../AngleMolecules/traj.vcf traj.vtf -i ../AngleMolecules/traj.data

To see the results of the following examples, either use visualization software
or follow the arrows in the figure within the manual

1) Using all three options in one command (shows default behaviour):

  Selected traj.vtf full.vtf -cx 0 0.2 0.6 0.8 -cy 0 0.2 0.5 0.9 -sc 2 -m 1 1 0

  This is equivalent to running the following three commands in sequence:

  a) constrain coordinates: no -cz option means, z-coordinate isn't constrained;
     -cx option means x-coordinates of saved beads must be 0-20% or 60-80% of
     the simulation box size in x-direction; simultaneously, the -cy means
     y-coordinates are similarly constrained
    Selected traj.vtf cx_cy.vtf -cx 0 0.2 0.6 0.8 -cy 0 0.2 0.5 0.9
  b) move all beads by the vector (<x box size>, <y box size>, 0)
    Selected cx_cy.vtf cx_cy-m.vtf -m 1 1 0
  c) scale all coordinates by dividing them by 2
    Selected cx_cy-m.vtf cx_cy-m-sc.vtf -sc 2

2) First move beads and then scale coordinates - the first command uses
   cx_cy.vtf from 1):

  Selected cx_cy.vtf cx_cy-sc.vtf -sc 2
  Selected cx_cy-sc.vtf cx_cy-sc-m.vtf -m 1 1 0

3) First scale coordinates and then constrain the saved ones - the first

  Selected traj.vtf sc2.vtf -sc 2
  Selected sc2.vtf sc2-cx_cy.vtf -cx 0 0.2 0.6 0.8 -cy 0 0.2 0.5 0.9
