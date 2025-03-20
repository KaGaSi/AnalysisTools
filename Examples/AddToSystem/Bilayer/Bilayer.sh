#!/bin/bash

###############################################################################
# This script creates a bilayer spanning yz plane of the simulation box and
# adds salt into the water phase at one side of the bilayer. One corner of the
# bilayer is composed of different molecules. The total number of particles
# correspond to particle density of three, i.e., typical density for a
# dissipative particle dynamics simulation.
#
# The first layer of the bilayer is done in one go, while the second layer is
# created in three steps. The first and second step create the larger part of
# the layer from two rectangles specified by -cy and -cz options. The third step
# fills in the 'notch'.
#
# Then, water with ions is added to one side of the bilayer, while pure water is
# added to the second side.
###############################################################################

# relative path to the AddToSystem utility
Add="../../../build/bin/AddToSystem"

# create first layer (400 A5B1 molecules)
sed "s/NUMBER/200/" A5B1.FIELD > FIELD
${Add} - FIELD 1.vtf -cx 8.9 9.0 --no-rotate --head --real

# second layer - create the 'notched' part composed of 390 A5B1 molecules
sed "s/NUMBER/160/" A5B1.FIELD > FIELD
${Add} 1.vtf FIELD 2.vtf -cx 8 8.1 -cy 2 10 -a 180 0 0 --add --head --real
sed "s/NUMBER/30/" A5B1.FIELD > FIELD
${Add} 2.vtf FIELD 1.vtf -cx 8 8.1 -cy 0 2 -cz 2 10 -a 180 0 0 --add --head --real
# second layer - fill the 'notch' with 10 E5D1 molecules
${Add} 1.vtf E5D1.FIELD 2.vtf -cx 8 8.1 -cy 0 2 -cz 0 2 -a 180 0 0 --head --add --real
# add only water to one side
sed "s/WATER/2100/" W.FIELD | sed "s/ION/0/" > FIELD
${Add} 2.vtf FIELD 1.vtf -cx 12 20 --add --real
# add water and ions to the other side
sed "s/WATER/1300/" W.FIELD | sed "s/ION/100/" > FIELD
${Add} 1.vtf FIELD Bilayer.vtf -cx 0 5 --add --real

# remove extraneous files
rm {1,2}.vtf FIELD
