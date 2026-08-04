#!/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

###############################################################################
# This script creates a bilayer spanning yz plane of the simulation box and adds
# salt into the water phase at one side of the bilayer. One corner of the
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
#
# Every step is given a seed (-s) so that the script is reproducible; remove
# the -s options to get a different bilayer on every run.
#
# While the packaged Bilayer.jpg shows a snapshot of the constructed bilayer,
# the newly generaged Bilayer.vtf structure/coordinate file may be viewed
# through, e.g., vmd.
###############################################################################

# path to the executable - default assumes AddToSystem is in the $PATH
bin="AddToSystem"

surf1="${ROOT}/A5B1.FIELD" # FIELD file defining the 1st surfactant molecule
surf2="${ROOT}/E5D1.FIELD" # FIELD file defining the 2nd surfactant molecule
solvent="${ROOT}/W.FIELD" # FIELD file defining the water and salt beads
count=0 # name the files in every step as <count>.vtf and <count>.FIELD

# create first layer from scratch (200 A5B1 molecules)
count=$((count+1))
sed "s/NUMBER/200/" ${surf1} > ${count}.FIELD
${bin} - ${count}.FIELD ${count}.vtf -cx 8.9 9.0 --no-rotate --head --real -s ${count}

# second layer - create the 'notched' part composed of 190 A5B1 molecules
# a) bigger rectangle
count=$((count+1))
sed "s/NUMBER/160/" ${surf1} > ${count}.FIELD
${bin} $((count-1)).vtf ${count}.FIELD ${count}.vtf -cx 8 8.1 -cy 2 10 -a 0 0 180 --add --head --real -s ${count}
# b) smaller rectangle
count=$((count+1))
sed "s/NUMBER/30/" ${surf1} > ${count}.FIELD
${bin} $((count-1)).vtf ${count}.FIELD ${count}.vtf -cx 8 8.1 -cy 0 2 -cz 2 10 -a 0 0 180 --add --head --real -s ${count}
# second layer - fill the 'notch' with 10 E5D1 molecules
count=$((count+1))
${bin} $((count-1)).vtf ${surf2} ${count}.vtf -cx 8 8.1 -cy 0 2 -cz 0 2 -a 0 0 180 --head --add --real -s ${count}
# add only water to one side
count=$((count+1))
sed "s/WATER/2100/" ${solvent} | sed "s/ION/0/" > ${count}.FIELD
${bin} $((count-1)).vtf ${count}.FIELD ${count}.vtf -cx 12 20 --add --real -s ${count}
# add water and ions to the other side
count=$((count+1))
sed "s/WATER/1300/" ${solvent} | sed "s/ION/100/" > ${count}.FIELD
${bin} $((count-1)).vtf ${count}.FIELD ${count}.vtf -cx 0 5 --add --real -s ${count}

cp ${count}.vtf Bilayer.vtf
