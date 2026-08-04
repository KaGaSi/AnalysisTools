#!/bin/bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/../../func.sh"

###############################################################################
# This scripts creates a wire-like aggregate spanning z-axis of the simulation
# box. The aggregates is composed of 20 layers of molecules, each layer
# comprises of 6 A4B6 molecules (defined in the A4B6.FIELD file) arranged in a
# 6-point star (i.e., 60 degrees apart). The simulation box size is defined at
# the first line of the A4B6.FIELD file (it's chosen so the layers are 0.8 units
# apart which is fairly reasonable for dissipative particle dynamics)
#
# Should a different structure be wanted, the number of layers can be changed in
# the outer for loop (it is 20 by default), while the number of molecules per
# layer can be changed in the inner for loop (it is 6 by default).
#
# The distance between the layers is specified by the cz1 and cz2 variables as
# lower and upper bound, respectively, for z-axis; the range is small so that
# all molecules of the layer have basically the same coordinate. By default, the
# distance between layers is 0.05 in relative units, meaning 0.8 in box units
# (the box's z length is 16).
#
# The angles between molecules are governed by the ${angle} variable; by
# default, it is 60 degrees between molecules in the one layer, while the layers
# are rotated relative to the privous one by 25 degrees, creating a screw
# structure.
###############################################################################

# path to the executable - default assumes AddToSystem is in the $PATH
bin="AddToSystem"
in_field1="${ROOT}/A4B6.FIELD"
in_field2="${ROOT}/W.FIELD"

# make 20 layers, each a 6-pointed star, stacking them in z-axis direction
for (( i=0; i<20; i++ )); do # go over whole z- coordinate
  # constraint in z-axis; i.e., what layer is being build
  cz1=$(calc "${i}*0.05" 2)
  cz2=$(calc "(${i}+1)*0.05" 2)
  for (( j=0; j<6; j++ )); do # create one layer
    # add z- to x- and y-axis constraints
    constraint="-cx 0.49 0.51 -cy 0.49 0.51 -cz ${cz1} ${cz2} --head"
    if [[ ${i} == 0 && ${j} == 0 ]]; then # first molecule creates a new file
      ${bin} - ${in_field1} new.data ${constraint} --no-rotate -s $((i*j+1))
    else # other molecules are added to an existing file
      # in each layer, molecules are 60° apart; layers are shifted by 25°
      angle=$(( j * 60 + i * 25 ))
      ${bin} old.data ${in_field1} new.data ${constraint} -a 0 0 ${angle} --add -s $((i*j+2))
    fi
    # move new file to be used in the next cycle as an input file
    mv {new,old}.data
  done
done
# add solvent outside the aggregate; corresponds to overall number density 3
${bin} old.data ${in_field2} Wire.data -ld 0.5 -bt A B --add -o Wire.vtf -s 111
# remove temporary file
rm old.data
