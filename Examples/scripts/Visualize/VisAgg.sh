#!/bin/bash

################################################################################
# A bash script to visualize aggregates via vmd.
#
# Specified agg and vtf files as well as classic linux utilities are used to
# generate a tcl script, using it with vmd to render a tga picture. The tcl
# script is then deleted.
#
# Note that it was designed for the provided VisAgg* files, so to use it with
# different files, some parts of the script may need changes (e.g., the position
# of the box via translate/scale vmd commands or the modselect lines in the vmd
# representation).
################################################################################

################################################################################
# file names and other variables
################################################################################
name=VisAgg
tcl=${name}.tcl # tcl script for vmd
agg=${name}.agg # aggregate file (generated via Aggregates)
vtf=${name}.vtf # structure and coordinate file
pic=${name}.tga
res=300 # resolution (size of vmd canvas in pixels)

##############################################################################
# start of the tcl script for vmd
##############################################################################
{
  echo "package require pbctools"
  echo "color Display Background white"
  echo "pbc box"
  echo "set mol 0"
  echo "set rep 0"
  echo "mol modselect \${rep} \${mol} none"
  echo "display resetview"
  echo "translate by -0.574 -0.565 0"
  echo "axes location Off"
} > ${tcl}
##############################################################################
# find number of aggregates for each timestep in the agg file and pick the
# timestep to visualize
##############################################################################
lines1=$(grep -n "Step:" "${agg}" | awk -F ":" '{print $1}')
lines1=(${lines1}) # array with those line numbers
i=1 # choose the timestep; for loop can be used to go over all steps in vtf
i1=$(( i - 1 )) # number line for the following '... Step: ...'
l1=$(( ${lines1[${i1}]} + 2 )) # first line with an aggregate
l2=$(( ${lines1[${i}]} - 0 )) # last line with an aggregate
##############################################################################
# loop over all aggregate lines in step ${i}, creating vmd representation (or
# more) for each aggregate
##############################################################################
n=-1 # ColorID
for (( j=l1; j<l2; j++ )); do
  line=$(head -n ${j} ${agg} | tail -n 1) # read the appropriate line
  line="resid ${line##*:}" # remove <size> : from the line
  ############################################################################
  # print all aggregates
  ############################################################################
  n=$(( n + 1)) # ColorID
  if (( n > 32 )); then # there are 32 colours defined in vmd
    n=0
  fi
  {
    # whol aggregate via small balls
    echo "set rep [expr \$rep + 1]"
    echo "mol addrep \${mol}"
    echo "mol modselect \${rep} \${mol} ${line}"
    echo "mol modstyle  \${rep} \${mol} cpk 0.5 0.3"
    echo "mol modcolor  \${rep} \${mol} ColorID ${n}"
    # accentuate hydrophobic core via big balls
    echo "set rep [expr \$rep + 1]"
    echo "mol addrep \${mol}"
    echo "mol modselect \${rep} \${mol} ${line} and name B"
    echo "mol modstyle  \${rep} \${mol} cpk 1.4 0.3"
    echo "mol modcolor  \${rep} \${mol} ColorID ${n}"
  } >> ${tcl}
done
{
  echo "render TachyonInternal ${pic}"
  echo "exit"
} >> ${tcl}
##############################################################################
# run vmd with the new script and delete the script
##############################################################################
vmd ${vtf} -e ${tcl} -size ${res} ${res} -dispdev text
rm ${tcl}
