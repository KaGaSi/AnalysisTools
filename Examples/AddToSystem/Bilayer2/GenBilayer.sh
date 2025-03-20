#!/bin/bash

################################################################################
# Created on 18/06/2024 for real simulation case                               #
#                                                                              #
# Uses: AnalysisTools commit 78872cecba3895259f4076023d1ca8386ed443b0          #
#       round() and run_err() bash functions from Scripting/shell/func.sh      #
#                                                                              #
# Script to create a bilayer from three kinds of molecules defined in three    #
# FIELD files (and supplemented by a FIELD file containing solvent and         #
# counterion beads). The bilayer is created in the middle of the z-direction   #
# of the box.                                                                  #
#                                                                              #
# The numbers of all beads/molecules as well as the box size are read from the #
# FIELD files. The script uses AddToSystem utility to first create a bilayer   #
# from the first type of molecule using the first FIELD file, then adds the    #
# molecules from the second FIELD, then from the last one and finally adds     #
# counterions and solvent beads, hopefully so the beads aren't inside the      #
# bilayer. There's also a 'reasonable' gap between the two leaflets of the     #
# bilayer.                                                                     #
#                                                                              #
# Each of the bilayers is created on a square grid; any molecules in excess of #
# the number for a square grid is than placed randomly in the bilayer.         #
#                                                                              #
# Note the generation is slow as it runs the AddToSystem utility to place each #
# individual molecule - presumably, GenLayers will take over one day...        #
################################################################################

source "${HOME}/Scripting/shell/func.sh"

Add=~/AnalysisTools/build/bin/AddToSystem

out=Bilayer.vtf

# test that all used FIELD files exist #{{{
file=( "FA_C16.FIELD" "BTAC.FIELD" "FA_C18.FIELD" "bulk.FIELD")
for i in "${file[@]}"; do
  if ! [ -f "${i}" ]; then
    echo "${i} file does not exist!"
    exit 1
  fi
done #}}}

# find box dimensions from a FIELD file #{{{
# box=$(head -n 1 "${files[0]}" | awk '{print $1}')
boxz=$(head -n 1 "${file[0]}" | awk '{print $3}')
# position of the 'c' bead to place the bilayer is in the box centre
czlo1=$(round "(${boxz}/2-0.4)/${boxz}" 5)
czlo2=$(round "(${boxz}/2-0.3)/${boxz}" 5)
# echo ${czlo1} ${czlo2}
czhi1=$(round "(${boxz}/2+0.3)/${boxz}" 5)
czhi2=$(round "(${boxz}/2+0.4)/${boxz}" 5)
# echo ${czhi1} ${czhi2} #}}}

################################################################################
# create FA_C16 bilayer
################################################################################
#{{{
number=294
field=${file[0]}
echo "##################################################"
echo "# Creating bilayer from ${field}"
echo "##################################################"

nx=$(round "sqrt(${number})" 0)
gapx=$(round "1/${nx}" 5)
ny=$(round "${nx}" 0)
gapy=$(round "1/${ny}" 5)

opt="--head --no-rotate --add -cz ${czhi1} ${czhi2}"
count=0
for (( i=0; i<nx; i++ )); do
  cx1=$(round "${gapx}*(${i}+0.5*1.4)" 5)
  cx2=$(round "${gapx}*(${i}+0.5*1.5)" 5)
  for (( j=0; j<ny; j++ )); do
    count=$(( count+1 ))
    if (( count <= number )); then
      cy1=$(round "${gapy}*(${j}+0.5*1.4)" 5)
      cy2=$(round "${gapy}*(${j}+0.5*1.5)" 5)
      constrain="-cx ${cx1} ${cx2} -cy ${cy1} ${cy2}"
      if [[ ${i} == 0 ]] && [[ ${j} == 0 ]]; then
        run_err "${Add} - ${field} new.vtf ${opt} ${constrain}"
      else
        run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
      fi
      mv {new,old}.vtf
    fi
  done
done
constrain="-bt c cc -ld 0.5"
for (( i=count; i<number; i++)) {
  run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
  mv {new,old}.vtf
}
echo "# ...finished 1st layer"

opt="--head -a 0 0 180 --add -cz ${czlo1} ${czlo2}"
count=0
for (( i=0; i<nx; i++ )); do
  cx1=$(round "${gapx}*(${i}+0.5*0.4)" 5)
  cx2=$(round "${gapx}*(${i}+0.5*0.5)" 5)
  for (( j=0; j<ny; j++ )); do
    count=$(( count+1 ))
    if (( count <= number )); then
      cy1=$(round "${gapy}*(${j}+0.5*0.4)" 5)
      cy2=$(round "${gapy}*(${j}+0.5*0.5)" 5)
      constrain="-cx ${cx1} ${cx2} -cy ${cy1} ${cy2}"
      run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
      mv {new,old}.vtf
    fi
  done
done
constrain="-bt c cc -ld 0.5"
for (( i=count; i<number; i++)) {
  run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
  mv {new,old}.vtf
}
echo "# ...finished 2nd layer"
#}}}
################################################################################
# add BTAC molecules
################################################################################
#{{{
number=293
field="${file[1]}"
echo "##################################################"
echo "# Creating bilayer from ${field}"
echo "##################################################"

nx=$(round "sqrt(${number})" 0)
gapx=$(round "1/${nx}" 5)
ny=$(round "${nx}" 0)
gapy=$(round "1/${ny}" 5)

opt="--head --no-rotate --add -cz ${czhi1} ${czhi2}"
count=0
for (( i=0; i<nx; i++ )); do
  cx1=$(round "${gapx}*(${i}+0.5*0.4)" 5)
  cx2=$(round "${gapx}*(${i}+0.5*0.5)" 5)
  for (( j=0; j<ny; j++ )); do
    count=$(( count+1 ))
    if (( count <= number )); then
      cy1=$(round "${gapy}*(${j}+0.5*0.4)" 5)
      cy2=$(round "${gapy}*(${j}+0.5*0.5)" 5)
      constrain="-cx ${cx1} ${cx2} -cy ${cy1} ${cy2}"
      run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
      mv {new,old}.vtf
    fi
  done
done
constrain="-bt c cc -ld 0.5"
for (( i=count; i<number; i++)) {
  run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
  mv {new,old}.vtf
}
echo "# ...finished 1st layer"

opt="--head -a 0 0 180 --add -cz ${czlo1} ${czlo2}"
count=0
for (( i=0; i<nx; i++ )); do
  cx1=$(round "${gapx}*(${i}+0.5*1.4)" 5)
  cx2=$(round "${gapx}*(${i}+0.5*1.5)" 5)
  for (( j=0; j<ny; j++ )); do
    count=$(( count+1 ))
    if (( count <= number )); then
      cy1=$(round "${gapy}*(${j}+0.5*1.4)" 5)
      cy2=$(round "${gapy}*(${j}+0.5*1.5)" 5)
      constrain="-cx ${cx1} ${cx2} -cy ${cy1} ${cy2}"
      run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
      mv {new,old}.vtf
    fi
  done
done
constrain="-bt c cc -ld 0.5"
for (( i=count; i<number; i++)) {
  run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
  mv {new,old}.vtf
}
echo "# ...finished 2nd layer"
#}}}
################################################################################
# add FA_C18 molecules
################################################################################
#{{{
number=614
field="${file[2]}"
echo "##################################################"
echo "# Creating bilayer from ${field}"
echo "##################################################"

nx=$(round "sqrt(${number})" 0)
gapx=$(round "1/${nx}" 5)
ny=$(round "${nx}" 0)
gapy=$(round "1/${ny}" 5)

opt="--head --no-rotate --add -cz ${czhi1} ${czhi2}"
count=0
for (( i=0; i<nx; i++ )); do
  cx1=$(round "${gapx}*(${i}+0.5*0.9)" 5)
  cx2=$(round "${gapx}*(${i}+0.5*1.1)" 5)
  for (( j=0; j<ny; j++ )); do
    count=$(( count+1 ))
    if (( count <= number )); then
      cy1=$(round "${gapy}*(${j}+0.5*0.9)" 5)
      cy2=$(round "${gapy}*(${j}+0.5*1.1)" 5)
      constrain="-cx ${cx1} ${cx2} -cy ${cy1} ${cy2}"
      run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
      mv {new,old}.vtf
    fi
  done
done
constrain="-bt c cc -ld 0.5"
for (( i=count; i<number; i++)) {
  run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
  mv {new,old}.vtf
}
echo "# ...finished 1st layer"

opt="--head -a 0 0 180 --add -cz ${czlo1} ${czlo2}"
count=0
for (( i=0; i<nx; i++ )); do
  cx1=$(round "${gapx}*(${i}+0.5*0.9)" 5)
  cx2=$(round "${gapx}*(${i}+0.5*1.1)" 5)
  for (( j=0; j<ny; j++ )); do
    count=$(( count+1 ))
    if (( count <= number )); then
      cy1=$(round "${gapy}*(${j}+0.5*0.9)" 5)
      cy2=$(round "${gapy}*(${j}+0.5*1.1)" 5)
      constrain="-cx ${cx1} ${cx2} -cy ${cy1} ${cy2}"
      run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
      mv {new,old}.vtf
    fi
  done
done
constrain="-bt c cc -ld 0.5"
for (( i=count; i<number; i++)) {
  run_err "${Add} old.vtf ${field} new.vtf ${opt} ${constrain}"
  mv {new,old}.vtf
}
echo "# ...finished 2nd layer"
#}}}
################################################################################
# add water and counterions to both sides of the bilayer
################################################################################
#{{{
echo "##################################################"
echo "# Adding solvent into the bulk from ${field}"
echo "##################################################"
field=${file[3]}
# find line number of 'indexed' line
n=$(grep -m 1 indexed old.vtf -n | awk -F ':' '{print $1}')
# grep all z-coordinates and find highes and lowest
n=$(( n+2 ))
loz=$(tail -n +${n} old.vtf | awk '{print $4}' | sort -n | head -n 1)
hiz=$(tail -n +${n} old.vtf | awk '{print $4}' | sort -n | tail -n 1)
# add a bit to the z-coordinates so the solvent cannot lie on top of any bead
loz=$(round "${loz} - 0.3")
hiz=$(round "${hiz} + 0.3")
# add the water - 1st side
opt="--real --add -cz 0 ${loz}"
run_err "${Add} old.vtf ${field} new.vtf ${opt}"
# 2nd side
opt="--real --add -cz ${hiz} ${boxz}"
run_err "${Add} new.vtf bulk.FIELD old.vtf ${opt}"
mv old.vtf ${out}
rm new.vtf
#}}}
