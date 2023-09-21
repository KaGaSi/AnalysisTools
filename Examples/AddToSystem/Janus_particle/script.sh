#!/bin/bash

bin=~/AnalysisTools/build/bin

# 1) create a seed molecule in the middle of a smaller box
#      The system is saved in 1.vcf and 1.vsf files; its box size is 10x5x5,
#      and the -cx/-cy/-cz options along with -gc option ensure placement of
#      of the molecule in the middle of the simulation box. The --no-rotate
#      option ensure the molecule is oriented exactly as defined by the FIELD-1
#      input file.
f1=1
$bin/AddToSystem -- $f1.v{s,c}f -b 10 5 5 -cx 4.99 5.01 -cy 2.49 2.51 -cz 2.49 2.51 -gc -f FIELD-1 --no-rotate
# 2) add first half-shells around the molecule
#      The beads added from FIELD-2 are placed at the distance of at least 0.2
#      and at most 2 from the seed molecule, and the -cx ensures the new beads
#      surround only part of the molecule.
f1=1
f2=2
$bin/AddToSystem $f1.vcf $f2.v{s,c}f -i $f1.vsf -ld 0.2 -hd 2 -bt A -cx 0 4 -f FIELD-2
# 3) add second half-shells around the molecule and save the result into a FIELD-like file as a single molecule
#      Again, the new beads are constrained to be at the distance of at least
#      0.2 and at most 2 from the seed molecule, while the -cx ensures their
#      placement around the still bare part of the molecule. The resulting
#      Janus structure is saved into mol.FIELD file as a single molecule (but
#      with no bonds) to allow placing the structure into other systems.
f1=2
f2=3
$bin/AddToSystem $f1.vcf $f2.v{s,c}f -i $f1.vsf -ld 0.2 -hd 2 -bt A -cx 4.2 10 -f FIELD-3 -f_out mol.FIELD 1
# 4) populate a box with structures from 3) via for loop to ensure they do not overlap.
#      At last, a new system is created in 4.vcf and 4.vsf files with box size
#      20x20x20 that contains 10 of the Janus structures. The structures are
#      added one by one with the first being placed completely at random and
#      the remaining nine added using a for loop at a distance of at least 6
#      from the centre of other structures to ensure no overlap between them.
f2=a
$bin/AddToSystem -- $f2.v{s,c}f -b 20 20 20 -gc -f mol.FIELD
for (( i=0 ; i<9 ; i++ )); do
  f1=a
  f2=b
  $bin/AddToSystem $f1.vcf $f2.v{s,c}f -i $f1.vsf -ld 6 -bt A B C -gc -f mol.FIELD
  mv $f2.vsf $f1.vsf
  mv $f2.vcf $f1.vcf
done
f2=4
mv $f1.vsf $f2.vsf
mv $f1.vcf $f2.vcf
