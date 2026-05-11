This example builds a CTAC/FA_C16/water system using a molecule library.
The library directory (./library/) is detected automatically because it
contains list_molecules.txt.  All bead/bond/angle types and DPD interaction
parameters are read from the library; system.inp only specifies how many of
each molecule to include.

Compared to example 01, this approach also includes angle potentials (from
list_angles.txt) and uses the published ML-parameterised interaction matrix
rather than hand-tuned values.

Important: GenParams does not auto-add counterions.  The 100 Cl- counterions
for CTAC must be listed as a separate 'molecule Cl 100' line.

Run:

  GenParams system.inp out.FIELD -db ./library/

out.FIELD will contain the full topology (beads, bonds, angles) for all four
molecule types, plus the complete pairwise DPD interaction table.
