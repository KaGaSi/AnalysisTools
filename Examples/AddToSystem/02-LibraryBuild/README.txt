This example shows how to build a system from scratch using a molecule library
(-lib/-mol/-ntot) and then extend it by adding more molecules to the resulting
coordinate file.

The library in ../../library/ is used directly; no per-example copy is needed.

Molecules used:
  FA_C16 - 9-bead coarse-grained fatty alcohol
  water  - single-bead solvent (bead type H2O)

Target: 50 FA_C16 molecules (450 beads) plus water to reach 12000 beads in a
20x20x20 box, i.e. 11550 water beads.

1) build the system (coordinates and topology):
AddToSystem - system.vtf -lib ../../library/ -mol FA_C16 50 \
  -b 20 20 20 -ntot 12000 water -sys system.nfo

  '-' (instead of input coordinate file) ... build system from scratch
  -lib ... relative path to the library of molecules
  -b ... system's orthogonal box size
  -mol ... molecule name from the library and its count
  -ntot ... add enough 'water' beads to reach the target count (11550 here)
  -sys ... write the list of molecule type names and counts for use in 2)

2) add new molecules to the existing system:
AddToSystem system.vtf system2.vtf -lib ../../library/ -mol CTAC 1% \
  -sys system.nfo system2.nfo

  -mol ... same as before but a percentage of system's bead number instead of
           number of molecules (leading to 0.01*12000 beads, so 12 molecules)
  -sys ... input system (required when using -lib with pre-existing system) and
           output system descriptions
  Note that --add is not used, so new beads will be exchanged for existing
    'water' beads, keeping total bead count constant.
  Note that the CTAC Cl- counterion is part considered part of the molecule
    (check via Info utility).

Run all steps:

  ./run.sh
