This example shows how to build a system from scratch using a molecule library
(-lib/-mol/-ntot), generate a FIELD file with DPD interactions using Info, and
then extend the system by adding more molecules to an existing coordinate file.

The library in ../../library/ is used directly; no per-example copy is needed.

Molecules used:
  FA_C16 — 9-bead coarse-grained palmitic acid
  water  — single-bead solvent (bead type H2O)

Target: 50 FA_C16 molecules (450 beads) + water to fill a 20×20×20 box at
density 3 (12000 beads total → 11550 water beads).

Step 1: build the system (coordinates + topology)

  AddToSystem - -lib ../../library/ -mol FA_C16 50 -ntot 12000 water \
                -b 20 20 20 system.vtf -sys system.sys

  -ntot adds enough water beads (floor((12000-450)/1) = 11550) to reach
  the target count.  Because the fill type is added first internally,
  all 11550 H2O beads are placed with random coordinates.

  -sys system.sys writes a system info file listing molecule type names
  and counts.  When building from scratch ('-'), a single -sys argument
  is output-only.  The file is needed in step 3.

Step 2: add more molecules to the existing system

  AddToSystem system.vtf -lib ../../library/ -mol FA_C16 10 system2.vtf \
              --add -sys system.sys

  When -lib is used with an existing coordinate file (anything other than
  '-'), -sys is required: AddToSystem needs to know the molecule type names
  for the beads already in the system.  Here system.sys (written in step 1)
  provides that mapping so that 10 more FA_C16 molecules can be placed.

Run all steps:

  ./run.sh
