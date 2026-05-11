This example builds a palmitic acid (FA_C16) / water system using a small
FIELD-file database.  Each molecule type is defined by a separate .FIELD file
in the molecules/ directory, one instance per file (nummols 1).  DPD interaction
parameters are given explicitly in system.inp.

The 'fill' keyword in system.inp computes the number of water beads needed to
reach the target density (3) in the 20x20x20 box after placing the 50 FA_C16
molecules.

Run:

  GenParams system.inp out.FIELD -db ./molecules/

out.FIELD will contain the full molecule topology followed by an 'interactions'
block with all pairwise DPD parameters.
