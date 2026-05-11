GenParams assembles force field parameters for a DPD simulation system and
writes them as a FIELD file and/or a LAMMPS parameter file.  Molecules are
loaded from a database directory given by -db (default: ./), which is
auto-detected by the presence of list_molecules.txt:

  library directory  — contains list_molecules.txt; molecule topology and
                       interaction parameters are read from structured files
  FIELD database     — no list_molecules.txt; each molecule is defined by a
                       <name>.FIELD file in the usual DL_MESO format

The examples below use two molecule types drawn from a real CG-DPD parameterisation
(Anderson et al., JCP 147, 094503, 2017; J. Chem. Theory Comput. 14, 2633, 2018):

  water   — single H2O bead (represents two water molecules)
  FA_C16  — 9-bead coarse-grained palmitic acid (C16:0)
  CTAC    — 9-bead cetyltrimethylammonium cation (in examples 02 and 03)
  Cl      — chloride counterion (required alongside CTAC)

01-FieldDatabase/
  Classic mode: one <name>.FIELD file per molecule type, DPD interaction
  parameters specified explicitly in the input file.

02-Library/
  Library mode: topology and interactions come entirely from the library;
  the input file only lists molecule counts.  Note that counterions (Cl)
  must be listed explicitly — they are not added automatically.

03-LammpsParams/
  Extends example 02: adds -lmp to write pair_coeff/bond_coeff/angle_coeff
  lines ready to paste into a LAMMPS input script.
