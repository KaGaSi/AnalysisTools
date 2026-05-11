This example extends 02-Library by also writing a LAMMPS parameter file.
The -lmp flag produces pair_coeff, bond_coeff, and angle_coeff lines that
can be included directly in a LAMMPS input script.

The library from ../02-Library/library/ is reused; no separate copy is needed.

Run:

  GenParams system.inp out.FIELD -db ../02-Library/library/ -lmp params.in

out.FIELD    — DL_MESO FIELD file with full topology and interaction table
params.in    — LAMMPS coeff lines, e.g.:
                 pair_coeff  1  1 dpd 25.00000 ${gamma} 1.00000 # H2O - H2O
                 pair_coeff  1  2 dpd 14.50000 ${gamma} 0.99000 # H2O - CH2OH
                 ...
                 bond_coeff  1 75.0 0.29
                 ...
                 angle_coeff 1 2.5 180.0

The ${gamma} placeholder is left for substitution in the LAMMPS script since
gamma is typically set globally (e.g. 'variable gamma equal 4.5').
