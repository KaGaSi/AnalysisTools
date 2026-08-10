#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Info examples for generating and analysing all supported formats
################################################################################

# 1) generate other file formats from the System.data LAMMPS data format that
#    contains most information
Info "${ROOT}/System.data" -o System.vsf
Info "${ROOT}/System.data" -o System.vtf
Info "${ROOT}/System.data" -o System.lammpstrj
Info "${ROOT}/System.data" -o System.xyz
Info "${ROOT}/System.data" -o System.FIELD

# 2) take the newly generated files (as well as itp and pdb files Info cannot
#    write) and analyse them, highlighting differences between the file types
Info System.vsf -c System.vtf --verbose > vsf.nfo
Info System.vtf -c System.vtf --verbose > vtf.nfo
Info System.lammpstrj --verbose > ltrj.nfo
Info System.xyz --verbose > xyz.nfo
Info System.FIELD --verbose > field.nfo
Info "${ROOT}/System.itp" --verbose > itp.nfo
Info "${ROOT}/System.pdb" --verbose > pdb.nfo
