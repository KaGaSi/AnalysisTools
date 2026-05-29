#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Info examples for generating and analysing all supported formats
################################################################################
INFO="${ROOT}/../../../build/bin/Info"

# 1) generate other file formats from the System.data LAMMPS data format that
#    contains most information
"${INFO}" "${ROOT}/System.data" -o System.vsf
"${INFO}" "${ROOT}/System.data" -o System.vtf
"${INFO}" "${ROOT}/System.data" -o System.lammpstrj
"${INFO}" "${ROOT}/System.data" -o System.xyz
"${INFO}" "${ROOT}/System.data" -o System.FIELD

# 2) take the newly generated files (as well as itp and pdb files Info cannot
#    write) and analyse them, highlighting differences between the file types
"${INFO}" System.vsf -c System.vtf --verbose > vtf.nfo
"${INFO}" System.vtf -c System.vtf --verbose > vtf.nfo
"${INFO}" System.lammpstrj --verbose > ltrj.nfo
"${INFO}" System.xyz --verbose > xyz.nfo
"${INFO}" System.FIELD --verbose > field.nfo
"${INFO}" "${ROOT}/System.itp" --verbose > itp.nfo
"${INFO}" "${ROOT}/System.pdb" --verbose > pdb.nfo
