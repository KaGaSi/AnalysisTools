#!/bin/env bash

################################################################################
# Build both a FIELD file and a LAMMPS parameter file for the same system.
# Reuses the molecule library from ../02-Library/library/
################################################################################

bin=GenParams

${bin} system.inp out.FIELD -db ../02-Library/library/ -lmp params.in
