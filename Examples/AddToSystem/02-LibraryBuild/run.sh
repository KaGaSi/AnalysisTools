#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Build a FA_C16/CTAC/water system from a molecule library
################################################################################
lib="${ROOT}/../../library/"
# path to the executable - default assumes AddToSystem is in the $PATH
bin="AddToSystem"

# 1) build system from scratch
${bin} - system.vtf -lib ${lib} -mol FA_C16 50 -ntot 12000 water \
  -b 20 20 20 -sys system.nfo

# 2) add 1% more CTAC molecules to the existing system.
${bin} system.vtf system2.vtf -lib ${lib} -mol CTAC 1% \
  -sys system.nfo system2.nfo
