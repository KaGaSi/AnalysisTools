#!/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# Build a FA_C16 / water system from a molecule library
################################################################################
lib="${ROOT}/../../library/"
bin="AddToSystem"

# Step 1: build system from scratch; -sys writes the system info
${bin} - sys.vtf -lib ${lib} -mol FA_C16 50 -ntot 12000 water \
  -b 20 20 20 -sys sys.nfo

# Step 2: add 1% more FA_C16 molecules to the existing system.
${bin} sys.vtf sys2.vtf -lib ${lib} -mol FA_C16 1% -sys sys.nfo sys2.nfo
