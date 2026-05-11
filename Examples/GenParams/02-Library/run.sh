#!/bin/env bash

################################################################################
# Build a FIELD file for a CTAC/FA_C16/water system using a molecule library.
# list_molecules.txt in ./library/ triggers automatic library-mode detection.
################################################################################

bin=GenParams

${bin} system.inp out.FIELD -db ./library/
