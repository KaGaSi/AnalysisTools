#!/bin/env bash

################################################################################
# Build a FIELD file for 50 FA_C16 molecules in water.
# Molecule definitions are loaded from ./molecules/<name>.FIELD
################################################################################

bin=GenParams

${bin} system.inp out.FIELD -db ./molecules/
