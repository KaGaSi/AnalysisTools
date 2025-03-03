#!/bin/bash

bin=../../build/bin/

${bin}/Selected ../AngleMolecules/traj.vcf traj.vtf -i ../AngleMolecules/traj.data

# 1) Using all three options in one command (shows default behaviour):
${bin}/Selected traj.vtf full.vtf -cx 0 0.2 0.6 0.8 -cy 0 0.2 0.5 0.9 -sc 2 -m 1 1 0

${bin}/Selected traj.vtf cx_cy.vtf -cx 0 0.2 0.6 0.8 -cy 0 0.2 0.5 0.9
${bin}/Selected cx_cy.vtf cx_cy-m.vtf -m 1 1 0
${bin}/Selected cx_cy-m.vtf cx_cy-m-sc.vtf -sc 2

# 2) First move beads and then scale coordinates - the first command uses
#    cx_cy.vtf from 1):
${bin}/Selected cx_cy.vtf cx_cy-sc.vtf -sc 2
${bin}/Selected cx_cy-sc.vtf cx_cy-sc-m.vtf -m 1 1 0

# 3) First scale coordinates and then constrain the saved ones - the first
${bin}/Selected traj.vtf sc2.vtf -sc 2
${bin}/Selected sc2.vtf sc2-cx_cy.vtf -cx 0 0.2 0.6 0.8 -cy 0 0.2 0.5 0.9
