#ifndef COORDINATES_H
#define COORDINATES_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"
#include <stdbool.h>

// put given vector into range <0,BoxLength)
vec3d RestorePBC(const vec3d coor, const vec3d BoxLength);
// remove pbc for a single molecule by joining it
void RemovePBCMolecule(int mol_id, SYSTEM *System);
// wrap coordinates into simulation box and/or join molecules
void WrapJoinCoordinates(SYSTEM *System, const bool wrap, const bool join);
// distance between two beads; in the range <-BoxLength/2,BoxLength/2)
vec3d Distance(const double id1[3], const double id2[3],
               const vec3d BoxLength);
// calculate centre of mass for a list of beads
void CentreOfMass(const int n, const int *list,
                  const SYSTEM System, double com[3]);
// calculate geometric centre for a list of beads
void GeomCentre(const int n, const int *list, const BEAD *Bead, double gc[3]);
// calculate gyration tensor eigenvalues (shape descriptors)
vec3d Gyration(const int n, const int *list, SYSTEM *System);

#endif
