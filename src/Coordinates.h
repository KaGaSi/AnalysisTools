#ifndef COORDINATES_H
#define COORDINATES_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"
#include <stdbool.h>

// put given vector into range <0,BoxLength)
vec3d RestorePBC(const vec3d coor, const vec3d BoxLength);
// convert BeadCoor positions to/from OrthoLength-scaled fractional coordinates
void CoorToFractional(SYSTEM *System);
void CoorFromFractional(SYSTEM *System);
// remove pbc for a single molecule by joining it
void RemovePBCMolecule(int mol_id, SYSTEM *System);
// wrap coordinates into simulation box and/or join molecules
void WrapJoinCoordinates(SYSTEM *System, const bool wrap, const bool join);
// physical minimum-image distance for orthogonal and triclinic boxes
vec3d DistancePBC(const vec3d r1, const vec3d r2, const BOX *box);
// distance between two beads; in the range <-BoxLength/2,BoxLength/2)
vec3d Distance(const vec3d id1, const vec3d id2, const vec3d BoxLength);
// calculate centre of mass for a list of beads
vec3d CentreOfMass(const int n, const int *list, const SYSTEM System);
// calculate geometric centre for a list of beads
vec3d GeomCentre(const int n, const int *list, const BEAD *Bead);
// calculate gyration tensor eigenvalues (shape descriptors)
vec3d Gyration(const int n, const int *list, SYSTEM *System);

#endif
