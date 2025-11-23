#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#define _POSIX_C_SOURCE 200809L

#include "pcg32.h"
#include "Arrays.h"
#include "Structs.h"
#include "Pairs.h"
#include "Globals.h"
#include "MathUtils.h"
#include "System.h"
#include "Errors.h"
#include "General.h"
#include "Options.h"
#include "ReadWrite.h"
#include "Aggregates.h"
#include <complex.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_poly.h> // solve cubic equation
#include <gsl/gsl_eigen.h>  // Jacobi method
#include <gsl/gsl_matrix.h> //
#include <gsl/gsl_vector.h> //

// Helper functions for dealing with SYSTEM structure
// identify bead type based on name
int FindBeadType(const char *name, const SYSTEM System);
// identify molecule type based on name
int FindMoleculeName(const char *name, const SYSTEM System);
// identify molecule type based on name only or on other parameters too
int FindMoleculeType(const SYSTEM Sys1, const MOLECULETYPE mt,
                     const SYSTEM Sys2, const int mode, const bool name);
// Helper functions for manipulating coordinates
// wrap coordinates into simulation box and/or join molecules
void WrapJoinCoordinates(SYSTEM *System, const bool wrap, const bool join);
// distance between two beads; in the range <-BoxLength/2,BoxLength/2)
vec3d Distance(const double id1[3], const double id2[3],
              const double BoxLength[3]);
// calculate centre of mass for a list of beads
void CentreOfMass(const int n, const int *list,
                  const SYSTEM System, double gc[3]);
// calculate geometric centre for a list of beads
void GeomCentre(const int n, const int *list, const BEAD *Bead, double gc[3]);
// identify input coordinate and structure files
bool InputCoorStruct(const int argc, char **argv, SYS_FILES *f);
int StructureFileType(const char *name);
int CoordinateFileType(const char *name);
int FileType(const char *name);
// create a cell-linked list
void LinkedList_old(const SYSTEM System, int **Head, int **Link,
                const double cell_size, int n_cells[3], int Dc[27][3]);
int SelectCell1_old(const int c1[3], const int n_cells[3]);
int SelectCell2_old(const int c1[3], const int n_cells[3],
                const int Dc[27][3], const int n);
// calculate gyration tensor and various shape descriptors
void Gyration(const int n, const int *list, SYSTEM *System, double eigen[3]);
// TODO: redo
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      const int contacts, ArrNDi *contact);
void RemovePBCAggregates(const double distance, const AGGREGATE *Aggregate,
                         SYSTEM *System);
bool UseStep(const COMMON_OPT opt, int step);
#endif
