#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"
#include "System.h"
#include "Errors.h"
#include "General.h"
#include "Options.h"
#include "Read.h"
#include "Write.h"
// #include <dirent.h>
// #include <sys/stat.h>

// Helper functions for dealing with SYSTEM structure
// identify bead type based on name
int FindBeadType(char name[], SYSTEM System);
// identify molecule type based on name
int FindMoleculeName(char name[], SYSTEM System);
// identify molecule type based on name only or on other parameters too
int FindMoleculeType(SYSTEM Sys1, MOLECULETYPE mt, SYSTEM Sys2, int mode);

// Helper functions for manipulating coordinates
// wrap coordinates into simulation box and/or join molecules
void WrapJoinCoordinates(SYSTEM *System, bool wrap, bool join);
// distance between two beads; in the range <-BoxLength/2,BoxLength/2)
void Distance(double id1[3], double id2[3],
              double BoxLength[3], double out[3]);
// calculate centre of mass for a list of beads
void CentreOfMass(int n, int list[], SYSTEM System, double gc[3]);
// calculate geometric centre for a list of beads
void GeomCentre(int n, int *list, BEAD *Bead, double gc[3]);
void ToFractionalCoor(SYSTEM *System);
void FromFractionalCoor(SYSTEM *System);

// get bead indices for bonds/angles/dihedrals (with some error checking)
int *BondIndices(SYSTEM System, int mol, int bond);
int *AngleIndices(SYSTEM System, int mol, int angle);
int *DihedralIndices(SYSTEM System, int mol, int dihed);
int *ImproperIndices(SYSTEM System, int mol, int dihed);

// identify input coordinate and structure files
// bool InputCoorStruct(int argc, char *argv[], char coor[], int *coor_type,
//                      char struc[], int *struc_type);
bool InputCoorStruct(int argc, char *argv[], SYS_FILES *f);
// identify type of provided structure file (mode=0: input, mode=1 output file)
int StructureFileType(char name[]);
int CoordinateFileType(char name[]);
int FileType(char name[]);

// create a cell-linked list
void LinkedList(SYSTEM System, int **Head, int **Link, double cell_size,
                int n_cells[3], int Dc[14][3]);
int SelectCell1(int c1[3], int n_cells[3]);
int SelectCell2(int c1[3], int n_cells[3], int Dc[14][3], int n);

// calculate gyration tensor and various shape descriptors
void Gyration(int n, int *list, COUNT Counts, BEADTYPE *BeadType,
              BEAD **Bead, double eigen[3]);

// TODO: redo
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      int contacts, int **contact);

// void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNT Counts,
//                          VECTOR BoxLength, BEADTYPE *BeadType, BEAD *Bead,
//                          MOLECULETYPE *MoleculeType, MOLECULE *Molecule);
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate,
                         SYSTEM *System);

bool UseStep(COMMON_OPT opt, int step);
#endif
