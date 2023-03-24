#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#define EXTENSION 16

#include "Errors.h"
#include "General.h"
#include "Options.h"
#include "Read.h"
#include "Structs.h"
#include "Write.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>

#define VSF_FILE 1
#define VCF_FILE 2
#define XYZ_FILE 3
#define LDATA_FILE 4
#define LTRJ_FILE 5
#define FIELD_FILE 6
#define CONFIG_FILE 7

// Helper functions for dealing with SYSTEM structure
// fill in some SYSTEM stutff
void FillMoleculeTypeBType(MOLECULETYPE *MoleculeType);
void FillMoleculeTypeChargeMass(MOLECULETYPE *MoleculeType,
                                BEADTYPE BeadType[]);
void FillBeadTypeIndex(SYSTEM *System);
void FillMoleculeTypeIndex(SYSTEM *System);
void CountBondAngleDihedralImproper(SYSTEM *System);
// Appends _# to bead/molecule types with the same name
void RenameBeadTypes(SYSTEM *System);
void RenameMoleculeTypes(SYSTEM *System);
// get bead indices for bonds/angles/dihedrals (with some error checking)
int *BondIndices(SYSTEM System, int mol, int bond);
int *AngleIndices(SYSTEM System, int mol, int angle);
int *DihedralIndices(SYSTEM System, int mol, int dihed);
int *ImproperIndices(SYSTEM System, int mol, int dihed);
// enrich molecule types with information from a second System structure
void ChangeMolecules(SYSTEM *Sys_orig, SYSTEM Sys_add, bool beads, bool name);
// test whether two bead types are identical
bool SameBeadType(BEADTYPE bt_1, BEADTYPE bt_2);
// create new bead/molecule type, realloc'ing the appropriate array
void NewBeadType(BEADTYPE *BeadType[], int *number_of_types, char name[],
                 double charge, double mass, double radius);
void NewMolType(MOLECULETYPE *MoleculeType[], int *n_types, char name[],
                int n_beads, int n_bonds, int n_angles, int n_dihedrals,
                int n_impropers);
// identify bead type based on name
int FindBeadType(char name[], SYSTEM System);
// identify molecule type based on name
int FindMoleculeName(char name[], SYSTEM System);
// identify molecule type based on name only or on other parameters too
int FindMoleculeType(SYSTEM Sys1, MOLECULETYPE mt, SYSTEM Sys2,
                     int mode, bool name);
// cleanse System by removing molecule/bead types with .Number=0, etc.
void PruneSystem(SYSTEM *System);
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box);
// copy molecule type
MOLECULETYPE CopyMoleculeType(MOLECULETYPE mt_old);
MOLECULETYPE CopyMoleculeTypeEssentials(MOLECULETYPE mt_old);
// check that the System struct doesn't contain an error
void CheckSystem(SYSTEM System, char file[]);

// Helper functions for manipulating coordinates
// wrap coordinates into simulation box and/or join molecules
void WrapJoinCoordinates(SYSTEM *System, bool wrap, bool join);
// distance between two beads; in the range <-BoxLength/2,BoxLength/2)
VECTOR Distance(VECTOR id1, VECTOR id2, VECTOR BoxLength);

// identify input coordinate and structure files
bool InputCoorStruct(int argc, char *argv[], char coor[], int *coor_type,
                     char struc[], int *struc_type);

// create a cell-linked list
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead, int **Head,
                int **Link, double cell_size, INTVECTOR *n_cells, int *Dcx,
                int *Dcy, int *Dcz);

// verbose output (print various structures and some such)
void VerboseOutput(SYSTEM System);
void PrintCount(COUNT Count);
void PrintBeadType(SYSTEM System);
void PrintMoleculeType(SYSTEM System);
void PrintMolecule(SYSTEM System);
void PrintBead(SYSTEM System);
void PrintBondType(SYSTEM System);
void PrintAngleType(SYSTEM System);
void PrintDihedralType(SYSTEM System);
void PrintImproperType(SYSTEM System);
// TODO: use SYSTEM
void PrintBondTypes(COUNTS Counts, PARAMS *bond_type);
// TODO: use SYSTEM
void PrintAngleTypes(COUNTS Counts, PARAMS *angle_type);
void PrintBox(BOX Box);
void PrintByline(FILE *ptr, int argc, char *argv[]);
void PrintStep(int *count_coor, int start, bool silent);

// calculate gyration tensor and various shape descriptors
VECTOR Gyration(int n, int *list, COUNTS Counts, BEADTYPE *BeadType,
                BEAD **Bead);

// memory-freeing functions
void FreeSystem(SYSTEM *System);
void FreeMoleculeType(MOLECULETYPE *MoleculeType);
void FreeMoleculeTypeEssentials(MOLECULETYPE *MoleculeType);

// TODO redo
void EvaluateContacts(COUNTS *Counts, AGGREGATE **Aggregate,
                      MOLECULE **Molecule, int contacts, int **contact);
void SortAggStruct(AGGREGATE **Aggregate, COUNTS Counts, MOLECULE *Molecule,
                   MOLECULETYPE *MoleculeType, BEAD **Bead,
                   BEADTYPE *BeadType);

#if 0 //{{{
// TODO redo
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNTS Counts,
                         VECTOR BoxLength, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule);
void PrintAggregate(COUNTS Counts, int *Index, MOLECULETYPE *MoleculeType,
                    MOLECULE *Molecule, BEAD *Bead, BEADTYPE *BeadType,
                    AGGREGATE *Aggregate);
void FreeAggregate(COUNTS Counts, AGGREGATE **Aggregate);
#endif //}}}
#endif
