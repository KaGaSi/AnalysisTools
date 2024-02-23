#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#include "Errors.h"
#include "General.h"
#include "Options.h"
#include "Read.h"
#include "Structs.h"
#include "Write.h"
// #include <dirent.h>
// #include <sys/stat.h>

// Helper functions for dealing with SYSTEM structure
// fill in some SYSTEM stutff
void FillMoleculeTypeBType(MOLECULETYPE *MoleculeType);
void ReFillMoleculeTypeBType(SYSTEM *System);
void FillMoleculeTypeChargeMass(MOLECULETYPE *MoleculeType,
                                BEADTYPE BeadType[]);
void FillBeadTypeIndex(SYSTEM *System);
void ReFillBeadTypeIndex(SYSTEM *System);
void FillMoleculeTypeIndex(SYSTEM *System);
void ReFillMoleculeTypeIndex(SYSTEM *System);
void FillIndexMol(SYSTEM *System);
void FillBondedUnbonded(SYSTEM *System);
void CountBondAngleDihedralImproper(SYSTEM *System);
void SortBonds(int (*bond)[3], int n);
void SortAngles(int (*angle)[4], int n);
void SortDihImp(int (*dihimp)[5], int n);
void FillSystemNonessentials(SYSTEM *System);
void FillInCoor(SYSTEM *System);
bool CalculateBoxData(BOX *Box, int mode);
// merge identical bead/molecule types
void MergeBeadTypes(SYSTEM *System, bool detailed);
void MergeMoleculeTypes(SYSTEM *System);
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
// TODO: CopySystem won't be static?
SYSTEM CopySystem(SYSTEM S_in);
// cleanse System by removing molecule/bead types with .Number=0, etc.
void PruneBondTypes(SYSTEM S_old, SYSTEM *System);
void PruneAngleTypes(SYSTEM S_old, SYSTEM *System);
void PruneDihedralTypes(SYSTEM S_old, SYSTEM *System);
void PruneImproperTypes(SYSTEM S_old, SYSTEM *System);
void PruneSystem(SYSTEM *System);
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box);
// copy molecule type
MOLECULETYPE CopyMoleculeType(MOLECULETYPE mt_old);
MOLECULETYPE CopyMoleculeTypeEssentials(MOLECULETYPE mt_old);
// check that the System struct doesn't contain an error
void CheckSystem(SYSTEM System, char file[]);
// simplify system for vtf output - remove stuff vtf does not support
void VtfSystem(SYSTEM *System);

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
// add/subtract Box.Low to/from coordinates
void AddLow(SYSTEM *System);
void SubtractLow(SYSTEM *System);

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
void PrintBondTypes(COUNT Counts, PARAMS *bond_type);
// TODO: use SYSTEM
void PrintAngleTypes(COUNT Counts, PARAMS *angle_type);
void PrintBox(BOX Box);
void PrintByline(char *file, int argc, char *argv[]);
void PrintStep(int *count_coor, int start, bool silent);

// calculate gyration tensor and various shape descriptors
void Gyration(int n, int *list, COUNT Counts, BEADTYPE *BeadType,
              BEAD **Bead, double eigen[3]);

// memory-freeing functions
void FreeSystem(SYSTEM *System);
void FreeMoleculeType(MOLECULETYPE *MoleculeType);
void FreeMoleculeTypeEssentials(MOLECULETYPE *MoleculeType);

// TODO redo
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      int contacts, int **contact);
void SortAggStruct(AGGREGATE *Aggregate, SYSTEM System);

// void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNT Counts,
//                          VECTOR BoxLength, BEADTYPE *BeadType, BEAD *Bead,
//                          MOLECULETYPE *MoleculeType, MOLECULE *Molecule);
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate,
                         SYSTEM *System);
void FreeAggregate(COUNT Count, AGGREGATE *Aggregate);

void PrintAggregate(SYSTEM System, AGGREGATE Aggregate[]);
#if 0 //{{{
// TODO redo
void PrintAggregate(COUNT Counts, int *Index, MOLECULETYPE *MoleculeType,
                    MOLECULE *Molecule, BEAD *Bead, BEADTYPE *BeadType,
                    AGGREGATE *Aggregate);
#endif //}}}
#endif
