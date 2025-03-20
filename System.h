#ifndef SYSTEM_H
#define SYSTEM_H

#include "AnalysisTools.h"

// fill in some SYSTEM arrays and some such
void FillMoleculeTypeBType(MOLECULETYPE *MoleculeType);
void ReFillMoleculeTypeBType(SYSTEM *System);
void FillMoleculeTypeChargeMass(MOLECULETYPE *MoleculeType,
                                BEADTYPE BeadType[]);
void FillBeadTypeIndex(SYSTEM *System);
void AllocFillBeadTypeIndex(SYSTEM *System);
void RefillBeadTypeIndex(SYSTEM *System);
void FillMoleculeTypeIndex(SYSTEM *System);
void ReFillMoleculeTypeIndex(SYSTEM *System);
// void FillIndexMol(SYSTEM *System);
void FillBondedUnbonded(SYSTEM *System);
void CountBondAngleDihedralImproper(SYSTEM *System);
void SortBonds(int (*bond)[3], int n);
void SortAngles(int (*angle)[4], int n);
void SortDihImp(int (*dihimp)[5], int n);
void SortAll(MOLECULETYPE *mt);
void FillSystemNonessentials(SYSTEM *System, bool bonds);
void FillInCoor(SYSTEM *System);
bool CalculateBoxData(BOX *Box, int mode);

// merge identical bead/molecule types
void MergeBeadTypes(SYSTEM *System, bool detailed);
void MergeMoleculeTypes(SYSTEM *System);

// Appends # to bead/molecule types with the same name
void RenameBeadTypes(SYSTEM *System);
void RenameMoleculeTypes(SYSTEM *System);

// test whether two bead types are identical
bool SameBeadType(const BEADTYPE bt_1, const BEADTYPE bt_2, const bool name);

// create new bead/molecule type, realloc'ing the appropriate array
void NewBeadType(BEADTYPE *BeadType[], int *number_of_types, char name[],
                 double charge, double mass, double radius);
void NewMolType(MOLECULETYPE *MoleculeType[], int *n_types, char name[],
                int n_beads, int n_bonds, int n_angles, int n_dihedrals,
                int n_impropers);

// copy System
MOLECULETYPE CopyMoleculeType(MOLECULETYPE mt_old);
MOLECULETYPE CopyMoleculeTypeEssentials(MOLECULETYPE mt_old);
SYSTEM CopySystem(SYSTEM S_in);

// cleanse System by removing molecule/bead types with .Number=0, etc.
void PruneSystem(SYSTEM *System);
// join two systems, possibly pruning it
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box, bool prune);

// check that the System struct doesn't contain an error
void CheckSystem(const SYSTEM System, const char *file);

// simplify system for vtf output - remove stuff vtf does not support
void VtfSystem(SYSTEM *System);

// finish filling System after essentials are filled
void FinishSystem(SYSTEM *System);

// enrich molecule types with information from a second System structure
void ChangeMolecules(SYSTEM *Sys_orig, SYSTEM Sys_add, bool name);

// add/subtract Box.Low to/from coordinates
void ChangeBoxByLow(SYSTEM *System, int sign);

void SortAggStruct(AGGREGATE *Aggregate, SYSTEM System);

// memory-freeing functions
void FreeSystem(SYSTEM *System);
void FreeMoleculeType(MOLECULETYPE *MoleculeType);
void FreeMoleculeTypeEssentials(MOLECULETYPE *MoleculeType);
void FreeAggregate(COUNT Count, AGGREGATE *Aggregate);

// realloc some System.*{,Coor} arrays
void ReallocBead(SYSTEM *System);
void ReallocBonded(SYSTEM *System);
void ReallocUnbonded(SYSTEM *System);
void ReallocMolecule(SYSTEM *System);
#endif
