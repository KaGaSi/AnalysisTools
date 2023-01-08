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

VECTOR ToFractional(VECTOR coor, BOX Box);
void ToFractionalCoor(SYSTEM *System);
VECTOR FromFractional(VECTOR coor, BOX Box);
void FromFractionalCoor(SYSTEM *System);

// get bead indices for bonds/angles/dihedrals (with some error checking)
int * BondIndices(SYSTEM System, int mol, int bond);
int * AngleIndices(SYSTEM System, int mol, int angle);
int * DihedralIndices(SYSTEM System, int mol, int dihed);
int * ImproperIndices(SYSTEM System, int mol, int dihed);

void ChangeMolecules(SYSTEM *Sys_orig, SYSTEM Sys_add, bool beads, bool name);
void ChangeMolecules_old(SYSTEM *Sys_orig, SYSTEM Sys_add, bool beads, bool name);

void CheckSystem(SYSTEM System, char file[]);

// InputCoor() //{{{
/**
 * \brief Function test input coordinate file is correct
 *
 * \param [out] vtf          is the coordinate file vtf or vcf?
 * \param [in]  file_coor    name of coordinate file
 * \param [in]  file_struct  name of structure file
 * \return false if file_coor has wrong extension, true otherwise
 */
bool InputCoor_old(bool *vtf, char *file_coor, char *file_struct); //}}}
int InputCoorStruct_old(int argc, char **argv,
                    char coor[], char vsf[], char lmp[], char field[]);
bool InputCoorStruct(int argc, char **argv, char coor[], int *coor_type,
                     char struc[], int *struc_type);

bool SameBeadType(BEADTYPE bt_1, BEADTYPE bt_2);
int FindBeadType(char name[], SYSTEM System);
int FindMoleculeName(char name[], SYSTEM System);
int FindMoleculeType_old(MOLECULETYPE mol, SYSTEM System, int mode, bool name);
int FindMoleculeType(SYSTEM Sys1, int mt, SYSTEM Sys2, int mode, bool name);

void FillSystemNonessentials(SYSTEM *System);
void FillMoleculeTypeBType(MOLECULETYPE *MoleculeType);
void FillMoleculeTypeChargeMass(MOLECULETYPE *MoleculeType,
                                BEADTYPE BeadType[]);
void FillBeadTypeIndex(SYSTEM *System);
void FillMoleculeTypeIndex(SYSTEM *System);
bool TriclinicCellData(BOX *Box, int mode);
void CountBondAngleDihedralImproper(SYSTEM *System);

void NewBeadType(BEADTYPE *BeadType[], int *number_of_types, char *name,
                 double charge, double mass, double radius);
void NewMolType(MOLECULETYPE *MoleculeType[], int *n_types, char *name,
                int n_beads, int n_bonds, int n_angles, int n_dihedrals,
                int n_impropers);

void RemovePBCMolecules(SYSTEM *System);
void RestorePBC(SYSTEM *System);

// Distancet() //{{{
/**
 * \brief Function to calculate distance vector between two beads.
 *
 * \param [in] id1         first coordinate vector
 * \param [in] id2         second coordinate vector
 * \param [in] BoxLength   dimensions of simulation box
 * \return distance vector between the two provided beads (without pbc)
 */
VECTOR Distance(VECTOR id1, VECTOR id2, VECTOR BoxLength); //}}}
// CentreOfMass() //{{{
/**
 * \brief Function to calculate centre of mass for a collection of beads.
 *
 * \param [in] n              number of beads
 * \param [in] list           list of bead ids (corresponding to indices in Bead
 * struct) \param [in] Bead           information about individual beads
 * (coordinates) \param [in] BeadType       information about beadtypes (masses)
 * \return coordinates of centre of mass of a given aggregate
 */
VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType); //}}}
// GeomCentre() //{{{
/**
 * \brief Function to calculate geometric centre for a collection of beads.
 *
 * \param [in] n              number of beads
 * \param [in] list           list of bead ids (corresponding to indices in Bead
 * struct) \param [in] Bead           information about individual beads
 * (coordinates) \return coordinates of geometric centre of a given aggregate
 */
VECTOR GeomCentre(int n, int *list, BEAD *Bead); //}}}
// Gyration() //{{{
/*
 * \brief Function calculating principal moments of the gyration tensor.
 *
 * \param [in] n             number of beads
 * \param [in] list          array of bead ids
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BoxLength     dimensions of simulation box
 * \param [in] BeadType      informationn about bead types
 * \param [in] Bead          informationn about individual beads
 * \return vector with principal moments of gyration tensor (sorted as x<y<z)
 */
VECTOR Gyration(int n, int *list, COUNTS Counts, BEADTYPE *BeadType,
                BEAD **Bead); //}}}

// LinkedList() //{{{
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead, int **Head,
                int **Link, double cell_size, INTVECTOR *n_cells, int *Dcx,
                int *Dcy, int *Dcz); //}}}

void SortBonds(int (*bond)[3], int num);
void SortAngles(int (*angle)[4], int num);
void SortDihImp(int (*dihedral)[5], int num);

SYSTEM CopySystem(SYSTEM Sys_in);
MOLECULETYPE CopyMoleculeType(MOLECULETYPE mt_old);
MOLECULETYPE CopyMoleculeTypeEssentials(MOLECULETYPE mt_old);
void PruneSystem(SYSTEM *System);
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box);

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

void FreeSystem(SYSTEM *System);
void FreeMoleculeType(MOLECULETYPE *MoleculeType);
void FreeMoleculeTypeEssentials(MOLECULETYPE *MoleculeType);

#if 0 //{{{
// TODO remove
// PrintAggregate() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintAggregate(COUNTS Counts, int *Index, MOLECULETYPE *MoleculeType,
                    MOLECULE *Molecule, BEAD *Bead, BEADTYPE *BeadType,
                    AGGREGATE *Aggregate); //}}}
void PrintBondTypes2(int number_of_bonds, PARAMS *bond_type);
// RestorePBC_old() //{{{
/**
 * \brief Function to restore pbc.
 *
 * \param [in]  Counts         numbers of beads, molecules, etc.
 * \param [in]  BoxLength      dimension of the simulation box
 * \param [out] Bead           information about individual beads (coordinates)
 */
void RestorePBC_old(COUNTS Counts, VECTOR BoxLength, BEAD **Bead); //}}}
// EvaluateContacts() //{{{
/**
 * \brief Function evaluating contacts for aggregate detection
 *
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [out] Aggregate     information about aggregates
 * \param [in]  Molecule      information about individual molecules
 * \param [in]  contacts      number of contacts for aggregate check
 * \param [in]  contact       2D array containing number of contacts between
 * molecules
 */
void EvaluateContacts(COUNTS *Counts, AGGREGATE **Aggregate,
                      MOLECULE **Molecule, int contacts, int **contact);
//}}}
// SortAggStruct() //{{{
/**
 * \brief Function to sort Aggregate struct.
 *
 * \param [out] Aggregate  Aggregate struct to sort
 * \param [in]  Counts     numbers of beads, molecules, etc.
 */
void SortAggStruct(AGGREGATE **Aggregate, COUNTS Counts, MOLECULE *Molecule,
                   MOLECULETYPE *MoleculeType, BEAD **Bead,
                   BEADTYPE *BeadType); //}}}
void PrintAngleTypes2(int number_of_angles, PARAMS *angle_type);
void PrintDihedralTypes2(int number_of_dihedrals, PARAMS *dihedral_type);
// VerboseOutput_old() //{{{
/**
 * \brief Function printing basic information about system if `-v` or `-V`
 * option is provided
 *
 * \param [in] input_vcf     .vcf coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BoxLength     dimension of the simulation box
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          informationn about individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      information about individual molecules
 */
void VerboseOutput_old(char *input_vcf, COUNTS Counts, VECTOR BoxLength,
                       BEADTYPE *BeadType, BEAD *Bead,
                       MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}
// VerboseOutput_oldish() //{{{
/**
 * \brief Function printing basic information about system if `-v` or `-V`
 * option is provided
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          informationn about individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      information about individual molecules
 */
void VerboseOutput_oldish(COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}
// FreeSystemInfo() //{{{
/**
 * \brief Free memory for all standard arrays and structures of arrays.
 *
 * \param [in]  Counts         numbers of beads, molecules, and types
 * \param [out] MoleculeType   information about molecule types
 * \param [out] Molecule       information about individual molecules
 * \param [out] BeadType       information about bead types
 * \param [out] Bead           information about individual beads
 */
void FreeSystemInfo(COUNTS Counts, MOLECULETYPE **MoleculeType,
                    MOLECULE **Molecule, BEADTYPE **BeadType, BEAD **Bead,
                    int **Index); //}}}
void FreeSystemInfo2(COUNTS Counts, MOLECULETYPE **MoleculeType,
                    MOLECULE **Molecule, int **Index_mol,
                    BEADTYPE **BeadType, BEAD **Bead, int **Index);
// PrintCounts_old()  //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts   numbers of beads, molecules, etc.
 */
void PrintCounts_old(COUNTS Counts);
//}}}
// PrintBeadType_old() //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts     numbers of beads, molecules, etc.
 * \param [in] BeadType   information about bead types
 */
void PrintBeadType_old(COUNTS Counts, BEADTYPE *BeadType); //}}}
void PrintBeadType2(int number, BEADTYPE *BeadType);
// PrintMoleculeType_old()  //{{{
/**
 * \brief Function printing MoleculeType structure.
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] MoleculeType  information about molecule types
 */
void PrintMoleculeType_old(COUNTS Counts, BEADTYPE *BeadType,
                       MOLECULETYPE *MoleculeType); //}}}
void PrintMoleculeType2(int number_of_types, BEADTYPE *BeadType,
                        MOLECULETYPE *MoleculeType);
// PrintMolecule_old() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule_old(int number_of_molecules, MOLECULETYPE *MoleculeType,
                   MOLECULE *Molecule, BEADTYPE *BeadType, BEAD *Bead); //}}}
// PrintBead_old() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead_old(COUNTS Counts, int *Index, BEADTYPE *BeadType, BEAD *Bead); //}}}
void PrintBead2(int number_of_beads, int *Index, BEADTYPE *BeadType,
                BEAD *Bead);
// FreeBead() //{{{
/**
 * \brief Free memory allocated for Bead struct array.
 *
 * \param [in]  number_of_beads   number of beads
 * \param [out] Bead              information about individual beads
 */
void FreeBead(int number_of_beads, BEAD **Bead); //}}}
// FreeMolecule() //{{{
/**
 * \brief Free memory allocated for Molecule struct array.
 *
 * \param [in]  number_of_molecules   number of molecules
 * \param [out] Molecule              information about individual molecules
 */
void FreeMolecule(int number_of_molecules, MOLECULE **Molecule); //}}}
// FreeMoleculeType_old() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  number_of_types  number of molecule types
 * \param [out] MoleculeType     information about molecule types
 */
void FreeMoleculeType_old(int number_of_types, MOLECULETYPE **MoleculeType); //}}}
// FindBeadType_old() //{{{
/** \brief Function to identify type of bead from its name
 *
 * \param [in]  name      bead name
 * \param [in]  Counts    numbers of beads, residues, etc.
 * \param [in]  BeadType  information about bead types
 * \return bead type id corresponding to index in BeadType struct (or -1 if
 * non-existent bead name)
 */
int FindBeadType_old(char *name, COUNTS Counts, BEADTYPE *BeadType); //}}}
int FindBeadType2(char *name, int types_of_beads, BEADTYPE *BeadType);
// FindMoleculeType_old() //{{{
/** \brief Function to identify type of molecule from its name
 *
 * \param [in]  name          molecule name
 * \param [in]  Counts        numbers of beads, residues, etc.
 * \param [in]  MoleculeType  information about bead types
 * \return molecule type      id corresponding to index in BeadType struct (or
 * -1 for non-existent molecule)
 */
int FindMoleculeType_old(char *name, COUNTS Counts,
                     MOLECULETYPE *MoleculeType); //}}}
int FindMoleculeType2(char *name, int number_of_types,
                      MOLECULETYPE *MoleculeType);
// RemovePBCMolecules_old() //{{{
/**
 * \brief Function to join all molecules.
 *
 * \param [in]  Counts         numbers of beads, molecules, etc.
 * \param [in]  BoxLength      dimension of the simulation box
 * \param [in]  BeadType       information about bead types
 * \param [out] Bead           information about individual beads (coordinates)
 * \param [in]  MoleculeType   information about molecule types
 * \param [in]  Molecule       information about individual molecules
 */
void RemovePBCMolecules_old(COUNTS Counts, VECTOR BoxLength, BEADTYPE *BeadType,
                            BEAD **Bead, MOLECULETYPE *MoleculeType,
                            MOLECULE *Molecule); //}}}
void RemovePBCMolecules2(COUNTS Counts, BOX Box, BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule);
void RemovePBCMolecules_new(COUNTS Counts, BOX Box, BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule);
void CopyBead_old(int number_of_beads, BEAD **b_out, BEAD *b_in, int mode);
void FreeSystem_old(SYSTEM *System);
void FillMolBTypes_old(int number_of_types, MOLECULETYPE **MoleculeType);
void FillMolMassCharge_old(int number_of_types, MOLECULETYPE **MoleculeType,
                       BEADTYPE *BeadType);
// FreeAggregate() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Aggregate   information about individual molecules
 */
void FreeAggregate(COUNTS Counts, AGGREGATE **Aggregate); //}}}
SYSTEM CopySystem_old(SYSTEM Sys_in);
void PruneSystem_old(SYSTEM *System);
// CopyBeadType() //{{{
/**
 * Function to copy BEADTYPE structure into a new one.
 */
void CopyBeadType(int number_of_types, BEADTYPE **bt_out, BEADTYPE *bt_in,
                  int mode); //}}}
// CopyMoleculeType() //{{{
/**
 * Function to copy MOLECULETYPE structure into a new one.
 */
void CopyMoleculeType_old(int number_of_types, MOLECULETYPE **mt_out,
                      MOLECULETYPE *mt_in, int mode); //}}}
// CopyMolecule() //{{{
/*
 * Function to copy a MOLECULE struct into a new one.
 */
void CopyMolecule(int number_of_molecules, MOLECULETYPE *mt, MOLECULE **m_out,
                  MOLECULE *m_in, int mode); //}}}
// RemovePBCAggregates() //{{{
/**
 * \brief Funcion to join all aggregates.
 *
 * \param [in]  distance       distance for closeness check (taken from agg
 * file) \param [in]  Aggregate      information about aggregates \param [in]
 * Counts         number of beads, molecu.es, etc. \param [in]  BoxLength
 * dimensions of the simulation box \param [in]  BeadType       information
 * about bead types \param [out] Bead           information about individual
 * beads (coordinates) \param [in]  MoleculeType   information about molecule
 * types \param [in]  Molecule       information about individual molecules
 */
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNTS Counts,
                         VECTOR BoxLength, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}
#endif //}}}
#endif
