/**
 * \file
 * \brief Functions common to all analysis utilities.
 */

#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include "General.h"
#include "Errors.h"
#include "Structs.h"
#include "Options.h"
#include "Read.h"
#include "Write.h"

// TransformMatrices()
void TriclinicCellData(BOX *Box);

void ToFractional(VECTOR *coor, BOX Box);
void ToFractionalCoor(int number_of_beads, BEAD **Bead, BOX Box);
VECTOR FromFractional(VECTOR coor, BOX Box);
void FromFractionalCoor(int number_of_beads, BEAD **Bead, BOX Box);

// InputCoor() //{{{
/**
 * \brief Function test input coordinate file is correct
 *
 * \param [out] vtf          is the coordinate file vtf or vcf?
 * \param [in]  file_coor    name of coordinate file
 * \param [in]  file_struct  name of structure file
 * \return false if file_coor has wrong extension, true otherwise
 */
bool InputCoor(bool *vtf, char *file_coor, char *file_struct); //}}}

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
// VerboseOutput() //{{{
/**
 * \brief Function printing basic information about system if `-v` or `-V`
 * option is provided
 *
 * \param [in] input_vcf     .vcf coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] Box           dimension and angles of the simulation box
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          informationn about individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      information about individual molecules
 */
void VerboseOutput(char *input_vcf, COUNTS Counts, BOX Box,
                   BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}

// PrintCounts()  //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts   numbers of beads, molecules, etc.
 */
void PrintCounts(COUNTS Counts);
//}}}

// PrintBeadType() //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts     numbers of beads, molecules, etc.
 * \param [in] BeadType   information about bead types
 */
void PrintBeadType(COUNTS Counts, BEADTYPE *BeadType); //}}}
void PrintBeadType2(int number, BEADTYPE *BeadType);

// PrintMoleculeTypeType()  //{{{
/**
 * \brief Function printing MoleculeType structure.
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] MoleculeType  information about molecule types
 */
void PrintMoleculeType(COUNTS Counts, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType); //}}}
void PrintMoleculeType2(int number_of_types, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType);

// PrintBead() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead(COUNTS Counts, int *Index, BEADTYPE *BeadType, BEAD *Bead); //}}}
void PrintBead2(int number_of_beads, int *Index, BEADTYPE *BeadType, BEAD *Bead);

// PrintMolecule() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule(int number_of_molecules,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                   BEADTYPE *BeadType, BEAD *Bead); //}}}

// PrintAggregate() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintAggregate(COUNTS Counts, int *Index,
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                    BEAD *Bead, BEADTYPE *BeadType, AGGREGATE *Aggregate); //}}}

// PrintBondTypes() //{{{
void PrintBondTypes(COUNTS Counts, PARAMS *bond_type);  //}}}
void PrintBondTypes2(int number_of_bonds, PARAMS *bond_type);

// PrintAngleTypes() //{{{
void PrintAngleTypes(COUNTS Counts, PARAMS *angle_type);  //}}}
void PrintAngleTypes2(int number_of_angles, PARAMS *angle_type);

void PrintDihedralTypes2(int number_of_dihedrals, PARAMS *dihedral_type);

// FindBeadType() //{{{
/** \brief Function to identify type of bead from its name
 *
 * \param [in]  name      bead name
 * \param [in]  Counts    numbers of beads, residues, etc.
 * \param [in]  BeadType  information about bead types
 * \return bead type id corresponding to index in BeadType struct (or -1 if non-existent bead name)
 */
int FindBeadType(char *name, COUNTS Counts, BEADTYPE *BeadType); //}}}
int FindBeadType2(char *name, int types_of_beads, BEADTYPE *BeadType);

// FindMoleculeType() //{{{
/** \brief Function to identify type of molecule from its name
 *
 * \param [in]  name          molecule name
 * \param [in]  Counts        numbers of beads, residues, etc.
 * \param [in]  MoleculeType  information about bead types
 * \return molecule type      id corresponding to index in BeadType struct (or -1 for non-existent molecule)
 */
int FindMoleculeType(char *name, COUNTS Counts, MOLECULETYPE *MoleculeType); //}}}
int FindMoleculeType2(char *name, int number_of_types, MOLECULETYPE *MoleculeType);

void FillMolBTypes(int number_of_types, MOLECULETYPE **MoleculeType);

// FillMolMassCharge() //{{{
/*
 * Function to calculate total mass and charge of molecules.
 */
void FillMolMassCharge(int number_of_types, MOLECULETYPE **MoleculeType,
                 BEADTYPE *BeadType); //}}}

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

// RemovePBCMolecules() //{{{
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
void RemovePBCMolecules_old(COUNTS Counts, VECTOR BoxLength,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}
// TODO: somehow generalise triclinic stuff
void RemovePBCMolecules(COUNTS Counts, BOX Box,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule);

// RemovePBCAggregates() //{{{
/**
 * \brief Funcion to join all aggregates.
 *
 * \param [in]  distance       distance for closeness check (taken from agg file)
 * \param [in]  Aggregate      information about aggregates
 * \param [in]  Counts         number of beads, molecu.es, etc.
 * \param [in]  BoxLength      dimensions of the simulation box
 * \param [in]  BeadType       information about bead types
 * \param [out] Bead           information about individual beads (coordinates)
 * \param [in]  MoleculeType   information about molecule types
 * \param [in]  Molecule       information about individual molecules
 */
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNTS Counts,
                         VECTOR BoxLength, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}

// RestorePBC_old() //{{{
/**
 * \brief Function to restore pbc.
 *
 * \param [in]  Counts         numbers of beads, molecules, etc.
 * \param [in]  BoxLength      dimension of the simulation box
 * \param [out] Bead           information about individual beads (coordinates)
 */
void RestorePBC_old(COUNTS Counts, VECTOR BoxLength, BEAD **Bead); //}}}
// RestorePBC() //{{{
/**
 * \brief Function to restore pbc.
 *
 * \param [in]  number_of_beads   number of beads in the system
 * \param [in]  BoxLength         dimension of the simulation box
 * \param [out] Bead              information about individual beads (coordinates)
 */
void RestorePBC(int number_of_beads, BOX Box, BEAD **Bead); //}}}

// CentreOfMass() //{{{
/**
 * \brief Function to calculate centre of mass for a collection of beads.
 *
 * \param [in] n              number of beads
 * \param [in] list           list of bead ids (corresponding to indices in Bead struct)
 * \param [in] Bead           information about individual beads (coordinates)
 * \param [in] BeadType       information about beadtypes (masses)
 * \return coordinates of centre of mass of a given aggregate
 */
VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType); //}}}

// GeomCentre() //{{{
/**
 * \brief Function to calculate geometric centre for a collection of beads.
 *
 * \param [in] n              number of beads
 * \param [in] list           list of bead ids (corresponding to indices in Bead struct)
 * \param [in] Bead           information about individual beads (coordinates)
 * \return coordinates of geometric centre of a given aggregate
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
VECTOR Gyration(int n, int *list, COUNTS Counts,
                BEADTYPE *BeadType, BEAD **Bead); //}}}

// EvaluateContacts() //{{{
/**
 * \brief Function evaluating contacts for aggregate detection
 *
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [out] Aggregate     information about aggregates
 * \param [in]  Molecule      information about individual molecules
 * \param [in]  contacts      number of contacts for aggregate check
 * \param [in]  contact       2D array containing number of contacts between molecules
 */
void EvaluateContacts(COUNTS *Counts, AGGREGATE **Aggregate,
                      MOLECULE **Molecule,
                      int contacts, int **contact);
//}}}

// SortAggStruct() //{{{
/**
 * \brief Function to sort Aggregate struct.
 *
 * \param [out] Aggregate  Aggregate struct to sort
 * \param [in]  Counts     numbers of beads, molecules, etc.
 */
void SortAggStruct(AGGREGATE **Aggregate, COUNTS Counts,
                   MOLECULE *Molecule, MOLECULETYPE *MoleculeType,
                   BEAD **Bead, BEADTYPE *BeadType); //}}}

// LinkedList() //{{{
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead,
                int **Head, int **Link, double cell_size, INTVECTOR *n_cells,
                int *Dcx, int *Dcy, int *Dcz); //}}}

// SortBonds() //{{{
/**
 * \brief Function to sort 2D array of bonds.
 *
 * \param [out] bond    2D array of bonds
 * \param [in]  length  number of bonds
 */
void SortBonds(int (*bond)[3], int length); //}}}

// SortAngles() //{{{
/**
 * \brief Function to sort 2D array of bonds.
 *
 * \param [out] angle   2D array of angles
 * \param [in]  length  number of angles
 */
void SortAngles(int (*angle)[4], int length); //}}}
void SortDihedrals(int (*dihedral)[5], int length);

// CopyBead() //{{{
/**
 * Function to copy BEAD structure into a new one.
 */
void CopyBead(int number_of_beads, BEAD **b_out, BEAD *b_in, int mode); //}}}

// CopyBeadType() //{{{
/**
 * Function to copy BEADTYPE structure into a new one.
 */
void CopyBeadType(int number_of_types, BEADTYPE **bt_out,
                  BEADTYPE *bt_in, int mode); //}}}

// CopyMoleculeType() //{{{
/**
 * Function to copy MOLECULETYPE structure into a new one.
 */
void CopyMoleculeType(int number_of_types, MOLECULETYPE **mt_out,
                      MOLECULETYPE *mt_in, int mode); //}}}

// CopyMolecule() //{{{
/*
 * Function to copy a MOLECULE struct into a new one.
 */
void CopyMolecule(int number_of_molecules, MOLECULETYPE *mt,
                  MOLECULE **m_out, MOLECULE *m_in, int mode); //}}}

// CopySystem() //{{{
/*
 * Function to copy the whole system - COUNTS, BEADTYPE, BEAD, MOLECULETYPE,
 * and MOLECULE structures and Index array.
 */
void CopySystem(COUNTS *Counts_out, COUNTS Counts_in,
                BEADTYPE **bt_out, BEADTYPE *bt_in,
                BEAD **bead_out, BEAD *bead_in, int **index_out, int *index_in,
                MOLECULETYPE **mt_out, MOLECULETYPE *m_in,
                MOLECULE **mol_out, MOLECULE *mol_in, int mode); //}}}

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

// FreeMoleculeType() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  number_of_types  number of molecule types
 * \param [out] MoleculeType     information about molecule types
 */
void FreeMoleculeType(int number_of_types, MOLECULETYPE **MoleculeType); //}}}

// FreeAggregate() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Aggregate   information about individual molecules
 */
void FreeAggregate(COUNTS Counts, AGGREGATE **Aggregate); //}}}

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
void FreeSystemInfo(COUNTS Counts, MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
                    BEADTYPE **BeadType, BEAD **Bead, int **Index); //}}}
#endif
