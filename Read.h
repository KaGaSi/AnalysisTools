/**
 * \file
 * \brief Functions reading files.
 */

#ifndef _READ_H_
#define _READ_H_

#include "AnalysisTools.h"

#define ERROR_LINE -1
#define BLANK_LINE 0
#define COMMENT_LINE 1
#define PBC_LINE 2
#define PBC_LINE_ANGLES 3
#define COOR_LINE_I 4
#define COOR_LINE_O 5
#define ATOM_LINE 6
#define BOND_LINE 7
#define TIME_LINE_I 8
#define TIME_LINE_O 9

// Functions to read vtf files //{{{

// Get the first pbc line from a vcf/vtf coordinate file.
void VtfReadPBC(char input_vcf[], BOX *Box);
// Read all information about the system from a vsf/vtf structure file.
void VtfReadStruct(char vsf_file[], bool detailed, COUNTS *Counts,
                   BEADTYPE *BeadType[], BEAD *Bead[], int *Index[],
                   MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[],
                   int *Index_mol[]);
// Read a single timestep from a vcf/vtf coordinate file
bool VtfReadTimestep(FILE *vcf, char vcf_file[], BOX *Box, COUNTS *Counts,
                     BEADTYPE BeadType[], BEAD *Bead[], int Index[],
                     MOLECULETYPE MoleculeType[], MOLECULE Molecule[],
                     int *InVcfFile[], int *file_line_count,
                     int step_count, char stuff[]);
// Discard a single timestep from a vcf/vtf coordinate file
bool VtfSkipTimestep(FILE *vcf, char vcf_file[],
                     int *file_line_count, int step_count);
bool VtfSkipCoorOrderedLine(FILE *fr);
// Find position of atom line keywords in the provided strtok'd line
int * VtfAtomLineValues(int words, char *split[]);
// functions checkingi validity of line types
int VtfCheckLineType(int words, char *split[], char *file, int line);
int VtfCheckCoorOrderedLine(int words, char *split[]);
int VtfCheckCoorIndexedLine(int words, char *split[]);
int VtfCheckCoordinateLine(int words, char *split[]);
int VtfCheckTimestepLine(int words, char *split[]);
int VtfCheckPbcLine(int words, char *split[]);
bool VtfCheckAtomLine(int words, char *split[]);
bool VtfCheckBondLine(int words, char *split[]);
 //}}}

// TODO will be changed - agg files
// ReadAggCommand() //{{{
/*
 * \brief Function reading Aggregate command from agg file.
 *
 * \param [in]  BeadType      information about bead types
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  input_coor    coordinate file
 * \param [in]  input_agg     aggregate file
 * \param [out] distance      <distance> parameter from Aggregate command
 * \param [out] contacts      <contacts> parameter from Aggregate command
 */
void ReadAggCommand(BEADTYPE *BeadType, COUNTS Counts,
                    char *input_coor, char *input_agg,
                    double *distance, int *contacts); //}}}
// SkipAgg() //{{{
/*
 * \brief Function to skip one timestep in coordinates file.
 *
 * \param [in] agg        pointer to the open agg file
 * \param [in] agg_file   agg file name
 */
void SkipAgg(FILE *agg, char *agg_file); //}}}
// ReadAggregates() //{{{
/*
 * \brief Function reading information about aggregates from `.agg` file
 *
 * \param [in]  fr            pointer to open aggregate file
 * \param [in]  agg_file      name of aggregate file
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  BeadType      information about bead types
 * \param [out] Bead          information about individual beads
 * \param [out] Aggregate     information about aggregates
 * \param [in]  MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 */
void ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts, AGGREGATE *Aggregate[],
                    BEADTYPE *BeadType, BEAD *Bead[],
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule[], int *Index); //}}}

// TODO will be changed - FIELD file
bool ReadFieldPbc(char *field, VECTOR *BoxLength);
void ReadFieldBeadType(char *field, COUNTS *Counts,
                       BEADTYPE *BeadType[], BEAD *Bead[]);
// ReadFieldMolecules() //{{{
void ReadFieldMolecules(char *field, COUNTS *Counts,
                        BEADTYPE *BeadType[], BEAD *Bead[],
                        MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[],
                        PARAMS *bond_type[], PARAMS *angle_type[],
                        PARAMS *dihedral_type[]); //}}}
// ReadField() //{{{
/*
 * \brief Function reading structure information from FIELD-like file
 *
 * \param [in]  field         input FIELD-like file
 * \param [out] BoxLength     simulation box size
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * \param [out] bond_types    information abount bond types
 * \param [out] angle_types   information abount angle types
 * */
void ReadField(char *field, VECTOR *BoxLength, COUNTS *Counts,
               BEADTYPE *BeadType[], BEAD *Bead[], int *Index[],
               MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[],
               PARAMS *bond_type[], PARAMS *angle_type[],
               PARAMS *dihedral_type[]); //}}}

// TODO will be changed - lammps data file
// ReadLmpData() //{{{
/*
 * Function reading all information from lammps data file
 *
 *  [in] data_field    input data file file
 * [out] bonds         number of bonds
 * [out] bond_type     information about bond types
 * [out] angles        number of angles
 * [out] angle_type    information about angle types
 * [out] BoxLength     simulation box size
 * [out] box_lo        minimum box coordinates
 * [out] Counts        numbers of beads, molecules, etc.
 * [out] BeadType      information about bead types
 * [out] Bead          informationn about individual beads
 * [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * [out] MoleculeType  information about molecule types
 * [out] Molecule      information about individual molecules
 * */
void ReadLmpData(char *data_file, int *bonds, PARAMS *bond_type[],
                 int *angles, PARAMS *angle_type[],
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE *BeadType[], BEAD *Bead[], int *Index[],
                 MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[]); //}}}

bool VtfSkipCoorOrderedLine(FILE *fr);

// helper functions //{{{
// Fill arrays of structures (BEADTYPE & MOLECULTYPE)
// NewBeadType() //{{{
/*
 * Function to add a new bead type to a BEADTYPE struct
 * (and increment the number of bead types).
 */
void NewBeadType(BEADTYPE *BeadType[], int *number_of_types, char *name,
                 double charge, double mass, double radius); //}}}
// NewMolType() //{{{
/*
 * Function to create a new molecule type in a MOLECULETYPE struct.
 */
void NewMolType(MOLECULETYPE *MoleculeType[], int *n_types, char *name,
                int n_beads, int n_bonds, int n_angles, int n_dihedrals); //}}}
// FillMolMass //{{{
/*
 * Function to calculate mass of all molecules. If at least one bead has
 * undefined mass, the mass of the molecule is also undefined.
 */
void FillMolMass(int number_of_types,
                 BEADTYPE *BeadType, MOLECULETYPE *MoleculeType[]); //}}}
// FillMolCharge //{{{
/*
 * Function to calculate charge of all molecules. If at least one bead has
 * undefined charge, the charge of the molecule is also undefined.
 */
void FillMolCharge(int number_of_types, BEADTYPE *BeadType,
                   MOLECULETYPE *MoleculeType[]); //}}}
// FillMolType //{{{
/*
 * Function to fill BType array and mass and charge for each molecule type.
 */
void FillMolType(int number_of_types, BEADTYPE *BeadType,
                 MOLECULETYPE *MoleculeType[]); //}}}
 //}}}
#endif
