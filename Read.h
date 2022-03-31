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

// GetPBC_old() //{{{
/*
 * \brief Function to get box dimensions from a vtf coordinate file.
 *
 * \param [in] coor_file   name of the coordinate file
 * \return vector with box dimensions
 */
VECTOR GetPBC_old(char *coor_file); //}}}
// GetPBC() //{{{
/*
 * \brief Function to get box dimensions from a vtf coordinate file.
 *
 * \param [in]  coor_file   name of the coordinate file
 * \param [out] Box         box dimensions and angles
 */
void VtfGetPBC(char *coor_file, BOX *Box); //}}}

// ReadAggCommand() //{{{
/**
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

bool VtfCheckTimestep(FILE *vcf, char *vcf_file, COUNTS *Counts,
                      BEADTYPE **BeadType, BEAD **Bead, int **Index,
                      MOLECULETYPE **MoleculeType, MOLECULE **Molecule);
bool VtfReadTimestep(FILE *vcf, char *vcf_file, BOX *Box, COUNTS *Counts,
                     BEADTYPE *BeadType, BEAD **Bead, int *Index,
                     MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                     int *file_line_count, int step_count);
bool VtfReadTimestep2(FILE *vcf, char *vcf_file, BOX *Box, COUNTS *Counts,
                     BEADTYPE *BeadType, BEAD **Bead, int *Index,
                     MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                     int *file_line_count, int step_count);
bool VtfSkipTimestep(FILE *vcf, char *vcf_file,
                     int *file_line_count, int step_count);
bool VtfSkipTimestep2(FILE *vcf, char *vcf_file,
                     int *file_line_count, int step_count);
void VtfReadStruct_old(char *vsf_file, bool detailed, COUNTS *Counts,
                      BEADTYPE **BeadType, BEAD **Bead, int **Index,
                      MOLECULETYPE **MoleculeType, MOLECULE **Molecule);
void VtfReadStruct(char *vsf_file, bool detailed, COUNTS *Counts,
                   BEADTYPE **BeadType, BEAD **Bead, int **Index,
                   MOLECULETYPE **MoleculeType, MOLECULE **Molecule);
void FullVtfRead(char *struct_file, char *vcf_file, bool detailed, bool vtf,
                 bool *indexed, int *struct_lines, BOX *Box, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule);
void FullVtfRead_new(char *struct_file, bool detailed, bool *indexed,
                 COUNTS *Counts, BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule);

// ReadStructure() //{{{
/**
 * \brief Function reading information from dl_meso FIELD and vsf
 * structure files.
 *
 * \param [in]  vsf_file      .vsf structure file
 * \param [in]  vcf_file      .vcf coordinate file
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * \return 'true' or 'false' for .vcf file with indexed or ordered
 * timesteps, respectively
 * */
bool ReadStructure(char *vsf_file, char *vcf_file, COUNTS *Counts,
                   BEADTYPE **BeadType, BEAD **Bead, int **Index,
                   MOLECULETYPE **MoleculeType, MOLECULE **Molecule); //}}}

int ReadVtfTimestepPreamble_old(bool *indexed, char *input_coor, FILE *vcf_file,
                            char **stuff, VECTOR *BoxLength, bool quit);
int ReadVtfTimestepPreamble(bool *indexed, char *input_coor, FILE *vcf_file,
                            char **stuff, BOX *Box, bool quit);
bool LastStep(FILE *vcf_file, FILE *agg_file);

// ReadCoordinates() //{{{
/**
 * \brief Function reading ordered coordinates from .vcf coordinate file.
 *
 * \param [in]  indexed    is the vcf indexed?
 * \param [in]  input_coor name of input coordinate file
 * \param [in]  vcf_file   pointer to the open coordinate file
 * \param [in]  Counts     numbers of beads, molecules, etc.
 * \param [in]  Index      bead indices between program and vsf
 * \param [out] Bead       coordinates of individual beads
 * \param [out] stuff      first line of a timestep
 */
void ReadCoordinates_old(bool indexed, char *input_coor, FILE *vcf_file, COUNTS Counts, int *Index, BEAD **Bead, char **stuff); //}}}

// ReadVcfCoordinates() //{{{
/**
 * \brief Function reading ordered coordinates from .vcf coordinate file.
 *
 * \param [in]  indexed    is the vcf indexed?
 * \param [in]  input_coor name of input coordinate file
 * \param [in]  vcf_file   pointer to the open coordinate file
 * \param [in]  Counts     numbers of beads, molecules, etc.
 * \param [in]  Index      bead indices between program and vsf
 * \param [out] Bead       coordinates of individual beads
 * \param [out] stuff      first line of a timestep
 */
void ReadVcfCoordinates_old(bool indexed, char *input_coor, FILE *vcf_file,
                        VECTOR *BoxLength, COUNTS Counts,
                        int *Index, BEAD **Bead, char **stuff); //}}}
void ReadVcfCoordinates(bool indexed, char *input_coor, FILE *vcf_file,
                        BOX *Box, COUNTS Counts,
                        int *Index, BEAD **Bead, char **stuff);

// SkipCoor() //{{{
/**
 * \brief Function to skip one timestep in coordinates file.
 *
 * \param [in]  vcf_file   file with vcf coordinates
 * \param [in]  Counts     number of beads in vcf file
 * \param [out] stuff      first line of a timestep
 * \return 1 if premature end of file or 0 for no error
 */
bool SkipCoor(FILE *vcf_file, COUNTS Counts, char **stuff); //}}}

// SkipAgg() //{{{
/**
 * \brief Function to skip one timestep in coordinates file.
 *
 * \param [in] agg        pointer to the open agg file
 * \param [in] agg_file   agg file name
 */
void SkipAgg(FILE *agg, char *agg_file); //}}}

// SkipVcfCoor() //{{{
/**
 * \brief Function to skip one timestep in coordinates file.
 *
 * \param [in]  vcf_file   file with vcf coordinates
 * \param [in]  input_coor name  of vcf_file
 * \param [in]  Counts     number of beads in vcf file
 * \param [out] stuff      first line of a timestep
 */
void SkipVcfCoor(FILE *vcf_file, char *input_coor,
                 COUNTS Counts, char **stuff); //}}}

// ReadAggregates() //{{{
/**
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
void ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts, AGGREGATE **Aggregate,
                    BEADTYPE *BeadType, BEAD **Bead,
                    MOLECULETYPE *MoleculeType, MOLECULE **Molecule, int *Index); //}}}

// ReadField() //{{{
/**
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
               BEADTYPE **BeadType, BEAD **Bead, int **Index,
               MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
               PARAMS **bond_type, PARAMS **angle_type,
               PARAMS **dihedral_type); //}}}

// ReadLmpData() //{{{
/**
 * \brief Function reading all information from lammps data file
 *
 * \param [in]  data_field    input data file file
 * \param [out] bonds         number of bonds
 * \param [out] bond_type     information about bond types
 * \param [out] angles        number of angles
 * \param [out] angle_type    information about angle types
 * \param [out] BoxLength     simulation box size
 * \param [out] box_lo        minimum box coordinates
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * */
void ReadLmpData(char *data_file, int *bonds, PARAMS **bond_type,
                 int *angles, PARAMS **angle_type,
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule); //}}}

// SkipCoorSteps() { //{{{
int SkipCoorSteps(FILE *vcf, char *input_coor, COUNTS Counts, int start, bool silent); //}}}

// SkipCoorAggSteps() { //{{{
int SkipCoorAggSteps(FILE *vcf, char *input_coor, FILE *agg, char *input_agg, COUNTS Counts, int start, bool silent); //}}}

int VtfCheckTimestepLine(int words, char split[SPL_STR][SPL_LEN]);

int VtfCheckPbcLine(int words, char split[SPL_STR][SPL_LEN], char *file);

bool VtfCheckAtomLine(int words, char split[SPL_STR][SPL_LEN],
                      char *file, int file_line_count);

bool VtfCheckBondLine(int words, char split[SPL_STR][SPL_LEN]);

bool VtfSkipCoorOrderedLine(FILE *fr);
int VtfCheckCoorOrderedLine(int words, char *split[SPL_STR]);
int VtfCheckCoorIndexedLine(int words, char *split[SPL_STR]);
int VtfCheckCoordinateLine(int words, char *split[SPL_STR]);
int VtfCheckLineType(int words, char split[SPL_STR][SPL_LEN], bool indexed,
                     char *file, int line);

// NewBeadType() //{{{
/*
 * Function to add a new bead type to a BEADTYPE struct
 * (and increment the number of bead types).
 */
void NewBeadType(BEADTYPE **BeadType, int *number_of_types, char *name,
                 double charge, double mass, double radius); //}}}

// NewMolType() //{{{
/*
 * Function to create a new molecule type in a MOLECULETYPE struct.
 */
void NewMolType(MOLECULETYPE **MoleculeType, int *n_types, char *name,
                int n_beads, int n_bonds, int n_angles, int n_dihedrals); //}}}

// FillMolMass //{{{
/*
 * Function to calculate mass of all molecules. If at least one bead has
 * undefined mass, the mass of the molecule is also undefined.
 */
void FillMolMass(int number_of_types,
                 BEADTYPE *BeadType, MOLECULETYPE **MoleculeType); //}}}

// FillMolCharge //{{{
/*
 * Function to calculate charge of all molecules. If at least one bead has
 * undefined charge, the charge of the molecule is also undefined.
 */
void FillMolCharge(int number_of_types, BEADTYPE *BeadType,
                   MOLECULETYPE **MoleculeType); //}}}

// FillMolType //{{{
/*
 * Function to fill BType array and mass and charge for each molecule type.
 */
void FillMolType(int number_of_types, BEADTYPE *BeadType,
                 MOLECULETYPE **MoleculeType); //}}}

// TODO not used
int VtfCountStructLines(bool vtf, char *input);
// TODO not to be used
void SkipVtfStructure(FILE *vcf, int struct_lines);
bool CheckVtfTimestepLine_old(int words, char split[SPL_STR][SPL_LEN]);
bool CheckVtfAtomLine_old(int words, char split[SPL_STR][SPL_LEN], char *error);
bool CheckVtfBondLine_old(int words, char split[SPL_STR][SPL_LEN], char *error);
bool CheckVtfCoordinateLine_old(int words, char split[SPL_STR][SPL_LEN],
                                bool indexed);
#endif
