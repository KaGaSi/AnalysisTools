#ifndef _READ_H_
#define _READ_H_

#include "AnalysisTools.h"

SYSTEM ReadStructure(int struct_type, char struct_file[],
                     int coor_type, char coor_file[], bool detailed,
                     bool vtf_coor_var, int pbc_xyz);
bool ReadTimestep(int coor_type, FILE *fr, char file[], SYSTEM *System,
                  int *line_count, bool vtf_var_coor);
bool SkipTimestep(int coor_type, FILE *f, char file1[], char file2[],
                  int *line_count);

#if 0 //{{{
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
void ReadAggCommand(BEADTYPE *BeadType, COUNT Counts,
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
void ReadAggregates(FILE *fr, char *agg_file, COUNT *Counts, AGGREGATE *Aggregate[],
                    BEADTYPE *BeadType, BEAD *Bead[],
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule[], int *Index); //}}}
#endif //}}}
#endif
