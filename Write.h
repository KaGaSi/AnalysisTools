/**
 * \file
 * \brief Functions writing to files.
 */

#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

// WriteCoorIndexed //{{{
/**
 * \brief Function writing indexed coordinates to a .vcf file.
 *
 * \param [in] vcf_file      name of output .vcf coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          coordinates of individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      coordinates of individual molecules
 * \param [in] stuff         array of chars containing comment line to place at the beginning
 */
void WriteCoorIndexed(FILE *vcf_file, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                      MOLECULETYPE *MoleculeType, MOLECULE *Molecule, char *stuff); //}}}

// WriteCoorXYZ() //{{{
/**
 * \brief Function for writing xyz coordinates
 *
 * \param [in] xyz_file      output .xyz coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          coordinates of individual beads
 */
void WriteCoorXYZ(FILE *xyz_file, COUNTS Counts,
                  BEADTYPE *BeadType, BEAD *Bead); //}}}

// WriteVsf() //{{{
/**
 * \brief Function writing vsf file
 *
 * \param [in] vsf_file      name of output .vsf structure file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          coordinates of individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      coordinates of individual molecules
 * \param [in] change        true or false if all molecules of a given type should contain the same beads as the first molecule of its kind in the vsf
 */
void WriteVsf(char *input_vsf, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
              MOLECULETYPE *MoleculeType, MOLECULE *Molecule, bool change); //}}}

// WriteAggregates() //{{{
/**
 * \brief Function writing agg file
 *
 * \param [in] step_count    current timestep
 * \param [in] agg_file      name of output .agg file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Bead          coordinates of individual beads
 * \param [in] Aggregates    information about aggregates
 */
void WriteAggregates(int step_count, char *agg_file, COUNTS Counts,
                     MOLECULETYPE *MoleculeType, BEAD *Bead, AGGREGATE *Aggregate); //}}}
#endif
