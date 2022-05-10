#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

// Append an indexed timestep to a vcf/vtf coordinate file
void VtfWriteCoorIndexed(FILE *vcf, char *stuff, int InFile[],
                         COUNTS Counts, BEAD *Bead, BOX Box);
// Append a timestep to an xyz file
void WriteCoorXYZ(FILE *xyz, COUNTS Counts, int InFile[],
                  BEADTYPE *BeadType, BEAD *Bead);
// Create a new vsf/vtf structure file
void WriteVsf(char *input_vsf, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
              MOLECULETYPE *MoleculeType, MOLECULE *Molecule, bool change);

// TODO will change
void WriteAggregates(int step_count, char *agg_file, COUNTS Counts,
                     MOLECULETYPE *MoleculeType, BEAD *Bead, AGGREGATE
                     *Aggregate);

// TODO will change
void WriteField(char *field, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                PARAMS *bond_type, PARAMS *angle_type, PARAMS *dihedral_type);

void PrintByline(FILE *ptr, int argc, char *argv[]);

// TODO remove
void WriteCoorIndexed(FILE *vcf_file, COUNTS Counts,
                      BEADTYPE *BeadType, BEAD *Bead,
                      MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                      char *stuff, BOX Box);
#endif
