#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

void WriteTimestep(int coor_type, char file[], SYSTEM System,
                   int count_coor, char stuff[], bool write[]);

void VtfWriteCoorIndexed(FILE *vcf, char stuff[], bool write[], SYSTEM System);
void XyzWriteCoor(FILE *xyz, bool write[], char *stuff, SYSTEM System);
void VtfWriteStruct(char file[], SYSTEM System, int type_def);
void WriteLmpData_old(SYSTEM System, char file_lmp[], bool srp, bool mass);
void WriteLmpData(SYSTEM System, char file_lmp[], bool srp, bool mass);
void LtrjWriteCoor(FILE *vcf, int step, bool write[], SYSTEM System);
void WriteField(SYSTEM System, char file_field[]);
void WriteConfig(SYSTEM System, char file_config[]);

// TODO will change
void WriteAggregates(int step_count, char *agg_file, COUNTS Counts,
                     MOLECULETYPE *MoleculeType, BEAD *Bead, AGGREGATE
                     *Aggregate);

#if 0
// TODO will change
void WriteField(char *field, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                PARAMS *bond_type, PARAMS *angle_type, PARAMS *dihedral_type);
#endif

// TODO remove
void WriteCoorIndexed(FILE *vcf_file, COUNTS Counts,
                      BEADTYPE *BeadType, BEAD *Bead,
                      MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                      char *stuff, BOX Box);
void WriteVsf_old(char *input_vsf, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
              MOLECULETYPE *MoleculeType, MOLECULE *Molecule, bool change);
void VtfWriteStruct_old(char file[], COUNTS Counts,
                    BEADTYPE BeadType[], BEAD Bead[],
                    MOLECULETYPE MoleculeType[], MOLECULE Molecule[]);
void VtfWriteCoorIndexed_old(FILE *vcf, char *stuff, int InFile[],
                         COUNTS Counts, BEAD *Bead, BOX Box);
void XyzWriteCoor_old(FILE *xyz, COUNTS Counts, int InFile[],
                  BEADTYPE *BeadType, BEAD *Bead);
#endif
