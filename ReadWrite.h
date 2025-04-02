#ifndef READWRITE_H
#define READWRITE_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

// General helper functions
void InitOutputCoorFile(const FILE_TYPE fout, const SYSTEM System,
                        const int argc, char **argv);
void CopyMoleculeTypeBeadsToMoleculeBeads(SYSTEM *System);
void FillAllMTypeStuff(SYSTEM *System, const int (*bond)[5], const int
                       (*angle)[5], const int (*dih)[5], const int (*imp)[5]);
void MinimizeMTypeStuffIds(SYSTEM *System);
void RemoveExtraTypes(SYSTEM *System);
void WriteBoxLengthAngles(FILE *fw, const BOX box);

// Reading functions
SYSTEM ReadStructure(const SYS_FILES f, const bool detailed);
bool ReadTimestep(const SYS_FILES f, FILE *fr,
                  SYSTEM *System, int *line_count);
bool SkipTimestep(const SYS_FILES f, FILE *fr, int *line_count);
int ReadAggregates(FILE *fr, const char *file, SYSTEM *System,
                   AGGREGATE *Aggregate, int *line_count);
bool SkipAggregates(FILE *fr, const char *file, int *line_count);

// writing functions
void WriteOutput(const SYSTEM System, const bool *write, FILE_TYPE fw,
                 const bool lmp_mass, const int vsf_def,
                 const int argc, char **argv);
void WriteOutputAll(const SYSTEM System, FILE_TYPE fw, const bool lmp_mass,
                    const int vsf_def, const int argc, char **argv);
void WriteTimestep(const FILE_TYPE f, const SYSTEM System, const int count_step,
                   const bool *write, const int argc, char **argv);
void WriteTimestepAll(FILE_TYPE f, SYSTEM System, int count_step,
                      int argc, char **argv);
void WriteStructure(FILE_TYPE f, const SYSTEM System, const int vsf_def_type,
                    const bool lmp_mass, const int argc, char **argv);
void WriteAggregates(const int step_count, const char *agg_file,
                     const SYSTEM System, const AGGREGATE *Aggregate);

// verbose output (print various structures and some such)
void VerboseOutput(const SYSTEM System);
void PrintCount(const COUNT Count);
void PrintBeadType(const SYSTEM System);
void PrintOneMolType(const SYSTEM System, const int n);
void PrintAllMolTypes(const SYSTEM System);
void Print1Molecule(const SYSTEM System, const int n);
void PrintMolecules(const SYSTEM System);
void PrintBead(const SYSTEM System);
void PrintBeadCoor(const SYSTEM System);
void PrintBondType(const SYSTEM System);
void PrintAngleType(const SYSTEM System);
void PrintDihedralType(const SYSTEM System);
void PrintImproperType(const SYSTEM System);
// TODO: use SYSTEM
void PrintBondTypes(const COUNT Counts, const PARAMS *bond_type);
// TODO: use SYSTEM
void PrintAngleTypes(const COUNT Counts, const PARAMS *angle_type);
void PrintBox(const BOX Box);
void PrintByline(const char *file, const int argc, char **argv);
void PrintStep(int *count_coor, const int start, const bool silent);
void PrintLastStep(const int coor, const int used, const bool silent);
void PrintAggregate(const SYSTEM System, const AGGREGATE *Aggregate);

void WriteFormatedDataLine(FILE *fw, const int columns, const double *data,
                           const int (*digits)[2]);
#endif
