#ifndef READWRITE_H
#define READWRITE_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"

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

// file type detection
int FileTypeFromString(const char *str);
bool InputCoorStruct(const int argc, char **argv, SYS_FILES *f);
int StructureFileType(const char *name);
int CoordinateFileType(const char *name);
int FileType(const char *name);

// output file header comment
void PrintByline(const char *file, const int argc, char **argv);
FILE * PrintBylineOpenFile(const char *f, const int argc, char **argv);

void WriteFormatedDataLine(FILE *fw, const int columns, const double *data,
                           const int (*digits)[2]);
#endif
