#ifndef READ_H
#define READ_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

// SYSTEM ReadStructure(int struct_type, char struct_file[], int coor_type,
//                      char coor_file[], bool detailed);
SYSTEM ReadStructure(SYS_FILES f, bool detailed);
bool ReadTimestep(SYS_FILES f, FILE *fr, SYSTEM *System, int *line_count);
bool SkipTimestep(SYS_FILES f, FILE *fr, int *line_count);
int ReadAggregates(FILE *fr, char file[], SYSTEM *System,
                   AGGREGATE Aggregate[], int *line_count);
#endif
