#ifndef READWRITELTRJ_H
#define READWRITELTRJ_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

int LtrjReadTimestep(FILE *fr, const char *file, SYSTEM *System,
                     int *line_count);
int LtrjSkipTimestep(FILE *fr, const char *file, int *line_count);
SYSTEM LtrjReadStruct(const char *file);
BOX LtrjReadPBC(const char *file);
void LtrjWriteCoor(FILE *fw, const int step,
                   const bool *write, const SYSTEM System);

#endif
