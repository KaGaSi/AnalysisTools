#ifndef READWRITEXYZ_H
#define READWRITEXYZ_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

int XyzReadTimestep(FILE *fr, const char *file,
                    SYSTEM *System, int *line_count);
SYSTEM XyzReadStruct(const char *file);
bool XyzSkipTimestep(FILE *fr, const char *file, int *line_count);

void XyzWriteCoor(FILE *fw, const bool *write, SYSTEM System);

#endif
