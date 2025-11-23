#ifndef READWRITELDATA_H
#define READWRITELDATA_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

int LmpDataReadTimestep(FILE *fr, const char *file,
                        SYSTEM *System, int *line_count);
SYSTEM LmpDataReadStruct(const char *file);
void WriteLmpData(const SYSTEM System, const char *file, const bool mass,
                  const int argc, char **argv);

#endif
