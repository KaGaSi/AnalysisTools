#ifndef READWRITEFIELD_H
#define READWRITEFIELD_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

SYSTEM FieldRead(const char *file);
void WriteField(const SYSTEM System, const char *file_field,
                const int argc, char **argv);

#endif

