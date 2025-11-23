#ifndef READWRITEGROMACS_H
#define READWRITEGROMACS_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"
SYSTEM ItpReadStruct(const char *file);
SYSTEM PdbReadStruct(const char *file);

#endif
