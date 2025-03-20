#ifndef READWRITEVTF_H
#define READWRITEVTF_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

SYSTEM VtfReadStruct(const char *file, const bool detailed);
int VtfReadTimestep(FILE *fr, const char *file,
                    SYSTEM *System, int *line_count);
int VtfSkipTimestep(FILE *fr, const char *file,
                    const char *vsf_file, int *line_count);
// read the first timestep to count beads in constant-size coordinate block
int VtfReadNumberOfBeads(const char *file);
// Get the first pbc line from a vcf/vsf/vtf coordinate file.
BOX VtfReadPBC(const char *file);

void VtfWriteStruct(char *file, SYSTEM System, int type_def,
                    const int argc, char **argv);
void VtfWriteCoorIndexed(FILE *fw, const bool *write, const SYSTEM System);
#endif
