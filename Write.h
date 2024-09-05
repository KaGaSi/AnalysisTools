#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

void InitCoorFile(FILE_TYPE fout, SYSTEM System, int argc, char *argv[]);
void WriteOutput(SYSTEM System, bool write[], FILE_TYPE fw,
                 bool lmp_mass, int vsf_def, int argc, char *argv[]);
void WriteTimestep(FILE_TYPE f, SYSTEM System, int count_step, bool write[],
                   int argc, char *argv[]);
void WriteStructure(FILE_TYPE f, SYSTEM System, int vsf_def_type,
                    bool lmp_mass, int argc, char *argv[]);
void WriteAggregates(int step_count, char *agg_file, SYSTEM System,
                     AGGREGATE *Aggregate);
#endif
