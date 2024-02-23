#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

void WriteTimestep(FILE_TYPE f, SYSTEM System, int count_step, bool write[]);
void WriteStructure(FILE_TYPE f, SYSTEM System, int vsf_def_type,
                    bool lmp_mass, int argc, char *argv[]);
void WriteAggregates(int step_count, char *agg_file, SYSTEM System,
                     AGGREGATE *Aggregate);
#endif
