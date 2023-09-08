#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

void WriteTimestep(int coor_type, char file[], SYSTEM System,
                   int count_step, bool write[]);
void WriteStructure(int struct_type, char file[], SYSTEM System,
                    int vsf_def_type, bool lmp_mass);
void WriteAggregates(int step_count, char *agg_file, SYSTEM System,
                     AGGREGATE *Aggregate);
#endif
