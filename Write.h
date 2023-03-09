#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

void WriteTimestep(int coor_type, char file[], SYSTEM System,
                   int count_coor, bool write[]);
void WriteStructure(int struct_type, char file[], SYSTEM System,
                    int vsf_def_type, bool lmp_mass);

#if 0 //{{{
// TODO will change
void WriteAggregates(int step_count, char *agg_file, COUNTS Counts,
                     MOLECULETYPE *MoleculeType, BEAD *Bead, AGGREGATE
                     *Aggregate);
#endif //}}}
#endif
