#ifndef AGGREGATES_H
#define AGGREGATES_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

typedef struct agg_choice{
  bool *x, *only, *m; // arrays for which molecules to use for -x, -only, and -m
  bool x_flag, only_flag, m_flag; // flags whethere -x, -only, and -m are used
  int range[2];       // -n; used with 1 & Count.Molecules if no -n
} AGG_CHOICE;

// determine what molecules belong to what aggregates
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System, OPT opt);
// based on options, should an aggregate be used for calculations?
bool UseAggregate(SYSTEM System, AGGREGATE *Aggregate, int id,
                  AGG_CHOICE agg, int *size, double *mass);
void AggChoiceOptions(const int argc, char **argv, AGG_CHOICE *opt,
                      SYSTEM System);
#endif
