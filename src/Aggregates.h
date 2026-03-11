#ifndef AGGREGATES_H
#define AGGREGATES_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"
#include "Arrays.h"
#include "PairHash.h"

// determine what molecules belong to what aggregates
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System,
                         OPT opt, bool *use_bt, ArrNDb *use_bt_pair);
// based on options, should an aggregate be used for calculations?
bool UseAggregate(SYSTEM System, AGGREGATE *Aggregate, int id,
                  AGG_PICKER agg, int *size, double *mass);
// detect -m, -x, -only, and -n options
void AggPickerOptions(const int argc, char **argv, AGG_PICKER *opt,
                      SYSTEM System);
// evaluate bead contacts to assign molecules to aggregates
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      const int contacts, PairHash *contact);
// remove PBC for aggregate molecules
void RemovePBCAggregates(const double distance, const AGGREGATE *Aggregate,
                         SYSTEM *System, const bool *use_bt);
#endif
