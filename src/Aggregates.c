#include "Aggregates.h"
#include "Errors.h"

// TODO: AggPickerOptions should be in Options.c, no?

// based on options, should an aggregate be used for calculations? //{{{
bool UseAggregate(SYSTEM System, AGGREGATE *Aggregate, int id,
                  AGG_PICKER agg, int *size, double *mass) {
  bool only = true;
  bool x = false;
  *size = 0;
  *mass = 0;
  for (int j = 0; j < Aggregate[id].nMolecules; j++) {
    MOLECULE *mol = &System.Molecule[Aggregate[id].Molecule[j]];
    int mtype = mol->Type;
    MOLECULETYPE *mt = &System.MoleculeType[mtype];
    if (agg.m[mtype]) {
      (*size)++;
      *mass += mt->Mass;
    }
    // if at least one unwanted molecule is present, don't use aggregate
    if (!agg.only[mtype]) {
      only = false;
      break;
    }
    // if at least one molecule isn't exluded, use aggregate
    if (!agg.x[mtype]) {
      x = true;
    }
  }
  if (*size == 0 || // no molecules remained in aggregate
      *size < agg.range[0] || *size > agg.range[1] || // -n: not in range
      !only || // -only: found molecule that weren't supposed to be in
      !x) { // -x: didn't find any un-excluded molecules
    return false;
  } else {
    return true;
  }
} //}}}
// detect -m, -x, -only, and -n options //{{{
void AggPickerOptions(const int argc, char **argv, AGG_PICKER *opt,
                      SYSTEM System) {
  COUNT *Count = &System.Count;
  // '-n' option
  opt->range[0] = 1;
  opt->range[1] = Count->Molecule;
  TwoNumbersOption(argc, argv, "-n", opt->range, 'i');
  if (opt->range[0] > opt->range[1]) {
    SwapInt(&opt->range[0], &opt->range[1]);
  }
  // '-m' option - define aggregate size as sum of those molecule types
  opt->m = calloc(Count->MoleculeType, sizeof *opt->m);
  opt->m_flag = true;
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->m, System)) {
    opt->m_flag = false;
    InitBoolArray(opt->m, Count->MoleculeType, true);
  }
  // '-only' option - use aggregates composed only of specified molecule types
  opt->only = calloc(Count->MoleculeType, sizeof *opt->only);
  opt->only_flag = true;
  if (!TypeOption(argc, argv, "-only", 'm', true, opt->only, System)) {
    opt->only_flag = false;
    InitBoolArray(opt->only, Count->MoleculeType, true);
  }
  // '-x' option - exclude aggregates composed only of specified molecule types
  opt->x = calloc(Count->MoleculeType, sizeof *opt->x);
  opt->x_flag = true;
  if (!TypeOption(argc, argv, "-x", 'm', true, opt->x, System)) {
    opt->x_flag = false;
  }
  // error checking //{{{
  // all molecule specified
  if (opt->x_flag) {
    bool overlap = true; // are all molecule types specified by -x?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (!opt->x[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("with all molecules listed, no aggregates would be detected");
      PrintErrorOption("-x");
      exit(1);
    }
  }
  // molecules specified by -m and -only do not overlap
  if (opt->m_flag && opt->only_flag) {
    bool overlap = false;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->m[i] && opt->only[i]) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      err_msg("for any aggregate to be used, at least one molecule "
              "must be specified in both options");
      PrintErrorOption("-m/-only");
      exit(1);
    }
  }
  // molecules specified by -only and -x must differ
  if (opt->only_flag && opt->x_flag) {
    bool overlap = true; // do the two array fully overlap?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->x[i] != opt->only[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("the lists of molecules must be different");
      PrintErrorOption("-x/-only");
      exit(1);
    }
  } //}}}
} //}}}
