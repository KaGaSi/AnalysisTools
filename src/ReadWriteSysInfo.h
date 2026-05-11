#ifndef READ_WRITE_SYS_INFO_H
#define READ_WRITE_SYS_INFO_H

#include "AnalysisTools.h"

// write system composition: "mol_name  n_molecules" line per molecule type
void WriteSysInfo(const char *filename, const SYSTEM *System);

// read system composition: assign MoleculeType[i].Name positionally from file.
// TODO: for now, ignores 1-word lines and doesn't use the molecule count
void ReadSysInfo(const char *filename, SYSTEM *System);

#endif
