#ifndef DEBUG_H
#define DEBUG_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"

// verbose output of the full System structure
void VerboseOutput(const SYSTEM System);
// print individual sub-structures
void PrintCount(const COUNT Count);
void PrintBeadType(const SYSTEM System);
void PrintOneMolType(const SYSTEM System, const int n);
void PrintAllMolTypes(const SYSTEM System);
void Print1Molecule(const SYSTEM System, const int n);
void PrintMolecules(const SYSTEM System);
void PrintBead(const SYSTEM System);
void PrintBeadCoor(const SYSTEM System);
void PrintBondType(const SYSTEM System);
void PrintAngleType(const SYSTEM System);
void PrintDihedralType(const SYSTEM System);
void PrintImproperType(const SYSTEM System);
// TODO: use SYSTEM
void PrintBondTypes(const COUNT Counts, const PARAMS *bond_type);
// TODO: use SYSTEM
void PrintAngleTypes(const COUNT Counts, const PARAMS *angle_type);
void PrintBox(const BOX Box);
// progress reporting
void PrintStep(int *count_coor, const int start, const bool silent);
void PrintLastStep(const int coor, const int used, const bool silent);
// aggregate debug output
void PrintAggregate(const SYSTEM System, const AGGREGATE *Aggregate);

#endif
