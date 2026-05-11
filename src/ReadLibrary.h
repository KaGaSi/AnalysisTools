#ifndef READ_LIBRARY_H
#define READ_LIBRARY_H

#include "AnalysisTools.h"

#define LIB_MAX_IDS 64
#define LIB_MAX_INTER 512

typedef struct {
  char id[16];
  int index; // index into System.BondType[]
} LIB_BOND_ID;

typedef struct {
  char id[16];
  int index; // index into System.AngleType[]
} LIB_ANGLE_ID;

typedef struct {
  char name1[BEAD_NAME];
  char name2[BEAD_NAME];
  double A, Rc, gamma;
} LIB_INTERACTION;

typedef struct {
  SYSTEM System;
  LIB_BOND_ID bond_id[LIB_MAX_IDS];
  int n_bond_ids;
  LIB_ANGLE_ID angle_id[LIB_MAX_IDS];
  int n_angle_ids;
  LIB_INTERACTION *inter;
  int n_inter;
} LIBRARY;

typedef struct {
  int n_beads;         // total bead count including counterion (-1 on error)
  char cion[MOL_NAME]; // counterion molecule name, empty string if none
} LIB_MOL_INFO;

// Read library metadata: all bead types, bond/angle types, and interactions
LIBRARY ReadLibrary(const char *lib_dir);

// Fill pot[i][j][{A,Rc,gamma}] from library interactions using bead type names
void FillPotFromLibrary(const LIBRARY *lib, const SYSTEM *System, ArrNDd *pot);

// Add n_mols copies of mol_name to lib->System
void ReadLibraryMolecule(const char *lib_dir, const char *mol_name,
                         int n_mols, LIBRARY *lib);

// Return n_beads and counterion name for mol_name from list_molecules.txt.
// n_beads includes counterion bead(s); cion is empty if has_counterion==0.
LIB_MOL_INFO LibraryMoleculeInfo(const char *lib_dir, const char *mol_name);

// Rename bead types in sys to library names using molecule type name matching.
// For each molecule type in sys whose name matches a library molecule file,
// each bead type used by that molecule is renamed to the library's bead type.
void RenameBeadTypesFromLibrary(SYSTEM *sys, const LIBRARY *lib,
                                const char *lib_dir);

void FreeLibrary(LIBRARY *lib);

#endif
