#ifndef READ_LIBRARY_H
#define READ_LIBRARY_H

#include "AnalysisTools.h"
#include <dirent.h>

#define LIB_MAX_IDS 64
#define LIB_MAX_INTER 512
#define LIB_MAX_BEADS 64
#define LIB_MAX_TOPO 256
#define LIB_MAX_CION 4

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

// One 'counterion <name> [count]' line
typedef struct {
  char name[MOL_NAME];
  int count;
} LIB_CION;

/*
 * One molecule file. Beads, bonds and angles are the molecule's own: a
 * counterion is a separate molecule named by the cion[] list, never beads
 * appended here.
 */
typedef struct {
  char name[MOL_NAME];
  bool bilayer;   // role: bilayer (true) or soluble (false)
  double M_w;     // 0 when the file gives none
  bool has_M_w;
  LIB_CION cion[LIB_MAX_CION];
  int n_cion;

  char bead_name[LIB_MAX_BEADS][BEAD_NAME];
  vec3d bead_pos[LIB_MAX_BEADS];
  int n_beads;

  char bond_id[LIB_MAX_TOPO][16];
  int bond_i[LIB_MAX_TOPO], bond_j[LIB_MAX_TOPO];
  int n_bonds;

  char angle_id[LIB_MAX_TOPO][16];
  int angle_i[LIB_MAX_TOPO], angle_j[LIB_MAX_TOPO], angle_k[LIB_MAX_TOPO];
  int n_angles;
} LIB_MOL;

// What callers need to size a system before building it
typedef struct {
  int n_beads;       // the molecule's own beads; -1 if there is no such file
  int n_beads_total; // plus one set of beads per declared counterion
  double M_w;        // 0 when the file gives none
  bool bilayer;
  LIB_CION cion[LIB_MAX_CION];
  int n_cion;
} LIB_MOL_INFO;

// Read library metadata: all bead types, bond/angle types, and interactions
LIBRARY ReadLibrary(const char *lib_dir);

// Fill pot[i][j][{A,Rc,gamma}] from library interactions using bead type names
void FillPotFromLibrary(const LIBRARY *lib, const SYSTEM *System, ArrNDd *pot);

// Read one molecule file. False only when the file does not exist; a file that
// exists but does not parse is fatal.
bool ReadLibraryMolFile(const char *lib_dir, const char *name, LIB_MOL *mol);

/*
 * Add n_mols copies of mol_name to lib->System. With with_cion, its declared
 * counterions come too, as molecules in their own right - n_mols * count of
 * each. Without it the caller places them, which is what the bilayer workflow
 * needs: a counterion belongs in the water, not in the leaflet its parent is
 * being constrained into.
 *
 * A counterion contributes its beads and its charge and nothing else: its own
 * counterion line is not followed, and its M_w is not consulted, because the
 * parent's M_w is the mass of the whole salt. Resolution is therefore always
 * exactly one level deep.
 */
void ReadLibraryMolecule(const char *lib_dir, const char *mol_name,
                         int n_mols, bool with_cion, LIBRARY *lib);

// Bead counts, mass, role and counterions of mol_name, from its own file
LIB_MOL_INFO LibraryMoleculeInfo(const char *lib_dir, const char *mol_name);

// Rename bead types in sys to library names using molecule type name matching.
// For each molecule type in sys whose name matches a library molecule file,
// each bead type used by that molecule is renamed to the library's bead type.
void RenameBeadTypesFromLibrary(SYSTEM *sys, const LIBRARY *lib,
                                const char *lib_dir);

void FreeLibrary(LIBRARY *lib);

#endif
