#ifndef READ_LIBRARY_H
#define READ_LIBRARY_H

// #include "AnalysisTools.h"
#include "Globals.h"
#include "Structs.h"
#include <dirent.h>

#define LIB_MAX_INTER 1024
// bond/angle type ID maximum length
#define LIB_ID 32

// defaults if missing from list_cross_interactions.txt or list_beadtypes.txt
// TODO: cross interactions should have some averages (arithmetic?)
#define LIB_DEF_GAMMA 4.5 // TODO: there should be gamma in the library files
#define LIB_DEF_A 25.0
#define LIB_DEF_RC 1.0
// structure for library bonds
typedef struct {
  char id[LIB_ID];
  int index; // index into System.BondType[]
} LIB_BOND_ID;
// structure for library angles
typedef struct {
  char id[LIB_ID];
  int index; // index into System.AngleType[]
} LIB_ANGLE_ID;
// structure for cross interactions
typedef struct {
  char name1[BEAD_NAME];
  char name2[BEAD_NAME];
  double A, Rc, gamma;
} LIB_INTERACTION;
// structure holding system information and bond/angle library information
typedef struct {
  SYSTEM System;
  int n_bond_ids;
  LIB_BOND_ID *bond_id;
  int n_angle_ids;
  LIB_ANGLE_ID *angle_id;
  int n_inter;
  LIB_INTERACTION *inter;
} LIBRARY;
// structure for library molecule's 'counterion <name> [count]' line
typedef struct {
  char name[MOL_NAME];
  int count;
} LIB_CION;
// structure for information about a single molecule frome <name>.txt file
typedef struct {
  char name[MOL_NAME];
  // TODO: the bilayer is Unilver case - delete?
  bool bilayer;   // role: bilayer (true) or soluble (false)
  double Mw;     // 0 when the file gives none
  bool has_Mw;
  LIB_CION *cion;
  int n_cion;

  char (*bead_name)[BEAD_NAME];
  vec3d *bead_pos;
  int n_beads;

  char (*bond_id)[LIB_ID];
  int *bond_i, *bond_j;
  int n_bonds;

  char (*angle_id)[LIB_ID];
  int *angle_i, *angle_j, *angle_k;
  int n_angles;
} LIB_MOL;

// What callers need to size a system before building it; cion is the caller's
// to FreeMolInfo()
// TODO: huh?
typedef struct {
  int n_beads;       // the molecule's own beads; -1 if there is no such file
  int n_beads_total; // plus one set of beads per declared counterion
  double M_w;        // 0 when the file gives none
  bool bilayer; // TODO: remove?
  int n_cion; // number of counterion lines
  LIB_CION *cion;
} LIB_MOL_INFO;
// Read library data: bead types, bond/angle types, and interactions
LIBRARY ReadLibrary(const char *lib_dir);
// Fill pot[i][j][{A,Rc,gamma}] from library interactions
void FillPotFromLibrary(const LIBRARY *lib, const SYSTEM *System, ArrNDd *pot);
// Append the DPD interactions to a FIELD file (WriteOutput() doesn't do that)
void AppendFieldInteractions(const char *file, const SYSTEM *System,
                             const LIBRARY *lib);
// read molecule file; returns false on nonexistent file, exit() on invalid file
bool ReadLibraryMolFile(const char *lib_dir, const char *name, LIB_MOL *mol);
// free structures; safe on both zeroed struct and one whose file is nonexistant
void FreeLibMol(LIB_MOL *mol);
void FreeMolInfo(LIB_MOL_INFO *info);

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

// fill bead counts, mass, and counterions of mol_name from its own file
LIB_MOL_INFO LibraryMoleculeInfo(const char *lib_dir, const char *mol_name);

// Rename bead types in sys to library names using molecule type name matching.
// For each molecule type in sys whose name matches a library molecule file,
// each bead type used by that molecule is renamed to the library's bead type.
void RenameBeadTypesFromLibrary(SYSTEM *sys, const LIBRARY *lib,
                                const char *lib_dir);

/*
 * What a utility's -lib does: read the library and let it name sys's bead
 * types and give them their masses, charges and radii. Call it right after
 * reading the structure, before any option that names a bead type.
 *
 * free_library true frees the library here and ignores lib (pass nullptr);
 * false hands it back in *lib, which the caller must then FreeLibrary().
 */
void ApplyLibraryToSystem(const char *lib_dir, SYSTEM *sys, LIBRARY *lib,
                          const bool free_library);

/*
 * Apply the two common options that change what a system's types are called:
 * -sys names the molecule types and -lib takes bead type names, masses,
 * charges and radii from a library (in that order, as the library matches
 * molecules by name). Does nothing when neither option is used.
 *
 * Call it right after ReadStructure() and before any option that names a type,
 * so that -bt and the like can use the library's names.
 *
 * free_library is passed through to ApplyLibraryToSystem(): true (what almost
 * every utility wants) frees the library here and ignores lib, false hands it
 * back in *lib for the caller to FreeLibrary(). Note that *lib is left alone
 * when -lib was not given at all, so a caller passing false should initialise
 * it and check commons.lib itself.
 */
void ApplyLibraryOptions(const COMMON_OPT commons, SYSTEM *sys, LIBRARY *lib,
                         const bool free_library);

void FreeLibrary(LIBRARY *lib);

#endif
