#include "ReadLibrary.h"
#include "General.h"

// TODO: I need to rewrite this

// static helpers
// build path to <library>/<file> //{{{
static void BuildPath(const char *lib_dir, const char *filename, char *out) {
  int len = strlen(lib_dir);
  if (len > 0 && lib_dir[len - 1] == '/') {
    snprintf(out, LINE, "%s%s", lib_dir, filename);
  } else {
    snprintf(out, LINE, "%s/%s", lib_dir, filename);
  }
} //}}}
// parse Fortran-style '0.19d0' exponent notation //{{{
static double ParseFortranDouble(const char *s) {
  char buf[64] = {0};
  for (int i = 0; i < 63 && s[i]; i++) {
    buf[i] = (s[i] == 'd' || s[i] == 'D') ? 'e' : s[i];
  }
  return strtod(buf, nullptr);
} //}}}
// ignore empty lines and comments and the first line with a number //{{{
static bool IgnoreLine (bool *found_count) {
  // empty and comment lines
  if (words == 0 || split[0][0] == '#') {
    return true;
  }
  // first valid line contains the muber of bond types (only check it exists)
  if (!*found_count) {
    long n;
    if (IsWholeNumber(split[0], &n)) {
      *found_count = true;
    }
    return true;
  }
  return false;
} //}}}
// read bead name from a single-bead molecule file; return true if n_beads==1
static bool read_single_bead_name(const char *lib_dir, const char *mol_name,
                                  char *bead_name_out) {
  char path[LINE], mol_file[LINE];
  snprintf(mol_file, LINE, "%s.txt", mol_name);
  BuildPath(lib_dir, mol_file, path);
  FILE *fr = fopen(path, "r");
  // skip non-existent molecule files
  if (!fr) {
    return false;
  }
  bool found_key = false,
       found_nbeads = false;
  long n_beads_file = -1;
  bead_name_out[0] = '\0';
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (words == 0 || split[0][0] == '#') {
      continue;
    }
    if (strncasecmp(split[0], "Bonds", 5) == 0 ||
        strncasecmp(split[0], "End",   3) == 0) {
      break;
    }
    if (!found_key) {
      found_key = true;
      continue;
    }
    if (!found_nbeads) {
      long n;
      if (IsWholeNumber(split[0], &n)) {
        n_beads_file = n;
        found_nbeads = true;
      }
      continue;
    }
    if (words >= 5 && bead_name_out[0] == '\0') {
      double x;
      if (!IsRealNumber(split[2], &x)) {
        continue;
      }
      s_strcpy(bead_name_out, split[1], BEAD_NAME);
    }
  }
  fclose(fr);
  if (n_beads_file == 1 && bead_name_out[0] != '\0') {
    return true;
  } else {
    return false;
  }
}

LIBRARY ReadLibrary(const char *lib_dir) { //{{{
  LIBRARY lib = {0};
  InitSystem(&lib.System);
  SYSTEM *Sys = &lib.System;
  COUNT *Count = &Sys->Count;
  // 1) read bead types and self-interactions from list_parameters.txt //{{{
  char path[LINE];
  double self_A[LIB_MAX_IDS] = {0},  // repulsion parameter
         self_Rc[LIB_MAX_IDS] = {0}; // bead radius
  int n_self = 0;
  BuildPath(lib_dir, "list_parameters.txt", path);
  FILE *fr = OpenFile(path, "r");
  bool found_count = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (IgnoreLine(&found_count)) {
      continue;
    }
    // line format for each bead type: index  bead_ID  mass  q  A  Rc
    double mass, charge, A, Rc;
    if (words < 6 ||
        !IsPosRealNumber(split[2], &mass) ||
        !IsRealNumber(split[3], &charge) ||
        !IsPosRealNumber(split[4], &A) ||
        !IsPosRealNumber(split[5], &Rc)) {
      continue;
    }
    NewBeadType(&Sys->BeadType, &Count->BeadType, split[1], charge, mass, Rc);
    if (n_self < LIB_MAX_IDS) {
      self_A[n_self]  = A;
      self_Rc[n_self] = Rc;
      n_self++;
    }
  }
  fclose(fr); //}}}
  // 2) read bond types from list_bonds.txt //{{{
  BuildPath(lib_dir, "list_bonds.txt", path);
  fr = OpenFile(path, "r");
  found_count = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (IgnoreLine(&found_count)) {
      continue;
    }
    // format: bond_ID  k  r0   (r0 may use Fortran 'd' exponent)
    // NOTE: k is lammps-style (k/2)
    double k = 0, r0 = 0;
    if (words < 3 || !IsPosRealNumber(split[1], &k)) {
      continue;
    }
    r0 = ParseFortranDouble(split[2]);
    if (r0 <= 0) {
      continue;
    }
    int idx = Count->BondType;
    Count->BondType++;
    Sys->BondType = s_realloc(Sys->BondType,
                              sizeof *Sys->BondType * Count->BondType);
    Sys->BondType[idx] = (PARAMS){2 * k, r0, 0, 0};
    if (lib.n_bond_ids < LIB_MAX_IDS) {
      s_strcpy(lib.bond_id[lib.n_bond_ids].id, split[0], 16);
      lib.bond_id[lib.n_bond_ids].index = idx;
      lib.n_bond_ids++;
    }
  }
  fclose(fr); //}}}
  // 3) read angle types from list_angles.txt //{{{
  BuildPath(lib_dir, "list_angles.txt", path);
  fr = OpenFile(path, "r");
  found_count = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (IgnoreLine(&found_count)) {
      continue;
    }
    // format: angle_ID  k  theta; k is lammps-style (k/2)
    double k = 0, theta = 0;
    if (words < 3 ||
        !IsPosRealNumber(split[1], &k) ||
        !IsPosRealNumber(split[2], &theta)) {
      continue;
    }
    int idx = Count->AngleType;
    Count->AngleType++;
    Sys->AngleType = s_realloc(Sys->AngleType,
                               sizeof *Sys->AngleType * Count->AngleType);
    Sys->AngleType[idx] = (PARAMS){2 * k, theta, 0, 0};
    if (lib.n_angle_ids < LIB_MAX_IDS) {
      s_strcpy(lib.angle_id[lib.n_angle_ids].id, split[0], 16);
      lib.angle_id[lib.n_angle_ids].index = idx;
      lib.n_angle_ids++;
    }
  }
  fclose(fr); //}}}
  // 4) build self-interactions array //{{{
  int n_bt = lib.System.Count.BeadType;
  int n_inter_alloc = n_bt + LIB_MAX_INTER;
  lib.inter = calloc(n_inter_alloc, sizeof *lib.inter);
  lib.n_inter = 0;
  for (int i = 0; i < n_self && i < n_bt; i++) {
    s_strcpy(lib.inter[lib.n_inter].name1, Sys->BeadType[i].Name, BEAD_NAME);
    s_strcpy(lib.inter[lib.n_inter].name2, Sys->BeadType[i].Name, BEAD_NAME);
    lib.inter[lib.n_inter].A = self_A[i];
    lib.inter[lib.n_inter].Rc = self_Rc[i];
    lib.inter[lib.n_inter].gamma = 4.5;
    lib.n_inter++;
  } //}}}
  // 5) read cross-interactions from list_cross_interactions.txt //{{{
  BuildPath(lib_dir, "list_cross_interactions.txt", path);
  fr = OpenFile(path, "r");
  found_count = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (IgnoreLine(&found_count)) {
      continue;
    }
    // format: index  name1  name2  A  Rc  [# comment]
    double A = 0, Rc = 0;
    if (words < 5 ||
        !IsPosRealNumber(split[3], &A) ||
        !IsPosRealNumber(split[4], &Rc)) {
      continue;
    }
    if (lib.n_inter >= n_inter_alloc) {
      n_inter_alloc = lib.n_inter + 64;
      lib.inter = s_realloc(lib.inter, n_inter_alloc * sizeof *lib.inter);
    }
    s_strcpy(lib.inter[lib.n_inter].name1, split[1], BEAD_NAME);
    s_strcpy(lib.inter[lib.n_inter].name2, split[2], BEAD_NAME);
    lib.inter[lib.n_inter].A = A;
    lib.inter[lib.n_inter].Rc = Rc;
    lib.inter[lib.n_inter].gamma = 4.5;
    lib.n_inter++;
  }
  fclose(fr); //}}}
  return lib;
} //}}}
// using bead type names, fill a potential array //{{{
void FillPotFromLibrary(const LIBRARY *lib, const SYSTEM *System, ArrNDd *pot) {
  int n_bt = System->Count.BeadType;
  for (int i = 0; i < n_bt; i++) {
    const char *n1 = System->BeadType[i].Name;
    for (int j = i; j < n_bt; j++) {
      const char *n2 = System->BeadType[j].Name;
      double A = 25.0,
             Rc = 1.0,
             gamma = 4.5;
      for (int k = 0; k < lib->n_inter; k++) {
        const LIB_INTERACTION *li = &lib->inter[k];
        if ((strcmp(li->name1, n1) == 0 && strcmp(li->name2, n2) == 0) ||
            (strcmp(li->name1, n2) == 0 && strcmp(li->name2, n1) == 0)) {
          A = li->A;
          Rc = li->Rc;
          gamma = li->gamma;
          break;
        }
      }
      SetArr3D(pot, i, j, 0, A);
      SetArr3D(pot, j, i, 0, A);
      SetArr3D(pot, i, j, 1, Rc);
      SetArr3D(pot, j, i, 1, Rc);
      SetArr3D(pot, i, j, 2, gamma);
      SetArr3D(pot, j, i, 2, gamma);
    }
  }
} //}}}
// add n_mols copies of mol_name to lib->System //{{{
void ReadLibraryMolecule(const char *lib_dir, const char *mol_name,
                         int n_mols, LIBRARY *lib) {
  SYSTEM *Sys = &lib->System;
  COUNT *Count = &Sys->Count;

  char mol_file[LINE], path[LINE];
  snprintf(mol_file, LINE, "%s.txt", mol_name);
  BuildPath(lib_dir, mol_file, path);
  FILE *fr = OpenFile(path, "r");
  // per-molecule temporary storage
  char bead_names[64][BEAD_NAME];
  vec3d bead_pos[64];
  char bond_ids_str[256][16];
  int bond_bi[256],
      bond_bj[256];
  char angle_ids_str[256][16];
  int angle_bi[256],
      angle_bj[256],
      angle_bk[256];
  int n_beads = 0,
      n_bonds = 0,
      n_angles = 0;
  // checks for section detection
  bool found_key = false,
       found_nbeads = false,
       in_bonds = false,
       in_angles = false,
       bonds_need_count = false,
       angles_need_count = false;
  // read the file in one loop: instead of separate bead type/bond/angle loops,
  // detect a keyword, keeping it until another keyword is encountered
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    // empty and comment lines
    if (words == 0 || split[0][0] == '#') {
      continue;
    }
    // detect section keyword //{{{
    if (strncasecmp(split[0], "Bonds", 5) == 0) {
      in_bonds = true;
      in_angles = false;
      bonds_need_count = true;
      continue;
    } else if (strncasecmp(split[0], "Angles", 6) == 0) {
      in_angles = true;
      in_bonds = false;
      angles_need_count = true;
      continue;
    } else if (strncasecmp(split[0], "End", 3) == 0) {
      break;
    } //}}}
    // skip the first non-comment non-empty line after a keyword  //{{{
    if (bonds_need_count) {
      bonds_need_count = false;
      continue;
    } else if (angles_need_count) {
      angles_need_count = false;
      continue;
    } //}}}
    // read bond parameters //{{{
    if (in_bonds) {
      // format: bond_ID  bead_i  bead_j
      if (words < 3) continue;
      if (n_bonds >= 256) {
        if (snprintf(ERROR_MSG, LINE, "molecule %s%s%s: more than 256 bonds",
                     ErrYellow(), mol_name, ErrRed()) < 0) {
          ErrorSnprintf();
        }
        PrintError();
        exit(1);
      }
      long bi, bj;
      if (!IsWholeNumber(split[1], &bi) ||
          !IsWholeNumber(split[2], &bj)) {
        continue;
      }
      if (bi < 1 || bj < 1) {
        if (snprintf(ERROR_MSG, LINE, "molecule %s%s%s: "
                     "bond bead index must be >= 1, got %s%ld%s and %s%ld%s",
                     ErrYellow(), mol_name, ErrRed(),
                     ErrYellow(), bi, ErrRed(),
                     ErrYellow(), bj, ErrRed()) < 0) {
          ErrorSnprintf();
        }
        PrintError();
        exit(1);
      }
      s_strcpy(bond_ids_str[n_bonds], split[0], 16);
      bond_bi[n_bonds] = bi;
      bond_bj[n_bonds] = bj;
      n_bonds++;
      continue;
    } //}}}
    // read andle parameters //{{{
    if (in_angles) {
      // format: angle_ID  bead_i  bead_j  bead_k
      if (words < 4) continue;
      if (n_angles >= 256) {
        if (snprintf(ERROR_MSG, LINE, "molecule %s%s%s: more than 256 angles",
                     ErrYellow(), mol_name, ErrRed()) < 0) {
          ErrorSnprintf();
        }
        PrintError();
        exit(1);
      }
      long bi, bj, bk;
      if (!IsWholeNumber(split[1], &bi) ||
          !IsWholeNumber(split[2], &bj) ||
          !IsWholeNumber(split[3], &bk)) {
        continue;
      }
      if (bi < 1 || bj < 1 || bk < 1) {
        if (snprintf(ERROR_MSG, LINE, "molecule %s%s%s: angle bead index "
                     "must be >= 1, got %s%ld%s, %s%ld%s/%s%ld%s",
                     ErrYellow(), mol_name, ErrRed(),
                     ErrYellow(), bi, ErrRed(),
                     ErrYellow(), bj, ErrRed(),
                     ErrYellow(), bk, ErrRed()) < 0) {
          ErrorSnprintf();
        }
        PrintError();
        exit(1);
      }
      s_strcpy(angle_ids_str[n_angles], split[0], 16);
      angle_bi[n_angles] = bi;
      angle_bj[n_angles] = bj;
      angle_bk[n_angles] = bk;
      n_angles++;
      continue;
    } //}}}
    // header section (no keyword found) //{{{
    if (!found_key) { // skip 'key value'
      found_key = true;
      continue;
    }
    if (!found_nbeads) { // skip 'number of beads' value
      long n;
      if (IsWholeNumber(split[0], &n)) {
        found_nbeads = true;
      }
      continue;
    }
    // bead line: index  bead_ID  x  y  z
    if (words >= 5) {
      if (n_beads >= 64) {
        if (snprintf(ERROR_MSG, LINE, "molecule %s%s%s: more than 64 beads",
                     ErrYellow(), mol_name, ErrRed()) < 0) {
          ErrorSnprintf();
        }
        PrintError();
        exit(1);
      }
      double x, y, z;
      if (!IsRealNumber(split[2], &x) ||
          !IsRealNumber(split[3], &y) ||
          !IsRealNumber(split[4], &z)) {
        continue;
      }
      s_strcpy(bead_names[n_beads], split[1], BEAD_NAME);
      bead_pos[n_beads] = (vec3d){.v = {x, y, z}};
      n_beads++;
    } //}}}
  }
  fclose(fr);
  // 1-bead molecules are added as free (unbonded) beads with no MoleculeType
  if (n_beads == 1) {
    int bt = FindBeadType(bead_names[0], *Sys);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "bead type '%s%s%s' not in library"
                   " (molecule %s%s%s)", ErrYellow(), bead_names[0], ErrRed(),
                   ErrYellow(), mol_name, ErrRed()) < 0) ErrorSnprintf();
      PrintError(); exit(1);
    }
    Sys->Bead = s_realloc(Sys->Bead,
                          (Count->Bead + n_mols) * sizeof *Sys->Bead);
    for (int m = 0; m < n_mols; m++) {
      int bid = Count->Bead + m;
      BEAD *bead = &Sys->Bead[bid];
      InitBead(bead);
      bead->Type = bt;
      bead->Molecule = -1;
      bead->InTimestep = true;
      bead->Position = bead_pos[0];
    }
    Count->Bead     += n_mols;
    Count->Unbonded += n_mols;
    Sys->BeadType[bt].Number += n_mols;
    return;
  }
  // Create MoleculeType entry
  int mt_idx = Count->MoleculeType;
  NewMolType(&Sys->MoleculeType, &Count->MoleculeType, (char *)mol_name,
             n_beads, n_bonds, n_angles, 0, 0);
  MOLECULETYPE *mt = &Sys->MoleculeType[mt_idx];
  mt->Number  = n_mols;
  mt->Named   = true;
  // Map bead names to BeadType indices
  for (int i = 0; i < n_beads; i++) {
    int bt = FindBeadType(bead_names[i], *Sys);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "bead type '%s%s%s' not in library"
                   "(molecule %s%s%s)", ErrYellow(), bead_names[i], ErrRed(),
                   ErrYellow(), mol_name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    mt->Bead[i] = bt;
  }
  // map bond ids to BondType indices
  for (int i = 0; i < n_bonds; i++) {
    mt->Bond[i][0] = bond_bi[i] - 1;
    mt->Bond[i][1] = bond_bj[i] - 1;
    int bt_idx = -1;
    for (int k = 0; k < lib->n_bond_ids; k++) {
      if (strcmp(lib->bond_id[k].id, bond_ids_str[i]) == 0) {
        bt_idx = lib->bond_id[k].index; break;
      }
    }
    mt->Bond[i][2] = bt_idx;
  }
  // map angle ids to AngleType indices
  for (int i = 0; i < n_angles; i++) {
    mt->Angle[i][0] = angle_bi[i] - 1;
    mt->Angle[i][1] = angle_bj[i] - 1;
    mt->Angle[i][2] = angle_bk[i] - 1;
    int at_idx = -1;
    for (int k = 0; k < lib->n_angle_ids; k++) {
      if (strcmp(lib->angle_id[k].id, angle_ids_str[i]) == 0) {
        at_idx = lib->angle_id[k].index; break;
      }
    }
    mt->Angle[i][3] = at_idx;
  }
  // grow Bead[] and Molecule[] arrays
  Sys->Bead = s_realloc(Sys->Bead, (Count->Bead + n_mols * n_beads) *
                        sizeof *Sys->Bead);
  Sys->Molecule = s_realloc(Sys->Molecule, (Count->Molecule + n_mols) *
                            sizeof *Sys->Molecule);
  for (int m = 0; m < n_mols; m++) {
    int mol_id = Count->Molecule + m;
    MOLECULE *mol = &Sys->Molecule[mol_id];
    InitMolecule(mol);
    mol->Type = mt_idx;
    mol->InTimestep = true;
    mol->Index = mol_id;
    mol->Bead = calloc(n_beads, sizeof *Sys->Molecule[mol_id].Bead);
    for (int b = 0; b < n_beads; b++) {
      int bid = Count->Bead + m * n_beads + b;
      BEAD *bead = &Sys->Bead[bid];
      InitBead(bead);
      bead->Type = mt->Bead[b];
      bead->Molecule = mol_id;
      bead->InTimestep = true;
      bead->Position = bead_pos[b];
      mol->Bead[b] = bid;
    }
  }
  Count->Molecule += n_mols;
  Count->Bead += n_mols * n_beads;
  Count->Bonded += n_mols * n_beads;
  for (int b = 0; b < n_beads; b++) {
    Sys->BeadType[mt->Bead[b]].Number += n_mols;
  }
} //}}}
// Read parent molecule and counterion together into a single MoleculeType.
// All beads per molecule are laid out contiguously: parent beads first, then
// counterion bead(s), so the LAMMPS data file writer sees them as one block. //{{{
void ReadLibraryMoleculeWithCion(const char *lib_dir, const char *mol_name,
                                  const char *cion_name, int n_mols,
                                  LIBRARY *lib) {
  SYSTEM *Sys = &lib->System;
  COUNT *Count = &Sys->Count;

  // --- read parent molecule file (identical parsing to ReadLibraryMolecule) ---
  char mol_file[LINE], path[LINE];
  snprintf(mol_file, LINE, "%s.txt", mol_name);
  BuildPath(lib_dir, mol_file, path);
  FILE *fr = OpenFile(path, "r");
  char bead_names[64][BEAD_NAME];
  vec3d bead_pos[64];
  char bond_ids_str[256][16];
  int bond_bi[256], bond_bj[256];
  char angle_ids_str[256][16];
  int angle_bi[256], angle_bj[256], angle_bk[256];
  int n_beads = 0, n_bonds = 0, n_angles = 0;
  bool found_key = false, found_nbeads = false,
       in_bonds = false, in_angles = false,
       bonds_need_count = false, angles_need_count = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (words == 0 || split[0][0] == '#') continue;
    if (strncasecmp(split[0], "Bonds",  5) == 0) {
      in_bonds = true; in_angles = false; bonds_need_count = true; continue;
    } else if (strncasecmp(split[0], "Angles", 6) == 0) {
      in_angles = true; in_bonds = false; angles_need_count = true; continue;
    } else if (strncasecmp(split[0], "End", 3) == 0) {
      break;
    }
    if (bonds_need_count)  { bonds_need_count  = false; continue; }
    if (angles_need_count) { angles_need_count = false; continue; }
    if (in_bonds) {
      if (words < 3) continue;
      long bi, bj;
      if (!IsWholeNumber(split[1], &bi) || !IsWholeNumber(split[2], &bj)) continue;
      s_strcpy(bond_ids_str[n_bonds], split[0], 16);
      bond_bi[n_bonds] = bi; bond_bj[n_bonds] = bj; n_bonds++;
      continue;
    }
    if (in_angles) {
      if (words < 4) continue;
      long bi, bj, bk;
      if (!IsWholeNumber(split[1], &bi) || !IsWholeNumber(split[2], &bj) ||
          !IsWholeNumber(split[3], &bk)) continue;
      s_strcpy(angle_ids_str[n_angles], split[0], 16);
      angle_bi[n_angles] = bi; angle_bj[n_angles] = bj;
      angle_bk[n_angles] = bk; n_angles++;
      continue;
    }
    if (!found_key)    { found_key    = true; continue; }
    if (!found_nbeads) {
      long n; if (IsWholeNumber(split[0], &n)) found_nbeads = true;
      continue;
    }
    if (words >= 5 && n_beads < 64) {
      double x, y, z;
      if (!IsRealNumber(split[2], &x) || !IsRealNumber(split[3], &y) ||
          !IsRealNumber(split[4], &z)) continue;
      s_strcpy(bead_names[n_beads], split[1], BEAD_NAME);
      bead_pos[n_beads] = (vec3d){.v = {x, y, z}};
      n_beads++;
    }
  }
  fclose(fr);

  // --- read counterion file (beads only; counterions have no bonds) ---
  char cion_bead_names[64][BEAD_NAME];
  vec3d cion_bead_pos[64];
  int n_cion = 0;
  snprintf(mol_file, LINE, "%s.txt", cion_name);
  BuildPath(lib_dir, mol_file, path);
  fr = OpenFile(path, "r");
  found_key = false; found_nbeads = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (words == 0 || split[0][0] == '#') continue;
    if (strncasecmp(split[0], "End",    3) == 0) break;
    if (strncasecmp(split[0], "Bonds",  5) == 0 ||
        strncasecmp(split[0], "Angles", 6) == 0) break;
    if (!found_key)    { found_key    = true; continue; }
    if (!found_nbeads) {
      long n; if (IsWholeNumber(split[0], &n)) found_nbeads = true;
      continue;
    }
    if (words >= 5 && n_cion < 64) {
      double x, y, z;
      if (IsRealNumber(split[2], &x) && IsRealNumber(split[3], &y) &&
          IsRealNumber(split[4], &z)) {
        s_strcpy(cion_bead_names[n_cion], split[1], BEAD_NAME);
        cion_bead_pos[n_cion] = (vec3d){.v = {x, y, z}};
        n_cion++;
      }
    }
  }
  fclose(fr);

  // --- build single MoleculeType with parent + counterion beads ---
  int total_beads = n_beads + n_cion;
  int mt_idx = Count->MoleculeType;
  NewMolType(&Sys->MoleculeType, &Count->MoleculeType, (char *)mol_name,
             total_beads, n_bonds, n_angles, 0, 0);
  MOLECULETYPE *mt = &Sys->MoleculeType[mt_idx];
  mt->Number = n_mols;
  mt->Named  = true;
  // parent bead types
  for (int i = 0; i < n_beads; i++) {
    int bt = FindBeadType(bead_names[i], *Sys);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "bead type '%s%s%s' not in library"
                   "(molecule %s%s%s)", ErrYellow(), bead_names[i], ErrRed(),
                   ErrYellow(), mol_name, ErrRed()) < 0) ErrorSnprintf();
      PrintError(); exit(1);
    }
    mt->Bead[i] = bt;
  }
  // counterion bead types
  for (int b = 0; b < n_cion; b++) {
    int bt = FindBeadType(cion_bead_names[b], *Sys);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "counterion bead type '%s%s%s' not in library",
                   ErrYellow(), cion_bead_names[b], ErrRed()) < 0) ErrorSnprintf();
      PrintError(); exit(1);
    }
    mt->Bead[n_beads + b] = bt;
  }
  // bond topology
  for (int i = 0; i < n_bonds; i++) {
    mt->Bond[i][0] = bond_bi[i] - 1;
    mt->Bond[i][1] = bond_bj[i] - 1;
    int bt_idx = -1;
    for (int k = 0; k < lib->n_bond_ids; k++) {
      if (strcmp(lib->bond_id[k].id, bond_ids_str[i]) == 0) {
        bt_idx = lib->bond_id[k].index; break;
      }
    }
    mt->Bond[i][2] = bt_idx;
  }
  // angle topology
  for (int i = 0; i < n_angles; i++) {
    mt->Angle[i][0] = angle_bi[i] - 1;
    mt->Angle[i][1] = angle_bj[i] - 1;
    mt->Angle[i][2] = angle_bk[i] - 1;
    int at_idx = -1;
    for (int k = 0; k < lib->n_angle_ids; k++) {
      if (strcmp(lib->angle_id[k].id, angle_ids_str[i]) == 0) {
        at_idx = lib->angle_id[k].index; break;
      }
    }
    mt->Angle[i][3] = at_idx;
  }
  // allocate Bead[] and Molecule[] arrays; each mol gets total_beads contiguous beads
  Sys->Bead = s_realloc(Sys->Bead,
                        (Count->Bead + n_mols * total_beads) * sizeof *Sys->Bead);
  Sys->Molecule = s_realloc(Sys->Molecule,
                            (Count->Molecule + n_mols) * sizeof *Sys->Molecule);
  for (int m = 0; m < n_mols; m++) {
    int mol_id = Count->Molecule + m;
    MOLECULE *mol = &Sys->Molecule[mol_id];
    InitMolecule(mol);
    mol->Type = mt_idx;
    mol->InTimestep = true;
    mol->Index = mol_id;
    mol->Bead = calloc(total_beads, sizeof *mol->Bead);
    // parent beads
    for (int b = 0; b < n_beads; b++) {
      int bid = Count->Bead + m * total_beads + b;
      BEAD *bead = &Sys->Bead[bid];
      InitBead(bead);
      bead->Type = mt->Bead[b];
      bead->Molecule = mol_id;
      bead->InTimestep = true;
      bead->Position = bead_pos[b];
      mol->Bead[b] = bid;
    }
    // counterion bead(s)
    for (int b = 0; b < n_cion; b++) {
      int bid = Count->Bead + m * total_beads + n_beads + b;
      BEAD *bead = &Sys->Bead[bid];
      InitBead(bead);
      bead->Type = mt->Bead[n_beads + b];
      bead->Molecule = mol_id;
      bead->InTimestep = true;
      bead->Position = cion_bead_pos[b];
      mol->Bead[n_beads + b] = bid;
    }
  }
  Count->Molecule += n_mols;
  Count->Bead     += n_mols * total_beads;
  Count->Bonded   += n_mols * total_beads;
  for (int b = 0; b < total_beads; b++) {
    Sys->BeadType[mt->Bead[b]].Number += n_mols;
  }
} //}}}
// get information about molecule type from list_molecules.txt file //{{{
LIB_MOL_INFO LibraryMoleculeInfo(const char *lib_dir, const char *mol_name) {
  LIB_MOL_INFO info = { .n_beads = -1, .cion = "" };
  char path[LINE];
  BuildPath(lib_dir, "list_molecules.txt", path);
  FILE *fr = OpenFile(path, "r");
  bool found_count = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (IgnoreLine(&found_count)) {
      continue;
    }
    // format: index  mol_name  Mw  n_beads  has_counterion  cion_name
    long nb, has_cion;
    if (words < 6 ||
        strcmp(split[1], mol_name) != 0) {
      continue;
    }
    if (IsWholeNumber(split[3], &nb)) {
      info.n_beads = nb;
    }
    if (IsWholeNumber(split[4], &has_cion) &&
        has_cion && strcmp(split[5], "None") != 0) {
      s_strcpy(info.cion, split[5], MOL_NAME);
    }
    break;
  }
  fclose(fr);
  return info;
} //}}}
// rename bead types in sys to library names //{{{
void RenameBeadTypesFromLibrary(SYSTEM *sys, const LIBRARY *lib,
                                const char *lib_dir) {
  char path[LINE];

  // pass 1) rename bead types used by named molecule types //{{{
  for (int mt = 0; mt < sys->Count.MoleculeType; mt++) {
    MOLECULETYPE *mt_sys = &sys->MoleculeType[mt];
    // skip unnamed molecules
    if (mt_sys->Name[0] == '\0') {
      continue;
    }
    char mol_file[LINE];
    snprintf(mol_file, LINE, "%s.txt", mt_sys->Name);
    BuildPath(lib_dir, mol_file, path);
    FILE *fr = fopen(path, "r");
    // skip non-existent molecule files
    if (!fr) {
      continue;
    }
    char bead_names[64][BEAD_NAME];
    int n_beads = 0;
    bool found_key = false, found_nbeads = false;
    // read name-based molecule file
    while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      // skip empty and comment lines
      if (words == 0 || split[0][0] == '#') {
        continue;
      }
      // end reading when a keyword is enctountered
      if (strncasecmp(split[0], "Bonds", 5) == 0 ||
          strncasecmp(split[0], "Angles", 6) == 0 ||
          strncasecmp(split[0], "End", 3) == 0) {
        break;
      }
      // skip 'key value' line
      if (!found_key) {
        found_key = true;
        continue;
      }
      // skip line with number of beads
      if (!found_nbeads) {
        long n;
        if (IsWholeNumber(split[0], &n)) {
          found_nbeads = true;
        }
        continue;
      }
      if (words >= 5 && n_beads < 64) {
        double x;
        if (!IsRealNumber(split[2], &x)) {
          continue;
        }
        s_strcpy(bead_names[n_beads], split[1], BEAD_NAME);
        n_beads++;
      }
    }
    fclose(fr);
    // if parent bead count doesn't match, try accounting for a merged counterion
    char cion_bead_names[64][BEAD_NAME];
    int n_cion = 0;
    if (n_beads != mt_sys->nBeads) {
      LIB_MOL_INFO info = LibraryMoleculeInfo(lib_dir, mt_sys->Name);
      if (info.cion[0] == '\0') continue;
      char cion_file[LINE], cion_path[LINE];
      snprintf(cion_file, LINE, "%s.txt", info.cion);
      BuildPath(lib_dir, cion_file, cion_path);
      FILE *cf = fopen(cion_path, "r");
      if (!cf) continue;
      bool ck = false, cn = false;
      while (ReadAndSplitLine(cf, SPL_STR, " \t\n")) {
        if (words == 0 || split[0][0] == '#') continue;
        if (strncasecmp(split[0], "Bonds",  5) == 0 ||
            strncasecmp(split[0], "Angles", 6) == 0 ||
            strncasecmp(split[0], "End",    3) == 0) break;
        if (!ck) { ck = true; continue; }
        if (!cn) { long nv; if (IsWholeNumber(split[0], &nv)) cn = true; continue; }
        if (words >= 5 && n_cion < 64) {
          double x;
          if (!IsRealNumber(split[2], &x)) continue;
          s_strcpy(cion_bead_names[n_cion], split[1], BEAD_NAME);
          n_cion++;
        }
      }
      fclose(cf);
      if (n_beads + n_cion != mt_sys->nBeads) continue;
    }
    // rename parent bead types
    for (int b = 0; b < n_beads; b++) {
      int sys_bt = mt_sys->Bead[b];
      int lib_bt = FindBeadType(bead_names[b], lib->System);
      if (lib_bt == -1) continue;
      BEADTYPE *lib_btp = &lib->System.BeadType[lib_bt];
      s_strcpy(sys->BeadType[sys_bt].Name, lib_btp->Name, BEAD_NAME);
      sys->BeadType[sys_bt].Charge = lib_btp->Charge;
      sys->BeadType[sys_bt].Mass   = lib_btp->Mass;
      sys->BeadType[sys_bt].Radius = lib_btp->Radius;
    }
    // rename counterion bead types
    for (int b = 0; b < n_cion; b++) {
      int sys_bt = mt_sys->Bead[n_beads + b];
      int lib_bt = FindBeadType(cion_bead_names[b], lib->System);
      if (lib_bt == -1) continue;
      BEADTYPE *lib_btp = &lib->System.BeadType[lib_bt];
      s_strcpy(sys->BeadType[sys_bt].Name, lib_btp->Name, BEAD_NAME);
      sys->BeadType[sys_bt].Charge = lib_btp->Charge;
      sys->BeadType[sys_bt].Mass   = lib_btp->Mass;
      sys->BeadType[sys_bt].Radius = lib_btp->Radius;
    }
  } //}}}
  // pass 2) rename free (unbonded) bead types
  // mark which bead type indices appear in at least one molecule type
  bool *in_mol = calloc(sys->Count.BeadType, sizeof *in_mol);
  if (!in_mol) {
    return;
  }
  for (int mt = 0; mt < sys->Count.MoleculeType; mt++) {
    MOLECULETYPE *m = &sys->MoleculeType[mt];
    for (int b = 0; b < m->nBeads; b++) {
      in_mol[m->Bead[b]] = true;
    }
  }

  // collect (bead_name, charge, expected_count) for single-bead library molecules:
  //   solvents (n_beads==1, no counterion): expected_count = -1 (match by charge only)
  //   counterions: expected_count = sum of parent molecule counts in sys
  //                (1 counterion bead per parent molecule)
  #define MAX_FREE_LIB 64
  char free_names[MAX_FREE_LIB][BEAD_NAME];
  char free_cion_mol[MAX_FREE_LIB][MOL_NAME]; // cion mol name (for count accumulation)
  double free_charges[MAX_FREE_LIB];
  int free_counts[MAX_FREE_LIB]; // -1 for solvents; sum of parent mols for counterions
  int n_free_lib = 0;

  BuildPath(lib_dir, "list_molecules.txt", path);
  FILE *f = fopen(path, "r");
  if (f) {
    bool found_count = false;

    while (ReadAndSplitLine(f, SPL_STR, " \t\n")) {
      if (words == 0 || split[0][0] == '#') continue;
      if (!found_count) {
        long n;
        if (IsWholeNumber(split[0], &n)) found_count = true;
        continue;
      }
      if (words < 6) continue;
      long nb, has_cion;
      if (!IsWholeNumber(split[3], &nb) || !IsWholeNumber(split[4], &has_cion)) continue;

      // single-bead solvent (no counterion, 1 bead): add once, count=-1
      if (!has_cion && nb == 1 && n_free_lib < MAX_FREE_LIB) {
        char bname[BEAD_NAME];
        if (read_single_bead_name(lib_dir, split[1], bname)) {
          int lib_bt = FindBeadType(bname, lib->System);
          if (lib_bt != -1) {
            s_strcpy(free_names[n_free_lib], bname, BEAD_NAME);
            free_cion_mol[n_free_lib][0] = '\0';
            free_charges[n_free_lib] = lib->System.BeadType[lib_bt].Charge;
            free_counts[n_free_lib++] = -1;
          }
        }
      }

      // counterion: accumulate parent molecule counts from sys
      if (has_cion && strcmp(split[5], "None") != 0) {
        const char *cion = split[5];
        // count parent molecules of this type present in sys
        int parent_count = 0;
        for (int mt = 0; mt < sys->Count.MoleculeType; mt++) {
          if (strcmp(sys->MoleculeType[mt].Name, split[1]) == 0)
            parent_count += sys->MoleculeType[mt].Number;
        }
        if (parent_count == 0) continue; // molecule not in this system
        // find existing entry for this cion or create one
        int idx = -1;
        for (int i = 0; i < n_free_lib; i++)
          if (strcmp(free_cion_mol[i], cion) == 0) { idx = i; break; }
        if (idx == -1 && n_free_lib < MAX_FREE_LIB) {
          char bname[BEAD_NAME];
          if (read_single_bead_name(lib_dir, cion, bname)) {
            int lib_bt = FindBeadType(bname, lib->System);
            if (lib_bt != -1) {
              idx = n_free_lib;
              s_strcpy(free_names[n_free_lib], bname, BEAD_NAME);
              s_strcpy(free_cion_mol[n_free_lib], cion, MOL_NAME);
              free_charges[n_free_lib] = lib->System.BeadType[lib_bt].Charge;
              free_counts[n_free_lib++] = 0;
            }
          }
        }
        if (idx != -1)
          free_counts[idx] += parent_count;
      }
    }
    fclose(f);
  }
  #undef MAX_FREE_LIB

  // for each free bead type in sys, match against the lib free list.
  // counterion entries (free_counts >= 0) are matched by charge+count;
  // solvent entries (free_counts == -1) are matched by charge alone.
  // charge+count match takes priority over charge-only.
  for (int bt = 0; bt < sys->Count.BeadType; bt++) {
    if (in_mol[bt]) continue;
    double q = sys->BeadType[bt].Charge;
    int cnt = sys->BeadType[bt].Number;
    // first try: charge + count (counterion entries)
    int match = -1;
    for (int i = 0; i < n_free_lib; i++) {
      if (free_counts[i] < 0) continue;
      if (fabs(free_charges[i] - q) < 0.01 && free_counts[i] == cnt) {
        if (match == -1) match = i;
        else { match = -1; break; } // ambiguous even with count
      }
    }
    // second try: charge only (solvent entries), if no counterion matched
    if (match == -1) {
      for (int i = 0; i < n_free_lib; i++) {
        if (free_counts[i] >= 0) continue;
        if (fabs(free_charges[i] - q) < 0.01) {
          if (match == -1) match = i;
          else { match = -1; break; }
        }
      }
    }
    if (match != -1) {
      int lib_bt = FindBeadType(free_names[match], lib->System);
      s_strcpy(sys->BeadType[bt].Name, free_names[match], BEAD_NAME);
      if (lib_bt != -1) {
        BEADTYPE *lib_btp = &lib->System.BeadType[lib_bt];
        sys->BeadType[bt].Charge = lib_btp->Charge;
        sys->BeadType[bt].Mass   = lib_btp->Mass;
        sys->BeadType[bt].Radius = lib_btp->Radius;
      }
    }
  }
  free(in_mol);
} //}}}
void FreeLibrary(LIBRARY *lib) { //{{{
  FreeSystem(&lib->System);
  free(lib->inter);
  lib->inter  = nullptr;
  lib->n_inter = 0;
} //}}}
