#include "ReadLibrary.h"
#include "General.h"

/*
 * The library grammar, used by every file here:
 *
 *   key value [value ...]     scalars, before the first block
 *   <keyword>                 opens a block
 *     idx col [col ...]       rows of that block, until the next keyword or EOF
 *
 * A line holding a single field is always a block header: every scalar carries
 * a value and every row carries an index plus at least one column. Comments run
 * from '#' to end of line, blank lines are insignificant, and the leading index
 * on a row is for the reader's benefit only - position in the block is the
 * index, so nothing here reads the printed one.
 */

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
// classify the line now in split[]/words; returns false for lines to skip //{{{
static bool LibLine(char *block, size_t block_size) {
  if (words == 0 || split[0][0] == '#') {
    return false;
  }
  if (words == 1) {
    s_strcpy(block, split[0], block_size);
    return false;
  }
  return true;
} //}}}
/*
 * One library candidate for a free (unbonded) bead type: the bead type name it
 * would get, the molecule that names it as a counterion (empty for a solvent),
 * the library charge, and the expected bead count - -1 when there is none,
 * which is what tells the two kinds apart when matching.
 */
typedef struct {
  char name[BEAD_NAME];
  char cion_of[MOL_NAME];
  double charge;
  int count;
} FREE_CAND;

// fatal: <file> is malformed. _Noreturn so the callers' parsed-value locals are
// not flagged as possibly-uninitialised on the failing branch. //{{{
_Noreturn static void LibError(const char *name, const char *what) {
  if (snprintf(ERROR_MSG, LINE, "library file %s%s%s: %s",
               ErrYellow(), name, ErrRed(), what) < 0) {
    ErrorSnprintf();
  }
  PrintError();
  exit(1);
} //}}}
// a bond/angle ID must fit, as a cut one could collide with another //{{{
static void CheckIdLen(const char *file, const char *id) {
  if (strlen(id) >= LIB_ID) {
    LibError(file, "bond/angle ID longer than LIB_ID characters");
  }
} //}}}

LIBRARY ReadLibrary(const char *lib_dir) { //{{{
  LIBRARY lib = {0};
  InitSystem(&lib.System);
  SYSTEM *Sys = &lib.System;
  COUNT *Count = &Sys->Count;
  char path[LINE], block[32];
  // 1) bead types and self-interactions from list_parameters.txt //{{{
  double *self_A = nullptr,  // repulsion parameter
         *self_Rc = nullptr; // bead radius
  int n_self = 0;
  BuildPath(lib_dir, "list_parameters.txt", path);
  FILE *fr = OpenFile(path, "r");
  block[0] = '\0';
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (!LibLine(block, sizeof block)) {
      continue;
    }
    if (strcmp(block, "bead_types") != 0) {
      continue;
    }
    // idx  bead_ID  mass  q  A  rc  [source ...]
    double mass, charge, A, Rc;
    if (words < 6 ||
        !IsPosRealNumber(split[2], &mass) ||
        !IsRealNumber(split[3], &charge) ||
        !IsPosRealNumber(split[4], &A) ||
        !IsPosRealNumber(split[5], &Rc)) {
      LibError("list_parameters.txt", "bead_types row must be "
               "'idx bead_ID mass q A rc [source]'");
    }
    NewBeadType(&Sys->BeadType, &Count->BeadType, split[1], charge, mass, Rc);
    self_A = s_realloc(self_A, (n_self + 1) * sizeof *self_A);
    self_Rc = s_realloc(self_Rc, (n_self + 1) * sizeof *self_Rc);
    self_A[n_self] = A;
    self_Rc[n_self] = Rc;
    n_self++;
  }
  fclose(fr); //}}}
  // 2) bond types from list_bonds.txt //{{{
  BuildPath(lib_dir, "list_bonds.txt", path);
  fr = OpenFile(path, "r");
  block[0] = '\0';
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (!LibLine(block, sizeof block)) {
      continue;
    }
    if (strcmp(block, "bond_types") != 0) {
      continue;
    }
    /*
     * idx  bond_ID  k_bond  r0
     *
     * Both the library and PARAMS.a use k of U = k (r - r0)^2 / 2, so k is
     * stored as read; it is lammps that folds the half into its K, and the
     * data writer is where that halving belongs.
     */
    double k = 0, r0 = 0;
    if (words < 4 ||
        !IsPosRealNumber(split[2], &k) ||
        !IsPosRealNumber(split[3], &r0)) {
      LibError("list_bonds.txt", "bond_types row must be "
               "'idx bond_ID k_bond r0'");
    }
    int idx = Count->BondType;
    Count->BondType++;
    Sys->BondType = s_realloc(Sys->BondType,
                              sizeof *Sys->BondType * Count->BondType);
    Sys->BondType[idx] = (PARAMS){k, r0, 0, 0};
    CheckIdLen("list_bonds.txt", split[1]);
    lib.bond_id = s_realloc(lib.bond_id,
                            (lib.n_bond_ids + 1) * sizeof *lib.bond_id);
    s_strcpy(lib.bond_id[lib.n_bond_ids].id, split[1], LIB_ID);
    lib.bond_id[lib.n_bond_ids].index = idx;
    lib.n_bond_ids++;
  }
  fclose(fr); //}}}
  // 3) angle types from list_angles.txt //{{{
  BuildPath(lib_dir, "list_angles.txt", path);
  fr = OpenFile(path, "r");
  block[0] = '\0';
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (!LibLine(block, sizeof block)) {
      continue;
    }
    if (strcmp(block, "angle_types") != 0) {
      continue;
    }
    // idx  angle_ID  k_angle  theta0    (k as read, as for bonds above)
    double k = 0, theta = 0;
    if (words < 4 ||
        !IsPosRealNumber(split[2], &k) ||
        !IsPosRealNumber(split[3], &theta)) {
      LibError("list_angles.txt", "angle_types row must be "
               "'idx angle_ID k_angle theta0'");
    }
    int idx = Count->AngleType;
    Count->AngleType++;
    Sys->AngleType = s_realloc(Sys->AngleType,
                               sizeof *Sys->AngleType * Count->AngleType);
    Sys->AngleType[idx] = (PARAMS){k, theta, 0, 0};
    CheckIdLen("list_angles.txt", split[1]);
    lib.angle_id = s_realloc(lib.angle_id,
                             (lib.n_angle_ids + 1) * sizeof *lib.angle_id);
    s_strcpy(lib.angle_id[lib.n_angle_ids].id, split[1], LIB_ID);
    lib.angle_id[lib.n_angle_ids].index = idx;
    lib.n_angle_ids++;
  }
  fclose(fr); //}}}
  // 4) build self-interactions array //{{{
  int n_bt = lib.System.Count.BeadType;
  int n_inter_alloc = n_bt + LIB_MAX_INTER;
  lib.inter = calloc(n_inter_alloc, sizeof *lib.inter);
  if (!lib.inter) {
    ErrorAlloc("lib.inter");
  }
  lib.n_inter = 0;
  for (int i = 0; i < n_self && i < n_bt; i++) {
    s_strcpy(lib.inter[lib.n_inter].name1, Sys->BeadType[i].Name, BEAD_NAME);
    s_strcpy(lib.inter[lib.n_inter].name2, Sys->BeadType[i].Name, BEAD_NAME);
    lib.inter[lib.n_inter].A = self_A[i];
    lib.inter[lib.n_inter].Rc = self_Rc[i];
    lib.inter[lib.n_inter].gamma = LIB_GAMMA;
    lib.n_inter++;
  }
  free(self_A);
  free(self_Rc); //}}}
  // 5) cross interactions from list_cross_interactions.txt //{{{
  BuildPath(lib_dir, "list_cross_interactions.txt", path);
  fr = OpenFile(path, "r");
  block[0] = '\0';
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (!LibLine(block, sizeof block)) {
      continue;
    }
    if (strcmp(block, "interactions") != 0) {
      continue;
    }
    // idx  bead_i  bead_j  A_ij  rc_ij  [source ...]
    double A = 0, Rc = 0;
    if (words < 5 ||
        !IsPosRealNumber(split[3], &A) ||
        !IsPosRealNumber(split[4], &Rc)) {
      LibError("list_cross_interactions.txt", "interactions row must be "
               "'idx bead_i bead_j A_ij rc_ij [source]'");
    }
    if (lib.n_inter >= n_inter_alloc) {
      n_inter_alloc = lib.n_inter + 64;
      lib.inter = s_realloc(lib.inter, n_inter_alloc * sizeof *lib.inter);
    }
    s_strcpy(lib.inter[lib.n_inter].name1, split[1], BEAD_NAME);
    s_strcpy(lib.inter[lib.n_inter].name2, split[2], BEAD_NAME);
    lib.inter[lib.n_inter].A = A;
    lib.inter[lib.n_inter].Rc = Rc;
    lib.inter[lib.n_inter].gamma = LIB_GAMMA;
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
      double A = LIB_DEFAULT_A,
             Rc = LIB_DEFAULT_RC,
             gamma = LIB_GAMMA;
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
// append the interactions block to a FIELD file //{{{
void AppendFieldInteractions(const char *file, const SYSTEM *System,
                             const LIBRARY *lib) {
  int n_bt = System->Count.BeadType;
  if (n_bt == 0) {
    return;
  }
  ArrNDd *pot = CreateArr3Dd(n_bt, n_bt, 3);
  FillPotFromLibrary(lib, System, pot);
  FILE *fw = OpenFile(file, "a");
  int n = n_bt * (n_bt - 1) / 2 + n_bt;
  fprintf(fw, "interactions %d <a_ij> <r_c> <gamma>\n", n);
  for (int i = 0; i < n_bt; i++) {
    for (int j = i; j < n_bt; j++) {
      fprintf(fw, "%10s %10s dpd %lf %lf %lf\n",
              System->BeadType[i].Name, System->BeadType[j].Name,
              GetArr3D(pot, i, j, 0), GetArr3D(pot, i, j, 1),
              GetArr3D(pot, i, j, 2));
    }
  }
  fclose(fw);
  FreeArrND(pot);
} //}}}
// read one molecule file //{{{
bool ReadLibraryMolFile(const char *lib_dir, const char *name, LIB_MOL *mol) {
  char mol_file[LINE], path[LINE];
  snprintf(mol_file, LINE, "%s.txt", name);
  BuildPath(lib_dir, mol_file, path);
  // zeroed before the file is even opened, so that a caller can FreeLibMol()
  // whatever it declared, whether the read found anything or not
  *mol = (LIB_MOL){0};
  FILE *fr = fopen(path, "r");
  if (!fr) {
    return false;
  }
  s_strcpy(mol->name, name, MOL_NAME);
  char block[32] = "";
  bool seen_role = false;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (!LibLine(block, sizeof block)) {
      continue;
    }
    // scalars, before the first block //{{{
    if (block[0] == '\0') {
      if (strcmp(split[0], "name") == 0) {
        s_strcpy(mol->name, split[1], MOL_NAME);
      } else if (strcmp(split[0], "role") == 0) {
        seen_role = true;
        if (strcmp(split[1], "bilayer") == 0) {
          mol->bilayer = true;
        } else if (strcmp(split[1], "soluble") == 0) {
          mol->bilayer = false;
        } else {
          LibError(mol_file, "role must be 'bilayer' or 'soluble'");
        }
      } else if (strcmp(split[0], "M_w") == 0) {
        double mw;
        if (!IsPosRealNumber(split[1], &mw)) {
          LibError(mol_file, "M_w must be a positive number");
        }
        mol->M_w = mw;
        mol->has_M_w = true;
      } else if (strcmp(split[0], "counterion") == 0) {
        mol->cion = s_realloc(mol->cion, (mol->n_cion + 1) * sizeof *mol->cion);
        LIB_CION *c = &mol->cion[mol->n_cion];
        s_strcpy(c->name, split[1], MOL_NAME);
        c->count = 1;
        if (words >= 3) {
          long n;
          if (!IsWholeNumber(split[2], &n) || n < 1) {
            LibError(mol_file, "counterion count must be a positive integer");
          }
          c->count = n;
        }
        mol->n_cion++;
      }
      // label and reference are for people; anything else is ignored so the
      // format can grow without old binaries refusing to load new files
      continue;
    } //}}}
    // beads: idx  bead_ID  x  y  z //{{{
    if (strcmp(block, "beads") == 0) {
      double x, y, z;
      if (words < 5 ||
          !IsRealNumber(split[2], &x) ||
          !IsRealNumber(split[3], &y) ||
          !IsRealNumber(split[4], &z)) {
        LibError(mol_file, "beads row must be 'idx bead_ID x y z'");
      }
      mol->bead_name = s_realloc(mol->bead_name,
                                 (mol->n_beads + 1) * sizeof *mol->bead_name);
      mol->bead_pos = s_realloc(mol->bead_pos,
                                (mol->n_beads + 1) * sizeof *mol->bead_pos);
      s_strcpy(mol->bead_name[mol->n_beads], split[1], BEAD_NAME);
      mol->bead_pos[mol->n_beads] = (vec3d){.v = {x, y, z}};
      mol->n_beads++;
      continue;
    } //}}}
    // bonds: idx  bond_ID  i  j //{{{
    if (strcmp(block, "bonds") == 0) {
      long bi, bj;
      if (words < 4 ||
          !IsWholeNumber(split[2], &bi) ||
          !IsWholeNumber(split[3], &bj)) {
        LibError(mol_file, "bonds row must be 'idx bond_ID i j'");
      }
      if (bi < 1 || bj < 1) {
        LibError(mol_file, "bond bead indices are 1-based");
      }
      CheckIdLen(mol_file, split[1]);
      int n = mol->n_bonds + 1;
      mol->bond_id = s_realloc(mol->bond_id, n * sizeof *mol->bond_id);
      mol->bond_i = s_realloc(mol->bond_i, n * sizeof *mol->bond_i);
      mol->bond_j = s_realloc(mol->bond_j, n * sizeof *mol->bond_j);
      s_strcpy(mol->bond_id[mol->n_bonds], split[1], LIB_ID);
      mol->bond_i[mol->n_bonds] = bi;
      mol->bond_j[mol->n_bonds] = bj;
      mol->n_bonds++;
      continue;
    } //}}}
    // angles: idx  angle_ID  i  j  k //{{{
    if (strcmp(block, "angles") == 0) {
      long bi, bj, bk;
      if (words < 5 ||
          !IsWholeNumber(split[2], &bi) ||
          !IsWholeNumber(split[3], &bj) ||
          !IsWholeNumber(split[4], &bk)) {
        LibError(mol_file, "angles row must be 'idx angle_ID i j k'");
      }
      if (bi < 1 || bj < 1 || bk < 1) {
        LibError(mol_file, "angle bead indices are 1-based");
      }
      CheckIdLen(mol_file, split[1]);
      int n = mol->n_angles + 1;
      mol->angle_id = s_realloc(mol->angle_id, n * sizeof *mol->angle_id);
      mol->angle_i = s_realloc(mol->angle_i, n * sizeof *mol->angle_i);
      mol->angle_j = s_realloc(mol->angle_j, n * sizeof *mol->angle_j);
      mol->angle_k = s_realloc(mol->angle_k, n * sizeof *mol->angle_k);
      s_strcpy(mol->angle_id[mol->n_angles], split[1], LIB_ID);
      mol->angle_i[mol->n_angles] = bi;
      mol->angle_j[mol->n_angles] = bj;
      mol->angle_k[mol->n_angles] = bk;
      mol->n_angles++;
      continue;
    } //}}}
    LibError(mol_file, "unknown block; expected beads, bonds or angles");
  }
  fclose(fr);
  if (mol->n_beads == 0) {
    LibError(mol_file, "no beads block");
  }
  if (!seen_role) {
    LibError(mol_file, "no role line");
  }
  // topology indices must land inside the bead list
  for (int i = 0; i < mol->n_bonds; i++) {
    if (mol->bond_i[i] > mol->n_beads || mol->bond_j[i] > mol->n_beads) {
      LibError(mol_file, "bond refers to a bead beyond the beads block");
    }
  }
  for (int i = 0; i < mol->n_angles; i++) {
    if (mol->angle_i[i] > mol->n_beads ||
        mol->angle_j[i] > mol->n_beads ||
        mol->angle_k[i] > mol->n_beads) {
      LibError(mol_file, "angle refers to a bead beyond the beads block");
    }
  }
  return true;
} //}}}
// release one molecule's arrays //{{{
void FreeLibMol(LIB_MOL *mol) {
  free(mol->cion);
  free(mol->bead_name);
  free(mol->bead_pos);
  free(mol->bond_id);
  free(mol->bond_i);
  free(mol->bond_j);
  free(mol->angle_id);
  free(mol->angle_i);
  free(mol->angle_j);
  free(mol->angle_k);
  *mol = (LIB_MOL){0};
} //}}}
// release the counterion list handed out by LibraryMoleculeInfo() //{{{
void FreeMolInfo(LIB_MOL_INFO *info) {
  free(info->cion);
  info->cion = nullptr;
  info->n_cion = 0;
} //}}}
// add n_mols copies of an already-read molecule to lib->System //{{{
static void AddMolecule(const LIB_MOL *mol, int n_mols, LIBRARY *lib) {
  SYSTEM *Sys = &lib->System;
  COUNT *Count = &Sys->Count;
  int n_beads = mol->n_beads;
  // single-bead molecules are free (unbonded) beads with no MoleculeType, which
  // is also what a monoatomic counterion should become
  if (n_beads == 1) {
    int bt = FindBeadType(mol->bead_name[0], *Sys);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "bead type '%s%s%s' not in library"
                   " (molecule %s%s%s)", ErrYellow(), mol->bead_name[0],
                   ErrRed(), ErrYellow(), mol->name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
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
      bead->Position = mol->bead_pos[0];
    }
    Count->Bead += n_mols;
    Count->Unbonded += n_mols;
    Sys->BeadType[bt].Number += n_mols;
    return;
  }
  int mt_idx = Count->MoleculeType;
  NewMolType(&Sys->MoleculeType, &Count->MoleculeType, (char *)mol->name,
             n_beads, mol->n_bonds, mol->n_angles, 0, 0);
  MOLECULETYPE *mt = &Sys->MoleculeType[mt_idx];
  mt->Number = n_mols;
  mt->Named = true;
  for (int i = 0; i < n_beads; i++) {
    int bt = FindBeadType(mol->bead_name[i], *Sys);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "bead type '%s%s%s' not in library"
                   " (molecule %s%s%s)", ErrYellow(), mol->bead_name[i],
                   ErrRed(), ErrYellow(), mol->name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    mt->Bead[i] = bt;
  }
  for (int i = 0; i < mol->n_bonds; i++) {
    mt->Bond[i][0] = mol->bond_i[i] - 1;
    mt->Bond[i][1] = mol->bond_j[i] - 1;
    int bt_idx = -1;
    for (int k = 0; k < lib->n_bond_ids; k++) {
      if (strcmp(lib->bond_id[k].id, mol->bond_id[i]) == 0) {
        bt_idx = lib->bond_id[k].index;
        break;
      }
    }
    if (bt_idx == -1) {
      if (snprintf(ERROR_MSG, LINE, "bond type '%s%s%s' not in list_bonds.txt"
                   " (molecule %s%s%s)", ErrYellow(), mol->bond_id[i],
                   ErrRed(), ErrYellow(), mol->name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    mt->Bond[i][2] = bt_idx;
  }
  for (int i = 0; i < mol->n_angles; i++) {
    mt->Angle[i][0] = mol->angle_i[i] - 1;
    mt->Angle[i][1] = mol->angle_j[i] - 1;
    mt->Angle[i][2] = mol->angle_k[i] - 1;
    int at_idx = -1;
    for (int k = 0; k < lib->n_angle_ids; k++) {
      if (strcmp(lib->angle_id[k].id, mol->angle_id[i]) == 0) {
        at_idx = lib->angle_id[k].index;
        break;
      }
    }
    if (at_idx == -1) {
      if (snprintf(ERROR_MSG, LINE, "angle type '%s%s%s' not in list_angles.txt"
                   " (molecule %s%s%s)", ErrYellow(), mol->angle_id[i],
                   ErrRed(), ErrYellow(), mol->name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    mt->Angle[i][3] = at_idx;
  }
  Sys->Bead = s_realloc(Sys->Bead, (Count->Bead + n_mols * n_beads) *
                        sizeof *Sys->Bead);
  Sys->Molecule = s_realloc(Sys->Molecule, (Count->Molecule + n_mols) *
                            sizeof *Sys->Molecule);
  for (int m = 0; m < n_mols; m++) {
    int mol_id = Count->Molecule + m;
    MOLECULE *new = &Sys->Molecule[mol_id];
    InitMolecule(new);
    new->Type = mt_idx;
    new->InTimestep = true;
    new->Index = mol_id;
    new->Bead = calloc(n_beads, sizeof *new->Bead);
    if (!new->Bead) {
      ErrorAlloc("Molecule.Bead");
    }
    for (int b = 0; b < n_beads; b++) {
      int bid = Count->Bead + m * n_beads + b;
      BEAD *bead = &Sys->Bead[bid];
      InitBead(bead);
      bead->Type = mt->Bead[b];
      bead->Molecule = mol_id;
      bead->InTimestep = true;
      bead->Position = mol->bead_pos[b];
      new->Bead[b] = bid;
    }
  }
  Count->Molecule += n_mols;
  Count->Bead += n_mols * n_beads;
  Count->Bonded += n_mols * n_beads;
  for (int b = 0; b < n_beads; b++) {
    Sys->BeadType[mt->Bead[b]].Number += n_mols;
  }
} //}}}
// add a molecule and, separately, its counterions //{{{
void ReadLibraryMolecule(const char *lib_dir, const char *mol_name,
                         int n_mols, bool with_cion, LIBRARY *lib) {
  LIB_MOL mol;
  if (!ReadLibraryMolFile(lib_dir, mol_name, &mol)) {
    if (snprintf(ERROR_MSG, LINE, "no library file for molecule %s%s%s",
                 ErrYellow(), mol_name, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  }
  AddMolecule(&mol, n_mols, lib);
  /*
   * The counterion is a molecule of its own, not beads appended to the parent.
   * Its own counterion line is deliberately not followed and its M_w is not
   * read: it contributes beads and charge only, so resolution stops here.
   */
  if (!with_cion) {
    FreeLibMol(&mol);
    return;
  }
  for (int i = 0; i < mol.n_cion; i++) {
    LIB_MOL cion;
    if (!ReadLibraryMolFile(lib_dir, mol.cion[i].name, &cion)) {
      if (snprintf(ERROR_MSG, LINE, "no library file for counterion %s%s%s"
                   " (molecule %s%s%s)", ErrYellow(), mol.cion[i].name,
                   ErrRed(), ErrYellow(), mol_name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    AddMolecule(&cion, n_mols * mol.cion[i].count, lib);
    FreeLibMol(&cion);
  }
  FreeLibMol(&mol);
} //}}}
// bead counts, mass, role and counterions of one molecule //{{{
LIB_MOL_INFO LibraryMoleculeInfo(const char *lib_dir, const char *mol_name) {
  LIB_MOL_INFO info = {.n_beads = -1, .n_beads_total = -1};
  LIB_MOL mol;
  if (!ReadLibraryMolFile(lib_dir, mol_name, &mol)) {
    return info;
  }
  info.n_beads = mol.n_beads;
  info.n_beads_total = mol.n_beads;
  info.M_w = mol.M_w;
  info.bilayer = mol.bilayer;
  info.n_cion = mol.n_cion;
  if (info.n_cion > 0) {
    info.cion = s_realloc(nullptr, info.n_cion * sizeof *info.cion);
  }
  for (int i = 0; i < mol.n_cion; i++) {
    info.cion[i] = mol.cion[i];
    LIB_MOL cion;
    if (!ReadLibraryMolFile(lib_dir, mol.cion[i].name, &cion)) {
      if (snprintf(ERROR_MSG, LINE, "no library file for counterion %s%s%s"
                   " (molecule %s%s%s)", ErrYellow(), mol.cion[i].name,
                   ErrRed(), ErrYellow(), mol_name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    info.n_beads_total += cion.n_beads * mol.cion[i].count;
    FreeLibMol(&cion);
  }
  FreeLibMol(&mol);
  return info;
} //}}}
// rename bead types in sys to library names //{{{
void RenameBeadTypesFromLibrary(SYSTEM *sys, const LIBRARY *lib,
                                const char *lib_dir) {
  // pass 1) bead types used by named molecule types //{{{
  for (int mt = 0; mt < sys->Count.MoleculeType; mt++) {
    MOLECULETYPE *mt_sys = &sys->MoleculeType[mt];
    if (mt_sys->Name[0] == '\0') {
      continue;
    }
    LIB_MOL mol;
    if (!ReadLibraryMolFile(lib_dir, mt_sys->Name, &mol)) {
      continue;
    }
    // counterions are their own molecules now, so the bead counts either match
    // or this is not the molecule the name suggests
    if (mol.n_beads != mt_sys->nBeads) {
      FreeLibMol(&mol);
      continue;
    }
    for (int b = 0; b < mol.n_beads; b++) {
      int sys_bt = mt_sys->Bead[b];
      int lib_bt = FindBeadType(mol.bead_name[b], lib->System);
      if (lib_bt == -1) {
        continue;
      }
      BEADTYPE *lib_btp = &lib->System.BeadType[lib_bt];
      s_strcpy(sys->BeadType[sys_bt].Name, lib_btp->Name, BEAD_NAME);
      sys->BeadType[sys_bt].Charge = lib_btp->Charge;
      sys->BeadType[sys_bt].Mass = lib_btp->Mass;
      sys->BeadType[sys_bt].Radius = lib_btp->Radius;
    }
    FreeLibMol(&mol);
  } //}}}
  // pass 2) free (unbonded) bead types //{{{
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
  /*
   * A LAMMPS data file has no bead names, so free beads are matched back to the
   * library by charge, and by count where the count is predictable. Walking the
   * directory is what replaced list_molecules.txt: every single-bead molecule
   * is a candidate, and every counterion of a molecule present in sys has an
   * expected count of (parent molecules) * (counterion count).
   */
  FREE_CAND *cand = nullptr;
  int n_free_lib = 0;

  DIR *dir = opendir(lib_dir);
  if (dir) {
    struct dirent *ent;
    while ((ent = readdir(dir)) != nullptr) {
      size_t len = strlen(ent->d_name);
      if (len < 5 || strcmp(ent->d_name + len - 4, ".txt") != 0) {
        continue;
      }
      if (strncmp(ent->d_name, "list_", 5) == 0) {
        continue;
      }
      char name[MOL_NAME];
      s_strcpy(name, ent->d_name, MOL_NAME);
      name[len - 4] = '\0';
      LIB_MOL mol;
      if (!ReadLibraryMolFile(lib_dir, name, &mol)) {
        continue;
      }
      /*
       * Any single-bead molecule is a free bead in the system, including one
       * that declares a counterion of its own - Na is a listed salt whose Na+
       * is loose in the water. It has no predictable count, so it can only be
       * matched on charge; leaving it out entirely is what made a Na+ bead get
       * renamed to the identically-charged H+.
       */
      if (mol.n_beads == 1) {
        int lib_bt = FindBeadType(mol.bead_name[0], lib->System);
        if (lib_bt != -1) {
          cand = s_realloc(cand, (n_free_lib + 1) * sizeof *cand);
          s_strcpy(cand[n_free_lib].name, mol.bead_name[0], BEAD_NAME);
          cand[n_free_lib].cion_of[0] = '\0';
          cand[n_free_lib].charge = lib->System.BeadType[lib_bt].Charge;
          cand[n_free_lib].count = -1;
          n_free_lib++;
        }
      }
      // this molecule's counterions, if it is in sys at all
      int parent_count = 0;
      for (int mt = 0; mt < sys->Count.MoleculeType; mt++) {
        if (strcmp(sys->MoleculeType[mt].Name, mol.name) == 0) {
          parent_count += sys->MoleculeType[mt].Number;
        }
      }
      if (parent_count == 0) {
        FreeLibMol(&mol);
        continue;
      }
      for (int c = 0; c < mol.n_cion; c++) {
        LIB_MOL cion;
        if (!ReadLibraryMolFile(lib_dir, mol.cion[c].name, &cion)) {
          continue;
        }
        // only a monoatomic counterion becomes a free bead
        if (cion.n_beads != 1) {
          FreeLibMol(&cion);
          continue;
        }
        int idx = -1;
        for (int i = 0; i < n_free_lib; i++) {
          if (strcmp(cand[i].cion_of, mol.cion[c].name) == 0) {
            idx = i;
            break;
          }
        }
        if (idx == -1) {
          int lib_bt = FindBeadType(cion.bead_name[0], lib->System);
          if (lib_bt != -1) {
            idx = n_free_lib;
            cand = s_realloc(cand, (n_free_lib + 1) * sizeof *cand);
            s_strcpy(cand[n_free_lib].name, cion.bead_name[0], BEAD_NAME);
            s_strcpy(cand[n_free_lib].cion_of, mol.cion[c].name, MOL_NAME);
            cand[n_free_lib].charge = lib->System.BeadType[lib_bt].Charge;
            cand[n_free_lib].count = 0;
            n_free_lib++;
          }
        }
        if (idx != -1) {
          cand[idx].count += parent_count * mol.cion[c].count;
        }
        FreeLibMol(&cion);
      }
      FreeLibMol(&mol);
    }
    closedir(dir);
  }

  // charge+count (counterions) wins over charge alone (solvents); an ambiguous
  // match renames nothing, because a wrong name is worse than no name
  for (int bt = 0; bt < sys->Count.BeadType; bt++) {
    if (in_mol[bt]) {
      continue;
    }
    double q = sys->BeadType[bt].Charge;
    int cnt = sys->BeadType[bt].Number;
    int match = -1;
    for (int i = 0; i < n_free_lib; i++) {
      if (cand[i].count < 0) {
        continue;
      }
      if (fabs(cand[i].charge - q) < 0.01 && cand[i].count == cnt) {
        if (match == -1) {
          match = i;
        } else {
          match = -1;
          break;
        }
      }
    }
    if (match == -1) {
      for (int i = 0; i < n_free_lib; i++) {
        if (cand[i].count >= 0) {
          continue;
        }
        if (fabs(cand[i].charge - q) < 0.01) {
          if (match == -1) {
            match = i;
          } else {
            match = -1;
            break;
          }
        }
      }
    }
    if (match != -1) {
      int lib_bt = FindBeadType(cand[match].name, lib->System);
      s_strcpy(sys->BeadType[bt].Name, cand[match].name, BEAD_NAME);
      if (lib_bt != -1) {
        BEADTYPE *lib_btp = &lib->System.BeadType[lib_bt];
        sys->BeadType[bt].Charge = lib_btp->Charge;
        sys->BeadType[bt].Mass = lib_btp->Mass;
        sys->BeadType[bt].Radius = lib_btp->Radius;
      }
    }
  }
  free(cand);
  free(in_mol); //}}}
} //}}}
void FreeLibrary(LIBRARY *lib) { //{{{
  FreeSystem(&lib->System);
  free(lib->inter);
  lib->inter = nullptr;
  lib->n_inter = 0;
  free(lib->bond_id);
  lib->bond_id = nullptr;
  lib->n_bond_ids = 0;
  free(lib->angle_id);
  lib->angle_id = nullptr;
  lib->n_angle_ids = 0;
} //}}}
