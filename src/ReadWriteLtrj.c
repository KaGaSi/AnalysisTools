#include "ReadWriteLtrj.h"
#include "System.h"
#include "Errors.h"

// TODO: sometimes, lammpstrj generated through lammps has ids>line in preamble
//       ...or is the norm if low-ids aren't saved?

// TODO: somehow optional writing of velocity, force, extra...
//       probably requires flags from outside ReadWrite files.

// maximum number of variables in 'ITEM: ATOM' line
static const int MAX_VAR = 22;

/*
 * Functions to read lammpstrj file (dump style custom) as a coordinate file via
 * LtrjReadTimestep() and LtrjSkipTimestep() and as a structure file via
 * LtrjReadStruct()
 */
// read timestep preamble, excluding 'ITEM: ATOMS' line
static int LtrjReadTimestepPreamble(FILE *fr, const char *file, BOX *box,
                                    int *line_count);
// test if next line is 'ITEM: TIMESTEP', then skip the section
static int LtrjSkipItemTimestep(FILE *fr, const char *file,
                                int *line_count);
// test if next line is 'ITEM: NUMBER OF ATOMS', then read the section
static int LtrjReadNumberOfAtoms(FILE *fr, const char *file,
                                 int *line_count);
// test if next line is 'ITEM: BOX BOUNDS', then read the section
static int LtrjReadPBCSection(FILE *fr, const char *file,
                              BOX *box, int *line_count);
// check if words & split contain 'ITEM: TIMESTEP' line
static bool LtrjCheckTimestepLine();
// check if words & split contain 'ITEM: NUMBER OF ATOMS' line
static bool LtrjCheckNumberAtomsLine();
// check if words & split contain 'ITEM: BOX BOUNDS ...' line
static bool LtrjCheckPbcLine();
// read 'ITEM: ATOMS ...' line, defining what variables are in which columns
static int LtrjReadAtomsLine(FILE *fr, const char *file, int *var_pos,
                             char vars[MAX_VAR][10], int unknown[6],
                             int *line_count);
// read an atom coordinate line
static int LtrjReadCoorLine(FILE *fr, BEAD *b, int b_count,
                            const int *var, int cols, int unknown[6]);
// fill a helper array with possible variables in 'ITEM: ATOMS ...' line
static void LtrjFillAtomVariables(char var[MAX_VAR][10]);
static void AssignPosVelForce(const BEAD in, BEAD *b);

// Use the first lammpstrj timestep as a definition of system composition //{{{
SYSTEM LtrjReadStruct(const char *file) {
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // read preamble, getting number of beads and box dimensions
  Count->Bead = LtrjReadTimestepPreamble(fr, file, &Sys.Box, &line_count);
  Count->BeadCoor = Count->Bead;
  Count->Unbonded = Count->Bead; // lammpstrj contains no bond information
  Count->UnbondedCoor = Count->Bead;
  Sys.Bead = s_realloc(Sys.Bead, sizeof *Sys.Bead * Count->Bead);
  Sys.BeadCoor = s_realloc(Sys.BeadCoor, sizeof *Sys.BeadCoor * Count->Bead);
  // read ITEM: ATOMS line & find positions of varaibles in a coordinate line
  int position[MAX_VAR];
  char var[MAX_VAR][10];
  int unknown[6];
  int cols = LtrjReadAtomsLine(fr, file, position, var, unknown, &line_count);
  // error - incorrect 'ITEM: ATOMS ...' line //{{{
  if (cols < 0) {
    err_msg("wrong 'ITEM: ATOMS' line in the first timestep");
    if (position[0] == -1) {
      char err[LINE];
      s_strcpy(err, ERROR_MSG, LINE);
      if (snprintf(ERROR_MSG, LINE, "%s (missing 'id' keyword)", err) < 0) {
        ErrorSnprintf();
      }
    }
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  } //}}}
  // read coordinate lines //{{{
  for (int i = 0; i < Count->Bead; i++) {
    BEAD line;
    InitBead(&line);
    line_count++;
    // read & check the coordinate line validity //{{{
    if (LtrjReadCoorLine(fr, &line, Sys.Count.Bead, position,
                         cols, unknown) < 0) {
      err_msg("invalid atom line (or not enough atom lines)");
      PrintErrorFileLine(file, line_count);
      exit(1);
    } //}}}
    int id = line.Type - 1;
    BEAD *b = &Sys.Bead[id];
    InitBead(b);
    AssignPosVelForce(line, b);
    if (b->InTimestep == true) {
      err_msg("multiple atoms with the same id");
      PrintErrorFileLine(file, line_count);
      exit(1);
    }
    b->InTimestep = true;
    Sys.BeadCoor[i] = id;
    /* find if the bead type exists based on 'element' variable
     *  (if 'element' is missing, all beads are of the same one type)
     */
    bool new = true;
    for (int j = 0; j < Count->BeadType; j++) {
      BEADTYPE *bt = &Sys.BeadType[j];
      if ((position[1] != -1 && strcmp(split[position[1]], bt->Name) == 0) ||
          (position[11] != -1 && strcmp(split[position[11]], bt->Name) == 0)) {
        bt->Number++;
        b->Type = j;
        new = false;
        break;
      }
    }
    if (new) { // create a new type
      int type = Count->BeadType;
      char name[BEAD_NAME];
      s_strcpy(name, "b0", BEAD_NAME);
      if (position[1] != -1) { // 'element' variable is present
        s_strcpy(name, split[position[1]], BEAD_NAME);
      } else if (position[11] != -1) { // 'type' variable is present
        s_strcpy(name, split[position[11]], BEAD_NAME);
      }
      NewBeadType(&Sys.BeadType, &Count->BeadType, name, CHARGE, MASS, RADIUS);
      BEADTYPE *bt_new = &Sys.BeadType[type];
      bt_new->Number = 1;
      b->Type = type;
    }
    // assign charge from 'q' column if present (split[] still holds this line)
    if (position[21] != -1) {
      double q;
      if (IsRealNumber(split[position[21]], &q)) {
        Sys.BeadType[b->Type].Charge = q;
      }
    }
  } //}}}
  fclose(fr);
  FillSystemNonessentials(&Sys, false);  // false for has_bonds
  // AllocFillBeadTypeIndex(&Sys);
  CheckSystem(Sys, file);
  ChangeBoxByLow(&Sys, -1);
  return Sys;
} //}}}
// Read a single timestep from lammpstrj file //{{{
int LtrjReadTimestep(FILE *fr, const char *file, SYSTEM *System,
                     int *line_count) {
  // set 'not in timestep' to all beads //{{{
  for (int i = 0; i < System->Count.Bead; i++) {
    System->Bead[i].InTimestep = false;
  } //}}}
  // set 'not in timestep' to all molecules //{{{
  System->Count.MoleculeCoor = 0;
  for (int i = 0; i < System->Count.Molecule; i++) {
    System->Molecule[i].InTimestep = false;
  } //}}}
  System->Count.BeadCoor = LtrjReadTimestepPreamble(fr, file, &System->Box,
                                                    line_count);
  if (System->Count.BeadCoor < 0) {
    return System->Count.BeadCoor;
  }
  // read ITEM: ATOMS line & find positions of varaibles in a coordinate line
  int position[MAX_VAR];
  char vars[MAX_VAR][10];
  int unknown[6];
  int cols = LtrjReadAtomsLine(fr, file, position, vars, unknown, line_count);
  if (cols < 0) { // error already printed by LtrjReadAtomsLine()
    return cols;
  }
  // read atom lines //{{{
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    BEAD line;
    (*line_count)++;
    if (LtrjReadCoorLine(fr, &line, System->Count.Bead, position,
                         cols, unknown) < 0) {
      err_msg("invalid atom line (or not enough atom lines)");
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
    int id = line.Type - 1;
    BEAD *b = &System->Bead[id];
    AssignPosVelForce(line, b);
    if (b->InTimestep) {
      err_msg("multiple atoms with the same id");
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
    b->InTimestep = true;
    System->BeadCoor[i] = id;
    if (b->Molecule != -1) {
      if (!System->Molecule[b->Molecule].InTimestep) {
        System->MoleculeCoor[System->Count.MoleculeCoor] = b->Molecule;
        System->Count.MoleculeCoor++;
      }
      System->Molecule[b->Molecule].InTimestep = true;
    }
  } //}}}
  // convert scaled (fractional) coordinates to Cartesian if needed
  // scaled vars: xs/ys/zs (12-14)
  //              xsu/ysu/zsu (18-20)
  //              xu/yu/zu (15-17) are Cartesian
  //              HUH???
  bool scaled = (position[18] != -1) ||
                (position[12] != -1 && position[15] == -1);
  if (scaled) {
    const BOX *box = &System->Box;
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      BEAD *b = &System->Bead[id];
      double sx = b->Position.v[0];
      double sy = b->Position.v[1];
      double sz = b->Position.v[2];
      b->Position.v[0] = box->transform[0][0] * sx +
                         box->transform[0][1] * sy +
                         box->transform[0][2] * sz;
      b->Position.v[1] = box->transform[1][1] * sy +
                         box->transform[1][2] * sz;
      b->Position.v[2] = box->transform[2][2] * sz;
    }
  }
  ChangeBoxByLow(System, -1);
  FillInCoor(System);
  return 1;
} //}}}
// TODO: skip lines based on the number of atoms?
// Skip a single timestep from lammpstrj file //{{{
int LtrjSkipTimestep(FILE *fr, const char *file, int *line_count) {
  /* read until two 'ITEM: TIMESTEP' lines are found
   *   the first should be the first line read, but who cares...
   *   the second is the beginning of the next timestep
   */
  fpos_t position;
  for (int i = 0; i < 2; i++) {
    do {
      fgetpos(fr, &position);
      (*line_count)++;
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        if (i == 0) {
          return -2; // before the first ITEM: TIMESTEP, so error
        } else {
          return 1; // there were some valid coor lines, so no error
        }
      }
    } while (words < 2 || strcmp(split[0], "ITEM:") != 0 ||
             strcmp(split[1], "TIMESTEP") != 0);
  }
  fsetpos(fr, &position); // restore the second 'ITEM: TIMESTEP' line
  (*line_count)--;        // the 'ITEM: TIMESTEP' will be re-read
  return 1;
} //}}}
// read pbc from the preamble of the first timestep //{{{
BOX LtrjReadPBC(const char *file) {
  BOX box = InitBox;
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  int test = LtrjSkipItemTimestep(fr, file, &line_count);
  if (test < 0) {
    if (test == -2) {
      ErrorEOF(file, "wrong 'ITEM: TIMESTEP' section");
    }
    exit(1);
  }
  test = LtrjReadNumberOfAtoms(fr, file, &line_count);
  if (test < 0) {
    exit(1);
  }
  test = LtrjReadPBCSection(fr, file, &box, &line_count);
  if (test < 0) {
    exit(1);
  }
  fclose(fr);
  return box;
} //}}}
// Helper functions for lammpstrj files
// LtrjReadTimestepPreamble() //{{{
static int LtrjReadTimestepPreamble(FILE *fr, const char *file, BOX *box,
                                    int *line_count) {
  int test = LtrjSkipItemTimestep(fr, file, line_count);
  if (test < 0) {
    return test;
  }
  int count = LtrjReadNumberOfAtoms(fr, file, line_count);
  if (count < 0) {
    return count;
  }
  test = LtrjReadPBCSection(fr, file, box, line_count);
  if (test < 0) {
    return test;
  }
  return count;
} //}}}
// LtrjSkipItemTimestep() //{{{
static int LtrjSkipItemTimestep(FILE *fr, const char *file,
                                int *line_count) {
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2; // proper eof - before the first line of a timestep
  }
  if (!LtrjCheckTimestepLine()) {
    err_msg("missing 'ITEM: TIMESTEP' line");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  // skip the timestep-counting line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing timestep in 'ITEM: TIMESTEP' section");
    return -2;
  }
  return 1;
} //}}}
// LtrjReadNumberOfAtoms() //{{{
static int LtrjReadNumberOfAtoms(FILE *fr, const char *file,
                                 int *line_count) {
  // read until 'ITEM: NUMBER OF ATOMS' line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing 'ITEM: NUMBER OF ATOMS' section");
    return -2;
  }
  if (!LtrjCheckNumberAtomsLine()) {
    err_msg("missing 'ITEM: NUMBER OF ATOMS' line");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  // read next line, i.e., the number of atoms
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing number of atoms");
    return -2;
  }
  long val = -1;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("number of atoms must be a non-zero whole number");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  return val;
} //}}}
// LtrjReadPBCSection() //{{{
static int LtrjReadPBCSection(FILE *fr, const char *file, BOX *box,
                              int *line_count) {
  // 1) read until 'ITEM:' line to find box type
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing 'ITEM: BOX BOUNDS' section");
    return -2;
  }
  if (!LtrjCheckPbcLine()) {
    err_msg("missing 'ITEM: BOX BOUNDS' line");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  // 2) read box dimensions
  if (strcmp(split[3], "pp") == 0 ||
      strcmp(split[3], "ff") == 0) { // orthogonal box
    vec3d bounds[2];
    for (int dd = 0; dd < 3; dd++) {
      (*line_count)++;
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'ITEM: BOX BOUNDS' section");
        return -2;
      }
      if (words < 2 ||
          !IsRealNumber(split[0], &bounds[0].v[dd]) ||
          !IsRealNumber(split[1], &bounds[1].v[dd]) ||
          bounds[1].v[dd] <= bounds[0].v[dd]) {
        err_msg("wrong line in 'ITEM: BOX BOUNDS' section");
        PrintErrorFileLine(file, *line_count);
        return -1;
      }
      box->OrthoLength.v[dd] = bounds[1].v[dd] - bounds[0].v[dd];
      box->Low.v[dd] = bounds[0].v[dd];
    }
    CalculateBoxData(box, 1);
  } else if (strcmp(split[3], "xy") == 0) { // triclinic box
    vec3d bounds[2], tilt;
    for (int dd = 0; dd < 3; dd++) {
      (*line_count)++;
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'ITEM: BOX BOUNDS' section");
        return -2;
      }
      if (words < 3 ||
          !IsRealNumber(split[0], &bounds[0].v[dd]) ||
          !IsRealNumber(split[1], &bounds[1].v[dd]) ||
          bounds[1].v[dd] <= bounds[0].v[dd] ||
          !IsRealNumber(split[2], &tilt.v[dd])) {
        err_msg("wrong pbc line");
        PrintErrorFileLine(file, *line_count);
        return -1;
      }
    }
    // see https://docs.lammps.org/Howto_triclinic.html
    vec3d from_bound[2];
    from_bound[0].x = Min3(0, tilt.v[0],
                           Min3(0, tilt.v[1], tilt.v[0] + tilt.v[1]));
    from_bound[1].x = Max3(0, tilt.v[0],
                           Max3(0, tilt.v[1], tilt.v[0] + tilt.v[1]));
    from_bound[0].y = Min3(0, 0, tilt.v[2]);
    from_bound[1].y = Max3(0, 0, tilt.v[2]);
    from_bound[0].z = 0;
    from_bound[1].z = 0;
    for (int dd = 0; dd < 3; dd++) {
      box->OrthoLength.v[dd] = (bounds[1].v[dd] - from_bound[1].v[dd]) -
                               (bounds[0].v[dd] - from_bound[0].v[dd]);
      box->Low.v[dd] = bounds[0].v[dd] - from_bound[0].v[dd];
    }
    box->transform[0][1] = tilt.v[0];
    box->transform[0][2] = tilt.v[1];
    box->transform[1][2] = tilt.v[2];
    CalculateBoxData(box, 1);
  } else { // not '<...> <...> <...> pp/xy ...' line
    err_msg("wrong ITEM: BOX BOUNDS line");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  return 1;
} //}}}
static bool LtrjCheckTimestepLine() { //{{{
  if (words >= 2 && strcmp(split[0], "ITEM:") == 0 &&
      strcmp(split[1], "TIMESTEP") == 0) {
    return true;
  } else {
    return false;
  }
} //}}}
static bool LtrjCheckNumberAtomsLine() { //{{{
  if (words >= 4 && strcmp(split[0], "ITEM:") == 0 &&
      strcmp(split[1], "NUMBER") == 0 && strcmp(split[2], "OF") == 0 &&
      strcmp(split[3], "ATOMS") == 0) {
    return true;
  } else {
    return false;
  }
} //}}}
static bool LtrjCheckPbcLine() { //{{{
  if (words >= 6 && strcmp(split[0], "ITEM:") == 0 &&
      strcmp(split[1], "BOX") == 0 && strcmp(split[2], "BOUNDS") == 0) {
    return true;
  } else {
    return false;
  }
} //}}}
// LtrjReadAtomsLine() //{{{
static int LtrjReadAtomsLine(FILE *fr, const char *file, int *var_pos,
                             char vars[MAX_VAR][10], int unknown[6],
                             int *line_count) {
  // generate array with possible variable names
  LtrjFillAtomVariables(vars);
  // initialize before any error return so callers never see garbage
  InitIntArray(var_pos, MAX_VAR, -1); // id, element, r[3], v[3], f[3], type
  InitIntArray(unknown, 6, -1);
  // read ITEM: ATOMS line //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing 'ITEM: ATOMS' line");
    return -2;
  }
  // error: line must be 'ITEMS: ATOMS <at least one more>'
  if (words < 3 || strcmp(split[0], "ITEM:") != 0 ||
      strcmp(split[1], "ATOMS") != 0) {
    err_msg("wrong 'ITEM: ATOMS ...' line");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }                                    //}}}
  int cols = -1;
  int count_unknown = 0;
  for (int i = 2; i < words; i++) {
    bool known = false;
    for (int j = 0; j < MAX_VAR; j++) {
      if (strcmp(split[i], vars[j]) == 0) {
        // column index: word position - 2 words (ITEMS: ATOMS)
        var_pos[j] = i - 2;
        // column count: column index + 1 for being count, not index
        cols = i - 2 + 1;
        known = true;
        break;
      }
    }
    if (!known && count_unknown < 6) {
      unknown[count_unknown] = i - 2;
      count_unknown++;
    }
  }
  // 'id' is mandatory
  if (var_pos[0] == -1) {
    err_msg("missing 'id' keyword in 'ITEM: ATOMS' line");
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  cols = words - 2; // count even the unknown columns
  return cols;
} //}}}
// LtrjReadCoorLine() //{{{
static int LtrjReadCoorLine(FILE *fr, BEAD *b, int b_count,
                            const int *var, int cols, int unknown[6]) {
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  InitBead(b);
  long id;
  // with regards to Position - unwrapped (xu/xsu) overwrite wrapped (x/xs)
  // atom ids are 1-based in a lammpstrj file
  if (words < cols || !IsNaturalNumber(split[var[0]], &id) || id > b_count ||
      (var[ 2] != -1 && !IsRealNumber(split[var[ 2]], &b->Position.v[0])) ||
      (var[ 3] != -1 && !IsRealNumber(split[var[ 3]], &b->Position.v[1])) ||
      (var[ 4] != -1 && !IsRealNumber(split[var[ 4]], &b->Position.v[2])) ||
      (var[ 5] != -1 && !IsRealNumber(split[var[ 5]], &b->Velocity.v[0])) ||
      (var[ 6] != -1 && !IsRealNumber(split[var[ 6]], &b->Velocity.v[1])) ||
      (var[ 7] != -1 && !IsRealNumber(split[var[ 7]], &b->Velocity.v[2])) ||
      (var[ 8] != -1 && !IsRealNumber(split[var[ 8]], &b->Force.v[0])) ||
      (var[ 9] != -1 && !IsRealNumber(split[var[ 9]], &b->Force.v[1])) ||
      (var[10] != -1 && !IsRealNumber(split[var[10]], &b->Force.v[2])) ||
      (var[12] != -1 && !IsRealNumber(split[var[12]], &b->Position.v[0])) ||
      (var[13] != -1 && !IsRealNumber(split[var[13]], &b->Position.v[1])) ||
      (var[14] != -1 && !IsRealNumber(split[var[14]], &b->Position.v[2])) ||
      (var[15] != -1 && !IsRealNumber(split[var[15]], &b->Position.v[0])) ||
      (var[16] != -1 && !IsRealNumber(split[var[16]], &b->Position.v[1])) ||
      (var[17] != -1 && !IsRealNumber(split[var[17]], &b->Position.v[2])) ||
      (var[18] != -1 && !IsRealNumber(split[var[18]], &b->Position.v[0])) ||
      (var[19] != -1 && !IsRealNumber(split[var[19]], &b->Position.v[1])) ||
      (var[20] != -1 && !IsRealNumber(split[var[20]], &b->Position.v[2]))) {
    return -1;
  }
  b->Type = id; // this will then be used to assign proper type to this bead

  for (int i = 0; i < 6; i++) {
    if (unknown[i] != -1) {
      IsRealNumber(split[unknown[i]], &b->Extra[i]);
    }
  }
  return 1;
} //}}}
static void LtrjFillAtomVariables(char var[MAX_VAR][10]) { //{{{
  s_strcpy(var[ 0], "id", 10);
  s_strcpy(var[ 1], "element", 10);
  s_strcpy(var[ 2], "x", 10);
  s_strcpy(var[ 3], "y", 10);
  s_strcpy(var[ 4], "z", 10);
  s_strcpy(var[ 5], "vx", 10);
  s_strcpy(var[ 6], "vy", 10);
  s_strcpy(var[ 7], "vz", 10);
  s_strcpy(var[ 8], "fx", 10);
  s_strcpy(var[ 9], "fy", 10);
  s_strcpy(var[10], "fz", 10);
  s_strcpy(var[11], "type", 10);
  s_strcpy(var[12], "xs", 10);
  s_strcpy(var[13], "ys", 10);
  s_strcpy(var[14], "zs", 10);
  s_strcpy(var[15], "xu", 10);
  s_strcpy(var[16], "yu", 10);
  s_strcpy(var[17], "zu", 10);
  s_strcpy(var[18], "xsu", 10);
  s_strcpy(var[19], "ysu", 10);
  s_strcpy(var[20], "zsu", 10);
  s_strcpy(var[21], "q",   10);
} //}}}
static void AssignPosVelForce(const BEAD in, BEAD *b) { //{{{
  for (int dd = 0; dd < 3; dd++) {
    b->Position.v[dd] = in.Position.v[dd];
    b->Velocity.v[dd] = in.Velocity.v[dd];
    b->Force.v[dd] = in.Force.v[dd];
  }
  for (int dd = 0; dd < 6; dd++) {
    b->Extra[dd] = in.Extra[dd];
  }
} //}}}

// LtrjWriteCoor() //{{{
void LtrjWriteCoor(FILE *fw, const int step,
                   const bool *write, const SYSTEM System) {
  // find out number of beads to save and if velocity/force should be saved
  int count_write = 0;
  bool vel = false; // TODO: define outside?
  bool force = false; // TODO: define outside?
  bool extra = false; // TODO: define outside?
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *b = &System.Bead[id];
    if (write[id]) {
      count_write++;
      if (b->Velocity.v[0] != 0 ||
          b->Velocity.v[1] != 0 ||
          b->Velocity.v[2] != 0) {
        vel = true;
      }
      if (b->Force.v[0] != 0 ||
          b->Force.v[1] != 0 ||
          b->Force.v[2] != 0) {
        force = true;
      }
    }
  }
  // print the step
  if (count_write > 0) {
    const BOX *box = &System.Box;
    fprintf(fw, "ITEM: TIMESTEP\n%d\n", step);
    fprintf(fw, "ITEM: NUMBER OF ATOMS\n%d\n", count_write);
    if (box->Volume == -1) {
      err_msg("unspecified box dimensions");
      PrintWarning();
    }
    // orthogonal box
    if (fabs(box->alpha -90) < 1e-3 &&
        fabs(box->beta  -90) < 1e-3 &&
        fabs(box->gamma -90) < 1e-3) {
      fprintf(fw, "ITEM: BOX BOUNDS pp pp pp\n");
      fprintf(fw, "%lf %lf\n", box->Low.x, box->Length.x + box->Low.x);
      fprintf(fw, "%lf %lf\n", box->Low.y, box->Length.y + box->Low.y);
      fprintf(fw, "%lf %lf\n", box->Low.z, box->Length.z + box->Low.z);
    } else {
      double lxy = box->transform[0][1];
      double lxz = box->transform[0][2];
      double lyz = box->transform[1][2];
      double lxyz = lxy + lxz;
      fprintf(fw, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
      fprintf(fw, "%lf %lf %lf\n",
              box->Low.x + Min3(0, lxy, Min3(0, lxz, lxyz)),
              box->Low.x + box->OrthoLength.x + Max3(0, lxy, Max3(0, lxz, lxyz)),
              lxy);
      fprintf(fw, "%lf %lf %lf\n",
              box->Low.y + Min3(0, 0.0, lyz),
              box->Low.y + box->OrthoLength.y + Max3(0, 0.0, lyz),
              lxz);
      fprintf(fw, "%lf %lf %lf\n",
              box->Low.z,
              box->Low.z + box->OrthoLength.z,
              lyz);
    }
    fprintf(fw, "ITEM: ATOMS id element x y z");
    if (vel) {
      fprintf(fw, " vx vy vz");
    }
    if (force) {
      fprintf(fw, " fx fy fz");
    }
    if (extra) {
      for (int dd = 0; dd < 6; dd++) {
        fprintf(fw, " extra[%d]", dd+1);
      }
    }
    // fprintf(fw, " mol");
    putc('\n', fw);
    for (int i = 0; i < System.Count.BeadCoor; i++) {
      int id = System.BeadCoor[i];
      BEAD *b = &System.Bead[id];
      if (write[id]) {
        int type = b->Type;
        fprintf(fw, "%8d %8s", id + 1, System.BeadType[type].Name);
        for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %8.4f", b->Position.v[dd] + box->Low.v[dd]);
        }
        if (vel) {
          for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %8.4f", b->Velocity.v[dd]);
          }
        }
        if (force) {
          for (int dd = 0; dd < 3; dd++) {
            fprintf(fw, " %8.4f", b->Force.v[dd]);
          }
        }
        if (extra) {
          for (int dd = 0; dd < 6; dd++) {
            fprintf(fw, " %8.4f", b->Extra[dd]);
          }
        }
        // fprintf(fw, " %5d", b->Molecule);
        putc('\n', fw);
      }
    }
  } else {
    err_msg("no beads to save");
    PrintWarning();
  }
} //}}}
