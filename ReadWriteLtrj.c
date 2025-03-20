#include "ReadWriteLtrj.h"
#include "General.h"

// maximum number of variables in 'ITEM: ATOM' line
static const int max_var = 12;

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
                             char vars[max_var][10], int *line_count);
// read an atom coordinate line
static int LtrjReadCoorLine(FILE *fr, BEAD *b, int b_count,
                            const int *var, int cols);
// fill a helper array with possible variables in 'ITEM: ATOMS ...' line
static void LtrjFillAtomVariables(char var[max_var][10]);
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
  int position[max_var];
  char var[max_var][10];
  int cols = LtrjReadAtomsLine(fr, file, position, var, &line_count);
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
    if (LtrjReadCoorLine(fr, &line, Sys.Count.Bead, position, cols) < 0) {
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
      if (position[11] != -1) { // 'type' variable is present
        s_strcpy(name, split[position[11]], BEAD_NAME);
      } else if (position[1] != -1) { // 'element' variable is present
        s_strcpy(name, split[position[1]], BEAD_NAME);
      }
      NewBeadType(&Sys.BeadType, &Count->BeadType, name, CHARGE, MASS, RADIUS);
      BEADTYPE *bt_new = &Sys.BeadType[type];
      bt_new->Number = 1;
      b->Type = type;
    }
  } //}}}
  fclose(fr);
  FillSystemNonessentials(&Sys, false);
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
  int position[max_var];
  char vars[max_var][10];
  int cols = LtrjReadAtomsLine(fr, file, position, vars, line_count);
  // read atom lines //{{{
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    BEAD line;
    (*line_count)++;
    if (LtrjReadCoorLine(fr, &line, System->Count.Bead, position, cols) < 0) {
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
  ChangeBoxByLow(System, -1);
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
    double bounds[2][3];
    for (int dd = 0; dd < 3; dd++) {
      (*line_count)++;
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'ITEM: BOX BOUNDS' section");
        return -2;
      }
      if (words < 2 ||
          !IsRealNumber(split[0], &bounds[0][dd]) ||
          !IsRealNumber(split[1], &bounds[1][dd]) ||
          bounds[1][dd] <= bounds[0][dd]) {
        err_msg("wrong line in 'ITEM: BOX BOUNDS' section");
        PrintErrorFileLine(file, *line_count);
        return -1;
      }
      box->OrthoLength[dd] = bounds[1][dd] - bounds[0][dd];
      box->Low[dd] = bounds[0][dd];
    }
    CalculateBoxData(box, 1);
  } else if (strcmp(split[3], "xy") == 0) { // triclinic box
    double bounds[2][3];
    double tilt[3];
    for (int dd = 0; dd < 3; dd++) {
      (*line_count)++;
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'ITEM: BOX BOUNDS' section");
        return -2;
      }
      if (words < 3 ||
          !IsRealNumber(split[0], &bounds[0][dd]) ||
          !IsRealNumber(split[1], &bounds[1][dd]) ||
          bounds[1][dd] <= bounds[0][dd] ||
          !IsRealNumber(split[2], &tilt[dd])) {
        err_msg("wrong pbc line");
        PrintErrorFileLine(file, *line_count);
        return -1;
      }
    }
    // see https://docs.lammps.org/Howto_triclinic.html
    double from_bound[2][3];
    from_bound[0][0] = Min3(0, tilt[0], Min3(0, tilt[1], tilt[0] + tilt[1]));
    from_bound[1][0] = Max3(0, tilt[0], Max3(0, tilt[1], tilt[0] + tilt[1]));
    from_bound[0][1] = Min3(0, 0, tilt[2]);
    from_bound[1][1] = Max3(0, 0, tilt[2]);
    from_bound[0][2] = 0;
    from_bound[1][2] = 0;
    for (int dd = 0; dd < 3; dd++) {
      box->OrthoLength[dd] = (bounds[1][dd] - from_bound[1][dd]) -
                             (bounds[0][dd] - from_bound[0][dd]);
    }
    box->transform[0][1] = tilt[0];
    box->transform[0][2] = tilt[1];
    box->transform[1][2] = tilt[2];
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
                             char vars[max_var][10], int *line_count) {
  // generate array with possible variable names
  LtrjFillAtomVariables(vars);
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
  InitIntArray(var_pos, max_var, -1); // id, element, r[3], v[3], f[3], type
  int cols = -1;
  for (int i = 2; i < words; i++) {
    for (int j = 0; j < max_var; j++) {
      if (strcmp(split[i], vars[j]) == 0) {
        var_pos[j] = i - 2;
        cols = i - 2 + 1;
        break;
      }
    }
  }
  if (cols < 0 || var_pos[0] == -1) {
    err_msg("wrong 'ITEM: ATOMS' line");
    return -1;
  }
  return cols;
} //}}}
// LtrjReadCoorLine() //{{{
static int LtrjReadCoorLine(FILE *fr, BEAD *b, int b_count,
                            const int *var, int cols) {
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  InitBead(b);
  long id;
  if (words < cols || !IsWholeNumber(split[var[0]], &id) || id > b_count ||
      (var[ 2] != -1 && !IsRealNumber(split[var[ 2]], &b->Position[0])) ||
      (var[ 3] != -1 && !IsRealNumber(split[var[ 3]], &b->Position[1])) ||
      (var[ 4] != -1 && !IsRealNumber(split[var[ 4]], &b->Position[2])) ||
      (var[ 5] != -1 && !IsRealNumber(split[var[ 5]], &b->Velocity[0])) ||
      (var[ 6] != -1 && !IsRealNumber(split[var[ 6]], &b->Velocity[1])) ||
      (var[ 7] != -1 && !IsRealNumber(split[var[ 7]], &b->Velocity[2])) ||
      (var[ 8] != -1 && !IsRealNumber(split[var[ 8]], &b->Force[0])) ||
      (var[ 9] != -1 && !IsRealNumber(split[var[ 9]], &b->Force[1])) ||
      (var[10] != -1 && !IsRealNumber(split[var[10]], &b->Force[2]))) {
    return -1;
  }
  b->Type = id;
  return 1;
} //}}}
static void LtrjFillAtomVariables(char var[max_var][10]) { //{{{
  s_strcpy(var[0], "id", 10);
  s_strcpy(var[1], "element", 10);
  s_strcpy(var[2], "x", 10);
  s_strcpy(var[3], "y", 10);
  s_strcpy(var[4], "z", 10);
  s_strcpy(var[5], "vx", 10);
  s_strcpy(var[6], "vy", 10);
  s_strcpy(var[7], "vz", 10);
  s_strcpy(var[8], "fx", 10);
  s_strcpy(var[9], "fy", 10);
  s_strcpy(var[10], "fz", 10);
  s_strcpy(var[11], "type", 10);
} //}}}
static void AssignPosVelForce(const BEAD in, BEAD *b) { //{{{
  for (int dd = 0; dd < 3; dd++) {
    b->Position[dd] = in.Position[dd];
    b->Velocity[dd] = in.Velocity[dd];
    b->Force[dd] = in.Force[dd];
  }
} //}}}

// LtrjWriteCoor() //{{{
void LtrjWriteCoor(FILE *fw, const int step,
                   const bool *write, const SYSTEM System) {
  // find out number of beads to save and if velocity/force should be saved
  int count_write = 0;
  bool vel = false;
  bool force = false;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *b = &System.Bead[id];
    if (write[id]) {
      count_write++;
      if (b->Velocity[0] != 0 ||
          b->Velocity[1] != 0 ||
          b->Velocity[2] != 0) {
        vel = true;
      }
      if (b->Force[0] != 0 ||
          b->Force[1] != 0 ||
          b->Force[2] != 0) {
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
    if (box->alpha == 90 && box->beta == 90 && box->gamma == 90) {
      fprintf(fw, "ITEM: BOX BOUNDS pp pp pp\n");
      fprintf(fw, "%lf %lf\n", box->Low[0], box->Length[0]+box->Low[0]);
      fprintf(fw, "%lf %lf\n", box->Low[1], box->Length[1]+box->Low[1]);
      fprintf(fw, "%lf %lf\n", box->Low[2], box->Length[2]+box->Low[2]);
    } else {
      fprintf(fw, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding[0], box->transform[0][1]);
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding[1], box->transform[0][2]);
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding[2], box->transform[1][2]);
    }
    fprintf(fw, "ITEM: ATOMS id element x y z");
    if (vel) {
      fprintf(fw, " vx vy vz");
    }
    if (force) {
      fprintf(fw, " fx fy fz");
    }
    // fprintf(fw, " mol");
    putc('\n', fw);
    for (int i = 0; i < System.Count.BeadCoor; i++) {
      int id = System.BeadCoor[i];
      BEAD *b = &System.Bead[id];
      if (write[id]) {
        int type = b->Type;
        fprintf(fw, "%8d %8s %8.4f %8.4f %8.4f", id + 1,
                System.BeadType[type].Name,
                b->Position[0]+box->Low[0],
                b->Position[1]+box->Low[1],
                b->Position[2]+box->Low[2]);
        if (vel) {
          for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %8.4f", b->Velocity[dd]);
          }
        }
        if (force) {
          for (int dd = 0; dd < 3; dd++) {
            fprintf(fw, " %8.4f", b->Force[dd]);
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
