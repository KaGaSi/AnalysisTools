#include "Read.h"
#include "AnalysisTools.h"
#include "Errors.h"
#include "Structs.h"
#include "System.h"
#include "Write.h"
#include <string.h>

// TODO: use FillInCoor()

// TODO: LmpDataReadDihedralCoeffs() and LmpDataReadImproperCoeffs() should read
//       up to three numbers, not assuming any format of the potential

// STATIC DEFINITIONS
// variables defining line types //{{{
static const int ERROR_LINE = -1;
static const int BLANK_LINE = 0;
static const int COMMENT_LINE = 1;
static const int PBC_LINE = 2;
static const int PBC_LINE_ANGLES = 3;
static const int COOR_LINE = 4;
static const int COOR_LINE_O = 5;
static const int ATOM_LINE = 6;
static const int BOND_LINE = 7;
static const int TIME_LINE = 8;
static const int TIME_LINE_O = 9; //}}}
/*
 * Functions to read lammpstrj file (dump style custom) as a coordinate file via
 * LtrjReadTimestep() and LtrjSkipTimestep() and as a structure file via
 * LtrjReadStruct()
 */ //{{{
static int LtrjReadTimestep(FILE *fr, char file[], SYSTEM *System,
                            int *line_count);
static int LtrjSkipTimestep(FILE *fr, char file[], int *line_count);
static SYSTEM LtrjReadStruct(char file[]);
static BOX LtrjReadPBC(char file[]);
// Helper functions for lammpstrj files
// read timestep preamble, excluding 'ITEM: ATOMS' line
static int LtrjReadTimestepPreamble(FILE *fr, char file[], BOX *box,
                                    int *line_count);
// test if next line is 'ITEM: TIMESTEP', then skip the section
static int LtrjSkipItemTimestep(FILE *fr, char file[], int *line_count);
// test if next line is 'ITEM: NUMBER OF ATOMS', then read the section
static int LtrjReadNumberOfAtoms(FILE *fr, char file[], int *line_count);
// test if next line is 'ITEM: BOX BOUNDS', then read the section
static int LtrjReadPBCSection(FILE *fr, char file[], BOX *box, int *line_count);
// check if words & split contain 'ITEM: TIMESTEP' line
static bool LtrjCheckTimestepLine();
// check if words & split contain 'ITEM: NUMBER OF ATOMS' line
static bool LtrjCheckNumberAtomsLine();
// check if words & split contain 'ITEM: BOX BOUNDS ...' line
static bool LtrjCheckPbcLine();
// read 'ITEM: ATOMS ...' line, defining what variables are in which columns
static int LtrjReadAtomsLine(FILE *fr, char file[], int n, int *var_pos,
                             char vars[n][10], int *line_count);
// read an atom coordinate line
static int LtrjReadCoorLine(FILE *fr, BEAD *b, int b_count,
                            const int *var, int cols);
// fill a helper array with possible variables in 'ITEM: ATOMS ...' line
static void LtrjFillAtomVariables(int n, char var[n][10]); //}}}
/*
 * Functions to read lammps data file as a structure file (LmpDataReadStruct) or
 * a coordinate file (LmpDataReadTimestep)
 */ //{{{
static int LmpDataReadTimestep(FILE *fr, char file[], SYSTEM *System,
                               int *line_count);
static SYSTEM LmpDataReadStruct(char file[]);
// Helper functions for lmpdata file
// read header of the lammps data file
static int LmpDataReadHeader(FILE *fr, char file[], SYSTEM *System,
                             int *line_count);
// read body of the lammps data file
static void LmpDataReadBody(FILE *fr, char file[], SYSTEM *System,
                            int atom_types, int *line_count);
// functions to read the various sections in the lammps data file
static void LmpDataReadMasses(FILE *fr, char file[], BEADTYPE name_mass[],
                              int lmp_types, int *line_count);
// TODO: non-harmonic bonds
static void LmpDataReadBondCoeffs(FILE *fr, char file[], SYSTEM *System,
                                  int *line_count);
// TODO: non-harmonic angles
static void LmpDataReadAngleCoeffs(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count);
// TODO: non-harmonic (lammps type) angles
static void LmpDataReadDihedralCoeffs(FILE *fr, char file[], SYSTEM *System,
                                      int *line_count);
// TODO: non-cvff (lammps type) impropers
static void LmpDataReadImproperCoeffs(FILE *fr, char file[], SYSTEM *System,
                                      int *line_count);
static void LmpDataReadAtoms(FILE *fr, char file[], SYSTEM *System,
                             BEADTYPE name_mass[], int atom_types,
                             int *line_count, int mode);
static void LmpDataReadVelocities(FILE *fr, char file[], SYSTEM *System,
                                  int *line_count);
static void LmpDataReadBonds(FILE *fr, char file[], COUNT Count, int (*bond)[3],
                             int *line_count);
static void LmpDataReadAngles(FILE *fr, char file[], COUNT Count,
                              int (*angle)[4], int *line_count);
static void LmpDataReadDihedrals(FILE *fr, char file[], COUNT Count,
                                 int (*dihedral)[5], int *line_count);
static void LmpDataReadImpropers(FILE *fr, char file[], COUNT Count,
                                 int (*improper)[5], int *line_count);
//}}}
/*
 * Functions to read vsf/vtf file as a structure file via VtfReadStruct() and
 * vcf/vtf file as a coordinate file via VtfReadTimestep() and VtfSkipTimestep()
 */ //{{{
static SYSTEM VtfReadStruct(char file[], bool detailed);
static int VtfReadTimestep(FILE *fr, char file[],
                           SYSTEM *System, int *line_count);
static int VtfSkipTimestep(FILE *fr, char file[], char vsf_file[],
                           int *line_count);
// Helper functions for vtf/vsf/vcf files
// identify type of a line
static int VtfCheckLineType(char file[], int line_count);
// functions to check if a line is of a specific type
static int VtfCheckCoorOrderedLine(double coor[3]);
static int VtfCheckCoorIndexedLine(double coor[3], long *index);
static int VtfCheckCoordinateLine(double coor[3], long *index);
static int VtfCheckTimestepLine();
static int VtfCheckPbcLine();
static bool VtfCheckAtomLine();
static bool VtfCheckBondLine();
// Find position of keywords in an atom line
static int *VtfAtomLineValues();
// Get the first pbc line from a vcf/vsf/vtf coordinate file.
static BOX VtfReadPBC(char file[]);
// process pbc line from a vtf file
static bool VtfPbcLine(BOX *box, int ltype);
// read timestep preamble in a vtf coordinate file
static int VtfReadTimestepPreamble(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count);
// read variable-size indexed coordinate block in a vtf file
static int VtfReadCoorBlockIndexed(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count);
// read ordered coordinate block in a vtf file
static int VtfReadCoorBlockOrdered(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count);
// read the first timestep to count beads in constant-size coordinate block
static int VtfReadNumberOfBeads(char file[]);
//}}}
/*
 * Function to read dl_meso FIELD-like file as a structure file
 */ //{{{
static SYSTEM FieldRead(char file[]);
// Helper functions for FIELD-like file
// read Species section
static void FieldReadSpecies(char file[], SYSTEM *System);
// read Molecules section
static void FieldReadMolecules(char file[], SYSTEM *System);
//}}}
/*
 * Functions to read xyz file as a coordinate file via XyzReadTimestep() and
 * XyzSkipTimestep() and as a structure file via XyzReadStruct()
 */ //{{{
static int XyzReadTimestep(FILE *fr, char file[], SYSTEM *System,
                           int *line_count);
static bool XyzSkipTimestep(FILE *fr, char file[], int *line_count);
static SYSTEM XyzReadStruct(char file[]);
// Helper functions for xyz file
static bool XyzCheckCoorLine(double coor[3]);
//}}}
/*
 * General helper functions
 */ //{{{
/*
 * Copy .Bead array from MoleculeType to Molecule. It assumes the number of
 * molecules is the same as the number of molecule types, i.e., the function is
 * to be used before molecule types are properly identified.
 */
static void CopyMoleculeTypeBeadsToMoleculeBeads(SYSTEM *System);
// fill MoleculeType[].Bond arrays with data from a separate bond array
static void FillMoleculeTypeBonds(SYSTEM *System, int (*bond)[3], int n);
// fill MoleculeType[].Angle arrays with data from a separate angle array
static void FillMoleculeTypeAngles(SYSTEM *System, int (*angle)[4], int n);
// fill MoleculeType[].Dihedral arrays with data from a separate dihedral array
static void FillMoleculeTypeDihedral(SYSTEM *System, int (*dihedral)[5], int n);
// fill MoleculeType[].Improper arrays with data from a separate improper array
static void FillMoleculeTypeImproper(SYSTEM *System, int (*improper)[5], int n);
// remove molecule/bead types with 0 molecules/beads
static void RemoveExtraTypes(SYSTEM *System);
//}}}

// STATIC IMPLEMENTATIONS
// lammpstrj //{{{
// Use the first lammpstrj timestep as a definition of system composition //{{{
static SYSTEM LtrjReadStruct(char file[]) {
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
  int max_vars = 11, position[max_vars];
  char var[max_vars][10];
  int cols = LtrjReadAtomsLine(fr, file, max_vars, position, var, &line_count);
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
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = line.Position[dd];
      b->Velocity[dd] = line.Velocity[dd];
      b->Force[dd] = line.Force[dd];
    }
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
      if (position[1] == -1 || strcmp(split[position[1]], bt->Name) == 0) {
        bt->Number++;
        b->Type = j;
        new = false;
        break;
      }
    }
    if (new) { // create a new type
      int type = Count->BeadType;
      if (position[1] != -1) { // 'element' variable is present
        NewBeadType(&Sys.BeadType, &Count->BeadType, split[position[1]], CHARGE,
                    MASS, RADIUS);
      } else { // 'element' variable is missing
        NewBeadType(&Sys.BeadType, &Count->BeadType, "b0", CHARGE, MASS,
                    RADIUS);
      }
      BEADTYPE *bt_new = &Sys.BeadType[type];
      bt_new->Number = 1;
      b->Type = type;
    }
  } //}}}
  fclose(fr);
  FillSystemNonessentials(&Sys);
  for (int i = 0; i < Count->Bead; i++) {
    Sys.UnbondedCoor[i] = i;
  }
  // AllocFillBeadTypeIndex(&Sys);
  CheckSystem(Sys, file);
  ChangeBoxByLow(&Sys, -1);
  return Sys;
} //}}}
// Read a single timestep from lammpstrj file //{{{
static int LtrjReadTimestep(FILE *fr, char file[], SYSTEM *System,
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
  int max_vars = 11, position[max_vars];
  char vars[max_vars][10];
  int cols = LtrjReadAtomsLine(fr, file, max_vars, position, vars, line_count);
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
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = line.Position[dd];
      b->Velocity[dd] = line.Velocity[dd];
      b->Force[dd] = line.Force[dd];
    }
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
static int LtrjSkipTimestep(FILE *fr, char file[], int *line_count) {
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
static BOX LtrjReadPBC(char file[]) {
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
static int LtrjReadTimestepPreamble(FILE *fr, char file[], BOX *box,
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
static int LtrjSkipItemTimestep(FILE *fr, char file[], int *line_count) { //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    // ErrorEOF(file); // proper eof - before the first line of a timestep
    return -2;
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
static int LtrjReadNumberOfAtoms(FILE *fr, char file[], int *line_count) { //{{{
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
static int LtrjReadPBCSection(FILE *fr, char file[], BOX *box,
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
    double bounds[2][3], tilt[3];
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
static int LtrjReadAtomsLine(FILE *fr, char file[], int max_vars, int *var_pos,
                             char var[max_vars][10], int *line_count) {
  // check for the correct maximum number of entries //{{{
  if (max_vars != 11) {
    err_msg("CODING: there should be at most 11 entries in 'ITEM: ATOMS' line");
    PrintError();
    exit(1);
  } //}}}
  // generate array with possible variable names
  LtrjFillAtomVariables(max_vars, var);
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
  InitIntArray(var_pos, max_vars, -1); // id, x, y, z, vx, vy, vz, fx, fy, fz
  int cols = -1;
  for (int i = 2; i < words; i++) {
    for (int j = 0; j < max_vars; j++) {
      if (strcmp(split[i], var[j]) == 0) {
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
static void LtrjFillAtomVariables(int n, char var[n][10]) { //{{{
  // check for the correct maximum number of entries //{{{
  if (n != 11) {
    err_msg("CODING: there should be at most 11 entries in ITEM: ATOMS line");
    PrintError();
    exit(1);
  } //}}}
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
} //}}}
  //}}}
// lmpdata //{{{
static SYSTEM LmpDataReadStruct(char file[]) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  FILE *fr = OpenFile(file, "r");
  int line_count = 0;
  int lmp_types = LmpDataReadHeader(fr, file, &System, &line_count);
  LmpDataReadBody(fr, file, &System, lmp_types, &line_count);
  fclose(fr);
  CalculateBoxData(&System.Box, 1);
  RemoveExtraTypes(&System);
  MergeBeadTypes(&System, true);
  MergeMoleculeTypes(&System);
  FillSystemNonessentials(&System);
  // take the data file as a coordinate file as well
  Count->BeadCoor = Count->Bead;
  Count->UnbondedCoor = Count->Unbonded;
  Count->BondedCoor = Count->Bonded;
  System.BeadCoor = s_realloc(System.BeadCoor,
                              Count->Bead * sizeof *System.BeadCoor);
  int c_unbonded = 0, c_bonded = 0;
  for (int i = 0; i < Count->Bead; i++) {
    System.BeadCoor[i] = i;
    if (System.Bead[i].Molecule == -1) {
      System.UnbondedCoor[c_unbonded] = i;
      c_unbonded++;
    } else {
      System.BondedCoor[c_bonded] = i;
      c_bonded++;
    }
  }
  // make molecule indices go from 0 (just to avoid large numbers) //{{{
  // 1) find the lowest molecule index
  int min_id = 1e7;
  for (int i = 0; i < System.Count.Molecule; i++) {
    if (System.Molecule[i].Index < min_id) {
      min_id = System.Molecule[i].Index;
    }
  }
  // 2) subtract the lowest molecule index from molecule indices
  for (int i = 0; i < System.Count.Molecule; i++) {
    System.Molecule[i].Index -= min_id;
  } //}}}
  if (Count->Molecule) {
    System.MoleculeCoor = s_realloc(System.MoleculeCoor, Count->Molecule *
                                    sizeof *System.MoleculeCoor);
  }
  Count->MoleculeCoor = Count->Molecule;
  for (int i = 0; i < Count->MoleculeCoor; i++) {
    System.MoleculeCoor[i] = i;
  }
  CheckSystem(System, file);
  return System;
} //}}}
// TODO: make Atoms section type detection (i.e., # <type>) into a function
// TODO: Atoms & Velocities section can be switched
// TODO: add Forces section reading
// read lammps data file as a coordinate file //{{{
static int LmpDataReadTimestep(FILE *fr, char file[], SYSTEM *System,
                               int *line_count) {
  COUNT *Count = &System->Count;
  // set 'in timestep' to all beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    System->Bead[i].InTimestep = true;
  } //}}}
  // set 'in timestep' to all molecules //{{{
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].InTimestep = false;
  } //}}}
  Count->MoleculeCoor = Count->Molecule;
  for (int i = 0; i < Count->MoleculeCoor; i++) {
    System->MoleculeCoor[i] = i;
  }

  // ignore the first line (comment) //{{{
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "\0");
    exit(1);
  } //}}}
  *line_count = 1;
  // read numer of atoms & box size //{{{
  // read until the first capital-letter-starting line
  fpos_t position;
  do {
    (*line_count)++;
    fgetpos(fr, &position);
    // read a line
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete lammps data file header");
      exit(1);
    }
    long val;
    // <int> atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      if (Count->Bead != val) {
        err_msg("coordinate lammps data file must contain all beads");
        PrintErrorFile(file, "\0", "\0");
        exit(1);
      } //}}}
      // <double> <double> xlo xhi //{{{
    } else if (words > 3 && strcmp(split[2], "xlo") == 0 &&
               strcmp(split[3], "xhi") == 0) {
      double xlo, xhi;
      if (!IsRealNumber(split[0], &xlo) || !IsRealNumber(split[1], &xhi)) {
        goto error;
      }
      System->Box.Low[0] = xlo;
      System->Box.OrthoLength[0] = xhi - xlo; //}}}
      // <double> <double> ylo yhi //{{{
    } else if (words > 3 && strcmp(split[2], "ylo") == 0 &&
               strcmp(split[3], "yhi") == 0) {
      double ylo, yhi;
      if (!IsRealNumber(split[0], &ylo) || !IsRealNumber(split[1], &yhi)) {
        goto error;
      }
      System->Box.Low[1] = ylo;
      System->Box.OrthoLength[1] = yhi - ylo; //}}}
      // <double> <double> zlo zhi //{{{
    } else if (words > 3 && strcmp(split[2], "zlo") == 0 &&
               strcmp(split[3], "zhi") == 0) {
      double zlo, zhi;
      if (!IsRealNumber(split[0], &zlo) || !IsRealNumber(split[1], &zhi)) {
        goto error;
      }
      System->Box.Low[2] = zlo;
      System->Box.OrthoLength[2] = zhi - zlo; //}}}
      // <double> <double> <double> xy xz yz //{{{
    } else if (words > 5 && strcmp(split[3], "xy") == 0 &&
               strcmp(split[4], "xz") == 0 && strcmp(split[5], "yz") == 0) {
      double xy, xz, yz;
      if (!IsRealNumber(split[0], &xy) || !IsRealNumber(split[1], &xz) ||
          !IsRealNumber(split[2], &yz)) {
        goto error;
      }
      System->Box.transform[0][1] = xy;
      System->Box.transform[0][2] = xz;
      System->Box.transform[1][2] = yz;
    } //}}}
  } while (words == 0 || split[0][0] < 'A' || split[0][0] > 'Z');
  //  return file pointer to before the first capital-letter-starting line
  fsetpos(fr, &position);
  (*line_count)--;
  CalculateBoxData(&System->Box, 1); //}}}
  // find 'Atoms' section (and skip the next blank line) //{{{
  do {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
  } while (words == 0 || strcmp(split[0], "Atoms") != 0);
  int mode = 0;
  if (words > 2 && split[1][0] == '#') {
    if (strcmp(split[2], "full") == 0) {
      mode = 0; // for 'Atoms # bond|angle|molecular|charge'
    } else if (strcmp(split[2], "bond") == 0 ||
               strcmp(split[2], "angle") == 0 ||
               strcmp(split[2], "molecular") == 0 ||
               strcmp(split[2], "charge") == 0) {
      mode = 1; // for 'Atoms # bond|angle|molecular'
    } else if (strcmp(split[2], "atomic") == 0) {
      mode = 2; // for 'Atoms # atomic'
    } else {    // for 'Atoms # <something else>'
      snprintf(ERROR_MSG, LINE, "unsupported atom_style (Atoms section) "
               "%s%s; assuming 'full'", ErrYellow(), split[2]);
      PrintWarnFile(file, "\0", "\0");
    }
  }
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  } //}}}
  // read atom lines //{{{
  Count->BeadCoor = Count->Bead;
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    long id;
    double pos[3];
    if (mode == 0) {
      // 'Atoms # full' or just 'Atoms':
      // <bead id> <mol id> <bead type id> <charge> <coordinates>
      if (words < 7 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsRealNumber(split[4], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[5], &pos[1]) || //
          !IsRealNumber(split[6], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    } else if (mode == 1 || mode == 3) {
      // 'Atoms # bond|angle|molecular':
      // <bead id> <mol id> <bead type id> <coordinates>
      // or 'Atoms # charge':
      // <bead id> <bead type id> <charge> <coordinates>
      if (words < 6 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsRealNumber(split[3], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[4], &pos[1]) || //
          !IsRealNumber(split[5], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    } else { // mode == 2
      // 'Atoms # atomic': <bead id> <bead type id> <coordinates>
      if (words < 5 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsRealNumber(split[2], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[3], &pos[1]) || //
          !IsRealNumber(split[4], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    }
    id--; // in lammps data file, ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = pos[dd] - System->Box.Low[dd];
    }
    System->BeadCoor[i] = id;
  } //}}}
  // find 'Velocities' section (and skip the next blank line) //{{{
  do {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ChangeBoxByLow(System, -1);
      FillInCoor(System);
      return 1; // Velocities section is not mandatory
    }
  } while (words == 0 || strcmp(split[0], "Velocities") != 0);
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ChangeBoxByLow(System, -1);
    FillInCoor(System);
    return -2;
  } //}}}
  // read velocity lines //{{{
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    long id;
    double vel[3];
    // error - wrong line //{{{
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        id > Count->Bead ||                // bead index
        !IsRealNumber(split[1], &vel[0]) || // bead velocities
        !IsRealNumber(split[2], &vel[1]) || //
        !IsRealNumber(split[3], &vel[2])) { //
      err_msg("wrong line in Velocities section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }     //}}}
    id--; // in lammps data file, ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Velocity[dd] = vel[dd];
    }
  } //}}}
  ChangeBoxByLow(System, -1);
  FillInCoor(System);
  return 1;
  error: // unrecognised line //{{{
    err_msg("unrecognised line in the file header");
    PrintErrorFileLine(file, *line_count);
    exit(1); //}}}
} //}}}
// read header //{{{
static int LmpDataReadHeader(FILE *fr, char file[], SYSTEM *System,
                             int *line_count) {
  COUNT *Count = &System->Count;
  // ignore the first line (comment) //{{{
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "\0");
    exit(1);
  } //}}}
  *line_count = 1;
  // read until a line starting with capital letter //{{{
  double atom_types = 0;
  fpos_t position;
  do {
    (*line_count)++;
    fgetpos(fr, &position);
    // read a line
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete lammps data file header");
      exit(1);
    }
    // evaluate the line
    long val;
    // <int> atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      Count->Bead = val;
      System->Bead = s_realloc(System->Bead,
                               Count->Bead * sizeof *System->Bead);
      /* Count->Bead = Count->BeadType as each bead can have different charge,
       * so at first, each bead will be its own type
       */
      Count->BeadType = val;
      System->BeadType = s_realloc(System->BeadType,
                                   Count->BeadType * sizeof *System->BeadType);
      for (int i = 0; i < Count->Bead; i++) {
        InitBead(&System->Bead[i]);
        InitBeadType(&System->BeadType[i]);
      } //}}}
    // <int> bonds //{{{
    } else if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Bond = val; //}}}
    // <int> angles //{{{
    } else if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Angle = val; //}}}
    // <int> dihedrals //{{{
    } else if (words > 1 && strcmp(split[1], "dihedrals") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Dihedral = val; //}}}
    // <int> impropers //{{{
    } else if (words > 1 && strcmp(split[1], "impropers") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Improper = val; //}}}
    // <int> atom types //{{{
    } else if (words > 2 && strcmp(split[1], "atom") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      atom_types = val; //}}}
    // <int> bond types //{{{
    } else if (words > 2 && strcmp(split[1], "bond") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (strcmp(split[0], "???") != 0) {
        if (!IsWholeNumber(split[0], &val)) {
          goto error;
        }
        Count->BondType = val;
        if (Count->BondType > 0) {
          System->BondType = s_realloc(System->BondType, Count->BondType *
                                       sizeof *System->BondType);
          for (int i = 0; i < Count->BondType; i++) {
            System->BondType[i] = InitParams;
          }
        }
      }
      //}}}
    // <int> angle types //{{{
    } else if (words > 2 && strcmp(split[1], "angle") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->AngleType = val;
      if (Count->AngleType > 0) {
        System->AngleType = s_realloc(System->AngleType, Count->AngleType *
                                      sizeof *System->AngleType);
        for (int i = 0; i < Count->AngleType; i++) {
          System->AngleType[i] = InitParams;
        }
      }
      //}}}
    // <int> dihedral types //{{{
    } else if (words > 2 && strcmp(split[1], "dihedral") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->DihedralType = val;
      if (Count->DihedralType > 0) {
        System->DihedralType = s_realloc(System->DihedralType,
                                         Count->DihedralType *
                                         sizeof *System->DihedralType);
        for (int i = 0; i < Count->DihedralType; i++) {
          System->DihedralType[i] = InitParams;
        }
      }
      //}}}
    // <int> improper types //{{{
    } else if (words > 2 && strcmp(split[1], "improper") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->ImproperType = val;
      if (Count->ImproperType > 0) {
        System->ImproperType = s_realloc(System->ImproperType,
                                         Count->ImproperType *
                                         sizeof *System->ImproperType);
        for (int i = 0; i < Count->ImproperType; i++) {
          System->ImproperType[i] = InitParams;
        }
      }
      //}}}
    // <double> <double> xlo xhi //{{{
    } else if (words > 3 && strcmp(split[2], "xlo") == 0 &&
               strcmp(split[3], "xhi") == 0) {
      double xlo, xhi;
      if (!IsRealNumber(split[0], &xlo) || !IsRealNumber(split[1], &xhi)) {
        goto error;
      }
      System->Box.Low[0] = xlo;
      System->Box.OrthoLength[0] = xhi - xlo; //}}}
    // <double> <double> ylo yhi //{{{
    } else if (words > 3 && strcmp(split[2], "ylo") == 0 &&
               strcmp(split[3], "yhi") == 0) {
      double ylo, yhi;
      if (!IsRealNumber(split[0], &ylo) || !IsRealNumber(split[1], &yhi)) {
        goto error;
      }
      System->Box.Low[1] = ylo;
      System->Box.OrthoLength[1] = yhi - ylo; //}}}
    // <double> <double> zlo zhi //{{{
    } else if (words > 3 && strcmp(split[2], "zlo") == 0 &&
               strcmp(split[3], "zhi") == 0) {
      double zlo, zhi;
      if (!IsRealNumber(split[0], &zlo) || !IsRealNumber(split[1], &zhi)) {
        goto error;
      }
      System->Box.Low[2] = zlo;
      System->Box.OrthoLength[2] = zhi - zlo; //}}}
    // <double> <double> <double> xy xz yz //{{{
    } else if (words > 5 && strcmp(split[3], "xy") == 0 &&
               strcmp(split[4], "xz") == 0 && strcmp(split[5], "yz") == 0) {
      double xy, xz, yz;
      if (!IsRealNumber(split[0], &xy) || !IsRealNumber(split[1], &xz) ||
          !IsRealNumber(split[2], &yz)) {
        goto error;
      }
      System->Box.transform[0][1] = xy;
      System->Box.transform[0][2] = xz;
      System->Box.transform[1][2] = yz;
    }                                                             //}}}
  } while (words == 0 || split[0][0] < 'A' || split[0][0] > 'Z'); //}}}
  //  return file pointer to before the first capital-letter-starting line
  fsetpos(fr, &position);
  (*line_count)--;
  // check all necessary parts were read //{{{
  if (Count->Bead == 0 || Count->BeadType == 0) {
    err_msg("missing atom or atom types in the file header");
    PrintErrorFile(file, "\0", "\0");
    putc('\n', stderr);
    exit(1);
  }
  if (System->Box.OrthoLength[0] == -1 ||
      System->Box.OrthoLength[1] == -1 ||
      System->Box.OrthoLength[2] == -1) {
    err_msg("missing box size in the file header");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
  }
  //}}}
  return atom_types;
error: // unrecognised line //{{{
  err_msg("unrecognised line in the file header");
  PrintErrorFileLine(file, *line_count);
  exit(1); //}}}
} //}}}
// read body //{{{
static void LmpDataReadBody(FILE *fr, char file[], SYSTEM *System,
                            int atom_types, int *line_count) {
  COUNT *Count = &System->Count;
  BEADTYPE *name_mass = calloc(atom_types, sizeof *name_mass);
  // create arrays for bonds/angles/dihedrals/impropers //{{{
  int(*bond)[3], (*angle)[4], (*dihedral)[5], (*improper)[5];
  if (Count->Bond > 0) {
    bond = calloc(Count->Bond, sizeof *bond);
  } else {
    bond = calloc(1, sizeof *bond);
  }
  if (Count->Angle > 0) {
    angle = calloc(Count->Angle, sizeof *angle);
  } else {
    angle = calloc(1, sizeof *angle);
  }
  if (Count->Dihedral > 0) {
    dihedral = calloc(Count->Dihedral, sizeof *dihedral);
  } else {
    dihedral = calloc(1, sizeof *dihedral);
  }
  if (Count->Improper > 0) {
    improper = calloc(Count->Improper, sizeof *improper);
  } else {
    improper = calloc(1, sizeof *improper);
  } //}}}
  // variables to check all sections are present
  bool atoms = false,                // Atoms section
       bonds[2] = {false, false},     // Bond Coeffs & Bonds sections
       angles[2] = {false, false},    // Angle Coeffs & Angles sections
       dihedrals[2] = {false, false}, // Dihedral Coeffs & Dihedrals sections
       impropers[2] = {false, false}; // Improper Coeffs & Impropers sections
  // read the rest of the file //{{{
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    (*line_count)++;
    // evaluate the line and read the appropriate section
    if (words > 0 && strcmp(split[0], "Masses") == 0) {
      LmpDataReadMasses(fr, file, name_mass, atom_types, line_count);
    } else if (words > 1 && strcmp(split[0], "Bond") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      bonds[0] = true;
      LmpDataReadBondCoeffs(fr, file, System, line_count);
    } else if (words > 1 && strcmp(split[0], "Angle") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      angles[0] = true;
      LmpDataReadAngleCoeffs(fr, file, System, line_count);
    } else if (words > 1 && strcmp(split[0], "Dihedral") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      dihedrals[0] = true;
      LmpDataReadDihedralCoeffs(fr, file, System, line_count);
    } else if (words > 1 && strcmp(split[0], "Improper") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      impropers[0] = true;
      LmpDataReadImproperCoeffs(fr, file, System, line_count);
    } else if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      atoms = true;
      int mode = 0; // for 'Atoms' or 'Atoms # full'
      if (words > 2 && split[1][0] == '#') {
        if (strcmp(split[2], "full") == 0) {
          mode = 0; // for 'Atoms # bond|angle|molecular'
        } else if (strcmp(split[2], "bond") == 0 ||
                   strcmp(split[2], "angle") == 0 ||
                   strcmp(split[2], "molecular") == 0) {
          mode = 1; // for 'Atoms # bond|angle|molecular'
        } else if (strcmp(split[2], "atomic") == 0) {
          mode = 2; // for 'Atoms # atomic'
        } else if (strcmp(split[2], "charge") == 0) {
          mode = 3; // for 'Atoms # charge'
        } else {    // for 'Atoms # <something else>'
          snprintf(ERROR_MSG, LINE, "unsupported atom_style (Atoms section) "
                   "%s%s; assuming 'full'", ErrYellow(), split[2]);
          PrintWarnFile(file, "\0", "\0");
        }
      }
      LmpDataReadAtoms(fr, file, System, name_mass, atom_types, line_count,
                       mode);
    } else if (words > 0 && strcmp(split[0], "Velocities") == 0) {
      LmpDataReadVelocities(fr, file, System, line_count);
    } else if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      bonds[1] = true;
      LmpDataReadBonds(fr, file, *Count, bond, line_count);
    } else if (words > 0 && strcmp(split[0], "Angles") == 0) {
      angles[1] = true;
      LmpDataReadAngles(fr, file, *Count, angle, line_count);
    } else if (words > 0 && strcmp(split[0], "Dihedrals") == 0) {
      dihedrals[1] = true;
      LmpDataReadDihedrals(fr, file, *Count, dihedral, line_count);
    } else if (words > 0 && strcmp(split[0], "Impropers") == 0) {
      impropers[1] = true;
      LmpDataReadImpropers(fr, file, *Count, improper, line_count);
    }
  } //}}}
  free(name_mass);
  // errors/warnings - missing sections //{{{
  if (!atoms) {
    err_msg("Missing Atoms section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Bond > 0 && !bonds[0]) {
    err_msg("Missing Bond Coeffs; setting all bond parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->BondType; i++) {
      System->BondType[i].a = 0;
      System->BondType[i].b = 0;
      System->BondType[i].c = 0;
    }
  }
  if (Count->Bond > 0 && !bonds[1]) {
    err_msg("Missing Bonds section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Angle > 0 && !angles[0]) {
    err_msg("Missing Angle Coeffs; setting all angle parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->AngleType; i++) {
      System->AngleType[i].a = 0;
      System->AngleType[i].b = 0;
      System->AngleType[i].c = 0;
    }
  }
  if (Count->Angle > 0 && !angles[1]) {
    err_msg("Missing Angles section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Dihedral > 0 && !dihedrals[0]) {
    err_msg("Missing Dihedral Coeffs; setting all dihedral parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->DihedralType; i++) {
      System->DihedralType[i].a = 0;
      System->DihedralType[i].b = 0;
      System->DihedralType[i].c = 0;
    }
  }
  if (Count->Dihedral > 0 && !dihedrals[1]) {
    err_msg("Missing Dihedrals section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Improper > 0 && !impropers[0]) {
    err_msg("Missing Improper Coeffs; setting all improper parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->ImproperType; i++) {
      System->ImproperType[i].a = 0;
      System->ImproperType[i].b = 0;
      System->ImproperType[i].c = 0;
    }
  }
  if (Count->Improper > 0 && !impropers[1]) {
    err_msg("Missing Impropers section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  } //}}}
  Count->MoleculeType = Count->Molecule;
  /*
   * Sort beads ascendingly in each MoleculeType (i.e., in each molecule);
   * necessary because bonds (etc.) are added according to their ascending bead
   * indexes, and no connction to the MoleculeType[].Bead array is made
   * TODO: I should rather fill Molecule[].Bead array and only later fill
   *       MoleculeType[].Bead based on connectivity of the molecule in
   *       question. Somehow...
   */
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    SortArrayInt(mt->Bead, mt->nBeads, 0);
  }
  CopyMoleculeTypeBeadsToMoleculeBeads(System);
  FillMoleculeTypeBonds(System, bond, Count->Bond);
  FillMoleculeTypeAngles(System, angle, Count->Angle);
  FillMoleculeTypeDihedral(System, dihedral, Count->Dihedral);
  FillMoleculeTypeImproper(System, improper, Count->Improper);
  free(bond);
  free(angle);
  free(dihedral);
  free(improper);
} //}}}
// read Masses section //{{{
static void LmpDataReadMasses(FILE *fr, char file[], BEADTYPE name_mass[],
                              int atom_types, int *line_count) {
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Masses' section");
    exit(1);
  }
  for (int i = 0; i < atom_types; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Masses' section");
      exit(1);
    }
    long type;
    double mass;
    if (words < 2 || !IsNaturalNumber(split[0], &type) || type > atom_types ||
        (!IsPosRealNumber(split[1], &mass) && strcmp(split[1], "???") != 0)) {
      err_msg("wrong line in Masses section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    type--;
    if (strcmp(split[1], "???") == 0) {
      err_msg("undefined mass (\'???\' in Masses section)");
      PrintWarnFileLine(file, *line_count);
      name_mass[type].Mass = MASS;
    } else {
      name_mass[type].Mass = mass;
    }
    // check for bead type name: # <name>
    if (words > 3 && split[2][0] == '#') {
      strncpy(name_mass[type].Name, split[3], BEAD_NAME);
      name_mass[type].Name[BEAD_NAME-1] = '\0'; // ensure null-termination
    } else {
      snprintf(name_mass[type].Name, BEAD_NAME, "b%d", i);
    }
  }
} //}}}
// read Bond Coeffs section //{{{
// TODO: non-harmonic bonds
static void LmpDataReadBondCoeffs(FILE *fr, char file[], SYSTEM *System,
                                  int *line_count) {
  COUNT *Count = &System->Count;
  // error - no bond types //{{{
  if (Count->BondType == 0) {
    err_msg("Bond Coeffs in a file with no bond types");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    return;
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Bond Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->BondType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Bond Coeffs' section");
      exit(1);
    }
    long type;
    double a, b;
    // error - wrong line //{{{
    if (words < 3 ||
        !IsNaturalNumber(split[0], &type) ||
        type > Count->BondType ||
        !IsRealNumber(split[1], &a) || a < 0 ||
        !IsRealNumber(split[2], &b) || b < 0) {
      err_msg("wrong line in Bond Coeffs section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }                                     //}}}
    System->BondType[type - 1].a = a * 2; // lammps uses k/2 as spring constant
    System->BondType[type - 1].b = b;
  }
} //}}}
// read Angle Coeffs section //{{{
// TODO: non-harmonic angles
static void LmpDataReadAngleCoeffs(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count) {
  COUNT *Count = &System->Count;
  // error - no angle types //{{{
  if (Count->AngleType == 0) {
    err_msg("Angle Coeffs in a file with no angle types");
    PrintWarnFile(file, "\0", "\0");
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Angle Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->AngleType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Angle Coeffs' section");
      exit(1);
    }
    long type;
    double a, b;
    // error - wrong line //{{{
    if (words < 3 ||
        !IsNaturalNumber(split[0], &type) ||
        type > Count->AngleType ||
        !IsRealNumber(split[1], &a) || a < 0 ||
        !IsRealNumber(split[2], &b) || b < 0) {
      err_msg("wrong line in Angle Coeffs section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }                                      //}}}
    System->AngleType[type - 1].a = a * 2; // lammps uses k/2 as spring constant
    System->AngleType[type - 1].b = b;
  }
} //}}}
// read Dihedral Coeffs section //{{{
// TODO: non-harmonic (lammps type) angles
static void LmpDataReadDihedralCoeffs(FILE *fr, char file[], SYSTEM *System,
                                      int *line_count) {
  COUNT *Count = &System->Count;
  // error - no dihedral types //{{{
  if (Count->DihedralType == 0) {
    err_msg("Dihedral Coeffs in a file with no dihedral types");
    PrintWarnFile(file, "\0", "\0");
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Dihedral Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->DihedralType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Dihedral Coeffs' section");
      exit(1);
    }
    long type;
    double a;
    long b, c;
    // error - wrong line //{{{
    // see https://docs.lammps.org/dihedral_harmonic.html for details
    if (words < 4 || !IsNaturalNumber(split[0], &type) ||
        type > Count->DihedralType ||
        !IsRealNumber(split[1], &a) ||
        (!IsIntegerNumber(split[2], &b) && labs(b) != 1) ||
        !IsWholeNumber(split[3], &c)) {
      err_msg("wrong line in Dihedral Coeffs section "
                        "(for now, only harmonic type accepted)");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    System->DihedralType[type-1].a = a * 2;
    System->DihedralType[type-1].b = b;
    System->DihedralType[type-1].c = c;
  }
} //}}}
// read Improper Coeffs section //{{{
// TODO: non-cvff (lammps type) impropers
static void LmpDataReadImproperCoeffs(FILE *fr, char file[], SYSTEM *System,
                                      int *line_count) {
  COUNT *Count = &System->Count;
  // error - no improper types //{{{
  if (Count->ImproperType == 0) {
    err_msg("Improper Coeffs in a file with no improper types");
    PrintWarnFile(file, "\0", "\0");
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Improper Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->ImproperType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Improper Coeffs' section");
      exit(1);
    }
    long type;
    double a;
    long b, c;
    // error - wrong line //{{{
    // see https://docs.lammps.org/improper_cvff.html for details
    if (words < 4 || !IsNaturalNumber(split[0], &type) ||
        type > Count->DihedralType ||
        !IsRealNumber(split[1], &a) ||
        (!IsIntegerNumber(split[2], &b) && labs(b) != 1) ||
        !IsWholeNumber(split[3], &c)) {
      err_msg("wrong line in Improper Coeffs (only cvff type accepted)");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    System->ImproperType[type - 1].a = a * 2;
    System->ImproperType[type - 1].b = b;
    System->ImproperType[type - 1].c = c;
  }
} //}}}
// read Atoms section //{{{
// mode: 0..full; 1..angle/bond; 2..atomic; 3..molecular; 4..charge
static void LmpDataReadAtoms(FILE *fr, char file[], SYSTEM *System,
                             BEADTYPE name_mass[], int atom_types,
                             int *line_count, int mode) {
  COUNT *Count = &System->Count;
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Atoms' section");
    exit(1);
  }
  bool warned = false; // to warn of undefined charge only once
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Atoms' section");
      exit(1);
    }
    long id, resid, type;
    double q;
    double pos[3];
    // check the line is correct & read/assign values //{{{
    // 'Atoms # full': <bead id> <mol id> <bead type id> <charge> <coordinates>
    if (mode == 0) {
      if (words < 7 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                   // bead index
          !IsIntegerNumber(split[1], &resid) || // molecule index
          !IsWholeNumber(split[2], &type) ||    // bead type (mass-defined)
          (!IsRealNumber(split[3], &q) &&       // bead charge
           strcmp(split[3], "???") != 0) ||     //
          !IsRealNumber(split[4], &pos[0]) ||    // Cartesean coordinates
          !IsRealNumber(split[5], &pos[1]) ||    //
          !IsRealNumber(split[6], &pos[2])) {    //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    } else if (mode == 1) {
      // 'Atoms # bond|angle|molecular':
      // <bead id> <mol id> <bead type id> <coordinates>
      if (words < 6 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                   // bead index
          !IsIntegerNumber(split[1], &resid) || // molecule index
          !IsWholeNumber(split[2], &type) ||    // bead type (mass-defined)
          !IsRealNumber(split[3], &pos[0]) ||    // Cartesean coordinates
          !IsRealNumber(split[4], &pos[1]) ||    //
          !IsRealNumber(split[5], &pos[2])) {    //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
      q = CHARGE;
    } else if (mode == 2) {
      // 'Atoms # atomic': <bead id> <bead type id> <coordinates>
      if (words < 5 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsWholeNumber(split[1], &type) || // bead type (mass-defined)
          !IsRealNumber(split[2], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[3], &pos[1]) || //
          !IsRealNumber(split[4], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
      resid = -1;
      q = CHARGE;
    } else {
      // 'Atoms # charge': <bead id> <bead type id> <charge> <coordinates>
      if (words < 6 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsWholeNumber(split[1], &type) || // bead type (mass-defined)
          (!IsRealNumber(split[2], &q) &&    // bead charge
           strcmp(split[3], "???") != 0) ||  //
          !IsRealNumber(split[3], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[4], &pos[1]) || //
          !IsRealNumber(split[5], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
      resid = -1;
    }
    //}}}
    id--;   // in lammps data file, bead ids start from 1
    type--; // in lammps data file, bead type ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = pos[dd];
    }
    b->InTimestep = true;
    b->Type = id;
    /*
     * bead type is same as bead index but with lmp-defined mass, i.e., there
     * is the same number of beads and bead types; the single-bead bead types
     * will be merged later
     */
    BEADTYPE *bt = &System->BeadType[id];
    if (type > atom_types) {
      snprintf(ERROR_MSG, LINE, "bead type is too high; id is %s%ld%s in "
               "a file with %s%d%s atom types", ErrRed(), type, ErrYellow(),
               ErrRed(), atom_types, ErrYellow());
      PrintErrorFile(file, "\0", "\0");
      exit(1);
    }
    bt->Mass = name_mass[type].Mass;
    if (strcmp(split[3], "???") == 0) {
      if (!warned) {
        err_msg("undefined charge (\'???\' in Atoms section)");
        PrintWarnFile(file, "\0", "\0");
        putc('\n', stderr);
        warned = true;
      }
      bt->Charge = CHARGE;
    } else {
      bt->Charge = q;
    }
    bt->Number = 1;
    s_strcpy(bt->Name, name_mass[type].Name, BEAD_NAME);
    if (resid >= 0) { // bead in a molecule //{{{
      Count->Bonded++;
      b->Molecule = resid;
      // resid ids may be discontinuous, so define Molecule for all possible ids
      if (resid > Count->HighestResid) {
        System->MoleculeType = s_realloc(System->MoleculeType, (resid + 1) *
                                         sizeof *System->MoleculeType);
        System->Molecule = s_realloc(System->Molecule, (resid + 1) *
                                     sizeof *System->Molecule);
        for (int j = (Count->HighestResid + 1); j <= resid; j++) {
          InitMoleculeType(&System->MoleculeType[j]);
          InitMolecule(&System->Molecule[j]);
        }
        Count->HighestResid = resid; // goes from 0
        Count->Molecule = resid + 1;
      }
      // does the molecule type (i.e., that molecule id) exist?
      MOLECULETYPE *mt_resid = &System->MoleculeType[resid];
      if (mt_resid->Number == 0) { // no, create new type
        char name[MOL_NAME] = "m";
        for (int j = 0; j < words; j++) {
          if (split[j][0] == '#' && words > j) {
            snprintf(name, MOL_NAME, "%s", split[j+1]);
            mt_resid->Flag = true;
            break;
          }
        }
        s_strcpy(mt_resid->Name, name, MOL_NAME);
        mt_resid->Number = 1;
        mt_resid->nBeads = 1;
        mt_resid->Bead = malloc(sizeof *mt_resid->Bead);
        mt_resid->Bead[0] = id;
      } else { // yes, assign the atom to that type
        int bead = mt_resid->nBeads;
        mt_resid->nBeads++;
        mt_resid->Bead = s_realloc(mt_resid->Bead, mt_resid->nBeads *
                                   sizeof *mt_resid->Bead);
        mt_resid->Bead[bead] = id;
      }
    } else {
      Count->Unbonded++;
    } //}}}
  }
} //}}}
// read Velocities section //{{{
static void LmpDataReadVelocities(FILE *fr, char file[], SYSTEM *System,
                                  int *line_count) {
  COUNT *Count = &System->Count;
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Velocities' section");
    exit(1);
  }
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Velocities' section");
      exit(1);
    }
    long id;
    double vel[3];
    // error - wrong line //{{{
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        id > Count->Bead ||                // bead index
        !IsRealNumber(split[1], &vel[0]) || // bead velocities
        !IsRealNumber(split[2], &vel[1]) || //
        !IsRealNumber(split[3], &vel[2])) { //
      err_msg("wrong line in Velocities section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    id--;
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Velocity[dd] = vel[dd];
    }
  }
} //}}}
// read Bonds section //{{{
static void LmpDataReadBonds(FILE *fr, char file[], COUNT Count, int (*bond)[3],
                             int *line_count) {
  bool *found = calloc(Count.Bond, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Bonds' section");
    exit(1);
  }
  bool warned = false;
  for (int i = 0; i < Count.Bond; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Bonds' section");
      exit(1);
    }
    long id, type, b_id[2];
    // errors //{{{
    // not <int> <int> <int> <int> line
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        (!IsNaturalNumber(split[1], &type) && strcmp(split[1], "???") != 0) ||
        !IsNaturalNumber(split[2], &b_id[0]) ||
        !IsNaturalNumber(split[3], &b_id[1])) {
      err_msg("wrong line in Bonds section");
      goto error;
    }
    if (id > Count.Bond) {
      err_msg("Bond index is too high");
      goto error;
    }
    if (type > Count.BondType) {
      err_msg("Bond type is too high");
      goto error;
    }
    if (b_id[0] > Count.Bead || b_id[1] > Count.Bead) {
      err_msg("Bead index in a bond is too high");
      goto error;
    }
    if (b_id[0] == b_id[1]) {
      err_msg("Same beads in a bond");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated bond index");
      goto error;
    } //}}}
    bond[id - 1][0] = b_id[0] - 1;
    bond[id - 1][1] = b_id[1] - 1;
    if (strcmp(split[1], "???") == 0) {
      if (!warned) {
        err_msg("undefined bond type (\'???\' in Bonds section)");
        PrintWarnFile(file, "\0", "\0");
        warned = true;
      }
      bond[id - 1][2] = -1;
    } else {
      bond[id - 1][2] = type - 1;
    }
    found[id - 1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
// read Angles section //{{{
static void LmpDataReadAngles(FILE *fr, char file[], COUNT Count,
                              int (*angle)[4], int *line_count) {
  bool *found = calloc(Count.Angle, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Angles' section");
    exit(1);
  }
  for (int i = 0; i < Count.Angle; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Angles' section");
      exit(1);
    }
    long id, type, a_id[3];
    // errors //{{{
    // not <int> <int> <int> <int> line
    if (words < 5 || !IsNaturalNumber(split[0], &id) ||
        (!IsNaturalNumber(split[1], &type) && strcmp(split[1], "???") != 0) ||
        !IsNaturalNumber(split[2], &a_id[0]) ||
        !IsNaturalNumber(split[3], &a_id[1]) ||
        !IsNaturalNumber(split[4], &a_id[2])) { //
      err_msg("wrong line in Angles section");
      goto error;
    }
    if (id > Count.Angle) {
      err_msg("Angle index is too high");
      goto error;
    }
    if (type > Count.AngleType) {
      err_msg("Angle type is too high");
      goto error;
    }
    if (a_id[0] > Count.Bead || a_id[1] > Count.Bead || a_id[2] > Count.Bead) {
      err_msg("Bead index in an angle is too high");
      goto error;
    }
    if (a_id[0] == a_id[1] || a_id[0] == a_id[2] || a_id[1] == a_id[2]) {
      err_msg("Same beads in an angle");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated angle index");
      goto error;
    } //}}}
    angle[id - 1][0] = a_id[0] - 1;
    angle[id - 1][1] = a_id[1] - 1;
    angle[id - 1][2] = a_id[2] - 1;
    angle[id - 1][3] = type - 1;
    found[id - 1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
// read Dihedrals section //{{{
static void LmpDataReadDihedrals(FILE *fr, char file[], COUNT Count,
                                 int (*dihedral)[5], int *line_count) {
  bool *found = calloc(Count.Dihedral, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Dihedrals' section");
    exit(1);
  }
  for (int i = 0; i < Count.Dihedral; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Dihedrals' section");
      exit(1);
    }
    long id, type, d_id[4];
    // errors //{{{
    // not <int> <int> <int> <int> line
    if (words < 6 || !IsNaturalNumber(split[0], &id) ||
        !IsNaturalNumber(split[1], &type) ||
        !IsNaturalNumber(split[2], &d_id[0]) ||
        !IsNaturalNumber(split[3], &d_id[1]) ||
        !IsNaturalNumber(split[4], &d_id[2]) ||
        !IsNaturalNumber(split[5], &d_id[3])) {
      err_msg("wrong line in Dihedrals section");
      goto error;
    }
    if (id > Count.Dihedral) {
      err_msg("Dihedral index is too high");
      goto error;
    }
    if (type > Count.DihedralType) {
      err_msg("Dihedral type is too high");
      goto error;
    }
    if (d_id[0] > Count.Bead ||
        d_id[1] > Count.Bead ||
        d_id[2] > Count.Bead ||
        d_id[3] > Count.Bead) {
      err_msg("Bead index in a dihedral is too high");
      goto error;
    }
    if (d_id[0] == d_id[1] || d_id[0] == d_id[2] || d_id[0] == d_id[3] ||
        d_id[1] == d_id[2] || d_id[1] == d_id[3] || d_id[2] == d_id[3]) {
      err_msg("Same beads in a single dihedral");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated dihedral index");
      goto error;
    } //}}}
    dihedral[id - 1][0] = d_id[0] - 1;
    dihedral[id - 1][1] = d_id[1] - 1;
    dihedral[id - 1][2] = d_id[2] - 1;
    dihedral[id - 1][3] = d_id[3] - 1;
    dihedral[id - 1][4] = type - 1;
    found[id - 1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
// read Impropers section //{{{
static void LmpDataReadImpropers(FILE *fr, char file[], COUNT Count,
                                 int (*improper)[5], int *line_count) {
  bool *found = calloc(Count.Angle, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Impropers' section");
    exit(1);
  }
  for (int i = 0; i < Count.Improper; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Impropers' section");
      exit(1);
    }
    long id, type, i_id[4];
    // errors //{{{
    // not 6x<int>
    if (words < 6 || !IsNaturalNumber(split[0], &id) ||
        !IsNaturalNumber(split[1], &type) ||
        !IsNaturalNumber(split[2], &i_id[0]) ||
        !IsNaturalNumber(split[3], &i_id[1]) ||
        !IsNaturalNumber(split[4], &i_id[2]) ||
        !IsNaturalNumber(split[5], &i_id[3])) {
      err_msg("wrong line in Impropers section");
      goto error;
    }
    if (id > Count.Improper) {
      err_msg("Improper index is too high");
      goto error;
    }
    if (type > Count.ImproperType) {
      err_msg("Improper type is too high");
      goto error;
    }
    if (i_id[0] > Count.Bead || i_id[1] > Count.Bead || i_id[2] > Count.Bead ||
        i_id[3] > Count.Bead) {
      err_msg("Bead index in a improper is too high");
      goto error;
    }
    if (i_id[0] == i_id[1] || i_id[0] == i_id[2] || i_id[0] == i_id[3] ||
        i_id[1] == i_id[2] || i_id[1] == i_id[3] || i_id[2] == i_id[3]) {
      err_msg("Same beads in a single improper");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated improper index");
      goto error;
    } //}}}
    improper[id - 1][0] = i_id[0] - 1;
    improper[id - 1][1] = i_id[1] - 1;
    improper[id - 1][2] = i_id[2] - 1;
    improper[id - 1][3] = i_id[3] - 1;
    improper[id - 1][4] = type - 1;
    found[id - 1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
  //}}}
// vtf //{{{
// VtfReadStruct() //{{{
/*
 * Read system information from vsf/vtf structure file. It can recognize bead
 * and molecule types based either on name only or on all information (name,
 * mass, charge, and radius for bead types; bead order, bonds, angles, and
 * dihedrals for molecule types).
 */
static SYSTEM VtfReadStruct(char file[], bool detailed) {
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  FILE *fr = OpenFile(file, "r");
  // define variables and structures and arrays //{{{
  int line_count = 0,   // total number of lines
      count_atoms = 0,  // number of 'atom <id>'
      default_atom = 0, // line number of the first 'atom default' line
      count_bonds = 0;  // number of bonds
  bool warned = false; // has 'a[tom] default' line warning been already issued?
  //}}}
  // read struct_file line by line, counting atoms and bonds //{{{
  /*
   * Test all the lines are correct, counting them and printing possible
   * warnings. For atom lines, take the highest <id>+1 as the number of atoms
   * (ids start from 0), save line number of the first 'atom default' line
   * (if present), and count molecules (well, get highest resid). For bond
   * lines, just count them. Ignore pbc, blank, or comment lines, stop reading
   * on timestep line, and exit with error for any other line (including a valid
   * coordinate line as they cannot be inside a structure block)
   */
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    line_count++;
    // read line
    int ltype = VtfCheckLineType(file, line_count);
    if (ltype == ATOM_LINE) { //{{{
      if (strcmp(split[1], "default") == 0) { // 'a[tom] default' line
        if (default_atom != 0 && !warned) { // warn of multiple defaults //{{{
          warned = true; // warn only once
          snprintf(ERROR_MSG, LINE, "multiple 'a[tom] default' lines %s, using"
                   " line %s%d%s as the default line%s\n", ErrCyan(),
                   ErrYellow(), default_atom, ErrCyan(), ErrColourReset());
          PrintWarnFile(file, "\0", "\0"); //}}}
        } else { // save line number of the first 'atom default' line
          default_atom = line_count;
        }
      } else { // 'a[tom] <id>' line
        int id = atoi(split[1]);
        count_atoms++;
        // highest bead index? (corresponds to the number of beads in vsf)
        if (id >= Count->Bead) {
          Count->Bead = id + 1; // +1 as bead ids start from 0 in vsf
        }
        // save values from the 'a[tom] <id>' line
        int *value = VtfAtomLineValues();
        // is the bead in a molecule?
        if (value[5] > -1) {
          int resid = atoi(split[value[5]]);
          if (resid > Count->HighestResid) {
            Count->HighestResid = resid;
          }
        }
      } //}}}
    } else if (ltype == BOND_LINE) {
      count_bonds++;
    } else if (ltype == TIME_LINE || ltype == TIME_LINE_O) {
      break;
    } else if (ltype == COOR_LINE || ltype == COOR_LINE_O) {
      err_msg("encountered a coordinate-like line inside the structure block ");
      PrintErrorFileLine(file, line_count);
      exit(1);
    } else if (ltype != BLANK_LINE && ltype != COMMENT_LINE &&
               ltype != PBC_LINE && ltype != PBC_LINE_ANGLES) { // e)
      // proper error message already established in VtfCheckLineType()
      PrintErrorFileLine(file, line_count);
      exit(1);
    }
  }
  fclose(fr); //}}}
  // error - no default line and too few atom lines //{{{
  if (default_atom == 0 && count_atoms != Count->Bead) {
    int undefined = Count->Bead - count_atoms;
    snprintf(ERROR_MSG, LINE, "not all beads defined ('atom default' line is "
             "omitted); %s%d%s bead(s) undefined", ErrYellow(), undefined,
             ErrRed());
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  } //}}}
  // assume there are as many molecules as resids (which start from 0, hence +1)
  Count->Molecule = Count->HighestResid + 1;
  // allocate & initialize structures //{{{
  Sys.BeadType = s_realloc(Sys.BeadType, sizeof *Sys.BeadType * Count->Bead);
  Sys.Bead = s_realloc(Sys.Bead, sizeof *Sys.Bead * Count->Bead);
  if (Count->Molecule > 0) {
    Sys.MoleculeType = s_realloc(Sys.MoleculeType,
                                 sizeof *Sys.MoleculeType * Count->Molecule);
    Sys.Molecule = s_realloc(Sys.Molecule, sizeof *Sys.Molecule *
                             Count->Molecule);
    Sys.MoleculeCoor = s_realloc(Sys.MoleculeCoor,
                                 Count->Molecule * sizeof *Sys.MoleculeCoor);
  }
  for (int i = 0; i < Count->Molecule; i++) {
    InitMoleculeType(&Sys.MoleculeType[i]);
    InitMolecule(&Sys.Molecule[i]);
  }
  for (int i = 0; i < Count->Bead; i++) {
    InitBeadType(&Sys.BeadType[i]);
    InitBead(&Sys.Bead[i]);
  } //}}}
  // array to save bonds - if there are any //{{{
  int (*bond)[3] = NULL;
  if (count_bonds != 0) {
    bond = calloc(count_bonds, sizeof *bond);
    bond[0][2] = -1; // no bond types in a vtf file
  } //}}}
  // default type - may not be used
  BEADTYPE bt_def;
  InitBeadType(&bt_def);

  fr = OpenFile(file, "r");
  // go through the file again to save atom and bond info //{{{
  /*
   * a)tom line: save bead information
   *   i) save 'atom default' line info into separate bt_def struct
   *   ii) atom <id> line: save into separate Sys.BeadType & Sys.Bead
   *   iii) resid keyword: save into separate Sys.MoleculeType & Sys.Molecule
   * b)ond line: save bonded bead indices into a bond struct for later use
   */
  count_bonds = 0;
  for (int line = 1; line <= line_count; line++) {
    // read line
    ReadAndSplitLine(fr, SPL_STR, " \t\n");
    if (words > 1 && split[0][0] == 'a') { // a)
      if (line == default_atom) { // 'a[tom] default' line //{{{
        // save values for the default bead type
        int *value = VtfAtomLineValues();
        if (value[0] != -1) { // TODO: should be always true (for now)
          snprintf(bt_def.Name, BEAD_NAME, "%s", split[value[0]]);
          bt_def.Name[BEAD_NAME-1] = '\0'; // ensure null-termination
        }
        if (value[1] != -1) {
          bt_def.Mass = atof(split[value[1]]);
        }
        if (value[2] != -1) {
          bt_def.Charge = atof(split[value[2]]);
        }
        if (value[3] != -1) {
          bt_def.Radius = atof(split[value[3]]);
        }
        //}}}
      } else if (split[1][0] != 'd') { // 'a[tom] <id>' line //{{{
        int id = atoi(split[1]);
        // warning - repeated atom line //{{{
        if (id < Count->Bead && Sys.BeadType[id].Number != 0) {
          count_atoms--; // this counted lines, not checked ids
          err_msg("atom defined multiple times; ignoring this line");
          PrintWarnFileLine(file, line);
          continue; // go to reading the next line
        } //}}}
        // save values from the 'a[tom] <id>' line
        int *value = VtfAtomLineValues();
        if (value[0] != -1) { // TODO: should never happen (for now)
          snprintf(Sys.BeadType[id].Name, BEAD_NAME, "%s", split[value[0]]);
          Sys.BeadType[id].Name[BEAD_NAME-1] = '\0'; // ensure null-termination
        }
        if (value[1] != -1) {
          Sys.BeadType[id].Mass = atof(split[value[1]]);
        }
        if (value[2] != -1) {
          Sys.BeadType[id].Charge = atof(split[value[2]]);
        }
        if (value[3] != -1) {
          Sys.BeadType[id].Radius = atof(split[value[3]]);
        }
        Sys.BeadType[id].Number = 1;
        Sys.Bead[id].Type = id;
        // is the bead in a molecule?
        if (value[5] > -1) {
          int resid = atoi(split[value[5]]);
          MOLECULETYPE *mt_resid = &Sys.MoleculeType[resid];
          if (mt_resid->Number == 0) { // new molecule type
            if (value[4] == -1) {
              s_strcpy(mt_resid->Name, "m", MOL_NAME);
              mt_resid->Flag = false;
            } else {
              strncpy(mt_resid->Name, split[value[4]], MOL_NAME);
              mt_resid->Name[MOL_NAME-1] = '\0';
              mt_resid->Flag = true;
            }
            mt_resid->Name[MOL_NAME-1] = '\0'; // null-terminate!
            mt_resid->Number = 1;
            mt_resid->nBeads = 1;
            mt_resid->Bead = malloc(sizeof *mt_resid->Bead);
            mt_resid->Bead[0] = id; // bead type = bead index
          } else {                  // not new moleclue type
            int bead = mt_resid->nBeads;
            mt_resid->nBeads++;
            mt_resid->Bead = s_realloc(mt_resid->Bead, sizeof *mt_resid->Bead *
                                       mt_resid->nBeads);
            mt_resid->Bead[bead] = id; // bead type = bead index
            if (value[4] != -1 && strcmp(mt_resid->Name, "m") == 0) {
              strncpy(mt_resid->Name, split[value[4]], MOL_NAME);
              mt_resid->Name[MOL_NAME-1] = '\0';
              mt_resid->Flag = true;
            }
          }
          Sys.Bead[id].Molecule = resid;
        }
      } //}}}
    } else if (words > 1 && split[0][0] == 'b') { // b)
      long val; //{{{
      if (words == 2) { // case 'bond <id>:<id>'
        char *index[SPL_STR];
        SplitLine(SPL_STR, index, split[1], ":");
        IsIntegerNumber(index[0], &val);
        bond[count_bonds][0] = val;
        IsIntegerNumber(index[1], &val);
        bond[count_bonds][1] = val;
      } else { // case 'bond <id>: <id>'
        // IsIntegerNumber() ignores the trailing ':'
        IsIntegerNumber(split[1], &val);
        bond[count_bonds][0] = val;
        IsIntegerNumber(split[2], &val);
        bond[count_bonds][1] = val;
      }
      bond[count_bonds][2] = -1;
      // assure index1<index2 (may be unnecessary, but definitely won't hurt)
      if (bond[count_bonds][0] > bond[count_bonds][1]) {
        SwapInt(&bond[count_bonds][0], &bond[count_bonds][1]);
      }
      count_bonds++; //}}}
    }
  } //}}}
  fclose(fr);
  Sys.BeadCoor = s_realloc(Sys.BeadCoor, Count->Bead * sizeof *Sys.BeadCoor);
  Count->BeadType = Count->Bead;
  Count->MoleculeType = Count->Molecule;
  // assign atom default to default beads & count bonded/unbonded beads //{{{
  // find first unused bead type and make it the default
  int def = -1;
  for (int i = 0; i < Count->Bead; i++) {
    if (Sys.BeadType[i].Number == 0) {
      def = i;
      Sys.BeadType[def] = bt_def;
      Sys.BeadType[def].Number = Count->Bead - count_atoms;
      Sys.Bead[def].Type = def;
      break;
    }
  }
  // add the default beads to their proper type & count Bonded/Unbonded
  for (int i = 0; i < Count->Bead; i++) {
    if (Sys.BeadType[i].Number == 0) { // default bead?
      Sys.Bead[i].Type = def;
    }
    if (Sys.Bead[i].Molecule == -1) { // is 'i' in a molecule?
      Count->Unbonded++;              // unbonded bead
    } else {
      Count->Bonded++; // bonded bead
    }
  }
  // just check that it counts the beads correctly
  if (Count->Bead != (Count->Unbonded + Count->Bonded)) {
    err_msg("wrong with bead counting; should never happen!");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  } //}}}
  CopyMoleculeTypeBeadsToMoleculeBeads(&Sys);
  FillMoleculeTypeBonds(&Sys, bond, count_bonds);
  if (count_bonds != 0) {
    free(bond);
  }
  RemoveExtraTypes(&Sys);
  MergeBeadTypes(&Sys, detailed);
  MergeMoleculeTypes(&Sys);
  FillSystemNonessentials(&Sys);
  CheckSystem(Sys, file);
  Sys.Box = VtfReadPBC(file);
  return Sys;
} //}}}
// VtfReadTimestep() //{{{
static int VtfReadTimestep(FILE *fr, char file[],
                           SYSTEM *System, int *line_count) {
  int timestep = VtfReadTimestepPreamble(fr, file, System, line_count);
  if (timestep < 0) {
    return timestep;
  }
  // set 'not in timestep' to all beads //{{{
  for (int i = 0; i < System->Count.Bead; i++) {
    System->Bead[i].InTimestep = false;
  } //}}}
  // set 'not in timestep' to all molecules //{{{
  for (int i = 0; i < System->Count.Molecule; i++) {
    System->Molecule[i].InTimestep = false;
  }
  //}}}
  // read coordinates
  COUNT *Count = &System->Count;
  Count->UnbondedCoor = 0;
  Count->BondedCoor = 0;
  if (timestep == TIME_LINE) {
    int test = VtfReadCoorBlockIndexed(fr, file, System, line_count);
    if (test < 0) {
      return test;
    }
  } else {
    int test = VtfReadCoorBlockOrdered(fr, file, System, line_count);
    if (test < 0) {
      return test;
    }
  }
  return 1; // coordinates read properly
} //}}}
// VtfSkipTimestep() //{{{
/*
 * Discard a single timestep from a vcf/vtf coordinate file. It assumes the
 * coordinates are ordered lines (i.e., three numbers); if a timestep line is
 * missing, it prints only a warning (and skips the ensuing coordinate lines).
 * Returns false only if a line cannot be read during preamble (e.g., on eof).
 */
static int VtfSkipTimestep(FILE *fr, char file[], char vsf_file[],
                           int *line_count) {
  // skip preamble - i.e., read until the first number-starting line
  fpos_t position;
  do {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
  } while (words == 0 || split[0][0] < '0' || split[0][0] > '9');
  // skip coordinate section - i.e., read until first non-number-starting line
  do {
    fgetpos(fr, &position);
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return 1;
    }
  } while (words > 0 && split[0][0] >= '0' && split[0][0] <= '9');
  (*line_count)--;
  fsetpos(fr, &position); // return to before that non-number-starting line
  return 1;
} //}}}
static int VtfCheckLineType(char file[], int line_count) { //{{{
  ERROR_MSG[0] = '\0'; // clear error message array
  // blank line
  if (words == 0) {
    return BLANK_LINE;
  }
  // comment line
  if (split[0][0] == '#') {
    return COMMENT_LINE;
  }
  // coordinate line
  double trash[3];
  long trash2;
  int test = VtfCheckCoordinateLine(trash, &trash2);
  if (test != ERROR_LINE) {
    return test;
  }
  // timestep line
  test = VtfCheckTimestepLine();
  if (test != ERROR_LINE) {
    return test;
  }
  // pbc line
  test = VtfCheckPbcLine();
  if (test != ERROR_LINE) {
    return test;
  }
  // atom line (vsf)
  if (VtfCheckAtomLine()) {
    return ATOM_LINE;
  }
  // bond line (vsf)
  if (VtfCheckBondLine()) {
    return BOND_LINE;
  }
  if (ERROR_MSG[0] == '\0') {
    err_msg("unrecognised line");
  }
  return ERROR_LINE;
} //}}}
// VtfAtomLineValues() //{{{
/*
 * Returns 'array' where subscripts correspond to n-th split[] for:
 * 0..name, 1..mass, 2..charge, 3..radius, 4..resame, 5..resid. If not present,
 * the corresponding element holds - 1.
 */
static int *VtfAtomLineValues() {
  static int value[6];
  for (int i = 0; i < 6; i++) {
    value[i] = -1;
  }
  for (int i = 2; i < words; i += 2) {
    if (split[i][0] == 'n') {
      value[0] = i + 1;
    } else if (split[i][0] == 'm') {
      value[1] = i + 1;
    } else if (split[i][0] == 'q' || strcmp(split[i], "charge") == 0) {
      value[2] = i + 1;
    } else if (split[i][0] == 'r' && strncmp(split[i], "res", 3) != 0) {
      value[3] = i + 1;
    } else if (strncmp(split[i], "resname", 3) == 0 &&
               strcmp(split[i], "resid") != 0) {
      value[4] = i + 1;
    } else if (strcmp(split[i], "resid") == 0) {
      value[5] = i + 1;
    }
  }
  return value;
} //}}}
static int VtfCheckCoorOrderedLine(double coor[3]) { //{{{
  if (words >= 3 &&
      IsRealNumber(split[0], &coor[0]) &&
      IsRealNumber(split[1], &coor[1]) &&
      IsRealNumber(split[2], &coor[2])) {
    return COOR_LINE_O;
  }
  return ERROR_LINE;
} //}}}
static int VtfCheckCoorIndexedLine(double coor[3], long *index) { //{{{
  // indexed line (may also be ordered)
  if (words >= 4 && IsWholeNumber(split[0], index) &&
      IsRealNumber(split[1], &coor[0]) &&
      IsRealNumber(split[2], &coor[1]) &&
      IsRealNumber(split[3], &coor[2])) {
    return COOR_LINE;
  } else {
    return ERROR_LINE;
  }
} //}}}
static int VtfCheckCoordinateLine(double coor[3], long *index) { //{{{
  // indexed line (may also be ordered)
  if (VtfCheckCoorIndexedLine(coor, index) == COOR_LINE) {
    return COOR_LINE;
  }
  // definitely ordered line
  if (VtfCheckCoorOrderedLine(coor) == COOR_LINE_O) {
    return COOR_LINE_O;
  }
  return ERROR_LINE;
} //}}}
static int VtfCheckTimestepLine() { //{{{
  // there are several possibilities how the timestep line can look
  /* ordered timestep:
   *   1) 't[imestep]'
   *   2) 't[imestep] o[rdered] ...'
   *   3) 'o[rdered] ...'
   */
  if ((words == 1 && split[0][0] == 't') ||                      // 1)
      (words > 1 && split[0][0] == 't' && split[1][0] == 'o') || // 2)
      split[0][0] == 'o') {                                      // 3)
    return TIME_LINE_O;
  }
  /* indexed timestep:
   *   1) 't[imestep] i[ndexed] ...'
   *   2) 'i[ndexed] ...'
   */
  if ((words > 1 && split[0][0] == 't' && split[1][0] == 'i') || // 1)
      split[0][0] == 'i') {                                      // 2)
    return TIME_LINE;
  }
  return ERROR_LINE; // not a timestep line
} //}}}
static int VtfCheckPbcLine() { //{{{
  // valid line: pbc <x> <y> <z> [<alpha> <beta> <gamm>]
  // unrecognised line
  double val;
  if (words < 4 || strcmp(split[0], "pbc") != 0 ||
      !IsPosRealNumber(split[1], &val) || !IsPosRealNumber(split[2], &val) ||
      !IsPosRealNumber(split[3], &val)) {
    return ERROR_LINE;
  } else if (words > 6 && IsPosRealNumber(split[4], &val) &&
             IsPosRealNumber(split[5], &val) &&
             IsPosRealNumber(split[6], &val)) {
    return PBC_LINE_ANGLES;
  } else {
    return PBC_LINE;
  }
} //}}}
static bool VtfCheckAtomLine() { //{{{
  long val_i;
  double val_d;
  // error - line not starting with a[tom] default/<id> //{{{
  if (split[0][0] != 'a' || (strcmp(split[1], "default") != 0 &&
                             !IsIntegerNumber(split[1], &val_i))) {
    return false;
  } //}}}
  // error - odd number of strings //{{{
  if ((words % 2) != 0) {
    err_msg("atom line with odd number of strings ");
    return false;
  } //}}}
  // check <keyword> <value> pairs
  bool name = false;
  for (int i = 2; i < words; i += 2) {
    // is n[ame] keyword present?
    if (split[i][0] == 'n') {
      name = true;
    }
    // error - resid not followed by non-negative integer //{{{
    if (strcmp(split[i], "resid") == 0 &&
        !IsWholeNumber(split[i + 1], &val_i)) {
      err_msg("atom line: 'resid' not followed by natural number");
      return false;
    } //}}}
    // error - charge|q //{{{
    if ((strcmp(split[i], "charge") == 0 || split[i][0] == 'q') &&
        !IsRealNumber(split[i + 1], &val_d)) {
      err_msg("atom line: 'charge|q' not followed by real number ");
      return false; //}}}
    // error - r[adius] not followed by positive number //{{{
    } else if (split[i][0] == 'r' &&               // possible r[adius]
               strncmp(split[i], "res", 3) != 0 && // it's not resid or resname
               !IsPosRealNumber(split[i + 1], &val_d)) {
      err_msg("atom line: 'r[adius]]]' not followed by positive real number ");
      return false; //}}}
    // error - m[ass] not followed by positive number //{{{
    } else if (split[i][0] == 'm' && !IsPosRealNumber(split[i + 1], &val_d)) {
      err_msg("atom line: 'm[ass]' not followed by positive real number");
      return false;
    } //}}}
  }
  // error - missing the mandatory n[ame] keyword //{{{
  if (!name) { // TODO: should eventually disappear
    err_msg("atom line: missing 'n[ame]' keyword ");
    return false;
  } //}}}
  return true; // valid atom line
} //}}}
static bool VtfCheckBondLine() { //{{{
  long val_i;
  // valid line: 'b[ond] <id>:[ ]<id> anything'
  // error - only one string or missing 'b[ond]' keyword
  if (words < 2 || split[0][0] != 'b') {
    return false;
  }
  // two strings - assume '<int>:<int>' and test it
  if (words == 2) {
    // ':'-split the second string
    char *index[SPL_STR], string[SPL_LEN];
    s_strcpy(string, split[1], SPL_LEN);
    int strings = SplitLine(2, index, string, ":");
    // int strings = SplitLine_old(index, string, ":");
    if (strings != 2 || !IsIntegerNumber(index[0], &val_i) || val_i < 0 ||
        !IsIntegerNumber(index[1], &val_i) || val_i < 0) {
      err_msg("bond line: only 'b[ond] <int>:<int>' or "
              "'b[ond] <int>: <int>' is valid (for now)");
      return false;
    }
  }
  // more than two strings - assume '<int>: <int>' and test it
  if (words > 2) {
    if (!IsIntegerNumber(split[1], &val_i) || val_i < 0 ||
        !IsIntegerNumber(split[2], &val_i) || val_i < 0) {
      err_msg("bond line: only 'b[ond] <int>:<int>' or "
              "'b[ond] <int>: <int>' is valid (for now)");
      return false;
    }
  }
  return true; // valid bond line
} //}}}
static BOX VtfReadPBC(char file[]) { //{{{
  BOX Box = InitBox;
  FILE *fr = OpenFile(file, "r");
  int line_count = 0;
  // read input file line by line
  while (true) {
    line_count++;
    // read line & split it via whitespace //{{{
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      break;
    } //}}}
    int ltype = VtfCheckLineType(file, line_count);
    // pbc line //{{{
    if (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES) {
      VtfPbcLine(&Box, ltype);
      break; //}}}
      // coordinate line - return from function //{{{
    } else if (ltype == COOR_LINE || ltype == COOR_LINE_O) {
      break;
    } //}}}
  };
  fclose(fr);
  return Box;
} //}}}
static bool VtfPbcLine(BOX *box, int ltype) { //{{{
  box->Length[0] = atof(split[1]);
  box->Length[1] = atof(split[2]);
  box->Length[2] = atof(split[3]);
  if (ltype == PBC_LINE_ANGLES) {
    box->alpha = atof(split[4]);
    box->beta = atof(split[5]);
    box->gamma = atof(split[6]);
  } else { // definitely orthogonal box
    box->alpha = 90;
    box->beta = 90;
    box->gamma = 90;
  }
  if (!CalculateBoxData(box, 0)) {
    return false;
  } else {
    return true;
  }
} //}}}
// read timestep preamble in a vtf file //{{{
static int VtfReadTimestepPreamble(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count) {
  fpos_t position; // to save file position
  int ltype = -1, timestep = -1;
  while (ltype != COOR_LINE && ltype != COOR_LINE_O) {
    // save file pointer position for when it's the first coordinate line
    fgetpos(fr, &position);
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    ltype = VtfCheckLineType(file, *line_count);
    if (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES) {
      BOX box = InitBox;
      if (!VtfPbcLine(&box, ltype)) {
        err_msg("ignoring invalid pbc line");
        PrintErrorFileLine(file, *line_count);
      } else {
        System->Box = box;
      }
    } else if (ltype == TIME_LINE || ltype == TIME_LINE_O) {
      // warn: multiple timestep lines //{{{
      if (timestep != -1) {
        err_msg("extra 'timestep' line in a vtf timestep");
        PrintWarnFileLine(file, *line_count);
      } //}}}
      timestep = ltype;
    } else if (ltype == ERROR_LINE) {
      err_msg("ignoring unrecognised line in a timestep preamble");
      PrintWarnFileLine(file, *line_count);
    }
  }
  // error: missing 'timestep' line //{{{
  if (timestep == -1) {
    err_msg("missing 'timestep' line");
    PrintErrorFile(file, "\0", "\0");
    return -1;
  } //}}}
  // restore file pointer to before the first coordinate line
  fsetpos(fr, &position);
  (*line_count)--; // the first coordinate line will be re-read
  return timestep;
} //}}}
// read variable-size indexed coordinate block in a vtf file //{{{
static int VtfReadCoorBlockIndexed(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count) {
  fpos_t position;
  COUNT *Count = &System->Count;
  Count->BeadCoor = 0;
  // read the first coordinate line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE, "premature end of file; "
             "no beads in an indexed coordinate block");
    PrintErrorFile(file, "\0", "\0");
    return -2;
  }
  double coordinate[3];
  long id;
  while (VtfCheckCoorIndexedLine(coordinate, &id) == COOR_LINE) {
    if (id > Count->Bead) {
      err_msg("bead index is too high");
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
    BEAD *bead_id = &System->Bead[id];
    if (bead_id->InTimestep) {
      err_msg("multiple beads with the same id; ignoring this line");
      PrintWarnFileLine(file, *line_count);
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        break; // 'proper' eof
      }
      continue;
    }
    bead_id->InTimestep = true;
    for (int dd = 0; dd < 3; dd++) {
      bead_id->Position[dd] = coordinate[dd];
    }
    double vel[3];
    if (words >= 7 &&
        IsRealNumber(split[4], &vel[0]) &&
        IsRealNumber(split[5], &vel[1]) &&
        IsRealNumber(split[6], &vel[2])) {
      for (int dd = 0; dd < 3; dd++) {
        bead_id->Velocity[dd] = vel[dd];
      }
    } else {
      bead_id->Velocity[0] = 0;
    }
    System->BeadCoor[Count->BeadCoor] = id;
    Count->BeadCoor++;
    if (Count->BeadCoor > Count->Bead) {
      err_msg("more beads in a variable-sized indexed "
              "coordinate block than in the whole system");
      return -1;
    }
    // read next line
    fgetpos(fr, &position);
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      break; // 'proper' eof
    }
  }
  FillInCoor(System);
  fsetpos(fr, &position);
  return 1;
} //}}}
// read ordered coordinate block in a vtf file //{{{
static int VtfReadCoorBlockOrdered(FILE *fr, char file[], SYSTEM *System,
                                   int *line_count) {
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      snprintf(ERROR_MSG, LINE,
               "premature end of ordered "
               "coordinate block (%s%d%s lines instead of %s%d%s)",
               ErrYellow(), i, ErrRed(), ErrYellow(), Count->Bead, ErrRed());
      PrintErrorFile(file, "\0", "\0");
      return -2;
    }
    double coor[3];
    if (VtfCheckCoorOrderedLine(coor) == ERROR_LINE) {
      snprintf(ERROR_MSG, LINE, "unrecognized line in constant-size coordinate"
               " block (%s%d%s-th line from %s%d%s)",
               ErrYellow(), i + 1, ErrRed(), ErrYellow(), Count->BeadCoor,
               ErrRed());
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
    BEAD *bead_i = &System->Bead[i];
    for (int dd = 0; dd < 3; dd++) {
      bead_i->Position[dd] = coor[dd];
    }
    double vel[3];
    if (words >= 6 &&
        IsRealNumber(split[3], &vel[0]) &&
        IsRealNumber(split[4], &vel[1]) &&
        IsRealNumber(split[5], &vel[2])) {
      for (int dd = 0; dd < 3; dd++) {
        bead_i->Velocity[dd] = vel[dd];
      }
    } else {
      for (int dd = 0; dd < 3; dd++) {
        bead_i->Velocity[dd] = 0;
      }
    }
    System->BeadCoor[i] = i;
  }
  FillInCoor(System);
  return 1;
} //}}}
int VtfReadNumberOfBeads(char file[]) { //{{{
  int ltype = -1, line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // read until vcf timestep line line
  while (ltype != TIME_LINE && ltype != TIME_LINE_O) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    ltype = VtfCheckLineType(file, line_count);
  }
  // test whether next line(s) are pbc
  fpos_t position;
  do {
    line_count++;
    fgetpos(fr, &position);
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    ltype = VtfCheckLineType(file, line_count);
  } while (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES);
  fsetpos(fr, &position);
  line_count--;
  // read until the first non-coordinate line
  int nbeads = -1; // that first non-coordinate line is also counted
  do {
    line_count++;
    nbeads++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      fclose(fr);
      return nbeads;
    }
    ltype = VtfCheckLineType(file, line_count);
  } while (ltype == COOR_LINE || ltype == COOR_LINE_O);
  fclose(fr);
  return nbeads;
} //}}}
  //}}}
// FIELD-like //{{{
static SYSTEM FieldRead(char file[]) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  // read box dimensions, if present //{{{
  FILE *fr = OpenFile(file, "r");
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "empty file");
    exit(1);
  }
  double box[3];
  if (words > 2 && IsRealNumber(split[0], &box[0]) &&
                   IsRealNumber(split[1], &box[1]) &&
                   IsRealNumber(split[2], &box[2])) {
    System.Box.Length[0] = box[0];
    System.Box.Length[1] = box[1];
    System.Box.Length[2] = box[2];
    CalculateBoxData(&System.Box, 0);
  }
  fclose(fr); //}}}
  FieldReadSpecies(file, &System);
  // fill System.Bead & System.Unbonded //{{{
  if (Count->Unbonded > 0) {
    System.Bead = s_realloc(System.Bead, Count->Bead * sizeof *System.Bead);
    System.Unbonded = s_realloc(System.Unbonded,
                                sizeof *System.Unbonded * Count->Unbonded);
    int count = 0;
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = 0; j < System.BeadType[i].Number; j++) {
        InitBead(&System.Bead[count]);
        System.Bead[count].Type = i;
        System.Bead[count].Molecule = -1;
        System.Unbonded[count] = count;
        count++;
      }
    }
    if (count != Count->Unbonded) {
      err_msg("wrong number of unbonded beads; should never happen!");
      PrintErrorFile(file, "\0", "\0");
      exit(1);
    }
  } //}}}
  FieldReadMolecules(file, &System);
  RemoveExtraTypes(&System);
  MergeBeadTypes(&System, true);
  MergeMoleculeTypes(&System);
  if (Count->Bead == 0) {
    err_msg("No beads in the file");
    PrintWarnFile(file, "\0", "\0");
  } else {
    System.BeadCoor = s_realloc(System.BeadCoor,
                                Count->Bead * sizeof *System.BeadCoor);
  }
  if (Count->Molecule > 0) {
    System.MoleculeCoor = s_realloc(System.MoleculeCoor, Count->Molecule *
                                    sizeof *System.MoleculeCoor);
  }
  FillSystemNonessentials(&System);
  CheckSystem(System, file);
  return System;
} //}}}
static void FieldReadSpecies(char file[], SYSTEM *System) { //{{{
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // skip till line starting with 'Species' keyword //{{{
  bool test;
  while ((test = ReadAndSplitLine(fr, SPL_STR, " \t\n"))) {
    line_count++;
    if (words > 0 && strncasecmp(split[0], "species", 6) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    ErrorEOF(file, "missing 'Species' line");
    exit(1);
  } //}}}
  COUNT *Count = &System->Count;
  // read number of bead types //{{{
  long types;
  if (words < 2 || !IsWholeNumber(split[1], &types)) {
    err_msg("incorrect 'Species' keyword line");
    PrintErrorFileLine(file, line_count);
    exit(1);
  } //}}}
  // read bead types //{{{
  bool warned = false; // warn about undefined mass/charge only once
  for (int i = 0; i < types; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Species' section");
      exit(1);
    }
    double mass, charge;
    long number;
    // error - illegal 'species' line //{{{
    if (words < 4 ||
        (!IsPosRealNumber(split[1], &mass) && strcmp(split[1], "???") != 0) ||
        (!IsRealNumber(split[2], &charge) && strcmp(split[2], "???") != 0) ||
        !IsWholeNumber(split[3], &number)) {
      err_msg("incorrect 'Species' line");
      PrintErrorFileLine(file, line_count);
      exit(1);
    } //}}}
    if (strcmp(split[1], "???") == 0) {
      mass = MASS;
      if (!warned) {
        err_msg("undefined charge/mass (\'???\' in species section)");
        PrintWarnFile(file, "\0", "\0");
        putc('\n', stderr);
        warned = true;
      }
    }
    if (strcmp(split[2], "???") == 0) {
      charge = CHARGE;
      if (!warned) {
        err_msg("undefined charge/mass (\'???\' in species section)");
        PrintWarnFile(file, "\0", "\0");
        putc('\n', stderr);
        warned = true;
      }
    }
    NewBeadType(&System->BeadType, &Count->BeadType, split[0], charge, mass,
                HIGHNUM);
    System->BeadType[Count->BeadType-1].Number = number;
    Count->Bead += number;
    Count->Unbonded += number;
  } //}}}
  fclose(fr);
} //}}}
static void FieldReadMolecules(char file[], SYSTEM *System) { //{{{
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // skip till line starting with 'Molecule' keyword //{{{
  bool test;
  while ((test = ReadAndSplitLine(fr, SPL_STR, " \t\n"))) {
    line_count++;
    if (words > 0 && strncasecmp(split[0], "molecules", 7) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    ErrorEOF(file, "missing 'Molecules' line");
    exit(1);
  } //}}}
  // read number of types //{{{
  long val;
  if (words < 2 || !IsWholeNumber(split[1], &val)) {
    err_msg("incorrect 'Molecules' keyword line");
    PrintErrorFileLine(file, line_count);
    exit(1);
  } //}}}
  COUNT *Count = &System->Count;
  Count->MoleculeType = val;
  if (Count->MoleculeType > 0) {
    System->MoleculeType = s_realloc(System->MoleculeType,
                                     sizeof *System->MoleculeType *
                                     Count->MoleculeType);
    // read molecule types
    /*
     * Order of entries:
     *   1) <molecule name>
     *   2) nummol[s] <int>
     *   3) bead[s] <int>
     *      <int> lines: <bead name> <x> <y> <z>
     *   4) bond[s] <int> (optional)
     *      <int> lines: harm <id1> <id2> <k> <r_0>
     *   5) angle[s] <int> (optional)
     *      <int> lines: harm <id1> <id2> <id3> <k> <theta_0>
     *   6) dihed[rals] <int> (optional)
     *      <int> lines: harm <id1> <id2> <id3> <id4> <k>
     *   7) improp[ers] <int> (optional)
     *      <int> lines: harm <id1> <id2> <id3> <id4> <k>
     *   8) finish keyword
     */
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      InitMoleculeType(mt_i);
      // 1) name //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words == 0) {
        err_msg("missing molecule name");
        PrintErrorFileLine(file, line_count);
        exit(1);
      }
      strncpy(mt_i->Name, split[0], MOL_NAME);
      mt_i->Name[MOL_NAME - 1] = '\0'; // null-terminate! //}}}
      // 2) number of molecules //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words < 2 || strcasecmp(split[0], "nummols") != 0 ||
          // !IsNaturalNumber(split[1], &val)) {
          !IsIntegerNumber(split[1], &val) || val < 0) {
        err_msg("incorrect 'nummols' line in a molecule entry");
        PrintErrorFileLine(file, line_count);
        exit(1);
      } else if (val == 0) {
        do {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
        } while(strcmp(split[0], "finish") != 0);
        i--;
        Count->MoleculeType--;
        continue;
      }
      mt_i->Number = val; //}}}
      // 3) beads in the molecule //{{{
      // a) number of beads //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words < 2 || strncasecmp(split[0], "beads", 4) != 0 ||
          !IsNaturalNumber(split[1], &val)) {
        err_msg("incorrect 'beads' line in a molecule entry");
        PrintErrorFileLine(file, line_count);
        exit(1);
      }
      mt_i->nBeads = val; //}}}
      mt_i->Bead = malloc(sizeof *mt_i->Bead * mt_i->nBeads);
      double (*coor)[3] = malloc(sizeof *coor * mt_i->nBeads);
      // b) beads themselves //{{{
      for (int j = 0; j < mt_i->nBeads; j++) {
        line_count++;
        // read a line //{{{
        if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
          ErrorEOF(file, "incomplete 'Molecules' section");
          exit(1);
        } //}}}
        // error - incorrect line //{{{
        if (words < 4 ||
            !IsRealNumber(split[1], &coor[j][0]) ||
            !IsRealNumber(split[2], &coor[j][1]) ||
            !IsRealNumber(split[3], &coor[j][2])) {
          err_msg("incorrect bead coordinate line in a molecule entry");
          PrintErrorFileLine(file, line_count);
          exit(1);
        } //}}}
        int btype = FindBeadType(split[0], *System);
        // error - unknown type //{{{
        if (btype == -1) {
          err_msg("unknown bead in a molecule entry");
          PrintErrorFileLine(file, line_count);
          ErrorBeadType(split[0], *System);
          exit(1);
        } //}}}
        mt_i->Bead[j] = btype;
      }                                //}}}
      int count = Count->Bead,         // for filling Molecule[] & Bead[]
          mol_count = Count->Molecule; //
      Count->Molecule += mt_i->Number;
      Count->Bead += mt_i->Number * mt_i->nBeads;
      Count->Bonded += mt_i->Number * mt_i->nBeads;
      System->Molecule = s_realloc(System->Molecule,
                                   sizeof *System->Molecule * Count->Molecule);
      System->Bead = s_realloc(System->Bead,
                               sizeof *System->Bead * Count->Bead);
      System->Bonded = s_realloc(System->Bonded,
                                 sizeof System->Bonded * Count->Bonded);
      for (int j = count; j < Count->Bead; j++) {
        InitBead(&System->Bead[j]);
      }
      // c) fill Molecule[] & Bead[] //{{{
      for (int j = 0; j < mt_i->Number; j++) {
        MOLECULE *mol = &System->Molecule[mol_count];
        mol->Type = i;
        mol->Index = mol_count;
        mol->Bead = malloc(sizeof *mol->Bead * mt_i->nBeads);
        mol->InTimestep = true;
        for (int k = 0; k < mt_i->nBeads; k++) {
          BEAD *bead = &System->Bead[count];
          bead->Type = mt_i->Bead[k];
          bead->Molecule = mol_count;
          for (int dd = 0; dd < 3; dd++) {
            bead->Position[dd] = coor[k][dd];
          }
          System->BeadType[bead->Type].Number++;
          System->Bonded[count - Count->Unbonded] = count;
          mol->Bead[k] = count;
          count++;
        }
        mol_count++;
      }
      // pro forma checks //{{{
      if (count != Count->Bead) {
        err_msg("wrong number of beads; should never happen!");
        PrintErrorFile(file, "\0", "\0");
        exit(1);
      }
      if (mol_count != Count->Molecule) {
        err_msg("wrong number of molecules; should never happen!");
        PrintErrorFile(file, "\0", "\0");
        exit(1);
      } //}}}
      //}}}
      free(coor); //}}}
      // 4) bonds in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "bonds", 4) == 0) {
        // a) number of bonds //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          err_msg("incorrect 'bonds' line in a molecule entry");
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nBonds = val; //}}}
        mt_i->Bond = malloc(sizeof *mt_i->Bond * mt_i->nBonds);
        // b) bonds themselves & bond types //{{{
        // TODO: for now, only harmonic bonds are considered
        bool warned = false;
        for (int j = 0; j < mt_i->nBonds; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[2];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 3 ||
              !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              beads[0] == beads[1]) {
            err_msg("incorrect bond line in a molecule entry");
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads) {
            err_msg("bead index in a bond is too high");
            PrintErrorFileLine(file, line_count);
            exit(1);
          }                                //}}}
          mt_i->Bond[j][0] = beads[0] - 1; // in FIELD, indices start from 1
          mt_i->Bond[j][1] = beads[1] - 1; //
          // read up to three values for bond type
          for (int k = 3; k <= 5; k++) {
            double *ptr = NULL;
            if (k == 3) {
              ptr = &values.a;
            } else if (k == 4) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined bond parameter (\'???\') "
                         "in molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          }
          // assign bond type //{{{
          int bond_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this bond type already exists
            for (int k = 0; k < Count->BondType; k++) {
              if (System->BondType[k].a == values.a &&
                  System->BondType[k].b == values.b &&
                  System->BondType[k].c == values.c) {
                bond_type = k;
                break;
              }
            }
            // create a new bond type if necessary
            if (bond_type == -1) {
              bond_type = Count->BondType;
              Count->BondType++;
              System->BondType = s_realloc(System->BondType, Count->BondType *
                                           sizeof *System->BondType);
              System->BondType[bond_type] = InitParams;
              System->BondType[bond_type].a = values.a;
              System->BondType[bond_type].b = values.b;
              System->BondType[bond_type].c = values.c;
            }
          } //}}}
          mt_i->Bond[j][2] = bond_type;
        } //}}}
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for "
                 "molecule %s%s", ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
      // 5) angles in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "angles", 4) == 0) {
        // a) number of angles //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          snprintf(ERROR_MSG, LINE, "incorrect line in 'angles' line for "
                   "molecule %s%s", ErrYellow(), mt_i->Name);
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nAngles = val; //}}}
        mt_i->Angle = malloc(sizeof *mt_i->Angle * mt_i->nAngles);
        // b) angles themselves & angle types //{{{
        // TODO: for now, only harmonic angles are considered
        bool warned = false; // to warn of undefined angle type only once
        for (int j = 0; j < mt_i->nAngles; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[3];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 4 || !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              !IsNaturalNumber(split[3], &beads[2]) || beads[0] == beads[1] ||
              beads[0] == beads[2] || beads[1] == beads[2]) {
            snprintf(ERROR_MSG, LINE,
                     "incorrect line in an entry for "
                     "molecule %s%s",
                     ErrYellow(), mt_i->Name);
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
              beads[2] > mt_i->nBeads) {
            snprintf(ERROR_MSG, LINE,
                     "bead index in an angle is too high in "
                     "molecule %s%s",
                     ErrYellow(), mt_i->Name);
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          mt_i->Angle[j][0] = beads[0] - 1; // in FIELD, bead indices go from 1
          mt_i->Angle[j][1] = beads[1] - 1; //
          mt_i->Angle[j][2] = beads[2] - 1; //
          // read up to three values for angle type
          for (int k = 4; k <= 6; k++) {
            double *ptr = NULL;
            if (k == 4) {
              ptr = &values.a;
            } else if (k == 5) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined angle parameter (\'???\') "
                         " for molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          }
          // assign angle type //{{{
          int angle_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this angle type already exists
            for (int k = 0; k < Count->AngleType; k++) {
              if (System->AngleType[k].a == values.a &&
                  System->AngleType[k].b == values.b &&
                  System->AngleType[k].c == values.c) {
                angle_type = k;
                break;
              }
            }
            // create a new angle type if necessary
            if (angle_type == -1) {
              angle_type = Count->AngleType;
              Count->AngleType++;
              System->AngleType = s_realloc(System->AngleType,
                                            Count->AngleType *
                                            sizeof *System->AngleType);
              System->AngleType[angle_type] = InitParams;
              System->AngleType[angle_type].a = values.a;
              System->AngleType[angle_type].b = values.b;
              System->AngleType[angle_type].c = values.c;
            }
          } //}}}
          mt_i->Angle[j][3] = angle_type;
        } //}}}
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE,
                 "unrecognised line in an entry for molecule "
                 "%s%s",
                 ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
      // 6) dihedrals in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "dihedrals", 5) == 0) {
        // a) number of dihedrals //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          snprintf(ERROR_MSG, LINE, "incorrect 'dihedrals' line for molecule "
                   "%s%s", ErrYellow(), mt_i->Name);
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nDihedrals = val; //}}}
        mt_i->Dihedral = malloc(sizeof *mt_i->Dihedral * mt_i->nDihedrals);
        // b) dihedrals themselves & dihedral types //{{{
        // TODO: for now, only harmonic dihedrals are considered
        bool warned = false; // to warn of undefined type only once
        for (int j = 0; j < mt_i->nDihedrals; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[4];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 5 ||
              !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              !IsNaturalNumber(split[3], &beads[2]) ||
              !IsNaturalNumber(split[4], &beads[3]) ||
              beads[0] == beads[1] ||
              beads[0] == beads[2] ||
              beads[0] == beads[3] ||
              beads[1] == beads[2] ||
              beads[1] == beads[3] ||
              beads[2] == beads[3]) {
            snprintf(ERROR_MSG, LINE, "incorrect dihedral line for molecule "
                     "%s%s", ErrYellow(), mt_i->Name);
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
              beads[2] > mt_i->nBeads || beads[3] > mt_i->nBeads) {
            err_msg("bead index in a improper is too high");
            PrintErrorFileLine(file, line_count);
            exit(1);
          }                                    //}}}
          mt_i->Dihedral[j][0] = beads[0] - 1; // in FIELD, indices start from 1
          mt_i->Dihedral[j][1] = beads[1] - 1; //
          mt_i->Dihedral[j][2] = beads[2] - 1; //
          mt_i->Dihedral[j][3] = beads[3] - 1; //
          // read up to three values for dihedral type //{{{
          for (int k = 5; k <= 7; k++) {
            double *ptr = NULL;
            if (k == 5) {
              ptr = &values.a;
            } else if (k == 6) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined dihedral paramter "
                         "(\'???\') in molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          } //}}}
          // assign dihedral type //{{{
          int dihedral_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this dihedral type already exists
            for (int k = 0; k < Count->DihedralType; k++) {
              if (System->DihedralType[k].a == values.a &&
                  System->DihedralType[k].b == values.b &&
                  System->DihedralType[k].c == values.c) {
                dihedral_type = k;
                break;
              }
            }
            // create a new dihedral type if necessary
            if (dihedral_type == -1) {
              dihedral_type = Count->DihedralType;
              Count->DihedralType++;
              System->DihedralType = s_realloc(System->DihedralType,
                                               Count->DihedralType *
                                               sizeof *System->DihedralType);
              System->DihedralType[dihedral_type] = InitParams;
              System->DihedralType[dihedral_type] = values;
            }
          }
          mt_i->Dihedral[j][4] = dihedral_type; //}}}
        }                                       //}}}
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for molecule "
                 "%s%s", ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
      // 7) impropers in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "impropers", 6) == 0) {
        // a) number of impropers //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          snprintf(ERROR_MSG, LINE,
                   "incorrect 'impropers' line in an entry "
                   "for molecule %s%s",
                   ErrYellow(), mt_i->Name);
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nImpropers = val; //}}}
        mt_i->Improper = malloc(sizeof *mt_i->Improper * mt_i->nImpropers);
        // b) impropers themselves & improper types //{{{
        bool warned = false; // to warn of undefined type only once
        for (int j = 0; j < mt_i->nImpropers; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[4];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 5 ||
              !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              !IsNaturalNumber(split[3], &beads[2]) ||
              !IsNaturalNumber(split[4], &beads[3]) ||
              beads[0] == beads[1] ||
              beads[0] == beads[2] ||
              beads[0] == beads[3] ||
              beads[1] == beads[2] ||
              beads[1] == beads[3] ||
              beads[2] == beads[3]) {
            snprintf(ERROR_MSG, LINE, "incorrect improper line for molecule "
                     "%s%s", ErrYellow(), mt_i->Name);
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
              beads[2] > mt_i->nBeads || beads[3] > mt_i->nBeads) {
            err_msg("bead index in a improper is too high");
            PrintErrorFileLine(file, line_count);
            exit(1);
          }                                    //}}}
          mt_i->Improper[j][0] = beads[0] - 1; // in FIELD, indices start from 1
          mt_i->Improper[j][1] = beads[1] - 1; //
          mt_i->Improper[j][2] = beads[2] - 1; //
          mt_i->Improper[j][3] = beads[3] - 1; //
          // read up to three values for improper type //{{{
          for (int k = 5; k <= 7; k++) {
            double *ptr = NULL;
            if (k == 5) {
              ptr = &values.a;
            } else if (k == 6) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined improper parameter "
                         "(\'???\') in molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          } //}}}
          // assign improper type //{{{
          int improper_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this improper type already exists
            for (int k = 0; k < Count->ImproperType; k++) {
              if (System->ImproperType[k].a == values.a &&
                  System->ImproperType[k].b == values.b &&
                  System->ImproperType[k].c == values.c) {
                improper_type = k;
                break;
              }
            }
            // create a new improper type if necessary
            if (improper_type == -1) {
              improper_type = Count->ImproperType;
              Count->ImproperType++;
              System->ImproperType = s_realloc(System->ImproperType,
                                               Count->ImproperType *
                                               sizeof *System->ImproperType);
              System->ImproperType[improper_type] = InitParams;
              System->ImproperType[improper_type] = values;
            }
          }
          mt_i->Improper[j][4] = improper_type; //}}}
        }                                       //}}}
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for molecule "
                 "%s%s", ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
      // finish keyword //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words == 0 || strcasecmp(split[0], "finish") != 0) {
        snprintf(ERROR_MSG, LINE,
                 "missing 'finish' at the end of an entry for "
                 "molecule %s%s",
                 ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
    }
    Count->HighestResid = Count->Molecule;
  }
  fclose(fr);
} //}}}
  //}}}
// xyz //{{{
static SYSTEM XyzReadStruct(char file[]) { //{{{
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // read number of beads //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing number of beads");
    exit(1);
  }
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz file");
    PrintErrorFileLine(file, line_count);
    exit(1);
  } //}}}
  Count->Bead = val;
  Count->BeadCoor = val;
  Count->Unbonded = val;
  Count->UnbondedCoor = val;
  Sys.Bead = s_realloc(Sys.Bead, sizeof *Sys.Bead * Count->Bead);
  Sys.BeadCoor = s_realloc(Sys.BeadCoor, sizeof *Sys.BeadCoor * Count->Bead);
  // read next line //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing comment line");
    exit(1);
  } //}}}
  // read atoms //{{{
  for (int i = 0; i < Count->Bead; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete timestep");
      exit(1);
    }
    double coor[3];
    if (!XyzCheckCoorLine(coor)) {
      err_msg("wrong coordinate line");
      PrintErrorFileLine(file, line_count);
      exit(1);
    }
    BEAD *b = &Sys.Bead[i];
    InitBead(b);
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = coor[dd];
    }
    b->InTimestep = true;
    Sys.BeadCoor[i] = i;
    // assign bead type or create a new one
    bool new = true;
    for (int j = 0; j < Count->BeadType; j++) {
      BEADTYPE *bt = &Sys.BeadType[j];
      if (strcmp(split[0], bt->Name) == 0) {
        bt->Number++;
        b->Type = j;
        new = false;
        break;
      }
    }
    if (new) {
      int type = Count->BeadType;
      NewBeadType(&Sys.BeadType, &Count->BeadType, split[0], HIGHNUM, HIGHNUM, HIGHNUM);
      BEADTYPE *bt_new = &Sys.BeadType[type];
      bt_new->Number = 1;
      b->Type = type;
    }
  } //}}}
  fclose(fr);
  FillSystemNonessentials(&Sys);
  for (int i = 0; i < Count->Bead; i++) {
    Sys.UnbondedCoor[i] = i;
  }
  CheckSystem(Sys, file);
  return Sys;
} //}}}
// XYZReadTimestep() //{{{
static int XyzReadTimestep(FILE *fr, char file[], SYSTEM *System,
                           int *line_count) {
  System->Count.MoleculeCoor = 0;
  // read number of beads //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz timestep");
    PrintWarnFileLine(file, *line_count);
    return -1;
  } else if (val > System->Count.Bead) {
    snprintf(ERROR_MSG, LINE,
             "too many beads in the timestep (maximum number is %s%d%s)",
             ErrYellow(), System->Count.Bead, ErrRed());
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  System->Count.BeadCoor = val; //}}}
  // ignore next line //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing comment line");
    return -2;
  }
  (*line_count)++; //}}}
  // read coordinates //{{{
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      snprintf(ERROR_MSG, LINE, "premature end of file "
               "(%s%d%s coordinate lines instead of %s%d%s)", ErrYellow(),
               i, ErrRed(), ErrYellow(), System->Count.BeadCoor, ErrRed());
      PrintErrorFileLine(file, *line_count);
      return -2;
    }
    double coor[3];
    if (XyzCheckCoorLine(coor)) {
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[i].Position[dd] = coor[dd];
      }
      System->BeadCoor[i] = i;
      (*line_count)++;
      double vel[3];
      if (words > 6 &&
          IsRealNumber(split[4], &vel[0]) &&
          IsRealNumber(split[5], &vel[1]) &&
          IsRealNumber(split[6], &vel[2])) {
        for (int dd = 0; dd < 3; dd++) {
          System->Bead[i].Velocity[dd] = vel[dd];
        }
      } else {
        for (int dd = 0; dd < 3; dd++) {
          System->Bead[i].Velocity[dd] = 0;
        }
      }
    } else {
      err_msg("unrecognized line (not enough valid coordinate lines)");
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
  } //}}}
  FillInCoor(System);
  return true;
} //}}}
static bool XyzSkipTimestep(FILE *fr, char file[], int *line_count) { //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz timestep");
    PrintWarnFileLine(file, *line_count);
    return -1;
  }
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  for (int i = 0; i < val; i++) {
    (*line_count)++;
    char a;
    while ((a = getc(fr) != '\n') && a != EOF)
      ;
  }
  return true;
} //}}}
static bool XyzCheckCoorLine(double coor[3]) { //{{{
  if (words > 3 &&
      IsRealNumber(split[1], &coor[0]) &&
      IsRealNumber(split[2], &coor[1]) &&
      IsRealNumber(split[3], &coor[2])) {
    return true;
  } else {
    return false;
  }
} //}}}
  //}}}
// General helper functions //{{{
static void CopyMoleculeTypeBeadsToMoleculeBeads(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    MOLECULE *mol_i = &System->Molecule[i];
    if (mt_i->nBeads == 1 && !mt_i->Flag) { // remove 'fake' molecules
      mt_i->Number = 0;
      System->Bead[mt_i->Bead[0]].Molecule = -1;
      free(mt_i->Bead);
      Count->Bonded--;
      Count->Unbonded++;
    }
    if (mt_i->Number == 1) { // copy beads only if molecule with id 'i' exists
      mol_i->Type = i;
      mol_i->Index = i;
      mol_i->Bead = malloc(sizeof *mol_i->Bead * mt_i->nBeads);
      memcpy(mol_i->Bead, mt_i->Bead, sizeof *mt_i->Bead * mt_i->nBeads);
    }
  }
} //}}}
// FillMoleculeTypeBonds() //{{{
static void FillMoleculeTypeBonds(SYSTEM *System, int (*bond)[3], int num) {
  COUNT *Count = &System->Count;
  // fill MoleculeType[].Bond array with bead indices
  for (int i = 0; i < num; i++) {
    int id[2];
    id[0] = bond[i][0];
    id[1] = bond[i][1];
    int bond_type = bond[i][2], mol = System->Bead[id[0]].Molecule;
    // warning - bonded beads in different molecules (skip the bond)  //{{{
    if (mol != System->Bead[id[1]].Molecule || mol == -1) {
      err_msg("bonded beads not in the same molecule; discarding this bond");
      PrintWarning();
      fprintf(stderr, "%sBead (molecule):", ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[0], ErrCyan(),
              ErrYellow(), System->Bead[id[0]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s)%s\n", ErrYellow(), id[1], ErrCyan(),
              ErrYellow(), System->Bead[id[1]].Molecule, ErrCyan(),
              ErrColourReset());
      continue;
    } //}}}
    MOLECULETYPE *mt_mol = &System->MoleculeType[mol];
    int bond = mt_mol->nBonds;
    mt_mol->nBonds++;
    if (bond == 0) { // first bond in a molecule - allocate one memory space
      mt_mol->Bond = malloc(sizeof *mt_mol->Bond);
    } else { // subsequent bonds - add one memory space
      mt_mol->Bond = s_realloc(mt_mol->Bond,
                               sizeof *mt_mol->Bond * mt_mol->nBonds);
    }
    mt_mol->Bond[bond][0] = id[0];
    mt_mol->Bond[bond][1] = id[1];
    mt_mol->Bond[bond][2] = bond_type;
  }
  // make the MoleculeType[].Bond bead indices go from 0 to nBeads //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    // 1) find lowest and highest index in the bond
    if (mt_i->nBonds) {
      int lowest = 1e9, highest = 0;
      for (int j = 0; j < mt_i->nBonds; j++) {
        if (mt_i->Bond[j][0] < lowest) {
          lowest = mt_i->Bond[j][0];
        }
        if (mt_i->Bond[j][1] < lowest) {
          lowest = mt_i->Bond[j][1];
        }
        if (mt_i->Bond[j][1] > highest) {
          highest = mt_i->Bond[j][1] + 1;
        }
      }
      // 2) a) make lowest index of Bond[][] 0
      //    b) find which indices are present to detect dicontinuities
      int diff = highest - lowest + 1; // +1 to count both highest & lowest ids
      if (diff < 0) { // never triggers; just to get rid of compiler warning
        err_msg("something went worng with bond indices!");
        PrintError();
        exit(1);
      }
      bool *present = calloc(diff, sizeof *present);
      for (int j = 0; j < mt_i->nBonds; j++) {
        int *id0 = &mt_i->Bond[j][0], *id1 = &mt_i->Bond[j][1];
        *id0 -= lowest;       // a)
        *id1 -= lowest;       // a)
        present[*id0] = true; // b)
        present[*id1] = true; // b)
      }
      // 3) find by how much to decrease indices (in case of discontinuities)
      int *decrease = calloc(diff, sizeof *decrease);
      for (int j = 0; j < diff; j++) {
        if (!present[j]) {
          for (int k = j; k < diff; k++) {
            decrease[k]++;
          }
        }
      }
      // 4) remove the discontinuities
      for (int j = 0; j < mt_i->nBonds; j++) {
        int *id0 = &mt_i->Bond[j][0], *id1 = &mt_i->Bond[j][1];
        mt_i->Bond[j][0] -= decrease[*id0];
        mt_i->Bond[j][1] -= decrease[*id1];
      }
      free(present);
      free(decrease);
      // 5) warning - index is too high; shouldn't happen
      for (int j = 0; j < mt_i->nBonds; j++) {
        if (mt_i->Bond[j][0] >= mt_i->nBeads ||
            mt_i->Bond[j][1] >= mt_i->nBeads) {
          err_msg("something is wrong in bond indices; should never happen!");
          PrintError();
        }
      }
    }
  } //}}}
} //}}}
// FillMoleculeTypeAngles() //{{{
static void FillMoleculeTypeAngles(SYSTEM *System, int (*angle)[4], int num) {
  COUNT *Count = &System->Count;
  // fill MoleculeType[].Bond array with bead indices
  for (int i = 0; i < num; i++) {
    int id[3];
    id[0] = angle[i][0];
    id[1] = angle[i][1];
    id[2] = angle[i][2];
    int angle_type = angle[i][3], mol = System->Bead[id[0]].Molecule;
    // warning - beads in different molecules (skip the bond)  //{{{
    if (mol != System->Bead[id[1]].Molecule ||
        mol != System->Bead[id[2]].Molecule || mol == -1) {
      err_msg("Beads in angle not in the same molecule; discarding this bond");
      PrintWarning();
      fprintf(stderr, "%sBead (molecule):", ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[0], ErrCyan(),
              ErrYellow(), System->Bead[id[0]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[1], ErrCyan(),
              ErrYellow(), System->Bead[id[1]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s)\n", ErrYellow(), id[2], ErrCyan(),
              ErrYellow(), System->Bead[id[2]].Molecule, ErrColourReset());
      continue;
    } //}}}
    MOLECULETYPE *mt_mol = &System->MoleculeType[mol];
    int angle = mt_mol->nAngles;
    mt_mol->nAngles++;
    if (angle == 0) {
      mt_mol->Angle = malloc(sizeof *mt_mol->Angle);
    } else {
      mt_mol->Angle = s_realloc(mt_mol->Angle,
                                sizeof *mt_mol->Angle * mt_mol->nAngles);
    }
    mt_mol->Angle[angle][0] = id[0];
    mt_mol->Angle[angle][1] = id[1];
    mt_mol->Angle[angle][2] = id[2];
    mt_mol->Angle[angle][3] = angle_type;
  }
  // make the MoleculeType[].Bond bead indices go from 0 to nBeads
  for (int i = 0; i < Count->MoleculeType; i++) {
    int lowest = 1e9;
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nAngles; j++) {
      if (mt_i->Angle[j][0] < lowest) {
        lowest = mt_i->Angle[j][0];
      }
      if (mt_i->Angle[j][1] < lowest) {
        lowest = mt_i->Angle[j][1];
      }
      if (mt_i->Angle[j][2] < lowest) {
        lowest = mt_i->Angle[j][2];
      }
    }
    for (int j = 0; j < mt_i->nAngles; j++) {
      mt_i->Angle[j][0] -= lowest;
      mt_i->Angle[j][1] -= lowest;
      mt_i->Angle[j][2] -= lowest;
      // warning - too high an intramolecular bead index; shouldn't happen //{{{
      if (mt_i->Angle[j][0] > mt_i->nBeads ||
          mt_i->Angle[j][1] > mt_i->nBeads ||
          mt_i->Angle[j][2] > mt_i->nBeads) {
        err_msg("error in angle's bead indices; should never happen!");
        PrintWarning();
      } //}}}
    }
  }
} //}}}
// FillMoleculeTypeDihedral() //{{{
static void FillMoleculeTypeDihedral(SYSTEM *System, int (*dihedral)[5],
                                     int num) {
  COUNT *Count = &System->Count;
  // fill MoleculeType[].Bond array with bead indices
  for (int i = 0; i < num; i++) {
    int id[4];
    id[0] = dihedral[i][0];
    id[1] = dihedral[i][1];
    id[2] = dihedral[i][2];
    id[3] = dihedral[i][3];
    int angle_type = dihedral[i][4], mol = System->Bead[id[0]].Molecule;
    // warning - beads in different molecules (skip the bond)  //{{{
    if (mol != System->Bead[id[1]].Molecule ||
        mol != System->Bead[id[2]].Molecule ||
        mol != System->Bead[id[3]].Molecule || mol == -1) {
      err_msg("Beads in dihedral do not share molecule; discarding this angle");
      PrintWarning();
      fprintf(stderr, "%sBead (molecule):", ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[0], ErrCyan(),
              ErrYellow(), System->Bead[id[0]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[1], ErrCyan(),
              ErrYellow(), System->Bead[id[1]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[2], ErrCyan(),
              ErrYellow(), System->Bead[id[2]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s)\n", ErrYellow(), id[3], ErrCyan(),
              ErrYellow(), System->Bead[id[3]].Molecule, ErrColourReset());
      continue;
    } //}}}
    MOLECULETYPE *mt_mol = &System->MoleculeType[mol];
    int angle = mt_mol->nDihedrals;
    mt_mol->nDihedrals++;
    if (angle == 0) {
      mt_mol->Dihedral = malloc(sizeof *mt_mol->Dihedral);
    } else {
      mt_mol->Dihedral = s_realloc(mt_mol->Dihedral, sizeof *mt_mol->Dihedral *
                                   mt_mol->nDihedrals);
    }
    mt_mol->Dihedral[angle][0] = id[0];
    mt_mol->Dihedral[angle][1] = id[1];
    mt_mol->Dihedral[angle][2] = id[2];
    mt_mol->Dihedral[angle][3] = id[3];
    mt_mol->Dihedral[angle][4] = angle_type;
  }
  // make the MoleculeType[].Bond bead indices go from 0 to nBeads
  for (int i = 0; i < Count->MoleculeType; i++) {
    int lowest = 1e9;
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nDihedrals; j++) {
      if (mt_i->Dihedral[j][0] < lowest) {
        lowest = mt_i->Dihedral[j][0];
      }
      if (mt_i->Dihedral[j][1] < lowest) {
        lowest = mt_i->Dihedral[j][1];
      }
      if (mt_i->Dihedral[j][2] < lowest) {
        lowest = mt_i->Dihedral[j][2];
      }
      if (mt_i->Dihedral[j][3] < lowest) {
        lowest = mt_i->Dihedral[j][3];
      }
    }
    for (int j = 0; j < mt_i->nDihedrals; j++) {
      mt_i->Dihedral[j][0] -= lowest;
      mt_i->Dihedral[j][1] -= lowest;
      mt_i->Dihedral[j][2] -= lowest;
      mt_i->Dihedral[j][3] -= lowest;
      // warning - too high an intramolecular bead index; shouldn't happen //{{{
      if (mt_i->Dihedral[j][0] > mt_i->nBeads ||
          mt_i->Dihedral[j][1] > mt_i->nBeads ||
          mt_i->Dihedral[j][2] > mt_i->nBeads ||
          mt_i->Dihedral[j][3] > mt_i->nBeads) {
        err_msg("error with dihedral bead indices; should never happen!");
        PrintWarning();
      } //}}}
    }
  }
} //}}}
// FillMoleculeTypeImproper() //{{{
static void FillMoleculeTypeImproper(SYSTEM *System, int (*improper)[5],
                                     int num) {
  COUNT *Count = &System->Count;
  // fill MoleculeType[].Bond array with bead indices
  for (int i = 0; i < num; i++) {
    int id[4];
    id[0] = improper[i][0];
    id[1] = improper[i][1];
    id[2] = improper[i][2];
    id[3] = improper[i][3];
    int angle_type = improper[i][4], mol = System->Bead[id[0]].Molecule;
    // warning - beads in different molecules (skip the bond)  //{{{
    if (mol != System->Bead[id[1]].Molecule ||
        mol != System->Bead[id[2]].Molecule ||
        mol != System->Bead[id[3]].Molecule || mol == -1) {
      err_msg("Beads in improper do not share molecule; discarding this angle");
      PrintWarning();
      fprintf(stderr, "%sBead (molecule):", ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[0], ErrCyan(),
              ErrYellow(), System->Bead[id[0]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[1], ErrCyan(),
              ErrYellow(), System->Bead[id[1]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id[2], ErrCyan(),
              ErrYellow(), System->Bead[id[2]].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s)\n", ErrYellow(), id[3], ErrCyan(),
              ErrYellow(), System->Bead[id[3]].Molecule, ErrColourReset());
      continue;
    } //}}}
    MOLECULETYPE *mt_mol = &System->MoleculeType[mol];
    int angle = mt_mol->nImpropers;
    mt_mol->nImpropers++;
    if (angle == 0) {
      mt_mol->Improper = malloc(sizeof *mt_mol->Improper);
    } else {
      mt_mol->Improper = s_realloc(mt_mol->Improper, sizeof *mt_mol->Improper *
                                   mt_mol->nImpropers);
    }
    mt_mol->Improper[angle][0] = id[0];
    mt_mol->Improper[angle][1] = id[1];
    mt_mol->Improper[angle][2] = id[2];
    mt_mol->Improper[angle][3] = id[3];
    mt_mol->Improper[angle][4] = angle_type;
  }
  // make the MoleculeType[].Bond bead indices go from 0 to nBeads
  for (int i = 0; i < Count->MoleculeType; i++) {
    int lowest = 1e9;
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nImpropers; j++) {
      if (mt_i->Improper[j][0] < lowest) {
        lowest = mt_i->Improper[j][0];
      }
      if (mt_i->Improper[j][1] < lowest) {
        lowest = mt_i->Improper[j][1];
      }
      if (mt_i->Improper[j][2] < lowest) {
        lowest = mt_i->Improper[j][2];
      }
      if (mt_i->Improper[j][3] < lowest) {
        lowest = mt_i->Improper[j][3];
      }
    }
    for (int j = 0; j < mt_i->nImpropers; j++) {
      mt_i->Improper[j][0] -= lowest;
      mt_i->Improper[j][1] -= lowest;
      mt_i->Improper[j][2] -= lowest;
      mt_i->Improper[j][3] -= lowest;
      // warning - too high an intramolecular bead index; shouldn't happen //{{{
      if (mt_i->Improper[j][0] > mt_i->nBeads ||
          mt_i->Improper[j][1] > mt_i->nBeads ||
          mt_i->Improper[j][2] > mt_i->nBeads ||
          mt_i->Improper[j][3] > mt_i->nBeads) {
        err_msg("error with improper bead indices; should never happen!");
        PrintWarning();
      } //}}}
    }
  }
} //}}}
// RemoveExtraTypes() { //{{{
/*
 * Remove bead and molecule types with .Number=0. It assumes the allocated
 * memory for BeadType and MoleculeType arrays of structures correspond to the
 * number of beads and molecules, respectively (i.e., not to the number of
 * types).
 */
static void RemoveExtraTypes(SYSTEM *System) {
  COUNT *Count = &System->Count;
  // BeadType & Bead[].Type
  int count = 0;
  int *bt_old_to_new = malloc(sizeof *bt_old_to_new * Count->BeadType);
  for (int i = 0; i < Count->BeadType; i++) {
    if (System->BeadType[i].Number != 0) {
      int bt_id = count;
      count++;
      if (bt_id != i) {
        System->BeadType[bt_id] = System->BeadType[i];
      }
      bt_old_to_new[i] = bt_id;
    }
  }
  Count->BeadType = count;
  for (int i = 0; i < Count->Bead; i++) {
    int old_type = System->Bead[i].Type;
    System->Bead[i].Type = bt_old_to_new[old_type];
  }
  // sync MoleculeType[].Bead (i.e., bead types) with the new BeadType
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (System->MoleculeType[i].Number != 0) {
      for (int j = 0; j < System->MoleculeType[i].nBeads; j++) {
        int id = System->MoleculeType[i].Bead[j];
        System->MoleculeType[i].Bead[j] = bt_old_to_new[id];
      }
    }
  }
  free(bt_old_to_new);
  // MoleculeType & Molecule
  count = 0;
  Count->Molecule = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    if (mt_i->Number != 0) {
      Count->Molecule += mt_i->Number;
      int mt_id = count;
      count++;
      if (mt_id != i) {
        // MoleculeType struct
        MOLECULETYPE *mt_new = &System->MoleculeType[mt_id];
        *mt_new = *mt_i;
        mt_new->Bead = malloc(sizeof *mt_new->Bead * mt_new->nBeads);
        memcpy(mt_new->Bead, mt_i->Bead, sizeof *mt_i->Bead * mt_i->nBeads);
        free(mt_i->Bead);
        if (mt_new->nBonds > 0) {
          mt_new->Bond = malloc(sizeof *mt_new->Bond * mt_new->nBonds);
          memcpy(mt_new->Bond, mt_i->Bond, sizeof *mt_i->Bond * mt_i->nBonds);
          free(mt_i->Bond);
        }
        if (mt_new->nAngles > 0) {
          mt_new->Angle = malloc(sizeof *mt_new->Angle * mt_new->nAngles);
          memcpy(mt_new->Angle, mt_i->Angle,
                 sizeof *mt_i->Angle * mt_i->nAngles);
          free(mt_i->Angle);
        }
        if (mt_new->nDihedrals > 0) {
          mt_new->Dihedral = malloc(sizeof *mt_new->Dihedral *
                                    mt_new->nDihedrals);
          memcpy(mt_new->Dihedral, mt_i->Dihedral,
                 sizeof *mt_i->Dihedral * mt_i->nDihedrals);
          free(mt_i->Dihedral);
        }
        if (mt_new->nImpropers > 0) {
          mt_new->Improper = malloc(sizeof *mt_new->Improper *
                                    mt_new->nImpropers);
          memcpy(mt_new->Improper, mt_i->Improper,
                 sizeof *mt_i->Improper * mt_i->nImpropers);
          free(mt_i->Improper);
        }
        // Molecule struct
        MOLECULE *mol_new = &System->Molecule[mt_id];
        mol_new->Type = mt_id;
        mol_new->Index = i;
        mol_new->Bead = malloc(sizeof *mol_new->Bead * mt_new->nBeads);
        for (int j = 0; j < mt_new->nBeads; j++) {
          int id = System->Molecule[i].Bead[j];
          System->Molecule[mt_id].Bead[j] = id;
          System->Bead[id].Molecule = mt_id;
        }
        free(System->Molecule[i].Bead);
      }
    }
  }
  Count->MoleculeType = count;
} //}}}
  //}}}

// THE VISIBLE FUNCTIONS
/*
 * Read the provided structure file and extra info from the
 * provided coordinate file if necessary.
 */ //{{{
SYSTEM ReadStructure(SYS_FILES f, bool detailed) {
  SYSTEM System;
  switch (f.stru.type) {
    case VTF_FILE:
    case VSF_FILE:
      System = VtfReadStruct(f.stru.name, detailed);
      break;
    case XYZ_FILE:
      System = XyzReadStruct(f.stru.name);
      break;
    case LTRJ_FILE:
      System = LtrjReadStruct(f.stru.name);
      break;
    case LDATA_FILE:
      System = LmpDataReadStruct(f.stru.name);
      break;
    case FIELD_FILE:
      System = FieldRead(f.stru.name);
      break;
    default:
      err_msg("unspecified structure file; should never happen!");
      PrintError();
      exit(1);
  }
  // read extra stuff from coordinate file if necessary
  if (f.coor.type == LTRJ_FILE && System.Box.Volume == -1) {
    System.Box = LtrjReadPBC(f.coor.name);
  }
  if (f.coor.type == VCF_FILE) {
    System.Count.BeadCoor = VtfReadNumberOfBeads(f.coor.name);
    if (System.Count.BeadCoor < 0) {
      exit(1);
    }
    if (System.Box.Volume == -1) {
      System.Box = VtfReadPBC(f.coor.name);
    }
  }
  WarnChargedSystem(System, f.stru.name, "\0", "\0");
  // warn if missing box dimensions (unless it's a pbc-less file type)
  if (System.Box.Volume == -1 && f.stru.type != VSF_FILE &&
      f.stru.type != FIELD_FILE && f.stru.type != XYZ_FILE) {
    err_msg("unspecified box dimensions in structure definition");
    PrintWarnFile(f.stru.name, "\0", "\0");
  }
  return System;
} //}}} //}}}
/*
 * Read a single timestep from the provided coordinate file.
 */ //{{{
bool ReadTimestep(SYS_FILES f, FILE *fr, SYSTEM *System, int *line_count) {
  switch (f.coor.type) {
    case VTF_FILE:
    case VCF_FILE:
      if (VtfReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    case XYZ_FILE:
      if (XyzReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    case LTRJ_FILE:;
      int line = *line_count;
      if (LtrjReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      // skip this step if it's a first one that contain only zeroes
      // ...huh? Why would this be a thing?
      if (line == 0 && System->Count.BeadCoor == System->Count.Bead) {
        bool zeroes = true;
        for (int i = 0; i < System->Count.Bead; i++) {
          if (System->Bead[i].Position[0] != 0 ||
              System->Bead[i].Position[1] != 0 ||
              System->Bead[i].Position[2] != 0) {
            zeroes = false;
            break;
          }
        }
        if (zeroes &&
            LtrjReadTimestep(fr, f.coor.name, System, line_count) < 0) {
          return false;
        }
      }
      break;
    case LDATA_FILE:
      if (LmpDataReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    default:
      snprintf(ERROR_MSG, LINE, "no action specified for coor_type %s%d",
               ErrYellow(), f.coor.type);
      PrintError();
      exit(1);
  }
  return true;
} //}}}
/*
 * Skip a single timestep from the provided coordinate file.
 */ //{{{
bool SkipTimestep(SYS_FILES f, FILE *fr, int *line_count) {
  switch (f.coor.type) {
    case VTF_FILE:
    case VCF_FILE:
      if (VtfSkipTimestep(fr, f.coor.name, f.stru.name, line_count) < 0) {
        return false;
      }
      break;
    case XYZ_FILE:
      if (!XyzSkipTimestep(fr, f.coor.name, line_count)) {
        return false;
      }
      break;
    case LTRJ_FILE:
      if (LtrjSkipTimestep(fr, f.coor.name, line_count) < 0) {
        return false;
      }
      break;
    case LDATA_FILE:
      err_msg("lammps data file contains only one step; should never trigger!");
      PrintWarnFile(f.coor.name, "\0", "\0");
      return false;
    default:
      snprintf(ERROR_MSG, LINE, "no action specified for coor_type %s%d",
               ErrYellow(), f.coor.type);
      PrintError();
      exit(1);
  }
  return true;
} //}}}
// Read aggregates from a single timestep //{{{
int ReadAggregates(FILE *fr, char file[], SYSTEM *System, AGGREGATE Aggregate[],
                   int *line_count) {
  COUNT *Count = &System->Count;
  // read Step:/Last Step: line //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE,
             "premature end of %s%s%s file; "
             "stopping reading",
             ErrYellow(), file, ErrCyan());
    PrintWarnFile(file, "\0", "\0");
    return -2;
  }
  if (words < 1) {
    snprintf(ERROR_MSG, LINE,
             "blank line instead of Step:/Last Step: line; "
             "stopping reading %s%s",
             ErrYellow(), file);
    PrintWarnFileLine(file, *line_count);
    return -2;
  }
  //}}}
  if (strcasecmp(split[0], "Last") == 0) {
    return -1; // no error - just the end of the timesteps
  }
  // read number of aggregates //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE, "premature end of %s%s%s file; stopping reading",
             ErrYellow(), file, ErrCyan());
    PrintWarnFile(file, "\0", "\0");
    return -2;
  }
  long int val;
  if (words < 1 || !IsNaturalNumber(split[0], &val)) {
    snprintf(ERROR_MSG, LINE,
             "incorrect line with number of aggregates stopping reading %s%s",
             ErrYellow(), file);
    PrintWarnFileLine(file, *line_count);
    return -2;
  } //}}}
  Count->Aggregate = val;
  // go through all aggregates, filling in the molecules
  for (int i = 0; i < Count->Aggregate; i++) {
    fscanf(fr, "%d :", &Aggregate[i].nMolecules);
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol;
      fscanf(fr, "%d", &mol);
      mol--; // in agg file, the numbers correspond to vmd
      Aggregate[i].Molecule[j] = mol;
      System->Molecule[mol].Aggregate = i;
    }
    while (getc(fr) != '\n')
      ;
  }
  // fill the rest of aggregate info
  for (int i = 0; i < Count->Aggregate; i++) {
    int count = 0;
    double mass = 0;
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      MOLECULE *mol = &System->Molecule[Aggregate[i].Molecule[j]];
      MOLECULETYPE *mt = &System->MoleculeType[mol->Type];
      if (mt->Mass == MASS || mass == -1) {
        mass = -1;
      } else {
        mass += mt->Mass;
      }
      mol->Aggregate = i;
      for (int k = 0; k < mt->nBeads; k++) {
        Aggregate[i].Bead[count] = mol->Bead[k];
        count++;
      }
    }
    Aggregate[i].nBeads = count;
    if (mass == -1) {
      Aggregate[i].Mass = MASS; // unspecified mass
    } else {
      Aggregate[i].Mass = mass; // valid mass
    }
  }
  return 1;
} //}}}
