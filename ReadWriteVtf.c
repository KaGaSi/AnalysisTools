#include "ReadWriteVtf.h"
#include "General.h"
#include "ReadWrite.h"

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

static int VtfCheckLineType(const char *file, const int line_count);
static int VtfCheckCoorOrderedLine(double coor[3]);
static int VtfCheckCoorIndexedLine(double coor[3], long *index);
static int VtfCheckCoordinateLine(double coor[3], long *index);
static int VtfCheckTimestepLine();
static int VtfCheckPbcLine();
static bool VtfCheckAtomLine();
static bool VtfCheckBondLine();
// Find position of keywords in an atom line
static int *VtfAtomLineValues();
// process pbc line from a vtf file
static bool VtfPbcLine(BOX *box, const int ltype);
// read timestep preamble in a vtf coordinate file
static int VtfReadTimestepPreamble(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count);
// read variable-size indexed coordinate block in a vtf file
static int VtfReadCoorBlockIndexed(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count);
// read ordered coordinate block in a vtf file
static int VtfReadCoorBlockOrdered(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count);
static void SimplifyResid(SYSTEM *System);

// VtfReadStruct() //{{{
/*
 * Read system information from vsf/vtf structure file. It can recognize bead
 * and molecule types based either on name only or on all information (name,
 * mass, charge, and radius for bead types; bead order, bonds, angles, and
 * dihedrals for molecule types).
 */
SYSTEM VtfReadStruct(const char *file, const bool detailed) {
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  FILE *fr = OpenFile(file, "r");
  // define variables and structures and arrays //{{{
  int line_count = 0;   // total number of lines
  int count_atoms = 0;  // number of 'atom <id>'
  int default_atom = 0; // line number of the first 'atom default' line
  // int Bond = 0;  // number of bonds
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
                   " line %s%d%s as the default line%s", ErrCyan(),
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
      Count->Bond++;
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
             "omitted); %s%d%s bead(s) undefined",
             ErrYellow(), undefined, ErrRed());
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
  int (*bond)[5] = NULL;
  if (Count->Bond != 0) {
    bond = calloc(Count->Bond, sizeof *bond);
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
  Count->Bond = 0;
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
            if (value[4] == -1) { // 'resname' vtf keyword missing
              s_strcpy(mt_resid->Name, "m", MOL_NAME);
              mt_resid->Flag = false;
            } else { // 'resname' vtf keyword present
              s_strcpy(mt_resid->Name, split[value[4]], MOL_NAME);
              mt_resid->Flag = true;
            }
            mt_resid->Number = 1;
            mt_resid->nBeads = 1;
            mt_resid->Bead = malloc(sizeof *mt_resid->Bead);
            mt_resid->Bead[0] = id; // bead type = bead index
          } else { // not new moleclue type
            int bead = mt_resid->nBeads;
            mt_resid->nBeads++;
            mt_resid->Bead = s_realloc(mt_resid->Bead, sizeof *mt_resid->Bead *
                                       mt_resid->nBeads);
            mt_resid->Bead[bead] = id; // bead type = bead index
            if (value[4] != -1 && strcmp(mt_resid->Name, "m") == 0) {
              s_strcpy(mt_resid->Name, split[value[4]], MOL_NAME);
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
        bond[Count->Bond][0] = val;
        IsIntegerNumber(index[1], &val);
        bond[Count->Bond][1] = val;
      } else { // case 'bond <id>: <id>'
        // IsIntegerNumber() ignores the trailing ':'
        IsIntegerNumber(split[1], &val);
        bond[Count->Bond][0] = val;
        IsIntegerNumber(split[2], &val);
        bond[Count->Bond][1] = val;
      }
      bond[Count->Bond][2] = -1;
      // assure index1<index2 (may be unnecessary, but definitely won't hurt)
      if (bond[Count->Bond][0] > bond[Count->Bond][1]) {
        SwapInt(&bond[Count->Bond][0], &bond[Count->Bond][1]);
      }
      Count->Bond++; //}}}
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
  FillAllMTypeStuff(&Sys, bond, NULL, NULL, NULL);
  if (Count->Bond != 0) {
    free(bond);
  }
  RemoveExtraTypes(&Sys);
  MergeBeadTypes(&Sys, detailed);
  MergeMoleculeTypes(&Sys);
  FillSystemNonessentials(&Sys, true);
  CheckSystem(Sys, file);
  Sys.Box = VtfReadPBC(file);
  return Sys;
} //}}}
// VtfReadTimestep() //{{{
int VtfReadTimestep(FILE *fr, const char *file,
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
int VtfSkipTimestep(FILE *fr, const char *file,
                    const char *vsf_file, int *line_count) {
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
static int VtfCheckLineType(const char *file, const int line_count) { //{{{
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
BOX VtfReadPBC(const char *file) { //{{{
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
static bool VtfPbcLine(BOX *box, const int ltype) { //{{{
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
static int VtfReadTimestepPreamble(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count) {
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
static int VtfReadCoorBlockIndexed(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count) {
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
static int VtfReadCoorBlockOrdered(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      snprintf(ERROR_MSG, LINE, "premature end of ordered "
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
int VtfReadNumberOfBeads(const char *file) { //{{{
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

// VtfWriteStruct() //{{{
void WriteBeadTypeInfo(FILE *fw, const BEADTYPE bt) {
  fprintf(fw, " name %8s", bt.Name);
  if (bt.Mass != MASS && bt.Mass != HIGHNUM) {
    fprintf(fw, " mass %12f", bt.Mass);
  }
  if (bt.Charge != CHARGE && bt.Charge != HIGHNUM) {
    fprintf(fw, " charge %12f", bt.Charge);
  }
  if (bt.Radius != RADIUS && bt.Radius != HIGHNUM) {
    fprintf(fw, " radius %12f", bt.Radius);
  }
}
void VtfWriteStruct(char *file, SYSTEM System, int type_def,
                    const int argc, char **argv) {
  SimplifyResid(&System);
  PrintByline(file, argc, argv);
  FILE *fw = OpenFile(file, "a");
  COUNT *Count = &System.Count;
  // default bead type //{{{
  if (type_def == -1) {
    // find most common type of bead and make it default
    int *count = calloc(Count->BeadType, sizeof *count);
    for (int i = 0; i < Count->Bead; i++) {
      if (System.Bead[i].Molecule == -1) {
        int type = System.Bead[i].Type;
        count[type]++;
      }
    }
    int max = 0;
    for (int i = 0; i < Count->BeadType; i++) {
      if (count[i] > max) {
        max = count[i];
        type_def = i;
      }
    }
    free(count);
  } //}}}
  // print default bead type //{{{
  if (type_def != -1) {
    BEADTYPE *bt = &System.BeadType[type_def];
    fprintf(fw, "atom default");
    WriteBeadTypeInfo(fw, *bt);
    putc('\n', fw);
  } //}}}
  // print beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    int btype = System.Bead[i].Type;
    int mol = System.Bead[i].Molecule;
    // print beads that are non-default, in a molecule, or have the highest id
    bool print = false;
    if (btype != type_def || mol != -1 || i == (Count->Bead - 1)) {
      print = true;
    }
    BEADTYPE *bt = &System.BeadType[btype];
    if (print) {
      fprintf(fw, "atom %7d", i);
      WriteBeadTypeInfo(fw, *bt);
      if (mol != -1) {
        int mtype = System.Molecule[mol].Type;
        int id = System.Molecule[mol].Index;
        char name[8];
        if (snprintf(name, 8, "%s", System.MoleculeType[mtype].Name) < 0) {
          ErrorSnprintf();
        }
        fprintf(fw, " resname %10s", name);
        fprintf(fw, " resid %5d", id);
      }
      putc('\n', fw);
    }
  } //}}}
  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol = &System.Molecule[i];
    MOLECULETYPE *mt = &System.MoleculeType[mol->Type];
    if (mt->nBonds > 0) {
      fprintf(fw, "# resid %d\n", i + 1); // in VMD resid start with 1
      for (int j = 0; j < mt->nBonds; j++) {
        int id[2];
        for (int aa = 0; aa < 2; aa++) {
          id[aa] = mt->Bond[j][aa];
          id[aa] = mol->Bead[id[aa]];
        }
        fprintf(fw, "bond %6d: %6d\n", id[0], id[1]);
      }
    }
  } //}}}
  WriteBoxLengthAngles(fw, System.Box);
  putc('\n', fw);
  fclose(fw);
} //}}}
// VtfWriteCoorIndexed() //{{{
void VtfWriteCoorIndexed(FILE *fw, const bool *write, const SYSTEM System) {
  fprintf(fw, "indexed\n");
  // print box size if present //{{{
  const BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "pbc %lf %lf %lf", box->Length[0],
                                   box->Length[1],
                                   box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
      fprintf(fw, "    %.3f %.3f %.3f", box->alpha, box->beta, box->gamma);
    }
    putc('\n', fw);
  } //}}}
  bool none = true;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    if (write[id]) {
      none = false;
      fprintf(fw, "%8d %8.4f %8.4f %8.4f\n", id, bead->Position[0],
                                                 bead->Position[1],
                                                 bead->Position[2]);
    }
  }
  if (none) {
    err_msg("no beads to save");
    PrintWarning();
  }
} //}}}
static void SimplifyResid(SYSTEM *System) { //{{{
  int lowest = 1e6;
  for (int i = 0; i < System->Count.Molecule; i++) {
    if (System->Molecule[i].Index < lowest) {
      lowest = System->Molecule[i].Index;
    }
  }
  for (int i = 0; i < System->Count.Molecule; i++) {
    System->Molecule[i].Index += -lowest + 1; // start from 1
  }
} //}}}
