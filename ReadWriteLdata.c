#include "ReadWriteLdata.h"
#include "Errors.h"
#include "General.h"
#include "ReadWrite.h"
#include "System.h"

// TODO: LmpDataReadDihedralCoeffs() and LmpDataReadImproperCoeffs() should read
//       up to three numbers, not assuming any format of the potential

// Helper functions for lmpdata file
// increment line count and read the line
static int CountLineReadLine(int *line_count, FILE *fr, const char *file,
                             const char *msg, const int mode);
// read header of the lammps data file
static void LmpDataReadHeader(FILE *fr, const char *file,
                              SYSTEM *System, int *line_count);
// read body of the lammps data file
static void LmpDataReadBody(FILE *fr, const char *file, SYSTEM *System,
                            int *line_count);
// functions to read the various sections in the lammps data file
static void LmpDataReadMasses(FILE *fr, const char *file, BEADTYPE *name_mass,
                              const int atom_types, int *line_count);
// TODO: eventually non-harmonic (non-cvff) stuff; also then add to FIELD
static void ReadCoeff(FILE *fr, const char *file, SYSTEM *System,
                      int *line_count, const int type);
static void ReadBodySection(SYSTEM *System, FILE *fr,
                            const char *file, int *line_count);
static void LmpDataReadAtoms(FILE *fr, const char *file, SYSTEM *System,
                             const BEADTYPE *name_mass,
                             int *line_count, const int mode);
static void LmpDataReadBADISection(FILE *fr, const char *file,
                                   const COUNT Count, int (*array)[5],
                                   int *line_count, const int type);
static void GetBeadVelocity(SYSTEM *System, const char *file,
                            const int line_count);
static bool ReadPbc(BOX *box);
static int CheckAtomsMode(const int words, char **split, const char *file);
static void CheckAtomsLine(const int mode, const int beads, long *id,
                           long *resid, long *type, double *q, double pos[3],
                           const char *file, const int line_count);
static void WriteStuff(FILE *fw, const SYSTEM System, const int mol,
                       int *count, const int num, int (**arr)[5], const int n);

// read files
// TODO: Atoms & Velocities section can be switched
// TODO: add Forces section reading
SYSTEM LmpDataReadStruct(const char *file) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  FILE *fr = OpenFile(file, "r");
  int line_count = 0;
  LmpDataReadHeader(fr, file, &System, &line_count);
  LmpDataReadBody(fr, file, &System, &line_count);
  fclose(fr);
  CalculateBoxData(&System.Box, 1);
  // till RemoveExtraTypes(), each bead is its own bead type
  Count->BeadType = Count->Bead;
  RemoveExtraTypes(&System);
  MergeBeadTypes(&System, true);
  MergeMoleculeTypes(&System);
  FillSystemNonessentials(&System, true);
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
// read lammps data file as a coordinate file //{{{
int LmpDataReadTimestep(FILE *fr, const char *file,
                        SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // set 'in timestep' to all beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    System->Bead[i].InTimestep = true;
  } //}}}
  // set 'in timestep' to all molecules //{{{
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].InTimestep = false;
  } //}}}
  // fill MoleculeCoor array - all molecules are present //{{{
  Count->MoleculeCoor = Count->Molecule;
  for (int i = 0; i < Count->MoleculeCoor; i++) {
    System->MoleculeCoor[i] = i;
  } //}}}
  *line_count = 0;
  CountLineReadLine(line_count, fr, file, "empty file", 1);
  // read numer of atoms & box size //{{{
  // read until the first capital-letter-starting line
  fpos_t position;
  do {
    fgetpos(fr, &position);
    CountLineReadLine(line_count, fr, file,
                      "incomplete lammps data file header", 1);
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
    }
    if (!ReadPbc(&System->Box)) {
      goto error;
    }
  } while (words == 0 || split[0][0] < 'A' || split[0][0] > 'Z');
  //  return file pointer to before the first capital-letter-starting line
  fsetpos(fr, &position);
  (*line_count)--;
  CalculateBoxData(&System->Box, 1); //}}}
  // find 'Atoms' section (and skip the next blank line) //{{{
  do {
    if (CountLineReadLine(line_count, fr, file, "", -2) < 0) {
      return -2;
    }
  } while (words == 0 || strcmp(split[0], "Atoms") != 0);
  int mode = CheckAtomsMode(words, split, file);
  if (CountLineReadLine(line_count, fr, file, "", -2) < 0) {
    return -2;
  } //}}}
  // read atom lines //{{{
  Count->BeadCoor = Count->Bead;
  for (int i = 0; i < Count->Bead; i++) {
    if (CountLineReadLine(line_count, fr, file, "", -2) < 0) {
      return -2;
    }
    long id, resid, type;
    double pos[3], q;
    CheckAtomsLine(mode, Count->Bead, &id, &resid,
                   &type, &q, pos, file, *line_count);
    id--; // in lammps data file, ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = pos[dd] - System->Box.Low[dd];
    }
    System->BeadCoor[i] = id;
  } //}}}
  // find 'Velocities' section (and skip the next blank line) //{{{
  do {
    if (CountLineReadLine(line_count, fr, file, "", -2) < 0) {
      ChangeBoxByLow(System, -1);
      FillInCoor(System);
      return 1; // Velocities section is not mandatory
    }
  } while (words == 0 || strcmp(split[0], "Velocities") != 0);
  if (CountLineReadLine(line_count, fr, file, "", -2) < 0) {
    ChangeBoxByLow(System, -1);
    FillInCoor(System);
    return -2;
  } //}}}
  // read velocity lines //{{{
  for (int i = 0; i < Count->Bead; i++) {
    if (CountLineReadLine(line_count, fr, file, "", -2) < 0) {
      return -2;
    }
    GetBeadVelocity(System, file, *line_count);
  } //}}}
  ChangeBoxByLow(System, -1);
  FillInCoor(System);
  return 1;
  error: // unrecognised line //{{{
    err_msg("unrecognised line in the file header");
    PrintErrorFileLine(file, *line_count);
    exit(1); //}}}
} //}}}

// increment line count and read the line //{{{
static int CountLineReadLine(int *line_count, FILE *fr, const char *file,
                             const char *msg, const int mode) {
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (mode > 0) {
      err_msg(msg);
      ErrorEOF(file, ERROR_MSG);
      exit(1);
    } else {
      return mode;
    }
  }
  return 0;
} //}}}
// read header //{{{
static void LmpDataReadHeader(FILE *fr, const char *file,
                              SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  *line_count = 0;
  CountLineReadLine(line_count, fr, file, "empty file", 1);
  // read until a line starting with capital letter //{{{
  fpos_t position;
  do {
    fgetpos(fr, &position);
    CountLineReadLine(line_count, fr, file,
                      "incomplete lammps data file header", 1);
    // evaluate the line
    long val;
    // <int> atom types //{{{
    if (words > 2 && strcmp(split[1], "atom") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      Count->BeadType = val; //}}}
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
    // <int> atoms //{{{
    } else if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      Count->Bead = val;
      System->Bead = s_realloc(System->Bead,
                               Count->Bead * sizeof *System->Bead);
      /*
       * Allocate the same memory for System.BeadType as each bead can have
       * different charge, so at first, each bead will be its own type
       */
      System->BeadType = s_realloc(System->BeadType,
                                   Count->Bead * sizeof *System->BeadType);
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
      Count->Improper = val;
    } //}}}
    if (!ReadPbc(&System->Box)) {
      goto error;
    }
  } while (words == 0 || split[0][0] < 'A' || split[0][0] > 'Z'); //}}}
  // return file pointer to before the first capital-letter-starting line
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
  } //}}}
  return;
error: // unrecognised line //{{{
  err_msg("unrecognised line in the file header");
  PrintErrorFileLine(file, *line_count);
  exit(1); //}}}
} //}}}
// read body //{{{
static void LmpDataReadBody(FILE *fr, const char *file, SYSTEM *System,
                            int *line_count) {
  COUNT *Count = &System->Count;
  // BEADTYPE *name_mass = calloc(atom_types, sizeof *name_mass);
  BEADTYPE *name_mass = calloc(Count->BeadType, sizeof *name_mass);
  // create arrays for bonds/angles/dihedrals/impropers //{{{
  int (*bond)[5], (*angle)[5], (*dihedral)[5], (*improper)[5];
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
      LmpDataReadMasses(fr, file, name_mass, Count->BeadType, line_count);
    } else if (words > 1 && strcmp(split[0], "Bond") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      bonds[0] = true;
      ReadCoeff(fr, file, System, line_count, 0);
      // check the parameters are correct - TODO: add more than harm
      for (int i = 0; i < Count->BondType; i++) {
        PARAMS *p = &System->BondType[i];
        if (p->a < 0 || p->b < 0) {
          snprintf(ERROR_MSG, LINE, "wrong Bond Coeff: %s%d %lf %lf%s",
                   ErrYellow(), i + 1, p->a, p->b, ErrRed());
          PrintErrorFile(file, "", "");
          exit(1);
        }
      }
    } else if (words > 1 && strcmp(split[0], "Angle") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      angles[0] = true;
      // check the parameters are correct - TODO: add more than harm
      ReadCoeff(fr, file, System, line_count, 1);
      for (int i = 0; i < Count->AngleType; i++) {
        PARAMS *p = &System->AngleType[i];
        if (p->a < 0 || p->b < 0) {
          snprintf(ERROR_MSG, LINE, "wrong Angle Coeff: %s%d %lf %lf%s",
                   ErrYellow(), i + 1, p->a, p->b, ErrRed());
          PrintErrorFile(file, "", "");
          exit(1);
        }
      }
    } else if (words > 1 && strcmp(split[0], "Dihedral") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      dihedrals[0] = true;
      ReadCoeff(fr, file, System, line_count, 2);
      // check the parameters are correct - TODO: add more than harm
      for (int i = 0; i < Count->DihedralType; i++) {
        PARAMS *p = &System->DihedralType[i];
        if (p->a < 0 || fabs(p->b) != 1 || p->c < 0) {
          snprintf(ERROR_MSG, LINE, "wrong Dihedral Coeff: %s%d %lf %lf %lf%s",
                   ErrYellow(), i + 1, p->a, p->b, p->c, ErrRed());
          PrintErrorFile(file, "", "");
          exit(1);
        }
      }
    } else if (words > 1 && strcmp(split[0], "Improper") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      impropers[0] = true;
      ReadCoeff(fr, file, System, line_count, 3);
      // check the parameters are correct - TODO: add more than cvff
      for (int i = 0; i < Count->ImproperType; i++) {
        PARAMS *p = &System->ImproperType[i];
        if (p->a < 0 || fabs(p->b) != 1 || p->c < 0) {
          snprintf(ERROR_MSG, LINE, "wrong Improper Coeff: %s%d %lf %lf %lf%s",
                   ErrYellow(), i + 1, p->a, p->b, p->c, ErrRed());
          PrintErrorFile(file, "", "");
          exit(1);
        }
      }
    } else if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      atoms = true;
      int mode = CheckAtomsMode(words, split, file);
      LmpDataReadAtoms(fr, file, System, name_mass, line_count, mode);
    } else if (words > 0 && strcmp(split[0], "Velocities") == 0) {
      ReadBodySection(System, fr, file, line_count);
    } else if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      bonds[1] = true;
      LmpDataReadBADISection(fr, file, *Count, bond, line_count, 0);
    } else if (words > 0 && strcmp(split[0], "Angles") == 0) {
      angles[1] = true;
      LmpDataReadBADISection(fr, file, *Count, angle, line_count, 1);
    } else if (words > 0 && strcmp(split[0], "Dihedrals") == 0) {
      dihedrals[1] = true;
      // LmpDataReadDihedrals(fr, file, *Count, dihedral, line_count);
      LmpDataReadBADISection(fr, file, *Count, dihedral, line_count, 2);
    } else if (words > 0 && strcmp(split[0], "Impropers") == 0) {
      impropers[1] = true;
      // LmpDataReadImpropers(fr, file, *Count, improper, line_count);
      LmpDataReadBADISection(fr, file, *Count, improper, line_count, 3);
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
    SortArray(mt->Bead, mt->nBeads, 0, 'i');
  }
  CopyMoleculeTypeBeadsToMoleculeBeads(System);
  FillAllMTypeStuff(System, bond, angle, dihedral, improper);
  free(bond);
  free(angle);
  free(dihedral);
  free(improper);
} //}}}
// read Masses section //{{{
static void LmpDataReadMasses(FILE *fr, const char *file, BEADTYPE *name_mass,
                              const int atom_types, int *line_count) {
  // skip one line
  CountLineReadLine(line_count, fr, file, "incomplete 'Masses' section", 1);
  for (int i = 0; i < atom_types; i++) {
    CountLineReadLine(line_count, fr, file, "incomplete 'Masses' section", 1);
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
      s_strcpy(name_mass[type].Name, split[3], BEAD_NAME);
    } else {
      snprintf(name_mass[type].Name, BEAD_NAME, "b%d", i);
    }
  }
} //}}}
// read one Bond/Angle/Dihedral/Improper Coeffs section //{{{
static void ReadCoeff(FILE *fr, const char *file, SYSTEM *System,
                      int *line_count, const int type) {
  // assign proper variables based on Bond/Angle/Dihedral/Improper section //{{{
  PARAMS *CoeffType = NULL;
  int *count_type;
  char msg[2][10];
  if (type == 0) {
    CoeffType = System->BondType;
    count_type = &System->Count.BondType;
    s_strcpy(msg[0], "Bond", 10);
    s_strcpy(msg[1], "bond", 10);
  } else if (type == 1) {
    CoeffType = System->AngleType;
    count_type = &System->Count.AngleType;
    s_strcpy(msg[0], "Angle", 10);
    s_strcpy(msg[1], "angle", 10);
  } else if (type == 2) {
    CoeffType = System->DihedralType;
    count_type = &System->Count.DihedralType;
    s_strcpy(msg[0], "Dihedral", 10);
    s_strcpy(msg[1], "dihedral", 10);
  } else if (type == 3) {
    CoeffType = System->ImproperType;
    count_type = &System->Count.ImproperType;
    s_strcpy(msg[0], "Improper", 10);
    s_strcpy(msg[1], "improper", 10);
  } else {
    err_msg("ReadCoeff: type must be 0 to 3");
    PrintError();
    exit(1);
  } //}}}
  // error - no bond types
  if (*count_type == 0) {
    snprintf(ERROR_MSG, LINE, "%s Coeffs in a file with no %s types",
             msg[0], msg[1]);
    err_msg("");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    return;
  }
  // skip one line
  CountLineReadLine(line_count, fr, file, "incomplete 'Coeffs' section", 1);
  for (int i = 0; i < *count_type; i++) {
    CountLineReadLine(line_count, fr, file, "incomplete 'Coeffs' section", 1);
    long num;
    double a[4] = {0, 0, 0, 0};
    // error - wrong line
    // Coeff line must have format <coeff type id> <num> [3x<num>]
    if (words < 2 ||
        !IsNaturalNumber(split[0], &num) || num > *count_type ||
        !IsRealNumber(split[1], &a[0]) || a[0] < 0) {
      snprintf(ERROR_MSG, LINE, "wrong line in '%s Coeffs' section", msg[0]);
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    // optional extra two numbers
    for (int j = 2; j < words && j < 5; j++) {
      if (!IsRealNumber(split[j], &a[j-1])) {
        a[j-1] = 0;
        break;
      }
    }
    CoeffType[num-1].a = a[0] * 2; // lammps uses k/2 as spring constant
    CoeffType[num-1].b = a[1];
    CoeffType[num-1].c = a[2];
    CoeffType[num-1].d = a[3];
  }
} //}}}
// TODO: should this be somehow connected to reading atom section? I guess not
//       since the atom stuff needs more arguments - rename this therefore
// TODO: how is it with the reading timestep? Can I use this or not? Oh, no; I
//       guess it's just the GetBeadVelocity bit as timestep requires not
//       exiting on error which the CountLineReadLine bit does! Actually, it
//       doesn't anymore! So yeah, use it! Somehow
// ReadBodySection() //{{{
static void ReadBodySection(SYSTEM *System, FILE *fr,
                            const char *file, int *line_count) {
  COUNT *Count = &System->Count;
  // skip one line
  CountLineReadLine(line_count, fr, file, "incomplete 'Velocities' section", 1);
  for (int i = 0; i < Count->Bead; i++) {
    CountLineReadLine(line_count, fr, file,
                      "incomplete 'Velocities' section", 1);
    GetBeadVelocity(System, file, *line_count);
  }
} //}}}
// assign velocity (if the line is correct) to a bead //{{{
static void GetBeadVelocity(SYSTEM *System, const char *file,
                            const int line_count) {
  long id;
  double vel[3];
  // error - wrong line //{{{
  if (words < 4 || !IsNaturalNumber(split[0], &id) || id > System->Count.Bead ||
      !IsRealNumber(split[1], &vel[0]) || // bead velocities
      !IsRealNumber(split[2], &vel[1]) || //
      !IsRealNumber(split[3], &vel[2])) { //
    err_msg("wrong line in Velocities section");
    PrintErrorFileLine(file, line_count);
    exit(1);
  } //}}}
  id--;
  for (int dd = 0; dd < 3; dd++) {
    System->Bead[id].Velocity[dd] = vel[dd];
  }
} //}}}
// read Atoms section //{{{
// mode: 0..full; 1..angle/bond; 2..atomic; 3..molecular; 4..charge
static void LmpDataReadAtoms(FILE *fr, const char *file, SYSTEM *System,
                             const BEADTYPE *name_mass,
                             int *line_count, const int mode) {
  COUNT *Count = &System->Count;
  int atom_types = Count->BeadType;
  // skip one line
  CountLineReadLine(line_count, fr, file, "incomplete 'Atoms' section", 1);
  bool warned = false; // to warn of undefined charge only once
  for (int i = 0; i < Count->Bead; i++) {
    CountLineReadLine(line_count, fr, file, "incomplete 'Atoms' section", 1);
    long id, resid, type;
    double pos[3], q;
    CheckAtomsLine(mode, Count->Bead, &id, &resid,
                   &type, &q, pos, file, *line_count);
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
          if (split[j][0] == '#' && (j + 1) < words) {
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
// read Bonds/Angles/Dihedrals/Impropers section //{{{
static void LmpDataReadBADISection(FILE *fr, const char *file,
                                   const COUNT Count, int (*array)[5],
                                   int *line_count, const int type) {
  // assign type-based variables //{{{
  int count, count_type, // number of bonds/angles/etc. and their types
      num; // number of required bead ids
  char txt[10];
  if (type == 0) {
    count = Count.Bond;
    count_type = Count.BondType;
    num = 2;
    s_strcpy(txt, "Bond", 10);
  } else if (type == 1) {
    count = Count.Angle;
    count_type = Count.AngleType;
    num = 3;
    s_strcpy(txt, "Angle", 10);
  } else if (type == 2) {
    count = Count.Dihedral;
    count_type = Count.DihedralType;
    num = 4;
    s_strcpy(txt, "Dihedral", 10);
  } else if (type == 3) {
    count = Count.Improper;
    count_type = Count.ImproperType;
    s_strcpy(txt, "Improper", 10);
    num = 4;
  } else {
    err_msg("LmpDataReadSection(): type must be 0 to 3");
    PrintError();
    exit(1);
  }
  char msg[LINE];
  snprintf(msg, LINE, "incomplete '%ss' section", txt); //}}}
  // does given id exists? For error when two bond/angle/etc. ids encountered
  bool *exists = calloc(count, sizeof *exists);
  // skip one line
  CountLineReadLine(line_count, fr, file, msg, 1);
  bool warned = false; // warn only once about '???' entry
  for (int i = 0; i < count; i++) {
    CountLineReadLine(line_count, fr, file, msg, 1);
    long id, type, b_id[num];
    // errors //{{{
    // line mus be <id> <type id> <bead id> <bead id> [<bead id>] [<bead id>]
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        (!IsNaturalNumber(split[1], &type) && strcmp(split[1], "???") != 0)) {
      snprintf(ERROR_MSG, LINE, "wrong line in '%ss' section", txt);
      goto error;
    }
    for (int j = 0; j < num; j++) {
      if (!IsNaturalNumber(split[2+j], &b_id[j])) {
        snprintf(ERROR_MSG, LINE, "wrong line in '%ss' section", txt);
        goto error;
      }
    }
    if (id > count) {
      snprintf(ERROR_MSG, LINE, "%s index is too high", txt);
      goto error;
    }
    if (type > count_type) {
      snprintf(ERROR_MSG, LINE, "%s type is too high", txt);
      goto error;
    }
    for (int j = 0; j < num; j++) {
      if (b_id[j] > Count.Bead) { // lammps data ids start from 1, so not >=
        snprintf(ERROR_MSG, LINE, "%s - bead index is too high", txt);
        goto error;
      }
    }
    for (int j = 0; j < (num - 1); j++) {
      for (int k = (j + 1); k < num; k++) {
        if (b_id[j] == b_id[k]) {
          snprintf(ERROR_MSG, LINE, "%s - identical bead ids", txt);
          goto error;
        }
      }
    }
    if (exists[id-1]) {
      snprintf(ERROR_MSG, LINE, "%s - repeated index", txt);
      goto error;
    } //}}}
    for (int j = 0; j < num; j++) {
      array[id-1][j] = b_id[j] - 1;
    }
    if (strcmp(split[1], "???") == 0) {
      if (!warned) {
        snprintf(ERROR_MSG, LINE, "undefined type (\'???\' in '%s' section)",
                 txt);
        PrintWarnFile(file, "\0", "\0");
        warned = true;
      }
      array[id-1][num] = -1;
    } else {
      array[id-1][num] = type - 1;
    }
    exists[id-1] = true;
  }
  free(exists);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
static bool ReadPbc(BOX *box) { //{{{
  double lo, hi;
  char *a[3] = {"xlo", "ylo", "zlo"};
  char *b[3] = {"xhi", "yhi", "zhi"};
  for (int i = 0; i < 3; i++) {
    if (words > 3 && strcmp(split[2], a[i]) == 0 &&
                     strcmp(split[3], b[i]) == 0) {
      if (!IsRealNumber(split[0], &lo) || !IsRealNumber(split[1], &hi)) {
        return false;
      }
      box->Low[i] = lo;
      box->OrthoLength[i] = hi - lo;
    }
  }
  // <double> <double> <double> xy xz yz
  if (words > 5 && strcmp(split[3], "xy") == 0 &&
             strcmp(split[4], "xz") == 0 && strcmp(split[5], "yz") == 0) {
    double xy, xz, yz;
    if (!IsRealNumber(split[0], &xy) || !IsRealNumber(split[1], &xz) ||
        !IsRealNumber(split[2], &yz)) {
      return false;
    }
    box->transform[0][1] = xy;
    box->transform[0][2] = xz;
    box->transform[1][2] = yz;
  }
  return true;
} //}}}
// check what mode is the Atoms section //{{{
static int CheckAtomsMode(const int words, char **split, const char *file) {
  int mode = 0;
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
  return mode;
} //}}}
// CheckAtomsLine() //{{{
static void CheckAtomsLine(const int mode, const int beads, long *id,
                           long *resid, long *type, double *q, double pos[3],
                           const char *file, const int line_count) {
  // bead index is always the first number
  if (words == 0 || !IsNaturalNumber(split[0], id) || *id > beads) {
    goto error;
  }
  // assign line position based on mode //{{{
  /*
   * vals[0]..minimum number of 'words'
   * Position on the line:
   * vals[1].. molecule id
   *      2 .. bead type
   *      3 .. charge
   *      4 .. x coordinate
  */
  int vals[5] = {-1, -1, -1, -1, -1};
  if (mode == 0) {
    // 'Atoms # full': <bead id> <mol id> <bead type id> <charge> <coordinates>
    vals[0] = 7; // minimum 'words'
    vals[1] = 1; // mol id
    vals[2] = 2; // bead type
    vals[3] = 3; // charge
    vals[4] = 4; // x coordinate
  } else if (mode == 1) {
    // 'Atoms # bond|angle|molecular':
    // <bead id> <mol id> <bead type id> <coordinates>
    vals[0] = 6; // minimum 'words'
    vals[1] = 1; // mol id
    vals[2] = 2; // bead type
    vals[4] = 3; // x coordinate
  } else if (mode == 2) {
    // 'Atoms # atomic': <bead id> <bead type id> <coordinates>
    vals[0] = 5; // minimum 'words'
    vals[2] = 1; // bead type
    vals[4] = 2; // x coordinate
  } else if (mode == 3) {
    // 'Atoms # charge': <bead id> <bead type id> <charge> <coordinates>
    vals[0] = 6; // minimum 'words'
    vals[2] = 1; // bead type
    vals[3] = 2; // charge
    vals[4] = 3; // x coordinate
  } else {
    err_msg("CheckAtomsLine(): mode must be 0 to 3");
    PrintError();
    exit(1);
  } //}}}
  // check the line is correct
  if (words < vals[0] ||
      (vals[1] != -1 && !IsIntegerNumber(split[vals[1]], resid)) ||
      (vals[2] != -1 && !IsWholeNumber(split[vals[2]], type)) ||
      (vals[3] != -1 && !IsRealNumber(split[vals[3]], q) &&
       strcmp(split[3], "???") != 0) ||
      (vals[4] != -1 && !IsRealNumber(split[vals[4]+0], &pos[0])) ||
      (vals[4] != -1 && !IsRealNumber(split[vals[4]+1], &pos[1])) ||
      (vals[4] != -1 && !IsRealNumber(split[vals[4]+2], &pos[2]))) {
    goto error;
  }
  if (mode == 1 || mode == 2) {
    *q = CHARGE;
  }
  if (mode == 2 || mode == 3) {
    *resid = -1;
  }
  return;
error:
  err_msg("wrong line in Atoms section");
  PrintErrorFileLine(file, line_count);
  exit(1);
} //}}}

// write files
void WriteAllStuff(FILE *fw, const SYSTEM System, const int type) { //{{{
  int num = 0;
  if (type == 0) { // bond
    num = 3;
    fprintf(fw, "\nBonds\n\n");
  } else if (type == 1) { // angle
    num = 4;
    fprintf(fw, "\nAngles\n\n");
  } else if (type == 2) { // dihedral
    num = 5;
    fprintf(fw, "\nDihedrals\n\n");
  } else /* if (type == 0) */ { // improper
    num = 5;
    fprintf(fw, "\nImpropers\n\n");
  }
  int count = 0;
  for (int i = 0; i < System.Count.Molecule; i++) {
    int mtype = System.Molecule[i].Type;
    MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
    int (**arr)[num];
    int *n;
    if (type == 0) { // bond
      arr = &mt_i->Bond;
      n = &mt_i->nBonds;
    } else if (type == 1) { // angle
      arr = &mt_i->Angle;
      n = &mt_i->nAngles;
    } else if (type == 2) { // dihedral
      arr = &mt_i->Dihedral;
      n = &mt_i->nDihedrals;
    } else /* if (type == 0) */ { // improper
      arr = &mt_i->Improper;
      n = &mt_i->nImpropers;
    }
    for (int j = 0; j < *n; j++) {
      WriteStuff(fw, System, i, &count, num, arr, j);
    }
  }
} //}}}

// write single bond/angle/dihedral/improper 'stuff' //{{{
static void WriteStuff(FILE *fw, const SYSTEM System, const int mol,
                       int *count, const int num, int (**arr)[5], const int n) {
  (*count)++;
  bool in = true;
  for (int aa = 0; aa < (num - 1); aa++) {
    int id = (*arr)[n][aa];
    id = System.Molecule[mol].Bead[id];
    if (!System.Bead[id].InTimestep) {
      in = false;
      break;
    }
  }
  if (in) {
    fprintf(fw, "%7d", *count);
    if ((*arr)[n][num-1] != -1) {
      fprintf(fw, " %6d", (*arr)[n][num-1] + 1);
    } else {
      fprintf(fw, "   ???");
    }
    for (int aa = 0; aa < (num - 1); aa++) {
      int id = (*arr)[n][aa];
      id = System.Molecule[mol].Bead[id];
      fprintf(fw, " %5d", id + 1);
    }
    putc('\n', fw);
  }
} //}}}
// write all bond/angle/dihedral/improper 'stuff'
// WriteLmpData() //{{{
void WriteLmpData(const SYSTEM System, const char *file, const bool mass,
                  const int argc, char **argv) {
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "Created via AnalysisTools v%s ", VERSION);
  fprintf(fw, "(https://github.com/KaGaSi/AnalysisTools); Command: ");
  PrintCommand(fw, argc, argv);
  putc('\n', fw);
  const COUNT *Count = &System.Count;
  // create new SYSTEM structure to figure out bead types if mass == true //{{{
  int mass_types = 0;
  int *bt_masstype_to_old = calloc(Count->BeadType, sizeof *bt_masstype_to_old);
  int *bt_old_to_masstype = calloc(Count->BeadType, sizeof *bt_old_to_masstype);
  if (mass) {
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt_i = &System.BeadType[i];
      bool found = false;
      for (int j = 0; j < mass_types; j++) {
        BEADTYPE *bt_j = &System.BeadType[bt_masstype_to_old[j]];
        if (bt_i->Mass == bt_j->Mass) {
          // bt_masstype_to_old[j] = i;
          bt_old_to_masstype[i] = j;
          found = true;
          break;
        }
      }
      if (!found) {
        bt_masstype_to_old[mass_types] = i;
        bt_old_to_masstype[i] = mass_types;
        mass_types++;
      }
    }
  } else {
    mass_types = Count->BeadType;
    for (int i = 0; i < Count->BeadType; i++) {
      bt_old_to_masstype[i] = i;
      bt_masstype_to_old[i] = i;
    }
  } //}}}
  // print counts //{{{
  fprintf(fw, "%10d atoms\n", Count->Bead);
  fprintf(fw, "%10d bonds\n", Count->Bond);
  fprintf(fw, "%10d angles\n", Count->Angle);
  fprintf(fw, "%10d dihedrals\n", Count->Dihedral);
  fprintf(fw, "%10d impropers\n", Count->Improper);
  putc('\n', fw);
  // // add one atom type for extra (possibly srp)
  // fprintf(fw, "%10d atom types\n", Count->BeadType + 1);
  fprintf(fw, "%10d atom types\n", mass_types);
  if (Count->Bond > 0 && Count->BondType > 0) {
    fprintf(fw, "%10d bond types\n", Count->BondType);
  }
  if (Count->Angle > 0 && Count->AngleType > 0) {
    fprintf(fw, "%10d angle types\n", Count->AngleType);
  }
  if (Count->Dihedral > 0 && Count->DihedralType > 0) {
    fprintf(fw, "%10d dihedral types\n", Count->DihedralType);
  }
  if (Count->Improper > 0 && Count->ImproperType > 0) {
    fprintf(fw, "%10d improper types\n", Count->ImproperType);
  }
  putc('\n', fw); //}}}
  // print box size //{{{
  fprintf(fw, "%.3f", System.Box.Low[0]);
  fprintf(fw, " %.3f xlo xhi\n", System.Box.Low[0] + System.Box.OrthoLength[0]);
  fprintf(fw, "%.3f", System.Box.Low[1]);
  fprintf(fw, " %.3f ylo yhi\n", System.Box.Low[1] + System.Box.OrthoLength[1]);
  fprintf(fw, "%.3f", System.Box.Low[2]);
  fprintf(fw, " %.3f zlo zhi\n", System.Box.Low[2] + System.Box.OrthoLength[2]);
  if (fabs(System.Box.alpha - 90) > 0.00001 ||
      fabs(System.Box.beta - 90) > 0.00001 ||
      fabs(System.Box.gamma - 90) > 0.00001) {
    fprintf(fw, "%.3f %.3f %.3f xy xz yz\n", System.Box.transform[0][1],
                                             System.Box.transform[0][2],
                                             System.Box.transform[1][2]);
  }
  putc('\n', fw); //}}}
  // print bead type masses //{{{
  fprintf(fw, "Masses\n\n");
  for (int i = 0; i < mass_types; i++) {
    BEADTYPE *bt = &System.BeadType[bt_masstype_to_old[i]];
    fprintf(fw, "%5d", i + 1);
    if (bt->Mass == MASS) {
      fprintf(fw, " ???");
    } else {
      fprintf(fw, " %lf", bt->Mass);
    }
    fprintf(fw, " # %s\n", bt->Name);
  } //}}}
  // print various coeffs //{{{
  if (Count->BondType > 0) {
    fprintf(fw, "\nBond Coeffs\n\n");
    for (int i = 0; i < Count->BondType; i++) {
      fprintf(fw, "%5d %lf %lf\n", i + 1, System.BondType[i].a / 2,
                                          System.BondType[i].b);
    }
  }
  if (Count->AngleType > 0) {
    fprintf(fw, "\nAngle Coeffs\n\n");
    for (int i = 0; i < Count->AngleType; i++) {
      fprintf(fw, "%5d %lf %lf\n", i + 1, System.AngleType[i].a / 2,
                                          System.AngleType[i].b);
    }
  }
  if (Count->DihedralType > 0) {
    fprintf(fw, "\nDihedral Coeffs\n\n");
    for (int i = 0; i < Count->DihedralType; i++) {
      fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.DihedralType[i].a / 2,
                                              System.DihedralType[i].b,
                                              System.DihedralType[i].c);
    }
  }
  if (Count->ImproperType > 0) {
    fprintf(fw, "\nImproper Coeffs\n\n");
    for (int i = 0; i < Count->ImproperType; i++) {
      fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.ImproperType[i].a / 2,
                                              System.ImproperType[i].b,
                                              System.ImproperType[i].c);
    }
  } //}}}
  // print atoms //{{{
  // if there is 0 molecule index, saved indices will get +1
  // TODO: why would I need to go from 1?
  bool zero = false;
  // for (int i = 0; i < Count->Molecule; i++) {
  //   if (System.Molecule[i].Index == 0) {
  //     zero = true;
  //     break;
  //   }
  // }
  fprintf(fw, "\nAtoms # full\n\n");
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    // <bead id>
    fprintf(fw, "%7d", id + 1);
    // <molecule id (-1 for no molecule)>
    int mol = bead->Molecule;
    if (mol != -1) {
      int id = System.Molecule[mol].Index;
      if (zero) {
        id++;
      }
      fprintf(fw, " %5d", id);
    } else {
      fprintf(fw, " %5d", -1);
    }
    // <bead type id>
    fprintf(fw, " %5d", bt_old_to_masstype[bead->Type] + 1);
    // <charge> from original System as the charge can differ in lmp data file
    int type = bead->Type;
    double q = System.BeadType[type].Charge;
    if (q == CHARGE) {
      fprintf(fw, " %15s", "???");
    } else {
      fprintf(fw, " %15f", q);
    }
    // coordinates
    for (int dd = 0; dd < 3; dd++) {
      fprintf(fw, " %15f", bead->Position[dd] + System.Box.Low[dd]);
    }
    // molecule name
    if (mol != -1) {
      int type = System.Molecule[mol].Type;
      fprintf(fw, " # %s", System.MoleculeType[type].Name);
    }
    putc('\n', fw);
  } //}}}
  // print velocities (if at least one non-zero) //{{{
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    double (*vel)[3] = &System.Bead[id].Velocity;
    if (fabs((*vel)[0]) > 1e-5 ||
        fabs((*vel)[1]) > 1e-5 ||
        fabs((*vel)[2]) > 1e-5) {
      fprintf(fw, "\nVelocities\n\n");
      for (int j = 0; j < Count->BeadCoor; j++) {
        id = System.BeadCoor[j];
        vel = &System.Bead[id].Velocity;
        fprintf(fw, "%7d", id + 1);
        for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %15f", (*vel)[dd]);
        }
        putc('\n', fw);
      }
      break;
    }
  } //}}}
  // print bonds/angles/dihedrals/impropers
  if (Count->Bond > 0) {
    WriteAllStuff(fw, System, 0);
  }
  if (Count->Angle > 0) {
    WriteAllStuff(fw, System, 1);
  }
  if (Count->Dihedral > 0) {
    WriteAllStuff(fw, System, 2);
  }
  if (Count->Improper > 0) {
    WriteAllStuff(fw, System, 3);
  }
  free(bt_masstype_to_old);
  free(bt_old_to_masstype);
  fclose(fw);
} //}}}
