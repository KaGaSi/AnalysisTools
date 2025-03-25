#include "AnalysisTools.h"
#include "ReadWriteGromacs.h"
#include "Errors.h"
#include "General.h"
#include "Globals.h"
#include "ReadWrite.h"
#include "Structs.h"
#include "System.h"

static bool CountLineReadLine(int *line_count, FILE *fr,
                              const char *file, char *msg);
static bool ReadMoleculetype(const char *file, FILE *fr,
                             int *line_count, SYSTEM *System);
static void ReadAtoms(const char *file, FILE *fr,
                      int *line_count, SYSTEM *System);
static void ReadBonds(const char *file, FILE *fr,
                      int *line_count, SYSTEM *System);
static void ReadAngles(const char *file, FILE *fr,
                      int *line_count, SYSTEM *System);
static void ReadDihedrals(const char *file, FILE *fr,
                          int *line_count, SYSTEM *System);
// not used for now
static void ReadBADI(const char *file, FILE *fr, int *line_count,
                     SYSTEM *System, int (**arr)[5], const int type);
static SYSTEM CountInFile(const char *file);
static void CountOneSection(const char *file, FILE *fr, int *line_count,
                            int *count, int min_num);

// read files
SYSTEM ItpReadStruct(const char *file) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  FILE *fr = OpenFile(file, "r");
  int line_count = 0;
  while (true) { // read multiple molecules if present
    if (!ReadMoleculetype(file, fr, &line_count, &System)) {
      break;
    }
    ReadAtoms(file, fr, &line_count, &System);
    do {
      fpos_t saved_pos;
      fgetpos(fr, &saved_pos);
      if (!CountLineReadLine(&line_count, fr, file, "")) {
        line_count--;
        break;
      }
      if ((words > 0 && strcasecmp(split[0], "[bonds]") == 0) ||
          (words > 2 && strcasecmp(split[1], "bonds") == 0)) {
        ReadBonds(file, fr, &line_count, &System);
      } else if ((words > 0 && strcasecmp(split[0], "[angles]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "angles") == 0)) {
        ReadAngles(file, fr, &line_count, &System);
      } else if ((words > 0 && strcasecmp(split[0], "[dihedrals]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "dihedrals") == 0)) {
        ReadDihedrals(file, fr, &line_count, &System);
      } else if ((words > 0 && strcasecmp(split[0], "[moleculetype]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "moleculetype") == 0)) {
        fsetpos(fr, &saved_pos);
        line_count--;
        break;
      }
    } while (true);
    Count->Molecule++;
    System.Molecule = realloc(System.Molecule,
                              Count->Molecule * sizeof *System.Molecule);
    MOLECULE *mol = &System.Molecule[Count->Molecule-1];
    MOLECULETYPE *mt = &System.MoleculeType[Count->MoleculeType-1];
    mol->Index = Count->Molecule - 1;
    mol->Type = Count->MoleculeType - 1;
    mol->Bead = malloc(mt->nBeads * sizeof *mol->Bead);
    mol->InTimestep = true;
    for (int i = 0; i < mt->nBeads; i++) {
      int id = Count->Bead - mt->nBeads + i;
      mol->Bead[i] = id;
    }
    Count->HighestResid = Count->Molecule;
  }
  fclose(fr);

  Count->Bonded = Count->Bead;
  Count->BeadCoor = Count->Bead;
  ReallocBead(&System);
  for (int i = 0; i < Count->Bead; i++) {
    System.BeadCoor[i] = i;
  }
  if (Count->Molecule > 0) {
    System.MoleculeCoor = s_realloc(System.MoleculeCoor, Count->Molecule *
                                    sizeof *System.MoleculeCoor);
  }
  FillSystemNonessentials(&System, true);

  CheckSystem(System, file);

  return System;
} //}}}
SYSTEM PdbReadStruct(const char *file) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // 1) find pbc - 'CRYST*' line; assume orthogonal box //{{{
  while (CountLineReadLine(&line_count, fr, file, "missing CRYST (pbc) line")) {
    line_count++;
    if (words > 0 && strncasecmp(split[0], "CRYST", 5) == 0) {
      if (words < 4 || !IsPosRealNumber(split[1], &System.Box.Length[0]) ||
                       !IsPosRealNumber(split[1], &System.Box.Length[1]) ||
                       !IsPosRealNumber(split[1], &System.Box.Length[2])) {
        err_msg("wrong 'CRYST' line (requires CRYST* <x> <y> <z>)");
        PrintErrorFileLine(file, line_count);
        exit(1);
      }
      CalculateBoxData(&System.Box, 0);
      break;
    }
  } //}}}
  // 2) find first ATOM line //{{{
  while (CountLineReadLine(&line_count, fr, file, "missing Atom line(s)")) {
    line_count++;
    if (words > 0 && strcasecmp(split[0], "ATOM") == 0) {
      break;
    }
  } //}}}
  // 3) read ATOM lines - assume no blanks or commnets //{{{
  // 1st ATOM was already read
  do {
    if (words < 1 || strcasecmp(split[0], "ATOM") != 0) {
      break;
    }
    long b_id, m_id;
    double coor[3];
    if (words < 9 || !IsNaturalNumber(split[1], &b_id) ||
                     !IsNaturalNumber(split[5], &m_id) ||
                     !IsRealNumber(split[6], &coor[0]) ||
                     !IsRealNumber(split[7], &coor[1]) ||
                     !IsRealNumber(split[8], &coor[2])) {
      err_msg("wrong ATOM line");
      PrintErrorFileLine(file, line_count);
      exit(1);
    }
    System.Bead = realloc(System.Bead, b_id * sizeof *System.Bead);
    // molecule and bead id's start at 1 in pdb
    m_id--;
    b_id--;
    // new bead type?
    int bt_id = FindBeadType(split[2], System);
    if (bt_id == -1) {
      bt_id = Count->BeadType;
      NewBeadType(&System.BeadType, &Count->BeadType,
                  split[2], CHARGE, MASS, RADIUS);
    }
    // new molecule type?
    int mt_id = FindMoleculeName(split[3], System);
    if (mt_id == -1) {
      mt_id = Count->MoleculeType;
      NewMolType(&System.MoleculeType, &Count->MoleculeType,
                 split[3], 1, 0, 0, 0, 0);
      System.MoleculeType[mt_id].nBeads = 0;
      System.MoleculeType[mt_id].Number = 0;
    }
    MOLECULETYPE *mt = &System.MoleculeType[mt_id];
    BEADTYPE *bt = &System.BeadType[bt_id];
    BEAD *b = &System.Bead[b_id];
    // new molecule?
    if (m_id > (Count->Molecule - 1)) {
      mt->Number++;
    }
    bt->Number++;
    b->Molecule = m_id;
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = coor[dd];
    }
    b->InTimestep = true;
    b->Type = bt_id;
    if (mt->Number == 1) {
      if (mt->nBeads > 0) {
        mt->Bead = realloc(mt->Bead, (mt->nBeads + 1) * sizeof *mt->Bead);
      }
      mt->Bead[mt->nBeads] = bt_id;
      mt->nBeads++;
    }
    Count->Bead++;
    Count->Molecule = m_id + 1;

    line_count++;
  } while (CountLineReadLine(&line_count, fr, file, "")); //}}}
  fclose(fr);

  // allocatate Molecule & Molecule[].Bead arrays
  System.Molecule = realloc(System.Molecule,
                            Count->Molecule * sizeof *System.Molecule);
  int m_id = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    for (int j = 0; j < mt->Number; j++) {
      MOLECULE *mol = &System.Molecule[m_id];
      mol->Bead = calloc(mt->nBeads, sizeof *mol->Bead);
      mol->Type = i;
      mol->InTimestep = true;
      mol->Index = m_id;
      m_id++;
    }
  }

  int *per_mol = calloc(Count->Molecule, sizeof *per_mol);
  for (int i = 0; i < Count->Bead; i++) {
    m_id = System.Bead[i].Molecule;
    MOLECULE *mol = &System.Molecule[System.Bead[i].Molecule];
    mol->Bead[per_mol[m_id]] = i;
    per_mol[m_id]++;
  }
  free(per_mol);

  Count->HighestResid = Count->Molecule;

  Count->Bonded = Count->Bead;
  Count->BeadCoor = Count->Bead;
  ReallocBead(&System);
  for (int i = 0; i < Count->Bead; i++) {
    System.BeadCoor[i] = i;
  }
  if (Count->Molecule > 0) {
    System.MoleculeCoor = s_realloc(System.MoleculeCoor, Count->Molecule *
                                    sizeof *System.MoleculeCoor);
  }
  FillSystemNonessentials(&System, true);

  CheckSystem(System, file);

  // PrintCount(*Count);
  // PrintBeadType(System);
  // PrintAllMolTypes(System);
  // PrintMolecules(System);
  // PrintBead(System);
  return System;
} //}}}

// increment line count and read the line //{{{
// TODO: mode?
static bool CountLineReadLine(int *line_count, FILE *fr,
                              const char *file, char *msg) {
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (msg[0] != '\0') {
      ErrorEOF(file, msg);
      exit(1);
    } else {
      return false;
    }
  }
  return true;
} //}}}
// static void ReadMoleculetype() //{{{
static bool ReadMoleculetype(const char *file, FILE *fr,
                             int *line_count, SYSTEM *System) {
  COUNT *Count = &System->Count;
  // find [ moleculetype ]
  do {
    if (!CountLineReadLine(line_count, fr, file, "")) {
      if (Count->MoleculeType > 0) {
        (*line_count)--;
        return false;
      } else {
        ErrorEOF(file, "searching for [ moleculetype ] keyword");
        exit(1);
      }
    }
  } while ((words < 1 || strcasecmp(split[0], "[moleculetype]") != 0) &&
           (words < 3 || strcasecmp(split[1], "moleculetype") != 0));
  // skip comments (;-starting lines) and blank lines and read moleculetype line
  do {
    CountLineReadLine(line_count, fr, file, "reading [ moleculetype ] section");
  } while (words == 0 || split[0][0] == ';');
  // must be <name> <number>, but I don't use <number>...
  if (words < 1) {
    err_msg("wrong [ moleculetype ] line");
    PrintErrorFileLine(file, *line_count);
    exit(1);
  }
  // new molecule type with no beads, bonds, etc.
  NewMolType(&System->MoleculeType, &Count->MoleculeType, split[0],
             0, 0, 0, 0, 0);
  return true;
} //}}}
// static void ReadAtoms() //{{{
static void ReadAtoms(const char *file, FILE *fr,
                      int *line_count, SYSTEM *System) {
  COUNT *Count = &System->Count;
  // find [ atoms ]
  do {
    CountLineReadLine(line_count, fr, file, "searching for [ atoms ] keyword");
  } while ((words < 1 || strcasecmp(split[0], "[atoms]") != 0) &&
           (words < 3 || strcasecmp(split[1], "atoms") != 0));
  // skip comments (;-starting lines) and blank lines and read first atoms line
  do {
    CountLineReadLine(line_count, fr, file, "reading [ atoms ] section");
  } while (words == 0 || split[0][0] == ';');
  // read the atoms lines
  fpos_t saved_pos;
  fgetpos(fr, &saved_pos);
  while (true) {
    if (words < 8) { // here, because 1st atoms line was already read
      fsetpos(fr, &saved_pos);
      break;
    }
    // assume 2nd column is name - may change
    char *a_name = split[4];
    // get mass //{{{
    double m;
    if (!IsRealNumber(split[7], &m) || m <= 0) {
      err_msg("8th column should be mass");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    // add new bead type?
    int b_id = FindBeadType(a_name, *System);
    if (b_id == -1) {
      NewBeadType(&System->BeadType, &Count->BeadType,
                  a_name, CHARGE, m, RADIUS);
      b_id = Count->BeadType - 1;
    }
    System->BeadType[b_id].Number++;
    MOLECULETYPE *mt = &System->MoleculeType[Count->MoleculeType-1];
    mt->nBeads++;
    mt->Bead = s_realloc(mt->Bead, mt->nBeads * sizeof *mt->Bead);
    mt->Bead[mt->nBeads-1] = b_id;

    Count->Bead++;
    // realloc & fill Bead struct
    System->Bead = realloc(System->Bead, Count->Bead * sizeof *System->Bead);
    BEAD *bead = &System->Bead[Count->Bead-1];
    bead->Molecule = Count->MoleculeType - 1;
    bead->Type = b_id;
    bead->InTimestep = true;
    // read next line
    fgetpos(fr, &saved_pos);
    CountLineReadLine(line_count, fr, file, "reading [ atoms ] section");
  };
} //}}}
// static void ReadBonds() //{{{
static void ReadBonds(const char *file, FILE *fr,
                      int *line_count, SYSTEM *System) {
  COUNT *Count = &System->Count;
  MOLECULETYPE *mt = &System->MoleculeType[Count->MoleculeType-1];
  // skip comments (;-starting lines) and blank lines and read first bonds line
  do {
    CountLineReadLine(line_count, fr, file,
                      "reading [ bonds ] section");
  } while (words == 0 || split[0][0] == ';');
  // read the bonds lines
  fpos_t saved_pos;
  fgetpos(fr, &saved_pos);
  do {
    long int bond[2];
    if (words < 2 || !IsNaturalNumber(split[0], &bond[0]) ||
                     !IsNaturalNumber(split[1], &bond[1])) {
      fsetpos(fr, &saved_pos);
      (*line_count)--;
      break;
    }
    if (bond[0] > mt->nBeads || bond[1] > mt->nBeads) {
      err_msg("atom id in a bond is too high");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    mt->nBonds++;
    if (mt->nBonds == 1) {
      mt->Bond = malloc(sizeof *mt->Bond);
    } else {
      mt->Bond = realloc(mt->Bond, mt->nBonds * sizeof *mt->Bond);
    }
    mt->Bond[mt->nBonds-1][0] = bond[0] - 1;
    mt->Bond[mt->nBonds-1][1] = bond[1] - 1;
    mt->Bond[mt->nBonds-1][2] = -1;
    fgetpos(fr, &saved_pos);
  } while (CountLineReadLine(line_count, fr, file, ""));
} //}}}
// static void ReadAngles() //{{{
static void ReadAngles(const char *file, FILE *fr,
                       int *line_count, SYSTEM *System) {
  COUNT *Count = &System->Count;
  MOLECULETYPE *mt = &System->MoleculeType[Count->MoleculeType-1];
  // skip comments (;-starting lines) and blank lines and read first angles line
  do {
    CountLineReadLine(line_count, fr, file,
                      "reading [ angles ] section");
  } while (words == 0 || split[0][0] == ';');
  // read the angles lines
  fpos_t saved_pos;
  fgetpos(fr, &saved_pos);
  do {
    long int angle[3];
    if (words < 3 || !IsNaturalNumber(split[0], &angle[0]) ||
                     !IsNaturalNumber(split[1], &angle[1]) ||
                     !IsNaturalNumber(split[2], &angle[2])) {
      fsetpos(fr, &saved_pos);
      (*line_count)--;
      break;
    }
    if (angle[0] > mt->nBeads ||
        angle[1] > mt->nBeads ||
        angle[2] > mt->nBeads) {
      err_msg("atom id in a angle is too high");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    mt->nAngles++;
    if (mt->nAngles == 1) {
      mt->Angle = malloc(sizeof *mt->Angle);
    } else {
      mt->Angle = realloc(mt->Angle, mt->nAngles * sizeof *mt->Angle);
    }
    mt->Angle[mt->nAngles-1][0] = angle[0] - 1;
    mt->Angle[mt->nAngles-1][1] = angle[1] - 1;
    mt->Angle[mt->nAngles-1][2] = angle[2] - 1;
    mt->Angle[mt->nAngles-1][3] = -1;
    fgetpos(fr, &saved_pos);
  } while (CountLineReadLine(line_count, fr, file, ""));
} //}}}
// static void ReadDihedrals() //{{{
static void ReadDihedrals(const char *file, FILE *fr,
                          int *line_count, SYSTEM *System) {
  COUNT *Count = &System->Count;
  MOLECULETYPE *mt = &System->MoleculeType[Count->MoleculeType-1];
  // skip comments (;-starting lines) and blank lines and read first dihds line
  do {
    CountLineReadLine(line_count, fr, file,
                      "reading [ dihedrals ] section");
  } while (words == 0 || split[0][0] == ';');
  // read the dihedrals lines
  fpos_t saved_pos;
  fgetpos(fr, &saved_pos);
  do {
    long int dihed[4];
    if (words < 4 || !IsNaturalNumber(split[0], &dihed[0]) ||
                     !IsNaturalNumber(split[1], &dihed[1]) ||
                     !IsNaturalNumber(split[2], &dihed[2]) ||
                     !IsNaturalNumber(split[3], &dihed[3])) {
      fsetpos(fr, &saved_pos);
      (*line_count)--;
      break;
    }
    if (dihed[0] > mt->nBeads ||
        dihed[1] > mt->nBeads ||
        dihed[2] > mt->nBeads) {
      err_msg("atom id in a dihedral is too high");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    mt->nDihedrals++;
    if (mt->nDihedrals == 1) {
      mt->Dihedral = malloc(sizeof *mt->Dihedral);
    } else {
      mt->Dihedral = realloc(mt->Dihedral,
                             mt->nDihedrals * sizeof *mt->Dihedral);
    }
    mt->Dihedral[mt->nDihedrals-1][0] = dihed[0] - 1;
    mt->Dihedral[mt->nDihedrals-1][1] = dihed[1] - 1;
    mt->Dihedral[mt->nDihedrals-1][2] = dihed[2] - 1;
    mt->Dihedral[mt->nDihedrals-1][3] = dihed[3] - 1;
    mt->Dihedral[mt->nDihedrals-1][4] = -1;
    fgetpos(fr, &saved_pos);
  } while (CountLineReadLine(line_count, fr, file, ""));
} //}}}

// unused - for now
// static void ReadBADI() //{{{
static void ReadBADI(const char *file, FILE *fr, int *line_count,
                     SYSTEM *System, int (**arr)[5], const int type) {
  COUNT *Count = &System->Count;
  // assign type-based variables //{{{
  int count, // number of bonds/angles/etc. and their types
      num; // number of required bead ids
  char txt[10];
  if (type == 0) {
    count = Count->Bond;
    num = 2;
    s_strcpy(txt, "bonds", 10);
  } else if (type == 1) {
    count = Count->Angle;
    num = 3;
    s_strcpy(txt, "angles", 10);
  } else if (type == 2) {
    count = Count->Dihedral;
    num = 4;
    s_strcpy(txt, "dihedrals", 10);
  } else if (type == 3) {
    count = Count->Improper;
    s_strcpy(txt, "impropers", 10);
    num = 4;
  } else {
    err_msg("ReadBADI(): type must be 0 to 3");
    PrintError();
    exit(1);
  }
  char msg[LINE];
  snprintf(msg, LINE, "incomplete '[ %ss ]' section", txt); //}}}
  MOLECULETYPE *mt = &System->MoleculeType[Count->MoleculeType-1];
  // skip comments (;-starting lines) and blank lines and read first bonds line
  do {
    CountLineReadLine(line_count, fr, file, msg);
  } while (words == 0 || split[0][0] == ';');
  // read the bonds lines
  do {
    long int id[num];
    if (words < num) {
      break;
    }
    if (!IsNaturalNumber(split[0], &id[0]) ||
        !IsNaturalNumber(split[1], &id[1])) {
      break;
    }
    if (id[0] > mt->nBeads || id[1] > mt->nBeads) {
      err_msg("atom id is too high");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    mt->nBonds++;
    if (mt->nBonds == 1) {
      mt->Bond = malloc(sizeof *mt->Bond);
    } else {
      mt->Bond = realloc(mt->Bond, mt->nBonds * sizeof *mt->Bond);
    }
    mt->Bond[mt->nBonds-1][0] = id[0] - 1;
    mt->Bond[mt->nBonds-1][1] = id[1] - 1;
    mt->Bond[mt->nBonds-1][2] = -1;
  } while (CountLineReadLine(line_count, fr, file, ""));
} //}}}
// static void CountOneSection() //{{{
static void CountOneSection(const char *file, FILE *fr, int *line_count,
                            int *count, int min_num) {
  while (CountLineReadLine(line_count, fr, file, "")) {
    fpos_t saved_pos;
    // only check non-blank, non-comment lines
    if (words > 0 && split[0][0] != ';') {
      if (split[0][0] == '[') { // exit loop on [ <section> ]
        fsetpos(fr, &saved_pos); // return pointer to the previous line
        (*line_count)--; // pointer rewound one line, so decrement count
        break;
      } else if (words < min_num) { // but too few words
        // correct section name //{{{
        int length = 20;
        char txt[length];
        switch (min_num) {
          case 8:
            s_strcpy(txt, "atoms", length);
            break;
          case 3:
            s_strcpy(txt, "bonds", length);
            break;
          case 4:
            s_strcpy(txt, "angles", length);
            break;
          case 5:
            s_strcpy(txt, "dihedrals/impropers", length);
            break;
        } //}}}
        snprintf(ERROR_MSG, LINE, "wrong line in [ %s ] section", txt);
        PrintErrorFileLine(file, *line_count);
        exit(1);
      } else { // correct line
        (*count)++;
      }
    }
    // save position in file
    fgetpos(fr, &saved_pos);
  }
} //}}}
static SYSTEM CountInFile(const char *file) { //{{{
  int line_count = 0;
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  FILE *fr = OpenFile(file, "r");
  // read first line, erroring out on empty file
  while (CountLineReadLine(&line_count, fr, file, "")) {
    // if not new molecule type, read next line
    if ((words < 1 || strcasecmp(split[0], "[moleculetype]") != 0) &&
        (words < 3 || strcasecmp(split[1], "moleculetype") != 0)) {
      continue;
    }
    // Count->MoleculeType++;
    char name[MOL_NAME];
    while (CountLineReadLine(&line_count, fr, file,
                             "missing  moleculetype name")) {
      if (words > 0) {
        s_strcpy(name, split[0], MOL_NAME);
        break;
      }
    }
    // if new molecule type, search for [ atoms ] //{{{
    bool section = false;
    while (CountLineReadLine(&line_count, fr, file, "")) {
      // if not [ atoms ] line, read next line
      if ((words > 0 && strcasecmp(split[0], "[atoms]") == 0) ||
          (words > 2 && strcasecmp(split[1], "atoms") == 0)) {
        section = true;
        break;
      }
      // found [ moleculetype ] before [ atoms ] - error
      if ((words > 0 && strcasecmp(split[0], "[moleculetype]") == 0) ||
          (words > 2 && strcasecmp(split[1], "moleculetype") == 0)) {
        break;
      }
    }
    if (!section) {
      err_msg("missing [ atoms ] section");
      PrintErrorFile(file, "\0", "\0");
      exit(1);
    } //}}}
    int beads = 0,
        bonds = 0,
        angles = 0,
        dihedrals = 0,
        impropers = 0;
    CountOneSection(file, fr, &line_count, &beads, 8);
    fpos_t saved_pos;
    // count bonds/angles/dihedrals/impropers
    while (CountLineReadLine(&line_count, fr, file, "")) {
      if ((words > 0 && strcasecmp(split[0], "[bonds]") == 0) ||
          (words > 2 && strcasecmp(split[1], "bonds") == 0)) {
        CountOneSection(file, fr, &line_count, &bonds, 3);
      } else if ((words > 0 && strcasecmp(split[0], "[angles]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "angles") == 0)) {
        CountOneSection(file, fr, &line_count, &angles, 4);
      } else if ((words > 0 && strcasecmp(split[0], "[dihedrals]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "dihedrals") == 0)) {
        CountOneSection(file, fr, &line_count, &dihedrals, 5);
      } else if ((words > 0 && strcasecmp(split[0], "[impropers]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "impropers") == 0)) {
        CountOneSection(file, fr, &line_count, &impropers, 5);
      } else if ((words > 0 && strcasecmp(split[0], "[moleculetype]") == 0) ||
                 (words > 2 && strcasecmp(split[1], "moleculetype") == 0)) {
        // start of the next molecule
        break;
      }
      fgetpos(fr, &saved_pos);
    }
    fsetpos(fr, &saved_pos);
    line_count--;
    System.MoleculeType = realloc(System.MoleculeType, Count->MoleculeType *
                                  sizeof *System.MoleculeType);
    NewMolType(&System.MoleculeType, &Count->MoleculeType, name, beads, bonds,
               angles, dihedrals, 0);
    MOLECULETYPE *mt = &System.MoleculeType[Count->MoleculeType-1];
    for (int i = 0; i < mt->nBeads; i++) {
      mt->Bead[i] = Count->Bead + i;
    }

    line_count--;
    Count->Bead += beads;
    Count->Bond += bonds;
    Count->Angle += angles;
    Count->Dihedral += dihedrals;
    Count->Improper += impropers;
    printf("Molecule %d: %d %d %d %d %d\n",
           Count->MoleculeType, beads, bonds, angles, dihedrals, impropers);
  }
  if (Count->MoleculeType == 0) {
    line_count--;
    ErrorEOF(file, "missing [ moleculetype ] section");
    exit(1);
  }
  fclose(fr);

  return System;
} //}}}

// write file
