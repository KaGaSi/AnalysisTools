#include "AnalysisTools.h"
#include "ReadWriteGromacs.h"
#include "Errors.h"
#include "General.h"
#include "ReadWrite.h"
#include "Structs.h"
#include "System.h"
#include <strings.h>

static bool CountLineReadLine(int *line_count, FILE *fr,
                              const char *file, char *msg);
static bool ReadMoleculeType(const char *file, FILE *fr,
                             int *line_count, SYSTEM *System);
static void ReadAtoms(const char *file, FILE *fr,
                      int *line_count, SYSTEM *System);
static void ReadBonds(const char *file, FILE *fr, int *line_count,
                      const int mtype, SYSTEM *System);

SYSTEM GromacsReadStruct(const char *file) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  FILE *fr = OpenFile(file, "r");
  int line_count = 0;
  while (true) {
    if (!ReadMoleculeType(file, fr, &line_count, &System)) {
      break;
    }
    ReadAtoms(file, fr, &line_count, &System);
    do {
      if (!CountLineReadLine(&line_count, fr, file, "")) {
        line_count--;
        break;
      }
      // [ bonds ] section
      if ((words > 0 && strcasecmp(split[0], "[bonds]") == 0) ||
          (words > 2 && strcasecmp(split[1], "bonds") == 0)) {
        ReadBonds(file, fr, &line_count, Count->MoleculeType - 1, &System);
      }
    } while (true);
    printf("%d\n", line_count);
    printf("%d %s\n", words, split[0]);
    System.Molecule = realloc(System.Molecule,
                              (Count->Molecule + 1) * sizeof *System.Molecule);
    MOLECULE *mol = &System.Molecule[Count->Molecule];
    MOLECULETYPE *mt = &System.MoleculeType[Count->MoleculeType-1];
    mol->Index = Count->Molecule;
    mol->Type = Count->MoleculeType - 1;
    mol->Bead = malloc(mt->nBeads * sizeof *mol->Bead);
    for (int i = 0; i < mt->nBeads; i++) {
      mol->Bead[Count->Bead-mt->nBeads+i] = Count->Bead - mt->nBeads + i;
    }
    Count->HighestResid = Count->Molecule;
    Count->Molecule++;
  }
  fclose(fr);

  Count->Bonded = Count->Bead;
  Count->BeadCoor = Count->Bead;
  ReallocBead(&System);
  for (int i = 0; i < Count->Bead; i++) {
    System.BeadCoor[i] = i;
  }
  FillSystemNonessentials(&System, true);

  printf("%d\n", line_count);

  CheckSystem(System, file);

  PrintCount(*Count);
  PrintBeadType(System);
  PrintAllMolTypes(System);

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
// static void ReadMoleculeType() //{{{
static bool ReadMoleculeType(const char *file, FILE *fr,
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
    CountLineReadLine(line_count, fr, file, "reading [ moleculetype ] section");
  } while (words == 0 || split[0][0] == ';');
  // read the atoms lines
  while (true) {
    if (words < 8) { // here, because 1st atoms line was already read
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
    // what molecule the bead is in?
    int m_id = FindMoleculeName(split[3], *System);
    if (m_id == -1) {
      err_msg("atom not in any known molecule");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    MOLECULETYPE *mt = &System->MoleculeType[m_id];
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
    CountLineReadLine(line_count, fr, file, "reading [ atoms ] section");
  };
} //}}}
// static void ReadBonds() //{{{
static void ReadBonds(const char *file, FILE *fr, int *line_count,
                      const int mtype, SYSTEM *System) {
  COUNT *Count = &System->Count;
  MOLECULETYPE *mt = &System->MoleculeType[mtype];
  // skip comments (;-starting lines) and blank lines and read first bonds line
  do {
    CountLineReadLine(line_count, fr, file,
                      "reading [ moleculetype ] section");
  } while (words == 0 || split[0][0] == ';');
  // read the bonds lines
  do {
    long int bond[2];
    if (words < 2 || !IsNaturalNumber(split[0], &bond[0]) ||
                     !IsNaturalNumber(split[1], &bond[1])) {
      break;
    }
    if (bond[0] > mt->nBeads || bond[1] > mt->nBeads) {
      err_msg("atom id in a bond is too high");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    mt->nBonds++;
    mt->Bond = realloc(mt->Bond, mt->nBonds * sizeof *mt->Bond);
    mt->Bond[mt->nBonds-1][0] = bond[0] - 1;
    mt->Bond[mt->nBonds-1][1] = bond[1] - 1;
    mt->Bond[mt->nBonds-1][2] = -1;
  } while (CountLineReadLine(line_count, fr, file, ""));
} //}}}
