#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
This utility reads a library of molecules/beads and produces FIELD file from a \
different text file.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <lib dir> <out1.FIELD> "
          "<in.data> <in.lmp> <out2.FIELD>\n\n", cmd);
  fprintf(ptr, "   <input>      input text file defining the system\n");
  fprintf(ptr, "   <lib dir>    path to the library directory\n");
  fprintf(ptr, "   <out1.FIELD> output FIELD with species to add\n");
  fprintf(ptr, "   <in.data>    LAMMPS data file\n");
  fprintf(ptr, "   <in.lmp>     LAMMPS input script\n");
  fprintf(ptr, "   <out2.FIELD> output FIELD from <in.data> and <in.lmp>\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 4, all = common + 0, count = 0, req_arg = 6;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options

  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count arguments

  // command line options //{{{
  // <input> - file defining the system
  char in_file[LINE] = "";
  snprintf(in_file, LINE, "%s", argv[++count]);
  // <lib dir> - path to directory 
  char lib_dir[LINE] = "";
  snprintf(lib_dir, LINE, "%s", argv[++count]);
  // <out1.FIELD> - output field file
  char add_field[LINE] = "";
  snprintf(add_field, LINE, "%s", argv[++count]);
  int add_type = FIELD_FILE;
  // <in.data>
  char in_data[LINE] = "";
  snprintf(in_data, LINE, "%s", argv[++count]);
  int in_data_type = LDATA_FILE;
  // <in.lmp>
  char in_lmp[LINE] = "";
  snprintf(in_lmp, LINE, "%s", argv[++count]);
  // <out2.FIELD> - output field file with original system
  char orig_field[LINE] = "";
  snprintf(orig_field, LINE, "%s", argv[++count]);
  int orig_type = FIELD_FILE;
  // read options
  bool silent, verbose, rubbish;
  int trash[1]; // some stuff for unused things in options
  CommonOptions(argc, argv, LINE, &verbose, &silent, &rubbish, &rubbish,
                trash, trash, trash, trash); //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // read input file with system description //{{{
  char line[LINE], *split[SPL_STR];
  int words, line_count = 0;
  FILE *fr = OpenFile(in_file, "r");
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(in_file, "");
    exit(1);
  } //}}}
  // read project name //{{{
  line_count++;
  ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, "\n");
  char project_id[LINE] = "";
  if (snprintf(project_id, LINE, "%s", split[0]) < 0) {
    ErrorSnprintf();
  } //}}}
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(in_file, "");
    exit(1);
  } //}}}
  // read number of species //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(in_file, "");
    exit(1);
  }
  long int val;
  if (!IsWholeNumber(split[0], &val)) {
    strcpy(ERROR_MSG, "wrong number of species line");
    PrintErrorFileLine(in_file, line_count, split, words);
    exit(1);
  } //}}}
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  Count->MoleculeType = val;
  Count->Molecule = val;
  System.MoleculeType = realloc(System.MoleculeType,
                                 sizeof *System.MoleculeType *
                                 Count->MoleculeType);
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(in_file, "");
    exit(1);
  } //}}}
  // read all species //{{{
  count = 0; // count molecules
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    InitMoleculeType(mt);
    // read species line: <name> <number>
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(in_file, "");
      exit(1);
    }
    long int val;
    if (words < 2 || !IsWholeNumber(split[1], &val)) {
      strcpy(ERROR_MSG, "wrong species line");
      PrintErrorFileLine(in_file, line_count, split, words);
      exit(1);
    }
    mt->Number = val;
    snprintf(mt->Name, MOL_NAME, "%s", split[0]);
  }
  fclose(fr); //}}}
  //}}}

  // read library
  // 1) bond types //{{{
  char file[LINE] = "";
  if (snprintf(file, LINE, "%s/list_bonds.txt", lib_dir) < 0) {
    ErrorSnprintf();
  }
  fr = OpenFile(file, "r");
  line_count = 0;
  // skip 2 comments //{{{
  for (int i = 0; i < 2; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
  } //}}}
  // read number of bonds //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  }
  if (!IsWholeNumber(split[0], &val)) {
    strcpy(ERROR_MSG, "wrong number of bonds line");
    PrintErrorFileLine(file, line_count, split, words);
    exit(1);
  } //}}}
  Count->BondType = val;
  System.BondType = realloc(System.BondType,
                             sizeof *System.BondType * Count->BondType);
  char (*bond_name)[BEAD_NAME] = malloc(sizeof *bond_name *
                                        Count->BondType);
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  } //}}}
  // read the bond types
  for (int i = 0; i < Count->BondType; i++) {
    System.BondType[i] = InitParams;
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
    if (words < 3 || !IsRealNumber(split[1], &System.BondType[i].a) ||
                     !IsRealNumber(split[2], &System.BondType[i].b)) {
      strcpy(ERROR_MSG, "wrong bond line");
      PrintErrorFileLine(file, line_count, split, words);
      exit(1);
    }
    if (snprintf(bond_name[i], BEAD_NAME, "%s", split[0]) < 0) {
      ErrorSnprintf();
    }
    System.BondType[i].a *= 2; // assumes lmp-style spring constant, i.e., k/2
  }
  fclose(fr); //}}}
  // 2) angle types //{{{
  if (snprintf(file, LINE, "%s/list_angles.txt", lib_dir) < 0) {
    ErrorSnprintf();
  }
  fr = OpenFile(file, "r");
  line_count = 0;
  // skip 2 comments //{{{
  for (int i = 0; i < 2; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
  } //}}}
  // read number of angles //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  }
  if (!IsWholeNumber(split[0], &val)) {
    strcpy(ERROR_MSG, "wrong number of angles line");
    PrintErrorFileLine(file, line_count, split, words);
    exit(1);
  } //}}}
  Count->AngleType = val;
  System.AngleType = realloc(System.AngleType,
                              sizeof *System.AngleType * Count->AngleType);
  char (*angle_name)[BEAD_NAME] = malloc(sizeof *angle_name *
                                         Count->AngleType);
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  } //}}}
  // read the angle types
  for (int i = 0; i < Count->AngleType; i++) {
    System.AngleType[i] = InitParams;
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
    if (words < 3 || !IsRealNumber(split[1], &System.AngleType[i].a) ||
                     !IsRealNumber(split[2], &System.AngleType[i].b)) {
      strcpy(ERROR_MSG, "wrong angle line");
      PrintErrorFileLine(file, line_count, split, words);
      exit(1);
    }
    if (snprintf(angle_name[i], BEAD_NAME, "%s", split[0]) < 0) {
      ErrorSnprintf();
    }
    System.AngleType[i].a *= 2; // assumes lmp-style spring constant, i.e., k/2
  }
  fclose(fr); //}}}
  // 3) parameters for bead types //{{{
  if (snprintf(file, LINE, "%s/list_parameters.txt", lib_dir) < 0) {
    ErrorSnprintf();
  }
  fr = OpenFile(file, "r");
  line_count = 0;
  // skip 2 comments //{{{
  for (int i = 0; i < 2; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
  } //}}}
  // read number of beads //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  }
  if (!IsWholeNumber(split[0], &val)) {
    strcpy(ERROR_MSG, "wrong number of beads line");
    PrintErrorFileLine(file, line_count, split, words);
    exit(1);
  } //}}}
  Count->BeadType = val;
  System.BeadType = realloc(System.BeadType,
                             sizeof *System.BeadType * Count->BeadType);
  double *a_ii = calloc(Count->BeadType, sizeof *a_ii);
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  } //}}}
  // read the parameters
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt = &System.BeadType[i];
    InitBeadType(bt);
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
    if (words < 6 || !IsRealNumber(split[2], &bt->Mass) || // mass
                     !IsRealNumber(split[3], &bt->Charge) || // charge
                     !IsRealNumber(split[4], &a_ii[i]) || // radius
                     !IsRealNumber(split[5], &bt->Radius)) { // a_{ii}
      strcpy(ERROR_MSG, "wrong bead parameter line");
      PrintErrorFileLine(file, line_count, split, words);
      exit(1);
    }
    if (snprintf(bt->Name, BEAD_NAME, "%s", split[1]) < 0) {
      ErrorSnprintf();
    }
  }
  fclose(fr); //}}}
  // 4) cross parameters for bead types //{{{
  if (snprintf(file, LINE, "%s/list_cross_interactions.txt", lib_dir) < 0) {
    ErrorSnprintf();
  }
  fr = OpenFile(file, "r");
  line_count = 0;
  // skip 2 comments //{{{
  for (int i = 0; i < 2; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
  } //}}}
  // read number of bonds //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  }
  long int n_cross;
  if (!IsWholeNumber(split[0], &n_cross)) {
    strcpy(ERROR_MSG, "wrong number of beads line");
    PrintErrorFileLine(file, line_count, split, words);
    exit(1);
  } //}}}
  PARAMS *cross = calloc(n_cross, sizeof *cross);
  // skip comment //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    ErrorEOF(file, "");
    exit(1);
  } //}}}
  // read the parameters
  for (int i = 0; i < n_cross; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
    if (words < 5 || !IsRealNumber(split[3], &cross[i].c) || // a_{ij}
                     !IsRealNumber(split[4], &cross[i].d)) { // R_{ij}
      strcpy(ERROR_MSG, "wrong cross parameter line");
      PrintErrorFileLine(file, line_count, split, words);
      exit(1);
    }
    cross[i].a = -1;
    cross[i].b = -1;
    // identify the two bead names
    for (int j = 0; j < Count->BeadType; j++) {
      BEADTYPE *bt = &System.BeadType[j];
      if (cross[i].a == -1 && strcmp(bt->Name, split[1]) == 0) {
        cross[i].a = j;
      }
      if (cross[i].b == -1 && strcmp(bt->Name, split[2]) == 0) {
        cross[i].b = j;
      }
      if (cross[i].a != -1 && cross[i].b != -1) {
        break;
      }
    }
  }
  fclose(fr); //}}}
  // 5) using in_file's species, read their individual files //{{{
  System.Molecule = realloc(System.Molecule,
                             sizeof *System.Molecule * Count->Molecule);
  int count_mol_beads = 0;
  Count->Bead = 0;
  BEAD *mol_beads = malloc(sizeof *mol_beads);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (snprintf(file, LINE, "%s/%s.txt", lib_dir, mt->Name) < 0) {
      ErrorSnprintf();
    }
    fr = OpenFile(file, "r");
    line_count = 0;
    // skip 2 comments //{{{
    for (int j = 0; j < 2; j++) {
      line_count++;
      if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
        ErrorEOF(file, "");
        exit(1);
      }
    } //}}}
    // read species type: 0 - solvent, 1 - counterion, 2 - molecule //{{{
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
    long int type;
    if (!IsWholeNumber(split[0], &type) || type > 2) {
      strcpy(ERROR_MSG, "wrong type of species line");
      PrintErrorFileLine(file, line_count, split, words);
      exit(1);
    } //}}}
    // skip comment //{{{
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    } //}}}
    // read number of beads //{{{
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    }
    if (!IsNaturalNumber(split[0], &val)) {
      strcpy(ERROR_MSG, "wrong type of species line");
      PrintErrorFileLine(file, line_count, split, words);
      exit(1);
    } //}}}
    mt->nBeads = val;
    mt->Bead = calloc(mt->nBeads, sizeof *mt->Bead);
    Count->Bead += mt->Number * mt->nBeads;
    count_mol_beads += mt->nBeads;
    mol_beads = realloc(mol_beads, sizeof *mol_beads * count_mol_beads);
    // skip comment //{{{
    line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "");
      exit(1);
    } //}}}
    // read beads //{{{
    for (int j = 0; j < mt->nBeads; j++) {
      line_count++;
      if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
        ErrorEOF(file, "");
        exit(1);
      }
      int id = count_mol_beads - mt->nBeads + j;
      BEAD *b = &mol_beads[id];
      if (words < 5 || !IsRealNumber(split[2], &b->Position.x) ||
                       !IsRealNumber(split[3], &b->Position.y) ||
                       !IsRealNumber(split[4], &b->Position.z)) {
        strcpy(ERROR_MSG, "wrong bead line");
        PrintErrorFileLine(file, line_count, split, words);
        exit(1);
      }
      int k;
      for (k = 0; k < Count->BeadType; k++) {
        if (strcmp(System.BeadType[k].Name, split[1]) == 0) {
          break;
        }
      }
      mt->Bead[j] = k;
      System.BeadType[k].Number += mt->Number;
    } //}}}
    if (type == 2) { // continue reading if not monomer bead
      // skip comment //{{{
      line_count++;
      if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
        ErrorEOF(file, "");
        exit(1);
      } //}}}
      // read keyword 'Bonds' or 'Angles' //{{{
      line_count++;
      if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
        ErrorEOF(file, "");
        exit(1);
      } //}}}
      while (strcasecmp(split[0], "End") != 0) {
        if (strcasecmp(split[0], "Bonds") == 0) {
          // skip comment //{{{
          line_count++;
          if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
            ErrorEOF(file, "");
            exit(1);
          } //}}}
          // read number of bonds //{{{
          line_count++;
          if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
            ErrorEOF(file, "");
            exit(1);
          }
          if (!IsNaturalNumber(split[0], &val)) {
            strcpy(ERROR_MSG, "wrong number of bonds line");
            PrintErrorFileLine(file, line_count, split, words);
            exit(1);
          } //}}}
          mt->nBonds = val;
          mt->Bond = calloc(mt->nBonds, sizeof *mt->Bond);
          // skip comment //{{{
          line_count++;
          if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
            ErrorEOF(file, "");
            exit(1);
          } //}}}
          // read bonds //{{{
          for (int j = 0; j < mt->nBonds; j++) {
            // bond line: <bond name> <bead 1> <bead 2>
            line_count++;
            if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
              ErrorEOF(file, "");
              exit(1);
            }
            long int bond[2];
            if (words < 3 ||
                !IsNaturalNumber(split[1], &bond[0]) || bond[0] > mt->nBeads ||
                !IsNaturalNumber(split[2], &bond[1]) || bond[1] > mt->nBeads) {
              strcpy(ERROR_MSG, "wrong bond line");
              PrintErrorFileLine(file, line_count, split, words);
              exit(1);
            }
            mt->Bond[j][0] = bond[0] - 1;
            mt->Bond[j][1] = bond[1] - 1;
            mt->Bond[j][2] = -1;
            // identify bond type from the name
            for (int k = 0; k < Count->BondType; k++) {
              if (strcmp(split[0], bond_name[k]) == 0) {
                mt->Bond[j][2] = k;
              }
            }
          } //}}}
        } else if (strcasecmp(split[0], "Angles") == 0) {
          // skip comment //{{{
          line_count++;
          if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
            ErrorEOF(file, "");
            exit(1);
          } //}}}
          // read number of angles //{{{
          line_count++;
          if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
            ErrorEOF(file, "");
            exit(1);
          }
          if (!IsNaturalNumber(split[0], &val)) {
            strcpy(ERROR_MSG, "wrong number of angles line");
            PrintErrorFileLine(file, line_count, split, words);
            exit(1);
          } //}}}
          // skip comment //{{{
          line_count++;
          if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
            ErrorEOF(file, "");
            exit(1);
          } //}}}
          mt->nAngles = val;
          mt->Angle = calloc(mt->nAngles, sizeof *mt->Angle);
          // read angles //{{{
          for (int j = 0; j < mt->nAngles; j++) {
            // angle line: <bond name> <bead 1> <bead 2> <bead 3>
            line_count++;
            if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
              ErrorEOF(file, "");
              exit(1);
            }
            long int ang[3];
            if (words < 4 ||
                !IsNaturalNumber(split[1], &ang[0]) || ang[0] > mt->nBeads ||
                !IsNaturalNumber(split[2], &ang[1]) || ang[1] > mt->nBeads ||
                !IsNaturalNumber(split[3], &ang[2]) || ang[2] > mt->nBeads) {
              strcpy(ERROR_MSG, "wrong angle line");
              PrintErrorFileLine(file, line_count, split, words);
              exit(1);
            }
            mt->Angle[j][0] = ang[0] - 1;
            mt->Angle[j][1] = ang[1] - 1;
            mt->Angle[j][2] = ang[2] - 1;
            mt->Angle[j][3] = -1;
            // identify bond type from the name
            for (int k = 0; k < Count->AngleType; k++) {
              if (strcmp(split[0], angle_name[k]) == 0) {
                mt->Angle[j][3] = k;
              }
            }
          } //}}}
        }
        // skip comment //{{{
        line_count++;
        if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
          ErrorEOF(file, "");
          exit(1);
        } //}}}
        // read keyword 'Bonds' or 'Angles' //{{{
        line_count++;
        if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
          ErrorEOF(file, "");
          exit(1);
        } //}}}
      }
    }
    fclose(fr);
  } //}}}

  // fill in the rest of the System //{{{
  // fill Bead & Molecule structs
  int count_bead = 0, id = 0;
  Count->Molecule = 0;
  System.Bead = realloc(System.Bead, sizeof *System.Bead * Count->Bead);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (mt->nBeads != 1) {
      System.Molecule = realloc(System.Molecule, sizeof *System.Molecule *
                                 (Count->Molecule + mt->Number));
      Count->Bonded += mt->nBeads * mt->Number;
    } else {
      Count->Unbonded += mt->Number;
    }
    for (int j = 0; j < mt->Number; j++) {
      if (mt->nBeads == 1) {
        BEAD *b = &System.Bead[count_bead];
        b->InTimestep = true;
        b->Position = mol_beads[id].Position;
        b->Molecule = -1;
        b->Type = mt->Bead[0];
        count_bead++;
      } else {
        MOLECULE *mol = &System.Molecule[Count->Molecule];
        mol->Bead = calloc(mt->nBeads, sizeof *mol->Bead);
        for (int k = 0; k < mt->nBeads; k++) {
          BEAD *b = &System.Bead[count_bead];
          b->InTimestep = true;
          b->Position = mol_beads[id+k].Position;
          b->Molecule = Count->Molecule;
          b->Type = mt->Bead[k];
          mol->Bead[k] = count_bead;
          mol->Type = i;
          mol->Index = Count->Molecule;
          count_bead++;
        }
        Count->Molecule++;
      }
    }
    if (mt->nBeads == 1) {
      mt->Number = 0;
    }
    id += mt->nBeads;
  }
  Count->HighestResid = Count->Molecule;
  System.BeadCoor = realloc(System.BeadCoor,
                             Count->Bead * sizeof *System.BeadCoor);
  FillSystemNonessentials(&System); //}}}

  SYSTEM Library = CopySystem(System);
  PruneSystem(&System);
  MergeMoleculeTypes(&System);
  MergeBeadTypes(&System, false);
  CheckSystem(System, file);

  if (verbose) {
    VerboseOutput(System);
  }

  // write output FIELD file
  WriteStructure(add_type, add_field, System, -1, false);

  SYSTEM S_data = ReadStructure(in_data_type, in_data, 0, "\0",
                                false, false, false);

  // opent lammpst input script
  fr = OpenFile(in_lmp, "r");

  bool working = false,
       *named = calloc(S_data.Count.BeadType, sizeof *named);
  do {
    // find next 'pair_coeff' line
    do {
      if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
        ErrorEOF(in_file, "");
        exit(1);
      }
    } while (words == 0 || strcmp(split[0], "pair_coeff") != 0);
    // if the line is a valid line, read the bead name
    long int id[2];
    if (words >= 11 && strcmp(split[3], "dpd") == 0 &&
        IsNaturalNumber(split[1], &id[0]) &&
        IsNaturalNumber(split[2], &id[1])) {
      id[0]--; // lammps bead type ids start from 1
      id[1]--; //
      if (!named[id[0]]) {
        snprintf(S_data.BeadType[id[0]].Name, BEAD_NAME, "%s", split[8]);
        named[id[0]] = true;
      }
      if (!named[id[1]]) {
        snprintf(S_data.BeadType[id[1]].Name, BEAD_NAME, "%s", split[10]);
        named[id[1]] = true;
      }
    }
    working = false;
    for (int i = 0; i < S_data.Count.BeadType; i++) {
      if (!named[i]) {
        working = true;
      }
    }
  } while (working);
  fclose(fr);

  // write the system from data file into FIELD with proper bead names
  WriteStructure(orig_type, orig_field, S_data, 0, false);

  // free memory //{{{
  free(bond_name);
  free(angle_name);
  FreeSystem(&System);
  FreeSystem(&S_data);
  FreeSystem(&Library);
  free(a_ii);
  free(cross);
  free(mol_beads);
  free(named); //}}}

  return 0;
}
