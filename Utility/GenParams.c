#include "../src/AnalysisTools.h"
#include "../src/ReadLibrary.h"
#include <dirent.h>

bool file_exists(char *filename) {
  struct stat buffer;
  return (stat(filename, &buffer) == 0);
}

static bool dir_has_field_files(const char *dir) {
  DIR *d = opendir(dir);
  if (!d) return false;
  struct dirent *e;
  while ((e = readdir(d)) != NULL) {
    int len = strlen(e->d_name);
    if (len > 6 && strcasecmp(e->d_name + len - 6, ".FIELD") == 0) {
      closedir(d);
      return true;
    }
  }
  closedir(d);
  return false;
}

static void db_validate(const char *db_dir) {
  char path[LINE];
  snprintf(path, LINE, "%slist_molecules.txt", db_dir);
  if (!file_exists(path) && !dir_has_field_files(db_dir)) {
    if (snprintf(ERROR_MSG, LINE,
                 "directory '%s' contains neither list_molecules.txt nor *.FIELD files",
                 db_dir) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  }
}

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "GenParams assembles force field parameters for a molecular system and writes "
  "them in one or more output formats. System composition is read from a text "
  "file (<input>); molecule definitions are loaded from a directory given by "
  "-db (default: ./). The directory type is detected automatically: if it "
  "contains list_molecules.txt it is treated as a molecule library (interactions "
  "included); otherwise individual <name>.FIELD files are expected. -ldata reads "
  "bead/bond/angle types from a LAMMPS data file instead; -db can still be used "
  "alongside it for interaction parameters. -lmp writes LAMMPS "
  "pair/bond/angle_coeff lines. Both <input> and <output> are optional depending "
  "on the mode used.",

  "Usage: GenParams [<input>] [<output>] [options]",
  .args = 0, // both <input> and <output> are optional
  .all = 9, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"[<input>]", NULL, "input text file (optional with -ldata/-db)", OPT_ARG},
  {"[<output>]", NULL, "output structure file (optional with -ldata/-db)", OPT_ARG},
  {"-lmp", "<file>", "write LAMMPS params.in (pair/bond/angle_coeff lines)", OPT_EXTRA},
  {"-ldata", "<file>", "read bead/bond/angle types from a LAMMPS data file", OPT_EXTRA},
  {"-db", "<dir>", "molecule database directory; auto-detected as library (list_molecules.txt) or FIELD files (default: ./)", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, false, HelpDesc, opts);
  int count = 0;
  // [<input>] - optional input text file
  char in[LINE] = "";
  if (count + 1 < argc && argv[count+1][0] != '-') {
    s_strcpy(in, argv[++count], LINE);
  }
  // [<output>] - optional output structure file
  FILE_TYPE out;
  out.name[0] = '\0';
  out.type = -1;
  if (count + 1 < argc && argv[count+1][0] != '-') {
    s_strcpy(out.name, argv[++count], LINE);
    out.type = StructureFileType(out.name);
  }

  SYS_FILES trash = InitSysFiles; // unused
  COMMON_OPT commons = CommonOptions(argc, argv, trash);
  char lmp_file[LINE] = "";
  FileOption(argc, argv, "-lmp", lmp_file);
  char ldata_file[LINE] = "";
  FileOption(argc, argv, "-ldata", ldata_file);
  char db_dir[LINE] = "./";
  bool db_given = FileOption(argc, argv, "-db", db_dir);
  // ensure db_dir ends with /
  int db_len = strlen(db_dir);
  if (db_len > 0 && db_dir[db_len-1] != '/') {
    strncat(db_dir, "/", LINE - db_len - 1);
  }
  // auto-detect directory type; validate early if -db was explicit
  char list_mol_path[LINE];
  snprintf(list_mol_path, LINE, "%slist_molecules.txt", db_dir);
  bool use_lib = file_exists(list_mol_path);
  if (db_given) {
    db_validate(db_dir);
  } //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;

  LIBRARY lib = {0};
  bool lib_molecules = false; // true when -lib used without -ldata for molecule loading
  FILE *f = NULL;

  // read the input text file (first pass: system composition + default potential) //{{{
  int line_count = 0;
  double density = 0,
         dpd[3] = {25, 1, 4.5};
  if (in[0] != '\0') {
    f = OpenFile(in, "r");
    while (ReadAndSplitLine(f, 7, " \t\n")) {
      line_count++;
      if (words == 0 || split[0][0] == '#') { // blank or comment line
        continue;
      } else if (strncasecmp("potential", split[0], 3) == 0) { //{{{
        double val[3] = {-1, -1, -1};
        if (words < 4 || strncasecmp("dpd", split[1], 3) != 0 ||
            (words > 4 && !IsPosRealNumber(split[4], &val[0])) ||
            (words > 5 && !IsPosRealNumber(split[5], &val[1])) ||
            (words > 6 && !IsPosRealNumber(split[6], &val[2]))) {
          fclose(f);
          goto err_in;
        } else if (strcmp(split[2], "*") == 0 && strcmp(split[3], "*") == 0) {
          if (val[0] != -1) { dpd[0] = val[0]; }
          if (val[1] != -1) { dpd[1] = val[1]; }
          if (val[2] != -1) { dpd[2] = val[2]; }
        } //}}}
      } else if (strncasecmp("box", split[0], 3) == 0) { //{{{
        if (ldata_file[0] == '\0' && !use_lib) {
          if (words >= 4) {
            if (!IsPosRealNumber(split[1], &System.Box.Length.x) ||
                !IsPosRealNumber(split[2], &System.Box.Length.y) ||
                !IsPosRealNumber(split[3], &System.Box.Length.z)) {
              fclose(f);
              goto err_in;
            }
          } else if (words >= 2) {
            if (!IsPosRealNumber(split[1], &System.Box.Volume)) {
              fclose(f);
              goto err_in;
            }
            System.Box.Length.x = cbrt(System.Box.Volume);
            System.Box.Length.y = System.Box.Length.x;
            System.Box.Length.z = System.Box.Length.x;
          } else {
            fclose(f);
            goto err_in;
          }
          CalculateBoxData(&System.Box, 0);
        } //}}}
      } else if (strncasecmp("density", split[0], 3) == 0) { //{{{
        if (ldata_file[0] == '\0' && !use_lib) {
          if (words < 2 || !IsPosRealNumber(split[1], &density)) {
            fclose(f);
            goto err_in;
          }
        } //}}}
      } else if (strncasecmp("molecule", split[0], 3) == 0) { //{{{
        if (ldata_file[0] == '\0') {
          long int n_Sys;
          if (words < 3 || (strncasecmp("fill", split[2], 1) != 0 &&
                            !IsWholeNumber(split[2], &n_Sys))) {
            fclose(f);
            goto err_in;
          }
          if (!db_given) { // default ./: validate lazily on first molecule
            db_validate(db_dir);
            db_given = true; // only validate once
          }
          if (use_lib) {
            // library mode: load molecule from db_dir
            lib_molecules = true;
            if (!IsWholeNumber(split[2], &n_Sys) || n_Sys <= 0) { continue; }
            // save mol name before ReadLibrary overwrites global split[]
            char mol_name[MOL_NAME];
            s_strcpy(mol_name, split[1], MOL_NAME);
            if (lib.n_inter == 0) { // initialise library on first molecule
              lib = ReadLibrary(db_dir);
            }
            ReadLibraryMolecule(db_dir, mol_name, (int)n_Sys, &lib);
          } else {
            // -db: load molecule from FIELD file (original behaviour)
            if (strncasecmp("fill", split[2], 1) == 0) {
              n_Sys = System.Box.Volume * density - Count->Bead;
              if (n_Sys <= 0) {
                if (snprintf(ERROR_MSG, LINE,
                             "No beads to 'fill': %ld beads too many", -n_Sys) < 0) {
                  ErrorSnprintf();
                }
                PrintWarning();
              }
            } else if (n_Sys == 0) {
              continue;
            }
            char name[MOL_NAME];
            snprintf(name, MOL_NAME, "%s", split[1]);
            SYS_FILES fmol = InitSysFiles;
            if (snprintf(fmol.stru.name, LINE, "%s%s.FIELD", db_dir, name) < 0) {
              ErrorSnprintf();
            }
            fmol.stru.type = FIELD_FILE;
            if (file_exists(fmol.stru.name)) {
              SYSTEM Sys = ReadStructure(fmol, false);
              Sys.Count.BeadCoor = Sys.Count.Bead;
              for (int i = 0; i < Sys.Count.Bead; i++) {
                Sys.BeadCoor[i] = i;
              }
              FillInCoor(&Sys);
              for (int i = 0; i < n_Sys; i++) {
                ConcatenateSystems(&System, Sys, System.Box, false);
              }
              FreeSystem(&Sys);
            } else {
              if (snprintf(ERROR_MSG, LINE, "Missing file %s", fmol.stru.name) < 0) {
                ErrorSnprintf();
              }
              fclose(f);
              PrintError();
              exit(1);
            }
          }
        } //}}}
      } else {
        fclose(f);
        goto err_in;
      }
    }
    fclose(f);
  }
  // finalise the system
  if (ldata_file[0] != '\0') {
    // read system topology (bead/bond/angle types + coeffs) from LAMMPS data file
    SYS_FILES ldata_f = InitSysFiles;
    s_strcpy(ldata_f.stru.name, ldata_file, LINE);
    ldata_f.stru.type = LDATA_FILE;
    FreeSystem(&System);
    System = ReadStructure(ldata_f, false);
    System.Count.BeadCoor = System.Count.Bead;
    for (int i = 0; i < System.Count.Bead; i++) {
      System.BeadCoor[i] = i;
    }
    FillInCoor(&System);
    if (use_lib && lib.n_inter == 0) {
      lib = ReadLibrary(db_dir); // load library for interactions only
    }
  } else if (lib_molecules) {
    // -lib without -ldata: finalise system built from library molecules
    FillSystemNonessentials(&lib.System, true); // fix Bonded[], BeadType[].Index, etc.
    System = lib.System;
    Count = &System.Count;
    System.BeadCoor = s_realloc(System.BeadCoor,
                                Count->Bead * sizeof *System.BeadCoor);
    System.Count.BeadCoor = Count->Bead;
    for (int i = 0; i < Count->Bead; i++) {
      System.BeadCoor[i] = i;
    }
    for (int i = 0; i < Count->Molecule; i++) {
      System.Molecule[i].InTimestep = true;
    }
    PruneSystem(&System, NULL);
  } else {
    for (int i = 0; i < Count->Molecule; i++) {
      System.Molecule[i].InTimestep = true;
    }
    PruneSystem(&System, NULL);
  } //}}}

  // array for dpd parameters: a_ij, r_c, gamma //{{{
  ArrNDd *pot = CreateArr3Dd(Count->BeadType, Count->BeadType, 3);
  if (!pot) {
    ErrorAlloc("pot");
  }
  if (use_lib) {
    FillPotFromLibrary(&lib, &System, pot);
  } else {
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = 0; j < Count->BeadType; j++) {
        for (int aa = 0; aa < 3; aa++) {
          SetArr3D(pot, i, j, aa, dpd[aa]);
        }
      }
    }
  } //}}}

  // reread the file to apply specific potential overrides //{{{
  if (in[0] != '\0') {
    FILE *f2 = OpenFile(in, "r");
    while (ReadAndSplitLine(f2, SPL_STR, " \t\n")) {
      line_count++;
      if (words > 0 && strncasecmp("potential", split[0], 3) == 0) {
        double val[3] = {-1, -1, -1};
        if (words < 4 || strncasecmp("dpd", split[1], 3) != 0 ||
            (words > 4 && !IsPosRealNumber(split[4], &val[0])) ||
            (words > 5 && !IsPosRealNumber(split[5], &val[1])) ||
            (words > 6 && !IsPosRealNumber(split[6], &val[2]))) {
          fclose(f2);
          goto err_in;
        }
        int bt_1 = FindBeadType(split[2], System),
            bt_2 = FindBeadType(split[3], System);
        if (bt_1 != -1 && bt_2 != -1) {
          if (bt_1 > bt_2) {
            SwapInt(&bt_1, &bt_2);
          }
          for (int aa = 0; aa < 3; aa++) {
            if (val[aa] != -1) {
              SetArr3D(pot, bt_1, bt_2, aa, val[aa]);
            }
          }
        }
      }
    }
    fclose(f2);
  } //}}}

  if (commons.verbose) { //{{{
    VerboseOutput(System);
    fprintf(stdout, "Potentials:\n");
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = i; j < Count->BeadType; j++) {
        fprintf(stdout, "%10s %10s", System.BeadType[i].Name,
                                     System.BeadType[j].Name);
        for (int aa = 0; aa < 3; aa++) {
          fprintf(stdout, " %lf", GetArr3D(pot, i, j, aa));
        }
        putchar('\n');
      }
    }
  } //}}}

  // write FIELD output if <output> was given and -ldata not used //{{{
  if (out.name[0] != '\0' && ldata_file[0] == '\0') {
    WriteOutputAll(System, out, false, false, argc, argv);
    if (out.type == FIELD_FILE) {
      f = OpenFile(out.name, "a");
      int n = Count->BeadType * (Count->BeadType - 1) / 2 + Count->BeadType;
      fprintf(f, "interactions %d <a_ij> <r_c> <gamma>\n", n);
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = i; j < Count->BeadType; j++) {
          fprintf(f, "%10s %10s dpd", System.BeadType[i].Name,
                                      System.BeadType[j].Name);
          for (int aa = 0; aa < 3; aa++) {
            fprintf(f, " %lf", GetArr3D(pot, i, j, aa));
          }
          putc('\n', f);
        }
      }
      fclose(f);
    }
  } //}}}

  // write LAMMPS params.in if -lmp was specified //{{{
  if (lmp_file[0] != '\0') {
    f = OpenFile(lmp_file, "w");
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = i; j < Count->BeadType; j++) {
        fprintf(f, "pair_coeff %4d %4d dpd %9.5f ${gamma} %7.5f # %s - %s\n",
                i + 1, j + 1,
                GetArr3D(pot, i, j, 0), GetArr3D(pot, i, j, 1),
                System.BeadType[i].Name, System.BeadType[j].Name);
      }
    }
    if (Count->BondType > 0) {
      fputc('\n', f);
      for (int i = 0; i < Count->BondType; i++) {
        fprintf(f, "bond_coeff  %4d %g %g\n",
                i + 1, System.BondType[i].a, System.BondType[i].b);
      }
    }
    if (Count->AngleType > 0) {
      fputc('\n', f);
      for (int i = 0; i < Count->AngleType; i++) {
        fprintf(f, "angle_coeff %4d %g %g\n",
                i + 1, System.AngleType[i].a, System.AngleType[i].b);
      }
    }
    fclose(f);
  } //}}}

  FreeArrND(pot);
  if (!lib_molecules) { // lib.System was copied to System; avoid double-free
    FreeLibrary(&lib);
  }
  FreeSystem(&System);

  return 0;

  err_in:
    err_msg("wrong line");
    PrintErrorFileLine(in, line_count);
    exit(1);
}
