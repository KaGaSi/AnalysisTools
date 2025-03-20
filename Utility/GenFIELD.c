#include "../AnalysisTools.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>   // stat

bool file_exists (char *filename) {
  struct stat   buffer;
  return (stat (filename, &buffer) == 0);
}

void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
GenFIELD generates a FIELD file from the supplied database of molecules and \
input text file specifying the system.\
\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <output> [options]\n\n", cmd);
  fprintf(ptr, "<input>             input text file\n");
  fprintf(ptr, "<output>            output structure file\n");
  fprintf(ptr, "[options]\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 4, all = common, count = 0, req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "--verbose", "--silent", "--help", "--version");

  count = 0; // count arguments
  OPT *opt = opt_create();

  char in[LINE];
  s_strcpy(in, argv[++count], LINE);

  FILE_TYPE out;
  s_strcpy(out.name, argv[++count], LINE);
  out.type = StructureFileType(out.name);

  char db_dir[LINE];
  s_strcpy(db_dir, "./", LINE);

  PrintCommand(stdout, argc, argv);

  SYS_FILES trash = InitSysFiles;
  opt->c = CommonOptions(argc, argv, trash);

  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;

  // read the input text file //{{{
  FILE *f = OpenFile(in, "r");
  int line_count = 0;
  double density = 0,
         dpd[3] = {25, 1, 4.5};
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
        goto err_in;
      } else if (strcmp(split[2], "*") == 0 && strcmp(split[3], "*") == 0) {
        if (val[0] != -1) {
          dpd[0] = val[0];
        }
        if (val[1] != -1) {
          dpd[1] = val[1];
        }
        if (val[2] != -1) {
          dpd[2] = val[2];
        }
      } //}}}
    } else if (strncasecmp("box", split[0], 3) == 0) { //{{{
      if (words >= 4) {
        if (!IsPosRealNumber(split[1], &System.Box.Length[0]) ||
            !IsPosRealNumber(split[2], &System.Box.Length[1]) ||
            !IsPosRealNumber(split[3], &System.Box.Length[2])) {
          goto err_in;
        }
      } else if (words >= 2) {
        if (!IsPosRealNumber(split[1], &System.Box.Volume)) {
          goto err_in;
        }
        System.Box.Length[0] = cbrt(System.Box.Volume);
        System.Box.Length[1] = System.Box.Length[0];
        System.Box.Length[2] = System.Box.Length[0];
      } else {
        goto err_in;
      }
      CalculateBoxData(&System.Box, 0); //}}}
    } else if (strncasecmp("density", split[0], 3) == 0) { //{{{
      if (words < 2 || !IsPosRealNumber(split[1], &density)) {
        goto err_in;
      } //}}}
    } else if (strncasecmp("molecule", split[0], 3) == 0) { //{{{
      long int n_Sys;
      if (words < 3 || (strncasecmp("fill", split[2], 1) != 0 &&
                        !IsWholeNumber(split[2], &n_Sys))) {
        goto err_in;
      }
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
      // for (int i = 0; i < words; i++) {
      //   printf(" %s", split[i]);
      // }
      // putchar('\n');
      SYS_FILES f = InitSysFiles;
      if (snprintf(f.stru.name, LINE, "%s%s.FIELD", db_dir, name) < 0) {
        ErrorSnprintf();
      }
      // printf("%s%s%s\n", Yellow(), f.stru.name, ColourReset());
      f.stru.type = FIELD_FILE;
      if (file_exists(f.stru.name)) {
        SYSTEM Sys = ReadStructure(f, false);
        // CheckSystem(Sys, f.stru.name);
        Sys.Count.BeadCoor = Sys.Count.Bead;
        for (int i = 0; i < Sys.Count.Bead; i++) {
          Sys.BeadCoor[i] = i;
        }
        FillInCoor(&Sys);
        for (int i = 0; i < n_Sys; i++) {
          // printf("%d\n", i);
          ConcatenateSystems(&System, Sys, System.Box, false);
        }
        // printf("%sADDED %s%s\n", Magenta(), name, ColourReset());
        FreeSystem(&Sys);
      } else {
        if (snprintf(ERROR_MSG, LINE, "Missing file %s", f.stru.name) < 0) {
          ErrorSnprintf();
        }
        PrintError();
        exit(1);
      }
      //}}}
    } else {
      goto err_in;
    }
  }
  for (int i = 0; i < Count->Molecule; i++) {
    System.Molecule[i].InTimestep = true;
  }
  PruneSystem(&System);
  fclose(f); //}}}

  // array for dpd parameters //{{{
  // dpd potential: a_ij; r_c; gamma
  double (**pot)[3] = calloc(Count->BeadType, sizeof (*pot)[3]);
  // default: a_ij = 25, r_c = 1, gamma = 4.5
  for (int i = 0; i < Count->BeadType; i++) {
    pot[i] = calloc (Count->BeadType, sizeof **pot);
    for (int j = 0; j < Count->BeadType; j++) {
      for (int aa = 0; aa < 3; aa++) {
        pot[i][j][aa] = dpd[aa];
      }
    }
  } //}}}

  // reread the file to get potential parameters //{{{
  f = OpenFile(in, "r");
  while (ReadAndSplitLine(f, SPL_STR, " \t\n")) {
    line_count++;
    if (words > 0 && strncasecmp("potential", split[0], 3) == 0) {
      double val[3] = {-1, -1, -1};
      if (words < 4 || strncasecmp("dpd", split[1], 3) != 0 ||
          (words > 4 && !IsPosRealNumber(split[4], &val[0])) ||
          (words > 5 && !IsPosRealNumber(split[5], &val[1])) ||
          (words > 6 && !IsPosRealNumber(split[6], &val[2]))) {
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
            pot[bt_1][bt_2][aa] = val[aa];
          }
        }
      }
    }
  }
  fclose(f); //}}}

  if (opt->c.verbose) { //{{{
    VerboseOutput(System);
    fprintf(stdout, "Potentials:\n");
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = i; j < Count->BeadType; j++) {
        fprintf(stdout, "%10s %10s", System.BeadType[i].Name,
                                     System.BeadType[j].Name);
        for (int aa = 0; aa < 3; aa++) {
          fprintf(stdout, " %lf", pot[i][j][aa]);
        }
        putchar('\n');
      }
    }
  } //}}}

  // write the output file if required
  WriteOutputAll(System, out, false, false, argc, argv);

  // print parameters to FIELD output file //{{{
  if (out.type == FIELD_FILE) {
    f = OpenFile(out.name, "a");
    int n = Count->BeadType * (Count->BeadType - 1) / 2 + Count->BeadType;
    fprintf(f, "interactions %d <a_ij> <r_c> <gamma>\n", n);
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = i; j < Count->BeadType; j++) {
        fprintf(f, "%10s %10s dpd", System.BeadType[i].Name,
                                    System.BeadType[j].Name);
        for (int aa = 0; aa < 3; aa++) {
          fprintf(f, " %lf", pot[i][j][aa]);
        }
        putc('\n', f);
      }
    }
    fclose(f);
  } //}}}

  for (int i = 0; i < Count->BeadType; i++) {
    free(pot[i]);
  }
  free(pot);
  FreeSystem(&System);
  free(opt);

  return 0;

  err_in:
    err_msg("wrong line");
    PrintErrorFileLine(in, line_count);
    exit(1);
}
