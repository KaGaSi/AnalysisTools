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

// calculate index in an 1D array that simulates a 3D one /{{{
int id3D(int i1, int i2, int i3, int *size) {
  return (i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
} //}}}

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
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

  // define options //{{{
  int common = 3, all = common + 0, count = 0, req_arg = 2;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  // extra options
  OptionCheck(argc, argv, count, req_arg, common, all, option, true); //}}}

  count = 0; // count arguments
  OPT *opt = opt_create();

  char in[LINE];
  snprintf(in, LINE, "%s", argv[++count]);

  FILE_TYPE out;
  snprintf(out.name, LINE, "%s", argv[++count]);
  out.type = StructureFileType(out.name);

  char db_dir[LINE];
  snprintf(db_dir, LINE, "./");

  PrintCommand(stdout, argc, argv);

  SYS_FILES trash = InitSysFiles;
  opt->c = CommonOptions(argc, argv, LINE, trash);

  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;

  // read the input text file //{{{
  FILE *f = OpenFile(in, "r");
  int line_count = 0;
  double density = 0;
  while (ReadAndSplitLine(f, 4, " \t\n")) {
    line_count++;
    if (words == 0 || split[0][0] == '#' || // blank or comment line
        strncasecmp("potential", split[0], 3) == 0) {
      continue;
    }
    if (strncasecmp("box", split[0], 3) == 0) { //{{{
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
  PruneSystem(&System);
  fclose(f); //}}}

  // array for dpd parameters //{{{
  int arr_pot[3];
  arr_pot[0] = Count->BeadType;
  arr_pot[1] = Count->BeadType;
  arr_pot[2] = 3;
  int arr_pot_size = arr_pot[0] * arr_pot[1] * arr_pot[2];
  // dpd potential: a_ij; r_c; gamma
  double *pot = malloc(arr_pot_size * sizeof *pot);
  // default: a_ij = 25, r_c = 1, gamma = 4.5
  for (int i = 0; i < arr_pot[0]; i++) {
    for (int j = 0; j < arr_pot[1]; j++) {
      int id = id3D(i, j, 0, arr_pot);
      pot[id] = 25;
      id = id3D(i, j, 1, arr_pot);
      pot[id] = 1;
      id = id3D(i, j, 2, arr_pot);
      pot[id] = 4.5;
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
          int id = id3D(bt_1, bt_2, aa, arr_pot);
          if (val[aa] != -1) {
            pot[id] = val[aa];
          }
        }
      }
    }
  }
  fclose(f); //}}}

  if (opt->c.verbose) { //{{{
    VerboseOutput(System);
    fprintf(stdout, "Potentials:\n");
    for (int i = 0; i < arr_pot[0]; i++) {
      for (int j = i; j < arr_pot[1]; j++) {
        fprintf(stdout, "%10s %10s", System.BeadType[i].Name,
                                     System.BeadType[j].Name);
        for (int aa = 0; aa < 3; aa++) {
          int id = id3D(i, j, aa, arr_pot);
          fprintf(stdout, " %lf", pot[id]);
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
          int id = id3D(i, j, aa, arr_pot);
          fprintf(f, " %lf", pot[id]);
        }
        putc('\n', f);
      }
    }
    fclose(f);
  } //}}}

  FreeSystem(&System);
  free(pot);
  free(opt);

  return 0;

  err_in:
    strcpy(ERROR_MSG, "wrong line");
    PrintErrorFileLine(in, line_count);
    exit(1);
}