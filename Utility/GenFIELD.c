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

  // read the input text file
  FILE *fr = OpenFile(in, "r");
  int line_count = 0;
  double density = 0;
  while (ReadAndSplitLine(fr, 4, " \t\n")) {
    line_count++;
    if (words == 0) {
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
                        !IsNaturalNumber(split[2], &n_Sys))) {
        goto err_in;
      }
      if (strncasecmp("fill", split[2], 1) == 0) {
        n_Sys = System.Box.Volume * density - Count->Bead;
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
  fclose(fr);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // write the output file if required
  bool *write = malloc(sizeof *write * Count->Bead);
  InitBoolArray(write, Count->Bead, true);
  WriteOutput(System, write, out, false, false, argc, argv);
  free(write);

  FreeSystem(&System);
  free(opt);

  return 0;

  err_in:
    strcpy(ERROR_MSG, "wrong line");
    PrintErrorFileLine(in, line_count);
    exit(1);
}
