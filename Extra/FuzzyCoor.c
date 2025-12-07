#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "FuzzyCoor adds a random value from "
  "<-max, max> interval to coordinates, introducing noise.",

  "Usage: %s <input> <max> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 12, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", NULL, "input coordinate file"},
  {"<max>", NULL, "mandatory double argument"},
  {"<output>", NULL, "output coordinate file"},
  {"-d", "[3x]<axis>", "apply in <axis> direction(s) (default x y z)"},
  {NULL}
}; //}}}

// Help() //{{{
void Help_old(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
} //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  int dim[4]; // -d: [0]...1D/2D/3D; rest axes (0-2)
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <max>
  double max = 0;
  if (!IsRealNumber(argv[++count], &max)) {
    ErrorNaN("<max>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  // <output> - output coordinate file
  FILE_TYPE fout;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name);
  //}}}

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  char dust[LINE];
  opt.dim[0] = 3; // 3D
  if (FileOption(argc, argv, "-2D", dust)) {
    opt.dim[0] = 2; // 2D
  }
  // -d option //{{{
  // default: all axes
  opt.dim[0] = 3;
  for (int dd = 0; dd < opt.dim[0]; dd++) {
    opt.dim[dd+1] = dd;
  }
  // find if the function is there
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0) {
      i++;
      opt.dim[0] = 0;
      for (int j = 0; j < 3; j++) {
        if ((i + j + 1) > argc) {
          break;
        } else if (argv[i+j][0] == 'x') {
          opt.dim[0]++;
          opt.dim[opt.dim[0]] = 0;
        } else if (argv[i+j][0] == 'y') {
          opt.dim[0]++;
          opt.dim[opt.dim[0]] = 1;
        } else if (argv[i+j][0] == 'z') {
          opt.dim[0]++;
          opt.dim[opt.dim[0]] = 2;
        } else {
          err_msg("must be x, y, and/or z");
          PrintErrorOption("-d");
          exit(1);
        }
      }
    }
  } //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // seed random number generator
  srand(time(0));

  InitOutputCoorFile(fout, System, argc, argv);

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (fout.type == LDATA_FILE && count_used == 1) {
        err_msg("only one timestep can be saved to lammps data file");
        PrintWarnFile(fout.name, "\0", "\0");
        count_coor--;
        break;
      }
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      // go over all beads in the coordinate file
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i]; // bead index
        BEAD *b = &System.Bead[id];
        for (int dd = 0; dd < opt.dim[0]; dd++) {
          // random number <-max, max>
          double n = ((double)rand() / RAND_MAX) * 2.0 * max - max;
          b->Position.v[opt.dim[dd+1]] += n;
        }
      }
      bool *write = malloc(Count->BeadCoor * sizeof *write);
      // first assume all are saved/excluded...
      InitBoolArray(write, Count->BeadCoor, true);
      WriteTimestep(fout, System, count_coor, write, argc, argv);
      free(write);
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!commons.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  //}}}

  return 0;
}
