#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Utility description",

  "Usage: Template <input> <output> <double> [options]\n\n",
  .args = 3, // number of mandatory arguments
  .all = 17, // number of valid lines OptSpec (not counting last {NULL})
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
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<output>", NULL, "output file", OPT_ARG},
  {"<double>", NULL, "mandatory double argument", OPT_ARG},
  {"--bool", NULL, "bool option (not common ones)", OPT_EXTRA},
  {"-int1", "<int>", "1-integer option", OPT_EXTRA},
  {"-int2", "2x<int>", "2-integer option", OPT_EXTRA},
  {"-f", "<f> [int(s)]", "filename + integer(s) option", OPT_EXTRA},
  {"-mt", "<name(s)>", "use specified molecule type(s)", OPT_EXTRA},
  {"-bt", "<name(s)>", "use specified bead type(s)", OPT_EXTRA},
  {NULL}
}; //}}}

// Help() //{{{
void Help_old(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
} //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  bool bool_opt;     // --bool
  bool *mt;          // -mt
  bool *bt;          // -bt
  int int1;          // -int1
  int int2[2];       // -int2
  char f_file[LINE]; // -f (filename)
  int f_list[100],   // -f (list of numbers)
      f_num;         // -f (number of those numbers)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

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
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // <double>
  // use other Is*Number() from General.c for other possible numbers
  double mandatory_double = 0;
  if (!IsRealNumber(argv[++count], &mandatory_double)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  opt.bool_opt = BoolOption(argc, argv, "--bool");
  FileNumbersOption(argc, argv, 0, 100, "-int", opt.f_list,
                    &opt.f_num, opt.f_file, 'i');
  opt.int1 = -1;
  OneNumberOption(argc, argv, "-int1", &opt.int1, 'i');
  opt.int2[0] = -1;
  TwoNumbersOption(argc, argv, "-int2", opt.int2, 'i');
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  opt.mt = calloc(System.Count.MoleculeType, sizeof *opt.mt);
  if (!TypeOption(argc, argv, "-mt", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }
  opt.bt = calloc(System.Count.BeadType, sizeof *opt.bt);
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System)) {
    InitBoolArray(opt.bt, Count->MoleculeType, true);
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // print options
  fprintf(stdout, "--bool .. %s\n", opt.bool_opt ? "yes" : " no");
  if (opt.int1 != -1) {
    fprintf(stdout, "-int1 .. %d\n", opt.int1);
  }
  if (opt.int2[0] != -1) {
    fprintf(stdout, "-int2 .. %d %d\n", opt.int2[0], opt.int2[1]);
  }
  if (opt.f_file[0] != '\0') {
    fprintf(stdout, "-f .. %s", opt.f_file);
    for (int i = 0; i < opt.f_num; i++) {
      fprintf(stdout, " %d", opt.f_list[i]);
    }
    putchar('\n');
  }
  putchar('\n');

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
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      // go over all beads in the coordinate file
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i]; // bead index
        int type = System.Bead[id].Type;
        if (opt.bt[type]) { // use only specified bead types
          printf("Bead %d of type %s\n", id, System.BeadType[type].Name);
        }
      }
      // go over all molecules in the coordinate file
      for (int i = 0; i < Count->MoleculeCoor; i++) {
        int id = System.MoleculeCoor[i];
        int type = System.Molecule[id].Type;
        if (opt.mt[type]) { // use only specified molecule types
          printf("Molecule %d of type %s\n",
                 id, System.MoleculeType[type].Name);
        }
      }
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
  free(opt.mt);
  free(opt.bt);
  FreeSystem(&System);
  //}}}

  return 0;
}
