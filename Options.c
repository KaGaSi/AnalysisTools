#include "Options.h"
#include <stdarg.h>

// STATIC DECLARATIONs
static void SilentOption(int argc, char *argv[], bool *verbose, bool *silent);
static bool VersionOption(int argc, char *argv[]);

// STATIC IMPLEMENTATIONS
// option for output verbosity (--silent) //{{{
static void SilentOption(int argc, char *argv[], bool *verbose, bool *silent) {
  *silent = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--silent") == 0) {
      *verbose = false;
      *silent = true;
      break;
    }
  }
} //}}}
// print AnalysisTools version number (--version) //{{{
static bool VersionOption(int argc, char *argv[]) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fprintf(stdout, "AnalysisTools by Karel Šindelka (KaGaSi), version %s"
              " (released %s)\n", VERSION, DATE);
      fprintf(stdout, "Download at https://github.com/KaGaSi/AnalysisTools/"
              "releases/tag/v%s\n", VERSION);
      return true;
    }
  }
  return false;
} //}}}

// THE VISIBLE FUNCTIONS

// print version or help and exit (--version and --help options) //{{{
void HelpVersionOption(int argc, char *argv[]) {
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      // Help(argv[0], false);
      exit(0);
    }
  }
} //}}}

// version/help printing and initial check of provided options //{{{
int OptionCheck2(int argc, char *argv[], int req, int common, int all,
                 bool check_extra, char opt[all][OPT_LENGTH], ...) {
  snprintf(argv[0], LINE, "%s", BareCommand(argv[0]));
  // copy options to an array
  va_list args;
  va_start(args, opt);
  // Loop through the variable arguments, copy each string
  for (int i = 0; i < all; i++) {
    const char *src = va_arg(args, const char *);
    s_strcpy(opt[i], src, OPT_LENGTH);
  }
  va_end(args);
  // --version option?
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  // --help option?
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      Help(argv[0], false, common, opt);
      exit(0);
    }
  }
  // correct number of mandatory options?
  int count = 0;
  while ((count + 1) < argc &&
         // there may be '-' as a mandatory argument
         (argv[count+1][0] != '-' || strnlen(argv[count+1], LINE) == 1)) {
    count++;
  }
  if (count < req) {
    ErrorArgNumber(count, req);
    PrintCommand(stderr, argc, argv);
    Help(argv[0], true, common, opt);
    exit(1);
  }
  // all options exist?
  for (int i = (count+1); i < argc; i++) {
    bool valid = false;
    for (int j = 0; j < all; j++) {
      double value;
      if (argv[i][0] != '-' || // assumes an argument to some option
          IsRealNumber(argv[i], &value) || // assumes negative numeric argument
          strcmp(argv[i], opt[j]) == 0) { // an option
        valid = true;
        break;
      }
    }
    if (!valid) {
      ErrorOption(argv[i]);
      PrintCommand(stderr, argc, argv);
      Help(argv[0], true, common, opt);
      exit(1);
    }
  }
  // warn if extra arguments (between required ones and options)
  if (check_extra && req != count) {
    char extra[LINE] = "\0";
    for (int i = (req + 1); i <= count; i++) {
      char cpy[LINE];
      s_strcpy(cpy, extra, LINE);
      if (snprintf(extra, LINE, "%s %s", cpy, argv[i]) < 0) {
        ErrorSnprintf();
      }
    }
    if (snprintf(ERROR_MSG, LINE, "command line arguments%s%s%s have no effect",
                 ErrYellow(), extra, ErrCyan()) < 0) {
      ErrorSnprintf();
    }
    PrintWarning();
  }
  return count;
} //}}}
// version/help printing and initial check of provided options //{{{
int OptionCheck(int argc, char *argv[], int auto_c, int req, int common,
                int all, char opt[all][OPT_LENGTH], bool check_extra) {
  snprintf(argv[0], LINE, "%s", BareCommand(argv[0]));
  // check the manually specified and automatically counted numbers are the same
  if (auto_c != all) {
    snprintf(ERROR_MSG, LINE, "bug: wrong numbers of arguments in %s%s%s "
             "(counted %s%d%s but %s%d%s specified)",
             ErrYellow(), argv[0], ErrRed(), ErrYellow(), auto_c, ErrRed(),
             ErrYellow(), all, ErrRed());
    PrintError();
    exit(1);
  }
  // --version option?
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  // --help option?
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      Help(argv[0], false, common, opt);
      exit(0);
    }
  }
  // correct number of mandatory options?
  int count = 0;
  while ((count + 1) < argc &&
         (argv[count+1][0] != '-' || strnlen(argv[count+1], LINE) == 1)) {
    count++;
  }
  if (count < req) {
    ErrorArgNumber(count, req);
    PrintCommand(stderr, argc, argv);
    Help(argv[0], true, common, opt);
    exit(1);
  }
  // all options exist?
  for (int i = (count+1); i < argc; i++) {
    bool valid = false;
    for (int j = 0; j < all; j++) {
      double value;
      if (argv[i][0] != '-' || // assumes an argument to some option
          IsRealNumber(argv[i], &value) || // assumes negative numeric argument
          strcmp(argv[i], opt[j]) == 0) {
        valid = true;
        break;
      }
    }
    if (!valid) {
      ErrorOption(argv[i]);
      PrintCommand(stderr, argc, argv);
      Help(argv[0], true, common, opt);
      exit(1);
    }
  }
  // warn if extra arguments (between required ones and options)
  if (check_extra && req != count) {
    char extra[LINE] = "\0";
    for (int i = (req + 1); i <= count; i++) {
      char cpy[LINE];
      s_strcpy(cpy, extra, LINE);
      if (snprintf(extra, LINE, "%s %s", cpy, argv[i]) < 0) {
        ErrorSnprintf();
      }
    }
    if (snprintf(ERROR_MSG, LINE, "command line arguments%s%s%s have no effect",
                 ErrYellow(), extra, ErrCyan()) < 0) {
      ErrorSnprintf();
    }
    PrintWarning();
  }
  return count;
} //}}}

// print help for common options //{{{
void CommonHelp(bool error, int n, char option[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
  }
  for (int i = 0; i < n; i++) {
    if (strcmp(option[i], "-i") == 0) {
      fprintf(ptr, "  -i <name>         input structure file if different "
                   "than the coordinate file\n");
    } else if (strcmp(option[i], "-st") == 0) {
      fprintf(ptr, "  -st <int>         starting timestep for calculation\n");
    } else if (strcmp(option[i], "-e") == 0) {
      fprintf(ptr, "  -e <end>          ending timestep for calculation\n");
    } else if (strcmp(option[i], "-sk") == 0) {
      fprintf(ptr, "  -sk <int>         leave out every 'skip' steps\n");
    } else if (strcmp(option[i], "--variable") == 0) {
      fprintf(ptr, "  --variable        vtf coordinate file with indexed "
              "timesteps with varying number of beads\n");
    } else if (strcmp(option[i], "-ltrj") == 0) {
      fprintf(ptr, "  -ltrj <int>       does lammpstrj ids go from 0 or 1?\n");
    } else if (strcmp(option[i], "--verbose") == 0) {
      fprintf(ptr, "  --verbose         verbose output\n");
    } else if (strcmp(option[i], "--silent") == 0) {
      fprintf(ptr, "  --silent          no output (overrides --verbose)\n");
    } else if (strcmp(option[i], "--help") == 0) {
      fprintf(ptr, "  --help            print this help and exit\n");
    } else if (strcmp(option[i], "--version") == 0) {
      fprintf(ptr, "  --version         print version number and exit\n");
    } else {
      snprintf(ERROR_MSG, LINE, "unknown common option %s%s%s!", ErrYellow(),
               option[i], ErrRed());
      PrintError();
      exit(1);
    }
  }
} //}}}

// detect options common for most utilities //{{{
COMMON_OPT CommonOptions(int argc, char *argv[], int length, SYS_FILES f) {
  COMMON_OPT opt;
  opt.start = -1;
  opt.end = -1;
  opt.skip = 0;
  // -v option - verbose output
  opt.verbose = BoolOption(argc, argv, "--verbose");
  // --silent option - silent mode
  SilentOption(argc, argv, &opt.verbose, &opt.silent);
  // starting/ending timestep
  if (IntegerOption1(argc, argv, "-st", &opt.start)) {
    if (opt.start <= 0) {
      s_strcpy(ERROR_MSG, "positive number required", LINE);
      PrintErrorOption("-st");
      exit(1);
    }
  } else {
    opt.start = 1;
  }
  if (IntegerOption1(argc, argv, "-e", &opt.end) && opt.end <= 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption("-e");
    exit(1);
  }
  ErrorStartEnd(opt.start, opt.end);
  // number of timesteps to skip per one used
  if (IntegerOption1(argc, argv, "-sk", &opt.skip) && opt.skip <= 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption("-sk");
    exit(1);
  }
  opt.skip++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  if (f.coor.type == LDATA_FILE) {
    opt.start = 1;
    opt.skip = 1;
    opt.end = 1;
  }
  return opt;
} //}}}

// exclude specified molecule names (-x <mol name(s)>) //{{{
bool ExcludeOption(int argc, char *argv[], SYSTEM *System) {
  // set all molecules to use
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    System->MoleculeType[i].Flag = true;
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-x") == 0) {
      // wrong argument to -x option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG, "missing an argument "
                 "(or molecule name beginning with a dash)", LINE);
        PrintErrorOption("-x");
        exit(1);
      } //}}}
      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {
        int type = FindMoleculeName(argv[i+1+j], *System);
        if (type == -1) { // is it in vsf?
          snprintf(ERROR_MSG, LINE, "non-existent molecule %s%s",
                   ErrYellow(), argv[i+1+j]);
          PrintErrorOption("-x");
          return true;
        } else {
          // exclude that molecule
          System->MoleculeType[type].Flag = false;
        }
        j++;
      }
    }
  }
  return false;
} //}}}

// join aggregates, saving the coordinates (-j <filename>) //{{{
bool JoinCoorOption(int argc, char *argv[], int *coor_type, char file[]) {
  file[0] = '\0'; // no -j option
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      // wrong argument to -j option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG, "missiont file name "
                 "(or the file begins with a dash)", LINE);
        PrintErrorOption("-j");
        return true;
      } //}}}
      snprintf(file, LINE, "%s", argv[i+1]);
      int ext = 3;
      char extension[4][EXTENSION];
      s_strcpy(extension[0], ".vcf", EXTENSION);
      s_strcpy(extension[1], ".vtf", EXTENSION);
      s_strcpy(extension[2], ".xyz", EXTENSION);
      s_strcpy(extension[3], ".lammpstrj", EXTENSION);
      *coor_type = ErrorExtension(file, ext, extension);
    }
  }
  return false;
} //}}}

// tag bead types with true/false (if missing, set all to opposite) //{{{
bool BeadTypeOption(int argc, char *argv[], char opt[],
                    bool use, bool flag[], SYSTEM System) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int pos = i;
      while (++pos < argc && argv[pos][0] != '-') {
        int btype = FindBeadType(argv[pos], System);
        if (btype == -1) {
          if (snprintf(ERROR_MSG, LINE, "non-existent bead name %s%s",
                       ErrYellow(), argv[pos]) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption(opt);
          ErrorBeadType(argv[pos], System);
          exit(1);
        }
        flag[btype] = use;
      }
      if (pos == (i + 1)) {
        s_strcpy(ERROR_MSG, "at least one bead type is required", LINE);
        PrintErrorOption(opt);
        exit(1);
      }
      return true; // option is present
    }
  }
  return false; // option is not present
} //}}}
// tag molecule types with true/false (if missing, set all to opposite) //{{{
bool MoleculeTypeOption(int argc, char *argv[], char opt[],
                        bool use, bool flag[], SYSTEM System) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int pos = i;
      while (++pos < argc && argv[pos][0] != '-') {
        int btype = FindMoleculeName(argv[pos], System);
        if (btype == -1) {
          if (snprintf(ERROR_MSG, LINE, "non-existent molecule name %s%s",
                       ErrYellow(), argv[pos]) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption(opt);
          ErrorMoleculeType(argv[pos], System);
          exit(1);
        }
        flag[btype] = use;
      }
      if (pos == (i + 1)) {
        s_strcpy(ERROR_MSG, "at least one bead type is required", LINE);
        PrintErrorOption(opt);
        exit(1);
      }
      return true; // option is present
    }
  }
  return false; // option is not present
} // }}}

// general boolean option //{{{
bool BoolOption(int argc, char *argv[], char opt[]) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      return true;
    }
  }
  return false;
} // }}}

// general option with multiple integer arguments (up to 'max') //{{{
bool IntegerOption(int argc, char *argv[], int max,
                   char opt[], int *count, int *values) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0, // number of arguments
          arg = i+1+n; // argument ids in the command
      long val;
      while (arg < argc && IsIntegerNumber(argv[arg], &val)) {
        values[n] = val;
        n++;
        arg = i+1+n;
        // warning - too many numeric arguments
        if (n > max) {
          snprintf(ERROR_MSG, LINE, "too many arguments; only the first %d "
                   "used", max);
          PrintWarnOption(opt);
          *count = max;
          return true;
        }
      }
      if (n == 0) {
        s_strcpy(ERROR_MSG, "missing argument(s)", LINE);
        PrintErrorOption(opt);
        return true;
      }
      *count = n;
      return true;
    }
  }
  return false;
}
bool IntegerOption1(int argc, char *argv[], char opt[], int *value) {
  int count = 0;
  if (IntegerOption(argc, argv, 1, opt, &count, value)) {
    if (count != 1) {
      s_strcpy(ERROR_MSG, "single numeric argument required", LINE);
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
}
bool IntegerOption2(int argc, char *argv[], char opt[], int value[2]) {
  int count = 0;
  if (IntegerOption(argc, argv, 2, opt, &count, value)) {
    if (count != 2) {
      s_strcpy(ERROR_MSG, "two numeric arguments required", LINE);
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
} //}}}

// general option with multiple double arguments (up to 'max') //{{{
bool DoubleOption(int argc, char *argv[], int max,
                  char opt[], int *count, double values[max]) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read integers
      int arg = i+1+n;
      double val;
      while (arg < argc && IsRealNumber(argv[arg], &val)) {
        values[n] = val;
        n++;
        arg = i+1+n;
        // warning - too many numeric arguments
        if (n > max) {
          snprintf(ERROR_MSG, LINE, "too many arguments; only the first %d "
                   "used", max);
          PrintErrorOption(opt);
          *count = n;
          return true;
        }
      }
      if (n == 0) {
        s_strcpy(ERROR_MSG, "missing argument(s)", LINE);
        PrintErrorOption(opt);
        return true;
      }
      *count = n;
      return true;
    }
  }
  return false;
}
bool DoubleOption1(int argc, char *argv[], char opt[], double *value) {
  int count = 0;
  if (DoubleOption(argc, argv, 1, opt, &count, value)) {
    if (count != 1) {
      s_strcpy(ERROR_MSG, "single numeric argument required", LINE);
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
}
bool DoubleOption2(int argc, char *argv[], char opt[], double value[2]) {
  int count = 0;
  if (DoubleOption(argc, argv, 2, opt, &count, value)) {
    if (count != 2) {
      s_strcpy(ERROR_MSG, "two numeric arguments required", LINE);
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
}
bool DoubleOption3(int argc, char *argv[], char opt[], double value[3]) {
  int count = 0;
  if (DoubleOption(argc, argv, 3, opt, &count, value)) {
    if (count != 3) {
      s_strcpy(ERROR_MSG, "three numeric arguments required", LINE);
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
} //}}}

// general option with filename and integer(s) arguments //{{{
bool FileIntegerOption(int argc, char *argv[], int min, int max, char opt[],
                       int *values, int *count, char file[]) {
  int n = 0;
  *count = 0;
  file[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - no output file name
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG, "missing file name "
                 "(or the file name begins with a dash)", LINE);
        PrintErrorOption(opt);
        exit(1);
      }
      snprintf(file, LINE, "%s", argv[i+1]);
      // read integers
      if (max == 0) {
        return true;
      } else {
        while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
          // Error - non-numeric or missing argument
          long val;
          if (!IsIntegerNumber(argv[i+2+n], &val)) {
            s_strcpy(ERROR_MSG, "arguments must be non-negative whole numbers",
                     LINE);
            PrintErrorOption(opt);
            exit(1);
          }
          values[n] = val;
          n++;
          // warning - too many numeric arguments
          if (n > max) {
            snprintf(ERROR_MSG, LINE, "too many arguments; only the first %d "
                     "used", max);
            PrintErrorOption(opt);
            *count = n;
            exit(1);
          }
        }
        if (n < min) {
          s_strcpy(ERROR_MSG, "not enough numeric arguments", LINE);
          PrintErrorOption(opt);
          exit(1);
        }
      }
      *count = n;
      return true; // option present
    }
  }
  return false; // option not present
}
bool FileOption(int argc, char *argv[], char opt[], char file[]) {
  int trash;
  if (FileIntegerOption(argc, argv, 0, 0, opt, &trash, &trash, file)) {
    return true;
  }
  return false;
} //}}}

// general option with filename and double(s) arguments //{{{
bool FileDoubleOption(int argc, char *argv[], int max, char opt[],
                      double *values, int *count, char file[]) {
  int n = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - no output file name
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG, "missing file name "
                 "(or the file name begins with a dash)", LINE);
        PrintErrorOption(opt);
        exit(1);
      }
      snprintf(file, LINE, "%s", argv[i+1]);
      // read integers
      if (max == 0) {
        return true;
      } else {
        while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
          // Error - non-numeric or missing argument
          double val;
          if (!IsRealNumber(argv[i+2+n], &val)) {
            s_strcpy(ERROR_MSG, "arguments must be non-negative whole numbers",
                     LINE);
            PrintErrorOption(opt);
            exit(1);
          }
          values[n] = val;
          n++;
          // warning - too many numeric arguments
          if (n > max) {
            snprintf(ERROR_MSG, LINE, "too many arguments; only the first %d "
                     "used", max);
            PrintErrorOption(opt);
            *count = n;
            exit(1);
          }
        }
        if (n == 0) {
          s_strcpy(ERROR_MSG, "missing numeric argument(s)", LINE);
          PrintErrorOption(opt);
          exit(1);
        }
      }
      *count = n;
      return true; // option present
    }
  }
  return false; // option not present
} //}}}
