#include "Options.h"

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
int OptionCheck(int argc, char *argv[], int auto_c, int req,
                int common, int all, char opt[all][OPT_LENGTH]) {
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
         (argv[count+1][0] != '-' || strlen(argv[count+1]) == 1)) {
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
      if (argv[i][0] != '-' || // assumes an argument to some option
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
    } else if (strcmp(option[i], "--detailed") == 0) {
      fprintf(ptr, "  --detailed        use name as well as charge, mass, "
              "and radius to identfy bead types (vtf structure files only)\n");
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
void CommonOptions(int argc, char *argv[], int length, bool *verbose,
                   bool *silent, bool *detailed,
                   int *start, int *end, int *skip) {
  // -v option - verbose output
  *verbose = BoolOption(argc, argv, "--verbose");
  // --silent option - silent mode
  SilentOption(argc, argv, verbose, silent);
  // --detailed option - base bead types on name, charge, mass, and radius
  *detailed = BoolOption(argc, argv, "--detailed");
  // starting/ending timestep
  *start = 1;
  int trash; // number of values from IntegerOption(); unused
  IntegerOption1(argc, argv, "-st", start);
  *end = -1;
  IntegerOption1(argc, argv, "-e", end);
  ErrorStartEnd(*start, *end);
  // number of timesteps to skip per one used
  if (IntegerOption(argc, argv, 1, "-sk", &trash, skip)) {
    fprintf(stderr, "%sCommand: %s", ErrRed(), ErrColourReset());
    PrintCommand(stderr, argc, argv);
    exit(1);
  }
  (*skip)++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
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
        strcpy(ERROR_MSG, "missing an argument "
               "(or molecule name beginning with a dash)");
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
        strcpy(ERROR_MSG, "missiont file name "
               "(or the file begins with a dash)");
        PrintErrorOption("-j");
        return true;
      } //}}}
      snprintf(file, LINE, "%s", argv[i+1]);
      int ext = 3;
      char extension[4][EXTENSION];
      strcpy(extension[0], ".vcf");
      strcpy(extension[1], ".vtf");
      strcpy(extension[2], ".xyz");
      strcpy(extension[3], ".lammpstrj");
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
        strcpy(ERROR_MSG, "at least one bead type is required");
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
        strcpy(ERROR_MSG, "at least one bead type is required");
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
      int n = 0; // number of arguments
      // read integers
      int arg = i+1+n;
      while (arg < argc && argv[arg][0] != '-') {
        // Error - non-numeric or missing argument
        long val;
        if (!IsIntegerNumber(argv[arg], &val)) {
          strcpy(ERROR_MSG, "each argument must be non-negative whole number");
          PrintErrorOption(opt);
          exit(1);
        }
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
        strcpy(ERROR_MSG, "missing argument(s)");
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
      strcpy(ERROR_MSG, "single numeric argument required");
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
      strcpy(ERROR_MSG, "two numeric arguments required");
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
      while (arg < argc && argv[arg][0] != '-') {
        // Error - non-numeric or missing argument
        double val;
        if (!IsRealNumber(argv[arg], &val)) {
          strcpy(ERROR_MSG, "each argument must be non-negative whole number");
          PrintErrorOption(opt);
          return true;
        }
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
        strcpy(ERROR_MSG, "missing argument(s)");
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
      strcpy(ERROR_MSG, "single numeric argument required");
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
      strcpy(ERROR_MSG, "two numeric arguments required");
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
      strcpy(ERROR_MSG, "three numeric arguments required");
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
} //}}}

// general option with filename and integer(s) arguments //{{{
bool FileIntegerOption(int argc, char *argv[], int max, char opt[],
                       int *values, int *count, char file[]) {
  int n = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - no output file name
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        strcpy(ERROR_MSG, "missing file name "
               "(or the file name begins with a dash)");
        PrintErrorOption(opt);
        return false;
      }
      snprintf(file, LINE, "%s", argv[i+1]);
      // read integers
      if (max == 0) {
        return false;
      } else {
        while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
          // Error - non-numeric or missing argument
          long val;
          if (!IsIntegerNumber(argv[i+2+n], &val)) {
            strcpy(ERROR_MSG, "each argument must be non-negative whole number");
            PrintErrorOption(opt);
            return true;
          }
          values[n] = val;
          n++;
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
          strcpy(ERROR_MSG, "missing numeric argument(s)");
          PrintErrorOption(opt);
          return true;
        }
      }
    }
  }
  *count = n;
  return false;
}
bool FileOption(int argc, char *argv[], char opt[], char file[]) {
  int trash;
  if (FileIntegerOption(argc, argv, 0, opt, &trash, &trash, file)) {
    return true;
  }
  return false;
} //}}}

#if 0 //{{{
// TODO: remove
// general option with a filename argument //{{{
bool FileIntegerOption(int argc, char *argv[], char opt[], char *name, int length) {
  name[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // wrong argument to the option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        strcpy(ERROR_MSG, "missing file name "
               "(or the file name begins with a dash)");
        PrintErrorOption(opt);
        return true;
      }
      if (snprintf(name, LINE, "%s", argv[i+1]) < 0) {
        ErrorSnprintf();
      }
    }
  }
  return false;
} //}}}
// TODO: redo
// MoleculeTypeOption() //{{{
/**
 * Generic option for molecule type that can take one
 * argument. The option is an argument of this function.
 */
bool MoleculeTypeOption(int argc, char *argv[], char opt[], int *moltype,
                        COUNTS Counts, MOLECULETYPE **MoleculeType) {

  *moltype = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        strcpy(ERROR_MSG, "missing file name "
               "(or the file name begins with a dash)");
        PrintErrorOption(opt);
        return true;
      } //}}}
      *moltype = FindMoleculeType_old(argv[i+1], Counts, *MoleculeType);
      if (*moltype == -1) {
        snprintf(ERROR_MSG, LINE, "non-existent molecule %s%s",
                 ErrYellow(), argv[i+1]);
        PrintErrorOption(opt);
        ErrorMoleculeType_old(Counts, *MoleculeType);
        return true;
      }
    }
  }
  return false;
} //}}}
// MoleculeTypeOption2() //{{{
/**
 * Generic option for molecule types that can take multiple arguments. The
 * option is an argument of this function.
 */
bool MoleculeTypeOption2(int argc, char *argv[], char opt[], int *moltype,
                         COUNTS Counts, MOLECULETYPE **MoleculeType) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // set all moltypes not to be used
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        moltype[j] = 0;
      }
      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        strcpy(ERROR_MSG, "missing molecule name "
               "(or the name begins with a dash)");
        PrintErrorOption(opt);
        ErrorMoleculeType_old(Counts, *MoleculeType);
        return true;
      } //}}}
      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {
        int type = FindMoleculeType_old(argv[i+1+j], Counts, *MoleculeType);
        if (type == -1) { // is argv[i+1+j] in vsf?
          strcpy(ERROR_MSG, "non-existent molecule name %s%s"
                 ErrYellow(), argv[i+1+j]);
          PrintErrorOption(opt);
          ErrorMoleculeType_old(Counts, *MoleculeType);
          return true;
        }
        moltype[type] = 1;
        j++;
      }
    }
  }
  return false;
} //}}}
// MoleculeTypeIntOption() //{{{
/**
 * Generic option for a single molecule type followed by a single integer
 * number. The option is an argument of this function.
 */
bool MoleculeTypeIntOption(int argc, int i, char *argv[], char opt[],
                           int *moltype, int *value, COUNTS Counts,
                           MOLECULETYPE *MoleculeType) {
  *moltype = -1;
  if (strcmp(argv[i], opt) == 0) {
    // Error - missing or wrong arguments //{{{
    if ((i+2) >= argc || argv[i+1][0] == '-' || !IsInteger_old(argv[i+2])) {
      strcpy(ERROR_MSG, "two arguments required (<mol name> <int>)");
      PrintErrorOption(opt);
      return true;
    } //}}}
    *moltype = FindMoleculeType_old(argv[i+1], Counts, MoleculeType);
    if (*moltype == -1) {
      snprintf(ERROR_MSG, LINE, "non-existent molecule %s %s"
               ErrYellow(), argv[i+1]);
      PrintErrorOption(opt);
      ErrorMoleculeType_old(Counts, MoleculeType);
      return true;
    }
    *value = atoi(argv[i+2]);
  }
  return false;
} //}}}
#endif //}}}
