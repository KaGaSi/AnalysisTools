#include "Options.h"

// STATIC DECLARATIONs
// option for output verbosity (--silent)
static void SilentOption(int argc, char *argv[], bool *verbose, bool *silent);

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

// THE VISIBLE FUNCTIONS
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
      fprintf(ptr, "      -i <name>      input structure file if different "
                   "than the coordinate file\n");
    } else if (strcmp(option[i], "-st") == 0) {
      fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
    } else if (strcmp(option[i], "-e") == 0) {
      fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
    } else if (strcmp(option[i], "-sk") == 0) {
      fprintf(ptr, "      -sk <int>      leave out every 'skip' steps\n");
    } else if (strcmp(option[i], "--detailed") == 0) {
      fprintf(ptr, "      --detailed     use name as well as charge, mass, "
              "and radius to identfy bead types (vtf structure files only)\n");
    } else if (strcmp(option[i], "--variable") == 0) {
      fprintf(ptr, "      --variable     vtf coordinate file with indexed "
              "timesteps with varying number of beads\n");
    } else if (strcmp(option[i], "-pbc") == 0) {
      fprintf(ptr, "      -pbc <int>     position of pbc in xyz file's comment"
              " line (of the first number)\n");
    } else if (strcmp(option[i], "-v") == 0) {
      fprintf(ptr, "      -v             verbose output\n");
    } else if (strcmp(option[i], "--silent") == 0) {
      fprintf(ptr, "      --silent       no output "
              "(overrides verbose option)\n");
    } else if (strcmp(option[i], "-h") == 0) {
      fprintf(ptr, "      -h             print this help and exit\n");
    } else if (strcmp(option[i], "--version") == 0) {
      fprintf(ptr, "      --version      print version number and exit\n");
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
                   bool *silent, bool *detailed, bool *vtf_var_coor,
                   int *pbc_xyz, int *start, int *end, int *skip) {
  // -v option - verbose output
  *verbose = BoolOption(argc, argv, "-v");
  // --silent option - silent mode
  SilentOption(argc, argv, verbose, silent);
  // --detailed option - base bead types on name, charge, mass, and radius
  *detailed = BoolOption(argc, argv, "--detailed");
  // vtf timesteps with variable number of beads
  *vtf_var_coor = BoolOption(argc, argv, "--variable");
  // starting/ending timestep
  *start = 1;
  if (IntegerOption(argc, argv, "-st", start)) {
    exit(1);
  }
  *end = -1;
  if (IntegerOption(argc, argv, "-e", end)) {
    exit(1);
  }
  ErrorStartEnd(*start, *end);
  // number of timesteps to skip per one used
  if (IntegerOption(argc, argv, "-sk", skip)) {
    exit(1);
  }
  (*skip)++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  // position of the first number of pbc in xyz file
  if (IntegerOption(argc, argv, "-pbc", pbc_xyz)) {
    exit(1);
  }
  if (*pbc_xyz == 0) {
    strcpy(ERROR_MSG, "position must be a positive number");
    PrintErrorOption("-pbc");
  }
} //}}}

// print AnalysisTools version number (--version) //{{{
bool VersionOption(int argc, char *argv[]) {
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

// exclude specified molecule names (-x <mol name(s)>) //{{{
bool ExcludeOption(int argc, char *argv[], SYSTEM *System) {
  // set all molecules to use
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    System->MoleculeType[i].Use = true;
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
          ErrorMoleculeType(*System);
          return true;
        } else {
          // exclude that molecule
          (*System).MoleculeType[type].Use = false;
        }
        j++;
      }
    }
  }
  return false;
} //}}}

// TODO: allow others file types, not just vtf/vcf
// join aggregates, saving the coordinates (-j <filename>) //{{{
bool JoinCoorOption(int argc, char *argv[], char *joined_vcf) {
  joined_vcf[0] = '\0'; // no -j option
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      // wrong argument to -j option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        strcpy(ERROR_MSG, "missiont file name "
               "(or the file begins with a dash)");
        PrintErrorOption("-j");
        return true;
      } //}}}
      snprintf(joined_vcf, LINE, "%s", argv[i+1]);
      // test if <joined.vcf> filename ends with '.vcf' //{{{
      char *dot = strrchr(joined_vcf, '.');
      if (!dot || (strcmp(dot, ".vcf") && strcmp(dot, ".vtf"))) {
        strcpy(ERROR_MSG, "file name must have '.vcf' extension");
        PrintErrorOption("-j");
        return true;
      } //}}}
    }
  }
  return false;
} //}}}

// tag which bead types to use (if not present, set to specified value) //{{{
bool BeadTypeOption(int argc, char *argv[], char *opt,
                    bool use, bool flag[], SYSTEM *System) {
  // specify what bead types to use - either specified by 'opt' or use all
  int types = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      types = i; // positon of the '-bt' argument in command
      // <type names> - names of bead types to save
      while (++types < argc && argv[types][0] != '-') {
        int type = FindBeadType(argv[types], *System);
        if (type == -1) {
          snprintf(ERROR_MSG, LINE, "non-existent bead name %s%s",
                   ErrYellow(), argv[types]);
          PrintErrorOption(opt);
          ErrorBeadType(argv[types], *System);
          return true;
        }
        flag[type] = true;
      }
    }
  }
  if (types == -1) {
    for (int i = 0; i < System->Count.BeadType; i++) {
      flag[i] = use;
    }
  }
  return false;
} // }}}

// tag which molecule types to use (if not present, set to specified value) //{{{
bool MoleculeTypeOption(int argc, char *argv[], char *opt,
                        bool use, bool flag[], SYSTEM System) {
  // specify what molecule types to use - either specified by 'opt' or use all
  int types = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      types = i; // positon of the '-bt' argument in command
      // <type names> - names of bead types to save
      while (++types < argc && argv[types][0] != '-') {
        int type = FindMoleculeName(argv[types], System);
        if (type == -1) {
          snprintf(ERROR_MSG, LINE, "non-existent bead name %s%s",
                   ErrYellow(), argv[types]);
          PrintErrorOption(opt);
          ErrorMoleculeType(System);
          return true;
        }
        flag[type] = true;
      }
    }
  }
  if (types == -1) {
    for (int i = 0; i < System.Count.BeadType; i++) {
      flag[i] = use;
    }
  }
  return false;
} // }}}

// general boolean option //{{{
bool BoolOption(int argc, char *argv[], char *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      return true;
    }
  }
  return false;
} // }}}

// general option with one integer argument //{{{
bool IntegerOption(int argc, char *argv[], char *opt, int *value) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - missing argument
      if ((i+1) >= argc) {
        strcpy(ERROR_MSG, "missing numeric argument");
        PrintErrorOption(opt);
        return true;
      }
      // Error - non-numeric
      long val;
      if (!IsIntegerNumber(argv[i+1], &val)) {
        strcpy(ERROR_MSG, "argument must be non-negative whole number");
        PrintErrorOption(opt);
        return true;
      }
      *value = val;
    }
  }
  return false;
} //}}}

// general option with one double argument //{{{
bool DoubleOption(int argc, char *argv[], char *opt, double *value) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - missing argument
      if ((i+1) >= argc) {
        strcpy(ERROR_MSG, "missing numeric argument");
        PrintErrorOption(opt);
        return true;
      }
      // Error - non-numeric
      if (!IsPosRealNumber(argv[i+1], value)) {
        strcpy(ERROR_MSG, "argument must be positive number");
        PrintErrorOption(opt);
        return true;
      }
    }
  }
  return false;
} //}}}

// TODO: join with IntegerOption and add max number of arguments as a parameter
// general option with multiple integer arguments (up to 100) //{{{
bool MultiIntegerOption(int argc, char *argv[], char *opt,
                        int *count, int *values) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read integers
      int arg = i+1+n;
      while ((arg) < argc && argv[arg][0] != '-') {
        // Error - non-numeric or missing argument
        long val;
        if (!IsIntegerNumber(argv[arg], &val)) {
          strcpy(ERROR_MSG, "each argument must be non-negative whole number");
          PrintErrorOption(opt);
          return true;
        }
        values[n] = val;
        n++;
        arg = i+1+n;
        // warning - too many numeric arguments
        if (n == 100) {
          strcpy(ERROR_MSG, "too many arguments; only the first 100 used");
          PrintErrorOption(opt);
          *count = n;
          return true;
        }
      }
      *count = n;
    }
  }
  return false;
} //}}}

// TODO: join with DoubleOption and add max number of arguments as a parameter
// general option with multiple double arguments (up to 100) //{{{
bool MultiDoubleOption(int argc, char *argv[], char *opt,
                       int *count, double *values) {
  *count = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read doubles
      int arg = i+1+n;
      while (arg < argc && IsRealNumber(argv[arg], &values[n])) {
        values[n] = atof(argv[arg]);
        n++;
        arg = i + 1 + n;
        // warning - too many numeric arguments
        if (n == 100) {
          strcpy(ERROR_MSG, "too many arguments; only the first 100 used");
          PrintErrorOption(opt);
          *count = n;
          return true;
        }
      }
      *count = n;
    }
  }
  return false;
} //}}}

// general option with filename and integer(s) arguments //{{{
bool FileIntsOption(int argc, char *argv[], char *opt, int *values,
                    int *count, char *file) {
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
        if (n == 100) {
          strcpy(ERROR_MSG, "too many arguments; only the first 100 used");
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
  *count = n;
  return false;
} //}}}

// TODO: join with FileIntsOption(), adding max integers (incl. 0) as parameter
// general option with a filename argument //{{{
bool FileOption(int argc, char *argv[], char *opt, char *name, int length) {
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
        strcpy(ERROR_MSG, "something wrong with snprintf");
        PrintError();
        return true;
      }
    }
  }
  return false;
} //}}}

#if 0 //{{{
// TODO redo
// MoleculeTypeOption() //{{{
/**
 * Generic option for molecule type that can take one
 * argument. The option is an argument of this function.
 */
bool MoleculeTypeOption(int argc, char *argv[], char *opt, int *moltype,
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
bool MoleculeTypeOption2(int argc, char *argv[], char *opt, int *moltype,
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
bool MoleculeTypeIntOption(int argc, int i, char *argv[], char *opt,
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
