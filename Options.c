#include "Options.h"

// TODO: make the bool Functions() into ints - -1 for option not present; 0 for
//       present; 1 for error. Then, use switch(Function()) in utils

// CommonHelp() //{{{
/**
 * Function to print help for common options, either for `-h` help option
 * or program error.
 */
void CommonHelp(bool error) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
  }
  fprintf(ptr, "   [standard options]\n");
  fprintf(ptr, "      -i <name>     input vtf structure file \
(default: traject.vs)\n");
  fprintf(ptr, "      --detailed    use name as weel as charge, mass, \
and radius to identfy bead types\n");
  fprintf(ptr, "      -v            verbose output\n");
  fprintf(ptr, "      --silent      no output (overrides verbose option)\n");
  fprintf(ptr, "      -h            print this help and exit\n");
  fprintf(ptr, "      --version     print version number and exit\n");
} //}}}

// CommonOptions() //{{{
/**
 * Function for options common to most of the utilities.
 */
void CommonOptions(int argc, char *argv[], char vsf_file[], int length,
                   bool *verbose, bool *silent, bool *detailed) {

  // -i <name> option - input structure file //{{{
  // test if '-i' option is there
  char name[LINE] = {'\0'};
  if (FileOption(argc, argv, "-i", name, length)) {
    exit(1);
  }
  // copy the name if '-i' option is present
  if (name[0] != '\0') {
    snprintf(vsf_file, LINE, "%s", name);
  }
  // test if structure file ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(vsf_file, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}
  // -v option - verbose output
  *verbose = BoolOption(argc, argv, "-v");
  // --silent option - silent mode
  SilentOption(argc, argv, verbose, silent);
  // --detailed option - base bead types on name, charge, mass, and radius
  *detailed = BoolOption(argc, argv, "--detailed");
} //}}}

// SilentOption() //{{{
/**
 * Option to not print anything to stdout (or at least no system
 * definitions and no Step: #). Overrides verbose option. Argument: `--silent`
 */
void SilentOption(int argc, char **argv, bool *verbose, bool *silent) {

  *silent = false;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--silent") == 0) {
      *verbose = false;
      *silent = true;

      break;
    }
  }
} //}}}

// VersionOption() //{{{
/**
 * Option to print version number of the program suite.
 */
bool VersionOption(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fprintf(stdout, "AnalysisTools by Karel Šindelka (KaGaSi), version %s", VERSION);
      fprintf(stdout, " (released %s)\n", DATE);
      fprintf(stdout, "Download at https://github.com/KaGaSi/AnalysisTools/releases/tag/v%s\n", VERSION);
      return true;
    }
  }

  return false;
} //}}}

bool ExcludeOption(int argc, char *argv[], SYSTEM *System) { //{{{
  // set all molecules to use //{{{
  for (int i = 0; i < (*System).Count.MoleculeType; i++) {
    (*System).MoleculeType[i].Use = true;
  } //}}}
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-x") == 0) {
      // wrong argument to -x option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "-x");
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing an argument (or molecule name beginning with a dash)\n\n");
        ColourReset(STDERR_FILENO);
        exit(1);
      } //}}}
      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {
        int type = FindMoleculeType(argv[i+1+j], *System);
        if (type == -1) { // is it in vsf?
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "-x");
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - non-existent ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", argv[i+1+j]);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " molecule\n\n");
          ColourReset(STDERR_FILENO);
          ErrorMoleculeType(*System);
          return(true);
        } else {
          // exclude that molecule
          (*System).MoleculeType[type].Use = false;
        }
        j++;
      }
    }
  }

  return(false);
} //}}}

// JoinCoorOption() //{{{
/**
 * Option whether to join aggregates and save joined coordinates into a
 * specified file. Arguments: `-j <joined.vcf>`
 */
bool JoinCoorOption(int argc, char **argv, char *joined_vcf) {

  joined_vcf[0] = '\0'; // no -j option
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      // wrong argument to -j option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "-j");
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing output file name");
        fprintf(stderr, " (or the file name begins with '-')\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      } //}}}
      snprintf(joined_vcf, LINE, "%s", argv[i+1]);
      // test if <joined.vcf> filename ends with '.vcf' //{{{
      char *dot = strrchr(joined_vcf, '.');
      if (!dot || (strcmp(dot, ".vcf") && strcmp(dot, ".vtf"))) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "-j");
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - file ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", joined_vcf);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " must have .vcf ending\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      } //}}}
    }
  }

  return(false);
} //}}}

// BeadTypeOption() //{{{
/**
 * Option to choose which bead types to use for calculation. If the option
 * is absent, all bead types are switched to the specified 'bool use' value.
 */
bool BeadTypeOption(int argc, char **argv, char *opt,
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
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - non-existent ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", argv[types]);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " bead type\n\n");
          ColourReset(STDERR_FILENO);
          ErrorBeadType(*System);
          return(true);
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
  return(false);
} // }}}

// BoolOption() //{{{
/**
 * Function for any boolean option (i.e. without argument).
 */
bool BoolOption(int argc, char **argv, char *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      return(true);
    }
  }
  return(false);
} // }}}

// IntegerOption() //{{{
/**
 * Function for any option with integer argument.
 */
bool IntegerOption(int argc, char **argv, char *opt, int *value) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - missing argument
      if ((i+1) >= argc) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing numeric argument\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      }
      // Error - non-numeric
      if (!IsInteger_old(argv[i+1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - argument must be non-negative whole number\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      }
      *value = atoi(argv[i+1]);
    }
  }

  return(false);
} //}}}

// DoubleOption() //{{{
/**
 * Function for any option with double argument.
 */
bool DoubleOption(int argc, char **argv, char *opt, double *value) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - missing argument
      if ((i+1) >= argc) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing numeric argument\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      }
      // Error - non-numeric
      if (!IsPosReal_old(argv[i+1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - argument must be positive number\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      }
      *value = atof(argv[i+1]);
    }
  }

  return(false);
} //}}}

// TODO: join with IntegerOption and add max number of arguments as a parameter
// MultiIntegerOption() //{{{
/**
 * Function for any option with two or more (up to 100) integer arguments.
 */
bool MultiIntegerOption(int argc, char **argv, char *opt,
                        int *count, int *values) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read integers
      int arg = i+1+n;
      while ((arg) < argc && argv[arg][0] != '-') {
        // Error - non-numeric or missing argument
        if (!IsInteger_old(argv[arg])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - all argument(s) must be");
          fprintf(stderr, " non-negative whole number(s)\n\n");
          ColourReset(STDERR_FILENO);
          return true;
        }
        values[n] = atoi(argv[arg]);
        n++;
        arg = i+1+n;
        // warning - too many numeric arguments
        if (n == 100) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - too many arguments; only the first 100 used\n\n");
          ColourReset(STDERR_FILENO);
          *count = n;
          return true;
        }
      }
      *count = n;
    }
  }

  return false;
} //}}}

// MultiDoubleOption() //{{{
/**
 * Function for any option with two or more (up to 100) double arguments.
 */
bool MultiDoubleOption(int argc, char **argv, char *opt,
                       int *count, double *values) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read doubles
      int arg = i+1+n;
//    while ((arg) < argc && argv[arg][0] != '-') {
      // A = arg < argc; B = argv[arg][0] == '-'; C = IsReal(argv[arg])
      // A B C | we want | A and ((B and C) or (!B and C)) = A and C
      // ----------------|------------------------------------------
      // 1 1 1 | 1       | 1
      // 1 1 0 | 0       | 0
      // 1 0 1 | 1       | 1
      // 0 0 0 | 0       | 0
      // 0 1 1 | 0       | 0
      // 0 1 0 | 0       | 0
      while (arg < argc && IsReal_old(argv[arg])) { // see expression table up
        // Error - non-numeric argument
//      if (!IsPosReal(argv[arg])) {
//        ColourText(STDERR_FILENO, RED);
//        fprintf(stderr, "\nError: ");
//        ColourText(STDERR_FILENO, YELLOW);
//        fprintf(stderr, "%s", opt);
//        ColourText(STDERR_FILENO, RED);
//        fprintf(stderr, " - argument(s) must be positive number(s)\n\n");
//        ColourReset(STDERR_FILENO);
//        return true;
//      }
        values[n] = atof(argv[arg]);
        n++;
        arg = i+1+n;
        // warning - too many numeric arguments
        if (n == 100) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - too many arguments; only the first 100 used\n\n");
          ColourReset(STDERR_FILENO);
          *count = n;
          return true;
        }
      }
      *count = n;
    }
  }

  return false;
} //}}}

// FileIntsOptions() //{{{
/**
 * Function for any option with a file name and up to 100 integer
 * arguments.
 */
bool FileIntsOption(int argc, char **argv, char *opt, int *values,
                    int *count, char *file) {

  int n = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - no output file name
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing output file name");
        fprintf(stderr, " (or the file name begins with '-')\n\n");
        ColourReset(STDERR_FILENO);
        return false;
      }
      snprintf(file, LINE, "%s", argv[i+1]);
      // read integers
      while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
        // Error - non-numeric or missing argument
        if (!IsInteger_old(argv[i+2+n])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - argument must be non-negative whole number\n\n");
          ColourReset(STDERR_FILENO);
          return true;
        }
        values[n] = atoi(argv[i+2+n]);
        n++;
        // warning - too many numeric arguments
        if (n == 100) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - too many arguments; only the first 100 used\n\n");
          ColourReset(STDERR_FILENO);
          *count = n;
          return true;
        }
      }
      if (n == 0) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing numeric argument(s)\n\n");
        ColourReset(STDERR_FILENO);
        return true;
      }
    }
  }
  *count = n;

  return false;
} //}}}

// FileOption() //{{{
/**
 * Generic option for file name. The option is an argument of this function.
 */
bool FileOption(int argc, char **argv, char *opt,
                char *name, int length) {
  name[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // wrong argument to the option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing output file name");
        fprintf(stderr, " (or the file name begins with '-')\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      }
      snprintf(name, LINE, "%s", argv[i+1]);
    }
  }
  return(false);
} //}}}

// StartEndTime() //{{{
/**
 * Options for starting and ending timesteps.
 */
void StartEndTime(int argc, char **argv, int *start, int *end) {
  *start = 1;
  if (IntegerOption(argc, argv, "-st", start)) {
    exit(1);
  }
  *end = -1;
  if (IntegerOption(argc, argv, "-e", end)) {
    exit(1);
  }
  ErrorStartEnd(*start, *end);
} //}}}

#if 0
// TODO redo
// MoleculeTypeOption() //{{{
/**
 * Generic option for molecule type that can take one
 * argument. The option is an argument of this function.
 */
bool MoleculeTypeOption(int argc, char **argv, char *opt, int *moltype,
                        COUNTS Counts, MOLECULETYPE **MoleculeType) {

  *moltype = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing output file name");
        fprintf(stderr, " (or the file name begins with '-')\n\n");
        ColourReset(STDERR_FILENO);
        return(true);
      } //}}}
      *moltype = FindMoleculeType_old(argv[i+1], Counts, *MoleculeType);
      if (*moltype == -1) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - non-existent ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", argv[i+1]);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " molecule\n\n");
        ColourReset(STDERR_FILENO);
        ErrorMoleculeType_old(Counts, *MoleculeType);
        return(true);
      }
    }
  }

  return(false);
} //}}}
// MoleculeTypeOption2() //{{{
/**
 * Generic option for molecule types that can take multiple arguments. The
 * option is an argument of this function.
 */
bool MoleculeTypeOption2(int argc, char **argv, char *opt, int *moltype,
                         COUNTS Counts, MOLECULETYPE **MoleculeType) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // set all moltypes not to be used
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        moltype[j] = 0;
      }
      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", opt);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing molecule name");
        fprintf(stderr, " (or the name begins with '-')\n\n");
        ColourReset(STDERR_FILENO);
        ErrorMoleculeType_old(Counts, *MoleculeType);
        return true;
      } //}}}
      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {
        int type = FindMoleculeType_old(argv[i+1+j], Counts, *MoleculeType);
        if (type == -1) { // is argv[i+1+j] in vsf?
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - non-existent ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", argv[i+1+j]);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " molecule\n\n");
          ColourReset(STDERR_FILENO);
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
bool MoleculeTypeIntOption(int argc, int i, char **argv, char *opt,
                           int *moltype, int *value, COUNTS Counts,
                           MOLECULETYPE *MoleculeType) {
  *moltype = -1;
  if (strcmp(argv[i], opt) == 0) {
    // Error - missing or wrong arguments //{{{
    if ((i+2) >= argc || argv[i+1][0] == '-' || !IsInteger_old(argv[i+2])) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", opt);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - two arguments required (<mol name> <int>)");
      ColourReset(STDERR_FILENO);
      return true;
    } //}}}
    *moltype = FindMoleculeType_old(argv[i+1], Counts, MoleculeType);
    if (*moltype == -1) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", opt);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - non-existent ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", argv[i+1]);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " molecule\n\n");
      ColourReset(STDERR_FILENO);
      ErrorMoleculeType_old(Counts, MoleculeType);
      return true;
    }
    *value = atoi(argv[i+2]);
  }
  return false;
} //}}}
// TODO remove
// BeadTypeOption_old() //{{{
/**
 * Option to choose which bead types to use for calculation. If the option
 * is absent, all bead types are switched to the specified 'bool use' value.
 */
bool BeadTypeOption_old(int argc, char **argv, char *opt, bool use,
                    COUNTS Counts, BEADTYPE **BeadType) {

  // specify what bead types to use - either specified by 'opt' or use all
  int types = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      types = i; // positon of the '-bt' argument in command
      // <type names> - names of bead types to save
      while (++types < argc && argv[types][0] != '-') {
        int type = FindBeadType_old(argv[types], Counts, *BeadType);
        if (type == -1) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", opt);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - non-existent ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", argv[types]);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " bead type\n\n");
          ColourReset(STDERR_FILENO);
          ErrorBeadType_old(Counts, *BeadType);
          return(true);
        }
        (*BeadType)[type].Use = true;
      }
    }
  }
  if (types == -1) {
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      (*BeadType)[i].Use = use;
    }
  }

  return(false);
} // }}}
// ExcludeOption_old() //{{{
/**
 * Option to exclude specified molecule types from calculations. Gives
 * specified molecule types `Use = false` and the rest `Use = true`.
 * Arguments: `-x <name(s)>`
 */
bool ExcludeOption_old(int argc, char **argv, COUNTS Counts,
                   MOLECULETYPE **MoleculeType) {

  // set all molecules to use //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = true;
  } //}}}
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-x") == 0) {
      // wrong argument to -x option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "-x");
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing an argument (or molecule name beginning with a dash)\n\n");
        ColourReset(STDERR_FILENO);
        exit(1);
      } //}}}
      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {
        int type = FindMoleculeType_old(argv[i+1+j], Counts, *MoleculeType);
        if (type == -1) { // is it in vsf?
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "-x");
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - non-existent ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", argv[i+1+j]);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " molecule\n\n");
          ColourReset(STDERR_FILENO);
          ErrorMoleculeType_old(Counts, *MoleculeType);
          return(true);
        } else {
          // exclude that molecule
          (*MoleculeType)[type].Use = false;
        }
        j++;
      }
    }
  }

  return(false);
} //}}}
#endif
