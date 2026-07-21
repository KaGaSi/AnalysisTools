#include "Options.h"
#include "Arrays.h"
#include "Errors.h"
#include "General.h"
#include "Globals.h"
#include "System.h"
#include <stdio.h>
#include <string.h>

// STATIC DECLARATIONS
static void SilentOption(const int argc, char **argv,
                         bool *verbose, bool *silent);
static bool VersionOption(const int argc, char **argv);
// some repeated warnings/errors
static void MissingFilenameError(const int argc, char **argv,
                                 const char *opt, const int i);
static bool TooManyArgsWarn(const int max, const int n,
                            const char *opt, int *count);
static void ArgumentNumberErr(const int count, const int n, const char *opt);
static void ArgumentMissingErr(const int n, const char *opt);

// print help //{{{
void Help(const bool error, const struct HelpHelp help,
          const struct OptSpec *options) {
  // pick output stream
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
  }
  // compute max width
  size_t maxlen = 0;
  for (size_t i = 0; options[i].opt; i++) {
    // 256 is the same as used in OptionCheck to ensure the string are fine
    size_t len = strnlen(options[i].opt, 256);
    if (options[i].extra) {
      len += strnlen(options[i].extra, 256) + 1;
    }
    if (len > maxlen) {
      maxlen = len;
    }
  }
  int width = (int)maxlen + 5;
  // print stuff
  // a) description (for no error case)
  if (!error) {
    fprintf(ptr, "%s\n\n", help.description);
  }
  // b) usage
  fprintf(ptr, "%s\n\n", help.usage);
  // c) mandatory arguments
  for (size_t i = 0; options[i].opt; i++) {
    if (options[i].kind != OPT_ARG) {
      continue;
    }
    char line[LINE];
    if (options[i].extra) {
      snprintf(line, sizeof(line), "%s %s", options[i].opt, options[i].extra);
    } else {
      snprintf(line, sizeof(line), "%s", options[i].opt);
    }
    fprintf(ptr, "%-*s%s\n", width + 2, line, options[i].desc);
  }
  fprintf(ptr, "[options]\n");
  // d) extra arguments
  for (size_t i = 0; options[i].opt; i++) {
    if (options[i].kind != OPT_EXTRA) {
      continue;
    }
    char line[LINE];
    if (options[i].extra) {
      snprintf(line, sizeof(line), "%s %s", options[i].opt, options[i].extra);
    } else {
      snprintf(line, sizeof(line), "%s", options[i].opt);
    }
    fprintf(ptr, "  %-*s%s\n", width, line, options[i].desc);
  }
  // e) common arguments
  for (size_t i = 0; options[i].opt; i++) {
    if (options[i].kind != OPT_COMMON) {
      continue;
    }
    char line[LINE];
    if (options[i].extra) {
      snprintf(line, sizeof(line), "%s %s", options[i].opt, options[i].extra);
    } else {
      snprintf(line, sizeof(line), "%s", options[i].opt);
    }
    fprintf(ptr, "  %-*s%s\n", width, line, options[i].desc);
  }
} //}}}

// version/help printing and initial check of provided options //{{{
// TODO: check for multiple options; warn that the first one is used
int OptionCheck(const int argc, char **argv, const bool check_extra,
                const struct HelpHelp desc, const struct OptSpec *opts) {
  // simple check the opts struct is filled in properly
  // test for the {nullptr}
  if (opts[desc.all].opt) {
    err_msg("last opts[] must be {nullptr}; or wrong count desc.opts");
    PrintError();
    exit(1);
  }
  // test strings are null-terminated & count mandatory arguments
  for (size_t i = 0; opts[i].opt; i++) {
    // ...uses 256 as a reasaonable maximum length (actully, unnecessarily long)
    if ((opts[i].extra && strnlen(opts[i].extra, 256) == 256) ||
        strnlen(opts[i].opt, 256) == 256 ||
        strnlen(opts[i].desc, 256) == 256) {
      err_msg("unterminated string in opts[] (or longer than 256 characters)");
      PrintError();
      exit(1);
    }
    // .kind must be present
    if (!opts[i].opt && !opts[i].kind) {
      if (snprintf(ERROR_MSG, LINE, "missiong .kind in %s%s%s",
                   ErrYellow(), opts[i].opt, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
  }
  // --version option?
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  // --help option?
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      Help(false, desc, opts);
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
  if (count < desc.args) {
    ErrorArgNumber(count, desc.args);
    PrintCommand(stderr, argc, argv);
    Help(true, desc, opts);
    exit(1);
  }
  // all options exist?
  for (int i = (count+1); i < argc; i++) {
    bool valid = false;
    for (int j = 0; opts[j].opt; j++) {
      double value;
      // check if cli argument is valid
      if (argv[i][0] != '-' || // argument to an option
          IsRealNumber(argv[i], &value) || // negative number
          strcmp(argv[i], opts[j].opt) == 0) { // an option
        valid = true;
        break;
      }
    }
    if (!valid) {
      ErrorOption(argv[i]);
      PrintCommand(stderr, argc, argv);
      Help(false, desc, opts);
      exit(1);
    }
  }
  // warn if extra arguments (between required ones and options)
  if (check_extra && desc.args != count) {
    char extra[LINE] = "\0";
    for (int i = (desc.args + 1); i <= count; i++) {
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
void CommonHelp(const bool error, const int n,
                const char option[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
  }
  for (int i = 0; i < n; i++) {
    if (strcmp(option[i], "-i") == 0) {
      fprintf(ptr, "  -i <stru>         input structure file if different "
                   "than the coordinate file\n");
    } else if (strcmp(option[i], "-ft") == 0) {
      fprintf(ptr, "  -ft <type>        file type of coordinate file: "
                   "vtf/vsf/vcf, xyz, data, ltrj "
                   "(skips extension-based detection)\n");
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
COMMON_OPT CommonOptions(const int argc, char **argv, const SYS_FILES f) {
  COMMON_OPT opt;
  opt.start = 1;
  opt.end = -1;
  opt.skip = 0;
  // -v option - verbose output
  opt.verbose = BoolOption(argc, argv, COMMON_OPTS[C_VERBOSE].opt);
  // --silent option - silent mode
  SilentOption(argc, argv, &opt.verbose, &opt.silent);
  // starting/ending timestep
  if (OneNumberOption(argc, argv, COMMON_OPTS[C_ST].opt, &opt.start, 'i') &&
      opt.start <= 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption(COMMON_OPTS[C_ST].opt);
    exit(1);
  }
  if (OneNumberOption(argc, argv, COMMON_OPTS[C_E].opt, &opt.end, 'i') &&
      opt.end <= 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption(COMMON_OPTS[C_E].opt);
    exit(1);
  }
  if (opt.end != -1 && opt.start > opt.end) {
    snprintf(ERROR_MSG, LINE, "starting step (%s%d%s) lower than ending step "
             "(%s%d%s)", ErrYellow(), opt.start, ErrRed(),
             ErrYellow(), opt.end, ErrRed());
    PrintErrorOption("-st/-e");
    exit(1);
  }
  // number of timesteps to skip per one used
  if (OneNumberOption(argc, argv, COMMON_OPTS[C_SK].opt, &opt.skip, 'i') &&
      opt.skip < 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption(COMMON_OPTS[C_SK].opt);
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
// tag bead/molecule types true/false //{{{
bool TypeOption(const int argc, char **argv, const char opt[], const int mode,
                const bool use, bool *flag, const SYSTEM System) {
  char text[9];
  if (mode == 'b') {
    s_strcpy(text, "bead", 9);
  } else if (mode == 'm') {
    s_strcpy(text, "molecule", 9);
  } else {
    err_msg("TypeOption(): mode must be 'b' or 'm'");
    PrintError();
    exit(1);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // pointer to the proper function
      int (*func)(const char *, const SYSTEM);
      if (mode == 'b') {
        func = &FindBeadType;
      } else {
        func = &FindMoleculeName;
      }
      // go over all the beads
      int pos = i;
      while (++pos < argc && argv[pos][0] != '-') {
        int type = func(argv[pos], System);
        if (type == -1) {
          err_msg("non-existent name");
          PrintErrorOption(opt);
          if (mode == 'b') {
            ErrorBeadType(argv[pos], System);
          } else {
            ErrorMoleculeType(argv[pos], System);
          }
          exit(1);
        }
        flag[type] = use;
      }
      if (pos == (i + 1)) {
        if (snprintf(ERROR_MSG, LINE, "at least one %s type is required",
                     text) < 0) {
          ErrorSnprintf();
        }
        PrintErrorOption(opt);
        exit(1);
      }
      return true; // option is present
    }
  }
  return false; // option is not present
} //}}}
// tag bead/molecule type pairs true/false //{{{
bool TypeOptionPair(const int argc, char **argv, const char opt[],
                    const int mode, const bool use, ArrNDb *flag,
                    const SYSTEM System) {
  char text[9];
  if (mode == 'b') {
    s_strcpy(text, "bead", 9);
  } else if (mode == 'm') {
    s_strcpy(text, "molecule", 9);
  } else {
    err_msg("TypeOption(): mode must be 'b' or 'm'");
    PrintError();
    exit(1);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      // count type names
      int pos = i, count_name = 0;
      while (++pos < argc && argv[pos][0] != '-') {
        count_name++;
      }
      // errors //{{{
      // a) no name supplied
      if (count_name == 0) {
        if (snprintf(ERROR_MSG, LINE, "at least one %s type pair is required",
                     text) < 0) {
          ErrorSnprintf();
        }
        PrintErrorOption(opt);
        exit(1);
      // b) odd number of names supplied
      } else if (count_name == 0 || (count_name % 2) != 0) {
        if (snprintf(ERROR_MSG, LINE, "even number of %s types is required",
                     text) < 0) {
          ErrorSnprintf();
        }
        PrintErrorOption(opt);
        exit(1);
      } //}}}
      // pointer to the proper function
      int (*func)(const char *, const SYSTEM);
      if (mode == 'b') {
        func = &FindBeadType;
      } else {
        func = &FindMoleculeName;
      }
      // go over the names
      for (int j = 0; j < count_name; j += 2) {
        int type[2];
        for (int dd = 0; dd < 2; dd++) {
          int id = i + j + dd + 1;
          type[dd] = func(argv[id], System);
          if (type[dd] == -1) {
            err_msg("non-existent name");
            PrintErrorOption(opt);
            if (mode == 'b') {
              ErrorBeadType(argv[id], System);
            } else {
              ErrorMoleculeType(argv[id], System);
            }
            exit(1);
          }
        }
        SetArr2D(flag, type[0], type[1], use);
        SetArr2D(flag, type[1], type[0], use);
      }
      return true; // option is present
    }
  }
  return false; // option is not present
} //}}}

// general boolean option //{{{
bool BoolOption(const int argc, char **argv, const char *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      return true;
    }
  }
  return false;
} // }}}
// general option with multiple integer arguments (up to 'max') //{{{
bool NumbersOption(const int argc, char **argv, const int max, const char *opt,
                   int *count, void *values, const char type) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read integers
      int arg = i+1+n;
      while (arg < argc) {
        if (type == 'i') {
          long val;
          if (!IsIntegerNumber(argv[arg], &val)) {
            break;
          }
          if (TooManyArgsWarn(max, n + 1, opt, count)) {
            return true;
          }
          int *num = (int *)values;
          int *a = &num[n];
          *a = val;
        } else {
          double val;
          if (!IsRealNumber(argv[arg], &val)) {
            break;
          }
          if (TooManyArgsWarn(max, n + 1, opt, count)) {
            return true;
          }
          double *num = (double *)values;
          double *a = &num[n];
          *a = val;
        }
        n++;
        arg = i+1+n;
      }
      ArgumentMissingErr(n, opt);
      *count = n;
      return true;
    }
  }
  return false;
}
bool OneNumberOption(const int argc, char **argv,
                     const char *opt, void *value, const char type) {
  int count = 0;
  if (NumbersOption(argc, argv, 1, opt, &count, value, type)) {
    ArgumentNumberErr(count, 1, opt);
    return true; // option present
  }
  return false; // option not present
}
bool TwoNumbersOption(const int argc, char **argv,
                      const char *opt, void *value, const char type) {
  int count = 0;
  if (NumbersOption(argc, argv, 2, opt, &count, value, type)) {
    ArgumentNumberErr(count, 2, opt);
    return true; // option present
  }
  return false; // option not present
}
bool ThreeNumbersOption(const int argc, char **argv,
                        const char *opt, void *value, const char type) {
  int count = 0;
  if (NumbersOption(argc, argv, 3, opt, &count, value, type)) {
    ArgumentNumberErr(count, 3, opt);
    return true; // option present
  }
  return false; // option not present
} //}}}
// general option with filename and integer(s)/double(s) arguments //{{{
bool FileNumbersOption(const int argc, char **argv, const int min,
                       const int max, const char *opt, void *values,
                       int *count, char *file, const char type) {
  int n = 0;
  *count = 0;
  file[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      MissingFilenameError(argc, argv, opt, i);
      snprintf(file, LINE, "%s", argv[i+1]);
      // read numbers
      if (max == 0) {
        return true;
      } else {
        while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
          if (type == 'i') {
            long val;
            if (!IsIntegerNumber(argv[i+2+n], &val)) {
              err_msg("arguments must be non-negative numbers");
              goto error;
            }
            if (TooManyArgsWarn(max, n + 1, opt, count)) {
              return true;
            }
            int *num = (int *)values;
            int *a = &num[n];
            *a = val;
          } else {
            double val;
            if (!IsRealNumber(argv[i+2+n], &val)) {
              err_msg("arguments must be non-negative numbers");
              goto error;
            }
            if (TooManyArgsWarn(max, n + 1, opt, count)) {
              return true;
            }
            double *num = (double *)values;
            double *a = &num[n];
            *a = val;
          }
          n++;
        }
        if (n < min) {
          s_strcpy(ERROR_MSG, "not enough numeric arguments", LINE);
          goto error;
        }
      }
      *count = n;
      return true; // option present
    }
  }
  return false; // option not present
  error:
    PrintErrorOption(opt);
    exit(1);
} //}}}
// general option with filename //{{{
bool FileOption(const int argc, char **argv, const char *opt, char *file) {
  int trash;
  if (FileNumbersOption(argc, argv, 0, 0, opt, &trash, &trash, file, 'i')) {
    return true;
  }
  return false;
} //}}}

// option for output verbosity (--silent) //{{{
static void SilentOption(const int argc, char **argv,
                         bool *verbose, bool *silent) {
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
static bool VersionOption(const int argc, char **argv) {
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
// exit if file name is missing in an option where it's required //{{{
static void MissingFilenameError(const int argc, char **argv,
                                 const char *opt, const int i) {
  // Error - no output file name
  if ((i+1) >= argc || argv[i+1][0] == '-') {
    s_strcpy(ERROR_MSG, "missing file name "
             "(or the file name begins with a dash)", LINE);
    PrintErrorOption(opt);
    exit(1);
  }
} //}}}
// warn if too many arguments //{{{
static bool TooManyArgsWarn(const int max, const int n,
                            const char *opt, int *count) {
  if (n > max) {
    snprintf(ERROR_MSG, LINE, "too many arguments; only the first %d "
             "used", max);
    PrintErrorOption(opt);
    *count = max;
    return true;
  } else {
    return false;
  }
} //}}}
// exit if wrong number of arguments  //{{{
static void ArgumentNumberErr(const int count, const int n, const char *opt) {
  if (count != n) {
    snprintf(ERROR_MSG, LINE, "%d numeric argument(s) required", n);
    PrintErrorOption(opt);
    exit(1);
  }
} //}}}
// exit if missing arguments //{{{
static void ArgumentMissingErr(const int n, const char *opt) {
  if (n == 0) {
    s_strcpy(ERROR_MSG, "missing numeric argument(s)", LINE);
    PrintErrorOption(opt);
    exit(1);
  }
} //}}}
