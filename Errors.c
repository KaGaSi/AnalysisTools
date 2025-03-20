#include "Errors.h"
#include "General.h"
#include "Globals.h"

static void PrintLine(FILE *f, const char *colour1, const char *colour2);

// print 'WARNING - <ERROR_MSG>\n' in cyan //{{{
void PrintWarning() {
  fprintf(stderr, "\n  %sWARNING - %s%s\n",
          ErrCyan(), ERROR_MSG, ErrColourReset());
} //}}}
// print 'ERROR - <ERROR_MSG>\n' in red //{{{
void PrintError() {
  fprintf(stderr, "\n  %sERROR - %s%s\n",
          ErrRed(), ERROR_MSG, ErrColourReset());
} //}}}
// print 'ERROR <option> - <ERROR_MSG>\n' in red and yellow //{{{
void PrintErrorOption(const char *opt) {
  fprintf(stderr, "\n  %sERROR: %s%s%s - %s%s\n",
          ErrRed(), ErrYellow(), opt, ErrRed(), ERROR_MSG, ErrColourReset());
} //}}}
// print 'WARNING <option> - <ERROR_MSG>\n' in cyan and yellow //{{{
void PrintWarnOption(const char *opt) {
  fprintf(stderr, "\n  %sWARNING: %s%s%s - %s%s\n",
          ErrCyan(), ErrYellow(), opt, ErrCyan(), ERROR_MSG, ErrColourReset());
} //}}}
// print 'ERROR - <ERROR_MSG>\nFile <file(s)>\n' //{{{
void PrintErrorFile(const char *file1, const char *file2, const char *file3) {
  PrintError();
  ErrorPrintFile(file1, file2, file3);
  putc('\n', stderr);
} //}}}
// print 'WARNING - <ERROR_MSG>\nFile <file(s)>\n' //{{{
void PrintWarnFile(const char *file1, const char *file2, const char *file3) {
  PrintWarning();
  WarnPrintFile(file1, file2, file3);
  putc('\n', stderr);
} //}}}
// print 'ERROR - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>' //{{{
void PrintErrorFileLine(const char *file, const int count) {
  PrintError();
  ErrorPrintFile(file, "", "");
  fprintf(stderr, "%s, line %s%d%s:\n", ErrRed(), ErrYellow(), count, ErrRed());
  ErrorPrintLine();
} //}}}
// print 'WARNING - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>' //{{{
void PrintWarnFileLine(const char *file, const int count) {
  PrintWarning();
  WarnPrintFile(file, "", "");
  fprintf(stderr, "%s, line %s%d%s:\n", ErrCyan(),
          ErrYellow(), count, ErrCyan());
  WarnPrintLine();
} //}}}
// print 'File <name(s)>' in colours //{{{
void WarnPrintFile(const char *file1, const char *file2,const  char *file3) {
  fprintf(stderr, "%sFile(s) %s", ErrCyan(), ErrYellow());
  const char *f[3];
  f[0] = file1;
  f[1] = file2;
  f[2] = file3;
  int n_files = 0;
  for (int i = 0; i < 3; i++) {
    if (f[i][0] != '\0') {
      if (n_files > 0) {
        fprintf(stderr, "%s, %s", ErrCyan(), ErrYellow());
      }
      fprintf(stderr, "%s", f[i]);
      n_files++;
    }
  }
  fprintf(stderr, "%s", ErrColourReset());
}
void ErrorPrintFile(const char *file1, const char *file2, const char *file3) {
  fprintf(stderr, "%sFile(s) %s", ErrRed(), ErrYellow());
  const char *f[3];
  f[0] = file1;
  f[1] = file2;
  f[2] = file3;
  if (strcmp(f[0], f[1]) == 0) {
    f[1] = "";
  } else if (strcmp(f[0], f[2]) == 0) {
    f[2] = "";
  }
  if (strcmp(f[1], f[2]) == 0) {
    f[2] = "";
  }
  int n_files = 0;
  for (int i = 0; i < 3; i++) {
    if (f[i][0] != '\0') {
      if (n_files > 0) {
        fprintf(stderr, "%s, %s", ErrRed(), ErrYellow());
      }
      fprintf(stderr, "%s", f[i]);
      n_files++;
    }
  }
  fprintf(stderr, "%s", ErrColourReset());
}
//}}}
// print 'Line: <line>|(blank)' in given colours //{{{
static void PrintLine(FILE *f, const char *colour1, const char *colour2) {
  if (words == 0) {
    fprintf(f, "%s(blank)\n%s", Colour(f, colour1), Colour(f, C_RESET));
  } else {
    fputs(Colour(f, colour2), f);
    for (int i = 0; i < words; i++) {
      if (i != 0) {
        putc(' ', f);
      }
      fprintf(f, "%s", split[i]);
    }
    fputs(Colour(f, C_RESET), f);
    // print line break if the last split[] doesn't end with one
    if (split[words-1][strnlen(split[words-1], SPL_LEN) - 1] != '\n') {
      putc('\n', f);
    }
  }
}
void ErrorPrintLine() {
  PrintLine(stderr, RED, YELLOW);
}
void WarnPrintLine() {
  PrintLine(stderr, CYAN, YELLOW);
}
//}}}
// print 'ERROR - premature end of file\n<file>\n' //{{{
void ErrorEOF(const char *file, char *msg) {
  if (msg[0] != '\0') {
    if (strnlen(msg, LINE) >= (LINE - 24)) {
      msg[LINE-24] = '\0';
    }
    snprintf(ERROR_MSG, LINE, "%s; premature end of file", msg);
  } else {
    err_msg("premature end of file");
  }
  PrintError();
  ErrorPrintFile(file, "", "");
  putc('\n', stderr);
} //}}}
// snprintf - just to shut up compiler warnings; should never trigger //{{{
void ErrorSnprintf() {
  err_msg("something went wrong with snprintf()");
  PrintError();
  exit(1);
} //}}}
// wrong number of commandline arguments //{{{
void ErrorArgNumber(const int count, const int need) {
  err_msg("insufficient number of arguments");
  PrintError();
  fprintf(stderr, "%ssupplied: %s%d%s, needed: %s%d%s\n", ErrRed(),
          ErrYellow(), count, ErrRed(), ErrYellow(), need, ErrColourReset());
} //}}}
// wrong file extension //{{{
int ErrorExtension(const char *file, const int number,
                   const char extension[][EXTENSION]) {
  char *dot = strrchr(file, '.');
  for (int i = 0; i < number; i++) {
    if (dot && strcasecmp(dot, extension[i]) == 0) {
      return i;
    }
  }
  err_msg("incorrect file extension");
  PrintError();
  ErrorPrintFile(file, "", "");
  fprintf(stderr, "%s; allowed extensions:", ErrRed());
  for (int i = 0; i < (number - 1); i++) {
    fprintf(stderr, " %s%s%s,", ErrYellow(), extension[i], ErrRed());
  }
  fprintf(stderr, " %s%s%s\n", ErrYellow(), extension[number - 1],
          ErrColourReset());
  return -1;
} //}}}
void ErrorOption(const char *option) { //{{{
  err_msg("non-existent option");
  PrintErrorOption(option);
} //}}}
void ErrorNaN(const char *option) { //{{{
  err_msg("non-numeric argument");
  PrintErrorOption(option);
} //}}}
void ErrorBeadType(const char *name, const SYSTEM System) { //{{{
  fprintf(stderr, "%sillegal name: %s%s%s\n",
          ErrRed(), ErrYellow(), name, ErrRed());
  fprintf(stderr, "Possible names: %s%s\n", ErrYellow(),
          System.BeadType[0].Name);
  for (int i = 1; i < System.Count.BeadType; i++) {
    fprintf(stderr, "                %s\n", System.BeadType[i].Name);
  }
  fprintf(stderr, "%s\n", ErrColourReset());
} //}}}
void ErrorMoleculeType(const char *name, const SYSTEM System) { //{{{
  fprintf(stderr, "%sillegal name: %s%s%s\n",
          ErrRed(), ErrYellow(), name, ErrRed());
  fprintf(stderr, "Possible molecule names: %s%s\n", ErrYellow(),
          System.MoleculeType[0].Name);
  for (int i = 1; i < System.Count.MoleculeType; i++) {
    fprintf(stderr, "                         %s\n",
            System.MoleculeType[i].Name);
  }
  fprintf(stderr, "%s\n", ErrColourReset());
} //}}}
// warn the system is not electrically neutral //{{{
void WarnChargedSystem(const SYSTEM System, const char *file1,
                       const char *file2, char const *file3) {
  double charge = 0;
  for (int i = 0; i < System.Count.BeadType; i++) {
    BEADTYPE *bt = &System.BeadType[i];
    // do nothing if at least one bead type had undefined charge
    if (bt->Charge == CHARGE || bt->Charge == HIGHNUM) {
      return;
    }
    charge += bt->Charge * bt->Number;
  }
  if (fabs(charge) > 0.00001) {
    err_msg("system with net electric charge");
    PrintWarning();
    WarnPrintFile(file1, file2, file3);
    fprintf(stderr, "%s; %sq = %lf%s\n", ErrCyan(), ErrYellow(), charge,
            ErrColourReset());
  }
} //}}}
void ErrorStartEnd(const int start, const int end) { //{{{
  if (end != -1 && start > end) {
    snprintf(ERROR_MSG, LINE, "starting step (%s%d%s) lower than ending step "
             "(%s%d%s)", ErrYellow(), start, ErrRed(),
             ErrYellow(), end, ErrRed());
    PrintErrorOption("-st/-e");
    exit(1);
  }
} //}}}
// copy message to the ERROR_MSG array //{{{
void err_msg(const char *str) {
  s_strcpy(ERROR_MSG, str, LINE);
} //}}}
