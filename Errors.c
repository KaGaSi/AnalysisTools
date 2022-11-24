#include "Errors.h"
char ERROR_MSG[LINE];

// simple messages //{{{
// print 'WARNING - <ERROR_MSG>' in cyan //{{{
void PrintWarning() {
  fprintf(stderr, "\n  %sWARNING - %s%s\n", ErrCyan(),
          ERROR_MSG, ErrColourReset());
} //}}}
// print 'ERROR - <ERROR_MSG>' in red //{{{
void PrintError() {
  fprintf(stderr, "\n  %sERROR - %s%s\n", ErrRed(),
          ERROR_MSG, ErrColourReset());
} //}}}
// print 'WARNING: <option> - <ERROR_MSG>' in cyan and yellow //{{{
void PrintWarningOption(char *opt) {
  fprintf(stderr, "\n  %sWARNING: %s%s%s - %s%s\n", ErrCyan(),
          ErrYellow(), opt, ErrCyan(), ERROR_MSG, ErrColourReset());
} //}}}
// print 'ERROR: <option> - <ERROR_MSG>' in red and yellow //{{{
void PrintErrorOption(char *opt) {
  fprintf(stderr, "\n  %sERROR: %s%s%s - %s%s\n", ErrRed(),
          ErrYellow(), opt, ErrRed(), ERROR_MSG, ErrColourReset());
} //}}}
// print 'ERROR: - <ERROR_MSG>\nFile <file(s)>' //{{{
void PrintErrorFile(char file1[], char file2[], char file3[]) {
  PrintError();
  ErrorPrintFile(file1, file2, file3);
//putc('\n', stderr);
} //}}}
// print 'WARNING: - <ERROR_MSG>\nFile <file(s)>' //{{{
void PrintWarnFile(char file1[], char file2[], char file3[]) {
  PrintWarning();
  WarnPrintFile(file1, file2, file3);
//putc('\n', stderr);
} //}}}
// print 'ERROR: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>' //{{{
void PrintErrorFileLine(char file[], int count,
                        char *split[SPL_STR], int words) {
  PrintError();
  ErrorPrintFile(file, "\0", "\0");
  fprintf(stderr, "%s, line %s%d%s:\n", ErrRed(), ErrYellow(), count, ErrRed());
  ErrorPrintLine2(split, words);
} //}}}
// print 'WARNING: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>' //{{{
void PrintWarningFileLine(char file[], int count,
                          char *split[SPL_STR], int words) {
  PrintWarning();
  WarnPrintFile(file, "\0", "\0");
  fprintf(stderr, "%s, line %s%d%s:\n", ErrCyan(), ErrYellow(),
                                        count, ErrCyan());
  WarnPrintLine(split, words);
} //}}}
// 'File <name> (<optional name>)' in colours //{{{
void PrintFile(FILE *f, char file1[], char colour[]) {
  fprintf(f, "%sFile %s%s%s", Colour(f, colour), file1, Colour(f, YELLOW),
                              Colour(f, C_RESET));
}
void WarnPrintFile(char file1[], char file2[], char file3[]) {
  fprintf(stderr, "%sFile(s) %s", ErrCyan(), ErrYellow());
  char *f[3];
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
void ErrorPrintFile(char file1[], char file2[], char file3[]) {
  fprintf(stderr, "%sFile(s) %s", ErrRed(), ErrYellow());
  char *f[3];
  f[0] = file1;
  f[1] = file2;
  f[2] = file3;
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
void PrintLine2(FILE *f, char *split[SPL_STR], int words,
                char *colour1, char *colour2) {
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
    if (split[words-1][strlen(split[words-1])-1] != '\n') {
      putc('\n', f);
    }
  }
}
void ErrorPrintLine2(char *split[SPL_STR], int words) {
  PrintLine2(stderr, split, words, RED, YELLOW);
}
void WarnPrintLine(char *split[SPL_STR], int words) {
  PrintLine2(stderr, split, words, CYAN, YELLOW);
}
//}}}
 //}}}

// ErrorArgNumber() //{{{
/**
 * Error when insufficient number of arguments
 */
void ErrorArgNumber(int count, int need) {
  strcpy(ERROR_MSG, "insufficient number of arguments");
  PrintError();
  fprintf(stderr, "%ssupplied: %s%d%s, needed: %s%d%s\n", ErrRed(), ErrYellow(),
          count, ErrRed(), ErrYellow(), need, ErrColourReset());
} //}}}

// ErrorExtension() //{{{
/**
 * Error when missing or incorrect file extension
 */
int ErrorExtension(char *file, int number, char extension[][EXTENSION]) {
  char *dot = strrchr(file, '.');
  for (int i = 0; i < number; i++) {
    if (dot && strcmp(dot, extension[i]) == 0) {
      return i;
    }
  }
  strcpy(ERROR_MSG, "incorrect file extension");
  PrintError();
  ErrorPrintFile(file, "\0", "\0");
  fprintf(stderr, "%s; allowed extensions:", ErrRed());
  for (int i = 0; i < (number-1); i++) {
    fprintf(stderr, " %s%s%s,", ErrYellow(), extension[i], ErrRed());
  }
  fprintf(stderr, " %s%s%s\n", ErrYellow(),
          extension[number-1], ErrColourReset());
  return -1;
} //}}}

// ErrorNaN() //{{{
/**
 * Error when non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  strcpy(ERROR_MSG, "non-numeric argument");
  PrintErrorOption(option);
} //}}}

// ErrorOption() //{{{
/**
 * Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  strcpy(ERROR_MSG, "non-existent option");
  PrintErrorOption(option);
} //}}}

void ErrorBeadType(SYSTEM System) { //{{{
  fprintf(stderr, "     Possible bead names: %s\n", System.BeadType[0].Name);
  for (int i = 1; i < System.Count.BeadType; i++) {
    fprintf(stderr, "                          %s\n", System.BeadType[i].Name);
  }
  putc('\n', stderr);
} //}}}
void ErrorMoleculeType(SYSTEM System) { //{{{
  fprintf(stderr, "   Possible molecule names: %s\n",
          System.MoleculeType[0].Name);
  for (int i = 1; i < System.Count.MoleculeType; i++) {
    fprintf(stderr, "                        %s\n",
            System.MoleculeType[i].Name);
  }
  putc('\n', stderr);
} //}}}

// ErrorPrintLine() //{{{
/**
 * Print provided strings (array of strings generally created using
 * SplitLine()) to error output.
 */
void ErrorPrintLine(char split[SPL_STR][SPL_LEN], int words) {
  ColourChange(STDERR_FILENO, RED);
  if (words == 0) {
    fprintf(stderr, "   Blank line encountered");
  } else {
    fprintf(stderr, "   Wrong line:\n");
    ColourChange(STDERR_FILENO, YELLOW);
    for (int i = 0; i < words; i++) {
      if (i != 0) {
        putc(' ', stderr);
      }
      fprintf(stderr, "%s", split[i]);
    }
  }
  fprintf(stderr, "\n\n");
  ColourReset(STDERR_FILENO);
} //}}}

// WarnChargedSystem() //{{{
void WarnChargedSystem(SYSTEM System, char file1[], char file2[],
                       char file3[]) {
  double charge = 0;
  for (int i = 0; i < System.Count.BeadType; i++) {
    // do nothing if at least one bead type had undefined charge
    if (System.BeadType[i].Charge == CHARGE) {
      return;
    }
    charge += System.BeadType[i].Charge * System.BeadType[i].Number;
  }
  if (fabs(charge) > 0.00001) {
    strcpy(ERROR_MSG, "system with net electric charge");
    PrintWarning();
    WarnPrintFile(file1, file2, file3);
    fprintf(stderr, "%s; %sq = %lf%s\n", ErrCyan(),
            ErrYellow(), charge, ErrColourReset());
  }
} //}}}

// ErrorStartEnd() //{{{
/**
 * Error when ending timestep is higher than the starting one.
 */
void ErrorStartEnd(int start, int end) {
  if (end != -1 && start > end) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-st");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " and ");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-e");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - starting step (");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", start);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, ") is higher than ending step (");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", end);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, ")\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
} //}}}

// PrintLine() //{{{
/*
 * Print a single line from a file in a given colour (or 'blank line' in
 * different colour).
 */
void PrintLine(char split[SPL_STR][SPL_LEN], int words,
               int col_line, int col_blank) {
  if (words == 0) {
    fprintf(stderr, "Blank line");
    ColourReset(STDERR_FILENO);
  } else {
    for (int i = 0; i < words; i++) {
      if (i != 0) {
        putc(' ', stderr);
      }
      fputs(split[i], stderr);
    }
    ColourReset(STDERR_FILENO);
  }
  fputs("\n\n", stderr);
} //}}}

// PrintFileLine() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void PrintFileLine(char *file, int line,
                   char split[SPL_STR][SPL_LEN], int words) {
  ColourChange(STDERR_FILENO, RED);
  fputs("File ", stderr);
  PrintFile(stderr, file, YELLOW);
  ColourChange(STDERR_FILENO, RED);
  fputs(", line ", stderr);
  ColourChange(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", line);
  ColourChange(STDERR_FILENO, RED);
  fputs(":\n", stderr);
  ColourChange(STDERR_FILENO, YELLOW);
  for (int i = 0; i < words; i++) {
    if (i != 0) {
      putc(' ', stderr);
    }
    fputs(split[i], stderr);
  }
  fputs("\n\n", stderr);
  ColourReset(STDERR_FILENO);
} //}}}
// PrintFileLine2() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void PrintFileLine2(char *file, int line,
                   char *split[SPL_STR], int words) {
  ColourChange(STDERR_FILENO, RED);
  fputs("File ", stderr);
  PrintFile(stderr, file, YELLOW);
  ColourChange(STDERR_FILENO, RED);
  fputs(", line ", stderr);
  ColourChange(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", line);
  ColourChange(STDERR_FILENO, RED);
  fputs(":\n", stderr);
  ColourChange(STDERR_FILENO, YELLOW);
  for (int i = 0; i < words; i++) {
    if (i != 0) {
      putc(' ', stderr);
    }
    fputs(split[i], stderr);
  }
  fputs("\n\n", stderr);
  ColourReset(STDERR_FILENO);
} //}}}

// ErrorPrintFull() //{{{
void ErrorPrintFull(char *file, int line,
                    char split[SPL_STR][SPL_LEN], int words) {
  ErrorPrintError();
  PrintFileLine(file, line, split, words);
} //}}}
// ErrorPrintFull2() //{{{
void ErrorPrintFull2(char *file, int line,
                    char *split[SPL_STR], int words) {
  ErrorPrintError();
  PrintFileLine2(file, line, split, words);
} //}}}

// WarnStopReading() //{{{
/*
 * Warning when stopping file reading due to some error in the file.
 */
void WarnStopReading(char *vcf_file, int line_count, int step_count,
                     char split[SPL_STR][SPL_LEN], int words) {
  WarnPrintWarning();
  if (step_count-1 > 0) {
    ColourChange(STDERR_FILENO, CYAN);
    fputs("Last step read from file ", stderr);
    PrintFile(stderr, vcf_file, YELLOW);
    ColourChange(STDERR_FILENO, CYAN);
    fputs(": ", stderr);
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", step_count);
  } else {
    ColourChange(STDERR_FILENO, CYAN);
    fputs("No valid timestep in file ", stderr);
    PrintFile(stderr, vcf_file, YELLOW);
  }
  ColourChange(STDERR_FILENO, CYAN);
  fprintf(stderr, "; error at line ");
  ColourChange(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", line_count);
  ColourChange(STDERR_FILENO, CYAN);
  fputs(":\n", stderr);
//PrintLine(split, words, YELLOW, CYAN);
} //}}}
// WarnStopReading2() //{{{
/*
 * Warning when stopping file reading due to some error in the file.
 */
void WarnStopReading2(char *vcf_file, int line_count, int step_count,
                     char *split[SPL_STR], int words) {
  WarnPrintWarning();
  if (step_count-1 > 0) {
    ColourChange(STDERR_FILENO, CYAN);
    fputs("Last step read from file ", stderr);
//  PrintFile(vcf_file, YELLOW);
    ColourChange(STDERR_FILENO, CYAN);
    fputs(": ", stderr);
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", step_count);
  } else {
    ColourChange(STDERR_FILENO, CYAN);
    fputs("No valid timestep in file ", stderr);
//  PrintFile(vcf_file, YELLOW);
  }
  ColourChange(STDERR_FILENO, CYAN);
  fprintf(stderr, "; error at line ");
  ColourChange(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", line_count);
  ColourChange(STDERR_FILENO, CYAN);
  fputs(":\n", stderr);
//PrintLine2(split, words, YELLOW, CYAN);
} //}}}

// ErrorEOF() //{{{
void ErrorEOF(char file1[]) {
  strcpy(ERROR_MSG, "premature end of file");
  PrintError();
  ErrorPrintFile(file1, "\0", "\0");
  putc('\n', stderr);
  exit(1);
} //}}}

// TODO remove
// ErrorPrintError_old() //{{{
/*
 * Function to print print error keyword in red
 */
void ErrorPrintError_old() {
  ColourChange(STDERR_FILENO, RED);
  fprintf(stderr, "\nError: ");
  ColourReset(STDERR_FILENO);
} //}}}
// ErrorPrintError() //{{{
/*
 * Function to print print error keyword in red
 */
void ErrorPrintError() {
  ColourChange(STDERR_FILENO, RED);
  fprintf(stderr, "\n  ERROR - %s\n", ERROR_MSG);
  ColourReset(STDERR_FILENO);
} //}}}
// ErrorDiscard()  //{{{
/**
 * Error when number of starting step is higher then the total number of steps
 * in a coordinate file
 */
bool ErrorDiscard(int start, int step, char *file, FILE *coor) {
  int test;
  if ((test = getc(coor)) == EOF) {
    fflush(stdout);
    strcpy(ERROR_MSG, "starting timestep is too high");
    PrintError();
    ErrorPrintFile(file, "\0", "\0");
    fprintf(stderr, "%s, number of timesteps:%s%d%s\n", ErrRed(), ErrYellow(),
                                                        step, ErrColourReset());
    return true;
  } else {
    ungetc(test, coor);
    return false;
  }
} //}}}
// WarnPrintWarning() //{{{
void WarnPrintWarning() {
  fprintf(stderr, "%sWARNING - %s", ErrCyan(), ErrColourReset());
} //}}}
// FilePrintFile() //{{{
void FilePrintFile(char *file, char *colour) {
  ColourChange(STDERR_FILENO, colour);
  fputs("File ", stderr);
  PrintFile(stderr, file, YELLOW);
} //}}}
// WarnElNeutrality() //{{{
/**
 * Function to warn if the system is not electrically neutral.
 */
void WarnElNeutrality(COUNTS Counts, BEADTYPE *BeadType, char *file) {
  double charge = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    // do nothing if at least one bead type had undefined charge
    if (BeadType[i].Charge == CHARGE) {
      return;
    }
    charge += BeadType[i].Charge * BeadType[i].Number;
  }
  if (fabs(charge) > 0.00001) {
    strcpy(ERROR_MSG, "system with net electric charge");
    PrintWarning();
    WarnPrintFile(file, "\0", "\0");
    fprintf(stderr, "%s, %sq = %lf%s\n", ErrCyan(), ErrYellow(), charge,
                                         ErrColourReset());
  }
} //}}}
// ErrorBeadType_old() //{{{
/**
 * Error when non-existent bead is used.
 */
void ErrorBeadType_old(COUNTS Counts, BEADTYPE *BeadType) {
  fprintf(stderr, "       Possible bead names: %s\n", BeadType[0].Name);
  for (int i = 1; i < Counts.TypesOfBeads; i++) {
    fprintf(stderr, "                            %s\n", BeadType[i].Name);
  }
  putc('\n', stderr);
} //}}}
// ErrorMoleculeType_old() //{{{
/**
 * Error when non-existent molecule is used.
 */
void ErrorMoleculeType_old(COUNTS Counts, MOLECULETYPE *MoleculeType) {
  fprintf(stderr, "       Possible molecule names: %s\n", MoleculeType[0].Name);
  for (int i = 1; i < Counts.TypesOfMolecules; i++) {
    fprintf(stderr, "                            %s\n", MoleculeType[i].Name);
  }
  putc('\n', stderr);
} //}}}
