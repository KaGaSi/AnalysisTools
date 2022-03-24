#include "Errors.h"
char ERROR_MSG[LINE];

// ErrorArgNumber() //{{{
/**
 * Error when insufficient number of arguments
 */
void ErrorArgNumber(int count, int need) {
  ErrorPrintError_old();
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, "too few mandatory arguments (");
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", count);
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, " instead of ");
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", need);
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, ")\n\n");
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
    ErrorPrintError_old();
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", file);
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, " - starting timestep (");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", start);
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, ") is higher than the total number of steps (");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", step);
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, ")\n\n");
    ColourReset(STDERR_FILENO);
    return true;
  } else {
    ungetc(test, coor);
    return false;
  }
} //}}}

// ErrorExtension() //{{{
/**
 * Error when missing or incorrect file extension
 */
int ErrorExtension(char *file, int number, char extension[][5]) {
  char *dot = strrchr(file, '.');
  for (int i = 0; i < number; i++) {
    if (dot && strcmp(dot, extension[i]) == 0) {
      return i;
    }
  }
  ErrorPrintError_old();
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%s", file);
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, " does not have a correct extension (");
  for (int i = 0; i < number; i++) {
    fprintf(stderr, "'%s'", extension[i]);
    if (i < (number-1)) {
      fprintf(stderr, ", ");
    } else {
      fprintf(stderr, ")\n\n");
    }
  }
  ColourReset(STDERR_FILENO);
  return -1;
} //}}}

// ErrorFileOpen() //{{{
/**
 * Error when open file
 */
void ErrorFileOpen(char *file, char mode) {
  ErrorPrintError_old();
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%s", file);
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, " - cannot open for ");
  switch(mode) {
    case 'r':
      fprintf(stderr, "reading\n");
      break;
    case 'w':
      fprintf(stderr, "writing\n");
      break;
    case 'a':
      fprintf(stderr, "appending\n");
      break;
    default :
      fprintf(stderr, "Use only r(ead), w(rite), or a(ppend).\n\n");
  }
  putchar('\n');
  ColourReset(STDERR_FILENO);
} //}}}

// ErrorNaN() //{{{
/**
 * Error when non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  ErrorPrintError_old();
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%s", option);
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, " - non-numeric argument\n\n");
  ColourReset(STDERR_FILENO);
} //}}}

// ErrorOption() //{{{
/**
 * Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  ErrorPrintError_old();
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, "non-existent option ");
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%s", option);
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, "\n\n");
  ColourReset(STDERR_FILENO);
} //}}}

// ErrorBeadType() //{{{
/**
 * Error when non-existent bead is used.
 */
void ErrorBeadType(COUNTS Counts, BEADTYPE *BeadType) {
  fprintf(stderr, "       Possible bead names: %s\n", BeadType[0].Name);
  for (int i = 1; i < Counts.TypesOfBeads; i++) {
    fprintf(stderr, "                            %s\n", BeadType[i].Name);
  }
  putc('\n', stderr);
} //}}}

// ErrorMoleculeType() //{{{
/**
 * Error when non-existent molecule is used.
 */
void ErrorMoleculeType(COUNTS Counts, MOLECULETYPE *MoleculeType) {
  fprintf(stderr, "       Possible molecule names: %s\n", MoleculeType[0].Name);
  for (int i = 1; i < Counts.TypesOfMolecules; i++) {
    fprintf(stderr, "                            %s\n", MoleculeType[i].Name);
  }
  putc('\n', stderr);
} //}}}

// ErrorPrintLine() //{{{
/**
 * Print provided strings (array of strings generally created using
 * SplitLine()) to error output.
 */
void ErrorPrintLine(char split[SPL_STR][SPL_LEN], int words) {
  ColourText(STDERR_FILENO, RED);
  if (words == 0) {
    fprintf(stderr, "   Blank line encountered");
  } else {
    fprintf(stderr, "   Wrong line:\n");
    ColourText(STDERR_FILENO, YELLOW);
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
    WarnPrintWarning();
    WarnPrintFile(file);
    ColourText(STDERR_FILENO, CYAN);
    fputs(", ", stderr);
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "q = %lf\n", charge);
    ColourReset(STDERR_FILENO);
  }
} //}}}

// ErrorStartEnd() //{{{
/**
 * Error when ending timestep is higher than the starting one.
 */
void ErrorStartEnd(int start, int end) {
  if (end != -1 && start > end) {
    ErrorPrintError_old();
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-st");
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, " and ");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-e");
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, " - starting step (");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", start);
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, ") is higher than ending step (");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d", end);
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, ")\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
} //}}}

// ErrorPrintError_old() //{{{
/*
 * Function to print print error keyword in red
 */
void ErrorPrintError_old() {
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, "\nError: ");
  ColourReset(STDERR_FILENO);
} //}}}

// ErrorPrintError() //{{{
/*
 * Function to print print error keyword in red
 */
void ErrorPrintError() {
  ColourText(STDERR_FILENO, RED);
  fprintf(stderr, "\n  ERROR - %s\n", ERROR_MSG);
  ColourReset(STDERR_FILENO);
} //}}}

// WarnPrintWarning() //{{{
/*
 * Function to print print warning keyword in cyan
 */
void WarnPrintWarning() {
  ColourText(STDERR_FILENO, CYAN);
  fprintf(stderr, "\n  WARNING - %s\n", ERROR_MSG);
  ColourReset(STDERR_FILENO);
} //}}}

// WarnPrintFile() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void WarnPrintFile(char *file) {
  ColourText(STDERR_FILENO, CYAN);
  fputs("File ", stderr);
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%s", file);
  ColourReset(STDERR_FILENO);
} //}}}

// PrintFile() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void PrintFile(char *file, int colour) {
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%s", file);
  ColourReset(STDERR_FILENO);
} //}}}

// FilePrintFile() //{{{
void FilePrintFile(char *file, int colour) {
  ColourText(STDERR_FILENO, colour);
  fputs("File ", stderr);
  PrintFile(file, YELLOW);
} //}}}

// ErrorPrintFileLine() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void PrintFileLine(char *file, int line,
                   char split[SPL_STR][SPL_LEN], int words) {
  ColourText(STDERR_FILENO, RED);
  fputs("File ", stderr);
  PrintFile(file, YELLOW);
  ColourText(STDERR_FILENO, RED);
  fputs(", line ", stderr);
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", line);
  ColourText(STDERR_FILENO, RED);
  fputs(":\n", stderr);
  ColourText(STDERR_FILENO, YELLOW);
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

// WarnStopReading() //{{{
/*
 * Warning when stopping file reading due to some error in the file
 */
void WarnStopReading(char *vcf_file, int line_count,
                     char split[SPL_STR][SPL_LEN], int words) {
  WarnPrintWarning();
  ColourText(STDERR_FILENO, CYAN);
  fputs("Stopped reading file ", stderr);
  PrintFile(vcf_file, YELLOW);
  ColourText(STDERR_FILENO, CYAN);
  fputs(" before line ", stderr);
  ColourText(STDERR_FILENO, YELLOW);
  fprintf(stderr, "%d", line_count);
  ColourText(STDERR_FILENO, CYAN);
  fputs(":\n", stderr);
  PrintLine(split, words, YELLOW, CYAN);
} //}}}

// PrintLine() //{{{
void PrintLine(char split[SPL_STR][SPL_LEN], int words,
               int col_line, int col_blank) {
  if (words == 0) {
    ColourText(STDERR_FILENO, col_blank);
    fprintf(stderr, "Blank line");
    ColourReset(STDERR_FILENO);
  } else {
    ColourText(STDERR_FILENO, col_line);
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
