#include "Errors.h"
char ERROR_MSG[LINE];

// ErrorArgNumber() //{{{
/**
 * Error when insufficient number of arguments
 */
void ErrorArgNumber(int count, int need) {
  ErrorPrintError_old();
  RedText(STDERR_FILENO);
  fprintf(stderr, "too few mandatory arguments (");
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%d", count);
  RedText(STDERR_FILENO);
  fprintf(stderr, " instead of ");
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%d", need);
  RedText(STDERR_FILENO);
  fprintf(stderr, ")\n\n");
  ResetColour(STDERR_FILENO);
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
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - starting timestep (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", start);
    RedText(STDERR_FILENO);
    fprintf(stderr, ") is higher than the total number of steps (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", step);
    RedText(STDERR_FILENO);
    fprintf(stderr, ")\n\n");
    ResetColour(STDERR_FILENO);
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
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%s", file);
  RedText(STDERR_FILENO);
  fprintf(stderr, " does not have a correct extension (");
  for (int i = 0; i < number; i++) {
    fprintf(stderr, "'%s'", extension[i]);
    if (i < (number-1)) {
      fprintf(stderr, ", ");
    } else {
      fprintf(stderr, ")\n\n");
    }
  }
  ResetColour(STDERR_FILENO);
  return -1;
} //}}}

// ErrorFileOpen() //{{{
/**
 * Error when open file
 */
void ErrorFileOpen(char *file, char mode) {
  ErrorPrintError_old();
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%s", file);
  RedText(STDERR_FILENO);
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
  ResetColour(STDERR_FILENO);
} //}}}

// ErrorNaN() //{{{
/**
 * Error when non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  ErrorPrintError_old();
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%s", option);
  RedText(STDERR_FILENO);
  fprintf(stderr, " - non-numeric argument\n\n");
  ResetColour(STDERR_FILENO);
} //}}}

// ErrorOption() //{{{
/**
 * Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  ErrorPrintError_old();
  RedText(STDERR_FILENO);
  fprintf(stderr, "non-existent option ");
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%s", option);
  RedText(STDERR_FILENO);
  fprintf(stderr, "\n\n");
  ResetColour(STDERR_FILENO);
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
  RedText(STDERR_FILENO);
  if (words == 0) {
    fprintf(stderr, "   Blank line encountered");
  } else {
    fprintf(stderr, "   Wrong line:\n");
    YellowText(STDERR_FILENO);
    for (int i = 0; i < words; i++) {
      if (i != 0) {
        putc(' ', stderr);
      }
      fprintf(stderr, "%s", split[i]);
    }
  }
  fprintf(stderr, "\n\n");
  ResetColour(STDERR_FILENO);
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
    CyanText(STDERR_FILENO);
    fputs(", ", stderr);
    YellowText(STDERR_FILENO);
    fprintf(stderr, "q = %lf\n", charge);
    CyanText(STDERR_FILENO);
    ResetColour(STDERR_FILENO);
  }
} //}}}

// ErrorStartEnd() //{{{
/**
 * Error when ending timestep is higher than the starting one.
 */
void ErrorStartEnd(int start, int end) {
  if (end != -1 && start > end) {
    ErrorPrintError_old();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "-st");
    RedText(STDERR_FILENO);
    fprintf(stderr, " and ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "-e");
    RedText(STDERR_FILENO);
    fprintf(stderr, " - starting step (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", start);
    RedText(STDERR_FILENO);
    fprintf(stderr, ") is higher than ending step (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", end);
    RedText(STDERR_FILENO);
    fprintf(stderr, ")\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
} //}}}

// ErrorPrintError_old() //{{{
/*
 * Function to print print error keyword in red
 */
void ErrorPrintError_old() {
  RedText(STDERR_FILENO);
  fprintf(stderr, "\nError: ");
  ResetColour(STDERR_FILENO);
} //}}}

// ErrorPrintError() //{{{
/*
 * Function to print print error keyword in red
 */
void ErrorPrintError() {
  RedText(STDERR_FILENO);
  fprintf(stderr, "\n  ERROR - %s\n", ERROR_MSG);
  ResetColour(STDERR_FILENO);
} //}}}

// WarnPrintWarning() //{{{
/*
 * Function to print print warning keyword in cyan
 */
void WarnPrintWarning() {
  CyanText(STDERR_FILENO);
  fprintf(stderr, "\n  WARNING - %s\n", ERROR_MSG);
  ResetColour(STDERR_FILENO);
} //}}}

// WarnPrintFile() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void WarnPrintFile(char *file) {
  CyanText(STDERR_FILENO);
  fputs("File ", stderr);
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%s", file);
  ResetColour(STDERR_FILENO);
} //}}}

// PrintFile() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void PrintFile(char *file) {
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%s", file);
  ResetColour(STDERR_FILENO);
} //}}}

// ErrorPrintFileLine() //{{{
/*
 * Function to print file name and the line number in colour.
 */
void PrintFileLine(char *file, int line,
                   char split[SPL_STR][SPL_LEN], int words) {
  RedText(STDERR_FILENO);
  fputs("File ", stderr);
  PrintFile(file);
  RedText(STDERR_FILENO);
  fputs(", line ", stderr);
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%d", line);
  RedText(STDERR_FILENO);
  fputs(":\n", stderr);
  YellowText(STDERR_FILENO);
  for (int i = 0; i < words; i++) {
    if (i != 0) {
      putc(' ', stderr);
    }
    fputs(split[i], stderr);
  }
  fputs("\n\n", stderr);
  ResetColour(STDERR_FILENO);
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
void WarnStopReading(char *vcf_file, int line_count, int step_count) {
  WarnPrintWarning();
  WarnPrintFile(vcf_file);
  CyanText(STDERR_FILENO);
  fputs(" (line ", stderr);
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%d", line_count);
  CyanText(STDERR_FILENO);
  fputs("), step ", stderr);
  YellowText(STDERR_FILENO);
  fprintf(stderr, "%d", step_count);
  CyanText(STDERR_FILENO);
  fputs(" (finished reading)\n", stderr);
  ResetColour(STDERR_FILENO);
} //}}}
