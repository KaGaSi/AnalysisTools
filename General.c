#include "General.h"
#include "Errors.h"

char line[LINE], *split[SPL_STR];
int words;

static bool ReadLine(FILE *fr, char *line);

// convert string into a number if possible //{{{
/* Functions to test provided string and convert it to a number type. Note that
 * the conversion stops when it encounters an illegal character, so only the
 * beginning of the string must be a legal number of the given type.
 *
 * Example strings: 1) 02.2x & 2) x2.2
 *   IsReal() on 1) gives val=2.2 and returns success (i.e., true)
 *   IsInteger() on 1) gives val=2 and returns success (i.e., true)
 *   On 2), all functions return failure (i.e., false)
 */
bool IsRealNumber(char str[], double *val) {
  char *endptr = NULL;
  *val = strtod(str, &endptr);
  if (endptr == str) {
    return false;
  }
  return true;
}
bool IsPosRealNumber(char str[], double *val) {
  if (IsRealNumber(str, val) && *val > 0) {
    return true;
  } else {
    return false;
  }
}
bool IsIntegerNumber(char str[], long *val) {
  char *endptr = NULL;
  *val = strtol(str, &endptr, 0);
  if (endptr == str) {
    return false;
  }
  return true;
}
bool IsNaturalNumber(char str[], long *val) {
  if (IsIntegerNumber(str, val) && *val > 0) {
    return true;
  } else {
    return false;
  }
}
bool IsWholeNumber(char str[], long *val) {
  if (IsIntegerNumber(str, val) && *val >= 0) {
    return true;
  } else {
    return false;
  }
} //}}}
double VectorLength(VECTOR a) { //{{{
  return sqrt(SQR(a.x) + SQR(a.y) + SQR(a.z));
} //}}}
double Min3(double x, double y, double z) { //{{{
  double min;
  if (x > y) {
    if (y > z) {
      min = z;
    } else {
      min = y;
    }
  } else if (x > z) {
    min = z;
  } else {
    min = x;
  }
  return min;
} //}}}
double Max3(double x, double y, double z) { //{{{
  double max;
  if (x < y) {
    if (y < z) {
      max = z;
    } else {
      max = y;
    }
  } else if (x < z) {
    max = z;
  } else {
    max = x;
  }
  return max;
} //}}}
// swapping functions //{{{
void SwapInt(int *a, int *b) {
  int swap = *a;
  *a = *b;
  *b = swap;
}
void SwapDouble(double *a, double *b) {
  double swap = *a;
  *a = *b;
  *b = swap;
}
void SwapBool(bool *a, bool *b) {
  bool swap = *a;
  *a = *b;
  *b = swap;
} //}}}
// Bubble sort an array; mode = 0: ascendingly, mode = 1: descendingly //{{{
void SortArray(int *array, int length, int mode) {
  if (mode != 0 && mode != 1) {
    strcpy(ERROR_MSG, "SortArray(): use 0 or 1 for sorting mode");
    PrintError();
    exit(1);
  }
  for (int i = 0; i < (length - 1); i++) {
    bool done = true;
    for (int j = 0; j < (length - i - 1); j++) {
      if (mode == 0 && array[j] > array[j + 1]) {
        SwapInt(&array[j], &array[j + 1]);
        done = false;
      }
      if (mode == 1 && array[j] < array[j + 1]) {
        SwapInt(&array[j], &array[j + 1]);
        done = false;
      }
    }
    if (done)
      break;
  }
} //}}}
VECTOR SortVector(VECTOR in) { //{{{
  VECTOR out;
  if (in.x < in.y) {
    if (in.y < in.z) {
      out.x = in.x;
      out.y = in.y;
      out.z = in.z;
    } else if (in.x < in.z) {
      out.x = in.x;
      out.y = in.z;
      out.z = in.y;
    } else {
      out.x = in.z;
      out.y = in.x;
      out.z = in.y;
    }
  } else {
    if (in.x < in.z) {
      out.x = in.y;
      out.y = in.x;
      out.z = in.z;
    } else if (in.y < in.z) {
      out.x = in.y;
      out.y = in.z;
      out.z = in.x;
    } else {
      out.x = in.z;
      out.y = in.y;
      out.z = in.x;
    }
  }
  return out;
} //}}}
bool ReadLine(FILE *fr, char *line) { //{{{
  if (!fgets(line, LINE, fr)) {
    return false; // error/EOF
  }
  // if the line is too long, skip the rest of it
  if (strcspn(line, "\n") == (LINE - 1)) {
    int test;
    do {
      test = getc(fr);
    } while (test != '\n' && test != EOF);
    if (test == EOF) {
      return false;
    }
  }
  return true;
} //}}}
int SplitLine(int max_str, char *out[], char line[], const char delim[]) { //{{{
  // split into words separated by delimiters in delim array
  int words = 0;
  out[words] = strtok(line, delim); // first word
  while (words < max_str && out[words] != NULL) {
    words++; // start from 1, as the first split is already done
    out[words] = strtok(NULL, delim);
  }
  return words;
} //}}}
// ReadAndSplitLine() //{{{
bool ReadAndSplitLine(FILE *fr, int max_strings, const char delim[]) {
  if (!ReadLine(fr, line)) {
    return false;
  }
  words = SplitLine(max_strings, split, line, delim);
  return true;
} //}}}
char * BareCommand(char cmd[]) {
  strcpy(line, cmd);
  int words = SplitLine(SPL_STR, split, line, "/");
  return split[words - 1];
}
void PrintCommand(FILE *ptr, int argc, char *argv[]) { //{{{
  // command may contain whole path - print only the string behind last '/'
  fprintf(ptr, "%s%s", Colour(ptr, WHITE), argv[0]);
  // print the rest of the command
  for (int i = 1; i < argc; i++) {
    fprintf(ptr, " %s", argv[i]);
  }
  fprintf(ptr, "%s\n", Colour(ptr, C_RESET));
} //}}}
// changing the text colour (and making it bold) for cli output //{{{
char *Colour(FILE *f, char colour[]) {
  if (isatty(fileno(f))) {
    return colour;
  } else {
    return "";
  }
}
// colours for stderr
char *ErrRed() { return Colour(stderr, RED); }
char *ErrCyan() { return Colour(stderr, CYAN); }
char *ErrYellow() { return Colour(stderr, YELLOW); }
char *ErrColourReset() { return Colour(stderr, C_RESET); }
// colours for stdout
char *Red() { return Colour(stdout, RED); }
char *Cyan() { return Colour(stdout, CYAN); }
char *Yellow() { return Colour(stdout, YELLOW); }
char *Magenta() { return Colour(stdout, MAGENTA); }
char *Green() { return Colour(stdout, GREEN); }
char *White() { return Colour(stdout, WHITE); }
char *ColourReset() { return Colour(stdout, C_RESET); }
void ColourChange(int a, char *colour) {
  if (isatty(a)) {
    FILE *ptr;
    if (a == STDOUT_FILENO) {
      ptr = stdout;
    } else if (a == STDERR_FILENO) {
      ptr = stderr;
    } else {
      strcpy(ERROR_MSG, "ColourChange() - error that should never happen!");
      PrintError();
      exit(1);
    }
    fputs(colour, ptr);
  }
}
//}}}
FILE *OpenFile(char *file, char *mode) { //{{{
  FILE *ptr = fopen(file, mode);
  if (ptr == NULL) {
    snprintf(ERROR_MSG, LINE, "%sERROR - cannot open file %s%s%s", ErrRed(),
             ErrYellow(), file, ErrRed());
    perror(ERROR_MSG);
    fputs(Colour(stderr, C_RESET), stderr);
    exit(1);
  }
  return ptr;
} //}}}
// initialize arrays to specified value //{{{
void InitDoubleArray(double array[], int n, int val) {
  for (int i = 0; i < n; i++) {
    array[i] = val;
  }
}
void InitIntArray(int array[], int n, int val) {
  for (int i = 0; i < n; i++) {
    array[i] = val;
  }
}
void InitBoolArray(bool array[], int n, bool val) {
  for (int i = 0; i < n; i++) {
    array[i] = val;
  }
}
void InitVecArray(VECTOR array[], int n, bool val) {
  for (int i = 0; i < n; i++) {
    array[i].x = val;
    array[i].y = val;
    array[i].z = val;
  }
}
void InitLong2DArray(long *array[], int m, int n, long val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      array[i][j] = val;
    }
  }
}
void InitDouble2DArray(double *array[], int m, int n, double val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      array[i][j] = val;
    }
  }
} //}}}
// test whether two arrays are the same //{{{
bool SameArray(int arr_1[], int arr_2[], int n) {
  for (int i = 0; i < n; i++) {
    if (arr_1[i] != arr_2[i]) {
      return false;
    }
  }
  return true;
} //}}}
