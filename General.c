#include "General.h"
#include "Errors.h"

// convert string into number if possible //{{{
/* Functions to test provided string and convert it to a number type. Note that
 * the conversion stops when it encounters an illegal character, so only the
 * beginning of the string must be a legal number of the given type.
 *
 * Example strings: 1) 02.2x & 2) x2.2
 *   IsReal() on 1) gives val=2.2 and returns success (i.e., true)
 *   IsInteger() on 1) gives val=2 and returns success (i.e., true)
 *   On 2), all functions return failure (i.e., false)
 */
bool IsReal(char *str, double *val) {
  char *endptr = NULL;
  *val = strtod(str, &endptr);
  if (endptr == str) {
    return false;
  }
  return true;
}
bool IsPosReal(char *str, double *val) {
  if (IsReal(str, val) && *val > 0) {
    return true;
  } else {
    return false;
  }
}
bool IsInteger(char *str, long *val) {
  char *endptr = NULL;
  *val = strtol(str, &endptr, 0);
  if (endptr == str) {
    return false;
  }
  return true;
}
bool IsPosInteger(char *str, long *val) {
  if (IsInteger(str, val) && *val > 0) {
    return true;
  } else {
    return false;
  }
}
bool IsNatural(char *str, long *val) {
  if (IsInteger(str, val) && *val >= 0) {
    return true;
  } else {
    return false;
  }
}
 //}}}

// Length() //{{{
/*
 * Function to calculate vector length.
 */
double Length(VECTOR a) {
  double length = sqrt(SQR(a.x) + SQR(a.y) + SQR(a.z));
  return length;
} //}}}

// Min3() //{{{
/*
 * Function returning the lowest number from three floats.
 */
double Min3(double x, double y, double z) {

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

// Max3() //{{{
/*
 * Function returning the highest number from three floats.
 */
double Max3(double x, double y, double z) {

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

// SortVector() //{{{
/*
 * Function returning sorted numbers x < y < z.
 */
VECTOR SortVector(VECTOR in) {
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
}
 //}}}

// SortArray() //{{{
// Bubble sort an array; mode = 0: sort ascendingly, mode = 1: sort
// descendingly.
void SortArray(int *array, int length, int mode) {
  if (mode != 0 && mode != 1) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError - SortArray(): use 0 or 1 for sorting mode\n");
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  for (int i = 0 ; i < (length-1); i++) {
    bool done = true;
    for (int j = 0 ; j < (length-i-1); j++) {
      if (mode == 0 && array[j] > array[j+1]) {
        SwapInt(&array[j], &array[j+1]);
        done = false;
      }
      if (mode == 1 && array[j] < array[j+1]) {
        SwapInt(&array[j], &array[j+1]);
        done = false;
      }
    }
    if (done)
      break;
  }
} //}}}

// ReadLine() //{{{
bool ReadLine(FILE *fr, int max_char, char *line) {
  if (!fgets(line, LINE, fr)) {
    return false; // error/EOF
  }
  // if the line is too long, skip the rest of it
  if (strcspn(line, "\n") == (LINE-1)) {
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

// SplitLine() //{{{
/*
 * Function that splits the provided line into individual strings.
 */
int SplitLine(int max_str, char *out[], char *line, const char *delim) {
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
bool ReadAndSplitLine(FILE *fr, int max_char, char *line, int *words, char *out[],
                      int max_strings, const char *delim) {
  if (!ReadLine(fr, max_char, line)) {
    return false;
  }
  *words = SplitLine(max_strings, out, line, delim);
  return true;
} //}}}

// PrintCommand() //{{{
/*
 * Function to print the used command (at most 30 words).
 */
void PrintCommand(FILE *ptr, int argc, char *argv[]) {
  // command may contain whole path - print only the string behind last '/'
  char *split[SPL_STR], str[LINE];
  strcpy(str, argv[0]);
  int words = SplitLine(SPL_STR, split, str, "/");
  fprintf(ptr, " %s%s", Colour(ptr, WHITE), split[words-1]);
  // print the rest of the command
  for (int i = 1; i < argc; i++)
    fprintf(ptr, " %s", argv[i]);
  fprintf(ptr, "%s\n", Colour(ptr, C_RESET));
} //}}}

// changing the text colour (and making it bold) for cli output //{{{
char *Colour(FILE *f, char *colour) {
  if (isatty(fileno(f))) {
    return colour;
  } else {
    return "";
  }
}
void PrintColour(FILE *f, char *colour) {
  if (isatty(fileno(f))) {
    fputs(colour, f);
  }
}
// colours for stderr
char *ErrRed() {
  return Colour(stderr, RED);
}
char *ErrCyan() {
  return Colour(stderr, CYAN);
}
char *ErrYellow() {
  return Colour(stderr, YELLOW);
}
char *ErrColourReset() {
  return Colour(stderr, C_RESET);
}
// colours for stdout
char *Red() {
  return Colour(stdout, RED);
}
char *Cyan() {
  return Colour(stdout, CYAN);
}
char *Yellow() {
  return Colour(stdout, YELLOW);
}
char *Magenta() {
  return Colour(stdout, MAGENTA);
}
char *Green() {
  return Colour(stdout, GREEN);
}
char *White() {
  return Colour(stdout, WHITE);
}
char *ColourReset() {
  return Colour(stdout, C_RESET);
}
void ColourChange(int a, char *colour) {
  if (isatty(a)) {
    FILE *ptr;
    if (a == STDOUT_FILENO) {
      ptr = stdout;
    } else if (a == STDERR_FILENO) {
      ptr = stderr;
    }
    fputs(colour, ptr);
  }
}
 //}}}

// OpenFile() //{{{
FILE *OpenFile(char *file, char *mode) {
  FILE *ptr = fopen(file, mode);
  if (ptr == NULL) {
    strcpy(ERROR_MSG, "cannot open file");
    PrintError();
    ErrorPrintFile(file, "\0");
    fputs(Colour(stderr, RED), stderr);
    perror(" ");
    fputs(Colour(stderr, C_RESET), stderr);
    exit(1);
  }
  return ptr;
} //}}}

// TODO: remove
// IsReal_old() //{{{
/*
 * Function to test if provided string is a real number.
 */
bool IsReal_old(char *a) {
  // only one dot and scientific e can be present
  bool dot = false,
       sci_e = false;
  // wrong first character - can be minus, dot, or number
  if (a[0] != '-' && a[0] != '.' && (a[0] < '0' || a[0] > '9')) {
    return false;
  } else if (a[0] == '.') {
    dot = true;
  }
  // test the remaining characters - either digit, or dot (but only 1 in total)
  for (int i = 1; i < strlen(a); i++) {
    if (a[i] == '.') {
      if (dot) { // has there been a dot already?
        return false;
      } else {
        dot = true;
      }
    } else if (a[i] == 'e' || a[i] == 'E') { // scientific notation?
      // format must be: e/E[-/+]<int>
      if (sci_e) {
        return false;
      } else {
        sci_e = true;
      }
      if ((i+1) >= strlen(a)) {
        return false;
      } else if (a[i+1] == '-' || a[i+1] == '+') {
        i++; // skip the '-' sign in the for loop
        if ((i+1) >= strlen(a)) {
          return false;
        }
      }
    } else if (a[i] < '0' || a[i] > '9') {
      return false;
    }
  }
  return true;
} //}}}
// IsInteger_old() //{{{
/*
 * Function to test if provided string is a non-negative whole number.
 */
bool IsInteger_old(char *a) {
  // test the remaining characters - either digit, or dot (but only 1 in total)
  for (int i = 0; i < strlen(a); i++) {
    if (a[i] < '0' || a[i] > '9') {
      return false;
    }
  }
  return true;
} //}}}
// IsPosReal_old() //{{{
/*
 * Function to test if provided string is a positive real number.
 */
bool IsPosReal_old(char *a) {
  if (IsReal_old(a) && atof(a) > 0) {
    return true;
  } else {
    return false;
  }
} //}}}
// IsNatural_old()  //{{{
bool IsNatural_old(char *a) {
  if (IsInteger_old(a) && atof(a) >= 0) {
    return true;
  } else {
    return false;
  }
} //}}}
// SplitLine_old() //{{{
/*
 * Function that splits the provided line into individual strings.
 */
int SplitLine_old(char out[SPL_STR][SPL_LEN], char *line, const char *delim) {
  // trim whitespaces at the beginning and end of line
  strcpy(line, TrimLine(line));
  // split into words separated by delimiters in delim array
  char *split[SPL_STR];
  int words = 0;
  split[words] = strtok(line, delim); // first word
  while (words < (SPL_STR-1) && split[words] != NULL) {
    words++; // start from 1, as the first split is already done
    split[words] = strtok(NULL, delim);
  }
  if (words == 0) {
    out[0][0] = '\0';
  } else {
    // copy splits into the output array
    for (int i = 0; i < words; i++) {
      snprintf(out[i], SPL_LEN, "%s", split[i]);
    }
  }
  return words;
} //}}}
// TrimLine() //{{{
/*
 * Function to trim whitespace from the beginning and end of a string.
 */
char * TrimLine(char *line) {
  int length = strlen(line);
  static char trimmed[LINE];
  strcpy(trimmed, line);
  // 1) trailing whitespace
  while (length > 1 &&
         (trimmed[length-1] == ' ' ||
          trimmed[length-1] == '\n' ||
          trimmed[length-1] == '\r' ||
          trimmed[length-1] == '\t')) {
    trimmed[length-1] = '\0';
    length = strlen(trimmed);
  }
  // 2) preceding whitespace
  while (length >= 1 &&
         (trimmed[0] == ' ' ||
          trimmed[0] == '\n' ||
          trimmed[0] == '\r' ||
          trimmed[0] == '\t')) {
    for (int i = 0; i < length; i++) { // line[length] contains '\0'
      trimmed[i] = trimmed[i+1];
    }
    length = strlen(trimmed);
  }
  return trimmed;
} //}}}
