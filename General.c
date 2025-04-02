#include "General.h"
#include "Errors.h"

static void CountDigits(const double num, int digits[2]);
static void MaxDigits(const double num, int max_digits[2]);

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
bool IsRealNumber(const char *str, double *val) {
  char *endptr = NULL;
  *val = strtod(str, &endptr);
  if (endptr == str) {
    return false;
  }
  return true;
}
bool IsPosRealNumber(const char *str, double *val) {
  if (IsRealNumber(str, val) && *val > 0) {
    return true;
  } else {
    return false;
  }
}
bool IsIntegerNumber(const char *str, long *val) {
  char *endptr = NULL;
  *val = strtol(str, &endptr, 0);
  if (endptr == str) {
    return false;
  }
  return true;
}
bool IsNaturalNumber(const char *str, long *val) {
  if (IsIntegerNumber(str, val) && *val > 0) {
    return true;
  } else {
    return false;
  }
}
bool IsWholeNumber(const char *str, long *val) {
  if (IsIntegerNumber(str, val) && *val >= 0) {
    return true;
  } else {
    return false;
  }
} //}}}
// bubble sort int/double array ascendingly/descendingly //{{{
void SortArray(void *array, const int length, const int mode, const char type) {
  if (mode != 0 && mode != 1) {
    err_msg("SortArray*(): use 0 or 1 for sorting mode");
    PrintError();
    exit(1);
  }
  if (type != 'i' && type != 'd') {
    err_msg("SortArray(): use 'i' or 'd' for integer or double array");
    PrintError();
    exit(1);
  }
  for (int i = 0; i < (length - 1); i++) {
    bool done = true;
    for (int j = 0; j < (length - i - 1); j++) {
      if (type == 'i') {
        int *arr = (int *)array;
        int *a = &arr[j];
        int *b = &arr[j];
        if (mode == 0) {
          b = &arr[j+1];
        } else {
          a = &arr[j+1];
        }
        if (*a > *b) {
          SwapInt(a, b);
          done = false;
        }
      } else {
        double *arr = (double *)array;
        double *a = &arr[j];
        double *b = &arr[j];
        if (mode == 0) {
          b = &arr[j+1];
        } else {
          a = &arr[j+1];
        }
        if (*a > *b) {
          SwapDouble(a, b);
        }
      }
    }
    if (done)
      break;
  }
} //}}}
bool ReadLine(FILE *fr, char *line) { //{{{
  if (!fgets(line, LINE, fr)) {
    return false; // error/EOF
  }
  // if the line is too long, skip the rest of it
  size_t len = strcspn(line, "\n");
  if (len == (LINE - 1) && line[len] != '\n') {
    int ch;
    while ((ch = getc(fr)) != '\n' && ch != EOF)
      ;
  }
  return true;
} //}}}
// SplitLine() //{{{
int SplitLine(const int max_str, char **out, char *line, const char delim[]) {
  // split into words separated by delimiters in delim array
  int words = 0;
  out[words] = strtok(line, delim); // first word
  while (words < max_str && out[words] != NULL) {
    words++; // start from 1, as the first split is already done
    out[words] = strtok(NULL, delim);
  }
  return words;
} //}}}
bool ReadAndSplitLine(FILE *fr, const int max_str, const char delim[]) { //{{{
  if (!ReadLine(fr, line)) {
    return false;
  }
  words = SplitLine(max_str, split, line, delim);
  return true;
} //}}}
void WriteSplitLine(FILE *f) { //{{{
  for (int i = 0; i < words; i++) {
    fprintf(f, " %s", split[i]);
  }
  putc('\n', f);
} //}}}
const char * StripPath(const char cmd[]) { //{{{
  // Find the last occurrence of '/' in the path
  const char *command = strrchr(cmd, '/');
  if (command) { // '/' found
    return command + 1;
  } else { // '/' not found
    return cmd;
  }
} //}}}
void PrintCommand(FILE *ptr, const int argc, char **argv) { //{{{
  fprintf(ptr, "%s%s", Colour(ptr, WHITE), StripPath(argv[0]));
  // print the rest of the command
  for (int i = 1; i < argc; i++) {
    fprintf(ptr, " %s", argv[i]);
  }
  fprintf(ptr, "%s\n", Colour(ptr, C_RESET));
} //}}}
// changing the text colour (and making it bold) for cli output //{{{
void ColourChange(const int a, const char *colour) {
  if (isatty(a)) {
    FILE *ptr;
    if (a == STDOUT_FILENO) {
      ptr = stdout;
    } else if (a == STDERR_FILENO) {
      ptr = stderr;
    } else {
      err_msg("ColourChange() - error that should never happen!");
      PrintError();
      exit(1);
    }
    fputs(colour, ptr);
  }
} //}}}
FILE * OpenFile(const char *file, char *mode) { //{{{
  FILE *ptr = fopen(file, mode);
  if (ptr == NULL) {
    snprintf(ERROR_MSG, LINE, "%sERROR - cannot open file %s%s%s",
             ErrRed(), ErrYellow(), file, ErrRed());
    perror(ERROR_MSG);
    fputs(Colour(stderr, C_RESET), stderr);
    exit(1);
  }
  return ptr;
} //}}}
// initialize arrays to specified value //{{{
void InitDoubleArray(double *array, const int n, const double val) {
  for (int i = 0; i < n; i++) {
    array[i] = val;
  }
}
void InitIntArray(int *array, const int n, const int val) {
  for (int i = 0; i < n; i++) {
    array[i] = val;
  }
}
void InitBoolArray(bool *array, const int n, const bool val) {
  for (int i = 0; i < n; i++) {
    array[i] = val;
  }
}
void InitLong2DArray(long **array, const int m, const int n, const long val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      array[i][j] = val;
    }
  }
}
void InitDouble2DArray(double **array, const int m, const int n,
                       const double val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      array[i][j] = val;
    }
  }
}
void InitInt2DArray(int *array, const int m, const int n, const int val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      array[i * n + j] = val;  // Access using row-major order
    }
  }
}
void InitBool2DArray(bool **array, const int m, const int n, const bool val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      array[i][j] = val;
    }
  }
} //}}}
bool SameArrayInt(const int *arr_1, const int *arr_2, const int n) { //{{{
  for (int i = 0; i < n; i++) {
    if (arr_1[i] != arr_2[i]) {
      return false;
    }
  }
  return true;
} //}}}
void s_strcpy(char *dest, const char *src, const size_t dest_size) { //{{{
  if (dest == NULL || src == NULL || dest_size == 0) {
    fprintf(stderr, "s_strcpy error...");
    exit(1);
  }
  size_t i;
  for (i = 0; i < dest_size - 1 && src[i] != '\0'; i++) {
    dest[i] = src[i];
  }
  dest[i] = '\0';
} //}}}
void* s_realloc(void *ptr, size_t new_size) { //{{{
  if (new_size == 0) {
    fprintf(stderr, "realloc error\n");
    exit(1);
  }
  void *temp = realloc(ptr, new_size);
  if (temp == NULL) {
    fprintf(stderr, "realloc error\n");
    exit(1);
  }
  return temp;
} //}}}
// stuff to count digits and print correctly column width //{{{
static void CountDigits(const double num, int digits[2]) {
  int max_precision = 5;
  double frac_part, int_part;
  // Handle negative numbers
  int neg = 0;
  if (num < 0) {
    neg = 1;
  }
  double temp = fabs(num);
  // Split into integer and fractional parts
  frac_part = modf(temp, &int_part);
  // Count integer digits
  if (int_part == 0) {
    digits[0] = 1 + neg;
  } else {
    digits[0] = (int)log10(int_part) + 1 + neg;
  }
  // Count fractional digits dynamically
  digits[1] = 0;
  while (digits[1] < max_precision) {
    frac_part *= 10;
    frac_part = round(frac_part * 1e9) / 1e9; // roundto avoid precision errors
    if (frac_part < 1e-9) {
      break; // stop when there's no more fraction left
    }
    frac_part -= floor(frac_part);
    digits[1]++;
  }
  // +1 for '.'
  if (digits[1] > 0) {
    digits[0]++;
  }
}
static void MaxDigits(const double num, int max_digits[2]) {
  int digits[2];
  CountDigits(num, digits);
  for (int i = 0; i < 2; i++) {
    if (digits[i] > max_digits[i]) {
      max_digits[i] = digits[i];
    }
  }
}
void FillMaxDigits(const int columns, const int n,
                   double *data[n], int (*digits)[2]) {
  for (int i = 0; i < n; i++) {
    for (int col = 0; col < columns; col++) {
      MaxDigits(data[i][col], digits[col]);
    }
  }
}
void Fprintf1(FILE *f, const double value, const int digits[2]) {
  fprintf(f, " %*.*f", digits[0] + digits[1], digits[1], value);
}
//}}}
