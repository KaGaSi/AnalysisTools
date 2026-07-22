#include "General.h"
#include "Arrays.h"
#include "Errors.h"
#include <errno.h>

// TODO: check the digits stuff
//       also, some needs deleting when it's implemented everywhere

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
  char *endptr = nullptr;
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
  char *endptr = nullptr;
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
// read a line from a file //{{{
bool ReadLine(FILE *fr, char *line) {
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
// split a string by specified delimiters //{{{
int SplitLine(const int max_str, char **out, char *line, const char delim[]) {
  // split into words separated by delimiters in delim array
  int words = 0;
  out[words] = strtok(line, delim); // first word
  while (words < (max_str - 1) && out[words] != nullptr) {
    words++; // start from 1, as the first split is already done
    out[words] = strtok(nullptr, delim);
  }
  return words;
} //}}}
// read a line from file and split it into individual strings //{{{
bool ReadAndSplitLine(FILE *fr, const int max_str, const char delim[]) {
  if (!ReadLine(fr, line)) {
    return false;
  }
  words = SplitLine(max_str, split, line, delim);
  return true;
} //}}}
// increment linecount, read and split next line //{{{
bool CountLineReadLine(int *line_count, FILE *fr,
                       const char *file, const char *msg) {
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (msg[0] != '\0') {
      ErrorEOF(file, (char *)msg);
      exit(1);
    } else {
      return false;
    }
  }
  return true;
} //}}}
// write a line into a file //{{{
void WriteSplitLine(FILE *f) {
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
// TODO: OK, this is ugly: 1) exit on cosmetics - really? 2) errno clobbering
void ColourChange(const int a, const char *colour) {
  int saved_errno = errno;
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
    errno = saved_errno;
  }
} //}}}
FILE * OpenFile(const char *file, char *mode) { //{{{
  FILE *ptr = fopen(file, mode);
  if (ptr == nullptr) {
    snprintf(ERROR_MSG, LINE, "%sERROR - cannot open file %s%s%s",
             ErrRed(), ErrYellow(), file, ErrRed());
    perror(ERROR_MSG);
    fputs(Colour(stderr, C_RESET), stderr);
    exit(1);
  }
  return ptr;
} //}}}
// initialize arrays to specified value //{{{
void InitDoubleArray(double *arr, const int n, const double val) {
  for (int i = 0; i < n; i++) {
    arr[i] = val;
  }
}
void InitIntArray(int *arr, const int n, const int val) {
  for (int i = 0; i < n; i++) {
    arr[i] = val;
  }
}
void InitBoolArray(bool *arr, const int n, const bool val) {
  for (int i = 0; i < n; i++) {
    arr[i] = val;
  }
}
void InitLong2DArray(long **arr, const int m, const int n, const long val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      arr[i][j] = val;
    }
  }
}
void InitDouble2DArray(double *arr, const int m, const int n,
                       const double val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      arr[i*n+j] = val; // Access using row-major order
    }
  }
}
void InitInt2DArray(int *arr, const int m, const int n, const int val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      arr[i*n+j] = val; // Access using row-major order
    }
  }
}
void InitBool2DArray(bool **arr, const int m, const int n, const bool val) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      arr[i][j] = val;
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
  if (dest == nullptr || src == nullptr || dest_size == 0) {
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
    fprintf(stderr, "s_realloc: new size is 0\n");
    exit(1);
  }
  void *temp = realloc(ptr, new_size);
  if (temp == nullptr) {
    fprintf(stderr, "s_realloc: NULL returned\n");
    exit(1);
  }
  return temp;
} //}}}
// stuff to count digits and print correctly column width //{{{
void CountDigits(const double num, int digits[2]) {
  int max_precision = 6;
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
void FprintfRow(FILE *fw, int columns,
                const double value[columns], const int digits[columns][2]) {
  for (int col = 0; col < columns; col++) {
    Fprintf1(fw, value[col], digits[col]);
  }
  putc('\n', fw);
}
//}}}

// Count meaningful decimal digits of a number up to max_precision
static int CountDecimalDigits(double x, int max_precision) {
  // sanity checks
  if (max_precision < 1) {
    return 0;
  }
  if (max_precision > 17) {
    max_precision = 17; // maximum precision of doubles
  }
  // create format string, e.g., '%.6f' for max_precision == 6
  char fmt[16];
  snprintf(fmt, sizeof(fmt), "%%.%df", max_precision);
  fmt[15] = '\0';
  // convert number to a string with given maximum number of decimals
  char buf[64];
  snprintf(buf, sizeof(buf), fmt, x);
  buf[63] = '\0';
  // find decimal point
  char *dot = strchr(buf, '.');
  if (!dot) { // no decimal point -> no digits to count
    return 0;
  }
  // walk backwards from the end of the string until non-zero character
  char *end = buf + strlen(buf);
  char *p = end - 1;
  while (p > dot && *p == '0') {
    p--;
  }
  // discard the dot if no non-zero characters behind it were found
  if (*p == '.') {
    p--;
  }
  // return digit count: p-dot is number of digits between '.' and last non-zero
  if ((p - dot) > 0) {
    return p - dot;
  } else {
    return 0;
  }
}
// Determine per-column precision and width from actual data
void ComputeColumnWidths(const int nrows, const int ncols, ArrNDd *data,
                         int max_precision) {
  for (int i = 0; i < ncols; i++) {
    size_t index[] = {nrows, i};
    SetArrND(data, index, 0);
    index[0] = nrows + 1;
    SetArrND(data, index, 0);
  } // Pass 1: find max precision actually needed
  for (int i = 0; i < ncols; i++) {
    int col_prec = 0;
    for (int j = 0; j < nrows; j++) {
      size_t index[] = {j, i};
      int p = CountDecimalDigits(GetArrND(data, index), max_precision);
      if (p > col_prec) {
        col_prec = p;
      }
    }
    size_t index[] = {nrows + 1, i};
    SetArrND(data, index, col_prec);
  }
  // Pass 2: compute column widths using chosen per-column precision
  for (int i = 0; i < ncols; i++) {
    for (int j = 0; j < nrows; j++) {
      char buf[64];
      // snprintf(buf, sizeof buf, "%.*f", precisions[i][1], data[j][i]);
      size_t index[] = {j, i};
      size_t index2[] = {nrows + 1, i};
      snprintf(buf, sizeof buf, "%.*f", (int)GetArrND(data, index2),
                                        GetArrND(data, index));
      int w = (int)strlen(buf);
      index2[0] = nrows;
      if (w > GetArrND(data, index2)) {
        SetArrND(data, index2, w);
      }
    }
  }
}
// Print one formatted value
void PrintDataValue(FILE *fw, const int nrows, const int row, const int col,
                    const ArrNDd *data) {
    size_t index[] = {row, col};
    size_t index0[] = {nrows, col};
    size_t index1[] = {nrows + 1, col};
    fprintf(fw, " %*.*f", (int)GetArrND(data, index0),
                          (int)GetArrND(data, index1),
                          GetArrND(data, index));
}
// Print one formatted row
void PrintDataRow(FILE *fw, const int nrows, const int row, const int ncols,
                  const ArrNDd *data) {
  for (int i = 0; i < ncols; i++) {
    PrintDataValue(fw, nrows, row, i, data);
  }
  putc('\n', fw);
}
// Print all formatted rows
void PrintDataAll(FILE *fw, const int nrows, const int ncols,
                  const ArrNDd *data) {
  for (int j = 0; j < nrows; j++) {
    PrintDataRow(fw, nrows, j, ncols, data);
  }
}

void SkipLine(FILE *fr) {
  int ch;
  while ((ch = getc(fr)) != '\n' && ch != EOF)
    ;
}
