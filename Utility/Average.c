#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Average uses binning method to analyse data stored in a supplied file. It \
prints average, statistical error and estimate of integrated autocorrelation \
time (tau). Empty lines and lines beginning with '#' are skipped. The program \
prints to standart output 4 values: <n_blocks> <simple average> <statistical \
error> <estimate of tau>.\n\n");
    fprintf(stdout, "\
A way to get a reasonable tau estimate is to use a wide range of <n_blocks> \
and then plot <tau> as a function of <n_blocks>. Since the number of data \
points in a block has to be larger than tau (e.g., ten times larger), \
plotting <number of data lines>/10/<n_blocks> vs. <n_blocks> will produce an \
exponential function that intersects the plotted <tau>. A value of tau near \
the intersection (but to the left where the exponential is above <tau>) can \
be considered a safe estimate for tau.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <column> <discard> <n_blocks>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename\n");
  fprintf(ptr, "   <column>          column number in the file to analyse\n");
  fprintf(ptr, "   <discard>         number of data lines to discard from the file beginning\n");
  fprintf(ptr, "   <n_blocks>        number of blocks for binning\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -h             print this help and exit\n");
  fprintf(ptr, "      --version      print version number and exit\n");
} //}}}

int main ( int argc, char** argv ) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }

  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0 ) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - filename of input file //{{{
  char input[LINE];
  strcpy(input, argv[++count]); //}}}

  // <column> - column number to analyze //{{{
  // Error - non-numeric argument
  if (!IsInteger_old(argv[++count])) {
    ErrorNaN("<column>");
    Help(argv[0], true);
    exit(1);
  }
  int column = atoi(argv[count]); //}}}

  // <discard> - number of lines to discard from the file beginning //{{{
  // Error - non-numeric argument
  if (!IsInteger_old(argv[++count])) {
    ErrorNaN("<discard>");
    Help(argv[0], true);
    exit(1);
  }
  int discard = atoi(argv[count]); //}}}

  // <n_blocks> - number of blocks for binning //{{{
  // Error - non-numeric argument
  if (!IsInteger_old(argv[++count])) {
    ErrorNaN("<n_blocks>");
    Help(argv[0], true);
    exit(1);
  }
  int n_blocks = atoi(argv[count]); //}}}

  double *data = malloc(sizeof *data * 1);

  // read data from <input> file //{{{
  FILE *fr = OpenFile(input, "r");

  int test, lines = 0, all_lines = 0;
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);

    // count data lines, comment, and blank lines in case of error
    all_lines++;

    // get whole line - max 1000 chars //{{{
    char line[LINE];
    fgets(line, sizeof(line), fr); //}}}

    char split[SPL_STR][SPL_LEN], delim[8];
    strcpy(delim, " \t");
    int words = SplitLine_old(split, line, delim);

    // if not empty line or comment continue //{{{
    if (split[0][0] != '#' &&
        split[0][0] != '\n') {
      // error - insufficient number of columns
      // TODO: colours
      if (words < column) {
        ErrorPrintError_old();
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\033[1;33m%s\033[1;31m - too few columns", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      // number of data lines
      lines++;
      // save the value
      if (discard < lines) {
        count = lines - discard - 1;
        data = realloc(data, sizeof *data * (count + 1));
        data[count] = atof(split[column-1]);
      }
    } //}}}
  }
  fclose(fr); //}}}

  // error - <discard> is too large //{{{
  // TODO colours
  if (discard >= lines) {
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\033[1;31m<discard>\033[1;31m - \033[1;33m%d\033[1;31m is too high\n\n", discard);
    fprintf(stderr, "\033[0m");
    Help(argv[0], true);
    exit(1);
  } //}}}

  // variables //{{{
  // number of data points must be dividable by 'n_blocks'
  int remainder = (lines - discard) % n_blocks;
  // total number of data points to consider
  count = lines - discard - remainder;
  // number of data points per block
  int data_per_block = count / n_blocks;
  // allocate array for block averages
  long double *avg_block = malloc(sizeof *avg_block * n_blocks);
  memset(avg_block, 0, sizeof *avg_block * n_blocks);
  // overall averages
  long double avg_all[2] = {0}; // [0] for avg, [1] for avg of squares //}}}

  // calculate averages //{{{
  int k = remainder; // first datapoint to consider
  for (int i = 0; i < n_blocks; i++) {
    for (int j = 0; j < data_per_block; j++) {
      avg_all[0] += data[k];
      avg_all[1] += SQR(data[k]);
      avg_block[i] += data[k++];
    }
  }

  avg_all[0] /= count;
  avg_all[1] /= count;
  for (int i = 0; i < n_blocks; i++) {
    avg_block[i] /= data_per_block;
  } //}}}

  // standard deviation for block averages
  double block_stdev = 0;
  for (int i = 0; i < n_blocks; i++) {
    block_stdev += (avg_block[i] - avg_all[0])*(avg_block[i] - avg_all[0]);
  }
  block_stdev /= n_blocks - 1;
  // statistical error
  double error = sqrt(block_stdev / n_blocks);
  // approximate integrated autocorrelation time
  double tau_int = 0.5 * data_per_block * block_stdev / (avg_all[1] - SQR(avg_all[0]));

  // print number of blocks, average, statistical error, and estimate of tau
  fprintf(stdout, "%4d %Lf %lf %lf\n", n_blocks, avg_all[0], error, tau_int);

  free(data);
  free(avg_block);

  return 0;
}
