#include "../AnalysisTools.h"

// TODO: add -st/-e/-sk(?)
// TODO: possibility to specify multiple columns

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
        REWRITE!!!\
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
  fprintf(ptr, "   %s <input> <column>\n\n", cmd);

  fprintf(ptr, "   <input>            input filename\n");
  fprintf(ptr, "   <column>           column number in the file to analyse\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -b <int>        block mode - number of blocks "
          "to split data into\n");
  fprintf(ptr, "      -m <file> <int> moving mode - output filename and"
          "number of data points per moving average\n");
  CommonHelp(error, n, opt);
} //}}}

int main ( int argc, char** argv ) {

  // define options //{{{
  int common = 5, all = common + 2, count = 0,
      req_arg = 2; // TODO: will be only two (<input> and <column(s)>)
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-b");
  strcpy(option[count++], "-m");
  OptionCheck(argc, argv, req_arg, common, all, option);
  //}}}

  count = 0;

  // <input> - filename of input file
  char input[LINE];
  strcpy(input, argv[++count]);

  // <column> - column number to analyze
  long int column;
  if (!IsNaturalNumber(argv[++count], &column)) {
    ErrorNaN("<column>");
    Help(argv[0], true, common, option);
    exit(1);
  }

  // // <discard> - number of lines to discard from the file beginning
  // long int start = 0;
  // if (!IsWholeNumber(argv[++count], &start)) {
  //   ErrorNaN("<discard>");
  //   Help(argv[0], true, common, option);
  //   exit(1);
  // }
  // options before reading system data
  bool silent, rubbish2;
  int start = 1, end = -1, skip = 0, rubbish;
  CommonOptions(argc, argv, LINE, &rubbish2, &silent, &rubbish2, &rubbish2,
                &rubbish, &start, &end, &skip);
  start--; // discarded steps rather than starting step //TODO: change


  // -b option: use block method to get overall average (and stderr and tau)
  int n_blocks = -1;
  if (IntegerOption(argc, argv, 1, "-b", &rubbish, &n_blocks)) {
    exit(1);
  }
  // -m option: calculate moving average
  int moving = -1;
  char output[LINE] = "";
  if (FileIntegerOption(argc, argv, 1, "-m", &rubbish, &moving, output)) {
    exit(1);
  }
  if (moving == -1 && n_blocks == -1) {
    strcpy(ERROR_MSG, "either -m or -b option must be used");
    PrintError();
    Help(argv[0], true, common, option);
    exit(1);
  }
  if (moving != -1 && end != -1 && (end - start - moving) < 0) {
    snprintf(ERROR_MSG, LINE, "nothing to compute: %s%d%s-point moving "
             "average from the total of %s%d%s datapoints", ErrYellow(),
             moving, ErrRed(), ErrYellow(), end - start, ErrRed());
    PrintError();
    Help(argv[0], true, common, option);
  }

  // array to save the data; realloc'd on the fly
  double *data = malloc(sizeof *data);

  // read data from <input> file //{{{
  FILE *fr = OpenFile(input, "r");

  int data_lines = 0, line_count = 0;
  while (true) {

    line_count++;

    char line[LINE], *split[SPL_STR];
    int words;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      break;
    }
    // if not empty line or comment continue //{{{
    if (words > 0 && split[0][0] != '#') {
      // error - insufficient number of columns //{{{
      if (words < column) {
        snprintf(ERROR_MSG, LINE, "too few columns (%s%d%s instead of %s%ld%s)",
                 ErrYellow(), words, ErrRed(), ErrYellow(), column, ErrRed());
        PrintErrorFileLine(input, line_count, split, words);
        exit(1);
      } //}}}
      data_lines++;
      // save the value
      if (start < data_lines) {
        count = data_lines - start - 1;
        data = realloc(data, sizeof *data * (count + 1));
        data[count] = atof(split[column-1]);
      }
    } //}}}
    if (end == data_lines) {
      break;
    }
  }
  fclose(fr); //}}}

  // error - <discard> is too large //{{{
  if (start >= data_lines) {
    if (snprintf(ERROR_MSG, LINE, "number of lines to discard (%s%d%s) is "
                 "greater than the number of lines in %s%s%s", ErrYellow(),
                 start, ErrRed(), ErrYellow(), input, ErrRed()) < 0) {
      strcpy(ERROR_MSG, "something wrong with snprintf()");
      exit(1);
    }
    PrintError();
    exit(1);
  } //}}}

  // -b mode //{{{
  if (n_blocks > -1) {
    // variables
    // number of data points must be divisible by 'n_blocks'
    int remainder = (data_lines - start) % n_blocks;
    // total number of data points to consider
    count = data_lines - start - remainder;
    // number of data points per block
    int data_per_block = count / n_blocks;
    // allocate array for block averages
    long double *avg_block = malloc(sizeof *avg_block * n_blocks);
    memset(avg_block, 0, sizeof *avg_block * n_blocks);
    // overall averages
    long double avg_all[2] = {0}; // [0] for avg, [1] for avg of squares

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
    }

    // standard deviation for block averages
    double block_stdev = 0;
    for (int i = 0; i < n_blocks; i++) {
      block_stdev += (avg_block[i] - avg_all[0])*(avg_block[i] - avg_all[0]);
    }
    block_stdev /= n_blocks - 1;
    // statistical error
    double error = sqrt(block_stdev / n_blocks);
    // approximate integrated autocorrelation time
    double tau_int = 0.5 * data_per_block * block_stdev /
                     (avg_all[1] - SQR(avg_all[0]));

    // print number of blocks, average, statistical error, and estimate of tau
    fprintf(stdout, "%4d %Lf %lf %lf\n", n_blocks, avg_all[0], error, tau_int);

    free(avg_block);
  } //}}}

  // -m mode //{{{
  if (moving > -1) {
    // number of input datapoints
    count = data_lines - start;
    // number of output datapoints
    int count_out = count - (moving - 1);
    FILE *fw = OpenFile(output, "w");
    PrintByline(fw, argc, argv);
    for (int i = 0; i < count_out; i++) {
      double tmp = 0;
      for (int j = 0; j < moving; j++) {
        tmp += data[i+j];
      }
      tmp /= moving;
      double step = start + i + 1 + 0.5 * (moving - 1);
      fprintf(fw, "%lf %.3e\n", step, tmp);
    }
    fclose(fw);
  } //}}}
  free(data);

  return 0;
}
