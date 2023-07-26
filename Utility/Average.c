#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Average utility calculates averages for specified column(s) from the input \
file (ignoring empty lines and lines beginning with '#'), printing the results \
to an output file. It has three operation modes based on which of the three \
options is supplied:\n\
1) for -tau option, the utility uses binning method to calculate average, \
statistical error and an estimate of integrated autocorrelation \
time (tau), outputting the specified number of blocks and three values for \
each column used: <simple average> <error> <tau> \
(all on a single line). In this mode, Average appends to the \
output file instead of rewriting it. See the manual for a way to obtain \
a reasonable estimate of tau via rerunning the utility several times.\n\
2) for -b option, the binning method is used to calculate per-block averages.\n\
3) for -m option, the moving method is used to smoothen the input data.");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <output> <column(s)>\n\n", cmd);

  fprintf(ptr, "   <input>        input filename\n");
  fprintf(ptr, "   <output>       output filename\n");
  fprintf(ptr, "   <column(s)>    file column number(s) to analyse\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -tau <int>  estimate tau mode - "
          "number of blocks to split data into\n");
  fprintf(ptr, "      -b <int>    block mode - "
          "number of datapoints per block\n");
  fprintf(ptr, "      -m <int>    moving mode - "
          "number of data points per moving average\n");
  CommonHelp(error, n, opt);
} //}}}

int main ( int argc, char** argv ) {

  // define options //{{{
  int common = 5, all = common + 3, count = 0,
      req_arg = 3; // TODO: will be only two (<input> and <column(s)>)
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-tau");
  strcpy(option[count++], "-b");
  strcpy(option[count++], "-m");
  OptionCheck(argc, argv, req_arg, common, all, option);
  //}}}

  count = 0;

  // <input> - filename of input file
  char input[LINE];
  strcpy(input, argv[++count]);

  char output[LINE] = "";
  strcpy(output, argv[++count]);

  // <column> - column number(s) to analyze
  // TODO: warning if multiple times the same column number
  long int *column = malloc(sizeof *column);
  int col_count = 0;
  while (++count < argc && argv[count][0] != '-') {
    if (!IsNaturalNumber(argv[count], &column[col_count])) {
      ErrorNaN("<column>");
      Help(argv[0], true, common, option);
      exit(1);
    }
    col_count++;
    column = realloc(column, sizeof *column * (col_count + 1));
  }
  int col_max = 0;
  for (int i = 0; i < col_count; i++) {
    if (column[i] > col_max) {
      col_max = column[i];
    }
  }

  bool silent, rubbish2;
  int start = 1, end = -1, skip = 0, rubbish = 0;
  CommonOptions(argc, argv, LINE, &rubbish2, &silent, &rubbish2, &rubbish2,
                &rubbish, &start, &end, &skip);
  start--; // discarded steps rather than starting step //TODO: change

  // -tau option: use block method to get overall average (and stderr and tau)
  int n_blocks = -1;
  if (IntegerOption1(argc, argv, "-tau", &n_blocks)) {
    exit(1);
  }
  // -b option: calculate block averages
  int data_per_block = -1;
  if (IntegerOption1(argc, argv, "-b", &data_per_block)) {
    exit(1);
  }
  // -m option: calculate moving average
  int moving = -1;
  if (IntegerOption1(argc, argv, "-m", &moving)) {
    exit(1);
  }
  if (moving == -1 && n_blocks == -1 && data_per_block == -1) {
    strcpy(ERROR_MSG, "-tau, -b, or -m option must be used");
    PrintError();
    Help(argv[0], true, common, option);
    exit(1);
  }
  if ((moving != -1 && n_blocks != -1) ||
      (moving != -1 && data_per_block != -1) ||
      (n_blocks != -1 && data_per_block != -1)) {
    strcpy(ERROR_MSG, "only one of the -tau, -b, and -m option can be used");
    PrintError();
    Help(argv[0], true, common, option);
  }
  if (moving != -1 && end != -1 && (end - start - moving) < 0) {
    snprintf(ERROR_MSG, LINE, "nothing to compute: %s%d%s-point moving "
             "average from the total of %s%d%s datapoints", ErrYellow(),
             moving, ErrRed(), ErrYellow(), end - start, ErrRed());
    PrintError();
    Help(argv[0], true, common, option);
  }

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  // array to save the data; realloc'd on the fly
  double **data = malloc(sizeof *data * col_count);
  for (int i = 0; i < col_count; i++) {
    data[i] = malloc(sizeof *data[i]);
  }

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
      if (words < col_max) {
        snprintf(ERROR_MSG, LINE, "too few columns (%s%d%s instead "
                 "of %s%d%s); file reading finished", ErrYellow(), words,
                 ErrCyan(), ErrYellow(), col_max, ErrCyan());
        PrintWarningFileLine(input, line_count, split, words);
        break;
      } //}}}
      data_lines++;
      // save the value
      if (start < data_lines) {
        count = data_lines - start - 1;
        for (int i = 0; i < col_count; i++) {
          data[i] = realloc(data[i], sizeof *data[i] * (count + 1));
          data[i][count] = atof(split[column[i]-1]);
        }
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
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  } //}}}

  // -tau mode //{{{
  if (n_blocks > -1) {
    // variables
    // number of data points must be divisible by 'n_blocks'
    int remainder = (data_lines - start) % n_blocks;
    // total number of data points to consider
    count = data_lines - start - remainder;
    // number of data points per block
    int data_per_block = count / n_blocks;
    // block averages
    long double **avg_block = malloc(sizeof **avg_block * col_count);
    // overall averages
    long double **avg_all = malloc(sizeof **avg_all * col_count);
    for (int i = 0; i < col_count; i++) {
      avg_block[i] = calloc(n_blocks, sizeof *avg_block[i]);
      avg_all[i] = calloc(2, sizeof *avg_all[i]); // [0] avg, [1] avg of squares
    }

    int k = remainder; // first datapoint to consider
    for (int i = 0; i < n_blocks; i++) {
      for (int j = 0; j < data_per_block; j++) {
        for (int col = 0; col < col_count; col++) {
          avg_all[col][0] += data[col][k];
          avg_all[col][1] += SQR(data[col][k]);
          avg_block[col][i] += data[col][k];
        }
        k++;
      }
    }
    for (int col = 0; col < col_count; col++) {
      avg_all[col][0] /= count;
      avg_all[col][1] /= count;
      for (int i = 0; i < n_blocks; i++) {
        avg_block[col][i] /= data_per_block;
      }
    }

    double *error = calloc(col_count, sizeof *error),
           *tau_int = calloc(col_count, sizeof *tau_int);
    for (int col = 0; col < col_count; col++) {
      // standard deviation for block averages
      double block_stdev = 0;
      for (int i = 0; i < n_blocks; i++) {
        block_stdev += SQR(avg_block[col][i] - avg_all[col][0]);
      }
      block_stdev /= n_blocks - 1;
      // statistical error
      error[col] = sqrt(block_stdev / n_blocks);
      // approximate integrated autocorrelation time
      tau_int[col] = 0.5 * data_per_block * block_stdev /
                     (avg_all[0][1] - SQR(avg_all[0][0]));
    }

    // print number of blocks, average, statistical error, and estimate of tau
    FILE *fw = OpenFile(output, "a");
    fprintf(fw, " %6d", n_blocks);
    for (int col = 0; col < col_count; col++) {
      fprintf(fw, " %Lf", avg_all[col][0]);
      fprintf(fw, " %lf", error[0]);
      fprintf(fw, " %lf", tau_int[0]);
    }
    putc('\n', fw);
    fclose(fw);

    for (int i = 0; i < col_count; i++) {
      free(avg_block[i]);
      free(avg_all[i]);
    }
    free(avg_block);
    free(avg_all);
    free(error);
    free(tau_int);
  } //}}}

  // -b mode //{{{
  if (data_per_block > -1) {
    // variables
    // number of data points must be divisible by 'data_per_block'
    int remainder = (data_lines - start) % data_per_block;
    // total number of data points to consider
    count = data_lines - start - remainder;
    // number of blocks
    int blocks = count / data_per_block;
    // block averages
    long double **avg_block = malloc(sizeof **avg_block * col_count);
    for (int i = 0; i < col_count; i++) {
      avg_block[i] = calloc(blocks, sizeof *avg_block[i]);
    }

    int k = remainder; // first datapoint to consider
    for (int i = 0; i < blocks; i++) {
      for (int j = 0; j < data_per_block; j++) {
        for (int col = 0; col < col_count; col++) {
          avg_block[col][i] += data[col][k];
        }
        k++;
      }
    }
    for (int col = 0; col < col_count; col++) {
      for (int i = 0; i < blocks; i++) {
        avg_block[col][i] /= data_per_block;
      }
    }
    FILE *fw = OpenFile(output, "w");
    for (int i = 0; i < blocks; i++) {
      for (int col = 0; col < col_count; col++) {
        fprintf(fw, " %Le", avg_block[col][i]);
      }
      putc('\n', fw);
    }
    fclose(fw);

    for (int i = 0; i < col_count; i++) {
      free(avg_block[i]);
    }
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
      for (int col = 0; col < col_count; col++) {
        double tmp = 0;
        for (int j = 0; j < moving; j++) {
          tmp += data[col][i+j];
        }
        tmp /= moving;
        fprintf(fw, " %e", tmp);
      }
      putc('\n', fw);
    }
    fclose(fw);
  } //}}}
  for (int i = 0; i < col_count; i++) {
    free(data[i]);
  }
  free(data);
  free(column);

  return 0;
}
