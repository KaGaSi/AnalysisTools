#include "../AnalysisTools.h"

// TODO: remove <output> in favour of -tau/-b/-m <output> <int> options

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
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
3) for -m option, the moving method is used to smoothen the input data.\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <output> <column(s)>\n\n", cmd);

  fprintf(ptr, "<input>             input filename\n");
  fprintf(ptr, "<output>            output filename\n");
  fprintf(ptr, "<column(s)>         column number(s) to analyse\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -tau <int>        estimate tau mode - "
          "number of blocks to split data into\n");
  fprintf(ptr, "  -b <int>          block mode - "
          "number of datapoints per block\n");
  fprintf(ptr, "  -m <int>          moving mode - "
          "number of data points per moving average\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int tau, block, moving; // -tau -b -m
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main ( int argc, char** argv ) {

  // define options & check their validity
  int common = 5, all = common + 3, count = 0,
      req_arg = 3; // TODO: will be only two (<input> and <column(s)>)
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, false, option,
               "-st", "-e", "--silent", "--help",
               "--version", "-tau", "-b", "-m");

  count = 0;
  OPT *opt = opt_create();

  char fin[LINE];
  s_strcpy(fin, argv[++count], LINE);
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  // <column> - column number(s) to analyze
  // TODO: warning if multiple times the same column number
  long int *column = malloc(sizeof *column);
  int col_count = 0;
  while (++count < argc && argv[count][0] != '-') {
    if (!IsNaturalNumber(argv[count], &column[col_count])) {
      ErrorNaN("<column>");
      Help(StripPath(argv[0]), true, common, option);
      exit(1);
    }
    col_count++;
    column = s_realloc(column, sizeof *column * (col_count + 1));
  }
  int col_max = 0;
  for (int i = 0; i < col_count; i++) {
    if (column[i] > col_max) {
      col_max = column[i];
    }
  }

  SYS_FILES trash = InitSysFiles; // unused
  opt->c = CommonOptions(argc, argv, trash);
  opt->c.start--; // discarded steps rather than starting step //TODO: change

  // -tau option: use block method to get overall average (and std_err and tau)
  opt->tau = 0;
  OneNumberOption(argc, argv, "-tau", &opt->tau, 'i');
  // -b option: calculate block averages
  opt->block = 0;
  OneNumberOption(argc, argv, "-b", &opt->block, 'i');
  // -m option: calculate moving average
  opt->moving = 0;
  OneNumberOption(argc, argv, "-m", &opt->moving, 'i');
  // if (opt->moving == 0 && opt->tau == 0 && opt->block == 0) {
  //   err_msg("one of -tau, -b, or -m options must be used");
  //   PrintError();
  //   Help(StripPath(argv[0]), true, common, option);
  //   exit(1);
  // }
  if ((opt->moving != 0 && opt->tau != 0) ||
      (opt->moving != 0 && opt->block != 0) ||
      (opt->tau != 0 && opt->block != 0)) {
    err_msg("only one of the -tau, -b, and -m option can be used");
    PrintError();
    Help(StripPath(argv[0]), true, common, option);
  }
  if (opt->moving != -1 && opt->c.end != -1 &&
      (opt->c.end - opt->c.start - opt->moving) < 0) {
    snprintf(ERROR_MSG, LINE, "nothing to compute: %s%d%s-point moving "
             "average from the total of %s%d%s datapoints", ErrYellow(),
             opt->moving, ErrRed(), ErrYellow(),
             opt->c.end - opt->c.start, ErrRed());
    PrintError();
    Help(StripPath(argv[0]), true, common, option);
  }

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // array to save the data; realloc'd on the fly
  double **data = malloc(sizeof *data * col_count);
  for (int i = 0; i < col_count; i++) {
    data[i] = malloc(sizeof *data[i]);
  }

  // read data from <input> file //{{{
  FILE *fr = OpenFile(fin, "r");

  int data_lines = 0, line_count = 0;
  while (true) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      break;
    }
    // if not empty line or comment continue //{{{
    if (words > 0 && split[0][0] != '#') {
      // error - insufficient number of columns //{{{
      if (words < col_max) {
        snprintf(ERROR_MSG, LINE, "too few columns (%s%d%s instead "
                 "of %s%d%s); file reading finished", ErrYellow(), words,
                 ErrCyan(), ErrYellow(), col_max, ErrCyan());
        PrintWarnFileLine(fin, line_count);
        break;
      } //}}}
      data_lines++;
      // save the value
      if (opt->c.start < data_lines) {
        count = data_lines - opt->c.start - 1;
        for (int i = 0; i < col_count; i++) {
          data[i] = s_realloc(data[i], sizeof *data[i] * (count + 1));
          data[i][count] = atof(split[column[i]-1]);
        }
      }
    } //}}}
    if (opt->c.end == data_lines) {
      break;
    }
  }
  fclose(fr); //}}}

  // error - <discard> is too large //{{{
  if (opt->c.start >= data_lines) {
    if (snprintf(ERROR_MSG, LINE, "number of lines to discard (%s%d%s) is "
                 "greater than the number of lines in %s%s%s", ErrYellow(),
                 opt->c.start, ErrRed(), ErrYellow(), fin, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  } //}}}

  data_lines -= opt->c.start;

  // -tau mode //{{{
  if (opt->tau > 0) {
    // variables
    // number of data points must be divisible by 'n_blocks'
    int remainder = (data_lines - opt->c.start) % opt->tau;
    // total number of data points to consider
    count = data_lines - opt->c.start - remainder;
    // number of data points per block
    int data_per_block = count / opt->tau;
    // block averages
    long double **avg_block = malloc(sizeof **avg_block * col_count);
    // overall averages
    long double **avg_all = malloc(sizeof **avg_all * col_count);
    for (int i = 0; i < col_count; i++) {
      avg_block[i] = calloc(opt->tau, sizeof *avg_block[i]);
      avg_all[i] = calloc(2, sizeof *avg_all[i]); // [0] avg, [1] avg of squares
    }

    int k = remainder; // first datapoint to consider
    for (int i = 0; i < opt->tau; i++) {
      for (int j = 0; j < data_per_block; j++) {
        for (int col = 0; col < col_count; col++) {
          avg_all[col][0] += data[col][k];
          avg_all[col][1] += Square(data[col][k]);
          avg_block[col][i] += data[col][k];
        }
        k++;
      }
    }
    for (int col = 0; col < col_count; col++) {
      avg_all[col][0] /= count;
      avg_all[col][1] /= count;
      for (int i = 0; i < opt->tau; i++) {
        avg_block[col][i] /= data_per_block;
      }
    }

    double *error = calloc(col_count, sizeof *error),
           *tau_int = calloc(col_count, sizeof *tau_int);
    for (int col = 0; col < col_count; col++) {
      // standard deviation for block averages
      double block_stdev = 0;
      for (int i = 0; i < opt->tau; i++) {
        block_stdev += Square(avg_block[col][i] - avg_all[col][0]);
      }
      block_stdev /= opt->tau - 1;
      // statistical error
      error[col] = sqrt(block_stdev / opt->tau);
      // approximate integrated autocorrelation time
      tau_int[col] = 0.5 * data_per_block * block_stdev /
                     (avg_all[0][1] - Square(avg_all[0][0]));
    }

    // print number of blocks, average, statistical error, and estimate of tau
    FILE *fw = OpenFile(fout, "a");
    fprintf(fw, " %6d", opt->tau);
    for (int col = 0; col < col_count; col++) {
      fprintf(fw, " %Lf", avg_all[col][0]);
      fprintf(fw, " %lf", error[0]);
      fprintf(fw, " %lf", tau_int[0]);
    }
    putc('\n', fw);
    fclose(fw);

    printf("%d blocks with %d datapoints\n", opt->tau, data_per_block);

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
  if (opt->block > 0) {
    // variables
    // number of data points must be divisible by 'data_per_block'
    int remainder = (data_lines - opt->c.start) % opt->block;
    // total number of data points to consider
    count = data_lines - opt->c.start - remainder;
    // number of blocks
    int blocks = count / opt->block;
    // block averages
    long double **avg_block = malloc(sizeof **avg_block * col_count);
    for (int i = 0; i < col_count; i++) {
      avg_block[i] = calloc(blocks, sizeof *avg_block[i]);
    }

    int k = remainder; // first datapoint to consider
    for (int i = 0; i < blocks; i++) {
      for (int j = 0; j < opt->block; j++) {
        for (int col = 0; col < col_count; col++) {
          avg_block[col][i] += data[col][k];
        }
        k++;
      }
    }
    for (int col = 0; col < col_count; col++) {
      for (int i = 0; i < blocks; i++) {
        avg_block[col][i] /= opt->block;
      }
    }
    FILE *fw = OpenFile(fout, "w");
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
  if (opt->moving > 0) {
    PrintByline(fout, argc, argv);
    // number of input datapoints
    count = data_lines - opt->c.start;
    // number of output datapoints
    int count_out = count - (opt->moving - 1);
    FILE *fw = OpenFile(fout, "a");
    for (int i = 0; i < count_out; i++) {
      for (int col = 0; col < col_count; col++) {
        double tmp = 0;
        for (int j = 0; j < opt->moving; j++) {
          tmp += data[col][i+j];
        }
        tmp /= opt->moving;
        fprintf(fw, " %e", tmp);
      }
      putc('\n', fw);
    }
    fclose(fw);
  } //}}}

  // calculate averages and errors and print them //{{{
  // 1) average
  double *avg = calloc(col_count, sizeof *avg);
  for (int col = 0; col < col_count; col++) {
    for (int i = 0; i < data_lines; i++) {
      avg[col] += data[col][i];
    }
    avg[col] /=  data_lines;
  }
  // 2) error
  // a) calculate the sum of squared differences from the mean
  double *std_dev = calloc(col_count, sizeof *std_dev);
  for (int col = 0; col < col_count; col++) {
    for (int i = 0; i < data_lines; i++) {
      std_dev[col] += Square(data[col][i] - avg[col]);
    }
    std_dev[col] = sqrt(std_dev[col] / (data_lines - 1));
  }
  // b) calculate the sample standard deviation and standard error
  double *std_err = calloc(col_count, sizeof *std_err);
  for (int col = 0; col < col_count; col++) {
    // double sig = std_dev[col];
    std_err[col] = std_dev[col] / sqrt(data_lines);
  }
  // print averages and errors to standard output
  for (int col = 0; col < col_count; col++) {
    // fprintf(stdout, "%lf %lf %lf\n", avg[col], std_dev[col], std_err[col]);
    fprintf(stdout, "%lf %lf\n", avg[col], std_err[col]);
  } //}}}

  for (int i = 0; i < col_count; i++) {
    free(data[i]);
  }
  free(data);
  free(column);
  free(avg);
  free(opt);

  return 0;
}
