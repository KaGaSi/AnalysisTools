#include "../src/AnalysisTools.h"
// TODO: remove <output> in favour of -tau/-b/-m <output> <int> options
// TODO: test it's counting correctly the number of data lines

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Average utility calculates averages for specified column(s) from the input "
  "file (ignoring empty lines and lines beginning with '#'), printing "
  " the results to an output file. It has three operation modes based on which "
  " of the three options is supplied:\n"
  "1) for -tau option, the utility uses binning method to calculate average, "
  "statistical error and an estimate of integrated autocorrelation "
  "time (tau), outputting the specified number of blocks and three values for "
  "each column used: <simple average> <error> <tau> "
  "(all on a single line). In this mode, Average appends to the "
  "output file instead of rewriting it. See the manual for a way to obtain "
  "a reasonable estimate of tau via rerunning the utility several times.\n"
  "2) for -b option, the binning method is used to calculate per-block "
  "averages.\n"
  "3) for -m option, the moving method is used to smoothen the input data.",

  "Usage: Average <input> <output> <column(s)>",
  .args = 3, // number of mandatory arguments
  .all = 11, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  // COMMON_OPTS[C_SK], // TODO: this should be there too, right?
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", NULL, "input filename", OPT_ARG},
  {"<output>", NULL, "output filename", OPT_ARG},
  {"<column(s)>", NULL, "column number(s) to analyse", OPT_ARG},
  {"-tau", "<int>", "estimate tau mode - number of blocks to split data into", OPT_EXTRA},
  {"-b", "<int>", "block mode - number of datapoints per block", OPT_EXTRA},
  {"-m", "<int>", "moving mode - number of data points per moving average", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  int tau, block, moving; // -tau -b -m
}; //}}}

int main ( int argc, char** argv ) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input>
  char fin[LINE];
  s_strcpy(fin, argv[++count], LINE);
  // <output>
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // <column> - column number(s) to analyze
  // TODO: warning if multiple times the same column number
  long int *column = malloc(sizeof *column);
  if (!column) {
    ErrorAlloc("column");
  }
  int col_count = 0;
  while (++count < argc && argv[count][0] != '-') {
    if (!IsNaturalNumber(argv[count], &column[col_count])) {
      ErrorNaN("<column>");
      Help(true, HelpDesc, opts);
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
  COMMON_OPT commons = CommonOptions(argc, argv, trash);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // -tau option: use block method to get overall average (and std_err and tau)
  opt.tau = 0;
  OneNumberOption(argc, argv, "-tau", &opt.tau, 'i');
  // -b option: calculate block averages
  opt.block = 0;
  OneNumberOption(argc, argv, "-b", &opt.block, 'i');
  // -m option: calculate moving average
  opt.moving = 0;
  OneNumberOption(argc, argv, "-m", &opt.moving, 'i');
  if ((opt.moving == 0 && opt.tau == 0 && opt.block == 0) ||
      (opt.moving != 0 && opt.tau != 0) ||
      (opt.moving != 0 && opt.block != 0) ||
      (opt.tau != 0 && opt.block != 0)) {
    err_msg("exactly one of the -tau, -b, and -m option must be used");
    PrintError();
    Help(true, HelpDesc, opts);
  }
  if (opt.moving != -1 && commons.end != -1 &&
      (commons.end - commons.start - opt.moving) < 0) {
    snprintf(ERROR_MSG, LINE, "nothing to compute: %s%d%s-point moving "
             "average from the total of %s%d%s datapoints", ErrYellow(),
             opt.moving, ErrRed(), ErrYellow(),
             commons.end - commons.start, ErrRed());
    PrintError();
    Help(true, HelpDesc, opts);
  } //}}}

  // array to save the data; realloc'd on the fly
  double **data = malloc(sizeof *data * col_count);
  if (!data) {
    ErrorAlloc("data");
  }
  for (int i = 0; i < col_count; i++) {
    data[i] = malloc(sizeof *data[i]);
    if (!data[i]) {
      ErrorAlloc("data[i]");
    }
  }

  // read data from <input> file //{{{
  FILE *fr = OpenFile(fin, "r");

  int data_lines = 0, line_count = 0;
  while (true) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      break;
    }
    // if not empty line or comment continue
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
      if (commons.start < data_lines) {
        count = data_lines - commons.start - 1;
        for (int i = 0; i < col_count; i++) {
          data[i] = s_realloc(data[i], sizeof *data[i] * (count + 1));
          data[i][count] = atof(split[column[i]-1]);
        }
      }
    }
    if (commons.end == data_lines) {
      break;
    }
  }
  fclose(fr); //}}}

  // error - starting line is too large //{{{
  if (commons.start > data_lines) {
    if (snprintf(ERROR_MSG, LINE, "starting data line (%s%d%s) is "
                 "greater than the number of lines in %s%s%s", ErrYellow(),
                 commons.start, ErrRed(), ErrYellow(), fin, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  } //}}}

  data_lines -= commons.start;

  // -tau mode //{{{
  if (opt.tau > 0) {
    // variables
    // number of data points must be divisible by 'n_blocks'
    int remainder = (data_lines - commons.start) % opt.tau;
    // total number of data points to consider
    count = data_lines - commons.start - remainder;
    // number of data points per block
    int data_per_block = count / opt.tau;
    // block averages
    ArrNDld *avg_block = CreateArr2Dld(col_count, opt.tau);
    // overall averages
    ArrNDld *avg_all = CreateArr2Dld(col_count, 2);
    if (!avg_block || !avg_all) {
      ErrorAlloc("avg_block/avg_all");
    }

    int k = remainder; // first datapoint to consider
    for (int i = 0; i < opt.tau; i++) {
      for (int j = 0; j < data_per_block; j++) {
        for (int col = 0; col < col_count; col++) {
          AddArr2D(avg_all, col, 0, data[col][k]);
          AddArr2D(avg_all, col, 1, Square(data[col][k]));
          AddArr2D(avg_block, col, i, data[col][k]);
        }
        k++;
      }
    }
    for (int col = 0; col < col_count; col++) {
      SetArr2D(avg_all, col, 0, GetArr2D(avg_all, col, 0) / count);
      SetArr2D(avg_all, col, 1, GetArr2D(avg_all, col, 1) / count);
      for (int i = 0; i < opt.tau; i++) {
        double val = GetArr2D(avg_block, col, i) / data_per_block;
        SetArr2D(avg_block, col, i, val);
      }
    }

    double *error = calloc(col_count, sizeof *error),
           *tau_int = calloc(col_count, sizeof *tau_int);
    if (!error || !tau_int) {
      ErrorAlloc("error/tau_int");
    }
    for (int col = 0; col < col_count; col++) {
      // standard deviation for block averages
      double block_stdev = 0;
      for (int i = 0; i < opt.tau; i++) {
        double val = GetArr2D(avg_block, col, i) - GetArr2D(avg_all, col, 0);
        block_stdev += Square(val);
      }
      block_stdev /= opt.tau - 1;
      // statistical error
      error[col] = sqrt(block_stdev / opt.tau);
      // approximate integrated autocorrelation time
      double val = GetArr2D(avg_all, 0, 1) - Square(GetArr2D(avg_all, 0, 0));
      tau_int[col] = 0.5 * data_per_block * block_stdev / val;
    }

    // print number of blocks, average, statistical error, and estimate of tau
    FILE *fw = OpenFile(fout, "a");
    fprintf(fw, " %6d", opt.tau);
    for (int col = 0; col < col_count; col++) {
      fprintf(fw, " %Lf", GetArr2D(avg_all, col, 0));
      fprintf(fw, " %lf", error[0]);
      fprintf(fw, " %lf", tau_int[0]);
    }
    putc('\n', fw);
    fclose(fw);

    printf("%d blocks with %d datapoints\n", opt.tau, data_per_block);

    FreeArrND(avg_block);
    FreeArrND(avg_all);
    free(error);
    free(tau_int);
  } //}}}

  // -b mode //{{{
  if (opt.block > 0) {
    // variables
    // number of data points must be divisible by 'data_per_block'
    int remainder = (data_lines - commons.start) % opt.block;
    // total number of data points to consider
    count = data_lines - commons.start - remainder;
    // number of blocks
    int blocks = count / opt.block;
    // block averages
    ArrNDld *avg_block = CreateArr2Dld(col_count, blocks);
    if (!avg_block) {
      ErrorAlloc("avg_block");
    }

    int k = remainder; // first datapoint to consider
    for (int i = 0; i < blocks; i++) {
      for (int j = 0; j < opt.block; j++) {
        for (int col = 0; col < col_count; col++) {
          AddArr2D(avg_block, col, i, data[col][k]);
        }
        k++;
      }
    }
    for (int col = 0; col < col_count; col++) {
      for (int i = 0; i < blocks; i++) {
        SetArr2D(avg_block, col, i, GetArr2D(avg_block, col, i) / opt.block);
      }
    }
    FILE *fw = OpenFile(fout, "w");
    for (int i = 0; i < blocks; i++) {
      for (int col = 0; col < col_count; col++) {
        fprintf(fw, " %Le", GetArr2D(avg_block, col, i));
      }
      putc('\n', fw);
    }
    fclose(fw);

    FreeArrND(avg_block);
  } //}}}

  // -m mode //{{{
  if (opt.moving > 0) {
    FILE *fw = PrintBylineOpenFile(fout, argc, argv);
    // number of input datapoints
    count = data_lines - commons.start;
    // number of output datapoints
    int count_out = count - (opt.moving - 1);
    for (int i = 0; i < count_out; i++) {
      for (int col = 0; col < col_count; col++) {
        double tmp = 0;
        for (int j = 0; j < opt.moving; j++) {
          tmp += data[col][i+j];
        }
        tmp /= opt.moving;
        fprintf(fw, " %e", tmp);
      }
      putc('\n', fw);
    }
    fclose(fw);
  } //}}}

  // calculate averages and errors and print them //{{{
  // 1) average
  double *avg = calloc(col_count, sizeof *avg);
  if (!avg) {
    ErrorAlloc("avg");
  }
  for (int col = 0; col < col_count; col++) {
    for (int i = 0; i < data_lines; i++) {
      avg[col] += data[col][i];
    }
    avg[col] /=  data_lines;
  }
  // 2) error
  // a) calculate the sum of squared differences from the mean
  double *std_dev = calloc(col_count, sizeof *std_dev);
  if (!std_dev) {
    ErrorAlloc("std_dev");
  }
  for (int col = 0; col < col_count; col++) {
    for (int i = 0; i < data_lines; i++) {
      std_dev[col] += Square(data[col][i] - avg[col]);
    }
    std_dev[col] = sqrt(std_dev[col] / (data_lines - 1));
  }
  // b) calculate the sample standard deviation and standard error
  double *std_err = calloc(col_count, sizeof *std_err);
  if (!std_err) {
    ErrorAlloc("std_der");
  }
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

  return 0;
}
