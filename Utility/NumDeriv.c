#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "TBA\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <x> <f(x)> <mode>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename\n");
  fprintf(ptr, "   <x>               column number of x values\n");
  fprintf(ptr, "   <f(x)>            column number of f(x) values\n");
  fprintf(ptr, "   <mode>            0=central, 1=left, 2=right\n");
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

  PrintCommand(stdout, argc, argv);

  // <x> //{{{
  // Error - non-numeric argument
  if (!IsInteger(argv[++count])) {
    ErrorNaN("<x>");
    Help(argv[0], true);
    exit(1);
  }
  int x = atoi(argv[count]); //}}}

  // <f(x)> //{{{
  // Error - non-numeric argument
  if (!IsInteger(argv[++count])) {
    ErrorNaN("<f(x)>");
    Help(argv[0], true);
    exit(1);
  }
  int f_x = atoi(argv[count]); //}}}

  // <mode> //{{{
  // Error - non-numeric argument
  if (!IsInteger(argv[++count])) {
    ErrorNaN("<n_blocks>");
    Help(argv[0], true);
    exit(1);
  }
  int mode = atoi(argv[count]); //}}}

  // get number of data lines in <input> //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  }
  int data_points = 0; // data lines only
  int lines = 0; // all lines
  while (true) {
    // get & split line
    char line[LINE], split[SPL_STR][SPL_LEN], delim[8];
    fgets(line, sizeof(line), fr);
    strcpy(delim, " \t");
    int words = SplitLine(split, line, delim);
    // exit on end of file
    if(feof(fr)) {
      break;
    }
    lines++;
    // if not empty line or comment continue
    if (split[0][0] != '#' &&
        split[0][0] != '\n') {
      // error - insufficient number of columns //{{{
      if (words < x || words < f_x) {
        ErrorPrintError_old();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - too few columns");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      data_points++;
    }
  }
  fclose(fr); //}}}

  // read data from <input> file //{{{
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  }
  count = 0;
  double **data = calloc(data_points, sizeof(double *));
  for (int i = 0; i < data_points; i++) {
    data[i] = calloc(2, sizeof(double));
  }
  for (int i = 0; i < lines; i++) {
    char line[LINE], split[SPL_STR][SPL_LEN], delim[8];
    fgets(line, sizeof(line), fr);
    strcpy(delim, " \t");
    SplitLine(split, line, delim);
    if (split[0][0] != '#' &&
        split[0][0] != '\n') {
      data[count][0] = atof(split[x-1]);
      data[count][1] = atof(split[f_x-1]);
      count++;
    }
  }
  fclose(fr); //}}}

  // calculate numerical derivative //{{{
  for (int i = 0; i < data_points; i++) {
    if (mode == 0) { // (f(x+a)-f(x-a))/2a
      if (i == (data_points-1)) {
        break;
      }
      if (i > 0) {
        printf("%lf %lf\n", data[i][0],
                            (data[i+1][1]-data[i-1][1])/(data[i+1][0]-data[i-1][0]));
      }
    } else if (mode == 1) { // (f(x)-f(x-a))/a
      if (i > 0) {
        printf("%lf %lf\n", data[i-1][0]+(data[i][0]-data[i-1][0])/2,
                            (data[i][1]-data[i-1][1])/(data[i][0]-data[i-1][0]));
      }
    } else { // (f(x+a)-f(x))/a
      if (i == (data_points-1)) {
        break;
      }
      printf("%lf %lf\n", data[i][0]+(data[i+1][0]-data[i][0])/2,
                          (data[i+1][1]-data[i][1])/(data[i+1][0]-data[i][0]));
    }
  } //}}}

  return 0;
}
