#include "../AnalysisTools.h"

// TODO: write function ReadXyzCoordinates()

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Config_from_xyz generates dl_meso CONFIG file from given step of a xyz file. \
If the given timestep is larger than the number of steps the coordinate file, \
the last step is used. Note that box dimensions must be added manually to the \
resulting CONFIG file, as xyz file does not contain box size.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.xyz> [options]\n\n", cmd);

  fprintf(ptr, "   <input.xyz>   filename (xyz format)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -st <step>   timestep for creating CONFIG \
(default: last)\n");
  fprintf(ptr, "      --silent     no output\n");
  fprintf(ptr, "      -h           print this help and exit\n");
  fprintf(ptr, "      --version    print version number and exit\n");
} //}}}

int main(int argc, char *argv[]) {

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
  int req_args = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
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
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input.xyz> - input coordinate file (must end with .xyz) //{{{
  char input_xyz[LINE] = "";
  snprintf(input_xyz, LINE, "%s", argv[++count]);

  // test if <input> filename ends with '.xyz' (required by VMD)
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".xyz");
  if (ErrorExtension(input_xyz, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  bool verbose, silent;
  SilentOption(argc, argv, &verbose, &silent);

  // timestep to create CONFIG file from //{{{
  int timestep = -1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // get number of beads from xyz file //{{{
  FILE *xyz = OpenFile(input_xyz, "r");
  int beads;
  if (fscanf(xyz, "%d", &beads) != 1) {
    // TODO: correct colours etc.
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "cannot read number of beads from \033[1;33m%s\033[1;31m\n\n", input_xyz);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  fclose(xyz); //}}}

  // main loop - read timesteps till the chosen one //{{{
  xyz = OpenFile(input_xyz, "r");

  fpos_t pos; // for saving pointer position in xyz file
  int test;
  count = 0;

  while ((test = getc(xyz)) != EOF && count != timestep) {
    ungetc(test, xyz);

    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    // skip remainder of number-of-beads line
    while (getc(xyz) != '\n')
      ;
    // skip comment line
    while (getc(xyz) != '\n')
      ;

    // save pointer position in file
    fgetpos(xyz, &pos);

    // skip timestep
    for (int i = 0; i < beads; i++) {
      while (getc(xyz) != '\n')
        ;
    }
  } //}}}

  // print config step number? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Config Step: %d\n", count);
  } //}}}

  // read the coordinates from the chosen timestep and create CONFIG //{{{
  // restore pointer position in xyz file
  fsetpos(xyz, &pos);
  FILE *config = OpenFile("CONFIG", "w");

  // print first stuff to CONFIG
  fprintf(config, "CONFIG from %s\n0 1\n", input_xyz);
  fprintf(config, "x 0 0\n");
  fprintf(config, "0 y 0\n");
  fprintf(config, "0 0 z\n");

  for (int i = 0; i < beads; i++) {
    char line[LINE];
    fgets(line, sizeof(line), xyz);
    char split[SPL_STR][SPL_LEN], delim[8];
    strcpy(delim, " \t");
    int words = SplitLine(split, line, delim);

    // error - less then four whitespace-separated strings //{{{
    if (words < 4) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", input_xyz);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - not enough columns in timestep %d\n", count);
      ColourReset(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}

    // test if split[1-3] are doubles //{{{
    for (int j = 1; j < 4; j++) {
      if (!IsReal(split[j])) {
        // TODO colours
        ErrorPrintError_old();
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\033[1;33m%s\033[1;31m", input_xyz);
        fprintf(stderr, " - non-numeric coordinate in %d. timestep\n\n", count);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
    } //}}}

    // print bead to CONFIG
    fprintf(config, "%s %d\n", split[0], i+1);
    fprintf(config, "%lf %lf %lf\n", atof(split[1]),
                                     atof(split[2]),
                                     atof(split[3]));
  }
  fclose(config);
  fclose(xyz); //}}}

  return 0;
}
