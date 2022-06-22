#include "../AnalysisTools.h"
void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info prints information about the provided system and can also create a new \
vtf structure file. The new structure file can differ from the original one if \
--detailed switch is used. \n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> [options]\n\n", cmd);
  fprintf(ptr, "   <input>   vtf input structure file\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --detailed        differentiate bead types not just \
by names\n");
  fprintf(ptr, "      -vsf <file.vsf>   create a new vsf structure file\n");
  fprintf(ptr, "      -c <file>         input coordinate file\n");
  fprintf(ptr, "      -f <file>         input FIELD-like file for extra \
structural information\n");
  fprintf(ptr, "      -v                more verbose output\n");
  fprintf(ptr, "      -h                print this help and exit\n");
  fprintf(ptr, "      --version         print version number and exit\n");
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
        strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-vsf") != 0 &&
        strcmp(argv[i], "-c") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // <input> - input structure file (must end with .vsf or .vtf) //{{{
  char input_vsf[LINE] = "";
  snprintf(input_vsf, LINE, "%s", argv[++count]);
  // test if <input> filename ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data
  bool verbose = BoolOption(argc, argv, "-v");
  bool detailed = BoolOption(argc, argv, "--detailed");

  // -vsf option //{{{
  char output_vsf[LINE] = "\0";
  if (FileOption(argc, argv, "-vsf", output_vsf, LINE)) {
    exit(1);
  }
  if (output_vsf[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".vsf");
    if (ErrorExtension(output_vsf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // -c option //{{{
  char input_coor[LINE] = "\0";
  if (FileOption(argc, argv, "-c", input_coor, LINE)) {
    exit(1);
  }
  if (input_coor[0] != '\0') {
    ext = 2;
    strcpy(extension[0], ".vcf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(input_coor, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // -f option //{{{
  char input_field[LINE] = "\0";
  if (FileOption(argc, argv, "-f", input_field, LINE)) {
    exit(1);
  } //}}}

  // read information from input file(s)
  SYSTEM System = VtfReadStruct(input_vsf, detailed);
  if (input_field[0] != '\0') {
    printf("%sfield: %s%s\n", Green(), input_field, ColourReset());
    SYSTEM field = FieldReadFull(input_field);
    PrintCount(field.Count);
  } else {
    printf("%sNO FIELD%s\n", Magenta(), ColourReset());
  }
  if (input_coor[0] != '\0') {
    FILE *coor = OpenFile(input_coor, "r");
    char stuff[LINE];
    int step_count = 0, file_line_count = 0;
    VtfReadTimestep(coor, input_coor, input_vsf, &System, &file_line_count,
                    step_count, stuff);
    fclose(coor);
    PrintBox(System.Box);
  }

  // print information
  VerboseOutput(System);
  if (verbose) { // -v option
    fprintf(stdout, "\nInformation about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(System);
  }

  if (output_vsf[0] != '\0') {
    VtfWriteStruct(output_vsf, System);
  }

  FreeSystem(&System);

  return 0;
}
