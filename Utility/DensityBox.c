#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
DensityBox utility calculates number \
density for all bead types in the direction of specified axis (x, y, or z).\
The utility works properly only for orthogonal boxes that do not change \
size \n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <axis> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file\n");
  fprintf(ptr, "   <width>    width of a single bin\n");
  fprintf(ptr, "   <output>   output density file (automatic ending \
'<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>     calculate along x, y, or z axis\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  CommonHelp(error);
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
  int req_args = 4; //}}}

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
        strcmp(argv[i], "-vs_in") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-x") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char in_coor[LINE] = "",
       in_vsf[LINE] = "",
       in_lmp[LINE] = "",
       in_field[LINE] = "",
       *struct_file;
  int coor_type, struct_type = 0;
  snprintf(in_coor, LINE, "%s", argv[++count]);
  // coor_type: 1..vcf, 2..xyz, 3..lammpstrj
  coor_type = InputCoorStruct_old(argc, argv, in_coor, in_vsf, in_lmp, in_field);
  if (coor_type == -1) {
    exit(1);
  }
  // struct_type: 1..vsf, 2..lmp, 3..FIELD
  if (in_vsf[0] != '\0') {
    struct_type = 1;
    struct_file = in_vsf;
  } else if (in_lmp[0] != '\0') {
    struct_type = 2;
    struct_file = in_lmp;
  } else if (in_field[0] != '\0') {
    struct_type = 3;
    struct_file = in_field;
  } else {
    struct_file = "\0";
  }
  // error if no structure file specified (except when input is xyz)
  if (struct_type == 0 && coor_type != 2) {
    strcpy(ERROR_MSG, "missing input structure file; \
acceptable only for xyz input coordinate file");
    PrintError();
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal_old(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - filename with bead densities //{{{
  char output_rho[LINE] = "";
  snprintf(output_rho, LINE, "%s", argv[++count]); //}}}

  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {
    axis = argv[count][0];
    // Error - not x/y/z
    if (axis != 'x' && axis != 'y' && axis != 'z') {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "<axis>");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - use 'x', 'y', or 'z'\n\n");
      ColourReset(STDERR_FILENO);
      Help(argv[0], true);
      exit(1);
    }
    // add <axis>.rho to the output filename
    char str[LINE];
    output_rho[LINE-5] = '\0'; // ensure output_rho isn't too long
    snprintf(str, LINE, "%s%c.rho", output_rho, axis);
    strcpy(output_rho, str);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read input data //{{{
  SYSTEM System;
  switch(struct_type) {
    case 0: // xyz
      System = XyzReadStruct(in_coor);
      break;
    case 1: // vsf/vtf
      System = VtfReadStruct(in_vsf, detailed);
      break;
    case 2: // lmp
      System = LmpDataRead(in_lmp);
      break;
    case 3: // field
      System = FieldRead(in_field);
      break;
  }
  // pbc from vcf/vtf coordinate file
  // TODO: do I need pbc now? VtfSkipTimestep should read them (and so should
  //       LmpReadCoor)
  if (coor_type == 1) {
    VtfReadPBC(in_coor, &System.Box);
  }
  if (!TriclinicCellData(&System.Box, 1)) {
    strcpy(ERROR_MSG, "wrong pbc data");
    PrintError();
    exit(1);
  }
  WarnChargedSystem(System, in_vsf, in_field, in_lmp);
  // warn if missing box dimensions
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "unspecified box dimensions");
    PrintWarning();
    WarnPrintFile(struct_file, in_coor, "\0");
    putc('\n', stderr);
  } //}}}

  COUNT *Count = &System.Count;

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, &System)) {
    exit(1);
  } //}}}

  // write initial stuff to output density file //{{{
  FILE *out = OpenFile(output_rho, "w");
  PrintByline(out, argc, argv);
  // print bead type names to output file
  fprintf(out, "# columns: (1) distance;");
  for (int i = 0; i < Count->BeadType; i++) {
    fprintf(out, " (%d) %s", i+2, System.BeadType[i].Name);
    if (i != (Count->BeadType-1)) {
      putc(';', out);
    }
  }
  putc('\n', out);
  fclose(out); //}}}

  // number of bins //{{{
  // TODO: what about the box size? Read pbc from wherever
//double max_dist;
//switch(axis) {
//  case 'x':
//    max_dist = 2 * System.Box.Length.x;
//    break;
//  case 'y':
//    max_dist = 2 * System.Box.Length.y;
//    break;
//  case 'z':
//    max_dist = 2 * System.Box.Length.z;
//    break;
//}
//int bins = ceil(max_dist / width); //}}}
  int bins = ceil(200 / width); //}}}

// TODO: sizeof ...argh!
  // allocate memory for density arrays //{{{
//long int **rho = malloc(Count->BeadType*sizeof(long int *));
//long int **rho_2 = malloc(Count->BeadType*sizeof(long int *));
  long int **rho = malloc(Count->BeadType * sizeof **rho);
  long int **rho_2 = malloc(Count->BeadType * sizeof **rho_2);
  for (int i = 0; i < Count->BeadType; i++) {
//  rho[i] = calloc(bins, sizeof(long int));
//  rho_2[i] = calloc(bins, sizeof(long int));
    rho[i] = calloc(bins, sizeof *rho[i]);
    rho_2[i] = calloc(bins, sizeof *rho_2[i]);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  // open input coordinate file
  FILE *coor = OpenFile(in_coor, "r");

  // skip first (start-1) steps //{{{
  int file_line_count = 0; // count lines in the coordinate file
  int count_coor = 0; // count timesteps from the beginning
  for (int i = 1; i < start; i++) {
    count_coor++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count_coor);
    } //}}}
    if (!SkipTimestep(coor_type, coor, in_coor, in_vsf, &file_line_count)) {
      strcpy(ERROR_MSG, "premature end of file");
      PrintErrorFile(in_coor, "\0", "\0");
      putc('\n', stderr);
      exit(1);
    }
  }
  // print number of discarded steps?
  if (!silent && isatty(STDOUT_FILENO)) {
    fflush(stdout);
    fprintf(stdout, "\rDiscarded %d steps        \n", count_coor);
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps in the main loop
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count++;
    count_coor++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_coor);
    } //}}}
    if (!ReadTimestep(coor_type, coor, in_coor,
                      &System, &file_line_count, stuff)) {
      count_coor--;
      break;
    }
    RestorePBC(&System);

  // TODO: sizeof ...argh!
    // allocate memory for temporary density arrays //{{{
//  int **temp_rho = malloc(Count->BeadType*sizeof(int *));
    int **temp_rho = malloc(Count->BeadType * sizeof **temp_rho);
    for (int i = 0; i < Count->BeadType; i++) {
//    temp_rho[i] = calloc(bins,sizeof(int));
      temp_rho[i] = calloc(bins, sizeof temp_rho[i]);
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Count->BeadCoor; i++) {
      bool use = true;
      int id = System.BeadCoor[i];
      BEAD *bead = &System.Bead[id];
      int mol = bead->Molecule;
      if (mol != -1) { // do not use excluded molecules (-x option)
        int mtype = System.Molecule[mol].Type;
        use = System.MoleculeType[mtype].Use;
      }
      if (use) {
        if (axis == 'x') {
          int j = bead->Position.x / width;
          temp_rho[bead->Type][j]++;
        } else if (axis == 'y') {
          int j = bead->Position.y / width;
          temp_rho[bead->Type][j]++;
        } else {
          int j = bead->Position.z / width;
          temp_rho[bead->Type][j]++;
        }
      }
    } //}}}

    // add from temporary density array to global density arrays //{{{
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < bins; k++) {
        rho[j][k] += temp_rho[j][k];
        rho_2[j][k] += SQR(temp_rho[j][k]);
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Count->BeadType; i++) {
      free(temp_rho[i]);
    }
    free(temp_rho); //}}}

    if (count_coor == end) {
      break;
    }
  }
  fclose(coor);
  // print last step?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_coor);
  } //}}}

  // write densities to output file(s) //{{{
  out = OpenFile(output_rho, "a");

  // calculate rdf
  double volume = width;
  if (axis == 'x') {
    volume *= System.Box.Length.y * System.Box.Length.z;
  } else if (axis == 'y') {
    volume *= System.Box.Length.x * System.Box.Length.z;
  } else {
    volume *= System.Box.Length.x * System.Box.Length.y;
  }
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.3f", width*(2*i+1)/2);
    for (int j = 0; j < Count->BeadType; j++) {
      double temp_rho = rho[j][i] / (volume * count);
      // print average value to output file
      fprintf(out, " %10f", temp_rho);
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  free(stuff);
  for (int i = 0; i < Count->BeadType; i++) {
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
