#include "../AnalysisTools.h"

// TODO: possible changing box size: make bins' width variable, keeping their
//       number, and work in relative coordinates (relative to instantaneous
//       dimensions) throughout the code

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
DensityBox utility calculates number \
density for all bead types in the direction of all axes (x, y, and z).\
The utility works properly only for orthogonal boxes that do not change \
size.\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <width> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single bin\n");
  fprintf(ptr, "<output>            output density files (automatic ending "
          "-x.rho, -y.rho, and -z.rho)\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -x <name(s)>      exclude specified molecule(s)\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 10, all = common + 1, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  // strcpy(option[count++], "--variable"); // TODO: makes no sense, I think
  strcpy(option[count++], "-pbc");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-x");
  OptionCheck(argc, argv, req_arg, common, all, option);
  //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate (and structure) file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type = -1, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <outputt> - filename
  char output_rho[LINE] = "";
  snprintf(output_rho, LINE, "%s", argv[++count]);
  output_rho[LINE-6] = '\0'; // for adding -<axis>.rho

  // options before reading system data
  bool silent, verbose, detailed;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &start, &end, &skip);

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file,
                                coor_type, coor_file, detailed, pbc_xyz);
  COUNT *Count = &System.Count;
  BOX *box = &System.Box;

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, &System)) {
    exit(1);
  } //}}}

  // number of bins //{{{
  if (box->Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintErrorFile(coor_file, struct_file, "\0");
    exit(1);
  }
  VECTOR bin;
  // TODO: *3 to assume box change of at most thrice as big
  //       probably change from width to number of bins per box?
  bin.x = ceil(box->Length.x / width) * 3;
  bin.y = ceil(box->Length.y / width) * 3;
  bin.z = ceil(box->Length.z / width) * 3; //}}}

  // allocate memory for arrays //{{{
  // just check if the bead type is at all present in the calculation
  bool *n_beads = calloc(Count->BeadType, sizeof *n_beads);
  // TODO: LONGINTVECTOR instead long int [3]?
  long int ***rho = malloc(3 * sizeof ***rho);
  for (int i = 0; i < 3; i++) {
    rho[i] = malloc(Count->BeadType * sizeof **rho);
  }
  for (int j = 0; j < Count->BeadType; j++) {
    rho[0][j] = calloc(bin.x, sizeof *rho[0][j]);
    rho[1][j] = calloc(bin.y, sizeof *rho[1][j]);
    rho[2][j] = calloc(bin.z, sizeof *rho[2][j]);
  } //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  // main loop //{{{
  FILE *fr = OpenFile(coor_file, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, start, silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (count_coor >= start && (count_coor <= end || end == -1) &&
        ((count_coor - start) % skip) == 0) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(coor_type, fr, coor_file, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
      // allocate memory for temporary density arrays
      int ***temp_rho = malloc(3 * sizeof(int **));
      for (int i = 0; i < 3; i++) {
        temp_rho[i] = malloc(Count->BeadType * sizeof(int *));
      }
      for (int i = 0; i < Count->BeadType; i++) {
        temp_rho[0][i] = calloc(bin.x, sizeof *temp_rho[0][i]);
        temp_rho[1][i] = calloc(bin.y, sizeof *temp_rho[1][i]);
        temp_rho[2][i] = calloc(bin.z, sizeof *temp_rho[2][i]);
      }

      // calculate densities //{{{
      for (int i = 0; i < Count->BeadCoor; i++) {
        use = true;
        int id = System.BeadCoor[i];
        BEAD *bead = &System.Bead[id];
        int mol = bead->Molecule;
        if (mol != -1) { // do not use excluded molecules (-x option)
          int mtype = System.Molecule[mol].Type;
          use = System.MoleculeType[mtype].Flag;
        }
        if (use) {
          n_beads[bead->Type] = true;
          // x direction
          int j = bead->Position.x / width;
          temp_rho[0][bead->Type][j]++;
          // y direction
          j = bead->Position.y / width;
          temp_rho[1][bead->Type][j]++;
          // z direction
          j = bead->Position.z / width;
          temp_rho[2][bead->Type][j]++;
        }
      } //}}}
      // add from temporary density arrays to global density arrays
      for (int j = 0; j < Count->BeadType; j++) {
        // x direction
        for (int k = 0; k < bin.x-1; k++) {
          rho[0][j][k] += temp_rho[0][j][k];
        }
        // y direction
        for (int k = 0; k < bin.y-1; k++) {
          rho[1][j][k] += temp_rho[1][j][k];
        }
        // z direction
        for (int k = 0; k < bin.z-1; k++) {
          rho[2][j][k] += temp_rho[2][j][k];
        }
      }
      // free temporary density array
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < Count->BeadType; j++) {
          free(temp_rho[i][j]);
        }
        free(temp_rho[i]);
      }
      free(temp_rho);
    } else {
      if (!SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count)) {
        count_coor--;
        break;
      }
    }
    if (count_coor == end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // write densities to output file(s) //{{{
  for (int ax = 0; ax < 3; ax++) {
    // axis-based variables
    double volume = width, size = -1;
    int n = 0; // number of bins
    char axis;
    if (ax == 0) {
      axis = 'x';
      size = box->Length.x;
      volume *= box->Length.y * box->Length.z;
      n = bin.x;
    } else if (ax == 1) {
      axis = 'y';
      size = box->Length.y;
      volume *= box->Length.x * box->Length.z;
      n = bin.y;
    } else {
      axis = 'z';
      size = box->Length.z;
      volume *= box->Length.x * box->Length.y;
      n = bin.z;
    }
    char file[LINE]; // filename <output>-<axis>.rho
    if (snprintf(file, LINE, "%s-%c.rho", output_rho, axis) < 0) {
      ErrorSnprintf();
    }
    // write initial stuff to output density file
    PrintByline(file, argc, argv);
    FILE *fw = OpenFile(file, "a");
    // print bead type names to output file
    fprintf(fw, "# columns: (1) distance");
    count = 1;
    for (int i = 0; i < Count->BeadType; i++) {
      if (n_beads[i]) {
        count++;
        fprintf(fw, "; (%d) %s", count, System.BeadType[i].Name);
      }
    }
    putc('\n', fw);
    // write rdf
    for (int i = 0; i < (n - 1); i++) {
      double dist = width * (2 * i + 1) / 2;
      if (dist > size) { // write only til the max box size
        break;
      }
      fprintf(fw, "%7.3f", dist); // absolute distance
      for (int j = 0; j < Count->BeadType; j++) {
        if (n_beads[j]){
          double temp_rho = rho[ax][j][i] / (volume * count_used);
          fprintf(fw, " %10f", temp_rho);
        }
      }
      putc('\n',fw);
    }
    fclose(fw);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      free(rho[i][j]);
    }
    free(rho[i]);
  }
  free(rho);
  free(n_beads); //}}}

  return 0;
}
