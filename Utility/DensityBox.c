#include "../AnalysisTools.h"

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
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  int common = 12;
  char option[common][OPT_LENGTH];
  strcpy(option[ 0], "-st");
  strcpy(option[ 1], "-e");
  strcpy(option[ 2], "-sk");
  strcpy(option[ 3], "-i");
  strcpy(option[ 4], "--variable");
  strcpy(option[ 5], "-pbc");
  strcpy(option[ 6], "-ltrj");
  strcpy(option[ 7], "--detailed");
  strcpy(option[ 8], "-v");
  strcpy(option[ 9], "--silent");
  strcpy(option[10], "--help");
  strcpy(option[11], "--version");
  CommonHelp(error, common, option);
} //}}}

int main(int argc, char *argv[]) {

  // --help/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
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
    if (argv[i][0] == '-' && strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-st") != 0 && strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 && strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "--variable") != 0 && strcmp(argv[i], "-pbc") != 0 &&
        strcmp(argv[i], "-ltrj") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-v") != 0 && strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "--help") && strcmp(argv[i], "--version") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

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
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - filename with bead densities //{{{
  char output_rho[LINE] = "";
  snprintf(output_rho, LINE, "%s", argv[++count]); //}}}

  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {
    axis = argv[count][0];
    // Error - not x/y/z
    if (axis != 'x' && axis != 'y' && axis != 'z') {
      strcpy(ERROR_MSG, "<axis> must be 'x', 'y', or 'z'");
      PrintError();
      Help(argv[0], true);
      exit(1);
    }
    // add <axis>.rho to the output filename
    char str[LINE];
    output_rho[LINE-5] = '\0'; // ensure output_rho isn't too long
    if (snprintf(str, LINE, "%s%c.rho", output_rho, axis) < 0) {
      strcpy(ERROR_MSG, "something wrong with snprintf()");
      exit(1);
    }
    strcpy(output_rho, str);
  } //}}}

  // options before reading system data
  bool silent, verbose, detailed, vtf_var;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1, ltrj_start_id = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &ltrj_start_id, &start, &end, &skip);

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz, &ltrj_start_id);

  COUNT *Count = &System.Count;

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, &System)) {
    exit(1);
  } //}}}

  // number of bins //{{{
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintErrorFile(coor_file, struct_file, "\0");
    exit(1);
  }
  double max_dist;
  switch(axis) {
   case 'x':
     max_dist = System.Box.Length.x;
     break;
   case 'y':
     max_dist = System.Box.Length.y;
     break;
   case 'z':
     max_dist = System.Box.Length.z;
     break;
  }
  int bins = ceil(max_dist / width) + 10; //}}}

  // allocate memory for arrays //{{{
  long int **rho = malloc(Count->BeadType * sizeof **rho);
  long int **rho_2 = malloc(Count->BeadType * sizeof **rho_2);
  long int *n_beads = calloc(Count->BeadType, sizeof *n_beads);
  for (int i = 0; i < Count->BeadType; i++) {
    rho[i] = calloc(bins, sizeof *rho[i]);
    rho_2[i] = calloc(bins, sizeof *rho_2[i]);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  // main loop //{{{
  FILE *fr = OpenFile(coor_file, "r");
  int count_coor = 0,     // count steps in the vcf file
      count_used = 0,    // count steps in output file
      line_count = 0;     // count lines in the vcf file
  count = 0; // count timesteps in the main loop
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    PrintStep(&count_coor, start, silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (count_coor >= start && (count_coor <= end || end == -1) &&
        ((count_coor - start) % skip) == 0) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(coor_type, fr, coor_file, &System, &line_count,
                        ltrj_start_id, vtf_var)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
      // allocate memory for temporary density arrays
      int **temp_rho = malloc(Count->BeadType * sizeof(int *));
      // WTF? ...the next line produces error!
      // int **temp_rho = malloc(Count->BeadType * sizeof **temp_rho);
      for (int i = 0; i < Count->BeadType; i++) {
        temp_rho[i] = calloc(bins, sizeof *temp_rho[i]);
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
          n_beads[bead->Type]++;
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
      // add from temporary density arrays to global density arrays
      for (int j = 0; j < Count->BeadType; j++) {
        for (int k = 0; k < bins-1; k++) {
          rho[j][k] += temp_rho[j][k];
          rho_2[j][k] += SQR(temp_rho[j][k]);
        }
      }
      // free temporary density array
      for (int i = 0; i < Count->BeadType; i++) {
        free(temp_rho[i]);
      }
      free(temp_rho);
    } else {
      if (!SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count)) {
        count_coor--;
        break;
      }
    }
    // TODO: isn't this superfluous, considering the 'use' above?
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
  // write initial stuff to output density file //{{{
  FILE *fw = OpenFile(output_rho, "w");
  PrintByline(fw, argc, argv);
  // print bead type names to output file
  fprintf(fw, "# columns: (1) distance");
  count = 1;
  for (int i = 0; i < Count->BeadType; i++) {
    if (n_beads[i] != 0) {
      count++;
      fprintf(fw, "; (%d) %s", count, System.BeadType[i].Name);
    }
  }
  putc('\n', fw); //}}}

  // calculate rdf
  double volume = width;
  if (axis == 'x') {
    volume *= System.Box.Length.y * System.Box.Length.z;
  } else if (axis == 'y') {
    volume *= System.Box.Length.x * System.Box.Length.z;
  } else {
    volume *= System.Box.Length.x * System.Box.Length.y;
  }
  for (int i = 0; i < bins-1; i++) {
    fprintf(fw, "%7.3f", width*(2*i+1)/2);
    for (int j = 0; j < Count->BeadType; j++) {
      if (n_beads[j] > 0 ){
        double temp_rho = rho[j][i] / (volume * count_used);
        fprintf(fw, " %10f", temp_rho);
      }
    }
    putc('\n',fw);
  }
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  free(stuff);
  for (int i = 0; i < Count->BeadType; i++) {
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho);
  free(n_beads); //}}}

  return 0;
}
