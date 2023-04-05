#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Density2D utility calculates number \
density for all bead types for a plane perpendicular to the specified axis \
direction (x, y, or z).\
\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <axis> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename \
(either vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single bin \
(in the two planar directions)\n");
  fprintf(ptr, "   <output>          output density file \
(automatic ending '<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>            calculate in plane perpendicular to x, y, \
or z axis\n");
  fprintf(ptr, "   <options>\n");
//fprintf(ptr, "      -n <int>       number of bins to average\n");
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
    if (argv[i][0] == '-' && strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-st") != 0 && strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 && strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "--variable") != 0 && strcmp(argv[i], "-pbc") != 0 &&
        strcmp(argv[i], "-ltrj") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-v") != 0 && strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "--help") != 0 && strcmp(argv[i], "--version") != 0) {
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

  // <output.rho> - filename with bead densities
  char output_rho[LINE];
  strcpy(output_rho, argv[++count]);

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

  // <axis> - x, y, or z //{{{
  // TODO: WTF? Why while loop when it's just one character?
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {
    axis = argv[count][0];
    if (axis != 'x' && axis != 'y' && axis != 'z') {
      strcpy(ERROR_MSG, "<axis> must be 'x', 'y', or 'z'");
      PrintError();
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // full output file name //{{{
  char str[LINE];
  if (snprintf(str, LINE, "%s%c.rho", output_rho, axis) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    exit(1);
  }
  strcpy(output_rho, str); //}}}

  // number of bins //{{{
  double max_dist[2];
  if (axis == 'x') {
    max_dist[0] = System.Box.Length.y;
    max_dist[1] = System.Box.Length.z;
  } else if (axis == 'y') {
    max_dist[0] = System.Box.Length.x;
    max_dist[1] = System.Box.Length.z;
  } else {
    max_dist[0] = System.Box.Length.x;
    max_dist[1] = System.Box.Length.y;
  }
  int bins[2];
  bins[0] = ceil(max_dist[0] / width);
  bins[1] = ceil(max_dist[1] / width); //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  // allocate memory for density arrays //{{{
  long int ***rho = malloc(Count->BeadType * sizeof(long int **));
  long int ***rho_2 = malloc(Count->BeadType * sizeof(long int **));
  for (int i = 0; i < Count->BeadType; i++) {
    rho[i] = calloc(bins[0], sizeof(long int *));
    rho_2[i] = calloc(bins[0], sizeof(long int *));
    for (int j = 0; j < bins[0]; j++) {
      rho[i][j] = calloc(bins[1], sizeof(long int));
      rho_2[i][j] = calloc(bins[1], sizeof(long int));
    }
  } //}}}

  // main loop //{{{
  FILE *fr = OpenFile(coor_file, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, start, silent);
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

      // allocate memory for temporary density arrays //{{{
      int ***temp_rho = malloc(Count->BeadType * sizeof(int **));
      for (int i = 0; i < Count->BeadType; i++) {
        temp_rho[i] = calloc(bins[0], sizeof(int *));
        for (int j = 0; j < bins[0]; j++) {
          temp_rho[i][j] = calloc(bins[1], sizeof(int));
        }
      } //}}}

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
          if (axis == 'x') {
            int j = bead->Position.y / width;
            int k = bead->Position.z / width;
            temp_rho[bead->Type][j][k]++;
          } else if (axis == 'y') {
            int j = bead->Position.x / width;
            int k = bead->Position.z / width;
            temp_rho[bead->Type][j][k]++;
          } else {
            int j = bead->Position.x / width;
            int k = bead->Position.y / width;
            temp_rho[bead->Type][j][k]++;
          }
        }
      } //}}}

      // add from temporary density array to global density arrays //{{{
      for (int j = 0; j < Count->BeadType; j++) {
        for (int k = 0; k < bins[0]; k++) {
          for (int l = 0; l < bins[1]; l++) {
            rho[j][k][l] += temp_rho[j][k][l];
            rho_2[j][k][l] += SQR(temp_rho[j][k][l]);
          }
        }
      } //}}}

      // free temporary density array //{{{
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = 0; j < bins[0]; j++) {
          free(temp_rho[i][j]);
        }
        free(temp_rho[i]);
      }
      free(temp_rho); //}}}
    } else {
      if (!SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count)) {
        count_coor--;
        break;
      }
    }

    if (end == count_coor) {
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
  FILE *fw = OpenFile(output_rho, "w");
  PrintByline(fw, argc, argv);

  // print bead type names to output file //{{{
  if (axis == 'x') {
    fprintf(fw, "# columns: (1) y; (2) z;");
  } else if (axis == 'y') {
    fprintf(fw, "# columns: (1) x; (2) z;");
  } else {
    fprintf(fw, "# columns: (1) x; (2) y;");
  }
  for (int i = 0; i < Count->BeadType; i++) {
    fprintf(fw, " (%d) %s", i+3, System.BeadType[i].Name);
    if (i != (Count->BeadType-1)) {
      putc(';', fw);
    }
  }
  putc('\n', fw);
//for (int i = 0; i < Count->BeadType; i++) {
//  fprintf(out, " %d: %s", 4*i+2, BeadType[i].Name);
//}
//fprintf(out, "\n# for each molecule type: rdp | stderr | rnp | stderr\n"); //}}}

  // calculate rdf
  double volume = SQR(width);
  if (axis == 'x') {
    volume *= System.Box.Length.x;
  } else if (axis == 'y') {
    volume *= System.Box.Length.y;
  } else {
    volume *= System.Box.Length.z;
  }
  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      fprintf(fw, " %7.3f", width * (i + 0.5)); // TODO: 2D, check Surface
      fprintf(fw, " %7.3f", width * (j + 0.5)); // TODO: 2D, check Surface
      for (int k = 0; k < Count->BeadType; k++) {
        double temp_rho = 0;
        temp_rho += rho[k][i][j] / (volume * count);
        // print average value to output file
        fprintf(fw, " %10f", temp_rho);
      }
      putc('\n',fw);
    }
    putc('\n',fw);
  }

  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = 0; j < bins[0]; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
