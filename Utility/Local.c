#include "../AnalysisTools.h"

// TODO: -e X --last leads to leaving at step X and saving that as the last;
//       is that okay?
// TODO: make some in_struct akin to in_coor (in InputCoorStruct() function)

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "TEXT TO BE ADDED\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <axis> <with> <output> [options]\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file"
          " (vcf, vtf, lammpstrj, or xyz format)\n");
  fprintf(ptr, "   <dimensions>      1 or 2 axis for 1D or 2D calculation"
          " (e.g., x or xy)\n");
  fprintf(ptr, "   <width>           width of a single bin"
          " (in the 1 or both dimensions)\n");
  fprintf(ptr, "   <output>          output file with local variables\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -sk <skip>     leave out every 'skip' steps\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps \
(if --last option is used, save also the last timestep)\n");
  fprintf(ptr, "      --last         use only the last step \
(-st/-e options are ignored; -n option is not)\n");
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
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
    count++;
  }
  // possible to omit <bead name(s)> if '--reverse' is used
  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && strcmp(argv[i], "-l_in") != 0 &&
        strcmp(argv[i], "-v") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "--silent") != 0 && strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 && strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 && strcmp(argv[i], "-sk") != 0 &&
        strcmp(argv[i], "-n") != 0 && strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "--last") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char in_coor[LINE] = "", in_vsf[LINE] = "", in_lmp[LINE] = "",
       in_field[LINE] = "", *struct_file;
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
  if (struct_type == 0 && coor_type != 2 && coor_type != 3) {
    strcpy(ERROR_MSG, "missing input structure file; \
acceptable only for xyz input coordinate file");
    PrintError();
    exit(1);
  } //}}}

  // <axis> //{{{
  count++;
  bool err = false;
  int dim[3] = {-1, -1, -1}; // [0]..1/2 for 1D/2D, [0..1]..0, 1, 2 for x, y, z
  if (argv[count][0] == 'x') {
    dim[1] = 0;
  } else if (argv[count][0] == 'y') {
    dim[1] = 1;
  } else if (argv[count][0] == 'z') {
    dim[1] = 2;
  } else {
    err = true;
  }
  if (argv[count][1] == '\0') {
    dim[0] = 1;
  } else {
    dim[0] = 2;
    if (argv[count][1] == 'x') {
      dim[2] = 0;
    } else if (argv[count][1] == 'y') {
      dim[2] = 1;
    } else if (argv[count][1] == 'z') {
      dim[2] = 2;
  } else {
    err = true;
    }
  }
  if (dim[1] == dim[2]) {
    err = true;
  }
  if (err) {
    strcpy(ERROR_MSG, "<axis> argument must contain axis labels (x, y, z)");
    PrintError();
    Help(argv[0], true);
  } //}}}

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  double width[2];
  if (!IsPosRealNumber(argv[++count], &width[0])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  width[1] = width[0]; // TODO: for now
  //}}}

  // <output> - output text file
  char out_file[LINE] = "";
  snprintf(out_file, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed);
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
    exit(1);
  }
  skip++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  // should output coordinates be joined?
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // use only the last step?
  bool last = BoolOption(argc, argv, "--last");
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // TODO: make reading various structure types into a function in Read.c
  // read input data //{{{
  SYSTEM System;
  switch (struct_type) {
  case 0: // xyz or lammpstrj
    if (coor_type == 2) {
      System = XyzReadStruct(in_coor);
    } else {
      System = LtrjReadStruct(in_coor);
    }
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
  WarnChargedSystem(System, in_vsf, in_field, in_lmp);
  // warn if missing box dimensions
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "unspecified box dimensions");
    PrintWarning();
    WarnPrintFile(struct_file, in_coor, "\0");
    putc('\n', stderr);
  } //}}}

  // number of bins //{{{
  // TODO: variable box size? For now, assume at most twice the size
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    // PrintErrorFile(); // TODO add file(s)
    exit(1);
  }
  BOX *box = &System.Box;
  int bin[2] = {-1, 1};
  for (int i = 0; i <= dim[0]; i++) {
    if (dim[i+1] == 0) {
      bin[i] = 2 * box->Length.x / width[i] + 1;
    } else if (dim[i+1] == 1) {
      bin[i] = 2 * box->Length.y / width[i] + 1;
    } else if (dim[i+1] == 2) {
      bin[i] = 2 * box->Length.z / width[i] + 1;
    }
  }
  if (bin[1] == 1) { // 1D case
    width[1] = 2 * Max3(box->Length.x, box->Length.y, box->Length.z) + 1;
  } //}}}

// TODO: warn 0 mass/charge/radius/etc.
  // '-n' option - specify timestep ids //{{{
  int n_opt_save[100] = {0}, n_opt_number = -1;
  if (MultiIntegerOption(argc, argv, "-n", &n_opt_number, n_opt_save)) {
    exit(1);
  }
  SortArray(n_opt_save, n_opt_number, 0); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  // print initial stuff to output vcf file //{{{
  FILE *out = OpenFile(out_file, "w");
  PrintByline(out, argc, argv);
  fclose(out); //}}}

  // main loop //{{{
  FILE *coor = OpenFile(in_coor, "r");
  fpos_t position1, position2; // two file pointers for finding the last step
  int n_opt_count = 0,         // count saved steps if -n option is used
      count_coor = 0,          // count steps in the vcf file
      count_used = 0,         // count steps in output file
      file_line_count = 0;     // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  // when 1D is calculated, only Temp[i][0] is used
  long **bead_count = calloc(bin[0], sizeof **bead_count); // count beads
  double **Temp = calloc(bin[0], sizeof **Temp); // temperature
  double **Ekin = calloc(bin[0], sizeof **Ekin); // kinetic energy
  for (int i = 0; i < bin[0]; i++) {
    bead_count[i] = calloc(bin[1], sizeof *bead_count);
    Temp[i] = calloc(bin[1], sizeof *Temp);
    Ekin[i] = calloc(bin[1], sizeof *Ekin);
  }
  FILE *t_out = OpenFile("t.txt", "w");
  while (true) {
    count_coor++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (last) {
        fprintf(stdout, "\rDiscarding step: %d", count_coor);
      } else {
        if (count_coor == start) {
          fprintf(stdout, "\rStarting step: %d\n", start);
        }
        fprintf(stdout, "\rStep: %d", count_coor);
      }
      fflush(stdout);
    } //}}}
    // decide whether this timestep is to be used //{{{
    bool use = false;
    /* no -n option - use if timestep
     *    1) is between start (-st option) and end (-e option)
     *    and
     *    2) isn't skipped (-sk option); skipping starts counting with 'start'
     */
    if (n_opt_number == -1) {
      if ((count_coor >= start && (count_coor <= end || end == -1)) && // 1)
          ((count_coor - start) % skip) == 0) {                        // 2)
        use = true;
      } else {
        use = false;
      }
      // definitely not use, if --last option is used
      if (last) {
        use = false;
      }
      // -n option is used - save the timestep if it's in the list
    } else if (n_opt_count < n_opt_number &&
               n_opt_save[n_opt_count] == count_coor) {
      use = true;
      n_opt_count++;
    }
    //}}}
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (!ReadTimestep(coor_type, coor, in_coor, &System, &file_line_count,
                        stuff)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false); // restore pbc
      // calculate the local properties
      InitLong2DArray(bead_count, bin[0], bin[1], 0);
      InitDouble2DArray(Temp, bin[0], bin[1], 0);
      InitDouble2DArray(Ekin, bin[0], bin[1], 0);
      double temperature = 0;
      for (int i = 0; i < System.Count.BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        int pos[2] = {0, 0};
        for (int j = 0; j < 2; j++) {
          if (dim[j+1] == 0) {
            pos[j] = b->Position.x / width[j];
          } else if (dim[j+1] == 1) {
            pos[j] = b->Position.y / width[j];
          } else if (dim[j+1] == 2) {
            pos[j] = b->Position.z / width[j];
          }
        }
        bead_count[pos[0]][pos[1]]++;
        double vel2 = SQR(b->Velocity.x) +
                      SQR(b->Velocity.y) +
                      SQR(b->Velocity.z);
        Temp[pos[0]][pos[1]] += bt->Mass * vel2;
        Ekin[pos[0]][pos[1]] += 0.5 * bt->Mass * vel2;
        temperature += bt->Mass * vel2;
      }
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          Temp[i][j] /= 3 * bead_count[i][j];
          // Ekin[i][j] /= bead_count[i][j];
        }
      }
      temperature /= 3 * System.Count.BeadCoor;
      fprintf(t_out, "%8d %8.4f\n", count_coor, temperature);
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(coor_type, coor, in_coor, in_vsf, &file_line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // save file position (last two because of --last) //{{{
    if ((count_coor % 2) == 0) {
      fgetpos(coor, &position1);
    } else {
      fgetpos(coor, &position2);
    } //}}}
    // decide whether to exit the main loop //{{{
    /* break the loop if
     *    1) all timesteps in the -n option are saved (and --last isn't used)
     *    or
     *    2) end timestep was reached (-e option)
     */
    if ((n_opt_count == n_opt_number && !last) || // 1)
        count_coor == end) {                      // 2)
      break;
    } //}}}
  }   //}}}
  fclose(t_out);
  // if --last option is used, read & save the last timestep //{{{
  if (last) {
    // restore file pointer
    if ((count_coor % 2) == 1) {
      fsetpos(coor, &position1);
    } else {
      fsetpos(coor, &position2);
    }
    ReadTimestep(coor_type, coor, in_coor, &System, &file_line_count, stuff);
    count_used++;
    WrapJoinCoordinates(&System, true, false); // restore pbc
  } //}}}
  fclose(coor);
  if (count_coor == 0) { // error - input file without a valid timestep //{{{
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(in_coor, "\0", "\0");
    fputc('\n', stderr); //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG, "no coordinates written (starting timestep higher \
than the number of timestep)");
    PrintWarning(); //}}}
  } else if (!silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
    fflush(stdout);
  } //}}}

  // free memory
  for (int i = 0; i < bin[1]; i++) {
    free(bead_count[i]);
    free(Temp[i]);
  }
  FreeSystem(&System);
  free(stuff);

  return 0;
}
