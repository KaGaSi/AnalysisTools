#include "../AnalysisTools.h"
// TODO: pbc_xyz in ReadTimestep

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
SelectedVcf creates new <output.vcf> file (and possibly xyz file) from \
<input> containing only selected bead types. Periodic boundary conditions \
can be either stripped away or applied (which happens first if both \
'--join' and '--wrap' options are used). Also, specified molecules can be \
excluded. However, AnalysisTools utilities can only read coordinate files \
containing all beads of any given type, so the usefulness is very limited \
(for, e.g., visualization using vmd).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <output> <bead(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file \
(vcf, vtf or xyz format)\n");
  fprintf(ptr, "   <output.vcf>      output coordinate file (vcf format)\n");
  fprintf(ptr, "   <bead(s)>         names of bead types to save"
          " (optional if '--reverse' used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --reverse      reverse <bead name(s)>, i.e., exclude the"
          " specified bead types (use all if no <bead names> are present)\n");
  fprintf(ptr, "      --join         join molecules (remove pbc)\n");
  fprintf(ptr, "      --wrap         wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -sk <skip>     leave out every 'skip' steps\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps"
          "(if --last option is used, save also the last timestep)\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      --last         use only the last step"
          "(-st/-e options are ignored; -n option is not)\n");
  fprintf(ptr, "      -pbc <int>         position of pbc in xyz file's comment"
          " line (of the first number)\n");
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
  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
    count++;
  }
  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "--reverse");
  // possible to omit <bead name(s)> if '--reverse' is used
  if (count < (req_args - 1) || (count == (req_args - 1) && !reverse)) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && strcmp(argv[i], "-l_in") != 0 &&
        strcmp(argv[i], "-i") != 0 && strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--detailed") != 0 && strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "--silent") != 0 && strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--reverse") != 0 && strcmp(argv[i], "--join") != 0 &&
        strcmp(argv[i], "--wrap") != 0 && strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 && strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-x") != 0 && strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-pbc") != 0 && strcmp(argv[i], "--last") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char in_coor[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(in_coor, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, in_coor, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // <output> - output vcf file //{{{
  char out_coor[LINE] = "";
  snprintf(out_coor, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vcf'
  int ext = 3;
  char extension[3][EXTENSION];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".xyz");
  strcpy(extension[2], ".lammpstrj");
  int coor_out_type = ErrorExtension(out_coor, ext, extension);
  if (coor_out_type == 0) {
    coor_out_type = VCF_FILE;
  } else if (coor_out_type == 1) {
    coor_out_type = XYZ_FILE;
  } else if (coor_out_type == 2) {
    coor_out_type = LTRJ_FILE;
  } else if (coor_out_type == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed);
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
    exit(1);
  }
  skip++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  bool join = BoolOption(argc, argv, "--join");
  bool wrap = BoolOption(argc, argv, "--wrap");
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  bool last = BoolOption(argc, argv, "--last");
  // position of the first number of pbc in xyz file
  int pbc_xyz = -1;
  if (IntegerOption(argc, argv, "-pbc", &pbc_xyz)) {
    exit(1);
  }
  if (pbc_xyz == 0) {
    strcpy(ERROR_MSG, "position must be a positive number");
    PrintErrorOption("-pbc");
  } //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, detailed, pbc_xyz);
  bool variable_coor = false;
  // read number of beads in a non-variable coordinate file
  if (!variable_coor) {
    printf("%d ", System.Count.BeadCoor);
    System.Count.BeadCoor = CoorReadNumberOfBeads(coor_type, in_coor);
    printf("%d\n", System.Count.BeadCoor);
  }
  // in case of lammpstrj file, find if atom ids start at 0
  int start_id = -1;
  if (coor_type == LTRJ_FILE) {
    start_id = LtrjLowIndex(in_coor);
  }

  // <bead names> - names of bead types to save //{{{
  bool *write = calloc(System.Count.Bead, sizeof *write),
       *write_bt = calloc(System.Count.BeadType, sizeof *write_bt);
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], System);
    if (type == -1) {
      strcpy(ERROR_MSG, "non-existent bead name");
      PrintErrorFile(struct_file, in_coor, "\0");
      ErrorBeadType(argv[count], System);
      exit(1);
    }
    write_bt[type] = true;
  }
  // if '--reverse' is used, switch Write bools for all bead types
  if (reverse) {
    for (int i = 0; i < System.Count.BeadType; i++) {
      if (write_bt[i]) {
        write_bt[i] = false;
      } else {
        write_bt[i] = true;
      }
    }
  }
  for (int i = 0; i < System.Count.Bead; i++) {
    int type = System.Bead[i].Type;
    if (write_bt[type]) {
      write[i] = true;
    }
  }
  free(write_bt); //}}}

  // '-x' option //{{{
  // TODO remove the bool flags from SYSTEM, i.e., also change ExcludeOption()
  if (ExcludeOption(argc, argv, &System)) {
    exit(1);
  }
  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    System.MoleculeType[i].Write = System.MoleculeType[i].Use;
  } //}}}

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
  FILE *out = OpenFile(out_coor, "w");
  if (coor_out_type == VCF_FILE) {
    PrintByline(out, argc, argv);
  }
  fclose(out); //}}}

  // main loop //{{{
  FILE *coor = OpenFile(in_coor, "r");
  fpos_t position1, position2; // two file pointers for finding the last step
  int n_opt_count = 0,         // count saved steps if -n option is used
      count_coor = 0,          // count steps in the vcf file
      count_saved = 0,         // count steps in output file
      file_line_count = 0;     // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
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
    // decide whether this timestep is to be saved //{{{
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
    } //}}}
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (!ReadTimestep(coor_type, coor, in_coor, &System, &file_line_count,
                        start_id, stuff)) {
        count_coor--;
        break;
      }
      count_saved++;
      WrapJoinCoordinates(&System, wrap, join);
      WriteTimestep(coor_out_type, out_coor, System, count_coor, stuff, write);
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(coor_type, coor, in_coor,
                        struct_file, &file_line_count)) {
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
  // if --last option is used, read & save the last timestep //{{{
  if (last) {
    // restore file pointer
    if ((count_coor % 2) == 1) {
      fsetpos(coor, &position1);
    } else {
      fsetpos(coor, &position2);
    }
    ReadTimestep(coor_type, coor, in_coor, &System, &file_line_count,
                 start_id, stuff);
    count_saved++;
    WrapJoinCoordinates(&System, wrap, join);
    WriteTimestep(coor_out_type, out_coor, System, count_coor, stuff, write);
  } //}}}
  fclose(coor);
  if (count_coor == 0) { // error - input file without a valid timestep //{{{
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(in_coor, "\0", "\0");
    fputc('\n', stderr); //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG, "no coordinates written (starting timestep higher"
           " than the number of timestep)");
    PrintWarning(); //}}}
  } else if (!silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (saved %d)\n", count_coor, count_saved);
    fflush(stdout);
  } //}}}

  // free memory
  FreeSystem(&System);
  free(stuff);
  free(write);

  return 0;
}
