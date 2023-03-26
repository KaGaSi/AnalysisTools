#include "../AnalysisTools.h"

// TODO: -x option not implemented (needs change in Write.c, I guess)

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

  fprintf(ptr, "   <input>           input coordinate file\n");
  fprintf(ptr, "   <output.vcf>      output coordinate file (vcf format)\n");
  fprintf(ptr, "   <bead(s)>         names of bead types to save"
          " (optional if '--reverse' used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --reverse      reverse <bead name(s)>, i.e., exclude the"
          " specified bead types (use all if no <bead names> are present)\n");
  fprintf(ptr, "      --join         join molecules (remove pbc)\n");
  fprintf(ptr, "      --wrap         wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps"
          "(--last overrides this option)\n");
  // fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      --last         use only the last step"
          "(-st/-e/-n options are ignored)\n");
  int common = 11;
  char option[common][OPT_LENGTH];
  strcpy(option[ 0], "-st");
  strcpy(option[ 1], "-e");
  strcpy(option[ 2], "-sk");
  strcpy(option[ 3], "-i");
  strcpy(option[ 4], "--variable");
  strcpy(option[ 5], "--detailed");
  strcpy(option[ 6], "-pbc");
  strcpy(option[ 7], "-v");
  strcpy(option[ 8], "--silent");
  strcpy(option[ 9], "--help");
  strcpy(option[10], "--version");
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
    if (argv[i][0] == '-' && strcmp(argv[i], "--reverse") != 0 &&
        strcmp(argv[i], "--join") != 0 && strcmp(argv[i], "--wrap") != 0 &&
        strcmp(argv[i], "-n") != 0 && // strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "--last") != 0 && strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 && strcmp(argv[i], "-sk") != 0 &&
        strcmp(argv[i], "-i") != 0 && strcmp(argv[i], "--variable") != 0 &&
        strcmp(argv[i], "--detailed") != 0 && strcmp(argv[i], "-pbc") != 0 &&
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
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
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

  // options before reading system data
  bool silent, verbose, detailed, vtf_var;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1, ltrj_start_id = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &ltrj_start_id, &start, &end, &skip);
  bool join = BoolOption(argc, argv, "--join");
  bool wrap = BoolOption(argc, argv, "--wrap");
  bool last = BoolOption(argc, argv, "--last");

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz, &ltrj_start_id);

  // <bead names> - names of bead types to save //{{{
  bool *write = calloc(System.Count.Bead, sizeof *write),
       *write_bt = calloc(System.Count.BeadType, sizeof *write_bt);
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], System);
    if (type == -1) {
      strcpy(ERROR_MSG, "non-existent bead name");
      PrintErrorFile(struct_file, coor_file, "\0");
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
  } //}}}

  // '-n' option - specify timestep ids //{{{
  int n_opt_save[100] = {0}, n_opt_number = -1;
  if (IntegerOption(argc, argv, 100, "-n", &n_opt_number, n_opt_save)) {
    exit(1);
  }
  SortArray(n_opt_save, n_opt_number, 0); //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  // print initial stuff to output vcf file //{{{
  FILE *out = OpenFile(out_coor, "w");
  if (coor_out_type == VCF_FILE) {
    PrintByline(out, argc, argv);
  }
  fclose(out); //}}}

  // main loop //{{{
  FILE *fr = OpenFile(coor_file, "r");
  // file pointers for finding the last valid step
  fpos_t *position = calloc(1, sizeof *position);
  // save line count at every fgetpos()
  int *bkp_line_count = calloc(1, sizeof *bkp_line_count);
  int n_opt_count = 0,    // count saved steps if -n option is used
      count_coor = 0,     // count steps in the vcf file
      count_saved = 0,    // count steps in output file
      line_count = 0;     // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, start, silent);
    position = realloc(position, count_coor * sizeof *position);
    fgetpos(fr, &position[count_coor-1]);
    bkp_line_count = realloc(bkp_line_count, count_coor *
                             sizeof *bkp_line_count);
    bkp_line_count[count_coor-1] = line_count;
    // decide whether this timestep is to be saved //{{{
    bool use = false;
    /* no -n option - use if timestep
     *    1) is between start (-st option) and end (-e option)
     *    and
     *    2) isn't skipped (-sk option); skipping starts counting with 'start'
     */
    if (n_opt_number == -1) {
      // definitely not use, if --last option is used
      if (last) {
        use = false;
      } else if (count_coor >= start && (count_coor <= end || end == -1) && // 1)
                 ((count_coor - start) % skip) == 0) {                        // 2)
        use = true;
      } else {
        use = false;
      }
      // -n option is used - save the timestep if it's in the list
    } else if (n_opt_count < n_opt_number &&
               n_opt_save[n_opt_count] == count_coor) {
      use = true;
      n_opt_count++;
    } //}}}
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (!ReadTimestep(coor_type, fr, coor_file, &System, &line_count,
                        ltrj_start_id, vtf_var)) {
        count_coor--;
        break;
      }
      count_saved++;
      WrapJoinCoordinates(&System, wrap, join);
      WriteTimestep(coor_out_type, out_coor, System, count_coor, write);
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(coor_type, fr, coor_file,
                        struct_file, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // ecide whether to exit the main loop //{{{
    /* break the loop if
     *    1) all timesteps in the -n option are saved (and --last isn't used)
     *    or
     *    2) end timestep was reached (-e option)
     */
    if (!last && // never break when --last is used
        (n_opt_count == n_opt_number || // 1)
        count_coor == end)) { // 2)
      break;
    } //}}}
  } //}}}
  // if --last option is used, read & save the last timestep //{{{
  if (last) {
    /* To through all saved file positions (last to first) and save a the first
     * valid step encountered.
     * Start at count_coor as the saved position is at the beginning of the last
     * timestep to be skipped, and count_coor-- is used beore quitting the while
     * loop.
     */
    for (int i = (count_coor); i >= 0; i--) {
      fsetpos(fr, &position[i]);
      line_count = bkp_line_count[i];
      if (ReadTimestep(coor_type, fr, coor_file, &System, &line_count,
                       ltrj_start_id, vtf_var)) {
        count_saved++;
        WrapJoinCoordinates(&System, wrap, join);
        WriteTimestep(coor_out_type, out_coor, System, count_coor, write);
        if (!silent) {
          if (isatty(STDOUT_FILENO)) {
            fprintf(stdout, "\r                          \r");
          }
          fprintf(stdout, "Saved Step: %d\n", i+1);
          fflush(stdout);
        }
        break;
      } else {
        snprintf(ERROR_MSG, LINE, "disregarding step %s%d%s", ErrYellow(),
                 i+1, ErrCyan());
        PrintWarnFile(coor_file, "\0", "\0");
      }
    } //}}}
  } else if (count_coor == 0) { // error - input file without a valid timestep //{{{
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(coor_file, "\0", "\0");
    fputc('\n', stderr); //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG, "no coordinates written (starting timestep higher"
           " than the total number of timesteps)");
    PrintWarning(); //}}}
  } else if (!silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (saved %d)\n", count_coor, count_saved);
    fflush(stdout);
  } //}}}
  fclose(fr);

  // free memory
  FreeSystem(&System);
  free(write);
  free(position);
  free(bkp_line_count);

  // // reading file from the end //{{{
  // FILE *fr = OpenFile(in_coor, "r");
  // fseek(fr, -1, SEEK_END);
  // char test = fgetc(fr);
  // if (test == '\n') {
  //   printf("empty line\n");
  //   fseek(fr, -2, SEEK_END);
  // }
  // test = fgetc(fr);
  // printf("%c", test);
  // do {
  //   fseek(fr, -2, SEEK_CUR);
  //   test = fgetc(fr);
  //   printf("%c", test);
  // } while (test != '\n');
  // fgets(in_coor, LINE, fr);
  // printf("%s", in_coor);
  // fclose(fr); //}}}

  return 0;
}
