#include "../AnalysisTools.h"

// TODO: -x option not implemented (needs change in Write.c, I guess)

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Selected creates new coordinate file in the extension-specified format \
that contains specified bead and/or molecule types. Periodic boundary \
conditions can be either stripped away or applied (which happens first if both \
'--join' and '--wrap' options are used).\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <output> <bead(s)> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output.vcf>        output coordinate file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -bt <bead type>   bead types to exclude\n");
  fprintf(ptr, "  -mt <mol type>    molecule types to exclude\n");
  fprintf(ptr, "  --reverse         reverse <bead name(s)>, i.e., save the"
          " specified bead types (-bt and/or -mt is needed)\n");
  fprintf(ptr, "  --join            join molecules (remove pbc)\n");
  fprintf(ptr, "  --wrap            wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "  -n <int(s)>       save only specified timesteps"
          "(--last overrides this option)\n");
  fprintf(ptr, "  --last            use only the last step"
          "(-st/-e/-n options are ignored)\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 10, all = common + 7, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "-pbc");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "--reverse");
  strcpy(option[count++], "--join");
  strcpy(option[count++], "--wrap");
  strcpy(option[count++], "-n");
  strcpy(option[count++], "--last");
  strcpy(option[count++], "-bt");
  strcpy(option[count++], "-mt");
  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate (and structure) file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // <output> - output coordinate file
  char coor_out_file[LINE] = "";
  snprintf(coor_out_file, LINE, "%s", argv[++count]);
  int coor_out_type = CoordinateFileType(coor_out_file, 1);

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &start, &end, &skip);
  bool reverse = BoolOption(argc, argv, "--reverse");
  bool join = BoolOption(argc, argv, "--join");
  bool wrap = BoolOption(argc, argv, "--wrap");
  bool last = BoolOption(argc, argv, "--last"); //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file,
                                coor_type, coor_file, detailed, pbc_xyz);

  // specify beads to save (possibly using -bt and/or -mt options) //{{{
  // -bt option - which bead types to exclude/use
  bool *write_bt = malloc(System.Count.BeadType * sizeof *write_bt);
  InitBoolArray(write_bt, System.Count.BeadType, true);
  if (reverse) { // save only specified bead types
    InitBoolArray(write_bt, System.Count.BeadType, false);
    BeadTypeOption(argc, argv, "-bt", true, write_bt, System);
  } else { // exclude specifed bead types
    InitBoolArray(write_bt, System.Count.BeadType, true);
    BeadTypeOption(argc, argv, "-bt", false, write_bt, System);
  }
  // -mt option - which molecule types to exclude/use
  bool *write_mt = malloc(System.Count.MoleculeType * sizeof *write_mt);
  InitBoolArray(write_mt, System.Count.MoleculeType, true);
  if (reverse) { // save only specified molecule types
    InitBoolArray(write_mt, System.Count.MoleculeType, false);
    MoleculeTypeOption(argc, argv, "-mt", true, write_mt, System);
  } else { // exclude specifed molecule types
    InitBoolArray(write_mt, System.Count.MoleculeType, true);
    MoleculeTypeOption(argc, argv, "-mt", false, write_mt, System);
  }
  // specify beads to save
  bool *write = malloc(System.Count.Bead * sizeof *write);
  InitBoolArray(write, System.Count.Bead, false);
  for (int i = 0; i < System.Count.Bead; i++) {
    int type = System.Bead[i].Type;
    if (write_bt[type]) {
      write[i] = true;
    }
  }
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    for (int j = 0; j < System.MoleculeType[i].Number; j++) {
      int mol = System.MoleculeType[i].Index[j];
      for (int k = 0; k < System.MoleculeType[i].nBeads; k++) {
        int id = System.Molecule[mol].Bead[k];
        if (write_mt[i]) {
          write[id] = true;
        } else {
          write[id] = false;
        }
      }
    }
  }
  free(write_mt);
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
  // ignore -st/-e/-sk when -n is used
  if (n_opt_number != -1) {
    start = 1;
    end = -1;
    skip = 1;
  }
  SortArray(n_opt_save, n_opt_number, 0); //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  // print initial stuff to output coordinate file //{{{
  if (coor_out_type == VCF_FILE) {
    FILE *out = OpenFile(coor_out_file, "w");
    PrintByline(out, argc, argv);
    fclose(out);
  } else if (coor_out_type == VTF_FILE) {
    WriteStructure(VSF_FILE, coor_out_file, System, -1, false);
    coor_out_type = VCF_FILE;
  } else {
    FILE *out = OpenFile(coor_out_file, "w");
    fclose(out);
  } //}}}

  // for lammps data as a coordinate file, only the one step is used
  if (coor_type == LDATA_FILE) {
    start = 1;
    skip = 1;
    end = 1;
  }

  FILE *fr = OpenFile(coor_file, "r");
  // main loop //{{{
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
      } else if (count_coor >= start &&              // 1)
                 (count_coor <= end || end == -1) && //
                 ((count_coor - start) % skip) == 0) { // 2)
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
      if (!ReadTimestep(coor_type, fr, coor_file, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_saved++;
      WrapJoinCoordinates(&System, wrap, join);
      WriteTimestep(coor_out_type, coor_out_file, System, count_coor, write);
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
      if (ReadTimestep(coor_type, fr, coor_file, &System, &line_count)) {
        count_saved++;
        WrapJoinCoordinates(&System, wrap, join);
        WriteTimestep(coor_out_type, coor_out_file, System, count_coor, write);
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
    PrintErrorFile(coor_file, "\0", "\0"); //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG, "no coordinates written (starting timestep is higher"
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
