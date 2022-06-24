#include "../AnalysisTools.h"

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
  fprintf(ptr, "   %s <input> <output.vcf> <bead(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>        input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <output.vcf>   output coordinate file (vcf format)\n");
  fprintf(ptr, "   <bead(s)>      names of bead types to save \
(optional if '--reverse' used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --reverse      reverse <bead name(s)>, i.e., exclude \
the specified bead types (use all if no <bead names> are present)\n");
  fprintf(ptr, "      --join         join molecules (remove pbc)\n");
  fprintf(ptr, "      --wrap         wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -sk <skip>     leave out every 'skip' steps\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps \
(if --last option is used, save also the last timestep)\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      -xyz <name>    output xyz file\n");
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
  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }
  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "--reverse");
  // possible to omit <bead name(s)> if '--reverse' is used
  if (count < (req_args-1) || (count == (req_args-1) && !reverse)) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "--reverse") != 0 &&
        strcmp(argv[i], "--join") != 0 &&
        strcmp(argv[i], "--wrap") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-xyz") != 0 &&
        strcmp(argv[i], "--last") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE] = "", input_vsf[LINE] = "";
  snprintf(input_coor, LINE, "%s", argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output.vcf> - output vcf file //{{{
  char output_vcf[LINE] = "";
  snprintf(output_vcf, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vcf'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, input_vsf, LINE, &verbose, &silent, &detailed);
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
    exit(1);
  }
  skip++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  // should output coordinates be joined?
  bool join = BoolOption(argc, argv, "--join");
  // should output coordinates be wrapped?
  bool wrap = BoolOption(argc, argv, "--wrap");
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // save into xyz file?
  char output_xyz[LINE] = "";
  if (FileOption(argc, argv, "-xyz", output_xyz, LINE)) {
    exit(1);
  }
  if (output_xyz[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".xyz");
    if (ErrorExtension(output_xyz, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  }
  // use only the last step?
  bool last = BoolOption(argc, argv, "--last");
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  SYSTEM System = VtfReadStruct(input_vsf, detailed);
  VtfReadPBC(input_coor, input_vsf, &System.Box);
  if (!TriclinicCellData(&System.Box)) {
    strcpy(ERROR_MSG, "wrong pbc data");
    PrintError();
    exit(1);
  } //}}}

  // <bead names> - names of bead types to save //{{{
  bool *write = calloc(System.Count.Bead, sizeof *write),
       *write_bt = calloc(System.Count.BeadType, sizeof *write_bt);
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], System);
    if (type == -1) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", input_coor);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - non-existent bead name ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s\n", argv[count]);
      ColourReset(STDERR_FILENO);
      ErrorBeadType(System);
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

  // '-n' option - specify bead ids //{{{
  int n_opt_save[100] = {0}, n_opt_number = -1;
  if (MultiIntegerOption(argc, argv, "-n", &n_opt_number, n_opt_save)) {
    exit(1);
  }
  SortArray(n_opt_save, n_opt_number, 0); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  fpos_t position1, position2;

  // print initial stuff to output vcf file
  FILE *out = OpenFile(output_vcf, "w");
  PrintByline(out, argc, argv);
  fclose(out);
  // make sure a new xyz file is created (if -xyz option is used)
  if (output_xyz[0] != '\0') {
    fclose(OpenFile(output_xyz, "w"));
  }

  // main loop //{{{
  int n_opt_count = 0, // count saved steps if -n option is used
      count_vcf = 0, // count steps in the vcf file
      file_line_count = 0; // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count_vcf++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (last) {
        fprintf(stdout, "\rDiscarding step: %d", count_vcf);
      } else {
        if (count_vcf == start) {
          fprintf(stdout, "\rStarting step: %d\n", start);
        }
        fprintf(stdout, "\rStep: %d", count_vcf);
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
      if ((count_vcf >= start && (count_vcf <= end || end == -1)) && // 1)
          ((count_vcf-start)%skip) == 0) { // 2)
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
               n_opt_save[n_opt_count] == count_vcf) {
      use = true;
      n_opt_count++;
    }
    //}}}
    // read and write the timestep, if it should be saved //{{{
    if (use) {
      if (!VtfReadTimestep(vcf, input_coor, input_vsf, &System,
                           &file_line_count, count_vcf, stuff)) {
        count_vcf--;
        break;
      }
      // transform coordinates into fractional ones for non-orthogonal box
      if (wrap || join) {
        ToFractionalCoor(&System);
      }
      if (wrap) { // wrap coordinates into the simulation box
        RestorePBC(&System);
      }
      if (join) { // join molecules by removing periodic boundary conditions
        RemovePBCMolecules(&System);
      }
      // transform back to 'normal' coordinates for non-orthogonal box
      if (wrap || join) {
        FromFractionalCoor(&System);
      }
      // write to output .vcf file
      out = OpenFile(output_vcf, "a");
      VtfWriteCoorIndexed(out, stuff, write, System);
      fclose(out);
      // write to xyz file?
      if (output_xyz[0] != '\0') {
        out = OpenFile(output_xyz, "a");
        XyzWriteCoor(out, write, System);
        fclose(out);
      }
      //}}}
    // skip the timestep, if it shouldn't be saved //{{{
    } else {
      if (!VtfSkipTimestep(vcf, input_coor, input_vsf, &file_line_count,
                           count_vcf)) {
        count_vcf--;
        break;
      }
    } //}}}
    // save file position (last two because of --last) //{{{
    if ((count_vcf%2) == 0) {
      fgetpos(vcf, &position1);
    } else {
      fgetpos(vcf, &position2);
    } //}}}
    // decide whether to exit the main loop //{{{
    /* break the loop if
     *    1) all timesteps in the -n option are saved (and --last isn't used)
     *    or
     *    2) end timestep was reached (-e option)
     */
    if ((n_opt_count == n_opt_number && !last) || // 1)
        count_vcf == end) { // 2)
      break;
    } //}}}
  } //}}}
  // if --last option is used, read & save the last timestep //{{{
  if (last) {
    if ((count_vcf%2) == 1) {
      fsetpos(vcf, &position1);
    } else {
      fsetpos(vcf, &position2);
    }
    VtfReadTimestep(vcf, input_coor, input_vsf, &System,
                    &file_line_count, count_vcf, stuff);
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(&System);
    // wrap and/or join molecules?
    if (wrap) {
      RestorePBC(&System);
    }
    if (join) {
      RemovePBCMolecules(&System);
    }
    // transform back to 'normal' coordinates for non-orthogonal box
    FromFractionalCoor(&System);
    // write to output .vcf file
    out = OpenFile(output_vcf, "a");
    VtfWriteCoorIndexed(out, stuff, write, System);
    fclose(out);
    // write to xyz file?
    if (output_xyz[0] != '\0') {
      out = OpenFile(output_xyz, "a");
      XyzWriteCoor(out, write, System);
      fclose(out);
    }
  //}}}
  // warn if no timesteps were actually written //{{{
  } else if (start > count_vcf) {
    strcpy(ERROR_MSG, "no coordinates written (starting timestep higher \
than the number of timestep)");
    PrintWarning();
  } //}}}
  fclose(vcf);
  // print last step count? //{{{
  // error - input coordinate file with no coordinates //{{{
  if (count_vcf == 0) {
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(input_coor, "\0");
    fputc('\n', stderr);
  //}}}
  } else if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
    fflush(stdout);
  } //}}}

  // free memory
  FreeSystem(&System);
  free(stuff);
  free(write);

  return 0;
}
