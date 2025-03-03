#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
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

  fprintf(ptr, "Usage: %s <input> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output coordinate file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -bt <bead type>   bead types to exclude\n");
  fprintf(ptr, "  -mt <mol type>    molecule types to exclude\n");
  fprintf(ptr, "  --keep            save only the specified types "
          " instead of excluding them");
  fprintf(ptr, "  --join            join molecules (remove pbc)\n");
  fprintf(ptr, "  --wrap            wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "  -n <int(s)>       save only specified timesteps"
          "(--last overrides this option)\n");
  fprintf(ptr, "  --last            use only the last step"
          "(-st/-e/-n options are ignored)\n");
  fprintf(ptr, "  -sc <float>       divide all coordinates by given value\n");
  fprintf(ptr, "  -m <3xfloat>      move all coordinates by given vector "
          "(-sc option is applied first)\n");
  fprintf(ptr, "  -cx 2×<float>     constrain x-coordinate to specified "
          "dimensions (in fraction of output box); multiple pairs possible\n");
  fprintf(ptr, "  -cy 2×<float>     constrain y-coordinate to specified "
          "dimensions (in fraction of output box); multiple pairs possible\n");
  fprintf(ptr, "  -cz 2×<float>     constrain z-coordinate to specified "
          "dimensions (in fraction of output box); multiple pairs possible\n");
  fprintf(ptr, "  --real            use real coordinates for "
          "-cx/-cy/-cz options instead of box fractions\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool bt, mt,             // -bt/-mt
       reverse,              // --reverse
       join, wrap, last,   // --join --wrap --last
       real;                 // --real
  int n_save[100], n_number; // -n
  double scale,              // -sc
         move[3];            // -m
  double ca[3][100];         // -cx/y/z ... slice(s)' coordinates
  int ca_count[3];           // -cx/y/z ... number of slices per axis
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

static void ScaleCoordinates(SYSTEM *System, const double scale) { //{{{
  if (scale != 1) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[id].Position[dd] /= scale;
      }
    }
  }
} //}}}
static void MoveCoordinates(SYSTEM *System, const double move[3]) { //{{{
  if (move[0] != 0 || move[1] != 0 || move[2] != 0) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[id].Position[dd] += move[dd];
      }
    }
  }
} //}}}
static void CopyWrite(const int num, bool *new, bool *old) { //{{{
  for (int i = 0; i < num; i++) {
    new[i] = old[i];
  }
} //}}}
// ConstrainCoordinates() //{{{
static void ConstrainCoordinates(SYSTEM *System, const OPT opt,
                                 bool *write_new, bool **write_orig) {
  if (opt.ca_count[0] == 0 && opt.ca_count[1] == 0 && opt.ca_count[2] == 0) {
    return;
  }
  // recalculate constraints if --real was not used
  double con[3][100];
  bool init[3] = {true, true, true};
  for (int dd = 0; dd < 3; dd++) {
    for (int i = 0; i < opt.ca_count[dd]; i++) {
      con[dd][i] = opt.ca[dd][i];
      init[dd] = false;
      if (!opt.real) {
        con[dd][i] *= System->Box.Length[dd];
      }
    }
  }
  COUNT *Count = &System->Count;
  // allocate array to store the original 'write'
  *write_orig = malloc(Count->Bead * sizeof *write_orig);
  // save the original 'write' array
  CopyWrite(Count->Bead, *write_orig, write_new);
  // go over beads in the coordinate file
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System->BeadCoor[i];
    // skip constraints for beads not to be written
    if (!(*write_orig)[id]) {
      continue;
    }
    // check -cx/-cy/-cz constraint
    double (*pos)[3] = &System->Bead[id].Position;
    bool save[3] = {init[0], init[1], init[2]};
    for (int dd = 0; dd < 3; dd++) {
      for (int j = 0; j < opt.ca_count[dd]; j+=2) {
        if ((*pos)[dd] >= con[dd][j] && (*pos)[dd] <= con[dd][j+1]) {
          save[dd] = true;
          break;
        }
      }
    }
    // printf("%d %d %d\n", save[0], save[1], save[2]);
    // if at least one axis constraint isn't met, don't save the bead
    if (!save[0] || !save[1] || !save[2]) {
      write_new[id] = false;
    }
  }
} //}}}
// make all the alterations and saves into one function
void TransformAndSave(SYSTEM *System, const OPT opt, FILE_TYPE f, bool write[],
                      const int count_coor, const int argc, char *argv[]) {
  bool *write2 = NULL; // used in case of -cx/-cy/-cz constraints
  ConstrainCoordinates(System, opt, write, &write2);
  ScaleCoordinates(System, opt.scale);
  MoveCoordinates(System, opt.move);
  WrapJoinCoordinates(System, opt.wrap, opt.join);
  WriteTimestep(f, *System, count_coor, write, argc, argv);
  if (opt.ca_count[0] > 0 ||
      opt.ca_count[1] > 0 ||
      opt.ca_count[2] > 0) {
    // restore the original 'write' array
    CopyWrite(System->Count.Bead, write, write2);
    // free the array alloc'd in ConstrainCoordinates()
    free(write2);
  }
}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 13, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version", "-bt", "-mt", "--keep", "--join", "--wrap", "-n",
              "--last", "-sc", "-m", "-cx", "-cy", "-cz", "--real");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <output> - output coordinate file
  FILE_TYPE fout;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  opt->reverse = BoolOption(argc, argv, "--reverse");
  opt->join = BoolOption(argc, argv, "--join");
  opt->wrap = BoolOption(argc, argv, "--wrap");
  opt->last = BoolOption(argc, argv, "--last");
  if (!OneNumberOption(argc, argv, "-sc", &opt->scale, 'd')) {
    opt->scale = 1;
  }
  opt->real = BoolOption(argc, argv, "--real");
  if (!ThreeNumbersOption(argc, argv, "-m", opt->move, 'd')) {
    for (int dd = 0; dd < 3; dd++) {
      opt->move[dd] = 0;
    }
  } else if (opt->real) {

  }
  // constraints (-cx/-cy/-cz) //{{{
  for (int dd = 0; dd < 3; dd++) {
    char option[10];
    if (dd == 0) {
      s_strcpy(option, "-cx", 10);
    } else if (dd == 1) {
      s_strcpy(option, "-cy", 10);
    } else {
      s_strcpy(option, "-cz", 10);
    }
    if (!NumbersOption(argc, argv, 100, option,
                       &opt->ca_count[dd], opt->ca[dd], 'd')) {
      opt->ca_count[dd] = 0;
    } else if ((opt->ca_count[dd] % 2) != 0) { // not even number of numbers
      goto err_constraint;
    }
    for (int i = 0; i < opt->ca_count[dd]; i+=2 ) {
      if (fabs(opt->ca[0][i] - opt->ca[dd][i+1]) < 0.0001) { // same numbers in a pair
        goto err_constraint;
      } else if (!opt->real && (opt->ca[dd][i] < 0 || opt->ca[dd][i] > 1 ||
                                opt->ca[dd][i+1] < 0 || opt->ca[dd][i+1] > 1)) {
        goto err_constraint;
      } else if (opt->ca[dd][i] > opt->ca[dd][i+1]) { // switch so [i] < [i+1]
        SwapDouble(&opt->ca[dd][i], &opt->ca[dd][i+1]);
      }
    }
  }
  // // -cy
  // if (!NumbersOption(argc, argv, 100, "-cy", &opt->c_count[1], opt->cy, 'd')) {
  //   opt->c_count[1] = 0;
  // } else if ((opt->c_count[1] % 2) != 0) { // not even number of numbers
  //   goto err_constraint;
  // }
  // for (int i = 0; i < opt->c_count[1]; i+=2 ) { // same numbers in a pair
  //   if (fabs(opt->cy[i] - opt->cy[i+1]) < 0.0001) {
  //     goto err_constraint;
  //   } else if (opt->cy[i] > opt->cy[i+1]) { // switch so [i] < [i+1]
  //     SwapDouble(&opt->cy[i], &opt->cy[i+1]);
  //   }
  // }
  // // -cz
  // if (!NumbersOption(argc, argv, 100, "-cz", &opt->c_count[2], opt->cz, 'd')) {
  //   opt->c_count[2] = 0;
  // } else if ((opt->c_count[2] % 2) != 0) { // not even number of numbers
  //   goto err_constraint;
  // }
  // for (int i = 0; i < opt->c_count[2]; i+=2 ) { // same numbers in a pair
  //   if (fabs(opt->cz[i] - opt->cz[i+1]) < 0.0001) {
  //     goto err_constraint;
  //   } else if (opt->cz[i] > opt->cz[i+1]) { // switch so [i] < [i+1]
  //     SwapDouble(&opt->cz[i], &opt->cz[i+1]);
  //   }
  // }
  // error - triggers only via goto command
  if (false) {
    err_constraint:
      err_msg("requires coordinate pairs (with two distinct numbers per pair) "
              "in units of box fraction, i.e., <0,1> (unless --real is used)");
      PrintErrorOption("-cx/-cy/-cz");
      exit(1);
  } //}}}
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  if (!opt->real) {
    for (int dd = 0; dd < 3; dd++) {
      opt->move[dd] *= System.Box.Length[dd];
    }
  }

  if (opt->join && Count->Molecule == 0) {
    err_msg("no molecules to join");
    PrintWarning();
  }

  // specify beads to save (possibly using -bt and/or -mt options) //{{{
  /*
   * reverse=true ... save only the specified species
   * reverse=false ... exclude the specified species
   *
   * First setting all bead/molecule types to !reverse and then adjusting this
   * if -bt/-mt options are present correctly specifies which
   * bead/molecule types to save
   */
  // auxiliary arrays holding which bead/molecule types to save
  bool *write_bt = malloc(Count->BeadType * sizeof *write_bt),
       *write_mt = malloc(Count->MoleculeType * sizeof *write_mt);
  // first assume all bead types are saved/excluded based on --reverse option...
  InitBoolArray(write_bt, Count->BeadType, !opt->reverse);
  // ... then adjust if -bt option is present
  opt->bt = TypeOption(argc, argv, "-bt", 'b', opt->reverse, write_bt, System);
  // first assume all molecule types are saved/excluded...
  InitBoolArray(write_mt, Count->MoleculeType, !opt->reverse);
  // ... then adjust if -mt option is present
  opt->mt = TypeOption(argc, argv, "-mt", 'm', opt->reverse, write_mt, System);
  // array for holding which beads to save/exclude
  bool *write = malloc(Count->Bead * sizeof *write);
  // first assume all are saved/excluded...
  InitBoolArray(write, Count->Bead, !opt->reverse);
  // then check possible -mt option...
  if (opt->mt) {
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < mt->Number; j++) {
        int mol = mt->Index[j];
        for (int k = 0; k < mt->nBeads; k++) {
          int id = System.Molecule[mol].Bead[k];
          if (write_mt[i] == opt->reverse) {
            write[id] = opt->reverse; // save/exclude based on --reverse
          }
        }
      }
    }
  }
  // ... and possible -bt option
  if (opt->bt) {
    for (int i = 0; i < Count->Bead; i++) {
      int type = System.Bead[i].Type;
      if (write_bt[type] == opt->reverse) {
        write[i] = opt->reverse;
      }
    }
  }
  // free the auxiliary arrays
  free(write_mt);
  free(write_bt); //}}}

  // '-n' option - specify timestep ids //{{{
  opt->n_number = -1;
  InitIntArray(opt->n_save, 100, 0);
  NumbersOption(argc, argv, 100, "-n", &opt->n_number, opt->n_save, 'i');
  // ignore -st/-e/-sk when -n is used
  if (opt->n_number != -1) {
    opt->c.start = 1;
    opt->c.end = -1;
    opt->c.skip = 1;
  }
  // SortArrayInt(opt->n_save, opt->n_number, 0);
  SortArray(opt->n_save, opt->n_number, 0, 'i'); //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // print initial stuff to output coordinate file //{{{
  if (fout.type == VCF_FILE) {
    PrintByline(fout.name, argc, argv);
  } else if (fout.type == VTF_FILE) {
    WriteStructure(fout, System, -1, false, argc, argv);
  } else { // ensure it's a new file
    FILE *out = OpenFile(fout.name, "w");
    fclose(out);
  } //}}}

  FILE *fr = OpenFile(in.coor.name, "r");
  // main loop //{{{
  // file pointers for finding the last valid step
  fpos_t *position = calloc(1, sizeof *position);
  // save line count at every fgetpos()
  int *bkp_line_count = calloc(1, sizeof *bkp_line_count);
  int n_opt_count = 0, // count saved steps if -n option is used
      count_coor = 0,  // count steps in the vcf file
      count_saved = 0, // count steps in output file
      line_count = 0;  // count lines in the coor file
  while (true) {
    if (opt->last) {
      if (!opt->c.silent && isatty(STDOUT_FILENO)) {
        count_coor++;
        fprintf(stdout, "\rDiscarding step: %d", count_coor);
      }
    } else {
      PrintStep(&count_coor, opt->c.start, opt->c.silent);
    }
    position = s_realloc(position, count_coor * sizeof *position);
    fgetpos(fr, &position[count_coor-1]);
    bkp_line_count = s_realloc(bkp_line_count, count_coor *
                               sizeof *bkp_line_count);
    bkp_line_count[count_coor-1] = line_count;
    // decide whether this timestep is to be saved //{{{
    bool use = false;
    // no -n option - possibly use only if --last not present
    if (opt->n_number == -1 && !opt->last && UseStep(opt->c, count_coor)) {
      use = true;
    // -n option is used - save the timestep if it's in the list
    } else if (n_opt_count < opt->n_number &&
               opt->n_save[n_opt_count] == count_coor) {
      use = true;
      n_opt_count++;
    } //}}}
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (fout.type == LDATA_FILE && count_saved == 1) {
        err_msg("only one timestep can be saved to lammps data file");
        PrintWarnFile(fout.name, "\0", "\0");
        count_coor--;
        break;
      }
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_saved++;
      TransformAndSave(&System, *opt, fout, write, count_coor, argc, argv);
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // decide whether to exit the main loop //{{{
    /* break the loop if
     *    1) all timesteps in the -n option are saved (and --last isn't used)
     *    or
     *    2) end timestep was reached (-e option)
     */
    if (!opt->last && // never break when --last is used
        (n_opt_count == opt->n_number || // 1)
         count_coor == opt->c.end)) {    // 2)
      break;
    } //}}}
  } //}}}
  // if --last option is used, read & save the last timestep //{{{
  if (opt->last) {
    /* To through all saved file positions (last to first) and save a the first
     * valid step encountered.
     * Start at count_coor as the saved position is at the beginning of the last
     * timestep to be skipped, and count_coor-- is used beore quitting the while
     * loop.
     */
    count_coor--; // decrement as the last step should be eof
    for (int i = (count_coor); i >= 0; i--) {
      fsetpos(fr, &position[i]);
      line_count = bkp_line_count[i];
      if (ReadTimestep(in, fr, &System, &line_count)) {
        TransformAndSave(&System, *opt, fout, write, count_coor, argc, argv);
        if (!opt->c.silent) {
          if (isatty(STDOUT_FILENO)) {
            fprintf(stdout, "\r                          \r");
          }
          fprintf(stdout, "Saved Step: %d\n", i+1);
          fflush(stdout);
        }
        break;
      } else {
        snprintf(ERROR_MSG, LINE, "disregarding step %s%d%s", ErrYellow(),
                 i + 1, ErrCyan());
        PrintWarnFile(in.coor.name, "\0", "\0");
      }
    } //}}}
  } else if (count_coor == 0) { // error - input file without a valid timestep //{{{
    remove(fout.name);
    err_msg("no valid timestep found");
    PrintErrorFile(in.coor.name, "\0", "\0"); //}}}
  } else if (opt->c.start > count_coor) { // warn if no timesteps were written //{{{
    remove(fout.name);
    err_msg("no coordinates written (starting timestep is higher "
            "than the total number of timesteps)");
    PrintWarning(); //}}}
  } else if (!opt->c.silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (saved %d)\n", count_coor, count_saved);
    fflush(stdout);
  } //}}}
  fclose(fr);

  // // reading file from the end //{{{
  // fr = OpenFile(in.coor.name, "r");
  // fseek(fr, 0, SEEK_END); // move to the end of the file
  // long pos = ftell(fr); // save the pointer position
  // pos--; // skip EOF character
  // pos--; // assume the last character in the file is '\n', so skip that too
  // char test;
  // line_count = 0;
  // while (pos >= 0) {
  //   fseek(fr, pos, SEEK_SET); // decrement pointer position and go there
  //   test = fgetc(fr);
  //   if (test == '\n' || pos == 0) {
  //     line_count++;
  //     if (pos == 0) {
  //       ungetc(test, fr);
  //     }
  //     ReadAndSplitLine(fr, SPL_STR, " \t\n");
  //     printf("%5d|%d: %s\n", line_count, words, split[0]);
  //   }
  //   pos--;
  // }
  // fclose(fr); //}}}

  // free memory
  FreeSystem(&System);
  free(write);
  free(position);
  free(bkp_line_count);
  free(opt);

  return 0;
}
