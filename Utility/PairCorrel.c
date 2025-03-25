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
PairCorrel utility calculates pair correlation function for specified \
bead types. All pairs of bead types (including same type pairs) are \
calculated - given A and B types, pcf between A-A, A-B and B-B are \
calculated.\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <width> <output> <bead(s)> ", cmd);
  fprintf(ptr, "[options]\n\n");

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a distribution bin\n");
  fprintf(ptr, "<output>            output file with pair correlation "
          "function(s)\n");
  fprintf(ptr, "<bead(s)>           bead name(s) for calculation "
          "(optional and ignored if '--all' is used)\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --all             use all bead types "
          "(overwrites <bead(s)>)\n");
  // fprintf(ptr, "  -m <max>          maximum distance for calculation\n");
  fprintf(ptr, "  -D2 <axis>        assume 2D system (e.g., slit) with "
          "non-periodic condition in <axis> direction");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool all; // --all
  int axis[3];
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 2, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, false, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--all", "-D2");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  } //}}}

  // <output> - filename with pcf(s)
  char fout_pcf[LINE] = "";
  s_strcpy(fout_pcf, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  opt->all = BoolOption(argc, argv, "--all");
  char str[LINE];
  if (FileOption(argc, argv, "-D2", str)) {
    if (str[0] == 'x') {
      opt->axis[0] = 1;
      opt->axis[1] = 2;
      opt->axis[2] = 0;
    } else if (str[0] == 'y') {
      opt->axis[0] = 0;
      opt->axis[1] = 2;
      opt->axis[2] = 1;
    } else if (str[0] == 'z') {
      opt->axis[0] = 0;
      opt->axis[1] = 1;
      opt->axis[2] = 2;
    } else {
      err_msg("requires argument 'x', 'y', or 'z'");
      PrintErrorOption("-D2");
      exit(1);
    }
  } else {
    opt->axis[0] = -1;
  }
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  double *box = System.Box.Length;

  // <bead(s)> - names of bead types to use //{{{
  if (opt->all) {
    for (int i = 0; i < Count->BeadType; i++) {
      System.BeadType[i].Flag = true;
    }
  } else {
    for (int i = 0; i < Count->BeadType; i++) {
      System.BeadType[i].Flag = false;
    }
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], System);
      if (type == -1) {
        ErrorBeadType(argv[count], System);
        exit(1);
      }
      if (System.BeadType[type].Flag) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
                 ErrYellow(), argv[count], ErrCyan());
        PrintWarning();
      }
      System.BeadType[type].Flag = true;
    }
    count--; // while always increments count at least once
    if (count < (req_arg + 1)) {
      err_msg("missing <bead(s)> or --all option");
      PrintError();
      PrintCommand(stderr, argc, argv);
      Help(StripPath(argv[0]), true, common, option);
      exit(1);
    }
  } //}}}

  // write initial stuff to output pcf file //{{{
  PrintByline(fout_pcf, argc, argv);
  FILE *out = OpenFile(fout_pcf, "a");
  fprintf(out, "# (1) distance");
  // print bead type names to output file //{{{
  count = 1;
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = i; j < Count->BeadType; j++) {
      if (System.BeadType[i].Flag && System.BeadType[j].Flag) {
        count++;
        fprintf(out, " (%d) %s-%s", count, System.BeadType[i].Name,
                                           System.BeadType[j].Name);
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  int bins;
  double max_dist;
  // TODO: not using 2D-dependent min/max because I want the brute force to
  //       finish at some time...
  // if (opt->axis[0] == -1) {
  //   bins = Max3(box[0], box[1], box[2]) / width;
  //   max_dist = 0.5 * Min3(box[0], box[1], box[2]);
  // } else {
  //   bins = Max3(box[opt->axis[0]], box[opt->axis[0]], box[opt->axis[1]]);
  //   bins /= width;
  //   max_dist = Min3(box[opt->axis[0]], box[opt->axis[0]], box[opt->axis[1]]);
  //   max_dist *= 0.5;
  // }
  bins = Max3(box[0], box[1], box[2]) / width;
  max_dist = 0.5 * Min3(box[0], box[1], box[2]);

  // allocate memory //{{{
  // array counting number of pairs
  long int **counter = calloc(Count->BeadType, sizeof *counter);
  long int *counter2 = calloc(Count->BeadType, sizeof *counter2);
  // array for counting individual used beads of each type
  // pair correlation function
  int ***pcf = malloc(Count->BeadType * sizeof **pcf);
  for (int i = 0; i < Count->BeadType; i++) {
    counter[i] = calloc(Count->BeadType, sizeof *counter[i]);
    pcf[i] = malloc(Count->BeadType * sizeof *pcf[i]);
    for (int j = 0; j < Count->BeadType; j++) {
      pcf[i][j] = calloc(bins, sizeof *pcf[i][j]);
    }
  } //}}}

  // // TODO: for cell-linked list
  // double cell_size = 5;
  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count timesteps from the beginning
      count_used = 0, // count steps used for calculation
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      // // TODO: trying cell-linked list (unsuccessfully) //{{{
      // // double cell_size = 3;
      // int n_cells[3], *Head, *Link, Dc[14][3];
      // LinkedList(System, &Head, &Link, cell_size, n_cells, Dc);
      // int c1[3];
      // for (c1[2] = 0; c1[2] < n_cells[2]; c1[2]++) {
      //   for (c1[1] = 0; c1[1] < n_cells[1]; c1[1]++) {
      //     for (c1[0] = 0; c1[0] < n_cells[0]; c1[0]++) {
      //       int cell1 = SelectCell1(c1, n_cells);
      //       // select first bead in the cell 'cell1'
      //       int i = Head[cell1];
      //       while (i != -1) {
      //         BEAD *b_i = &System.Bead[System.BeadCoor[i]];
      //         if (!System.BeadType[b_i->Type].Flag) {
      //           i = Link[i];
      //           continue;
      //         }
      //         for (int k = 0; k < 14; k++) {
      //           int cell2 = SelectCell2(c1, n_cells, Dc, k);
      //
      //           int j;
      //           if (cell1 == cell2) { // next bead in 'cell1'
      //             j = Link[i];
      //           } else { // first bead in 'cell2'
      //             j = Head[cell2];
      //           }
      //           while (j != -1) {
      //             BEAD *b_j = &System.Bead[System.BeadCoor[j]];
      //             if (!System.BeadType[b_j->Type].Flag) {
      //               j = Link[j];
      //               continue;
      //             }
      //             int btype_i = b_i->Type;
      //             int btype_j = b_j->Type;
      //             if (btype_i > btype_j) {
      //               SwapInt(&btype_i, &btype_j);
      //             }
      //             counter[btype_i][btype_j]++;
      //             // calculate distance between i and j beads
      //             double dist[3];
      //             Distance(b_i->Position, b_j->Position, box, dist);
      //             dist[0] = VectLength(dist);
      //             if (dist[0] < max_dist) {
      //               int l = dist[0] / width;
      //               pcf[btype_i][btype_j][l]++;
      //             }
      //             j = Link[j];
      //           }
      //         }
      //         i = Link[i];
      //       }
      //     }
      //   }
      // }
      // free(Head);
      // free(Link); //}}}
      for (int i = 0; i < Count->BeadCoor; i++) { //{{{
        int id_i = System.BeadCoor[i];
        BEAD *b_i = &System.Bead[id_i];
        if (!System.BeadType[b_i->Type].Flag) {
          continue;
        }
        for (int j = (i + 1); j < Count->BeadCoor; j++) {
          int id_j = System.BeadCoor[j];
          BEAD *b_j = &System.Bead[id_j];
          if (!System.BeadType[b_j->Type].Flag) {
            continue;
          }
          int btype_i = b_i->Type;
          int btype_j = b_j->Type;
          if (btype_i > btype_j) {
            SwapInt(&btype_i, &btype_j);
          }
          counter[btype_i][btype_j]++;
          double temp[2];
          // make the non-periodic coordinate 0
          if (opt->axis[0] != -1) {
            temp[0] = b_i->Position[opt->axis[2]];
            b_i->Position[opt->axis[2]] = 0;
            temp[1] = b_j->Position[opt->axis[2]];
            b_j->Position[opt->axis[2]] = 0;
          }
          double dist[3];
          Distance(b_i->Position, b_j->Position, box, dist);
          // return the non-periodic coordinate - just pro forma
          if (opt->axis[0] != -1) {
            b_i->Position[opt->axis[2]] = temp[0];
            b_j->Position[opt->axis[2]] = temp[1];
          }
          dist[0] = VectLength(dist);
          if (dist[0] < max_dist) {
            int l = dist[0] / width;
            pcf[btype_i][btype_j][l]++;
          }
        }
      } //}}}
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // write data to output file(s) //{{{
  out = OpenFile(fout_pcf, "a");

  // calculate pcf
  for (int j = 0; j < bins; j++) {
    if ((width * (j+1)) > max_dist) {
      break;
    }
    // calculate volume of every shell that will be averaged
    double shell;
    // radius of outer and inner sphere
    double rad[2] = {width * (j + 1), width * j};
    if (opt->axis[0] == -1) {
      // maximum radius of complete sphere
      double max_r = Min3(box[0], box[1], box[2]) / 2;
      // volume of outer and inner spheres
      double sphere[2] = {4.0 / 3 * Cube(rad[0]), 4.0 / 3 * Cube(rad[1])};
      // volume of outer and inner sphere's cut-off tops (0 for full sphere)
      double top[2] = {0, 0};
      if (rad[0] > max_r) { // is the outer sphere cut-off?
        // volume of one cut-off spherical top of the outer sphere
        top[0] = Square(rad[0] - max_r) * (2 * rad[0] + max_r) / 3;
        if (rad[1] > max_r) { // is the inner sphere cut-off?
          // volume of one cut-off spherical top of the inner sphere
          top[1] = Square(rad[1] - max_r) * (2 * rad[1] + max_r) / 3;
        }
      }
      // volume is outer sphere w/o its tops minus inner sphere w/o its tops
      shell = PI * (sphere[0] - 2 * top[0] - (sphere[1] - 2 * top[1]));
    } else {
      // maximum radius of complete circle
      double max_r = Min3(box[opt->axis[0]], box[opt->axis[1]], HIGHNUM) / 2;
      // area of outer and inner circles
      double circle[2] = {PI * Square(rad[0]), PI * Square(rad[1])};
      // area of outer and inner circle's cut-off tops (0 for full circle)
      double top[2] = {0, 0};
      if (rad[0] > max_r) { // is the outer circle cut-off?
        double phi[2] = {2 * acos(max_r / rad[0]), 2 * acos(max_r / rad[1])};
        // area of one cut-off top of the outer circle
        top[0] = Square(rad[0]) / 2 * (phi[0] - sin(phi[0]));
        if (rad[1] > max_r) { // is the inner circle cut-off?
          // area of one cut-off top of the inner circle
          top[1] = Square(rad[1]) / 2 * (phi[1] - sin(phi[1]));
        }
      }
      // area is outer circle w/o its tops minus inner circle w/o its tops
      shell = circle[0] - 2 * top[0] - (circle[1] - 2 * top[1]);
    }
    // write middle distance of the bin
    fprintf(out, "%8.5f", (rad[0] + rad[1]) / 2);
    // write pcf for all pairs for the given bin
    for (int k = 0; k < Count->BeadType; k++) {
      for (int l = k; l < Count->BeadType; l++) {
        if (System.BeadType[k].Flag && System.BeadType[l].Flag) {
          double norm_factor = System.Box.Volume / (counter[k][l] * shell);
          if (opt->axis[0] != -1) {
            norm_factor /= System.Box.Length[opt->axis[2]];
          }
          // // TODO: norm factor for linked-list via effective volume
          // double norm_factor = Cube(3 * cell_size) / (counter[k][l] * shell);
          fprintf(out, " %10f", pcf[k][l][j] * norm_factor);
        }
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      free(pcf[i][j]);
    }
    free(pcf[i]);
    free(counter[i]);
  }
  free(pcf);
  free(counter);
  free(counter2);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
