#include "../AnalysisTools.h"

/* Debug: print Head chain for the two problem beads if present */
void DumpCellChain(const char *tag, int cell, int *Head, int *Link, const SYSTEM *S) {
  printf("Dump %s cell %d: Head=%d\n", tag, cell, Head[cell]);
  for (int p = Head[cell]; p != -1; p = Link[p]) {
    int id = S->BeadCoor[p];
    printf("  list-index p=%d -> bead-id=%d pos=(%g,%g,%g)\n",
           p, id,
           S->Bead[id].Position.v[0],
           S->Bead[id].Position.v[1],
           S->Bead[id].Position.v[2]);
  }
}


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
          "non-periodic condition in <axis> direction\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool all; // --all
  int axis[3];
  bool *bt;
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// TODO: <bead(s)> mandatory to -bt opt (default - all)
// TODO: move to some library file
static inline void CorrectBTypeOrder(int *btype_i, int *btype_j) {
  if (btype_i > btype_j) {
    SwapInt(btype_i, btype_j);
  }
}

// calculate distance between i and j beads //{{{
void CalculatePCF(int id_i, int id_j, SYSTEM System,
                  ArrNDi *pcf, int bins, double max_dist, double width) {
  int i = System.BeadCoor[id_i];
  int j = System.BeadCoor[id_j];
  BEAD *b_i = &System.Bead[i];
  BEAD *b_j = &System.Bead[j];
  int btype_i = b_i->Type;
  int btype_j = b_j->Type;
  CorrectBTypeOrder(&btype_i, &btype_j);
  // calculate distance between the two beads
  vec3d d = Distance(b_i->Position.v, b_j->Position.v, System.Box.Length);
  double dist = VectLength(d);
  if (dist < max_dist) {
    int l = dist / width;
    if (l < bins) {
      AddArr3D(pcf, btype_i, btype_j, l, 1);
    }
  }
} //}}}

// callback function to call for each pair of particles
// i & j ... indices in System.BeadCoor array
typedef void (*pair_cb_t)(int i, int j, const SYSTEM System, void *userdata);
// callback function for checking if bead/mol should be used
// btype ... bead type
typedef bool (*check_cb_t)(int btype, void *userdata);
// bead-pair traversing function (linked list) //{{{
void TraverseLinkedListPairs(const SYSTEM System, const double cell_size,
                             pair_cb_t pair_callback, void *ud_pair,
                             check_cb_t check_callback, void *check_ud) {
  int *Head, *Link;
  vec3i n_cells = LinkedList(System, &Head, &Link, cell_size);

  // neighbour offsets: home cell + 13 half-shell
  const vec3i neighbour[14] = {
    { .v = { 0, 0, 0} },
    { .v = { 1, 0, 0} },
    { .v = { 1, 1, 0} },
    { .v = {-1, 1, 0} },
    { .v = { 0, 1, 0} },
    { .v = { 0, 0, 1} },
    { .v = {-1, 0, 1} },
    { .v = { 1, 0, 1} },
    { .v = {-1,-1, 1} },
    { .v = { 0,-1, 1} },
    { .v = { 1,-1, 1} },
    { .v = {-1, 1, 1} },
    { .v = { 0, 1, 1} },
    { .v = { 1, 1, 1} },
  };
  vec3i c1;
  for (c1.z = 0; c1.z < n_cells.z; c1.z++) {
    for (c1.y = 0; c1.y < n_cells.y; c1.y++) {
      for (c1.x = 0; c1.x < n_cells.x; c1.x++) {
        int cell1 = SelectCell1(c1, n_cells);
        int i = Head[cell1];
        while (i != -1) {
          int id_i = System.BeadCoor[i];
          if (!check_callback(System.Bead[id_i].Type, check_ud)) {
            i = Link[i];
            continue;
          }
          // loop over all 14 neighbour offsets (home + 13 neighbours)
          for (int k = 0; k < 14; k++) {
            int cell2 = SelectCell2(c1, n_cells, neighbour, k);

            int j = Head[cell2];
            while (j != -1) {
              int id_j = System.BeadCoor[j];
              if (!check_callback(System.Bead[id_j].Type, check_ud)) {
                j = Link[j];
                continue;
              }
              // avoid double-counting in home cell
              if (cell1 != cell2 || i < j) {
                pair_callback(i, j, System, ud_pair);
              }
              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  }
  free(Head);
  free(Link);
} //}}}
// bead-pair traversing function (brute O(N^2) nested loops) //{{{
void TraverseBrutePairs(const SYSTEM System,
                        pair_cb_t pair_callback, void *ud_pair,
                        check_cb_t check_callback, void *check_ud) {
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id_i = System.BeadCoor[i];
    if (!check_callback(System.Bead[id_i].Type, check_ud)) {
      continue;
    }
    for (int j = (i + 1); j < System.Count.BeadCoor; j++) {
      int id_j = System.BeadCoor[j];
      if (!check_callback(System.Bead[id_j].Type, check_ud)) {
        continue;
      }
      pair_callback(i, j, System, ud_pair);
    }
  }
} //}}}
bool CheckBeadType(int btype, OPT *opt) { //{{{
  if (!opt->bt[btype]) {
    return false;
  } else {
    return true;
  }
} //}}}

struct pcf_args {
  ArrNDi *pcf;
  int bins;
  double max_dist;
  double width;
};
void CalculatePCF_adaptor(int id_i, int id_j, const SYSTEM System, void *ud) {
  struct pcf_args *p = (struct pcf_args*)ud;
  CalculatePCF(id_i, id_j, System, p->pcf, p->bins, p->max_dist, p->width);
}

struct check_args {
  OPT *opt;
};
bool CheckBeadType_adaptor(int type, void *ud) {
  struct check_args *p = (struct check_args*)ud;
  return CheckBeadType(type, p->opt);
}

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
  opt->bt = calloc(Count->BeadType, sizeof *opt->bt);
  if (opt->all) {
    for (int i = 0; i < Count->BeadType; i++) {
      opt->bt[i] = true;
    }
  } else {
    for (int i = 0; i < Count->BeadType; i++) {
      opt->bt[i] = false;
    }
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], System);
      if (type == -1) {
        ErrorBeadType(argv[count], System);
        exit(1);
      }
      if (opt->bt[type]) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
                 ErrYellow(), argv[count], ErrCyan());
        PrintWarning();
      }
      opt->bt[type] = true;
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
  FILE *out = PrintBylineOpenFile(fout_pcf, argc, argv);
  fprintf(out, "# (1) distance");
  // print bead type names to output file
  count = 1;
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = i; j < Count->BeadType; j++) {
      if (CheckBeadType(i, opt) && CheckBeadType(j, opt)) {
        count++;
        fprintf(out, " (%d) %s-%s", count,
                System.BeadType[i].Name, System.BeadType[j].Name);
      }
    }
  }
  putc('\n', out);
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
  // TODO: shitty stuff - should be -D2 opt dependant or some such
  // max_dist = 0.5 * Min3(box[0], box[1], box[2]);
  max_dist = 0.5 * Max3(box[0], box[1], box[2]);
  max_dist = 10;
  double cell_size = max_dist;

  // pair correlation function
  ArrNDi *pcf = CreateArr3Di(Count->BeadType, Count->BeadType, bins);

  // main loop //{{{
  struct check_args check = { opt }; // for CheckBeadType_adaptor
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count timesteps from the beginning
      count_used = 0, // count steps used for calculation
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    if (UseStep(opt->c, count_coor)) {
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      struct pcf_args args = { pcf, bins, max_dist, width };
      if (true) {
        TraverseLinkedListPairs(System, cell_size, CalculatePCF_adaptor, &args,
                                CheckBeadType_adaptor, &check);
      } else {
        TraverseBrutePairs(System, CalculatePCF_adaptor, &args,
                           CheckBeadType_adaptor, &check);
      }
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
    // calculate volume of every shell that will be averaged, taking into
    // account corrections for truncated spheres
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
      shell = PI * (sphere[0] - 6 * top[0] - (sphere[1] - 6 * top[1]));
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
    // write data
    fprintf(out, "%8.5f", (rad[0] + rad[1]) / 2);
    // write pcf for all pairs for the given bin
    for (int k = 0; k < Count->BeadType; k++) {
      for (int l = k; l < Count->BeadType; l++) {
        BEADTYPE *bt_k = &System.BeadType[k];
        BEADTYPE *bt_l = &System.BeadType[k];
        if (CheckBeadType(k, opt) && CheckBeadType(l, opt)) {
          int pairs = bt_k->Number * (bt_l->Number - 1) / 2;
          double norm_factor = System.Box.Volume / (shell * pairs * count_used);
          if (opt->axis[0] != -1) {
            norm_factor /= System.Box.Length[opt->axis[2]];
          }
          fprintf(out, " %10f", GetArr3D(pcf, k, l, j) * norm_factor);
        }
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeArrND(pcf);
  free(opt->bt);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
