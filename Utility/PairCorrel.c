#include "../src/AnalysisTools.h"
// TODO: <bead(s)> mandatory to -bt opt (default - all)
// TODO: 2D version

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "PairCorrel utility calculates pair correlation function for specified "
  "bead types. All pairs of bead types (including same type pairs) are "
  "calculated - given A and B types, pcf between A-A, A-B and B-B are "
  "calculated.",

  "Usage: PairCorrel <input> <width> <output> <bead(s)> [options]",
  .args = 3, // number of mandatory arguments
  .all = 15, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<width>", NULL, "width of a distribution bin", OPT_ARG},
  {"<output>", NULL, "output file with pair correlation function(s)", OPT_ARG},
  {"<bead(s)>", NULL, "bead name(s) for calculation (optional and ignored if '--all' is used)", OPT_ARG},
  {"-d", "<dist>", "maximum distance for RDF calculation (default: 1/3 of the shortest box side length)", OPT_EXTRA},
  {"--all", NULL, "use all bead types (overwrites <bead(s)>)", OPT_EXTRA},
  // {"-D2", "<axis>", "assume 2D system (e.g., slit) with non-periodic condition in <axis> direction", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  bool all; // --all
  vec3i axis;
  bool *bt;
  double max_dist;
}; //}}}

// TODO: move to some library file
static inline void CorrectBTypeOrder(int *btype_i, int *btype_j) {
  if (*btype_i > *btype_j) {
    SwapInt(btype_i, btype_j);
  }
}

// functions to plug into traversal functions
// calculate PCF, i.e., distance between i and j beads //{{{
// the calculation itself
static void CalculatePCF(int id_i, int id_j, SYSTEM System,
                         ArrNDi *pcf, int bins, double max_dist, double width) {
  int i = System.BeadCoor[id_i];
  int j = System.BeadCoor[id_j];
  BEAD *b_i = &System.Bead[i];
  BEAD *b_j = &System.Bead[j];
  // calculate distance between the two beads
  vec3d d = Distance(b_i->Position.v, b_j->Position.v, System.Box.Length);
  double dist = VectLength(d);
  if (dist < max_dist) {
    int l = dist / width;
    if (l < bins) {
      int btype_i = b_i->Type;
      int btype_j = b_j->Type;
      CorrectBTypeOrder(&btype_i, &btype_j);
      AddArr3D(pcf, btype_i, btype_j, l, 1);
    }
  }
}
// structure for the callback function
struct pcf_args {
  ArrNDi *pcf;
  int bins;
  double max_dist;
  double width;
};
// adaptor for the CalculatePCF() function
static void CalculatePCF_adaptor(int id_i, int id_j,
                                 const SYSTEM System, void *ud) {
  struct pcf_args *p = (struct pcf_args*)ud;
  CalculatePCF(id_i, id_j, System, p->pcf, p->bins, p->max_dist, p->width);
} //}}}
// condition for using specified beads //{{{
// check based on supplied type (needed for writing to file)
static bool CheckBeadType(int btype, OPT opt) {
  if (!opt.bt[btype]) {
    return false;
  } else {
    return true;
  }
}
// check based on supplied System.BeedCoor id
static bool CheckBead(int id_i, SYSTEM System, OPT opt) {
  int i = System.BeadCoor[id_i];
  return CheckBeadType(System.Bead[i].Type, opt);
}
// structure for the callback function
struct check_args {
  OPT opt;
};
// adaptor for the CalculatePCF() function
static bool CheckBeadType_adaptor(int id_i, SYSTEM System, void *ud) {
  struct check_args *p = (struct check_args*)ud;
  return CheckBead(id_i, System, p->opt);
}
//}}}

void Calculation(SYSTEM *System, OPT opt, ArrNDi *pcf, int bins,
                 double width, double cell_size) {
  WrapJoinCoordinates(System, true, false);
  struct pcf_args args = { pcf, bins, opt.max_dist, width };
  struct check_args check = { opt };
  TraversePairs(*System, cell_size, CalculatePCF_adaptor, &args,
                CheckBeadType_adaptor, &check);
}
// structure for the callback function
struct user_data {
  OPT opt;
  ArrNDi *pcf;
  int bins;
  double width, cell_size;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, p->opt, p->pcf, p->bins, p->width, p->cell_size);
};

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, false, HelpDesc, opts);
  OPT opt;
  int count = 0;
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
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <output> - filename with pcf(s)
  char fout_pcf[LINE] = "";
  s_strcpy(fout_pcf, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }
  // use all bead types in the structure file?
  opt.all = BoolOption(argc, argv, "--all");
  // 2D calculation?
  char str[LINE];
  if (FileOption(argc, argv, "-D2", str)) {
    if (str[0] == 'x') {
      opt.axis.v[0] = 1;
      opt.axis.v[1] = 2;
      opt.axis.v[2] = 0;
    } else if (str[0] == 'y') {
      opt.axis.v[0] = 0;
      opt.axis.v[1] = 2;
      opt.axis.v[2] = 1;
    } else if (str[0] == 'z') {
      opt.axis.v[0] = 0;
      opt.axis.v[1] = 1;
      opt.axis.v[2] = 2;
    } else {
      err_msg("requires argument 'x', 'y', or 'z'");
      PrintErrorOption("-D2");
      exit(1);
    }
  } else {
    opt.axis.v[0] = -1;
  }
  //}}}

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  vec3d box = System.Box.Length;

  // <bead(s)> - names of bead types to use //{{{
  if (!(opt.bt = calloc(Count->BeadType, sizeof *opt.bt))) {
    ErrorAlloc("opt.bt");
  }
  if (opt.all) {
    for (int i = 0; i < Count->BeadType; i++) {
      opt.bt[i] = true;
    }
  } else {
    for (int i = 0; i < Count->BeadType; i++) {
      opt.bt[i] = false;
    }
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], System);
      if (type == -1) {
        ErrorBeadType(argv[count], System);
        exit(1);
      }
      if (opt.bt[type]) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
                 ErrYellow(), argv[count], ErrCyan());
        PrintWarning();
      }
      opt.bt[type] = true;
    }
    count--; // while always increments count at least once
    if (count < (HelpDesc.args + 1)) {
      err_msg("missing <bead(s)> or --all option");
      PrintError();
      PrintCommand(stderr, argc, argv);
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  int bins;
  opt.max_dist = Min3(box.x, box.y, box.z) / 3;
  if (OneNumberOption(argc, argv, "-d", &opt.max_dist, 'd') &&
      opt.max_dist <= 0) {
    err_msg("distance must be a positive number");
    ErrorOption("-d");
  }
  double cell_size = opt.max_dist;
  bins = opt.max_dist / width;

  // pair correlation function
  ArrNDi *pcf = CreateArr3Di(Count->BeadType, Count->BeadType, bins);
  if (!pcf) {
    ErrorAlloc("pcf");
  }

  STEP step = InitStep;
  struct user_data ud = { opt, pcf, bins, width, cell_size };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // write data to output file(s) //{{{
  // header
  FILE *out = PrintBylineOpenFile(fout_pcf, argc, argv);
  fprintf(out, "# (1) distance");
  // print bead type names to output file
  int ncols = 1;
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = i; j < Count->BeadType; j++) {
      if (CheckBeadType(i, opt) && CheckBeadType(j, opt)) {
        ncols++;
        fprintf(out, " (%d) %s-%s", ncols,
                System.BeadType[i].Name, System.BeadType[j].Name);
      }
    }
  }
  putc('\n', out);
  // collate data
  int nrows = bins;
  ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
  if (!data) {
    ErrorAlloc("data");
  }
  for (int i = 0; i < nrows; i++) {
    // calculate volume/surface of a shell  //{{{
    // account (somewhat) for truncated stuff (high max_dist)
    double shell;
    // radius of outer and inner sphere
    double rad[2] = {width * (i + 1), width * i};
    // 3D (sphere)
    if (opt.axis.v[0] == -1) {
      // maximum radius of complete sphere
      double max_r = Min3(box.x, box.y, box.z) / 2;
      // volume of outer and inner spheres
      double sphere[2] = {4.0 / 3 * Cube(rad[0]), 4.0 / 3 * Cube(rad[1])};
      // volume of outer and inner sphere's cut-off tops (0 for full sphere)
      double top[2] = {0, 0};
      if (rad[0] > max_r) { // is the outer sphere cut off?
        // volume of one cut-off spherical top of the outer sphere
        top[0] = Square(rad[0] - max_r) * (2 * rad[0] + max_r) / 3;
        if (rad[1] > max_r) { // is the inner sphere cut-off?
          // volume of one cut-off spherical top of the inner sphere
          top[1] = Square(rad[1] - max_r) * (2 * rad[1] + max_r) / 3;
        }
      }
      // volume is outer sphere w/o its tops minus inner sphere w/o its tops
      // assumes cubic box - hence the 6
      shell = PI * (sphere[0] - 6 * top[0] - (sphere[1] - 6 * top[1]));
    // 2D (circle)
    } else {
      // maximum radius of complete circle
      double max_r;
      if (box.v[opt.axis.v[0]] < box.v[opt.axis.v[1]]) {
        max_r = box.v[opt.axis.v[0]];
      } else {
        max_r = box.v[opt.axis.v[1]];
      }
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
      // assumes square - hance the 4
      shell = circle[0] - 4 * top[0] - (circle[1] - 4 * top[1]);
    } //}}}
    count = 0;
    SetArr2D(data, i, count++, (rad[0] + rad[1]) / 2);
    // fprintf(out, "%8.5f", (rad[0] + rad[1]) / 2);
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        if (CheckBeadType(j, opt) && CheckBeadType(k, opt)) {
          BEADTYPE *bt_j = &System.BeadType[j];
          BEADTYPE *bt_k = &System.BeadType[k];
          int pairs = bt_j->Number;
          if (j != k) {
            pairs *= bt_k->Number;
          } else {
            pairs *= (bt_k->Number - 1) / 2;
          }
          double norm_factor = System.Box.Volume / (shell * pairs * step.used);
          if (opt.axis.v[0] != -1) { // 2D ...TODO: implement
            norm_factor /= System.Box.Length.v[opt.axis.v[2]];
          }
          SetArr2D(data, i, count++, GetArr3D(pcf, j, k, i) * norm_factor);
        }
      }
    }
  }
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(out, nrows, ncols, data);
  FreeArrND(data);
  fclose(out); //}}}

  // free memory //{{{
  FreeArrND(pcf);
  free(opt.bt);
  FreeSystem(&System);
  //}}}

  return 0;
}

// backup - will the MainLoop() function work?
// // main loop //{{{
// FILE *fr = OpenFile(in.coor.name, "r");
// int count_coor = 0, // count timesteps from the beginning
//     count_used = 0, // count steps used for calculation
//     line_count = 0; // count lines in the vcf file
// while (true) {
//   PrintStep(&count_coor, commons.start, commons.silent);
//   if (UseStep(commons, count_coor)) {
//     if (!ReadTimestep(in, fr, &System, &line_count)) {
//       count_coor--;
//       break;
//     }
//     count_used++;
//     WrapJoinCoordinates(&System, true, false);
//     struct pcf_args args = { pcf, bins, opt.max_dist, width };
//     struct check_args check = { opt };
//     TraversePairs(System, cell_size, CalculatePCF_adaptor, &args,
//                   CheckBeadType_adaptor, &check);
//   } else {
//     if (!SkipTimestep(in, fr, &line_count)) {
//       count_coor--;
//       break;
//     }
//   }
//   // exit the main loop if reached user-specied end timestep
//   if (count_coor == commons.end) {
//     break;
//   }
// }
// fclose(fr);
// PrintLastStep(count_coor, count_used, commons.silent); //}}}
