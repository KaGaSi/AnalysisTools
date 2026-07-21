#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "PairCorrel utility calculates pair correlation function for specified "
  "bead types. By default, alll pairs of bead types (including same type "
  "pairs) are calculated, but --pairs option can be used to explicitly specify "
  "what bead type pairs to use.",

  "Usage: PairCorrel <input> <width> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 16, // number of valid lines OptSpec (not counting last {nullptr})
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
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<width>", nullptr, "width of a distribution bin", OPT_ARG},
  {"<output>", nullptr, "output file with pair correlation function(s)", OPT_ARG},
  {"-bt", "<name(s)>", "bead types to use (default: all); with --pairs, "
    "names are read as consecutive pairs A B [C D ...]", OPT_EXTRA},
  {"--pairs", nullptr, "-bt specifies explicit bead type pairs instead of "
    "individual types (requires -bt)", OPT_EXTRA},
  {"-d", "<dist>", "maximum distance for RDF calculation "
    "(default: 1/3 of the shortest box side length)", OPT_EXTRA},
  {"-D2", "<axis>", "assume 2D system (e.g., slit) with non-periodic "
    "condition in <axis> direction", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  bool pairs;      // --pairs
  vec3i axis;      // -D2
  bool *bt;        // -bt
  ArrNDb *bt_pair; // bt pairs (affected by --pairs)
  double max_dist; // -d
}; //}}}

// functions to plug into traversal functions
// calculate PCF, i.e., distance between i and j beads //{{{
// the calculation itself
static void CalculatePCF(int id_i, int id_j, SYSTEM System,
                         ArrNDi *pcf, ArrNDb *bt_pair,
                         int bins, double max_dist, double width, vec3i axis) {
  int i = System.BeadCoor[id_i];
  int j = System.BeadCoor[id_j];
  BEAD *b_i = &System.Bead[i];
  BEAD *b_j = &System.Bead[j];
  // calculate distance between the two beads
  double dist;
  if (axis.v[0] == -1) {
    vec3d d = DistancePBC(b_i->Position, b_j->Position, &System.Box);
    dist = VectLength(d);
  } else {
    // 2D: minimum-image distance in the two in-plane axes only
    double d2 = 0;
    for (int ax = 0; ax < 2; ax++) {
      double delta = b_i->Position.v[axis.v[ax]] - b_j->Position.v[axis.v[ax]];
      double len = System.Box.Length.v[axis.v[ax]];
      delta -= len * round(delta / len);
      d2 += Square(delta);
    }
    dist = sqrt(d2);
  }
  if (dist < max_dist) {
    int l = dist / width;
    if (l < bins) {
      int btype_i = b_i->Type;
      int btype_j = b_j->Type;
      SortPairAsc(&btype_i, &btype_j);
      if (GetArr2D(bt_pair, btype_i, btype_j)) {
        AddArr3D(pcf, btype_i, btype_j, l, 1);
      }
    }
  }
}
// structure for the callback function
struct pcf_args {
  ArrNDi *pcf;
  ArrNDb *bt_pair;
  int bins;
  double max_dist;
  double width;
  vec3i axis;
};
// adaptor for the CalculatePCF() function
static void CalculatePCF_adaptor(int id_i, int id_j,
                                 const SYSTEM System, void *ud) {
  struct pcf_args *p = (struct pcf_args*)ud;
  CalculatePCF(id_i, id_j, System, p->pcf, p->bt_pair,
               p->bins, p->max_dist, p->width, p->axis);
} //}}}
// condition for using specified beads //{{{
static bool CheckBead(int id_i, SYSTEM System, OPT opt) {
  int i = System.BeadCoor[id_i];
  return opt.bt[System.Bead[i].Type];
}
// adaptor for the callback function
struct check_args {
  OPT opt;
};
static bool CheckBeadType_adaptor(int id_i, SYSTEM System, void *ud) {
  struct check_args *p = (struct check_args*)ud;
  return CheckBead(id_i, System, p->opt);
}
//}}}

void Calculation(SYSTEM *System, OPT opt, ArrNDi *pcf, int bins,
                 double width, double cell_size) {
  WrapJoinCoordinates(System, true, false);
  struct pcf_args args = { pcf, opt.bt_pair, bins,
                           opt.max_dist, width, opt.axis };
  struct check_args check = { opt };
  if (opt.axis.v[0] == -1) {
    TraversePairs(*System, cell_size, CalculatePCF_adaptor, &args,
                  CheckBeadType_adaptor, &check);
  } else {
    // 2D: in-plane distances, so cells must not be binned along the
    // non-periodic axis (beads far apart along it can still be close in-plane)
    TraversePairs2D(*System, cell_size, opt.axis.v[2], CalculatePCF_adaptor,
                    &args, CheckBeadType_adaptor, &check);
  }
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
  opt.pairs = BoolOption(argc, argv, "--pairs");
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

  // -bt/--pairs options //{{{
  opt.bt = calloc(Count->BeadType, sizeof *opt.bt);
  opt.bt_pair = CreateArr2Db(Count->BeadType, Count->BeadType);
  if (!opt.bt || !opt.bt_pair) {
    ErrorAlloc("opt.bt/opt.bt_pair");
  }
  if (opt.pairs) {
    if (!TypeOptionPair(argc, argv, "-bt", 'b', true, opt.bt_pair, System)) {
      err_msg("option -bt is mandatory with --pairs");
      PrintErrorOption("--pairs");
      exit(1);
    }
    // derive per-type participation from the named pairs
    TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System);
  } else {
    if (TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System)) {
      // enable all pairs between selected types
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = 0; j < Count->BeadType; j++) {
          SetArr2D(opt.bt_pair, i, j, opt.bt[i] && opt.bt[j]);
        }
      }
    } else {
      // no -bt: use all bead types and all pairs
      InitBoolArray(opt.bt, Count->BeadType, true);
      FillArrND(opt.bt_pair, true);
    }
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  int bins;
  if (opt.axis.v[0] == -1) {
    opt.max_dist = Min3(box.x, box.y, box.z) / 3;
  } else {
    opt.max_dist = fmin(box.v[opt.axis.v[0]], box.v[opt.axis.v[1]]) / 3;
  }
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
      if (GetArr2D(opt.bt_pair, i, j)) {
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
      // maximum radius of fully inscribed circle
      double max_r = fmin(box.v[opt.axis.v[0]], box.v[opt.axis.v[1]]) / 2;
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
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        if (GetArr2D(opt.bt_pair, j, k)) {
          BEADTYPE *bt_j = &System.BeadType[j];
          BEADTYPE *bt_k = &System.BeadType[k];
          int pairs;
          if (j != k) {
            pairs = bt_j->Number * bt_k->Number;
          } else {
            pairs = bt_j->Number * (bt_j->Number - 1) / 2;
          }
          double norm_factor = System.Box.Volume / (shell * pairs * step.used);
          if (opt.axis.v[0] != -1) { // 2D: Volume/Length = in-plane area
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
  FreeArrND(opt.bt_pair);
  free(opt.bt);
  FreeSystem(&System);
  //}}}

  return 0;
}
