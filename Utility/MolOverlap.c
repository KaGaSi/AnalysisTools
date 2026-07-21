#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "MolOverlap calculates tail interdigitation in a bilayer. "
  "For each xy-grid column the midplane is the z centre-of-mass of the "
  "specified tail beads; molecules are assigned to leaflets by comparison "
  "with that midplane. Overlap is defined as "
  "lower_leaflet_inner_z - upper_leaflet_inner_z "
  "(positive = overlap, negative = gap). "
  "The distribution and average overlap are written to <output>. "
  "At least one <mol> <first> <last> trio must be given at the end of the "
  "command line (1-indexed bead positions, multiple trios allowed).",

  "Usage: MolOverlap <input> <width> <output> [options] "
  "<mol> <first> <last> [<mol> <first> <last> ...]",
  .args = 3,
  .all = 15,
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
  {"<width>", nullptr, "distribution bin width", OPT_ARG},
  {"<output>", nullptr, "output file", OPT_ARG},
  {"--joined", nullptr, "input coordinates are already joined", OPT_EXTRA},
  {"-g", "<real>", "xy grid cell width (default: 2.5)", OPT_EXTRA},
  {"-r", "<real>", "distribution half-range (default: 5.0)", OPT_EXTRA},
  {nullptr}
}; //}}}

// molecule-type + tail-bead range (0-indexed internally) //{{{
struct TRIO {
  int mt;    // molecule type index
  int first; // 0-indexed first tail bead position within molecule
  int last;  // 0-indexed last tail bead position within molecule
}; //}}}

// option struct //{{{
struct OPT {
  bool   join;
  double gridw;
}; //}}}

// data passed into the MainLoopCoor callback //{{{
struct calc_data {
  struct OPT    opt;
  struct TRIO  *trios;
  int           n_trios;
  double       *dist;       // distribution histogram [bins]
  double        total_ov;   // running sum of overlap values
  long          total_cnt;  // count of valid per-column measurements
  int           bins;       // number of distribution bins
  int           nx, ny;     // fixed grid dimensions
  double        width;      // bin width
  double        range;      // half-range
}; //}}}

static void Calculation(SYSTEM *System, struct calc_data *cd) {
  BOX *box = &System->Box;
  WrapJoinCoordinates(System, true, cd->opt.join);

  int ncells = cd->nx * cd->ny;

  double *mid_sum  = calloc(ncells, sizeof *mid_sum);
  int    *mid_cnt  = calloc(ncells, sizeof *mid_cnt);
  double *upper    = malloc(ncells * sizeof *upper);
  double *lower    = malloc(ncells * sizeof *lower);
  if (!mid_sum || !mid_cnt || !upper || !lower) {
    ErrorAlloc("MolOverlap grid arrays");
  }
  for (int i = 0; i < ncells; i++) {
    upper[i] =  HIGHNUM;
    lower[i] = -HIGHNUM;
  }

  // Pass 1: accumulate per-column z-COM for the midplane
  for (int t = 0; t < cd->n_trios; t++) {
    MOLECULETYPE *mtype = &System->MoleculeType[cd->trios[t].mt];
    int nb = cd->trios[t].last - cd->trios[t].first + 1;
    for (int mi = 0; mi < mtype->Number; mi++) {
      MOLECULE *mol = &System->Molecule[mtype->Index[mi]];
      double sx = 0, sy = 0, sz = 0;
      for (int b = cd->trios[t].first; b <= cd->trios[t].last; b++) {
        vec3d p = System->Bead[mol->Bead[b]].Position;
        sx += p.x; sy += p.y; sz += p.z;
      }
      int cx = (int)(fmod(sx / nb, box->Length.x) / cd->opt.gridw);
      int cy = (int)(fmod(sy / nb, box->Length.y) / cd->opt.gridw);
      if (cx < 0) cx = 0;
      if (cx >= cd->nx) cx = cd->nx - 1;
      if (cy < 0) cy = 0;
      if (cy >= cd->ny) cy = cd->ny - 1;
      mid_sum[cx * cd->ny + cy] += sz / nb;
      mid_cnt[cx * cd->ny + cy]++;
    }
  }

  // Finalise per-column midplane
  double *midplane = malloc(ncells * sizeof *midplane);
  bool   *has_mid  = calloc(ncells, sizeof *has_mid);
  if (!midplane || !has_mid) {
    ErrorAlloc("MolOverlap midplane");
  }
  for (int c = 0; c < ncells; c++) {
    if (mid_cnt[c] > 0) {
      midplane[c] = mid_sum[c] / mid_cnt[c];
      has_mid[c]  = true;
    }
  }
  free(mid_sum);
  free(mid_cnt);

  // Pass 2: assign each molecule to a leaflet, track inner-bead z
  for (int t = 0; t < cd->n_trios; t++) {
    MOLECULETYPE *mtype = &System->MoleculeType[cd->trios[t].mt];
    int nb = cd->trios[t].last - cd->trios[t].first + 1;
    for (int mi = 0; mi < mtype->Number; mi++) {
      MOLECULE *mol = &System->Molecule[mtype->Index[mi]];
      double sx = 0, sy = 0, sz = 0;
      double z_min =  HIGHNUM;
      double z_max = -HIGHNUM;
      for (int b = cd->trios[t].first; b <= cd->trios[t].last; b++) {
        vec3d p = System->Bead[mol->Bead[b]].Position;
        sx += p.x; sy += p.y; sz += p.z;
        if (p.z < z_min) z_min = p.z;
        if (p.z > z_max) z_max = p.z;
      }
      double comx = sx / nb, comy = sy / nb, comz = sz / nb;
      int cx = (int)(fmod(comx, box->Length.x) / cd->opt.gridw);
      int cy = (int)(fmod(comy, box->Length.y) / cd->opt.gridw);
      if (cx < 0) cx = 0;
      if (cx >= cd->nx) cx = cd->nx - 1;
      if (cy < 0) cy = 0;
      if (cy >= cd->ny) cy = cd->ny - 1;
      int cell = cx * cd->ny + cy;
      if (!has_mid[cell]) continue;

      if (comz > midplane[cell]) {
        // upper leaflet: innermost bead has smallest z
        if (z_min < upper[cell]) upper[cell] = z_min;
      } else {
        // lower leaflet: innermost bead has largest z
        if (z_max > lower[cell]) lower[cell] = z_max;
      }
    }
  }

  // Accumulate into distribution
  for (int c = 0; c < ncells; c++) {
    if (!has_mid[c] || upper[c] >= HIGHNUM || lower[c] <= -HIGHNUM) continue;
    double ov = lower[c] - upper[c];
    int k = (int)((ov + cd->range) / cd->width);
    if (k < 0) k = 0;
    if (k >= cd->bins) k = cd->bins - 1;
    cd->dist[k]++;
    cd->total_ov  += ov;
    cd->total_cnt++;
  }

  free(midplane);
  free(has_mid);
  free(upper);
  free(lower);
}

static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  (void)step;
  Calculation(System, (struct calc_data *)userdata);
}

int main(int argc, char *argv[]) {

  // mandatory positional arguments //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;

  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) exit(1);

  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  }

  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  COMMON_OPT commons = CommonOptions(argc, argv, in); //}}}

  // --joined option (opt.join == true -> needs joining)
  opt.join = !BoolOption(argc, argv, "--joined");

  // -g: xy grid cell width
  opt.gridw = 2.5;
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-g") == 0) {
      double v;
      if (!IsPosRealNumber(argv[i + 1], &v)) {
        err_msg("requires a positive real number");
        PrintErrorOption("-g");
        exit(1);
      }
      opt.gridw = v;
      break;
    }
  }

  // -r: distribution half-range
  double range = 5.0;
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-r") == 0) {
      double v;
      if (!IsPosRealNumber(argv[i + 1], &v)) {
        err_msg("argument to -r must be a positive real number");
        PrintErrorOption("-r");
        exit(1);
      }
      range = v;
      break;
    }
  }

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);

  // parse mandatory positional trios from end of argv //{{{
  // scan backward in groups of 3: <mol> <first> <last>
  int trio_start = argc;
  for (int i = argc - 3; i > count; i -= 3) {
    long v1, v2;
    if (FindMoleculeName(argv[i], System) >= 0 &&
        IsNaturalNumber(argv[i + 1], &v1) &&
        IsNaturalNumber(argv[i + 2], &v2)) {
      trio_start = i;
    } else {
      break;
    }
  }
  int n_trios = (argc - trio_start) / 3;
  if (n_trios < 1) {
    err_msg("at least one <mol> <first> <last> trio is required");
    Help(true, HelpDesc, opts);
    exit(1);
  }

  struct TRIO *trios = malloc(n_trios * sizeof *trios);
  if (!trios) ErrorAlloc("trios");

  for (int t = 0; t < n_trios; t++) {
    int base = trio_start + 3 * t;
    int mt = FindMoleculeName(argv[base], System);
    if (mt < 0) {
      ErrorMoleculeType(argv[base], System);
      exit(1);
    }
    long v1, v2;
    IsNaturalNumber(argv[base + 1], &v1);
    IsNaturalNumber(argv[base + 2], &v2);
    int first = (int)v1 - 1; // convert 1-indexed → 0-indexed
    int last  = (int)v2 - 1;
    int nb = System.MoleculeType[mt].nBeads;
    if (first >= nb || last >= nb) {
      char msg[LINE];
      snprintf(msg, LINE,
               "bead indices for '%s' out of range [1, %d]", argv[base], nb);
      err_msg(msg);
      exit(1);
    }
    if (first > last) {
      err_msg("<first> must not exceed <last>");
      exit(1);
    }
    trios[t].mt    = mt;
    trios[t].first = first;
    trios[t].last  = last;
  } //}}}

  if (commons.verbose) VerboseOutput(System);

  BOX *box = &System.Box;
  if (box->Volume == -1) {
    err_msg("missing box dimensions");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  }

  int nx = box->Length.x / opt.gridw;
  int ny = box->Length.y / opt.gridw;
  if (nx < 1) {
    nx = 1;
  }
  if (ny < 1) {
    ny = 1;
  }

  int bins = 2 * range / width;
  if (bins < 1) {
    bins = 1;
  }

  double *dist = calloc(bins, sizeof *dist);
  if (!dist) {
    ErrorAlloc("dist");
  }

  struct calc_data cd = {
    .opt       = opt,
    .trios     = trios,
    .n_trios   = n_trios,
    .dist      = dist,
    .total_ov  = 0.0,
    .total_cnt = 0,
    .bins      = bins,
    .nx        = nx,
    .ny        = ny,
    .width     = width,
    .range     = range,
  };

  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &cd);

  // write output //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# (1) overlap; (2) normalised count\n");

  double total = 0;
  for (int k = 0; k < bins; k++) total += cd.dist[k];

  ArrNDd *data = CreateArr2Dd(bins + 2, 2);
  if (!data) ErrorAlloc("data");
  for (int k = 0; k < bins; k++) {
    double ov_center = -range + width * (k + 0.5);
    double norm = total > 0 ? cd.dist[k] / total : 0;
    SetArr2D(data, k, 0, ov_center);
    SetArr2D(data, k, 1, norm);
  }
  ComputeColumnWidths(bins, 2, data, 6);
  PrintDataAll(fw, bins, 2, data);
  FreeArrND(data);

  double avg = cd.total_cnt > 0 ? cd.total_ov / cd.total_cnt : 0;
  fprintf(fw, "# Average overlap: %.6f\n", avg);
  fclose(fw); //}}}

  free(trios);
  free(dist);
  FreeSystem(&System);

  return 0;
}
