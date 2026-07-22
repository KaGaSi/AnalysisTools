#include "../src/AnalysisTools.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "BilayerInPlane calculates the in-plane 2D radial distribution function "
  "g2D(r) and the hexatic order parameter |Phi6| for specified head-group "
  "beads in each leaflet of a bilayer. Leaflets are assigned by comparing "
  "each head-bead z-position to the global midplane (z-COM of all head "
  "beads). Output: <output>-rdf.txt (per-leaflet g2D) and "
  "<output>-phi6.txt (per-leaflet |Phi6| distribution + average). "
  "At least one <mol> <bead> pair is required." ,

  "Usage: BilayerInPlane <input> <width> <output> <mol> <bead> "
  "[<mol> <bead> ...] [options]",
  .args = 5,
  .all = 17,
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
  {"<width>", nullptr, "bin width for g2D(r) and |Phi6| distribution", OPT_ARG},
  {"<output>", nullptr, "output base name (appends -rdf.txt and -phi6.txt)",
    OPT_ARG},
  {"<mol>", nullptr, "molecule name", OPT_ARG},
  {"<bead>", nullptr, "1-indexed head bead within molecule", OPT_ARG},
  {"--joined", nullptr, "input coordinates are already joined", OPT_EXTRA},
  {"-r", "<real>", "neighbour cutoff for |Phi6| in DPD units (default: 1.5)",
    OPT_EXTRA},
  {"-a", "<axis>", "bilayer normal axis: x, y, or z (default: z)", OPT_EXTRA},
  {nullptr}
}; //}}}

// head-group bead specification (one per mol/bead pair on argv) //{{{
struct HEADSPEC {
  int mt;   // molecule type index
  int bead; // 0-indexed bead position within molecule
}; //}}}

// option struct //{{{
struct OPT {
  bool join;    // --join
  double r_cut; // Phi6 neighbour cutoff
  int axis;     // -a  bilayer normal axis (0=x, 1=y, 2=z)
}; //}}}

// data passed into the MainLoopCoor callback //{{{
struct calc_data {
  struct OPT opt;
  struct HEADSPEC *specs;
  int n_specs,
      n_heads;   // total head beads (sum of mtype.Number for all specs)
  // per-frame temp arrays (pre-allocated, reset each frame)
  vec3d *h;
  int *leaflet;
  double *phi6r, *phi6i;
  int *nb_cnt;
  // accumulated per-leaflet data
  double *rdf[2];       // g2D histogram counts [rdf_bins]
  long n2_sum[2];       // sum N^2 per leaflet over frames (RDF normalisation)
  double *phi6_dist[2]; // |Phi6| histogram counts [phi6_bins]
  double  phi6_sum[2];  // sum |Phi6| per leaflet
  long phi6_cnt[2];     // count of valid |Phi6| values per leaflet
  // parameters
  int rdf_bins,
      phi6_bins;
  double rdf_max;
  double width;
}; //}}}

static void Calculation(SYSTEM *System, struct calc_data *cd) {
  WrapJoinCoordinates(System, true, cd->opt.join);

  int a  = cd->opt.axis;
  int p1 = (a + 1) % 3;
  int p2 = (a + 2) % 3;
  vec3d *length = &System->Box.Length;

  // Collect head bead positions and sum for midplane
  double z_sum = 0.0;
  int nh = 0;
  for (int s = 0; s < cd->n_specs; s++) {
    MOLECULETYPE *mtype = &System->MoleculeType[cd->specs[s].mt];
    for (int mi = 0; mi < mtype->Number; mi++) {
      MOLECULE *mol = &System->Molecule[mtype->Index[mi]];
      vec3d p = System->Bead[mol->Bead[cd->specs[s].bead]].Position;
      cd->h[nh] = p;
      z_sum += p.v[a];
      nh++;
    }
  }
  if (nh == 0) {
    return;
  }

  // Assign leaflets by global midplane
  double midplane = z_sum / nh;
  int n_leaf[2] = {0, 0};
  for (int i = 0; i < nh; i++) {
    if (cd->h[i].v[a] > midplane) {
      cd->leaflet[i] = 1;
    } else {
      cd->leaflet[i] = 0;
    }
    n_leaf[cd->leaflet[i]]++;
  }

  // Reset per-frame Phi6 accumulators
  memset(cd->phi6r,  0, nh * sizeof *cd->phi6r);
  memset(cd->phi6i,  0, nh * sizeof *cd->phi6i);
  memset(cd->nb_cnt, 0, nh * sizeof *cd->nb_cnt);

  // Main pair loop: accumulate g2D and Phi6 simultaneously
  for (int i = 0; i < nh; i++) {
    for (int j = (i + 1); j < nh; j++) {
      if (cd->leaflet[i] != cd->leaflet[j]) {
        continue;
      }

      double dx = cd->h[j].v[p1] - cd->h[i].v[p1];
      double dy = cd->h[j].v[p2] - cd->h[i].v[p2];
      dx -= round(dx / length->v[p1]) * length->v[p1]; // 2D minimum image
      dy -= round(dy / length->v[p2]) * length->v[p2];
      double r = sqrt(Square(dx) + Square(dy));
      int l = cd->leaflet[i];

      // g2D: bin unique pair
      if (r > 0 && r < cd->rdf_max) {
        int k = (int)(r / cd->width);
        if (k < cd->rdf_bins) {
          cd->rdf[l][k]++;
        }
      }

      // Phi6: accumulate if within neighbour cutoff
      if (r > 0 && r < cd->opt.r_cut) {
        double theta = atan2(dy, dx);
        double cos6 = cos(6.0 * theta);
        double sin6 = sin(6.0 * theta);
        // exp(6i*theta_ji) == exp(6i*theta_ij) because theta_ji = theta_ij + pi
        // and 6*pi is a multiple of 2*pi, so the contribution is symmetric
        cd->phi6r[i] += cos6;
        cd->phi6i[i] += sin6;
        cd->nb_cnt[i]++;
        cd->phi6r[j] += cos6;
        cd->phi6i[j] += sin6;
        cd->nb_cnt[j]++;
      }
    }
  }

  // Accumulate Phi6 distribution
  for (int i = 0; i < nh; i++) {
    if (cd->nb_cnt[i] == 0) {
      continue;
    }
    int l = cd->leaflet[i];
    double mag = sqrt(Square(cd->phi6r[i]) + Square(cd->phi6i[i]))
                 / cd->nb_cnt[i];
    int k = mag / cd->width;
    if (k >= cd->phi6_bins) {
      k = cd->phi6_bins - 1;
    }
    cd->phi6_dist[l][k]++;
    cd->phi6_sum[l] += mag;
    cd->phi6_cnt[l]++;
  }

  for (int l = 0; l < 2; l++) {
    cd->n2_sum[l] += Square(n_leaf[l]);
  }
}

static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  (void)step;
  Calculation(System, (struct calc_data *)userdata);
}

int main(int argc, char *argv[]) {

  // mandatory positional arguments //{{{
  OptionCheck(argc, argv, false, HelpDesc, opts);
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

  opt.r_cut = 1.5;
  if (OneNumberOption(argc, argv, "-r", &opt.r_cut, 'd') && opt.r_cut <= 0) {
    err_msg("requires a positive real number");
    PrintErrorOption("-r");
    exit(1);
  }

  opt.axis = 2;
  char axis_arg[LINE];
  if (FileOption(argc, argv, "-a", axis_arg)) {
    if (axis_arg[0] == 'x') {
      opt.axis = 0;
    } else if (axis_arg[0] == 'y') {
      opt.axis = 1;
    } else if (axis_arg[0] == 'z') {
      opt.axis = 2;
    } else {
      err_msg("requires argument 'x', 'y', or 'z'");
      PrintErrorOption("-a");
      exit(1);
    }
  }

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);

  // parse <mol> <bead> [<mol> <bead> ...] //{{{
  int spec_start = count + 1;
  int n_specs = 0;
  for (int i = spec_start; i < argc && argv[i][0] != '-'; i += 2) {
    if ((i + 1) >= argc) {
      break;
    }
    long v;
    if (FindMoleculeName(argv[i], System) >= 0 &&
        IsNaturalNumber(argv[i + 1], &v)) {
      n_specs++;
    } else {
      break;
    }
  }
  if (n_specs < 1) {
    err_msg("at least one <mol> <bead> pair is required");
    PrintError();
    exit(1);
  }

  struct HEADSPEC *specs = malloc(n_specs * sizeof *specs);
  if (!specs) {
    ErrorAlloc("specs");
  }

  int n_heads = 0;
  for (int s = 0; s < n_specs; s++) {
    int base = spec_start + 2 * s;
    int mt = FindMoleculeName(argv[base], System);
    if (mt < 0) {
      ErrorMoleculeType(argv[base], System);
      exit(1);
    }
    long v;
    IsNaturalNumber(argv[base+1], &v);
    int bead = v - 1; // 1-indexed → 0-indexed
    int nb = System.MoleculeType[mt].nBeads;
    if (bead < 0 || bead >= nb) {
      if (snprintf(ERROR_MSG, LINE, "wrong bead index %s%d%s for %s%s%s ; "
          "must be [1, %d]", ErrYellow(), bead + 1, ErrRed(),
          ErrYellow(), argv[base], ErrRed(), nb) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    specs[s].mt   = mt;
    specs[s].bead = bead;
    n_heads += System.MoleculeType[mt].Number;
  } //}}}

  if (commons.verbose) VerboseOutput(System);

  BOX *box = &System.Box;
  if (box->Volume == -1) {
    err_msg("missing box dimensions");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  }

  int p1 = (opt.axis + 1) % 3;
  int p2 = (opt.axis + 2) % 3;
  double l1 = box->Length.v[p1],
         l2 = box->Length.v[p2];
  double lmin = l2;
  if (l1 < l2) {
    lmin = l1;
  }
  double rdf_max = lmin / 2;
  int rdf_bins  = rdf_max / width;
  int phi6_bins = ceil(1.0 / width);
  if (rdf_bins  < 1) {
    rdf_bins  = 1;
  }
  if (phi6_bins < 1) {
    phi6_bins = 1;
  }

  // pre-allocate per-frame temp arrays
  vec3d  *h = malloc(n_heads * sizeof *h);
  int *leaflet = malloc(n_heads * sizeof *leaflet);
  double *phi6r = malloc(n_heads * sizeof *phi6r);
  double *phi6i = malloc(n_heads * sizeof *phi6i);
  int *nb_cnt = malloc(n_heads * sizeof *nb_cnt);
  if (!h || !leaflet || !phi6r || !phi6i || !nb_cnt) {
    ErrorAlloc("BilayerInPlane temp arrays");
  }

  // per-leaflet output arrays
  double *rdf0 = calloc(rdf_bins, sizeof *rdf0);
  double *rdf1 = calloc(rdf_bins, sizeof *rdf1);
  double *phi0 = calloc(phi6_bins, sizeof *phi0);
  double *phi1 = calloc(phi6_bins, sizeof *phi1);
  if (!rdf0 || !rdf1 || !phi0 || !phi1) {
    ErrorAlloc("BilayerInPlane output arrays");
  }

  struct calc_data cd = {
    .opt = opt,
    .specs = specs,
    .n_specs = n_specs,
    .n_heads = n_heads,
    .h = h,
    .leaflet = leaflet,
    .phi6r = phi6r,
    .phi6i = phi6i,
    .nb_cnt = nb_cnt,
    .rdf = { rdf0, rdf1 },
    .n2_sum = { 0, 0 },
    .phi6_dist = { phi0, phi1 },
    .phi6_sum = { 0.0, 0.0 },
    .phi6_cnt = { 0, 0 },
    .rdf_bins = rdf_bins,
    .phi6_bins = phi6_bins,
    .rdf_max = rdf_max,
    .width = width,
  };

  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &cd);

  // write g2D output //{{{
  char frdf[LINE+16];
  snprintf(frdf, sizeof frdf, "%s-rdf.txt", fout);
  FILE *fw = PrintBylineOpenFile(frdf, argc, argv);
  fprintf(fw, "# (1) r; (2) lower_leaflet; (3) upper_leaflet\n");

  // g(r) = count * Lx*Ly / (n2_sum * pi * r * dr)
  // (counting unique pairs i<j; factor of 2 absorbed into denominator)
  ArrNDd *data = CreateArr2Dd(rdf_bins + 2, 3);
  if (!data) ErrorAlloc("rdf data");
  for (int k = 0; k < rdf_bins; k++) {
    double r = (k + 0.5) * width;
    SetArr2D(data, k, 0, r);
    for (int l = 0; l < 2; l++) {
      double g = 0;
      if (cd.n2_sum[l] > 0) {
        g = cd.rdf[l][k] * l1 * l2 / (cd.n2_sum[l] * M_PI * r * width);
      } else {
        g = 0;
      }
      SetArr2D(data, k, l + 1, g);
    }
  }
  ComputeColumnWidths(rdf_bins, 3, data, 6);
  PrintDataAll(fw, rdf_bins, 3, data);
  FreeArrND(data);
  fclose(fw); //}}}

  // write Phi6 output //{{{
  char fphi6[LINE+16];
  snprintf(fphi6, sizeof fphi6, "%s-phi6.txt", fout);
  fw = PrintBylineOpenFile(fphi6, argc, argv);
  fprintf(fw, "# (1) |Phi6|; (2) lower_leaflet; (3) upper_leaflet\n");

  double phi6_total[2] = {0.0, 0.0};
  for (int l = 0; l < 2; l++)
    for (int k = 0; k < phi6_bins; k++)
      phi6_total[l] += cd.phi6_dist[l][k];

  data = CreateArr2Dd(phi6_bins + 2, 3);
  if (!data) {
    ErrorAlloc("phi6 data");
  }
  for (int k = 0; k < phi6_bins; k++) {
    double center = (k + 0.5) * width;
    SetArr2D(data, k, 0, center);
    for (int l = 0; l < 2; l++) {
      double norm;
      if (phi6_total[l] > 0) {
        norm = cd.phi6_dist[l][k] / phi6_total[l];
      } else {
        norm = 0;
      }
      SetArr2D(data, k, l + 1, norm);
    }
  }
  ComputeColumnWidths(phi6_bins, 3, data, 6);
  PrintDataAll(fw, phi6_bins, 3, data);
  FreeArrND(data);

  fprintf(fw, "# Average |Phi6|:");
  for (int l = 0; l < 2; l++) {
    double avg;
    if (cd.phi6_cnt[l] > 0) {
      avg = cd.phi6_sum[l] / cd.phi6_cnt[l];
    } else {
      avg = 0;
    }
    if (l == 0) {
      fprintf(fw, " lower");
    } else {
      fprintf(fw, " upper");
    }
    fprintf(fw, "=%.6f", avg);
  }
  putc('\n', fw);
  fclose(fw); //}}}

  free(specs);
  free(h);
  free(leaflet);
  free(phi6r);
  free(phi6i);
  free(nb_cnt);
  free(rdf0);
  free(rdf1);
  free(phi0);
  free(phi1);
  FreeSystem(&System);

  return 0;
}
