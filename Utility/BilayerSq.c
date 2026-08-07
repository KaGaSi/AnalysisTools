#include "../src/AnalysisTools.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "BilayerSq calculates the 2D in-plane static structure factor and "
  "its radially-averaged 1D form for each bilayer leaflet, then averages "
  "both leaflets. Leaflets are assigned by comparing each head-bead z-position "
  "to the global midplane (z-COM of all specified head beads). "
  "Output: <output>-2d.txt (qx, qy, S) and <output>-1d.txt (q, S). "
  "At least one <mol> <bead> pair is required after the mandatory arguments "
  "(bead index is 1-indexed).",

  "Usage: BilayerSq <input> <delta_q> <output> <mol> <bead> "
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
  {"<delta_q>", nullptr, "q-point spacing in DPD units", OPT_ARG},
  {"<output>", nullptr, "output base name (appends -2d.txt and -1d.txt)",
    OPT_ARG},
  {"<mol>", nullptr, "molecule name", OPT_ARG},
  {"<bead>", nullptr, "1-indexed head bead within molecule", OPT_ARG},
  {"--joined", nullptr, "input coordinates are already joined", OPT_EXTRA},
  {"-q", "<int>",  "max bins per axis (default: 50)", OPT_EXTRA},
  {"-a", "<axis>", "bilayer normal axis: x, y, or z (default: z)", OPT_EXTRA},
  {nullptr}
}; //}}}

// head-group bead specification //{{{
struct HEADSPEC {
  int mt;
  int bead; // 0-indexed within molecule
}; //}}}

// option struct //{{{
struct OPT {
  bool join;   // --joined
  int maxqbin; // -q
  int axis;    // -a: bilayer normal axis (0=x, 1=y, 2=z)
}; //}}}

// data passed into the MainLoopCoor callback //{{{
struct calc_data {
  struct OPT       opt;
  struct HEADSPEC *specs;
  int n_specs, n_heads;
  // per-frame temporaries (pre-allocated)
  vec3d *h;
  int *leaflet;
  double *sf_re, *sf_im; // cos/sin sums over one leaflet, size n_qpts
  // accumulated S(q): indexed [leaflet][iq_offset * (maxqbin+1) + jq]
  double *sqxy[2];
  long frame_count;
  // parameters
  int n_qpts; // (2*maxqbin+1) * (maxqbin+1)
  double delta_q;
}; //}}}

static void Calculation(SYSTEM *System, struct calc_data *cd) { //{{{
  WrapJoinCoordinates(System, true, cd->opt.join);

  int a  = cd->opt.axis;
  int p1 = (a + 1) % 3;
  int p2 = (a + 2) % 3;

  // collect head positions and assign leaflets by global midplane //{{{
  double z_sum = 0.0;
  int nh = 0;
  for (int i = 0; i < cd->n_specs; i++) {
    MOLECULETYPE *mtype = &System->MoleculeType[cd->specs[i].mt];
    for (int j = 0; j < mtype->Number; j++) {
      MOLECULE *mol = &System->Molecule[mtype->Index[j]];
      int id = mol->Bead[cd->specs[i].bead];
      vec3d p = System->Bead[id].Position;
      cd->h[nh] = p;
      z_sum += p.v[a];
      nh++;
    }
  }
  if (nh == 0) {
    return;
  } //}}}

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

  int mqb  = cd->opt.maxqbin;
  int njq  = mqb + 1; // number of jq values: 0..maxqbin
  double dq = cd->delta_q;

  for (int l = 0; l < 2; l++) {
    if (n_leaf[l] == 0) {
      continue;
    }

    memset(cd->sf_re, 0, cd->n_qpts * sizeof *cd->sf_re);
    memset(cd->sf_im, 0, cd->n_qpts * sizeof *cd->sf_im);

    // accumulate cos/sin structure factor sums  //{{{
    // outer loop over particles for cache-friendly access to sqxy.
    for (int i = 0; i < nh; i++) {
      if (cd->leaflet[i] != l) {
        continue;
      }
      double xi = cd->h[i].v[p1],
             yi = cd->h[i].v[p2];
      for (int iq = -mqb; iq <= mqb; iq++) {
        double qx_xi = iq * dq * xi;
        int ioff = (iq + mqb) * njq;
        for (int jq = 0; jq <= mqb; jq++) {
          double arg = qx_xi + jq * dq * yi;
          int id = ioff + jq;
          cd->sf_re[id] += cos(arg);
          cd->sf_im[id] += sin(arg);
        }
      }
    } //}}}

    // S(q) = |sfsum|^2 / N
    double inv_n = 1.0 / n_leaf[l];
    for (int k = 0; k < cd->n_qpts; k++) {
      double s = Square(cd->sf_re[k]) + Square(cd->sf_im[k]);
      cd->sqxy[l][k] += s * inv_n;
    }
  }

  cd->frame_count++;
} //}}}

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
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }

  double delta_q;
  if (!IsPosRealNumber(argv[++count], &delta_q)) {
    ErrorNaN("<delta_q>");
    Help(true, HelpDesc, opts);
    exit(1);
  }

  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  COMMON_OPT commons = CommonOptions(argc, argv, in); //}}}

  // --joined option (opt.join == true -> needs joining)
  opt.join = !BoolOption(argc, argv, "--joined");

  opt.maxqbin = 50;
  OneNumberOption(argc, argv, "-q", &opt.maxqbin, 'i');
  if (opt.maxqbin <= 0) {
    err_msg("requires a positive integer");
    PrintErrorOption("-q");
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
    if (i + 1 >= argc) {
      break;
    }
    long v;
    if (FindMoleculeName(argv[i], System) >= 0 &&
        IsNaturalNumber(argv[i+1], &v)) {
      n_specs++;
    } else {
      break;
    }
  }
  if (n_specs < 1) {
    err_msg("at least one <mol> <bead> pair is required");
    PrintError();
    Help(true, HelpDesc, opts);
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
    int bead = v - 1;
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
    specs[s].mt = mt;
    specs[s].bead = bead;
    n_heads += System.MoleculeType[mt].Number;
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  if (System.Box.Volume == -1) { //{{{
    err_msg("missing box dimensions");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  } //}}}

  int n_qpts = (2 * opt.maxqbin + 1) * (opt.maxqbin + 1);

  // pre-allocate per-frame arrays
  vec3d *h = malloc(n_heads * sizeof *h);
  int *leaflet = malloc(n_heads * sizeof *leaflet);
  double *sf_re = malloc(n_qpts  * sizeof *sf_re),
         *sf_im = malloc(n_qpts  * sizeof *sf_im);
  double *sq0 = calloc(n_qpts, sizeof *sq0),
         *sq1 = calloc(n_qpts, sizeof *sq1);
  if (!h || !leaflet || !sf_re || !sf_im || !sq0 || !sq1) {
    ErrorAlloc("h/leaflet/sf_re/sf_im/sq0/sq1");
  }

  struct calc_data cd = {
    .opt = opt,
    .specs = specs,
    .n_specs = n_specs,
    .n_heads = n_heads,
    .h = h,
    .leaflet = leaflet,
    .sf_re = sf_re,
    .sf_im = sf_im,
    .sqxy = { sq0, sq1 },
    .frame_count = 0,
    .n_qpts = n_qpts,
    .delta_q = delta_q,
  };

  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &cd);

  if (cd.frame_count == 0) {
    err_msg("no frames processed");
    PrintError();
    exit(1);
  }

  // average over frames and both leaflets: S_avg=(S_lower+S_upper)/(2*nframes)
  double scale = 1.0 / (2.0 * cd.frame_count);
  double *sq_avg = malloc(n_qpts * sizeof *sq_avg);
  if (!sq_avg) {
    ErrorAlloc("sq_avg");
  }
  for (int i = 0; i < n_qpts; i++) {
    sq_avg[i] = (sq0[i] + sq1[i]) * scale;
  }

  int njq = opt.maxqbin + 1;

  // write 2D S(q) //{{{
  char f2d[LINE+16];
  snprintf(f2d, sizeof f2d, "%s-2d.txt", fout);
  FILE *fw = PrintBylineOpenFile(f2d, argc, argv);
  fprintf(fw, "# (1) qx; (2) qy; (3) S(qx,qy)\n");
  for (int iq = -opt.maxqbin; iq <= opt.maxqbin; iq++) {
    double qx = iq * delta_q;
    int ioff = (iq + opt.maxqbin) * njq;
    for (int jq = 0; jq <= opt.maxqbin; jq++) {
      double qy = jq * delta_q;
      double s = sq_avg[ioff+jq];
      fprintf(fw, " %12.6f %12.6f %12.6f\n", qx, qy, s);
      // reconstruct negative-qy half by symmetry S(qx,qy) = S(-qx,-qy)
      if (jq > 0) {
        fprintf(fw, " %12.6f %12.6f %12.6f\n", -qx, -qy, s);
      }
    }
  }
  fclose(fw); //}}}

  // write 1D radially-averaged S(q) //{{{
  int n_1d = (int)(sqrt(2.0) * opt.maxqbin) + 2;
  double *sum1d = calloc(n_1d, sizeof *sum1d);
  int *cnt1d = calloc(n_1d, sizeof *cnt1d);
  if (!sum1d || !cnt1d) {
    ErrorAlloc("sum1d/cnt1d");
  }

  for (int i = -opt.maxqbin; i <= opt.maxqbin; i++) {
    int ioff = (i + opt.maxqbin) * njq;
    for (int j = 0; j <= opt.maxqbin; j++) {
      int k = (int)sqrt((Square(i) + Square(j)));
      if (k < n_1d) {
        sum1d[k] += sq_avg[ioff+j];
        cnt1d[k]++;
      }
    }
  }

  char f1d[LINE+16];
  snprintf(f1d, sizeof f1d, "%s-1d.txt", fout);
  fw = PrintBylineOpenFile(f1d, argc, argv);
  fprintf(fw, "# (1) q; (2) S(q)\n");

  // count non-empty bins for ArrNDd; skip k=0 (q=0 forward-scattering, S=N)
  int n_out = 0;
  for (int k = 1; k < n_1d; k++) {
    if (cnt1d[k] > 0) {
      n_out++;
    }
  }

  ArrNDd *data = CreateArr2Dd(n_out + 2, 2);
  if (!data) {
    ErrorAlloc("1d data");
  }
  int row = 0;
  for (int k = 1; k < n_1d; k++) {
    if (cnt1d[k] == 0) {
      continue;
    }
    SetArr2D(data, row, 0, (k + 0.5) * delta_q);
    SetArr2D(data, row, 1, sum1d[k] / cnt1d[k]);
    row++;
  }
  ComputeColumnWidths(n_out, 2, data, 6);
  PrintDataAll(fw, n_out, 2, data);
  FreeArrND(data);
  fclose(fw); //}}}

  // free all arrays //{{{
  free(specs);
  free(h);
  free(leaflet);
  free(sf_re);
  free(sf_im);
  free(sq0);
  free(sq1);
  free(sq_avg);
  free(sum1d);
  free(cnt1d);
  FreeSystem(&System); //}}}

  return 0;
}
