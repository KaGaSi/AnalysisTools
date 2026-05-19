#include "../src/AnalysisTools.h"

// TODO: -a option for axis (leave -z as default)
// TODO: change BilayerSq ... [options] <mol> <bead> [<mol> <bead>]
//       into BilayerSq ... <mol> <bead> [<mol> <bead>] [options]
//       also changes .args - 3->5
// TODO: why did I use those *{x,y,z} vars instead of vec3d?

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
  "At least one <mol> <bead> pair is required at the end of the command line "
  "(bead index is 1-indexed).",

  "Usage: BilayerSq <input> <delta_q> <output> [options] "
  "<mol> <bead> [<mol> <bead> ...]",
  .args = 3,
  .all = 14,
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
  {"<input>",   NULL,    "input coordinate file",                          OPT_ARG},
  {"<delta_q>", NULL,    "q-point spacing in DPD units    ",               OPT_ARG},
  {"<output>",  NULL,    "output base name (appends -2d.txt and -1d.txt)", OPT_ARG},
  {"--joined",  NULL,    "input coordinates are already joined",           OPT_EXTRA},
  {"-q",        "<int>", "max bins per axis (default: 50)",                OPT_EXTRA},
  {NULL}
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
}; //}}}

// data passed into the MainLoopCoor callback //{{{
struct calc_data {
  struct OPT       opt;
  struct HEADSPEC *specs;
  int n_specs, n_heads;
  // per-frame temporaries (pre-allocated)
  double *hx, *hy, *hz;
  int *leaflet;
  double *sfx, *sfy; // cos/sin sums over one leaflet, size n_qpts
  // accumulated S(q): indexed [leaflet][iq_offset * (maxqbin+1) + jq]
  double *sqxy[2];
  long frame_count;
  // parameters
  int n_qpts; // (2*maxqbin+1) * (maxqbin+1)
  double delta_q;
}; //}}}

static void Calculation(SYSTEM *System, struct calc_data *cd) { //{{{
  WrapJoinCoordinates(System, true, cd->opt.join);

  // collect head positions and assign leaflets by global midplane //{{{
  double z_sum = 0.0;
  int nh = 0;
  for (int i = 0; i < cd->n_specs; i++) {
    MOLECULETYPE *mtype = &System->MoleculeType[cd->specs[i].mt];
    for (int j = 0; j < mtype->Number; j++) {
      MOLECULE *mol = &System->Molecule[mtype->Index[j]];
      int id = mol->Bead[cd->specs[i].bead];
      vec3d p = System->Bead[id].Position;
      cd->hx[nh] = p.x;
      cd->hy[nh] = p.y;
      cd->hz[nh] = p.z;
      z_sum += p.z;
      nh++;
    }
  }
  if (nh == 0) {
    return;
  } //}}}

  double midplane = z_sum / nh;
  int n_leaf[2] = {0, 0};
  for (int i = 0; i < nh; i++) {
    if (cd->hz[i] > midplane) {
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

    memset(cd->sfx, 0, cd->n_qpts * sizeof *cd->sfx);
    memset(cd->sfy, 0, cd->n_qpts * sizeof *cd->sfy);

    // accumulate cos/sin structure factor sums  //{{{
    // outer loop over particles for cache-friendly access to sqxy.
    for (int i = 0; i < nh; i++) {
      if (cd->leaflet[i] != l) {
        continue;
      }
      double xi = cd->hx[i],
             yi = cd->hy[i];
      for (int iq = -mqb; iq <= mqb; iq++) {
        double qx_xi = iq * dq * xi;
        int ioff = (iq + mqb) * njq;
        for (int jq = 0; jq <= mqb; jq++) {
          double arg = qx_xi + jq * dq * yi;
          int id = ioff + jq;
          cd->sfx[id] += cos(arg);
          cd->sfy[id] += sin(arg);
        }
      }
    } //}}}

    // S(q) = |sfsum|^2 / N
    double inv_n = 1.0 / n_leaf[l];
    for (int k = 0; k < cd->n_qpts; k++) {
      double s = cd->sfx[k] * cd->sfx[k] + cd->sfy[k] * cd->sfy[k];
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
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;

  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) exit(1);

  double delta_q;
  if (!IsPosRealNumber(argv[++count], &delta_q)) {
    ErrorNaN("<delta_q>");
  }

  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  COMMON_OPT commons = CommonOptions(argc, argv, in); //}}}

  // --joined option (opt.join == true -> needs joining)
  opt.join = !BoolOption(argc, argv, "--joined");

  // -q option - maxqbin
  // TODO: why not OneNumberOption()?
  opt.maxqbin = 50;
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-q") == 0) {
      long v;
      if (!IsNaturalNumber(argv[i+1], &v)) {
        err_msg("argument must be a positive integer");
        PrintErrorOption("-q");
        exit(1);
      }
      opt.maxqbin = v;
      break;
    }
  }

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);

  // parse <mol> <bead> pairs from end of argv //{{{
  int spec_start = argc;
  for (int i = argc - 2; i > count; i -= 2) {
    long v;
    if (FindMoleculeName(argv[i], System) >= 0 &&
        IsNaturalNumber(argv[i+1], &v)) {
      spec_start = i;
    } else {
      break;
    }
  }
  int n_specs = (argc - spec_start) / 2;
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
      ErrorMoleculeType(argv[base], System); exit(1);
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
  // TODO: why not vec3d?
  double *hx = malloc(n_heads * sizeof *hx),
         *hy = malloc(n_heads * sizeof *hy),
         *hz = malloc(n_heads * sizeof *hz);
  int *leaflet = malloc(n_heads * sizeof *leaflet);
  double *sfx = malloc(n_qpts  * sizeof *sfx),
         *sfy = malloc(n_qpts  * sizeof *sfy);
  double *sq0 = calloc(n_qpts, sizeof *sq0),
         *sq1 = calloc(n_qpts, sizeof *sq1);
  if (!hx || !hy || !hz || !leaflet || !sfx || !sfy || !sq0 || !sq1) {
    ErrorAlloc("hx/hy/hz/leaflet/sfx/sfy/sq0/sq1");
  }

  struct calc_data cd = {
    .opt = opt,
    .specs = specs,
    .n_specs = n_specs,
    .n_heads = n_heads,
    .hx = hx,
    .hy = hy,
    .hz = hz,
    .leaflet = leaflet,
    .sfx = sfx,
    .sfy = sfy,
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
  for (int i = 0; i < n_qpts; i++)
    sq_avg[i] = (sq0[i] + sq1[i]) * scale;

  int njq = opt.maxqbin + 1;

  // write 2D S(q) //{{{
  char f2d[LINE + 16];
  snprintf(f2d, sizeof f2d, "%s-2d.txt", fout);
  FILE *fw = PrintBylineOpenFile(f2d, argc, argv);
  fprintf(fw, "# (1) qx; (2) qy; (3) S(qx,qy)\n");
  for (int iq = -opt.maxqbin; iq <= opt.maxqbin; iq++) {
    double qx = iq * delta_q;
    int ioff = (iq + opt.maxqbin) * njq;
    for (int jq = 0; jq <= opt.maxqbin; jq++) {
      double qy = jq * delta_q;
      double s  = sq_avg[ioff + jq];
      fprintf(fw, " %12.6f %12.6f %12.6f\n",  qx,  qy, s);
      // reconstruct negative-qy half by symmetry S(qx,qy) = S(-qx,-qy)
      if (jq > 0)
        fprintf(fw, " %12.6f %12.6f %12.6f\n", -qx, -qy, s);
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

  char f1d[LINE + 16];
  snprintf(f1d, sizeof f1d, "%s-1d.txt", fout);
  fw = PrintBylineOpenFile(f1d, argc, argv);
  fprintf(fw, "# (1) q; (2) S(q)\n");

  // count non-empty bins for ArrNDd; skip k=0 (q=0 forward-scattering, S=N)
  int n_out = 0;
  for (int k = 1; k < n_1d; k++) if (cnt1d[k] > 0) n_out++;

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
  free(hx);
  free(hy);
  free(hz);
  free(leaflet);
  free(sfx);
  free(sfy);
  free(sq0);
  free(sq1);
  free(sq_avg);
  free(sum1d);
  free(cnt1d);
  FreeSystem(&System); //}}}

  return 0;
}
