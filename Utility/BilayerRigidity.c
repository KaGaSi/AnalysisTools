#include "../src/AnalysisTools.h"
#include <gsl/gsl_multifit.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "BilayerRigidity computes the bending rigidity Kc*A_l (Kc times area per "
  "lipid, in DPD units) from the splay distribution P(si), where "
  "si = (nvec_j - nvec_i) / |r_ij| is the splay vector between neighbouring "
  "molecules. nvec is the orientation unit vector from the (averaged) "
  "tail-end position to the head bead. A two-step Gaussian fit of "
  "-2*ln P(si) vs si gives Kc*A_l as the quadratic coefficient, with "
  "uncertainty from the GSL covariance matrix. Both si.x and si.y are "
  "binned into the same histogram. "
  "Multiple <mol> <first> <last> trios with the same molecule name are "
  "allowed: their tail-end bead positions are averaged per molecule "
  "instance (use this for multi-tail lipids). At least one trio is "
  "required at the end of the command line (bead indices are 1-indexed; "
  "first = head bead, last = tail-end bead).",

  "Usage: BilayerRigidity <input> <width> <output> [options] "
  "<mol> <first> <last> [<mol> <first> <last> ...]",
  .args = 3,
  .all = 16,
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
  {"<input>",  NULL,    "input coordinate file",                                   OPT_ARG},
  {"<width>",  NULL,    "bin width for P(si) histogram (default in propscalc: 0.01)", OPT_ARG},
  {"<output>", NULL,    "output file",                                             OPT_ARG},
  {"--joined", NULL,    "input coordinates are already joined",                    OPT_EXTRA},
  {"-d",       "<real>","max 3D pair distance in DPD units (default: 2.0)",        OPT_EXTRA},
  {"-m",       "<real>","half-range of si histogram in DPD units (default: 2.0)", OPT_EXTRA},
  {"-t",       "<real>","min |nvec_z| tilt threshold (default: 0.5)",              OPT_EXTRA},
  {NULL}
}; //}}}

// Per-molecule-type specification.
// Multiple trios with the same mol name produce one MOL_SPEC with n_tails > 1;
// tail-end positions are averaged per molecule instance. //{{{
struct MOL_SPEC {
  int  mt;
  int  head;     // 0-indexed head bead (from first trio for this mt)
  int  n_tails;
  int *tails;    // 0-indexed tail-end beads [n_tails], heap-allocated
}; //}}}

// option struct //{{{
struct OPT {
  bool   join;
  double max_dist;
  double tilt_thresh;
}; //}}}

// data passed into the MainLoopCoor callback //{{{
struct calc_data {
  struct OPT       opt;
  struct MOL_SPEC *specs;
  int              n_specs;
  int              n_mols;
  double          *hx, *hy, *hz;
  double          *nvx, *nvy, *nvz;
  double          *psi_dist;
  long             n_pairs;
  int              n_bins;
  int              maxrbin;
  double           delta_psi;
}; //}}}

// Fit y = a + b*x + c*x^2 using GSL QR least-squares. //{{{
// c_err_out receives the standard error of the quadratic coefficient (may be NULL).
// rmsd_out receives sqrt(chisq/n) (may be NULL).
// Returns false on GSL failure or n < 3.
static bool poly2_fit(int n, const double *xx, const double *yy,
                      double *a_out, double *b_out, double *c_out,
                      double *c_err_out, double *rmsd_out) {
  if (n < 3) return false;

  gsl_matrix *X   = gsl_matrix_alloc(n, 3);
  gsl_vector *y   = gsl_vector_alloc(n);
  gsl_vector *c   = gsl_vector_alloc(3);
  gsl_matrix *cov = gsl_matrix_alloc(3, 3);
  gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(n, 3);

  for (int i = 0; i < n; i++) {
    gsl_matrix_set(X, i, 0, 1.0);
    gsl_matrix_set(X, i, 1, xx[i]);
    gsl_matrix_set(X, i, 2, xx[i] * xx[i]);
    gsl_vector_set(y, i, yy[i]);
  }

  double chisq;
  int status = gsl_multifit_linear(X, y, c, cov, &chisq, ws);

  bool ok = (status == GSL_SUCCESS);
  if (ok) {
    *a_out = gsl_vector_get(c, 0);
    *b_out = gsl_vector_get(c, 1);
    *c_out = gsl_vector_get(c, 2);
    if (c_err_out) *c_err_out = sqrt(gsl_matrix_get(cov, 2, 2));
    if (rmsd_out)  *rmsd_out  = sqrt(chisq / n);
  }

  gsl_multifit_linear_free(ws);
  gsl_matrix_free(cov);
  gsl_vector_free(c);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  return ok;
} //}}}

static void Calculation(SYSTEM *System, struct calc_data *cd) {
  WrapJoinCoordinates(System, true, cd->opt.join);

  // Collect head positions and orientation unit vectors.
  // For multi-tail specs, the tail-end position is the average of all tails.
  int nh = 0;
  for (int s = 0; s < cd->n_specs; s++) {
    MOLECULETYPE *mtype = &System->MoleculeType[cd->specs[s].mt];
    int head_b  = cd->specs[s].head;
    int n_tails = cd->specs[s].n_tails;
    for (int mi = 0; mi < mtype->Number; mi++) {
      MOLECULE *mol = &System->Molecule[mtype->Index[mi]];
      vec3d hpos = System->Bead[mol->Bead[head_b]].Position;
      double tx = 0, ty = 0, tz = 0;
      for (int k = 0; k < n_tails; k++) {
        vec3d tp = System->Bead[mol->Bead[cd->specs[s].tails[k]]].Position;
        tx += tp.x; ty += tp.y; tz += tp.z;
      }
      tx /= n_tails; ty /= n_tails; tz /= n_tails;
      double ox = hpos.x - tx;
      double oy = hpos.y - ty;
      double oz = hpos.z - tz;
      double len = sqrt(ox*ox + oy*oy + oz*oz);
      cd->hx[nh] = hpos.x;
      cd->hy[nh] = hpos.y;
      cd->hz[nh] = hpos.z;
      if (len > 1e-10) {
        cd->nvx[nh] = ox / len;
        cd->nvy[nh] = oy / len;
        cd->nvz[nh] = oz / len;
      } else {
        cd->nvx[nh] = 0.0;
        cd->nvy[nh] = 0.0;
        cd->nvz[nh] = 1.0;
      }
      nh++;
    }
  }

  double Lx = System->Box.Length.x;
  double Ly = System->Box.Length.y;
  double Lz = System->Box.Length.z;
  double max_dist    = cd->opt.max_dist;
  double tilt_thresh = cd->opt.tilt_thresh;
  double dp          = cd->delta_psi;
  int    maxrbin     = cd->maxrbin;

  for (int i = 0; i < nh - 1; i++) {
    if (fabs(cd->nvz[i]) < tilt_thresh) continue;
    for (int j = i + 1; j < nh; j++) {
      if (fabs(cd->nvz[j]) < tilt_thresh) continue;
      double dx = cd->hx[j] - cd->hx[i];
      double dy = cd->hy[j] - cd->hy[i];
      double dz = cd->hz[j] - cd->hz[i];
      dx -= round(dx / Lx) * Lx;
      dy -= round(dy / Ly) * Ly;
      dz -= round(dz / Lz) * Lz;
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      if (r > max_dist || r < 1e-10) continue;
      double six = (cd->nvx[j] - cd->nvx[i]) / r;
      double siy = (cd->nvy[j] - cd->nvy[i]) / r;
      long bx = lround(six / dp);
      long by = lround(siy / dp);
      if (labs(bx) <= (long)maxrbin) cd->psi_dist[(int)bx + maxrbin]++;
      if (labs(by) <= (long)maxrbin) cd->psi_dist[(int)by + maxrbin]++;
      cd->n_pairs++;
    }
  }
}

static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  (void)step;
  Calculation(System, (struct calc_data *)userdata);
}

int main(int argc, char *argv[]) {

  // mandatory positional arguments //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  int count = 0;

  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) exit(1);

  double delta_psi;
  if (!IsPosRealNumber(argv[++count], &delta_psi)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  }

  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  COMMON_OPT commons = CommonOptions(argc, argv, in); //}}}

  // --joined
  bool join = !BoolOption(argc, argv, "--joined");

  // -d: max pair distance
  double max_dist = 2.0;
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-d") == 0) {
      double v;
      if (!IsPosRealNumber(argv[i + 1], &v)) {
        err_msg("argument to -d must be a positive real number");
        PrintErrorOption("-d");
        exit(1);
      }
      max_dist = v;
      break;
    }
  }

  // -m: half-range of si histogram
  double max_si = 2.0;
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-m") == 0) {
      double v;
      if (!IsPosRealNumber(argv[i + 1], &v)) {
        err_msg("argument to -m must be a positive real number");
        PrintErrorOption("-m");
        exit(1);
      }
      max_si = v;
      break;
    }
  }

  // -t: tilt threshold for |nvec_z|
  double tilt_thresh = 0.5;
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-t") == 0) {
      double v;
      if (!IsPosRealNumber(argv[i + 1], &v)) {
        err_msg("argument to -t must be a positive real number");
        PrintErrorOption("-t");
        exit(1);
      }
      tilt_thresh = v;
      break;
    }
  }

  if (!commons.silent) PrintCommand(stdout, argc, argv);

  SYSTEM System = ReadStructure(in, false);

  // parse positional trios from end of argv //{{{
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
  int n_raw = (argc - trio_start) / 3;
  if (n_raw < 1) {
    err_msg("at least one <mol> <first> <last> trio is required");
    Help(true, HelpDesc, opts);
    exit(1);
  }

  // Build MOL_SPEC array: group trios by molecule type.
  // Multiple trios for the same mol type collect multiple tail-end beads. //{{{
  struct MOL_SPEC *specs = malloc(n_raw * sizeof *specs);
  if (!specs) ErrorAlloc("specs");
  // Allocate tails conservatively: at most n_raw tails per spec
  for (int s = 0; s < n_raw; s++) {
    specs[s].tails = malloc(n_raw * sizeof *specs[s].tails);
    if (!specs[s].tails) ErrorAlloc("specs[].tails");
    specs[s].n_tails = 0;
  }
  int n_specs = 0;

  for (int t = 0; t < n_raw; t++) {
    int base = trio_start + 3 * t;
    int mt = FindMoleculeName(argv[base], System);
    if (mt < 0) { ErrorMoleculeType(argv[base], System); exit(1); }
    long v1, v2;
    IsNaturalNumber(argv[base + 1], &v1);
    IsNaturalNumber(argv[base + 2], &v2);
    int first = (int)v1 - 1;
    int last  = (int)v2 - 1;
    int nb = System.MoleculeType[mt].nBeads;
    if (first >= nb || last >= nb) {
      char msg[LINE];
      snprintf(msg, LINE, "bead indices for '%s' out of range [1, %d]", argv[base], nb);
      err_msg(msg);
      exit(1);
    }
    if (first == last) {
      err_msg("<first> and <last> must specify different beads");
      exit(1);
    }

    // Find existing spec for this mt, or create a new one
    int found = -1;
    for (int s = 0; s < n_specs; s++) {
      if (specs[s].mt == mt) { found = s; break; }
    }
    if (found < 0) {
      specs[n_specs].mt      = mt;
      specs[n_specs].head    = first;
      specs[n_specs].n_tails = 0;
      found = n_specs++;
    } else if (specs[found].head != first) {
      char msg[LINE];
      snprintf(msg, LINE,
               "conflicting head bead for '%s': %d vs %d (all trios for the "
               "same molecule must share the same <first> bead)",
               argv[base], specs[found].head + 1, first + 1);
      err_msg(msg);
      exit(1);
    }
    specs[found].tails[specs[found].n_tails++] = last;
  } //}}}

  if (commons.verbose) VerboseOutput(System);

  BOX *box = &System.Box;
  if (box->Volume == -1) {
    err_msg("missing box dimensions");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  }

  int n_mols = 0;
  for (int s = 0; s < n_specs; s++)
    n_mols += System.MoleculeType[specs[s].mt].Number;

  int maxrbin = (int)round(max_si / delta_psi);
  int n_bins  = 2 * maxrbin + 1;

  double *hx  = malloc(n_mols * sizeof *hx);
  double *hy  = malloc(n_mols * sizeof *hy);
  double *hz  = malloc(n_mols * sizeof *hz);
  double *nvx = malloc(n_mols * sizeof *nvx);
  double *nvy = malloc(n_mols * sizeof *nvy);
  double *nvz = malloc(n_mols * sizeof *nvz);
  if (!hx || !hy || !hz || !nvx || !nvy || !nvz)
    ErrorAlloc("BilayerRigidity temp arrays");

  double *psi_dist = calloc(n_bins, sizeof *psi_dist);
  if (!psi_dist) ErrorAlloc("psi_dist");

  struct OPT opt = {
    .join        = join,
    .max_dist    = max_dist,
    .tilt_thresh = tilt_thresh,
  };
  struct calc_data cd = {
    .opt       = opt,
    .specs     = specs,
    .n_specs   = n_specs,
    .n_mols    = n_mols,
    .hx = hx, .hy = hy, .hz = hz,
    .nvx = nvx, .nvy = nvy, .nvz = nvz,
    .psi_dist  = psi_dist,
    .n_pairs   = 0,
    .n_bins    = n_bins,
    .maxrbin   = maxrbin,
    .delta_psi = delta_psi,
  };

  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &cd);

  // Normalise P(si) //{{{
  double norm = 0.0;
  for (int k = 0; k < n_bins; k++) norm += psi_dist[k];
  if (norm > 0.0)
    for (int k = 0; k < n_bins; k++) psi_dist[k] /= norm; //}}}

  // Two-step Gaussian fit via GSL //{{{
  double kcl = 0.0, kcl_err = 0.0, rmsd2 = 0.0;
  if (norm > 0.0) {
    double pmax = 0.0;
    for (int k = 0; k < n_bins; k++)
      if (psi_dist[k] > pmax) pmax = psi_dist[k];
    double cutval = pmax / 10.0;

    int xmin = -maxrbin, xmax = maxrbin;
    for (int i = -maxrbin; i <= 0; i++) {
      if (psi_dist[i + maxrbin]  < cutval) xmin =  i + 1;
      if (psi_dist[-i + maxrbin] < cutval) xmax = -i - 1;
    }

    int nfit1 = xmax - xmin + 1;
    double *xx1 = malloc(nfit1 * sizeof *xx1);
    double *yy1 = malloc(nfit1 * sizeof *yy1);
    if (!xx1 || !yy1) ErrorAlloc("fit1 arrays");
    int j1 = 0;
    for (int i = xmin; i <= xmax; i++) {
      double p = psi_dist[i + maxrbin];
      if (p <= 0.0) continue;
      xx1[j1] = i * delta_psi;
      yy1[j1] = -log(p);
      j1++;
    }
    double a1, b1, c1;
    double sigma = 0.0;
    if (poly2_fit(j1, xx1, yy1, &a1, &b1, &c1, NULL, NULL) && c1 > 0.0)
      sigma = sqrt(0.5 / c1);
    free(xx1); free(yy1);

    if (sigma > 0.0) {
      int xrange = (int)(sigma / delta_psi) + 1;
      if (xrange > maxrbin) xrange = maxrbin;
      int nfit2 = 2 * xrange + 1;
      double *xx2 = malloc(nfit2 * sizeof *xx2);
      double *yy2 = malloc(nfit2 * sizeof *yy2);
      if (!xx2 || !yy2) ErrorAlloc("fit2 arrays");
      int j2 = 0;
      for (int i = -xrange; i <= xrange; i++) {
        double p = psi_dist[i + maxrbin];
        if (p <= 0.0) continue;
        xx2[j2] = i * delta_psi;
        yy2[j2] = -2.0 * log(p);
        j2++;
      }
      double a2, b2, c2;
      if (poly2_fit(j2, xx2, yy2, &a2, &b2, &c2, &kcl_err, &rmsd2) && c2 > 0.0)
        kcl = c2;
      free(xx2); free(yy2);
    }
  } //}}}

  // write output //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# (1) si; (2) P(si)\n");

  ArrNDd *data = CreateArr2Dd(n_bins + 2, 2);
  if (!data) ErrorAlloc("data");
  for (int k = 0; k < n_bins; k++) {
    int bin_idx = k - maxrbin;
    SetArr2D(data, k, 0, bin_idx * delta_psi);
    SetArr2D(data, k, 1, psi_dist[k]);
  }
  ComputeColumnWidths(n_bins, 2, data, 6);
  PrintDataAll(fw, n_bins, 2, data);
  FreeArrND(data);

  fprintf(fw, "# Kc*A_l: %.6f +/- %.6f  rmsd: %.6f\n", kcl, kcl_err, rmsd2);
  fclose(fw); //}}}

  for (int s = 0; s < n_specs; s++) free(specs[s].tails);
  free(specs);
  free(hx); free(hy); free(hz);
  free(nvx); free(nvy); free(nvz);
  free(psi_dist);
  FreeSystem(&System);

  return 0;
}
