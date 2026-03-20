#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Structure calculates distribution of local 2D bond order orientation "
  "parameters (BOOPs) for 3- to 8-fold symmetry.",

  "Usage: Structure <input> <output> <double> [options]",
  .args = 3, // number of mandatory arguments
  .all = 12, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<width>", NULL, "width of a bin", OPT_ARG},
  {"<output>", NULL, "output file", OPT_ARG},
  {"-pb", "<file>", "save per-bead BOOPs from the last step (automatic ending -<symmetry>.txt)", OPT_EXTRA},
  {NULL}
}; //}}}

// calculate bond orientation order parameter for a single bead //{{{
void ComputeBOOP(SYSTEM System, int n, int n_sym, int sym[n_sym], ArrNDd *boop) {
  int max_neigh = sym[n_sym-1];
  int nearest[max_neigh]; // nearest neighbour id
  double min_sqdist[max_neigh]; // nearest neighbour's square distance
  for (int i = 0; i < max_neigh; i++) {
    nearest[i] = -1;
    min_sqdist[i] = HIGHNUM;
  }
  // find nearest neighbours
  BEAD *b = &System.Bead[n];
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    if (id == n) {
      continue;
    }
    vec3d d = Distance(b->Position, System.Bead[id].Position, System.Box.Length);
    double sq_r = SqVectLength(d);
    for (int j = 0; j < max_neigh; j++) {
      if (sq_r < min_sqdist[j]) {
        for (int k = (max_neigh - 1); k > j; k--) {
          min_sqdist[k] = min_sqdist[k-1];
          nearest[k] = nearest[k-1];
        }
        min_sqdist[j] = sq_r;
        nearest[j] = id;
        break;
      }
    }
  }
  for (int i = 0; i < max_neigh; i++) {
    if (nearest[i] == -1) {
      err_msg("Not enough neighbours! Should never happen!!");
      PrintError();
    }
  }

  // compute angles for all max_neigh neighbours for use in all the symmetries
  double theta[max_neigh];
  for (int j = 0; j < max_neigh; j++) {
    vec3d d = Distance(System.Bead[nearest[j]].Position, b->Position,
                       System.Box.Length);
    theta[j] = atan2(d.v[1], d.v[0]);
  }
  // calculate boop
  for (int i = 0; i < n_sym; i++) {
    double q_real = 0, q_imag = 0;
    for (int j = 0; j < sym[i]; j++) {
      q_real += cos(sym[i] * theta[j]);
      q_imag += sin(sym[i] * theta[j]);
    }
    q_real /= sym[i];
    q_imag /= sym[i];
    SetArr3D(boop, n, i, 0, q_real);
    SetArr3D(boop, n, i, 1, q_imag);
  }
} //}}}

// structure for options //{{{
struct OPT {
  char per_bead_file[LINE]; // -pb option
}; //}}}

// g_n pair callback - called once per pair (i<j) by TraversePairs //{{{
struct gn_args {
  ArrNDd *boop;
  ArrNDd *g_n;
  long int *g_n_counts;
  int n_sym;
  double r_max;
  double dr;
};
static void GnPair(int i, int j, const SYSTEM System, void *ud) {
  struct gn_args *p = (struct gn_args*)ud;
  int id_i = System.BeadCoor[i];
  int id_j = System.BeadCoor[j];
  vec3d dist = Distance(System.Bead[id_i].Position, System.Bead[id_j].Position,
                        System.Box.Length);
  double r_ij = VectLength(dist);
  if (r_ij >= p->r_max) {
    return;
  }
  int bin = r_ij / p->dr;
  for (int k = 0; k < p->n_sym; k++) {
    double re = GetArr3D(p->boop, id_i, k, 0) * GetArr3D(p->boop, id_j, k, 0);
    double im = GetArr3D(p->boop, id_i, k, 1) * GetArr3D(p->boop, id_j, k, 1);
    AddArr2D(p->g_n, k, bin, re + im);
  }
  p->g_n_counts[bin]++;
}
static bool AcceptAll(int i, const SYSTEM System, void *ud) {
  (void)i; (void)System; (void)ud;
  return true;
} //}}}

// per-step calculation //{{{
struct calc_data {
  ArrNDd *boop;
  ArrNDd *boop_distr_type;
  ArrNDd *g_n;
  long int *g_n_counts;
  int n_sym;
  int *sym;
  int bins;
  double width;
  double r_max;
  double dr;
};
static void Calculation(SYSTEM *System, STEP *step, void *userdata) {
  (void)step;
  struct calc_data *p = (struct calc_data*)userdata;
  COUNT *Count = &System->Count;

  FillArrND(p->boop, 0);
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System->BeadCoor[i];
    int type = System->Bead[id].Type;
    ComputeBOOP(*System, id, p->n_sym, p->sym, p->boop);
    for (int j = 0; j < p->n_sym; j++) {
      double res = sqrt(Square(GetArr3D(p->boop, id, j, 0)) +
                        Square(GetArr3D(p->boop, id, j, 1)));
      // error - should never happen but better safe than sorry //{{{
      // TODO: proper error with bead id and whatnot
      if (res < 0 || res > 1) {
        err_msg("boop must be <0,1>!");
        PrintError();
      } //}}}
      int k = res / p->width;
      if (k == p->bins) {
        k--;
      }
      AddArr3D(p->boop_distr_type, type, k, j, 1);
    }
  }

  // g_n correlation: TraversePairs visits each pair once (vs. original O(N^2) double loop)
  struct gn_args args = { p->boop, p->g_n, p->g_n_counts,
                          p->n_sym, p->r_max, p->dr };
  TraversePairs(*System, p->r_max, GnPair, &args, AcceptAll, NULL);
} //}}}

int main(int argc, char *argv[]) {

  // command line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <width> - width of the bin for boop distribution
  double width = 0;
  if (!IsRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!FileOption(argc, argv, "-pb", opt.per_bead_file)) {
    opt.per_bead_file[0] = '\0';
  }
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  BOX *Box = &System.Box;

  if (commons.verbose) {
    VerboseOutput(System);
  }

  int n_sym = 6; // number of symmetries to calculate
  int sym[n_sym];
  sym[0] = 3; // n-fold symmetry
  sym[1] = 4; //
  sym[2] = 5; //
  sym[3] = 6; //
  sym[4] = 7; //
  sym[5] = 8; //

  int bins = 1 / width;
  // TODO: this must be somehow made better...
  double r_max = Max3(Box->Length.x, Box->Length.y, Box->Length.z) / 2;
  double dr = 1;
  int bins_g_n = r_max / dr + 1;

  // arrays for boop distributions
  ArrNDd *boop_distr_type = CreateArr3Dd(Count->BeadType, bins, n_sym);
  ArrNDd *g_n = CreateArr2Dd(n_sym, bins_g_n);
  if (!boop_distr_type || !g_n) {
    err_msg("ArrNDd constructor failed (boop_distr_type/g_n)");
    PrintError();
    exit(1);
  }
  long int *g_n_counts = calloc(bins_g_n, sizeof *g_n_counts);
  ArrNDd *boop = CreateArr3Dd(Count->Bead, n_sym, 2);
  if (!boop) {
    err_msg("ArrNDd constructor failed (boop)");
    PrintError();
    exit(1);
  }
// main loop //{{{
  struct calc_data ud = { boop, boop_distr_type, g_n, g_n_counts,
                          n_sym, sym, bins, width, r_max, dr };
  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation, &ud);
  //}}}

  if (opt.per_bead_file[0] != '\0') {
    // open all per-bead files once, write all beads, close - not per bead
    FILE *pb_files[n_sym];
    for (int i = 0; i < n_sym; i++) {
      char pb_name[LINE];
      if (snprintf(pb_name, LINE, "%s-%d.txt", opt.per_bead_file, sym[i]) < 0) {
        ErrorSnprintf();
      }
      pb_files[i] = OpenFile(pb_name, "w");
    }
    for (int i = 0; i < Count->BeadCoor; i++) {
      int id = System.BeadCoor[i];
      for (int j = 0; j < n_sym; j++) {
        double res = sqrt(Square(GetArr3D(boop, id, j, 0)) +
                          Square(GetArr3D(boop, id, j, 1)));
        fprintf(pb_files[j], "%lf\n", res);
      }
    }
    for (int i = 0; i < n_sym; i++) {
      fclose(pb_files[i]);
    }
  }

  // TODO: print the per-bead type stuff
  // flattening per-type distribution array //{{{
  ArrNDd *boop_distr = CreateArr2Dd(bins, n_sym);
  if (!boop_distr) {
    err_msg("ArrNDd constructor failed (boop)");
    PrintError();
    exit(1);
  }
  // loop order: type -> bin -> sym matches boop_distr_type's [type][bin][sym] strides
  double *norm = calloc(n_sym, sizeof *norm);
  for (int k = 0; k < Count->BeadType; k++) {
    for (int i = 0; i < bins; i++) {
      for (int j = 0; j < n_sym; j++) {
        AddArr2D(boop_distr, i, j, GetArr3D(boop_distr_type, k, i, j));
      }
    }
  }
  for (int i = 0; i < bins; i++) {
    for (int j = 0; j < n_sym; j++) {
      norm[j] += GetArr2D(boop_distr, i, j);
    }
  } //}}}

  // print boop distribution //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  // print header
  count = 1;
  fprintf(fw, "# (%d) boop; symmetry:", count++);
  bool first = true;
  for (int i = 0; i < n_sym; i++) {
    if (!first) {
      fprintf(fw, ", ");
    }
    fprintf(fw, "(%d) %d-fold", count++, sym[i]);
    first = false;
  }
  putc('\n', fw);
  // print data
  // collate data //{{{
  int ncols = n_sym + 1;
  int nrows = bins;
  ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
  if (!data) {
    err_msg("ArrNDd constructor failed (data)");
    PrintError();
    exit(1);
  }
  for (int i = 0; i < bins; i++) {
    count = -1;
    SetArr2D(data, i, ++count, width * (2 * i + 1) / 2);
    for (int j = 0; j < n_sym; j++) {
      double val = GetArr2D(boop_distr, i, j) / norm[j];
      SetArr2D(data, i, ++count, val);
    }
  } //}}}
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(fw, nrows, ncols, data);
  FreeArrND(data);
  fclose(fw); //}}}

  // create file name TODO: needs something better, I guess
  // print g_n to file //{{{
  char tmp[LINE] = "",
       fout2[LINE] = "";
  s_strcpy(tmp, fout, strlen(fout) - 3); // assumes .txt ending...
  if (snprintf(fout2, LINE, "%s-corr.txt", tmp) < 0) {
    ErrorSnprintf();
  }
  // print headers
  fw = PrintBylineOpenFile(fout2, argc, argv);
  count = 1;
  fprintf(fw, "# (%d) distance; boop correlation:", count++);
  first = true;
  for (int i = 0; i < n_sym; i++) {
    if (!first) {
      fprintf(fw, ", ");
    }
    fprintf(fw, " (%d) %d-fold", count++, sym[i]);
    first = false;
  }
  putc('\n', fw);
  // collate data //{{{
  ncols = n_sym + 1;
  nrows = bins_g_n;
  data = CreateArr2Dd(nrows + 2, ncols);
  for (int i = 0; i < bins_g_n; i++) {
    count = -1;
    SetArr2D(data, i, ++count, dr * (2 * i + 1) / 2);
    for (int j = 0; j < n_sym; j++) {
      double val = 0;
      if (g_n_counts[i] > 0) {
        val = GetArr2D(g_n, j, i) / g_n_counts[i];
      }
      SetArr2D(data, i, ++count, val);
    }
  } //}}}
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(fw, nrows, ncols, data);
  FreeArrND(data);
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  FreeArrND(boop_distr_type);
  FreeArrND(g_n);
  FreeArrND(boop);
  FreeArrND(boop_distr);
  free(g_n_counts);
  free(norm);
  FreeSystem(&System);
  //}}}

  return 0;
}
