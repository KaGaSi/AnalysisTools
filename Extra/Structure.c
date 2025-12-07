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
  double min_dist[max_neigh]; // nearest neighbour's distance
  for (int i = 0; i < max_neigh; i++) {
    nearest[i] = -1;
    min_dist[i] = HIGHNUM;
  }
  // find nearest neighbours
  BEAD *b = &System.Bead[n];
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    if (id == n) {
      continue;
    }
    BEAD *b_i = &System.Bead[id];
    vec3d d = Distance(b->Position.v, b_i->Position.v, System.Box.Length);
    double r = VectLength(d);
    for (int j = 0; j < max_neigh; j++) {
      if (r < min_dist[j]) {
        for (int k = (max_neigh - 1); k > j; k--) {
          min_dist[k] = min_dist[k-1];
          nearest[k] = nearest[k-1];
        }
        min_dist[j] = r;
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

  // calculate boop
  for (int i = 0; i < n_sym; i++) {
    double q_real = 0;
    double q_imag = 0;
    for (int j = 0; j < sym[i]; j++) {
      BEAD *b2 = &System.Bead[nearest[j]];
      vec3d d = Distance(b2->Position.v, b->Position.v, System.Box.Length);
      double theta = atan2(d.v[1], d.v[0]);
      q_real += cos(sym[i] * theta);
      q_imag += sin(sym[i] * theta);
    }
    q_real /= sym[i];
    q_imag /= sym[i];
    SetArr3D(boop, n, i, 0, q_real);
    SetArr3D(boop, n, i, 1, q_imag);
  }
} //}}}

// Help() //{{{
void Help_old(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
} //}}}

// structure for options //{{{
struct OPT {
  char per_bead_file[LINE]; // -pb option
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
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
  double r_max = Max3(Box->Length[0], Box->Length[1], Box->Length[2]) / 2;
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

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  ArrNDd *boop = CreateArr3Dd(Count->Bead, n_sym, 2);
  if (!boop) {
    err_msg("ArrNDd constructor failed (boop)");
    PrintError();
    exit(1);
  }
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;

      FillArrND(boop, 0);
      // go over all beads in the coordinate file
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i]; // bead index
        int type = System.Bead[id].Type;
        ComputeBOOP(System, id, n_sym, sym, boop);
        for (int j = 0; j < n_sym; j++) {
          double res = sqrt(Square(GetArr3D(boop, id, j, 0)) +
                            Square(GetArr3D(boop, id, j, 1)));
          // error - should never happen but better safe than sorry //{{{
          // TODO: proper error with bead id and whatnot
          if (res < 0 || res > 1) {
            err_msg("boop must be <0,1>!");
            PrintError();
          } //}}}
          int k = res / width; // edge case for res == 1 (put into highest bin)
          if (k == bins) {
            k--; // edge case for res[j] -> 1
          }
          AddArr3D(boop_distr_type, type, k, j, res);
        }
      }

      for (int i = 0; i < Count->BeadCoor; i++) {
        for (int j = 0; j < Count->BeadCoor; j++) {
          if (i == j) {
            continue;
          }
          int id_i = System.BeadCoor[i];
          int id_j = System.BeadCoor[j];
          vec3d *pos_i = &System.Bead[id_i].Position;
          vec3d *pos_j = &System.Bead[id_j].Position;
          vec3d dist = Distance(pos_i->v, pos_j->v, System.Box.Length);
          double r_ij = VectLength(dist);
          if (r_ij >= r_max) {
            continue;
          }
          int bin = r_ij / dr;

          // boop are complex numbers
          // boop_i.real * boop_j.real + boop_i.imag * boop_j.imag
          for (int k = 0; k < n_sym; k++) {
            double re = GetArr3D(boop, id_i, k, 0) * GetArr3D(boop, id_j, k, 0);
            double im = GetArr3D(boop, id_i, k, 1) * GetArr3D(boop, id_j, k, 1);
            double correlation = re + im;
            AddArr2D(g_n, k, bin, correlation);
          }
          g_n_counts[bin]++;
        }
      } //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!commons.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  if (opt.per_bead_file[0] != '\0') {
    // new file
    for (int i = 0; i < n_sym; i++) {
      char fout[LINE] = "";
      if (snprintf(fout, LINE, "%s-%d.txt", opt.per_bead_file, sym[i]) < 0) {
        ErrorSnprintf();
      }
      FILE *fw = OpenFile(fout, "w");
      fclose(fw);
    }
    // print per-bead values
    for (int i = 0; i < Count->BeadCoor; i++) {
      int id = System.BeadCoor[i]; // bead index
      for (int j = 0; j < n_sym; j++) {
        double res = sqrt(Square(GetArr3D(boop, id, j, 0)) +
                          Square(GetArr3D(boop, id, j, 1)));

        char fout[LINE] = "";
        if (snprintf(fout, LINE, "%s-%d.txt", opt.per_bead_file, sym[j]) < 0) {
          ErrorSnprintf();
        }
        FILE *fw = OpenFile(fout, "a");
        fprintf(fw, "%lf\n", res);
        fclose(fw);
      }
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
  // normalisation factor for distribution
  double *norm = calloc(n_sym, sizeof *norm);
  for (int i = 0; i < bins; i++) {
    for (int j = 0; j < n_sym; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        AddArr2D(boop_distr, i, j, GetArr3D(boop_distr_type, k, i, j));
      }
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
      if (g_n_counts[i] > 0) {
        double val = GetArr2D(g_n, j, i) / g_n_counts[i];
        SetArr2D(data, i, ++count, val);
      } else {
        SetArr2D(data, i, ++count, 0);
      }
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
