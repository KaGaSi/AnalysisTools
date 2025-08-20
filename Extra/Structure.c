#include "../AnalysisTools.h"

// Syrex-created; not yet tested whatsoever...
// Finds the cutoff index for neighbours using the largest relative gap
// Inputs:
//   min_dist[]: sorted array of neighbour distances (ascending), length max_neigh
//   max_neigh: maximum number of neighbours considered
//   gap_threshold: minimum relative gap ratio to accept (e.g. 0.2)
//   fallback_neigh: fallback fixed neighbour count if no gap found
// Returns:
//   neighbour_count: number of neighbours selected for the first shell
int FindNeighbourShell(double min_dist[], int max_neigh,
                      double gap_threshold, int fallback_neigh) {
  int cutoff = fallback_neigh;
  double max_ratio = 0.0;
  for (int i = 0; i < max_neigh - 1; i++) {
    double gap = min_dist[i+1] - min_dist[i];
    if (min_dist[i] < 1e-12) continue; // avoid divide by zero
    double ratio = gap / min_dist[i];
    if (ratio > max_ratio) {
      max_ratio = ratio;
      cutoff = i + 1;
    }
  }
  if (max_ratio < gap_threshold) {
    cutoff = fallback_neigh; // no significant gap found
  }
  if (cutoff > max_neigh) cutoff = max_neigh; // safety clamp
  return cutoff;
}

// calculate bond orientation order parameter for a single bead //{{{
void ComputeBOOP(SYSTEM System, int n,
                 int n_sym, int sym[n_sym], double res[n_sym]) {
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
    double d[3];
    Distance(b->Position, b_i->Position, System.Box.Length, d);
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
      err_msg("Huh? Not enough neighbours! Should never happen!!");
      PrintError();
    }
  }

  // calculate boop
  double q_real[n_sym];
  double q_imag[n_sym];
  for (int i = 0; i < n_sym; i++) {
    q_real[i] = 0;
    q_imag[i] = 0;
  }

  for (int i = 0; i < n_sym; i++) {
    // Syrex-suggested 'Adaptive neighbour shell detection via gap in sorted distances'
    // double gap_threshold = 0.2;    // tweak this parameter
    // int fallback_neigh = sym[i];   // number of neighbours expected for symmetry i
    // int neigh_count = FindNeighbourShell(min_dist, max_neigh, gap_threshold, fallback_neigh);
    // ...if used, change next loop (go to neigh_count if found or
    // fallback_neigh if not)
    for (int j = 0; j < sym[i]; j++) {
      BEAD *b2 = &System.Bead[nearest[j]];
      double d[3];
      Distance(b2->Position, b->Position, System.Box.Length, d);
      double theta = atan2(d[1], d[0]);
      q_real[i] += cos(sym[i] * theta);
      q_imag[i] += sin(sym[i] * theta);
    }
    q_real[i] /= sym[i];
    q_imag[i] /= sym[i];
    res[i] = sqrt(Square(q_real[i]) + Square(q_imag[i]));
  }
} //}}}

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Calculate distribution of bond order orientation parameters for 3- to \
8-fold symmetry.\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <output> <double> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a bin\n");
  fprintf(ptr, "<output>            output file\n");
  fprintf(ptr, "[options]\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 0, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // mandatory options //{{{
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
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  }
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  //}}}

  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  if (opt->c.verbose) {
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

  // arrays for boop distribution //{{{
  double ***boop_distr_type = calloc(Count->BeadType, sizeof *boop_distr_type);
  for (int i = 0; i < Count->BeadType; i++) {
    boop_distr_type[i] = calloc(bins, sizeof *boop_distr_type[i]);
    for (int j = 0; j <  bins; j++) {
      boop_distr_type[i][j] = calloc(n_sym, sizeof *boop_distr_type[i][j]);
    }
  } //}}}

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      // go over all beads in the coordinate file
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i]; // bead index
        int type = System.Bead[id].Type;
        double res[n_sym];
        ComputeBOOP(System, id, n_sym, sym, res);
        for (int j = 0; j < n_sym; j++) {
          // error - should never happen but better safe than sorry //{{{
          // TODO: proper error with bead id and whatnot
          if (res[j] > 1 || res[j] < 0) {
            err_msg("boop must be <0,1>!");
            PrintError();
          } //}}}
          int k = bins - 1; // edge case for res[] == 1 (put into highest bin)
          if (res[j] < 1) {
            // if (j == 0) {
            //   printf("%lf %lf ...", res[j], width);
            // }
            k = res[j] / width;
            // if (j == 0) {
            //   printf(" %d\n", k);
            // }
          }
          boop_distr_type[type][k][j] += res[j];
        }
      }
      //}}}
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
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // flattening per-type distribution array
  // TODO: print the per-bead type stuff
  double **boop_distr = calloc(bins, sizeof *boop_distr);
  for (int i = 0; i < bins; i++) {
    boop_distr[i] = calloc(n_sym, sizeof *boop_distr[i]);
  }
  // normalisation factor for distribution
  double *norm = calloc(n_sym, sizeof *norm);
  for (int i = 0; i < bins; i++) {
    for (int j = 0; j < n_sym; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        boop_distr[i][j] += boop_distr_type[k][i][j];
      }
      norm[j] += boop_distr[i][j];
    }
  }

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
  // determine width of each column //{{{
  int columns = n_sym + 1;
  int digits[columns][2];
  InitInt2DArray((int *)digits, columns, 2, 0);
  double *data[bins]; // array for data
  for (int i = 0; i < bins; i++) {
    data[i] = calloc(columns, sizeof data[i]);
    count = -1;
    data[i][++count] = width * (2 * i + 1) / 2;
    for (int j = 0; j < n_sym; j++) {
      data[i][++count] = boop_distr[i][j] / norm[j];
    }
  }
  FillMaxDigits(columns, bins, data, digits); //}}}
  for (int i = 0; i < bins; i++) {
    for (int col = 0; col < columns; col++) {
      Fprintf1(fw, data[i][col], digits[col]);
    }
    putc('\n', fw);
    free(data[i]);
  }
  fclose(fw);

  // free memory - to make valgrind happy //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = 0; j < bins; j++) {
      free(boop_distr_type[i][j]);
      free(boop_distr[j]);
    }
    free(boop_distr_type[i]);
  }
  free(boop_distr_type);
  free(boop_distr);
  free(norm);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
