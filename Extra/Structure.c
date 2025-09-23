#include "../AnalysisTools.h"

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
    vec3 d = Distance(b->Position.v, b_i->Position.v, System.Box.Length);
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
  for (int i = 0; i < n_sym; i++) {
    double q_real = 0;
    double q_imag = 0;
    for (int j = 0; j < sym[i]; j++) {
      BEAD *b2 = &System.Bead[nearest[j]];
      vec3 d = Distance(b2->Position.v, b->Position.v, System.Box.Length);
      double theta = atan2(d.v[1], d.v[0]);
      q_real += cos(sym[i] * theta);
      q_imag += sin(sym[i] * theta);
    }
    q_real /= sym[i];
    q_imag /= sym[i];
    res[i] = sqrt(Square(q_real) + Square(q_imag));
  }
} //}}}

// // calculate bond orientation order parameter for a single bead //{{{
// // ...calculated results weighed by distance from centre particle - g(r)
// // based Gaussian thingy
// void ComputeBOOP(SYSTEM System, int n,
//                  int n_sym, int sym[n_sym], double res[n_sym]) {
//   double r_nn = 1.0; // TODO: make into opt - should be \approx 1st min in g(r)
//   double r_min = 1.2;
//   double sigma = (r_min - r_nn) / 2;
//   int max_neigh = sym[n_sym-1];
//   int nearest[max_neigh]; // nearest neighbour id
//   double min_dist[max_neigh]; // nearest neighbour's distance
//   for (int i = 0; i < max_neigh; i++) {
//     nearest[i] = -1;
//     min_dist[i] = HIGHNUM;
//   }
//   // find nearest neighbours
//   BEAD *b = &System.Bead[n];
//   for (int i = 0; i < System.Count.BeadCoor; i++) {
//     int id = System.BeadCoor[i];
//     if (id == n) {
//       continue;
//     }
//     BEAD *b_i = &System.Bead[id];
//     vec3 d = Distance(b->Position.v, b_i->Position.v, System.Box.Length);
//     double r = VectLength(d);
//     for (int j = 0; j < max_neigh; j++) {
//       if (r < min_dist[j]) {
//         for (int k = (max_neigh - 1); k > j; k--) {
//           min_dist[k] = min_dist[k-1];
//           nearest[k] = nearest[k-1];
//         }
//         min_dist[j] = r;
//         nearest[j] = id;
//         break;
//       }
//     }
//   }
//   for (int i = 0; i < max_neigh; i++) {
//     if (nearest[i] == -1) {
//       err_msg("Huh? Not enough neighbours! Should never happen!!");
//       PrintError();
//     }
//   }
//
//   // calculate boop
//   for (int i = 0; i < n_sym; i++) {
//     double q_real = 0;
//     double q_imag = 0;
//     double wsum = 0;
//     for (int j = 0; j < max_neigh; j++) {
//       BEAD *b2 = &System.Bead[nearest[j]];
//       vec3 d = Distance(b2->Position.v, b->Position.v, System.Box.Length);
//       double theta = atan2(d.v[1], d.v[0]);
//       double w = exp(- Square((VectLength(d) - r_nn) / sigma));
//       q_real += w * cos(sym[i] * theta);
//       q_imag += w * sin(sym[i] * theta);
//       wsum += w;
//       printf("%d: %lf\n", j, wsum);
//     }
//     printf("%lf\n", wsum);
//     if (wsum < 1e-12) {
//       res[i] = 0;
//     } else {
//       q_real /= wsum;
//       q_imag /= wsum;
//       res[i] = sqrt(Square(q_real) + Square(q_imag));
//     }
//   }
// } //}}}

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
          int k =  res[j] / width; // edge case for res[] == 1 (put into highest bin)
          if (k == bins) {
            k--; // edge case for res[j] -> 1
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
