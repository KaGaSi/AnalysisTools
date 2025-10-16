#include "../AnalysisTools.h"

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
    size_t index[3] = {n, i, 0};
    SetArrND(boop, index, q_real);
    index[2] = 1;
    SetArrND(boop, index, q_imag);
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
  BOX *Box = &System.Box;

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
  // TODO: this must be somehow made better...
  double r_max = Max3(Box->Length[0], Box->Length[1], Box->Length[2]) / 2;
  double dr = 1;
  int bins_g_n = r_max / dr + 1;

  // arrays for boop distributions
  size_t shape3D[3] = {Count->BeadType, bins, n_sym};
  ArrNDd *boop_distr_type = CreateArrNDd(3, shape3D);
  ArrNDd *g_n = CreateArr2Dd(n_sym, bins_g_n);
  long int *g_n_counts = calloc(bins_g_n, sizeof *g_n_counts);

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  shape3D[0] = Count->Bead;
  shape3D[1] = n_sym;
  shape3D[2] = 2;
  ArrNDd *boop = CreateArrNDd(3, shape3D);
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

      FillArrND(boop, 0);
      // go over all beads in the coordinate file
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i]; // bead index
        int type = System.Bead[id].Type;
        // double res[n_sym];
        ComputeBOOP(System, id, n_sym, sym, boop);
        for (int j = 0; j < n_sym; j++) {
          size_t id_0[3] = {id, j, 0};
          size_t id_1[3] = {id, j, 1};
          double res = sqrt(Square(GetArrND(boop, id_0)) +
                            Square(GetArrND(boop, id_1)));
          // error - should never happen but better safe than sorry //{{{
          // TODO: proper error with bead id and whatnot
          if (res > 1 || res < 0) {
            err_msg("boop must be <0,1>!");
            PrintError();
          } //}}}
          int k = res / width; // edge case for res == 1 (put into highest bin)
          if (k == bins) {
            k--; // edge case for res[j] -> 1
          }
          size_t index[3] = {type, k, j};
          AddArrND(boop_distr_type, index, res);
        }
      }

      for (int i = 0; i < Count->BeadCoor; i++) {
        for (int j = 0; j < Count->BeadCoor; j++) {
          if (i == j) {
            continue;
          }
          int id_i = System.BeadCoor[i];
          int id_j = System.BeadCoor[j];
          vec3 *pos_i = &System.Bead[id_i].Position;
          vec3 *pos_j = &System.Bead[id_j].Position;
          vec3 dist = Distance(pos_i->v, pos_j->v, System.Box.Length);
          double r_ij = VectLength(dist);
          if (r_ij >= r_max) {
            continue;
          }
          int bin = r_ij / dr;

          // boop are complex numbers
          // boop_i.real * boop_j.real + boop_i.imag * boop_j.imag
          for (int k = 0; k < n_sym; k++) {
            size_t id_i0[3] = {id_i, k, 0};
            size_t id_i1[3] = {id_i, k, 1};
            size_t id_j0[3] = {id_j, k, 0};
            size_t id_j1[3] = {id_j, k, 1};
            double corr_real = GetArrND(boop, id_i0) * GetArrND(boop, id_j0) +
                               GetArrND(boop, id_i1) * GetArrND(boop, id_j1);
            AddArr2D(g_n, k, bin, corr_real);
          }
          g_n_counts[bin]++;
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
  ArrNDd *boop_distr = CreateArr2Dd(bins, n_sym);
  // normalisation factor for distribution
  double *norm = calloc(n_sym, sizeof *norm);
  for (int i = 0; i < bins; i++) {
    for (int j = 0; j < n_sym; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        size_t index[] = {k, i, j};
        AddArr2D(boop_distr, i, j, GetArrND(boop_distr_type, index));
      }
      norm[j] += GetArr2D(boop_distr, i, j);
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
      data[i][++count] = GetArr2D(boop_distr, i, j) / norm[j];
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

  // print g_n to file
  // create file name TODO: needs something better, I guess
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
  // determine width of each column //{{{
  int columns2 = n_sym + 1;
  int digits2[columns2][2];
  InitInt2DArray((int *)digits2, columns2, 2, 0);
  double *data2[bins_g_n]; // array for data
  for (int i = 0; i < bins_g_n; i++) {
    data2[i] = calloc(columns2, sizeof data2[i]);
    count = -1;
    data2[i][++count] = dr * (2 * i + 1) / 2;
    for (int j = 0; j < n_sym; j++) {
      if (g_n_counts[i] > 0) {
        data2[i][++count] = GetArr2D(g_n, j, i) / g_n_counts[i];
      } else {
        data2[i][++count] = 0;
      }
    }
  }
  FillMaxDigits(columns2, bins_g_n, data2, digits2); //}}}
  for (int i = 0; i < bins_g_n; i++) {
    for (int col = 0; col < columns2; col++) {
      Fprintf1(fw, data2[i][col], digits2[col]);
    }
    putc('\n', fw);
    free(data2[i]);
  }
  fclose(fw);

  // free memory - to make valgrind happy //{{{
  FreeArrND(boop_distr_type);
  FreeArrND(g_n);
  FreeArrND(boop);
  FreeArrND(boop_distr);
  free(g_n_counts);
  free(norm);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
