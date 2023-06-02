#include "../AnalysisTools.h"

// TODO: add different units in file_extra 'potential' section?
// TODO: warning if missing parameters
// TODO: cutoff & CG into potential section of -fx <file>

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "TEXT TO BE ADDED\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output_gl> <output> [options]\n\n", cmd);

  fprintf(ptr, "   <input>        input coordinate file\n");
  fprintf(ptr, "   <width>        width of a single bin\n");
  fprintf(ptr, "   <output_gl>    output file with global observables\n");
  fprintf(ptr, "   <output>       output files with local variables "
               "(automatic endings -x.txt, -y.txt, and -z.txt\n");
  fprintf(ptr, "   [options]\n");
  // fprintf(ptr, "      --per-step  save a new file for each step"
  //         "(adds '-<step>.txt' to <output>)\n");
  fprintf(ptr, "      -fx <file>  file with extra information\n");
  fprintf(ptr, "      --pot       also calculate potential (and pressure via "
               "virial) ;only shifted Sutton-Chen for now\n");
  fprintf(ptr, "      --com       use centre of mass as the coordinate system "
               "centre; default: geometric centre");
  CommonHelp(error, n, opt);
} //}}}

typedef struct FF {
  double n, m, a, c, eps, rcut, rep1, rep2, dd1, dd2, CG;
} FF;

const static int max_params = 100;

// read potential parameters from the extra file //{{{
// shifted Sutton-Chen potential for one bead type //{{{
void ShiftedSCReadParameter();
void ShiftedSCParameters(FILE *fr, char file[], int *line_count, FF **ff,
                         SYSTEM System) {
  int el_params = 4, // number of parameters for the potential per-element
      gl_params = 2; // number of global parameters for the potential
  // global parameters
  char gl[gl_params][LINE];
  int count = 0;
  strcpy(gl[count++], "CG");
  strcpy(gl[count++], "cutoff");
  bool gl_done[gl_params];
  for (int i = 0; i < gl_params; i++) {
    gl_done[i] = false;
  }
  // per-element parameters
  char el[el_params][LINE];
  count = 0;
  strcpy(el[count++], "mn");
  strcpy(el[count++], "a");
  strcpy(el[count++], "c");
  strcpy(el[count++], "eps");
  bool el_done[el_params];
  for (int i = 0; i < el_params; i++) {
    el_done[i] = false;
  }

  char line[LINE], *split[SPL_STR];
  int words = 0;
  // read global parameters //{{{
  fpos_t position;
  fgetpos(fr, &position);
  while (ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
    // break the loop? //{{{
    if (gl_done[0] && gl_done[1] && // CG & cutoff keywords read
        words != 0 && // not a blank line
        strcasecmp(split[0], gl[0]) != 0 && // not extra CG line
        strcasecmp(split[0], gl[1]) != 0 && // not extra cutoff line
        split[0][0] != '#') { // not a comment
      break;
    } //}}}
    fgetpos(fr, &position); // may be a line with bead name
    (*line_count)++;
    double val;
    if (words > 1 && split[0][0] != '#') {
      // the second word must be a number
      if (!IsRealNumber(split[1], &val)) {
        goto error;
      } else if (strcasecmp(split[0], gl[0]) == 0) {
        // warn - multiple 'CG' lines
        if (gl_done[0]) {
          if (snprintf(ERROR_MSG, LINE, "%s%s%s parameter in 'ShiftedSC' is "
                       "specified multiple times; using this line",
                       ErrYellow(), gl[0], ErrCyan())) {
            strcpy(ERROR_MSG, "something wrong with snprintf()");
            PrintErrorFile(file, "\0", "\0");
            exit(1);
          }
          PrintWarningFileLine(file, *line_count, split, words);
        }
        ff[0][0].CG = val;
        gl_done[0] = true;
      } else if (strcasecmp(split[0], gl[1]) == 0) {
        // warn - multiple 'cutoff' lines
        if (gl_done[1]) {
          if (snprintf(ERROR_MSG, LINE, "%s%s%s parameter in 'ShiftedSC' is "
                       "specified multiple times; using this line",
                       ErrYellow(), gl[1], ErrCyan())) {
            strcpy(ERROR_MSG, "something wrong with snprintf()");
            PrintErrorFile(file, "\0", "\0");
            exit(1);
          }
          PrintWarningFileLine(file, *line_count, split, words);
        }
        ff[0][0].rcut = val;
        gl_done[1] = true;
      }
    }
  }
  // restore file pointer to the first bead type name
  fsetpos(fr, &position); //}}}
  // read parameters for all bead types //{{{
  while (strcasecmp(split[0], "finish") != 0) {
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'potential' section "
                     "(or missing 'finish' keyword)");
      exit(1);
    }
    (*line_count)++;
    int btype = FindBeadType(split[0], System);
    // skip section with an unknown bead type //{{{
    if (btype == -1) {
      snprintf(ERROR_MSG, LINE, "skipping unknown bead type %s%s%s "
               "in 'potential' section", ErrYellow(), split[0], ErrCyan());
      PrintWarningFileLine(file, *line_count, split, words);
      for (int i = 0; i < el_params; i++) {
        if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR,
                              " \t\n")) {
          ErrorEOF(file, "incomplete 'potential' section "
                         "(or missing 'finish' keyword)");
          exit(1);
        }
        (*line_count)++;
      }
      continue;
    } //}}}
    bool mn = false, a = false, c = false, eps = false;
    fgetpos(fr, &position);
    // read parameters for bead type 'btype' //{{{
    while (ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      // break the loop? //{{{
      if (mn && c && eps && a && // all keywords read
          words != 0 && // not a blank line
          strcasecmp(split[0], "mn") != 0 &&  // not extra keyword line
          strcasecmp(split[0], "c") != 0 &&   //
          strcasecmp(split[0], "eps") != 0 && //
          strcasecmp(split[0], "a") != 0 &&   //
          split[0][0] != '#') { // not a comment
        break;
      } //}}}
      fgetpos(fr, &position);
      (*line_count)++;
      if (words > 0) {
        if (strcasecmp(split[0], "mn") == 0) {
          long int val[2];
          if (words < 3 || !IsWholeNumber(split[1], &val[0]) ||
              !IsWholeNumber(split[2], &val[1])) {
            goto error;
          // m & n parameters
          } else if (mn) {
            snprintf(ERROR_MSG, LINE, "%sm%s and %sn%s parameters for %s%s%s "
                     "is specified multiple times; using this line",
                     ErrYellow(), ErrCyan(), ErrYellow(), ErrCyan(),
                     ErrYellow(), System.BeadType[btype].Name, ErrCyan());
            PrintWarningFileLine(file, *line_count, split, words);
          }
          if (val[0] > val[1]) {
            ff[btype][btype].m = val[0];
            ff[btype][btype].n = val[1];
          } else {
            ff[btype][btype].m = val[1];
            ff[btype][btype].n = val[0];
          }
          mn = true;
        // epsilon parameter
        } else if (strcasecmp(split[0], "eps") == 0) {
          double val;
          if (words < 2 || !IsRealNumber(split[1], &val)) {
            goto error;
          } else if (eps) {
            snprintf(ERROR_MSG, LINE, "%seps%s parameter for %s%s%s is "
                     "specified multiple times; using this line",
                     ErrYellow(), ErrCyan(), ErrYellow(),
                     System.BeadType[btype].Name, ErrCyan());
            PrintWarningFileLine(file, *line_count, split, words);
          }
          ff[btype][btype].eps = val;
          eps = true;
        // c parameter
        } else if (strcasecmp(split[0], "c") == 0) {
          double val;
          if (words < 2 || !IsRealNumber(split[1], &val)) {
            goto error;
          } else if (c) {
            snprintf(ERROR_MSG, LINE,
                     "%sc%s parameter for %s%s%s is specified "
                     "multiple times; using this line",
                     ErrYellow(), ErrCyan(), ErrYellow(),
                     System.BeadType[btype].Name, ErrCyan());
            PrintWarningFileLine(file, *line_count, split, words);
          }
          ff[btype][btype].c = val;
          c = true;
        // a parameter
        } else if (strcasecmp(split[0], "a") == 0) {
          double val;
          if (words < 2 || !IsRealNumber(split[1], &val)) {
            goto error;
          } else if (a) {
            snprintf(ERROR_MSG, LINE,
                     "%sa%s parameter for %s%s%s is specified "
                     "multiple times; using this line",
                     ErrYellow(), ErrCyan(), ErrYellow(),
                     System.BeadType[btype].Name, ErrCyan());
            PrintWarningFileLine(file, *line_count, split, words);
          }
          ff[btype][btype].a = val;
          a = true;
        // exit loop when 'finish' line encountered
        } else if (strcasecmp(split[0], "finish") == 0) {
          break;
        // ignore comment line
        } else if (split[0][0] != '#') {
          goto error;
        }
      }
    }
    fsetpos(fr, &position); //}}}
    // error if not all parameters specified //{{{
    if (!mn || !eps || !c || !a) {
      strcpy(ERROR_MSG, "premature end to 'potential' section");
      PrintErrorFileLine(file, *line_count, split, words);
      exit(1);
    } //}}}
  } //}}}
  return;
  error:
    strcpy(ERROR_MSG, "wrong line for 'ShiftedSC' potential");
    PrintErrorFileLine(file, *line_count, split, words);
    exit(1);
} //}}}
  //}}}

// calculate some prefactors for the potential //{{{
void SCPrefactors(SYSTEM System, double **par[max_params]) {
  // prefactors for shifted potential
  for (int i = 0; i < System.Count.BeadType; i++) {
    for (int j = 0; j < System.Count.BeadType; j++) {
      double m = par[0][i][j],
             n = par[1][i][j],
             a = par[2][i][j],
             eps = par[4][i][j],
             rc = par[5][i][j];
      // (a / r_c)^m * eps ...for repulsive part of SC
      par[6][i][j] = pow((a / rc), m) * eps;
      // m * eps * (a / r_c)^m / r_c ...for repulsive part of SC
      par[7][i][j] = par[6][i][j] * m / rc;
      // (a / r_c)^n ...for local density
      par[8][i][j] = pow((a / rc), n);
      // n * (a / r_c)^n / r_c ...for local density
      par[9][i][j] = par[8][i][j] * n / rc;
    }
  }
} //}}}

// calculate local density for Sutton-Chen potential //{{{
void LocalDensitySC(SYSTEM System, int bead1, int bead2,
                    double **par[max_params], double *density) {
  // rho = (a / r_ij)^n - (a / r_c)^n + (a / r_c)^n * n / r_c * (r_ij - r_c)
  BEAD *b_i = &System.Bead[bead1], *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type, bt_j = b_j->Type;
  // Sutton-Chen parameters for beads i and j
  double a = par[2][bt_i][bt_j],
         n = par[1][bt_i][bt_j],
         dd1 = par[8][bt_i][bt_j],
         dd2 = par[9][bt_i][bt_j],
         rc = par[5][bt_i][bt_j];
  if (rij < rc) {
    if (rij == 0) {
      strcpy(ERROR_MSG, "zero bead-bead distance!");
      PrintError();
      exit(1);
    }
    double tmp = pow(a / rij, n) - dd1 + dd2 * (rij - rc);
    density[bead1] += tmp;
    density[bead2] += tmp;
  }
} //}}}

// calculate potential, force and wirial for shifted Sutton-Chen potential //{{{
void ShiftedSuttonChen(SYSTEM System, int bead1, int bead2,
                       double **par[max_params], double *density,
                       double *potential, VECTOR *force,
                       double *Pxx, double *Pyy, double *Pzz, double *UCR) {
  BEAD *b_i = &System.Bead[bead1], *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type, bt_j = b_j->Type;
  // Sutton-Chen parameters
  double m = par[0][bt_i][bt_j],
         n = par[1][bt_i][bt_j],
         a = par[2][bt_i][bt_j],
         c_ii = par[3][bt_i][bt_i],
         c_jj = par[3][bt_j][bt_j],
         eps = par[4][bt_i][bt_j],
         eps_ii = par[4][bt_i][bt_i],
         eps_jj = par[4][bt_j][bt_j],
         rc = par[5][bt_i][bt_j],
         rep1 = par[6][bt_i][bt_j],
         rep2 = par[7][bt_i][bt_j],
         dd2 = par[9][bt_i][bt_j];
  if (rij < rc) {
    if (rij == 0) {
      strcpy(ERROR_MSG, "zero bead-bead distance!");
      PrintError();
      exit(1);
    }
    // helper value of (a / r_ij)^m & (a / r_ij)^n
    double a_div_r = a / rij;
    double a_div_r_m = pow(a_div_r, m),
           a_div_r_n = pow(a_div_r, n);
    // 1) repulsive energy
    // eps*(a/r_ij)^m-eps*(a/r_c)^m+eps*(a/r_c)^m*m/r_c*(r_ij-r_c)
    double tmp = eps * a_div_r_m - rep1 + rep2 * (rij - rc);
    // add to per-particle repulsive energy
    *UCR += tmp;
    potential[bead1] += tmp;
    potential[bead2] += tmp;
    // 2) repulsive virial
    // -m * eps * (a / r_ij)^m + m * eps * (a / r_c)^m / r_c * r_ij
    tmp = -m * eps * a_div_r_m + rep2 * rij;
    double w_rep = tmp;
    // 3) density-dependent virial
    // 0.5 * (c_i * eps_ii / sqrt(\rho_i) +
    //        c_j * eps_jj / sqrt(\rho_j) ...
    double c_eps_i = c_ii * eps_ii / sqrt(density[bead1]),
           c_eps_j = c_jj * eps_jj / sqrt(density[bead2]);
    tmp = 0.5 * (c_eps_i + c_eps_j);
    // ... * (n * (a / r_ij)^n - n * (a / r_c)^n / r_c * r_ij)
    tmp *= n * a_div_r_n - dd2 * rij;
    double w_dd = tmp;
    // 4) pairwise forces
    double fij = -(w_rep + w_dd) / SQR(rij);
    // force in axes' directions
    VECTOR f;
    f.x = fij * dist.x;
    f.y = fij * dist.y;
    f.z = fij * dist.z;
    // pressure tensor components
    *Pxx = f.x * dist.x;
    *Pyy = f.y * dist.y;
    *Pzz = f.z * dist.z;
    // per-particle force
    force[bead1].x += f.x;
    force[bead1].y += f.y;
    force[bead1].z += f.z;
    force[bead2].x -= f.x;
    force[bead2].y -= f.y;
    force[bead2].z -= f.z;
    // printf("%6d %6d", bead1, bead2);
    // printf(" %lf %lf %lf", dist.x, dist.y, dist.z);
    // printf(" %lf %lf %lf", w_dd, w_rep, fij);
    // printf(" %lf %lf %lf\n", f.x, f.y, f.z);
  }
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 12, all = common + 3, count = 0, req_arg = 2;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "--variable");
  strcpy(option[count++], "-pbc");
  strcpy(option[count++], "-ltrj");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-fx");
  strcpy(option[count++], "--pot");
  strcpy(option[count++], "--com");
  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type, struct_file,
                       &struct_type)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output_gl> & <output>
  char out_global[LINE] = "", out_local[LINE] = "";
  snprintf(out_global, LINE, "%s", argv[++count]);
  snprintf(out_local, LINE, "%s", argv[++count]);
  out_local[LINE - 6] = '\0'; // for adding -<axis>.rho

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);
  int trash;
  char file_extra[LINE] = "";
  if (FileIntegerOption(argc, argv, 0, "-fx", &trash, &trash, file_extra)) {
    exit(1);
  }
  // TODO: based on used potential; i.e., add to the potential reading function
  bool calculate_density = true; // for local density-dependent potentials
  // calculate pair-wise potential, etc.
  bool calculate_pairwise = BoolOption(argc, argv, "--pot");
  // use centre of mass instead of geometric centre
  bool use_com = BoolOption(argc, argv, "--com");
  //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz);

  COUNT *Count = &System.Count;
  BOX *box = &System.Box;

  // number of bins //{{{
  // error - box size must be specified //{{{
  if (box->Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintErrorFile(coor_file, struct_file, "\0"); // TODO add file(s)
    exit(1);
  } //}}}
  int tmp_bin[3];
  // TODO: *3 to assume box change of at most thrice as big
  //       remake to realloc on the fly
  tmp_bin[0] = ceil(box->Length.x / width);
  tmp_bin[1] = ceil(box->Length.y / width);
  tmp_bin[2] = ceil(box->Length.z / width);
  int n_bins = Max3(tmp_bin[0], tmp_bin[1], tmp_bin[2]) * 2;
  if ((n_bins % 2) == 1) {
    n_bins++;
  }
  int n_bins_orig = n_bins; //}}}

  FF **sc = malloc(sizeof **sc * Count->BeadType);
  for (int i = 0; i < Count->BeadType; i++) {
    sc[i] = calloc(Count->BeadType, sizeof *sc[i]);
  }

  // TODO: generalize the potential, etc.
  int CG = 1;
  char potential[LINE] = ""; // potential type
  double real_m = 1,         // mass
      real_E = 1,            // energy
      real_l = 1,            // length
      real_T = 1,            // temperature
      real_P = 1,            // pressure
      rcut = 1; // cut-off in real units (will be reduced by CG & real_l)
  // read extra information (if present) //{{{
  if (file_extra[0] != '\0') {
    FILE *fr = OpenFile(file_extra, "r");
    char line[LINE], *split[SPL_STR];
    int words, line_count = 0;
    while (ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      line_count++;
      double val;
      if (words > 1 && split[0][0] != '#') {
        if (strcasecmp("potential", split[0]) == 0) {
          // only read the section if potential should be calculated
          if (calculate_pairwise) {
            strcpy(potential, split[1]);
            if (strcasecmp(potential, "ShiftedSC") == 0) {
              ShiftedSCParameters(fr, file_extra, &line_count, sc, System);
              CG = sc[0][0].CG; // TODO: make it better?
              rcut = sc[0][0].rcut; // TODO: make it better?
            }
          } else {
            do {
              if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR,
                                    " \t\n")) {
                ErrorEOF(file_extra, "skipping 'potential' section");
                exit(1);
              }
              line_count++;
            } while (words == 0 || strcasecmp(split[0], "finish") != 0);
          }
        } else if (!IsRealNumber(split[1], &val)) {
          goto warning;
        } else if (strncasecmp("Energy", split[0], 2) == 0) {
          real_E = val;
        } else if (strcasecmp("Temperature", split[0]) == 0) {
          real_T = val;
        } else if (strcasecmp("Length", split[0]) == 0) {
          real_l = val;
        } else if (strcasecmp("Mass", split[0]) == 0) {
          real_m = val;
        } else if (strcasecmp("Pressure", split[0]) == 0) {
          real_P = val;
        } else {
          warning:
            strcpy(ERROR_MSG, "unrecognized line");
            PrintWarningFileLine(file_extra, line_count, split, words);
        }
      }
    }
    fclose(fr);
  } //}}}
  // TODO: still only shifted SC
  // cross-terms for potentials //{{{
  for (int i = 0; i < (Count->BeadType - 1); i++) {
    for (int j = (i + 1); j < Count->BeadType; j++) {
      sc[i][j].n = (sc[i][i].n + sc[j][j].n) / 2;
      sc[i][j].m = (sc[i][i].m + sc[j][j].m) / 2;
      sc[i][j].eps = sqrt(sc[i][i].eps * sc[j][j].eps);
      sc[i][j].c = -1; // no cross-term in SC
      sc[i][j].a = sqrt(sc[i][i].a * sc[j][j].a);
      // sc[i][j].rcut = rcut;
      sc[j][i].n = sc[i][j].n;
      sc[j][i].m = sc[i][j].m;
      sc[j][i].eps = sc[i][j].eps;
      sc[j][i].c = sc[i][j].c;
      sc[j][i].a = sc[i][j].a;
      // sc[j][i].rcut = rcut;
    }
  } //}}}
  // coarse-graining of bead mass and potentials //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Mass *= CUBE(CG) / real_m;
    for (int j = 0; j < Count->BeadType; j++) {
      sc[i][j].a *= CG / real_l;
      sc[i][j].eps *= CUBE(CG) / real_T; // TODO: always T?
      sc[i][j].rcut = rcut * CG / real_l;
    }
  } //}}}
  // warn if no potential info even if potential is to be calculated //{{{
  if (calculate_pairwise && potential[0] == '\0') {
    strcpy(ERROR_MSG, "missing information about potential");
    PrintWarnFile(file_extra, "\0", "\0");
  } //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  // warning - unspecified mass //{{{
  for (int i = 0; i < System.Count.BeadType; i++) {
    if (System.BeadType[i].Mass == MASS) {
      snprintf(ERROR_MSG, LINE, "unspecified mass (bead type %s%s%s)",
               ErrYellow(), System.BeadType[i].Name, ErrCyan());
      PrintWarnFile(struct_file, coor_file, "\0");
    }
  } //}}}

  // write header to output file //{{{
  FILE *fw = OpenFile(out_global, "w");
  PrintByline(fw, argc, argv);
  count = 0;
  fprintf(fw, "# (%d) timestep;", ++count);
  fprintf(fw, " (%d) T;", ++count);
  fprintf(fw, " (%d) KE;", ++count);
  if (calculate_pairwise) {
    fprintf(fw, " (%d) P(KE);", ++count);
    fprintf(fw, " (%d) P(W);", ++count);
    fprintf(fw, " (%d) P(W+KE);", ++count);
    fprintf(fw, " (%d) Pxx;", ++count);
    fprintf(fw, " (%d) Pyy;", ++count);
    fprintf(fw, " (%d) Pzz;", ++count);
    fprintf(fw, " (%d) W;", ++count);
    fprintf(fw, " (%d) uCR;", ++count);
    fprintf(fw, " (%d) uCD;", ++count);
  }
  putc('\n', fw);
  fclose(fw); //}}}

  // variables for collecting observables //{{{
  bool *bt_in_use = calloc(Count->BeadType, sizeof *bt_in_use);
  long double(*KEnergy)[3] = calloc(n_bins, sizeof *KEnergy),
             (*Press_vir)[3] = calloc(n_bins, sizeof *Press_vir);
  long int(**BeadCount)[3] = malloc(Count->BeadType * sizeof **BeadCount);
  for (int i = 0; i < Count->BeadType; i++) {
    BeadCount[i] = calloc(n_bins, sizeof *BeadCount[i]);
  }
  double temp_sum = 0, press_kin_sum = 0, virial_sum = 0, press_vir_sum = 0,
         centre_avg[3] = {0, 0, 0};
  //}}}

  double **ff_par[max_params];
  for (int i = 0; i < max_params; i++) {
    ff_par[i] = calloc(Count->BeadType, sizeof *ff_par);
    for (int j = 0; j < Count->BeadType; j++) {
      ff_par[i][j] = calloc(Count->BeadType, sizeof *ff_par);
    }
  }
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = i; j < Count->BeadType; j++) {
      count = 0;
      ff_par[ 0][i][j] = sc[i][j].m;
      ff_par[ 1][i][j] = sc[i][j].n;
      ff_par[ 2][i][j] = sc[i][j].a;
      ff_par[ 3][i][j] = sc[i][j].c;
      ff_par[ 4][i][j] = sc[i][j].eps;
      ff_par[ 5][i][j] = sc[i][j].rcut;
      ff_par[ 6][i][j] = sc[i][j].rep1;
      ff_par[ 7][i][j] = sc[i][j].rep2;
      ff_par[ 8][i][j] = sc[i][j].dd1;
      ff_par[ 9][i][j] = sc[i][j].dd2;
      ff_par[10][i][j] = sc[i][j].CG;
      for (int k = 0; k < 11; k++) {
        ff_par[k][j][i] = ff_par[k][i][j];
      }
    }
  }

  SCPrefactors(System, ff_par);

  // main loop //{{{
  FILE *coor = OpenFile(coor_file, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file

  while (true) {
    PrintStep(&count_coor, start, silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if ((count_coor >= start && (count_coor <= end || end == -1)) &&
        ((count_coor - start) % skip) == 0) {
      use = true;
    }
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (!ReadTimestep(coor_type, coor, coor_file, &System, &line_count,
                        vtf_var)) {
        count_coor--;
        break;
      }
      // some variables (may not be used); TODO: clear names
      double *loc_dens = calloc(Count->Bead, sizeof *loc_dens);
      double *uSC = calloc(Count->Bead, sizeof *uSC), Pxx = 0, Pyy = 0, Pzz = 0,
             energy_kin = 0, UCR = 0, UCD = 0; // TODO: remove?
      VECTOR *force = calloc(Count->Bead, sizeof *force);
      bool *test_used = calloc(Count->Bead, sizeof *test_used);

      count_used++;
      box = &System.Box;
      // TODO: this is for ml's stuff - coordinate centre is in the box centre
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        b->Position.x += box->Length.x / 2;
        b->Position.y += box->Length.y / 2;
        b->Position.z += box->Length.z / 2;
      }
      WrapJoinCoordinates(&System, true, false); // restore pbc
      VECTOR centre;
      if (use_com) {
        centre = CentreOfMass(Count->BeadCoor, System.BeadCoor, System);
      } else {
        centre = GeomCentre(Count->BeadCoor, System.BeadCoor, System.Bead);
      }
      centre_avg[0] += centre.x;
      centre_avg[1] += centre.y;
      centre_avg[2] += centre.z;
      if (calculate_pairwise) { //{{{
        // calculete forces
        INTVECTOR n_cells;
        int *Head, *Link, Dcx[14], Dcy[14], Dcz[14];
        // access beads in the linked list through BeadCoor[]
        // TODO: rcut - pass maximum of sc[][].rcut, I guess
        LinkedList(System, &Head, &Link, sc[0][0].rcut, &n_cells, Dcx, Dcy,
                   Dcz);
        // calculate density (if necessary) //{{{
        if (calculate_density) {
          for (int c1z = 0; c1z < n_cells.z; c1z++) {
            for (int c1y = 0; c1y < n_cells.y; c1y++) {
              for (int c1x = 0; c1x < n_cells.x; c1x++) {
                // select first cell
                int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;
                // select first bead in the cell 'cell1'
                int i = Head[cell1];
                while (i != -1) {
                  for (int k = 0; k < 14; k++) {
                    int c2x = c1x + Dcx[k]; //{{{
                    int c2y = c1y + Dcy[k];
                    int c2z = c1z + Dcz[k];
                    // periodic boundary conditions for cells //{{{
                    if (c2x >= n_cells.x)
                      c2x -= n_cells.x;
                    else if (c2x < 0)
                      c2x += n_cells.x;

                    if (c2y >= n_cells.y)
                      c2y -= n_cells.y;
                    else if (c2y < 0)
                      c2y += n_cells.y;

                    if (c2z >= n_cells.z)
                      c2z -= n_cells.z; //}}}

                    // select second cell
                    int cell2 = c2x + c2y * n_cells.x +
                                c2z * n_cells.x * n_cells.y; //}}}
                    // select bead in the cell 'cell2' //{{{
                    int j;
                    if (cell1 == cell2) { // next bead in 'cell1'
                      j = Link[i];
                    } else { // first bead in 'cell2'
                      j = Head[cell2];
                    } //}}}
                    while (j != -1) {
                      int id_i = System.BeadCoor[i], id_j = System.BeadCoor[j];
                      LocalDensitySC(System, id_i, id_j, ff_par, loc_dens);
                      j = Link[j];
                    }
                  }
                  i = Link[i];
                }
              }
            }
          }
        } //}}}
        // calculate pair-wise forces, energy, virial //{{{
        for (int c1z = 0; c1z < n_cells.z; c1z++) {
          for (int c1y = 0; c1y < n_cells.y; c1y++) {
            for (int c1x = 0; c1x < n_cells.x; c1x++) {
              // select first cell
              int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;
              // select first bead in the cell 'cell1'
              int i = Head[cell1];
              while (i != -1) {
                int id_i = System.BeadCoor[i];
                for (int k = 0; k < 14; k++) {
                  int c2x = c1x + Dcx[k]; //{{{
                  int c2y = c1y + Dcy[k];
                  int c2z = c1z + Dcz[k];
                  // periodic boundary conditions for cells //{{{
                  if (c2x >= n_cells.x)
                    c2x -= n_cells.x;
                  else if (c2x < 0)
                    c2x += n_cells.x;

                  if (c2y >= n_cells.y)
                    c2y -= n_cells.y;
                  else if (c2y < 0)
                    c2y += n_cells.y;

                  if (c2z >= n_cells.z)
                    c2z -= n_cells.z; //}}}

                  // select second cell
                  int cell2 =
                      c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}
                  // select bead in the cell 'cell2' //{{{
                  int j;
                  if (cell1 == cell2) { // next bead in 'cell1'
                    j = Link[i];
                  } else { // first bead in 'cell2'
                    j = Head[cell2];
                  } //}}}
                  while (j != -1) {
                    int id_j = System.BeadCoor[j];
                    // bin ids in the three directions //{{{
                    BEAD *b_i = &System.Bead[id_i];
                    BEAD *b_j = &System.Bead[id_j];
                    int bin_i[3], bin_j[3];
                    // hopefully, bin will always be positive
                    bin_i[0] =
                        (b_i->Position.x - centre.x) / width + n_bins_orig / 2;
                    bin_i[1] =
                        (b_i->Position.y - centre.y) / width + n_bins_orig / 2;
                    bin_i[2] =
                        (b_i->Position.z - centre.z) / width + n_bins_orig / 2;
                    // hopefully, bin will always be positive
                    bin_j[0] =
                        (b_j->Position.x - centre.x) / width + n_bins_orig / 2;
                    bin_j[1] =
                        (b_j->Position.y - centre.y) / width + n_bins_orig / 2;
                    bin_j[2] = (b_j->Position.z - centre.z) / width +
                               n_bins_orig / 2; //}}}
                    double pxx_ij = 0, pyy_ij = 0, pzz_ij = 0;
                    ShiftedSuttonChen(System, id_i, id_j, ff_par, loc_dens, uSC,
                                      force, &pxx_ij, &pyy_ij, &pzz_ij, &UCR);
                    Pxx += pxx_ij;
                    Pyy += pyy_ij;
                    Pzz += pzz_ij;
                    for (int k = 0; k < 3; k++) {
                      if (bin_i[k] == bin_j[k]) {
                        double press = pxx_ij + pyy_ij + pzz_ij;
                        Press_vir[bin_i[k]][k] += press;
                      }
                    }
                    j = Link[j];
                  }
                }
                i = Link[i];
              }
            }
          }
        } //}}}
        free(Head);
        free(Link);
      } //}}}
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        // kinetic enrgy
        double vel2 =
            SQR(b->Velocity.x) + SQR(b->Velocity.y) + SQR(b->Velocity.z);
        double ke = 0.5 * bt->Mass * vel2;
        energy_kin += ke;
        if (calculate_pairwise) {
          // Sutton-Chen potential
          double ui = sc[b->Type][b->Type].c * sc[b->Type][b->Type].eps *
                      sqrt(loc_dens[id]);
          UCD -= ui;
          uSC[id] = 0.5 * uSC[id] - ui;
          // copy previously calculated forces
          b->Force = force[id];
        }
        // bin id in the three directions
        int bin[3];
        // hopefully, bin will always be positive
        bin[0] = (b->Position.x - centre.x) / width + n_bins_orig / 2;
        bin[1] = (b->Position.y - centre.y) / width + n_bins_orig / 2;
        bin[2] = (b->Position.z - centre.z) / width + n_bins_orig / 2;

        // realloc memory if bin is too high //{{{
        if (Max3(bin[0], bin[1], bin[2]) >= n_bins) {
          printf("\nchanged: %d", n_bins);
          n_bins = Max3(bin[0], bin[1], bin[2]) + 1;
          printf(" -> %d\n", n_bins);
          KEnergy = realloc(KEnergy, n_bins * sizeof *KEnergy);
          for (int j = 0; j < Count->BeadType; j++) {
            BeadCount[j] = realloc(BeadCount[j], n_bins * sizeof *BeadCount[j]);
          }
        } //}}}
        // error - negative bin //{{{
        if (bin[0] < 0 || bin[1] < 0 || bin[2] < 0) {
          strcpy(ERROR_MSG, "negative bin id; should never happen!");
          PrintError();
          exit(1);
        } //}}}
        for (int j = 0; j < 3; j++) {
          BeadCount[b->Type][bin[j]][j]++;
          KEnergy[bin[j]][j] += ke;
        }
        // specify the given bead type is used
        bt_in_use[b->Type] = true;
      }

      // global observables
      double temperature = 2 * energy_kin / (3 * Count->BeadCoor - 3);
      // kinetic pressure (i.e., ideal gas contribution)
      double press_kin = 2 * energy_kin / box->Volume;
      // virial from the axes' contributions
      double virial = (Pxx + Pyy + Pzz) / 3;
      // virial contribution to pressure
      double press_vir = virial / box->Volume;
      // write global values to output file //{{{
      fw = OpenFile(out_global, "a");
      fprintf(fw, " %8d", count_coor);
      fprintf(fw, " %10.5e", temperature * real_T / CUBE(CG));
      fprintf(fw, " %10.5e", energy_kin * real_E / CUBE(CG));
      if (calculate_pairwise) {
        fprintf(fw, " %10.5e", press_kin * real_P);
        fprintf(fw, " %10.5e", press_vir * real_P);
        fprintf(fw, " %10.5e", (press_kin + press_vir) * real_P);
        fprintf(fw, " %10.5e", Pxx);
        fprintf(fw, " %10.5e", Pyy);
        fprintf(fw, " %10.5e", Pzz);
        fprintf(fw, " %10.5e", virial);
        fprintf(fw, " %10.5e", UCR / Count->BeadCoor * real_E / CUBE(CG));
        fprintf(fw, " %10.5e", UCD / Count->BeadCoor * real_E / CUBE(CG));
      }
      putc('\n', fw);
      fclose(fw); //}}}

      // add global observables to total sums
      temp_sum += temperature;
      press_kin_sum += press_kin;
      virial_sum += virial;
      press_vir_sum += press_vir;

      free(loc_dens);
      free(uSC);
      free(force);
      free(test_used);
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(coor_type, coor, coor_file, struct_file, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    if (count_coor == end) {
      break;
    }
  }
  fclose(coor);          //}}}
  if (count_coor == 0) { // error - input file without a valid timestep //{{{
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(coor_file, "\0", "\0");
    fputc('\n', stderr);           //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG,
           "no timestep used; "
           "starting timestep is higher than the number of timesteps");
    PrintWarning();     //}}}
  } else if (!silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
    fflush(stdout);
  } //}}}

  // average centre of mass/geometry
  centre_avg[0] /= count_used;
  centre_avg[1] /= count_used;
  centre_avg[2] /= count_used;

  // print profiles of observables //{{{
  // open the three files & print headers //{{{
  char f[3][LINE]; // complete filenames
  if (snprintf(f[0], LINE, "%s-x.rho", out_local) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    PrintError();
    exit(1);
  }
  if (snprintf(f[1], LINE, "%s-y.rho", out_local) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    PrintError();
    exit(1);
  }
  if (snprintf(f[2], LINE, "%s-z.rho", out_local) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    PrintError();
    exit(1);
  }
  FILE *fw2[3];
  for (int i = 0; i < 3; i++) {
    fw2[i] = OpenFile(f[i], "w");
    // print headers
    PrintByline(fw2[i], argc, argv);
    fprintf(fw2[i], "# columns: (1) distance");
    count = 1;
    fprintf(fw2[i], "; (%d) T", ++count);
    if (calculate_pairwise) {
      fprintf(fw2[i], "; (%d) P(KE)", ++count);
      fprintf(fw2[i], "; (%d) P(W)", ++count);
      fprintf(fw2[i], "; (%d) P", ++count);
    }
    for (int j = 0; j < Count->BeadType; j++) {
      if (bt_in_use[j]) {
        count++;
        fprintf(fw2[i], "; (%d) %s", count, System.BeadType[j].Name);
      }
    }
    putc('\n', fw2[i]);
  } //}}}
  // find first/last non-zero bin in every direction //{{{
  int non_zero_lo[3] = {-1, -1, -1};
  for (int i = 0; i < n_bins; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < 3; k++) {
        if (non_zero_lo[k] == -1 && BeadCount[j][i][k] > 0) {
          non_zero_lo[k] = i;
        }
      }
    }
  }
  int non_zero_hi[3] = {-1, -1, -1};
  for (int i = (n_bins - 1); i >= 0; i--) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < 3; k++) {
        if (non_zero_hi[k] == -1 && BeadCount[j][i][k] > 0) {
          non_zero_hi[k] = i;
        }
      }
    }
  } //}}}
  // helper variables //{{{
  long int(*BeadCount_tot)[3] = calloc(n_bins, sizeof *BeadCount_tot);
  for (int i = 0; i < n_bins; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < 3; k++) {
        BeadCount_tot[i][k] += BeadCount[j][i][k];
      }
    }
  }
  double volume[3];
  volume[0] = width * box->Length.y * box->Length.z;
  volume[1] = box->Length.x * width * box->Length.z;
  volume[2] = box->Length.x * box->Length.y * width;
  //}}}
  // write profiles //{{{
  for (int i = 0; i < (n_bins - 1); i++) {
    double dist = width * (2 * i + 1) / 2;
    for (int j = 0; j < 3; j++) {
      if (i >= non_zero_lo[j] && i <= non_zero_hi[j]) {
        double abs_dist =
            dist - centre_avg[j] - width * (2 * non_zero_lo[j] + 1) / 2;
        fprintf(fw2[j], "%7.3f", abs_dist * real_l); // absolute distance
        // temperature - don't use 'count_used' as both KEnergy & BeadCount_tot
        // are sums from all the counted steps
        double temp = 2 * KEnergy[i][j] / (3 * BeadCount_tot[i][j] - 3);
        fprintf(fw2[j], " %e", temp * real_T / CUBE(CG));
        if (calculate_pairwise) {
          // kinetic (ideal gas) pressure
          double press_kin = 2 * KEnergy[i][j] / (volume[j] * count_used);
          fprintf(fw2[j], " %e", press_kin * real_P);
          // virial pressure
          double press_vir = Press_vir[i][j] / (3 * volume[j] * count_used);
          fprintf(fw2[j], " %e", press_vir * real_P);
          // total pressure
          fprintf(fw2[j], " %e", (press_kin + press_vir) * real_P);
        }
        // density
        for (int k = 0; k < Count->BeadType; k++) {
          if (bt_in_use[k]) {
            double rho = BeadCount[k][i][j] / (volume[j] * count_used);
            rho *= CUBE(CG) / CUBE(real_l);
            fprintf(fw2[j], " %e", rho);
          }
        }
        putc('\n', fw2[j]);
      }
    }
  } //}}}
  fclose(fw2[0]);
  fclose(fw2[1]);
  fclose(fw2[2]); //}}}

  // free memory //{{{
  FreeSystem(&System);
  for (int i = 0; i < Count->BeadType; i++) {
    free(sc[i]);
  }
  free(sc);
  for (int i = 0; i < Count->BeadType; i++) {
    free(BeadCount[i]);
  }
  free(BeadCount);
  free(KEnergy);
  free(Press_vir);
  free(BeadCount_tot);
  free(bt_in_use);
  //}}}

  return 0;
}
