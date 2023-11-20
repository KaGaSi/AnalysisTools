#include "../AnalysisTools.h"

// TODO: add different units in file_extra 'potential' section?
// TODO: warning if missing parameters
// TODO: calculate_density - base it on used potential
// TODO: potentials into number codes so I don't need to use strcmp all the time

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "TEXT TO BE ADDED\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <width> <out_gl> <out_lc> "
          "[options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single bin\n");
  fprintf(ptr, "<out_gl>            output file with global observables\n");
  fprintf(ptr, "<output>            output files with local variables "
          "(three files with automatic endings -x.txt, -y.txt, and -z.txt\n");
  fprintf(ptr, "[options]\n");
  // fprintf(ptr, "      --per-step  save a new file for each step"
  //         "(adds '-<step>.txt' to <output>)\n");
  fprintf(ptr, "  -fx <file>        file with extra information\n");
  fprintf(ptr, "  --pot             also calculate potential (and pressure via "
          "virial) ;only shifted Sutton-Chen for now\n");
  fprintf(ptr, "  --com             use centre of mass as the coordinate system "
          "centre; default: simulation box centre\n");
  fprintf(ptr, "  --cog             use geometic centre as the coordinate system "
          "centre; default: simulation box centre\n");
  CommonHelp(error, n, opt);
} //}}}

const static int max_params = 100;

// read potential parameters from the extra file //{{{
// (shifted) Sutton-Chen potential for one bead type //{{{
void SuttonChenParameters(FILE *fr, char file[], int *line_count,
                          double **par[max_params], SYSTEM System) {
  // variables //{{{
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
  } //}}}
  // read global parameters //{{{
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    (*line_count)++;
    // break the loop when 'bead' or 'finish' keyword encountered //{{{
    if (words != 0 && (strcasecmp(split[0], "bead") == 0 ||
                       strcasecmp(split[0], "finish") == 0)) {
      break;
    } //}}}
    double val = 0;
    if (words > 1 && split[0][0] != '#') {
      // the second word must be a number
      if (!IsRealNumber(split[1], &val)) {
        goto error;
      } else {
        int i = 0;
        for (; i < gl_params; i++) {
          if (strcasecmp(split[0], gl[i]) == 0) {
            break;
          }
        }
        // unrecognized line
        if (i == gl_params) {
          goto error;
        }
        // warning - multiply defined parameter
        if (gl_done[i]) {
          if (snprintf(ERROR_MSG, LINE, "%s%s%s Sutton-Chen parameter is "
                       "specified multiple times; using this line",
                       ErrYellow(), gl[i], ErrCyan()) < 0) {
            ErrorSnprintf();
          }
          PrintWarnFileLine(file, *line_count);
        }
        gl_done[i] = true;
        switch (i) {
        case 0:
          par[10][0][0] = val; // CG
          break;
        case 1:
          par[5][0][0] = val; // rcut
          break;
        }
      }
    }
  } //}}}
  // read parameters for every bead types //{{{
  do { // first 'bead <name>' was already read
    int btype = FindBeadType(split[1], System);
    // skip section with an unknown bead type //{{{
    if (btype == -1) {
      if (snprintf(ERROR_MSG, LINE, "skipping unknown bead type "
                   "%s%s%s in the 'potential' section", ErrYellow(),
                   split[1], ErrCyan()) < 0) {
        ErrorSnprintf();
      }
      PrintWarnFileLine(file, *line_count);
      while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        (*line_count)++;
        if (words > 0 && (strcasecmp(split[0], "bead") == 0 ||
                          strcasecmp(split[0], "finish") == 0)) {
          break;
        }
      } //}}}
    // read parameters for bead type 'btype' //{{{
    } else {
      for (int i = 0; i < el_params; i++) {
        el_done[i] = false;
      }
      while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        (*line_count)++;
        // break loop when 'bead' or 'finish' keyword encountered
        if (words > 0 && (strcasecmp(split[0], "bead") == 0 ||
                          strcasecmp(split[0], "finish") == 0)) {
          break;
        }
        // ignore blanka and comment lines
        if (words == 0 || split[0][0] == '#') {
          continue;
        }
        double val = 0;
        // the second word must be a number
        if (!IsRealNumber(split[1], &val)) {
          goto error;
        } else {
          // identify keyword
          int i = 0;
          for (; i < el_params; i++) {
            if (strcasecmp(split[0], el[i]) == 0) {
              break;
            }
          }
          // unrecognized keyword
          if (i == el_params) {
            goto error;
          }
          // multiply defined keyword
          if (el_done[i]) {
            if (snprintf(ERROR_MSG, LINE, "%s%s%s parameter in 'ShiftedSC' "
                         "for %s%s%s is specified multiple times; "
                         "using this line", ErrYellow(), el[i],
                         ErrCyan(), ErrYellow(),
                         System.BeadType[btype].Name, ErrCyan()) < 0) {
              ErrorSnprintf();
            }
            PrintWarnFileLine(file, *line_count);
          }
          el_done[i] = true; // keyword 'i' found
          // assign value(s)
          switch (i) {
            case 0:;
              double val2;
              if (words < 3 || !IsPosRealNumber(split[2], &val2) || val <= 0) {
                goto error;
              }
              if (val > val2) {
                par[0][btype][btype] = val;  // m
                par[1][btype][btype] = val2; // n
              } else {
                par[0][btype][btype] = val2; // m
                par[1][btype][btype] = val;  // n
              }
              break;
            case 1:
              par[2][btype][btype] = val; // a
              break;
            case 2:
              par[3][btype][btype] = val; // c
              break;
            case 3:
              par[4][btype][btype] = val; // epsilon
              break;
          }
        }
      }
      // error if not all parameters specified
      if (!el_done[0] || !el_done[1] || !el_done[2] || !el_done[3]) {
        strcpy(ERROR_MSG, "premature end to 'potential' section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    }
  } while (strcasecmp(split[0], "finish") != 0); //}}}
  //}}}
  for (int i = 1; i < System.Count.BeadType; i++) {
    par[ 5][i][i] = par[ 5][0][0];
    par[10][i][i] = par[10][0][0];
  }
  return;
  error:
    strcpy(ERROR_MSG, "wrong line for Sutton-Chen potential");
    PrintErrorFileLine(file, *line_count);
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
      // (a / r_c)^m * eps ...for repulsive part of ShiftedSC
      par[6][i][j] = pow((a / rc), m) * eps; // not used, I think
      // m * eps * (a / r_c)^m / r_c ...for repulsive part of ShiftedSC
      par[7][i][j] = par[6][i][j] * m / rc;
      // (a / r_c)^n ...for local density
      par[8][i][j] = pow((a / rc), n);
      // n * (a / r_c)^n / r_c ...for local density
      par[9][i][j] = par[8][i][j] * n / rc;
    }
  }
} //}}}

// calculate local density for (shifted) Sutton-Chen potential //{{{
// shifted Sutton-Chen //{{{
void LocalDensityShiftedSC(SYSTEM System, int bead1, int bead2,
                    double **par[max_params], double *density) {
  // rho = (a / r_ij)^n - (a / r_c)^n + (a / r_c)^n * n / r_c * (r_ij - r_c)
  BEAD *b_i = &System.Bead[bead1], *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type, bt_j = b_j->Type;
  // Sutton-Chen parameters for beads i and j
  double a = par[2][bt_i][bt_j], n = par[1][bt_i][bt_j],
         dd1 = par[8][bt_i][bt_j], dd2 = par[9][bt_i][bt_j],
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
// Sutton-Chen //{{{
void LocalDensitySC(SYSTEM System, int bead1, int bead2,
                    double **par[max_params], double *density) {
  // rho = (a / r_ij)^n - (a / r_c)^n + (a / r_c)^n * n / r_c * (r_ij - r_c)
  BEAD *b_i = &System.Bead[bead1], *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type, bt_j = b_j->Type;
  // Sutton-Chen parameters for beads i and j
  double a = par[2][bt_i][bt_j], n = par[1][bt_i][bt_j],
         rc = par[5][bt_i][bt_j];
  if (rij < rc) {
    if (rij == 0) {
      strcpy(ERROR_MSG, "zero bead-bead distance!");
      PrintError();
      exit(1);
    }
    double tmp = pow(a / rij, n);
    density[bead1] += tmp;
    density[bead2] += tmp;
  }
} //}}}
 //}}}

// calculate potential, force and wirial for shifted Sutton-Chen potential //{{{
// shifted Sutton-Chen //{{{
void ShiftedSuttonChen(SYSTEM System, int bead1, int bead2,
                       double **par[max_params], double *density,
                       VECTOR *P_diag) {
  BEAD *b_i = &System.Bead[bead1], *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type, bt_j = b_j->Type;
  // Sutton-Chen parameters
  double m = par[0][bt_i][bt_j], n = par[1][bt_i][bt_j],
         a = par[2][bt_i][bt_j],
         c_ii = par[3][bt_i][bt_i], c_jj = par[3][bt_j][bt_j],
         eps_ij = par[4][bt_i][bt_j],
         eps_ii = par[4][bt_i][bt_i], eps_jj = par[4][bt_j][bt_j],
         rc = par[5][bt_i][bt_j],
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
    double a_div_r_m = pow(a_div_r, m), a_div_r_n = pow(a_div_r, n);
    // repulsive virial
    // -m * eps * (a / r_ij)^m + m * eps * (a / r_c)^m / r_c * r_ij
    double w_rep = -m * eps_ij * a_div_r_m + rep2 * rij;
    // density-dependent virial
    // 0.5 * (c_i * eps_ii / sqrt(\rho_i) + c_j * eps_jj / sqrt(\rho_j)) *
    // (n * (a / r_ij)^n - n * (a / r_c)^n / r_c * r_ij)
    double c_eps_i = c_ii * eps_ii / sqrt(density[bead1]),
           c_eps_j = c_jj * eps_jj / sqrt(density[bead2]);
    double w_dd = 0.5 * (c_eps_i + c_eps_j) * (n * a_div_r_n - dd2 * rij);
    // pairwise forces
    double fij = -(w_rep + w_dd) / SQR(rij);
    System.Bead[bead1].Force.x += fij * dist.x;
    System.Bead[bead1].Force.y += fij * dist.y;
    System.Bead[bead1].Force.z += fij * dist.z;
    System.Bead[bead2].Force.x -= fij * dist.x;
    System.Bead[bead2].Force.y -= fij * dist.y;
    System.Bead[bead2].Force.z -= fij * dist.z;
    // pressure tensor components
    P_diag->x = fij * SQR(dist.x);
    P_diag->y = fij * SQR(dist.y);
    P_diag->z = fij * SQR(dist.z);
  }
} //}}}
// Sutton-Chen //{{{
// TODO: doesn't work - at least against my lammps stuff
void SuttonChen(SYSTEM System, int bead1, int bead2, double **par[max_params],
                double *density, VECTOR *P_diag) {
  BEAD *b_i = &System.Bead[bead1], *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type, bt_j = b_j->Type;
  // Sutton-Chen parameters
  double m = par[0][bt_i][bt_j], n = par[1][bt_i][bt_j],
         a = par[2][bt_i][bt_j],
         c_ii = par[3][bt_i][bt_i], c_jj = par[3][bt_j][bt_j],
         eps_ij = par[4][bt_i][bt_j],
         eps_ii = par[4][bt_i][bt_i], eps_jj = par[4][bt_j][bt_j],
         rc = par[5][bt_i][bt_j];
  if (rij < rc) {
    if (rij == 0) {
      strcpy(ERROR_MSG, "zero bead-bead distance!");
      PrintError();
      exit(1);
    }
    // helper value of (a / r_ij)^m & (a / r_ij)^n
    double a_div_r = a / rij;
    double a_div_r_m = pow(a_div_r, m), a_div_r_n = pow(a_div_r, n);
    // repulsive virial
    // -m * eps * (a / r_ij)^m
    double w_rep = -m * eps_ij * a_div_r_m;
    // density-dependent virial
    // 0.5 * n * (c_i * eps_ii / sqrt(\rho_i) +
    //            c_j * eps_jj / sqrt(\rho_j)) * (a / r_ij)^n
    double c_eps_i = c_ii * eps_ii / sqrt(density[bead1]),
           c_eps_j = c_jj * eps_jj / sqrt(density[bead2]);
    double w_dd = 0.5 * n * (c_eps_i + c_eps_j) * a_div_r_n;
    // pairwise forces
    double fij = -(w_rep + w_dd) / SQR(rij);
    System.Bead[bead1].Force.x += fij * dist.x;
    System.Bead[bead1].Force.y += fij * dist.y;
    System.Bead[bead1].Force.z += fij * dist.z;
    System.Bead[bead2].Force.x -= fij * dist.x;
    System.Bead[bead2].Force.y -= fij * dist.y;
    System.Bead[bead2].Force.z -= fij * dist.z;
    // pressure tensor components
    P_diag->x = fij * SQR(dist.x);
    P_diag->y = fij * SQR(dist.y);
    P_diag->z = fij * SQR(dist.z);
  }
} //}}}
  //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 11, all = common + 4, count = 0, req_arg = 2;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  // strcpy(option[count++], "--variable"); // TODO: makes no sense, I think
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
  strcpy(option[count++], "--cog");
  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorOption("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output_gl> & <output>
  char out_global[LINE] = "", out_local[LINE] = "";
  snprintf(out_global, LINE, "%s", argv[++count]);
  snprintf(out_local, LINE, "%s", argv[++count]);
  out_local[LINE-6] = '\0'; // for adding -<axis>.rho

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &start, &end, &skip);
  char file_extra[LINE] = "";
  if (FileOption(argc, argv, "-fx", file_extra)) {
    exit(1);
  }
  // calculate pair-wise potential, etc.
  bool calculate_pairwise = BoolOption(argc, argv, "--pot");
  // file with potential parameters must be present if --pot is used
  if (calculate_pairwise && file_extra[0] == '\0') {
    strcpy(ERROR_MSG, "if --pot is used, -fx option identifying the potential "
           "is mandatory");
    PrintError();
    Help(argv[0], true, common, option);
    exit(1);
  }
  // use centre of mass (or centre of geometry) instead of box centre
  bool use_com = BoolOption(argc, argv, "--com");
  bool use_cog = BoolOption(argc, argv, "--cog");
  if (use_com && use_cog) {
    strcpy(ERROR_MSG, "both --com and --cog used; --com takes precedence");
    PrintWarning();
  } //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file,
                                coor_type, coor_file, detailed, pbc_xyz);
  COUNT *Count = &System.Count;
  BOX *box = &System.Box;

  // error - undefined mass when --com option is used //{{{
  bool mass = true;
  if (use_com) {
    for (int i = 0; i < Count->BeadType; i++) {
      if (System.BeadType[i].Mass == MASS) {
        mass = false;
        break;
      }
    }
    if (!mass) {
      strcpy(ERROR_MSG, "system contains beads with undefined mass; "
             "--com cannot be used");
      PrintError();
      exit(1);
    }
  } //}}}
  // array for parameters of the chosen potential
  double **ff_par[max_params];
  if (calculate_pairwise) {
    for (int i = 0; i < max_params; i++) {
      ff_par[i] = calloc(Count->BeadType, sizeof *ff_par);
      for (int j = 0; j < Count->BeadType; j++) {
        ff_par[i][j] = calloc(Count->BeadType, sizeof *ff_par);
      }
    }
  }

  char potential[LINE] = ""; // potential type
  double real_m = 1, // mass
         real_E = 1, // energy
         real_l = 1, // length
         real_T = 1, // temperature
         real_P = 1, // pressure
         CG = 1;     // coarse-graining level
  bool calculate_density = false; // for local density-dependent potentials
  // read extra information (if present) //{{{
  if (file_extra[0] != '\0') {
    FILE *fr = OpenFile(file_extra, "r");
    int line_count = 0;
    while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      line_count++;
      double val;
      if (words > 1 && split[0][0] != '#') {
        if (strcasecmp("potential", split[0]) == 0) {
          // only read the section if potential should be calculated
          if (calculate_pairwise) {
            strcpy(potential, split[1]);
            if (strcasecmp(potential, "ShiftedSC") == 0 ||
                strcasecmp(potential, "SuttonChen") == 0) {
              calculate_density = true;
              SuttonChenParameters(fr, file_extra, &line_count, ff_par, System);
            }
          } else {
            do {
              if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
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
          PrintWarnFileLine(file_extra, line_count);
        }
      }
    }
    fclose(fr);
    if (calculate_pairwise && potential[0] == '\0') {
      strcpy(ERROR_MSG, "missing potential information "
                   "(mandatory when --pot is used)");
      PrintErrorFile(file_extra, "\0", "\0");
      exit(1);
    }
  } //}}}
  // cross-terms for potentials //{{{
  if (calculate_pairwise) {
    if (strcasecmp(potential, "ShiftedSC") == 0 ||
        strcasecmp(potential, "SuttonChen") == 0) {
      CG = ff_par[10][0][0]; // CG level is the same for all beads
      double rcut = ff_par[5][0][0]; // cut-off is the same for all beads
      for (int i = 0; i < (Count->BeadType - 1); i++) {
        for (int j = (i + 1); j < Count->BeadType; j++) {
          ff_par[0][i][j] = (ff_par[0][i][i] + ff_par[0][j][j]) / 2; // m
          ff_par[1][i][j] = (ff_par[1][i][i] + ff_par[1][j][j]) / 2; // n
          ff_par[2][i][j] = sqrt(ff_par[2][i][i] * ff_par[2][j][j]); // a
          ff_par[3][i][j] = -1; // no cross-term in SC // c
          ff_par[4][i][j] = sqrt(ff_par[4][i][i] * ff_par[4][j][j]); // eps
          ff_par[5][i][j] = rcut; // rcut

          ff_par[0][j][i] = ff_par[0][i][j]; // m
          ff_par[1][j][i] = ff_par[1][i][j]; // n
          ff_par[2][j][i] = ff_par[2][i][j]; // a
          ff_par[3][j][i] = ff_par[3][i][j]; // c
          ff_par[4][j][i] = ff_par[4][i][j]; // eps
          ff_par[5][j][i] = ff_par[5][i][j]; // rcut
        }
      }
      // coarse-graining parameters of the potential //{{{
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = 0; j < Count->BeadType; j++) {
          ff_par[2][i][j] *= CG / real_l; // a
          ff_par[4][i][j] *= CUBE(CG) / real_T; // eps // TODO: always T?
          ff_par[5][i][j] *= CG / real_l; // rcut
        }
      }
    } //}}}
  } //}}}

  // coarse-graining of bead mass (or warn if unspecified mass)
  for (int i = 0; i < Count->BeadType; i++) {
    // warning - unspecified mass
    if (System.BeadType[i].Mass == MASS) {
      snprintf(ERROR_MSG, LINE, "unspecified mass (bead type %s%s%s)",
               ErrYellow(), System.BeadType[i].Name, ErrCyan());
      PrintWarnFile(struct_file, coor_file, "\0");
    } else {
      System.BeadType[i].Mass *= CUBE(CG) / real_m;
    }
  }

  if (verbose) {
    VerboseOutput(System);
  }

  // write header to output file //{{{
  PrintByline(out_global, argc, argv);
  FILE *fw = OpenFile(out_global, "a");
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
    // fprintf(fw, " (%d) uCR;", ++count);
    // fprintf(fw, " (%d) uCD;", ++count);
  }
  putc('\n', fw);
  fclose(fw); //}}}

  // variables for collecting observables //{{{
  bool *bt_in_use = calloc(Count->BeadType, sizeof *bt_in_use);
  long double (*KEnergy)[3] = NULL,
              (*Press_vir)[3] = NULL;
  long int (**BeadCount)[3] = NULL;
  // overall system-wide sums
  double sum_T = 0, // temperature
         sum_KE = 0, // kinetic energy
         sum_W = 0, // virial
         sum_Pxx = 0, sum_Pyy = 0, sum_Pzz = 0, // pressure tensor diagonal
         sum_PW = 0, sum_PKE = 0, // virial and KE pressure components
         sum_centre[3] = {0, 0, 0}; // centre of box/mass/geometry
  //}}}

  if (calculate_pairwise && strcasecmp(potential, "ShiftedSC") == 0) {
    SCPrefactors(System, ff_par);
  }

  // main loop //{{{
  FILE *coor = OpenFile(coor_file, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0, // count lines in the vcf file
      n_bins = 0; // number of bins is calculated at the first used timestep
  double *loc_dens = NULL;
  if (calculate_density) {
    loc_dens = calloc(Count->Bead, sizeof *loc_dens);
  }
  while (true) {
    PrintStep(&count_coor, start, silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if ((count_coor >= start && (count_coor <= end || end == -1)) &&
        ((count_coor - start) % skip) == 0) {
      use = true;
    }
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (!ReadTimestep(coor_type, coor, coor_file, &System, &line_count)) {
        count_coor--;
        break;
      }
      // define number of bins on the first used step //{{{
      if (count_coor == start) {
        int bin[3];
        bin[0] = ceil(box->Length.x / width);
        bin[1] = ceil(box->Length.y / width);
        bin[2] = ceil(box->Length.z / width);
        // make number of bins twice as large to be on the save side
        n_bins = Max3(bin[0], bin[1], bin[2]) * 2;
        if ((n_bins % 2) == 0) {
          n_bins++;
        }
        KEnergy = calloc(n_bins, sizeof *KEnergy),
        Press_vir = calloc(n_bins, sizeof *Press_vir);
        BeadCount = malloc(Count->BeadType * sizeof **BeadCount);
        for (int i = 0; i < Count->BeadType; i++) {
          BeadCount[i] = calloc(n_bins, sizeof *BeadCount[i]);
        }
      } //}}}
      count_used++;
      box = &System.Box;
      WrapJoinCoordinates(&System, true, false); // restore pbc
      // choose centre of the Cartesian coordinates //{{{
      VECTOR centre;
      if (use_com) {
        centre = CentreOfMass(Count->BeadCoor, System.BeadCoor, System);
      } else if (use_cog) {
        centre = GeomCentre(Count->BeadCoor, System.BeadCoor, System.Bead);
      } else {
        centre.x = box->Length.x / 2;
        centre.y = box->Length.y / 2;
        centre.z = box->Length.z / 2;
      } //}}}
      sum_centre[0] += centre.x;
      sum_centre[1] += centre.y;
      sum_centre[2] += centre.z;

      // calculate forces, pressure, etc. if desired //{{{
      for (int i = 0; i < Count->Bead; i++) {
        System.Bead[i].Force.x = 0;
        System.Bead[i].Force.y = 0;
        System.Bead[i].Force.z = 0;
      }
      double Pxx = 0, Pyy = 0, Pzz = 0;
      if (calculate_pairwise) {
        // calculete forces
        INTVECTOR n_cells;
        int *Head, *Link, Dcx[14], Dcy[14], Dcz[14];
        // access beads in the linked list through BeadCoor[]
        // TODO: rcut - pass maximum of ff_par[5][][], I guess
        LinkedList(System, &Head, &Link, ff_par[5][0][0], &n_cells, Dcx, Dcy, Dcz);
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
                      if (strcasecmp(potential, "ShiftedSC") == 0) {
                        LocalDensityShiftedSC(System, id_i, id_j, ff_par, loc_dens);
                      } else if (strcasecmp(potential, "SuttonChen") == 0) {
                        LocalDensitySC(System, id_i, id_j, ff_par, loc_dens);
                      }
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
                    bin_i[0] = (b_i->Position.x - centre.x) / width
                               + n_bins / 2;
                    bin_i[1] = (b_i->Position.y - centre.y) / width
                               + n_bins / 2;
                    bin_i[2] = (b_i->Position.z - centre.z) / width
                               + n_bins / 2;
                    bin_j[0] = (b_j->Position.x - centre.x) / width
                               + n_bins / 2;
                    bin_j[1] = (b_j->Position.y - centre.y) / width
                               + n_bins / 2;
                    bin_j[2] = (b_j->Position.z - centre.z) / width
                               + n_bins / 2; //}}}
                    VECTOR P_diag = {0, 0, 0};
                    if (strcasecmp(potential, "ShiftedSC") == 0) {
                      ShiftedSuttonChen(System, id_i, id_j, ff_par,
                                        loc_dens, &P_diag);
                    } else if (strcasecmp(potential, "SuttonChen") == 0) {
                      SuttonChen(System, id_i, id_j, ff_par, loc_dens, &P_diag);
                    }
                    Pxx += P_diag.x;
                    Pyy += P_diag.y;
                    Pzz += P_diag.z;
                    for (int k = 0; k < 3; k++) {
                      if (bin_i[k] == bin_j[k] &&
                          bin_i[k] >= 0 && bin_i[k] < n_bins) {
                        double press = P_diag.x + P_diag.y + P_diag.z;
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
        }
        free(Head);
        free(Link); //}}}
      } //}}}
      // calculate temperature, density, etc. //{{{
      double energy_kin = 0;
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        // kinetic enrgy
        double vel2 = SQR(b->Velocity.x) +
                      SQR(b->Velocity.y) +
                      SQR(b->Velocity.z);
        double ke = 0.5 * bt->Mass * vel2;
        energy_kin += ke;
        // bin id in the three directions
        int bin[3];
        bin[0] = (b->Position.x - centre.x) / width + n_bins / 2;
        bin[1] = (b->Position.y - centre.y) / width + n_bins / 2;
        bin[2] = (b->Position.z - centre.z) / width + n_bins / 2;
        for (int j = 0; j < 3; j++) {
          if (bin[j] >= 0 && bin[j] < n_bins) {
            BeadCount[b->Type][bin[j]][j]++;
            KEnergy[bin[j]][j] += ke;
          }
        }
        // specify the given bead type is used
        bt_in_use[b->Type] = true;
      } //}}}

      // global observables
      double temperature = 2 * energy_kin / (3 * Count->BeadCoor - 3),
             press_kin = 2 * energy_kin / box->Volume,
             virial = (Pxx + Pyy + Pzz) / 3,
             press_vir = virial / box->Volume;
      // write global values to output file //{{{
      fw = OpenFile(out_global, "a");
      fprintf(fw, " %5d", count_coor);
      fprintf(fw, " %e", temperature * real_T / CUBE(CG));
      fprintf(fw, " %e", energy_kin * real_E / CUBE(CG));
      if (calculate_pairwise) {
        fprintf(fw, " %e", press_kin * real_P);
        fprintf(fw, " %e", press_vir * real_P);
        fprintf(fw, " %e", (press_kin + press_vir) * real_P);
        fprintf(fw, " %e", Pxx);
        fprintf(fw, " %e", Pyy);
        fprintf(fw, " %e", Pzz);
        fprintf(fw, " %e", virial);
      }
      putc('\n', fw);
      fclose(fw); //}}}
      // add global observables to total sums //{{{
      sum_T += temperature;
      sum_KE += energy_kin;
      sum_PKE += press_kin;
      sum_W += virial;
      sum_PW += press_vir;
      sum_Pxx += Pxx;
      sum_Pyy += Pyy;
      sum_Pzz += Pzz; //}}}
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
  fclose(coor); //}}}
  if (count_coor == 0) { // error - input file without a valid timestep //{{{
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(coor_file, "\0", "\0");
    fputc('\n', stderr);           //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG, "no timestep used; "
           "starting timestep is higher than the number of timesteps");
    PrintWarning();     //}}}
  } else if (!silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
    fflush(stdout);
  } //}}}

  // print profiles of observables //{{{
  // open the three files & print headers //{{{
  char f[3][LINE]; // complete filenames
  if (snprintf(f[0], LINE, "%s-x.rho", out_local) < 0) {
    ErrorSnprintf();
  }
  if (snprintf(f[1], LINE, "%s-y.rho", out_local) < 0) {
    ErrorSnprintf();
  }
  if (snprintf(f[2], LINE, "%s-z.rho", out_local) < 0) {
    ErrorSnprintf();
  }
  FILE *fw2[3];
  for (int i = 0; i < 3; i++) {
    PrintByline(f[i], argc, argv);
    // print headers
    fw2[i] = OpenFile(f[i], "w");
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
  // find first bin to save in every direction //{{{
  // ...based on the geometric centre of the initial box and the averaged centre
  int first_bin[3] = {-1, -1, -1}, last_bin[3] = {-1, -1, -1};
  sum_centre[0] /= count_used;
  sum_centre[1] /= count_used;
  sum_centre[2] /= count_used;
  for (int k = 0; k < 3; k++) {
    first_bin[k] = n_bins / 2 - sum_centre[k] / width;
    last_bin[k] = n_bins / 2 + sum_centre[k] / width;
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
  for (int j = 0; j < 3; j++) {
    for (int i = first_bin[j]; i <= last_bin[j]; i++) {
      // absolute distance with 0 in the average centre of box/mass/geometry
      double abs_dist = width * (i + 0.5 - first_bin[j]) - sum_centre[j];
      fprintf(fw2[j], "%7.3f", abs_dist * real_l); // absolute distance
      // temperature - don't use 'count_used' as both KEnergy & BeadCount_tot
      //               are sums from all the counted steps
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
  } //}}}
  for (int i = 0; i < 3; i++) {
    fclose(fw2[i]);
  } //}}}

  // append overall averages to the 'out_global' file //{{{
  fw = OpenFile(out_global, "a");
  count = 1;
  fprintf(fw, "# (%d) T;", count++);
  fprintf(fw, " (%d) KE;", count++);
  if (calculate_pairwise) {
    fprintf(fw, " (%d) P(KE);", count++);
    fprintf(fw, " (%d) P(W);", count++);
    fprintf(fw, " (%d) P;", count++);
    fprintf(fw, " (%d) Pxx;", count++);
    fprintf(fw, " (%d) Pyy;", count++);
    fprintf(fw, " (%d) Pzz;", count++);
    fprintf(fw, " (%d) W;", count++);
  }
  putc('\n', fw);
  fprintf(fw, "#     ");
  fprintf(fw, " %e", sum_T * real_T / CUBE(CG) / count_used);
  fprintf(fw, " %e", sum_KE * real_E / CUBE(CG) / count_used);
  if (calculate_pairwise) {
    fprintf(fw, " %e", sum_PKE * real_P / count_used);
    fprintf(fw, " %e", sum_PW * real_P / count_used);
    fprintf(fw, " %e", (sum_PKE + sum_PW) * real_P / count_used);
    fprintf(fw, " %e", sum_Pxx / count_used);
    fprintf(fw, " %e", sum_Pyy / count_used);
    fprintf(fw, " %e", sum_Pzz / count_used);
    fprintf(fw, " %e", sum_W / count_used);
  }
  putc('\n', fw);
  fclose(fw); //}}}

  // free memory //{{{
  FreeSystem(&System);
  if (calculate_pairwise) {
    for (int i = 0; i < max_params; i++) {
      for (int j = 0; j < Count->BeadType; j++) {
        free(ff_par[i][j]);
      }
      free(ff_par[i]);
    }
  }
  for (int i = 0; i < Count->BeadType; i++) {
    free(BeadCount[i]);
  }
  free(BeadCount);
  free(KEnergy);
  free(Press_vir);
  free(BeadCount_tot);
  free(bt_in_use);
  if (calculate_density) {
    free(loc_dens);
  }
  //}}}

  return 0;
}
