#include "../AnalysisTools.h"

// TODO: reduced eps in SC - for now, the reduction is 230.5K (Al value)
//       Should I maybe use real units, so it can be reduced via red_E
//       Or should I add eps keyword to extra file (along with type of
//       potential, if I ever wanted to add one)
// TODO: stuff must be calculated relative to centre of mass becaus of the
//       variable-sized box otherwise which (presumably) gets bigger/smaller
//       in all directions, leaving the centre of mass roughly the same

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
  CommonHelp(error, n, opt);
} //}}}

typedef struct FF {
  double n, m,
         a, c, eps,
         rcut,
         rep1, rep2, dd1, dd2;
} FF;

// calculate some prefactors for the potential //{{{
void SCPrefactors(SYSTEM System, FF **sc) {
// prefactors for shifted potential
  for (int i = 0; i < System.Count.BeadType; i++) {
    for (int j = 0; j < System.Count.BeadType; j++) {
      // (a / r_c)^m * eps ...for repulsive part of SC
      sc[i][j].rep1 = pow((sc[i][j].a / sc[i][j].rcut), sc[i][j].m);
      sc[i][j].rep1 *= sc[i][j].eps;
      // m * eps * (a / r_c)^m / r_c ...for repulsive part of SC
      sc[i][j].rep2 = sc[i][j].rep1 * sc[i][j].m / sc[i][j].rcut;
      // (a / r_c)^n ...for local density
      sc[i][j].dd1 = pow((sc[i][j].a / sc[i][j].rcut), sc[i][j].n);
      // n * (a / r_c)^n / r_c ...for local density
      sc[i][j].dd2 = sc[i][j].dd1 * sc[i][j].n / sc[i][j].rcut;
    }
  }
} //}}}

// calculate local density for Sutton-Chen potential //{{{
void LocalDensitySC(SYSTEM System, int bead1, int bead2,
                    FF **sc, double *density) {
  BEAD *b_i = &System.Bead[bead1],
       *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type,
      bt_j = b_j->Type;
  if (rij < sc[bt_i][bt_j].rcut) {
    if (rij == 0) {
      strcpy(ERROR_MSG, "zero bead-bead distance!");
      PrintError();
      exit(1);
    }
    double a = sc[bt_i][bt_j].a;
    int n = sc[bt_i][bt_j].n;
    // density: (a / r_ij)^n ...
    double a_div_r = a / rij;
    double tmp = pow(a_div_r, n);
    // (shifted potential) ... - (a / r_c)^n ...
    tmp -= sc[bt_i][bt_j].dd1;
    // ... + (a / r_c)^n * n / r_c * (r_ij - r_c)
    tmp += sc[bt_i][bt_j].dd2 * (rij - sc[bt_i][bt_j].rcut);
    density[bead1] += tmp;
    density[bead2] += tmp;
    // density[id_i]++;
    // density[id_j]++;
  }
} //}}}

// calculate potential, force and wirial for shifted Sutton-Chen potential //{{{
void ShiftedSuttonChen(SYSTEM System, int bead1, int bead2, FF **sc,
                       double *density, double *potential, VECTOR *force,
                       double *Pxx, double *Pyy, double *Pzz) {
  BEAD *b_i = &System.Bead[bead1],
       *b_j = &System.Bead[bead2];
  VECTOR dist = Distance(b_i->Position, b_j->Position, System.Box.Length);
  double rij = VectorLength(dist);
  int bt_i = b_i->Type,
      bt_j = b_j->Type;
  if (rij < sc[bt_i][bt_j].rcut) {
    if (rij == 0) {
      strcpy(ERROR_MSG, "zero bead-bead distance!");
      PrintError();
      exit(1);
    }
    // helper value of (a / r_ij)^m & (a / r_ij)^n
    double a_div_r = sc[bt_i][bt_j].a / rij;
    double a_div_r_m = pow(a_div_r, sc[bt_i][bt_j].m),
           a_div_r_n = pow(a_div_r, sc[bt_i][bt_j].n);
    // 1) repulsive energy
    // eps * (a / r_ij)^m ...
    double tmp = sc[bt_i][bt_j].eps * a_div_r_m;
    // (shifted potential) ... - eps * (a / r_c)^m ...
    tmp -= sc[bt_i][bt_j].rep1;
    // ... + eps * (a / r_c)^m * m / r_c * (r_ij - r_c)
    tmp += sc[bt_i][bt_j].rep2 * (rij - sc[bt_i][bt_j].rcut);
    // add to per-particle repulsive energy
    potential[bead1] += tmp;
    potential[bead2] += tmp;
    // u_rep[id_i]++;
    // u_rep[id_j]++;
    // 2) repulsive virial
    // repulsive: -m * eps * (a / r_ij)^m ...
    tmp = -sc[bt_i][bt_j].m * sc[bt_i][bt_j].eps * a_div_r_m;
    // ... + m * eps * (a / r_c)^m / r_c * r_ij
    tmp += sc[bt_i][bt_j].rep2 * rij;
    double w_rep = tmp;
    // 3) density-dependent virial
    // 0.5 * (c_i * eps_ii / sqrt(\rho_i) +
    //        c_j * eps_jj / sqrt(\rho_j) ...
    double c_eps_i = sc[bt_i][bt_i].c * sc[bt_i][bt_i].eps /
                     sqrt(density[bead1]),
           c_eps_j = sc[bt_j][bt_j].c * sc[bt_j][bt_j].eps /
                     sqrt(density[bead2]);
    tmp = 0.5 * (c_eps_i + c_eps_j);
    // ... * (n * (a / r_ij)^n - n * (a / r_c)^n / r_c * r_ij)
    tmp *= sc[bt_i][bt_j].n * a_div_r_n - sc[bt_i][bt_j].dd2 * rij;
    double w_dd = tmp;
    // 4) pairwise forces
    double fij = -(w_rep + w_dd) / SQR(rij);
    // force in axes' directions
    VECTOR f;
    f.x = fij * dist.x;
    f.y = fij * dist.y;
    f.z = fij * dist.z;
    // pressure tensor components
    *Pxx += f.x * dist.x;
    *Pyy += f.y * dist.y;
    *Pzz += f.z * dist.z;
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

  int elements = 2; // Count->BeadType basically
  FF **sc = malloc(sizeof **sc * elements);
  for (int i = 0; i < elements; i++) {
    sc[i] = calloc(elements, sizeof *sc[i]);
  }
  // Al-Al
  sc[0][0].n = 6;
  sc[0][0].m = 8;
  sc[0][0].eps = 230.5;
  sc[0][0].c = 25.398;
  sc[0][0].a = 4.049;
  sc[0][0].rcut = 1;
  // Zr-Zr
  sc[1][1].n = 5;
  sc[1][1].m = 9;
  sc[1][1].eps = 4437.0;
  sc[1][1].c = 5.462;
  sc[1][1].a = 3.203;
  sc[0][0].rcut = 1;
  // Al-Zr
  sc[0][1].n = (sc[0][0].n + sc[1][1].n) / 2;
  sc[0][1].m = (sc[0][0].m + sc[1][1].m) / 2;
  sc[0][1].eps = sqrt(sc[0][0].eps * sc[1][1].eps);
  sc[0][1].c = -1; // no cross-term in SC
  sc[0][1].a = sqrt(sc[0][0].a * sc[1][1].a);
  sc[0][0].rcut = 1;
  // Zr-Al
  sc[1][0].n = sc[0][1].n;
  sc[1][0].m = sc[0][1].m;
  sc[1][0].eps = sc[0][1].eps;
  sc[1][0].c = sc[0][1].c;
  sc[1][0].a = sc[0][1].a;
  sc[0][0].rcut = 1;

  // define options //{{{
  int common = 12, all = common + 2, count = 0,
      req_arg = 2;
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
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output_gl> & <output>
  char out_global[LINE] = "", out_local[LINE] = "";
  snprintf(out_global, LINE, "%s", argv[++count]);
  snprintf(out_local, LINE, "%s", argv[++count]);
  out_local[LINE-6] = '\0'; // for adding -<axis>.rho

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int start= 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);
  int trash;
  char file_extra[LINE] = "";
  if (FileIntegerOption(argc, argv, 0, "-fx", &trash, &trash, file_extra)) {
    exit(1);
  }
  // TODO: add potential selection into the extra file; this will depend on it
  bool calculate_density = true; // for local density-dependent potentials
  // calculate pair-wise potential, etc.
  bool calculate_pairwise = BoolOption(argc, argv, "--pot");
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
  VECTOR tmp_bin;
  // TODO: *3 to assume box change of at most thrice as big
  //       remake to realloc on the fly
  tmp_bin.x = ceil(box->Length.x / width);
  tmp_bin.y = ceil(box->Length.y / width);
  tmp_bin.z = ceil(box->Length.z / width);
  double n_bins = Max3(tmp_bin.x, tmp_bin.y, tmp_bin.z); //}}}

  // coarse-graining & reducing parameters //{{{
  int CG = 1;
  // basis for reduced units
  double real_m = 1, // mass
         real_E = 1, // energy
         real_l = 1, // length
         real_T = 1, // temperature
         real_P = 1, // pressure
         rcut = 1;   // cut-off in real units (TODO: I guess?)
  // read extra information (if present) //{{{
  if (file_extra[0] != '\0') {
    FILE *fr = OpenFile(file_extra, "r");
    char line[LINE], *split[SPL_STR];
    int words, line_count = 0;
    while (ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      line_count++;
      double val;
      if (words > 1 && split[0][0] != '#') {
        if (!IsRealNumber(split[1], &val)) {
          goto warning;
        }
        if (strncasecmp("Energy", split[0], 2) == 0) {
          real_E = val;
        } else if (strcasecmp("Temperature", split[0]) == 0) {
          real_T = val;
        } else if (strcasecmp("Length", split[0]) == 0) {
          real_l = val;
        } else if (strcasecmp("Mass", split[0]) == 0) {
          real_m = val;
        } else if (strcasecmp("Pressure", split[0]) == 0) {
          real_P = val;
        } else if (strcasecmp("CG", split[0]) == 0) {
          CG = val;
        } else if (strcasecmp("rcut", split[0]) == 0) {
          rcut = val;
        } else {
          warning:
            strcpy(ERROR_MSG, "unrecognized line");
            PrintWarningFileLine(file_extra, line_count, split, words);
        }
      }
    }
    fclose(fr);
  } //}}}
  for (int i = 0; i < elements; i++) {
    for (int j = 0; j < elements; j++) {
      sc[i][j].eps *= CUBE(CG);
      sc[i][j].a *= CG;
      sc[i][j].rcut = rcut * CG;
    }
  }
  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Mass *= CUBE(CG);
  }
  for (int i = 0; i < elements; i++) {
    for (int j = 0; j < elements; j++) {
      sc[i][j].a /= real_l;
      // sc[i][j].eps /= real_E;
      sc[i][j].eps /= 230.5; // TODO: generalize somehow
      sc[i][j].rcut /= real_l;
    }
  }
  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Mass /= real_m;
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
  fprintf(fw, " (%d) P(KE);", ++count);
  fprintf(fw, " (%d) P(W);", ++count);
  fprintf(fw, " (%d) P(W+KE);", ++count);
  fprintf(fw, " (%d) Pxx;", ++count);
  fprintf(fw, " (%d) Pyy;", ++count);
  fprintf(fw, " (%d) Pzz;", ++count);
  fprintf(fw, " (%d) W;", ++count);
  putc('\n', fw);
  fclose(fw); //}}}

  SCPrefactors(System, sc);

  // variables for collecting observables //{{{
  LONGVECTOR *KEnergy = calloc(Count->Bead, sizeof *KEnergy);
  LONGINTVECTOR **BeadCount = malloc(Count->BeadType * sizeof **BeadCount);
  bool *bt_in_use = calloc(Count->BeadType, sizeof *bt_in_use);
  for (int i = 0; i < Count->BeadType; i++) {
    BeadCount[i] = calloc(n_bins, sizeof *BeadCount[i]);
  }
  double temp_sum = 0,
         press_kin_sum = 0,
         virial_sum = 0,
         press_vir_sum = 0;
  //}}}

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
      if (!ReadTimestep(coor_type, coor, coor_file, &System,
                        &line_count, vtf_var)) {
        count_coor--;
        break;
      }
      // some variables (may not be used); TODO: clear names
      double *loc_dens = calloc(Count->Bead, sizeof *loc_dens);
      double *uSC = calloc(Count->Bead, sizeof *uSC),
             Pxx = 0, Pyy = 0, Pzz = 0,
             energy_kin = 0;
      VECTOR *force = calloc(Count->Bead, sizeof *force);
      bool *test_used = calloc(Count->Bead, sizeof *test_used);

      count_used++;
      box = &System.Box;
      WrapJoinCoordinates(&System, true, false); // restore pbc
      if (calculate_pairwise) {
        // calculete forces
        INTVECTOR n_cells;
        int *Head, *Link, Dcx[14], Dcy[14], Dcz[14];
        // access beads in the linked list through BeadCoor[]
        // TODO: rcut - pass maximum of sc[][].rcut, I guess
        LinkedList(System, &Head, &Link, sc[0][0].rcut, &n_cells, Dcx, Dcy, Dcz);
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
                    int cell2 = c2x +
                                c2y * n_cells.x +
                                c2z * n_cells.x * n_cells.y; //}}}
                    // select bead in the cell 'cell2' //{{{
                    int j;
                    if (cell1 == cell2) { // next bead in 'cell1'
                      j = Link[i];
                    } else { // first bead in 'cell2'
                      j = Head[cell2];
                    } //}}}
                    while (j != -1) {
                      int id_i = System.BeadCoor[i],
                          id_j = System.BeadCoor[j];
                      LocalDensitySC(System, id_i, id_j, sc, loc_dens);
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
                  int cell2 = c2x +
                              c2y * n_cells.x +
                              c2z * n_cells.x * n_cells.y; //}}}
                  // select bead in the cell 'cell2' //{{{
                  int j;
                  if (cell1 == cell2) { // next bead in 'cell1'
                    j = Link[i];
                  } else { // first bead in 'cell2'
                    j = Head[cell2];
                  } //}}}
                  while (j != -1) {
                    int id_j = System.BeadCoor[j];
                    ShiftedSuttonChen(System, id_i, id_j, sc, loc_dens, uSC,
                                      force, &Pxx, &Pyy, &Pzz);
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
      }
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
        // Sutton-Chen potential
        double ui = sc[b->Type][b->Type].c * sc[b->Type][b->Type].eps *
                    sqrt(loc_dens[id]);
        uSC[id] = 0.5 * uSC[id] - ui;
        // copy previously calculated forces
        b->Force = force[id];
        // bin id in the three directions
        INTVECTOR bin;
        bin.x = b->Position.x / width;
        bin.y = b->Position.y / width;
        bin.z = b->Position.z / width;
        // realloc memory if bin is too high
        if (Max3(bin.x, bin.y, bin.z) >= n_bins) {
          n_bins += 10;
          KEnergy = realloc(KEnergy, n_bins * sizeof *KEnergy);
          for (int j = 0; j < Count->BeadType; j++) {
            BeadCount[j] = realloc(BeadCount[j], n_bins * sizeof *BeadCount[j]);
          }
        }
        // count beads in bins for density, temperature
        BeadCount[b->Type][bin.x].x++;
        BeadCount[b->Type][bin.y].y++;
        BeadCount[b->Type][bin.z].z++;
        // kinetic energy for temperature
        KEnergy[bin.x].x += ke;
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
      fprintf(fw, " %8.3e", temperature * real_T / CUBE(CG));
      fprintf(fw, " %10.3e", energy_kin * real_E / CUBE(CG));
      fprintf(fw, " %10.3e", press_kin * real_P);
      fprintf(fw, " %10.3e", press_vir * real_P);
      fprintf(fw, " %10.3e", (press_kin + press_vir) * real_P);
      fprintf(fw, " %10.3e", Pxx);
      fprintf(fw, " %10.3e", Pyy);
      fprintf(fw, " %10.3e", Pzz);
      fprintf(fw, " %10.3e", virial);
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
      if (!SkipTimestep(coor_type, coor, coor_file,
                        struct_file, &line_count)) {
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
    fputc('\n', stderr); //}}}
  } else if (start > count_coor) { // warn if no timesteps were written //{{{
    strcpy(ERROR_MSG, "no timestep used; "
           "starting timestep is higher than the number of timesteps");
    PrintWarning(); //}}}
  } else if (!silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
    fflush(stdout);
  } //}}}

  char fx[LINE], fy[LINE], fz[LINE];
  if (snprintf(fx, LINE, "%s-x.rho", out_local) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    PrintError();
    exit(1);
  }
  if (snprintf(fy, LINE, "%s-y.rho", out_local) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    PrintError();
    exit(1);
  }
  if (snprintf(fz, LINE, "%s-z.rho", out_local) < 0) {
    strcpy(ERROR_MSG, "something wrong with snprintf()");
    PrintError();
    exit(1);
  }
  FILE *fwx = OpenFile(fx, "w");
  FILE *fwy = OpenFile(fy, "w");
  FILE *fwz = OpenFile(fz, "w");
  PrintByline(fwx, argc, argv);
  PrintByline(fwy, argc, argv);
  PrintByline(fwz, argc, argv);
  fprintf(fwx, "# columns: (1) distance");
  fprintf(fwy, "# columns: (1) distance");
  fprintf(fwz, "# columns: (1) distance");
  count = 1;
  for (int i = 0; i < Count->BeadType; i++) {
    if (bt_in_use[i]) {
      count++;
      fprintf(fwx, "; (%d) %s", count, System.BeadType[i].Name);
      fprintf(fwy, "; (%d) %s", count, System.BeadType[i].Name);
      fprintf(fwz, "; (%d) %s", count, System.BeadType[i].Name);
    }
  }
  putc('\n', fwx);
  putc('\n', fwy);
  putc('\n', fwz);
  // write rdf
  VECTOR volume;
  volume.x = width * box->Length.y * box->Length.z;
  volume.y = box->Length.x * width * box->Length.z;
  volume.z = box->Length.x * box->Length.y * width;
  for (int i = 0; i < (n_bins - 1); i++) {
    double dist = width * (2 * i + 1) / 2;
    fprintf(fwx, "%7.3f", dist); // absolute distance
    fprintf(fwy, "%7.3f", dist); // absolute distance
    fprintf(fwz, "%7.3f", dist); // absolute distance
    for (int j = 0; j < Count->BeadType; j++) {
      if (bt_in_use[j]){
        VECTOR rho;
        rho.x = BeadCount[j][i].x / (volume.x * count_used);
        rho.y = BeadCount[j][i].y / (volume.y * count_used);
        rho.z = BeadCount[j][i].z / (volume.z * count_used);
        fprintf(fwx, " %10f", rho.x);
        fprintf(fwy, " %10f", rho.y);
        fprintf(fwz, " %10f", rho.z);
      }
    }
    putc('\n',fwx);
    putc('\n',fwy);
    putc('\n',fwz);
  }
  fclose(fwx);
  fclose(fwy);
  fclose(fwz);

  // free memory //{{{
  FreeSystem(&System);
  for (int i = 0; i < elements; i++) {
    free(sc[i]);
  }
  free(sc);
  for (int i = 0; i < Count->BeadType; i++) {
    free(BeadCount[i]);
  }
  free(BeadCount);
  free(KEnergy);
  free(bt_in_use);
  //}}}

  return 0;
}
