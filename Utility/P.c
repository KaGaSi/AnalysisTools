#include "../AnalysisTools.h"

typedef struct SC {
  double n, m,
         a, c, eps;
} SC;

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "TEXT TO BE ADDED\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out> [options]\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file"
          " (vcf, vtf, lammpstrj, or xyz format)\n");
  fprintf(ptr, "   <out>             output file\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --per-step     save a new file for each step"
          "(adds '-<step>.txt' to <out1>)\n");
  fprintf(ptr, "      -avg <n>       number of adjacent bins to average\n");
  fprintf(ptr, "      -fx <file>     file with extra information\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  int elements = 2; // Count->BeadType basically
  SC sc[elements][elements];
  // Al-Al
  sc[0][0].n = 6;
  sc[0][0].m = 8;
  sc[0][0].eps = 230.5;
  sc[0][0].c = 25.398;
  sc[0][0].a = 4.049;
  // Zr-Zr
  sc[1][1].n = 5;
  sc[1][1].m = 9;
  sc[1][1].eps = 4437.0;
  sc[1][1].c = 5.462;
  sc[1][1].a = 3.203;
  // Al-Zr
  sc[0][1].n = (sc[0][0].n + sc[1][1].n) / 2;
  sc[0][1].m = (sc[0][0].m + sc[1][1].m) / 2;
  sc[0][1].eps = sqrt(sc[0][0].eps * sc[1][1].eps);
  sc[0][1].c = -1; // no cross-term in SC
  sc[0][1].a = sqrt(sc[0][0].a * sc[1][1].a);
  // Zr-Al
  sc[1][0].n = sc[0][1].n;
  sc[1][0].m = sc[0][1].m;
  sc[1][0].eps = sc[0][1].eps;
  sc[1][0].c = sc[0][1].c;
  sc[1][0].a = sc[0][1].a;

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

  // <out1> & <out2>
  char out1_file[LINE] = "";
  snprintf(out1_file, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int start= 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);
  int trash;
  char file_extra[LINE] = "";
  if (FileIntegerOption(argc, argv, 0, "-fx", &trash, &trash, file_extra)) {
    exit(1);
  } //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz);

  COUNT *Count = &System.Count;

  bool *write = malloc(Count->Bead * sizeof *write);
  for (int i = 0; i < Count->Bead; i++) {
    write[i] = true;
  }

  int CG = 1;
  double rcut = 1;
  // basis for reduced units
  double real_m = 1, // mass
         real_E = 1, // energy
         real_l = 1, // length
         real_T = 1, // temperature
         real_P = 1; // pressure
  // read extra information (if present)
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
  }
  rcut *= CG;
  for (int i = 0; i < elements; i++) {
    for (int j = 0; j < elements; j++) {
      sc[i][j].eps *= CUBE(CG);
      sc[i][j].a *= CG;
    }
  }
  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Mass *= CUBE(CG);
  }
  // reduce parameters
  for (int i = 0; i < elements; i++) {
    for (int j = 0; j < elements; j++) {
      sc[i][j].a /= real_l;
      // sc[i][j].eps /= real_E;
      sc[i][j].eps /= 230.5;
    }
  }
  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Mass /= real_m;
  }
  rcut /= real_l;

  // calculate some prefactors for shifted potential //{{{
  // eps * (a / r_c)^m ...for repulsive part of SC
  double **rep1 = malloc(Count->BeadType * sizeof **rep1);
  // m * eps * (a / r_c)^m / r_c ...for repulsive part of SC
  double **rep2 = malloc(Count->BeadType * sizeof **rep2);
  // (a / r_c)^n ...for local density
  double **dd1 = malloc(Count->BeadType * sizeof **dd1);
  // n * (a / r_c)^n / r_c ...for local density
  double **dd2 = malloc(Count->BeadType * sizeof **dd2);
  for (int i = 0; i < Count->BeadType; i++) {
    rep1[i] = calloc(Count->BeadType, sizeof *rep1);
    rep2[i] = calloc(Count->BeadType, sizeof *rep2);
    dd1[i] = calloc(Count->BeadType, sizeof *dd1);
    dd2[i] = calloc(Count->BeadType, sizeof *dd2);
    for (int j = 0; j < Count->BeadType; j++) {
      rep1[i][j] = sc[i][j].eps * pow((sc[i][j].a / rcut), sc[i][j].m);
      rep2[i][j] = rep1[i][j] * sc[i][j].m / rcut;
      dd1[i][j] = pow((sc[i][j].a / rcut), sc[i][j].n);
      dd2[i][j] = dd1[i][j] * sc[i][j].n / rcut;
    }
  } //}}}

  BOX *box = &System.Box;
  if (box->Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintErrorFile(coor_file, struct_file, "\0"); // TODO add file(s)
    exit(1);
  }

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

  FILE *fw = OpenFile(out1_file, "w");
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
  fclose(fw);

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
      count_used++;
      box = &System.Box;
      WrapJoinCoordinates(&System, true, false); // restore pbc
      // calculate kinetic energy //{{{
      double energy_kin = 0;
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        double vel2 = SQR(b->Velocity.x) +
                      SQR(b->Velocity.y) +
                      SQR(b->Velocity.z);
        energy_kin += 0.5 * bt->Mass * vel2;
      }
      double temperature = 2 * energy_kin / (3 * Count->BeadCoor - 3);
      double press_kin = 2 * energy_kin / box->Volume;
      //}}}
      // calculete forces
      INTVECTOR n_cells;
      int *Head, *Link, Dcx[14], Dcy[14], Dcz[14];
      // access beads in the linked list through BeadCoor[]
      LinkedList(System, &Head, &Link, rcut, &n_cells, Dcx, Dcy, Dcz);
      // density //{{{
      double *density = calloc(Count->Bead, sizeof *density);
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
                  BEAD *b_i = &System.Bead[id_i],
                       *b_j = &System.Bead[id_j];
                  VECTOR dist = Distance(b_i->Position, b_j->Position, box->Length);
                  double rij = VectorLength(dist);
                  if (rij < rcut) {
                    if (rij == 0) {
                      strcpy(ERROR_MSG, "zero bead-bead distance!");
                      PrintError();
                      exit(1);
                    }
                    int bt_i = b_i->Type,
                        bt_j = b_j->Type;
                    double a = sc[bt_i][bt_j].a;
                    int n = sc[bt_i][bt_j].n;
                    // density: (a / r_ij)^n ...
                    double a_div_r = a / rij;
                    double tmp = pow(a_div_r, n);
                    // (shifted potential) ... - (a / r_c)^n ...
                    tmp -= dd1[bt_i][bt_j];
                    // ... + (a / r_c)^n * n / r_c * (r_ij - r_c)
                    tmp += dd2[bt_i][bt_j] * (rij - rcut);
                    density[id_i] += tmp;
                    density[id_j] += tmp;
                    // density[id_i]++;
                    // density[id_j]++;
                  }
                  j = Link[j];
                }
              }
              i = Link[i];
            }
          }
        }
      } //}}}
      // forces, energy, virial //{{{
      double *uSC = calloc(Count->Bead, sizeof *uSC),
             Pxx = 0, Pyy = 0, Pzz = 0;
      VECTOR *force = calloc(Count->Bead, sizeof *force);
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
                  BEAD *b_i = &System.Bead[id_i],
                       *b_j = &System.Bead[id_j];
                  VECTOR dist = Distance(b_i->Position, b_j->Position, box->Length);
                  double rij = VectorLength(dist);
                  if (rij < rcut) {
                    if (rij == 0) {
                      strcpy(ERROR_MSG, "zero bead-bead distance!");
                      PrintError();
                      exit(1);
                    }
                    int bt_i = b_i->Type,
                        bt_j = b_j->Type;
                    // helper value of (a / r_ij)^m & (a / r_ij)^n
                    double a_div_r = sc[bt_i][bt_j].a / rij;
                    double a_div_r_m = pow(a_div_r, sc[bt_i][bt_j].m),
                           a_div_r_n = pow(a_div_r, sc[bt_i][bt_j].n);
                    // 1) repulsive energy
                    // eps * (a / r_ij)^m ...
                    double tmp = sc[bt_i][bt_j].eps * a_div_r_m;
                    // (shifted potential) ... - eps * (a / r_c)^m ...
                    tmp -= rep1[bt_i][bt_j];
                    // ... + eps * (a / r_c)^m * m / r_c * (r_ij - r_c)
                    tmp += rep2[bt_i][bt_j] * (rij - rcut);
                    // add to per-particle repulsive energy
                    uSC[id_i] += tmp;
                    uSC[id_j] += tmp;
                    // u_rep[id_i]++;
                    // u_rep[id_j]++;
                    // 2) repulsive virial
                    // repulsive: -m * eps * (a / r_ij)^m ...
                    tmp = -sc[bt_i][bt_j].m * sc[bt_i][bt_j].eps * a_div_r_m;
                    // ... + m * eps * (a / r_c)^m / r_c * r_ij
                    tmp += rep2[bt_i][bt_j] * rij;
                    double w_rep = tmp;
                    // 3) density-dependent virial
                    // 0.5 * (c_i * eps_ii / sqrt(\rho_i) +
                    //        c_j * eps_jj / sqrt(\rho_j) ...
                    double c_eps_i = sc[bt_i][bt_i].c * sc[bt_i][bt_i].eps /
                                     sqrt(density[id_i]),
                           c_eps_j = sc[bt_j][bt_j].c * sc[bt_j][bt_j].eps /
                                     sqrt(density[id_j]);
                    tmp = 0.5 * (c_eps_i + c_eps_j);
                    // ... * (n * (a / r_ij)^n - n * (a / r_c)^n / r_c * r_ij)
                    tmp *= sc[bt_i][bt_j].n * a_div_r_n - dd2[bt_i][bt_j] * rij;
                    double w_dd = tmp;
                    // 4) pairwise forces
                    double fij = -(w_rep + w_dd) / SQR(rij);
                    // force in axes' directions
                    VECTOR f;
                    f.x = fij * dist.x;
                    f.y = fij * dist.y;
                    f.z = fij * dist.z;
                    // pressure tensor components
                    Pxx += f.x * dist.x;
                    Pyy += f.y * dist.y;
                    Pzz += f.z * dist.z;
                    // per-particle force
                    force[id_i].x += f.x;
                    force[id_i].y += f.y;
                    force[id_i].z += f.z;
                    force[id_j].x -= f.x;
                    force[id_j].y -= f.y;
                    force[id_j].z -= f.z;
                    printf("%6d %6d", id_i, id_j);
                    printf(" %lf %lf %lf", dist.x, dist.y, dist.z);
                    printf(" %lf %lf %lf", w_dd, w_rep, fij);
                    printf(" %lf %lf %lf\n", f.x, f.y, f.z);
                  }
                  j = Link[j];
                }
              }
              i = Link[i];
            }
          }
        }
      } //}}}
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        int bt = b->Type;
        double ui = sc[bt][bt].c * sc[bt][bt].eps * sqrt(density[id]);
        uSC[id] = 0.5 * uSC[id] - ui;
        b->Force = force[id];
      }
      double virial = (Pxx + Pyy + Pzz) / 3;
      double press_vir = virial / box->Volume;

      fw = OpenFile(out1_file, "a");
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
      fclose(fw);

      free(density);
      free(uSC);
      free(force);
      free(Head);
      free(Link);
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

  // free memory //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    free(rep1[i]);
    free(rep2[i]);
    free(dd1[i]);
    free(dd2[i]);
  }
  free(rep1);
  free(rep2);
  free(dd1);
  free(dd2);
  free(write);
  FreeSystem(&System);
  //}}}

  return 0;
}
