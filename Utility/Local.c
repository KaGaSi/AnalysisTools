#include "../AnalysisTools.h"

typedef struct SC {
  double n, m,
         a, c, eps;
} SC;

static void WriteData(FILE *fw, double **T, double **KE,
                      double **P, double **W) {
}

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "TEXT TO BE ADDED\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <axis> <with> <out1> <out2> [options]\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file"
          " (vcf, vtf, lammpstrj, or xyz format)\n");
  fprintf(ptr, "   <axis>            1 or 2 axis for 1D or 2D calculation"
          " (e.g., x or xy)\n");
  fprintf(ptr, "   <width>           width of a single bin"
          " (in the 1 or both dimensions)\n");
  fprintf(ptr, "   <out1>            output file with local variables\n");
  fprintf(ptr, "   <out2>            output file with per-timestep values\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --per-step     save a new file for each step"
          "(adds '-<step>.txt' to <out1>)\n");
  fprintf(ptr, "      -avg <n>       number of adjacent bins to average\n");
  fprintf(ptr, "      -fx <file>     file with extra information\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  char out_traj[LINE] = "out2.lammpstrj";
  FILE *fw = OpenFile(out_traj, "w");
  FILE *fw_data = OpenFile("data.txt", "w");
  fclose(fw);
  int CG = 4,
      elements = 2; // Count->BeadType basically
  double rcut = 12 * CG;
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
  for (int i = 0; i < elements; i++) {
    for (int j = 0; j < elements; j++) {
      sc[i][j].eps *= CUBE(CG);
      sc[i][j].a *= CG;
    }
  }

  // define options //{{{
  int common = 12, all = common + 3, count = 0,
      req_arg = 5;
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
  strcpy(option[count++], "--per-step");
  strcpy(option[count++], "-avg");
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

  // <axis> //{{{
  count++;
  bool err = false;
  int dim[3] = {-1, -1, -1}; // [0]..1/2 for 1D/2D, [1..2]..0, 1, 2 for x, y, z
  if (argv[count][0] == 'x') {
    dim[1] = 0;
  } else if (argv[count][0] == 'y') {
    dim[1] = 1;
  } else if (argv[count][0] == 'z') {
    dim[1] = 2;
  } else {
    err = true;
  }
  if (argv[count][1] == '\0') {
    dim[0] = 1;
  } else {
    dim[0] = 2;
    if (argv[count][1] == 'x') {
      dim[2] = 0;
    } else if (argv[count][1] == 'y') {
      dim[2] = 1;
    } else if (argv[count][1] == 'z') {
      dim[2] = 2;
  } else {
    err = true;
    }
  }
  if (dim[1] == dim[2]) {
    err = true;
  }
  if (err) {
    strcpy(ERROR_MSG, "<axis> argument must contain axis labels (x, y, z)");
    PrintError();
    Help(argv[0], true, common, option);
  } //}}}

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  double width[2];
  if (!IsPosRealNumber(argv[++count], &width[0])) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  }
  width[1] = width[0]; // TODO: for now
  //}}}

  // <out1> & <out2>
  char out1_file[LINE] = "", out2_file[LINE] = "";
  snprintf(out1_file, LINE, "%s", argv[++count]);
  snprintf(out2_file, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int start= 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);
  bool per_step = BoolOption(argc, argv, "--per-step");
  int arr[1] = {0}, trash;
  if (IntegerOption(argc, argv, 1, "-avg", &trash, arr)) {
    exit(1);
  }
  int avg = arr[0];
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

  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Mass *= CUBE(CG);
  }
  bool *write = malloc(Count->Bead * sizeof *write);
  for (int i = 0; i < Count->Bead; i++) {
    write[i] = true;
  }

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
        } else {
          warning:
            strcpy(ERROR_MSG, "unrecognized line");
            PrintWarningFileLine(file_extra, line_count, split, words);
        }
      }
    }
    fclose(fr);
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

  // number of bins //{{{
  // TODO: variable box size? For now, assume at most twice the size?
  //       nah, just make the bins variable width based on instantenous box size
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintErrorFile(coor_file, struct_file, "\0"); // TODO add file(s)
    exit(1);
  }
  BOX *box = &System.Box;
  int bin[2] = {-1, 1};
  for (int i = 0; i < dim[0]; i++) {
    if (dim[i+1] == 0) {
      bin[i] = box->Length.x / width[i] + 1;
    } else if (dim[i+1] == 1) {
      bin[i] = box->Length.y / width[i] + 1;
    } else if (dim[i+1] == 2) {
      bin[i] = box->Length.z / width[i] + 1;
    }
  }
  if (bin[1] == 1) { // 1D case // TODO: change?
    width[1] = 2 * Max3(box->Length.x, box->Length.y, box->Length.z) + 1;
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

  // print initial stuff to output file (if not -per-step) //{{{
  if (!per_step) {
    FILE *out = OpenFile(out1_file, "w");
    PrintByline(out, argc, argv);
    fclose(out);
  }
  fw = OpenFile(out2_file, "w");
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
  putc('\n', fw);
  fclose(fw);
  //}}}

  // main loop //{{{
  FILE *coor = OpenFile(coor_file, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file

  // when 1D is calculated, only Temp[i][0] is used
  long **bead_count = calloc(bin[0], sizeof **bead_count); // count beads
  double **KE = calloc(bin[0], sizeof **KE); // kinetic energy
  double **sum_KE = calloc(bin[0], sizeof **sum_KE);
  double **Temp = calloc(bin[0], sizeof **Temp); // temperature
  double **sum_Temp = calloc(bin[0], sizeof **sum_Temp);
  double **vir = calloc(bin[0], sizeof **vir); // virial
  double **sum_vir = calloc(bin[0], sizeof **sum_vir);
  double **sum_beadcount = calloc(bin[0], sizeof **sum_beadcount);
  for (int i = 0; i < bin[0]; i++) {
    bead_count[i] = calloc(bin[1], sizeof *bead_count);
    KE[i] = calloc(bin[1], sizeof *KE);
    sum_KE[i] = calloc(bin[1], sizeof *sum_KE);
    Temp[i] = calloc(bin[1], sizeof *Temp);
    sum_Temp[i] = calloc(bin[1], sizeof *sum_Temp);
    vir[i] = calloc(bin[1], sizeof *vir);
    sum_vir[i] = calloc(bin[1], sizeof *sum_vir);
    sum_beadcount[i] = calloc(bin[1], sizeof *sum_beadcount);
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
      if (!ReadTimestep(coor_type, coor, coor_file, &System,
                        &line_count, vtf_var)) {
        count_coor--;
        break;
      }
      count_used++;
      box = &System.Box;
      // define width based on instantenous box size //{{{
      for (int j = 0; j < 2; j++) {
        if (dim[j+1] == 0) {
          width[j] = box->Length.x / bin[j];
        } else if (dim[j+1] == 1) {
          width[j] = box->Length.y / bin[j];
        } else if (dim[j+1] == 2) {
          width[j] = box->Length.z / bin[j];
        }
      } //}}}
      WrapJoinCoordinates(&System, true, false); // restore pbc
      InitLong2DArray(bead_count, bin[0], bin[1], 0);
      InitDouble2DArray(Temp, bin[0], bin[1], 0);
      InitDouble2DArray(KE, bin[0], bin[1], 0);
      // calculate the local properties //{{{
      double temperature = 0, energy_kin = 0, virial = 0;
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        int pos[2] = {0, 0};
        for (int j = 0; j < 2; j++) {
          if (dim[j+1] == 0) {
            pos[j] = b->Position.x / width[j];
          } else if (dim[j+1] == 1) {
            pos[j] = b->Position.y / width[j];
          } else if (dim[j+1] == 2) {
            pos[j] = b->Position.z / width[j];
          }
        }
        // printf("width %lf %lf pos %3d %3d centre %7.3f %7.3f %7.3f "
        //        "box %7.3f %7.3f %7.3f\n", width[0], width[1],
        //        pos[0], pos[1], centre.x, centre.y, centre.z,
        //        System.Box.Length.x, System.Box.Length.y, System.Box.Length.z);
        // printf(" %5d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
        //        id, b->Position.x, b->Position.y, b->Position.z,
        //        b->Velocity.x, b->Velocity.y, b->Velocity.z,
        //        b->Force.x, b->Force.y, b->Force.z);
        bead_count[pos[0]][pos[1]]++;
        double vel2 = SQR(b->Velocity.x) +
                      SQR(b->Velocity.y) +
                      SQR(b->Velocity.z);
        KE[pos[0]][pos[1]] += 0.5 * bt->Mass * vel2;
        vir[pos[0]][pos[1]] += b->Force.x * b->Position.x +
                               b->Force.y * b->Position.y +
                               b->Force.z * b->Position.z;
        energy_kin += 0.5 * bt->Mass * vel2;
        virial += b->Force.x * b->Position.x +
                  b->Force.y * b->Position.y +
                  b->Force.z * b->Position.z;
      }
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          sum_KE[i][j] += KE[i][j];
          Temp[i][j] = 2 * KE[i][j] / 3;
          sum_Temp[i][j] += Temp[i][j];
          sum_vir[i][j] += vir[i][j];
          sum_beadcount[i][j] += bead_count[i][j];
        }
      }
      temperature = 2 * energy_kin / (3 * Count->BeadCoor - 3);
      double press_kin = 2 * energy_kin / (3 * box->Volume);
      double press_vir = virial / (3 * box->Volume);
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
      WriteTimestep(LTRJ_FILE, out_traj, System, count_coor, write);
      virial = (Pxx + Pyy + Pzz) / 3;
      // printf("\nP: %lf %lf %lf %lf\n", Pxx, Pyy, Pzz, virial);
      press_vir = virial / box->Volume;
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        fprintf(fw_data, "%d %lf %lf\n", id + 1, uSC[id], density[id]);
      }

      fw = OpenFile(out2_file, "a");
      fprintf(fw, " %8d", count_coor);
      fprintf(fw, " %8.3e", temperature * real_T / CUBE(CG));
      fprintf(fw, " %10.3e", energy_kin * real_E / CUBE(CG));
      fprintf(fw, " %10.3e", press_kin * real_P);
      fprintf(fw, " %10.3e", press_vir * real_P);
      fprintf(fw, " %10.3e", (press_kin + press_vir) * real_P);
      fprintf(fw, " %10.3e", Pxx);
      fprintf(fw, " %10.3e", Pyy);
      fprintf(fw, " %10.3e", Pzz);
      putc('\n', fw);
      fclose(fw);

      free(density);
      free(uSC);
      free(force);
      free(Head);
      free(Link);
      // int avg = 0; // use data from bins -avg to avg
      // TODO: add some averaging: sum up bins x-avg to x+avg, plot to x
      // save values (if using --pers-step) //{{{
      if (per_step) {
        char file[LINE];
        if (snprintf(file, LINE, "%s-%04d.txt", out1_file, count_coor) < 0) {
          strcpy(ERROR_MSG, "something wrong with snprintf()");
          PrintErrorFile(file, "\0", "\0");
          exit(1);
        }
        FILE *fw2 = OpenFile(file, "w");
        // 1D case //{{{
        if (dim[0] == 1) {
          double volume = (2 * avg + 1) * width[0];
          char *axis;
          if (dim[1] == 0) { //{{{
            axis = "x";
            volume *= System.Box.Length.y * System.Box.Length.z;
          } else if (dim[1] == 1) {
            axis = "y";
            volume *= System.Box.Length.x * System.Box.Length.z;
          } else {
            axis = "z";
            volume *= System.Box.Length.x * System.Box.Length.y;
          } //}}}
          count = 0;
          fprintf(fw2, "# (%d) %s", ++count, axis);
          fprintf(fw2, "; (%d) T", ++count);
          fprintf(fw2, "; (%d) KE", ++count);
          fprintf(fw2, "; (%d) P", ++count);
          fprintf(fw2, "; (%d) W", ++count);
          fprintf(fw2, "; (%d) avg beads in the bin", ++count);
          putc('\n', fw2);
          for (int i = 0; i < bin[0]; i++) {
            double avg_KE = 0, avg_T = 0, avg_vir = 0;
            int avg_count = 0;
            for (int j = (i - avg); j <= (i + avg); j++) {
              int k = j;
              if (k < 0) {
                k = bin[0] + k;
              } else if (k >= bin[0]) {
                k = k - bin[0];
              }
              avg_KE += KE[k][0];
              avg_T += Temp[k][0];
              avg_vir += sum_vir[k][0];
              avg_count += bead_count[k][0];
            }
            int boxes = (2 * avg + 1);
            fprintf(fw2, " %8.3f", width[0] * (2 * i + 1) / 2);
            fprintf(fw2, " %8.3f", avg_T / avg_count * real_T);
            fprintf(fw2, " %8.3f", avg_KE / boxes * real_E);
            fprintf(fw2, " %8.3f", (avg_T  + avg_vir / 3) / volume);
            fprintf(fw2, " %8.3f", avg_vir / avg_count);
            fprintf(fw2, " %8.3f", (double)(avg_count) / boxes);
            putc('\n', fw2);
          }
        //}}}
        // 2D case //{{{
        } else {
          double volume = SQR(2 * avg + 1) * width[0] * width[0];
          char *axis[2];
          if (dim[1] == 0) {
            axis[0] = "x";
            if (dim[2] == 1) {
              axis[1] = "y";
              volume *= System.Box.Length.z;
            } else {
              axis[1] = "z";
              volume *= System.Box.Length.y;
            }
          } else if (dim[1] == 1) {
            axis[0] = "y";
            if (dim[2] == 0) {
              axis[1] = "x";
              volume *= System.Box.Length.z;
            } else {
              axis[1] = "z";
              volume *= System.Box.Length.x;
            }
          } else {
            axis[0] = "z";
            if (dim[2] == 2) {
              axis[1] = "x";
              volume *= System.Box.Length.y;
            } else {
              axis[1] = "y";
              volume *= System.Box.Length.x;
            }
          }
          count = 0;
          fprintf(fw2, "# (%d) %s", ++count, axis[0]);
          fprintf(fw2, "; (%d) %s", ++count, axis[1]);
          fprintf(fw2, "; (%d) T", ++count);
          fprintf(fw2, "; (%d) KE", ++count);
          fprintf(fw2, "; (%d) P", ++count);
          fprintf(fw2, "; (%d) W", ++count);
          fprintf(fw2, "; (%d) avg beads in the bin", ++count);
          putc('\n', fw2);
          for (int i = 0; i < bin[0]; i++) {
            for (int j = 0; j < bin[1]; j++) {
              double avg_KE = 0, avg_T = 0, avg_vir = 0;
              int avg_count = 0;
              for (int k = (i-avg); k <= (i+avg); k++) {
                int m = k;
                if (m < 0) {
                  m = bin[0] + m;
                } else if (m >= bin[0]) {
                  m = m - bin[0];
                }
                for (int l = (j-avg); l <= (j+avg); l++) {
                  int n = l;
                  if (n < 0) {
                    n = bin[1] + n;
                  } else if (n >= bin[1]) {
                    n = n - bin[1];
                  }
                  avg_KE += KE[m][n];
                  avg_T += Temp[m][n];
                  avg_vir += vir[m][n];
                  avg_count += bead_count[m][n];
                }
              }
              int boxes = SQR(2 * avg + 1);
              fprintf(fw2, " %8.3f", width[0] * (2 * i + 1) / 2);
              fprintf(fw2, " %8.3f", width[1] * (2 * j + 1) / 2);
              fprintf(fw2, " %8.3f", avg_T / avg_count * real_T);
              fprintf(fw2, " %8.3f", avg_KE / boxes * real_E);
              fprintf(fw2, " %8.3f", (avg_T + avg_vir / 3) / volume);
              fprintf(fw2, " %8.3f", avg_vir / avg_count);
              fprintf(fw2, " %8.3f\n", (double)(avg_count) / boxes);
              putc('\n', fw2);
            }
            putc('\n', fw2);
          }
        }
        //}}}
        fclose(fw2);
      } //}}}
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

  // save average (if not using --pers-step) //{{{
  if (!per_step) {
    fw = OpenFile(out1_file, "a");
    // 1D case //{{{
    if (dim[0] == 1) {
      double volume = (2 * avg + 1) * width[0];
      char *axis;
      if (dim[1] == 0) {
        axis = "x";
        volume *= System.Box.Length.y * System.Box.Length.z;
      } else if (dim[1] == 1) {
        axis = "y";
        volume *= System.Box.Length.x * System.Box.Length.z;
      } else {
        axis = "z";
        volume *= System.Box.Length.x * System.Box.Length.y;
      }
      count = 0;
      fprintf(fw, "# (%d) %s", ++count, axis);
      fprintf(fw, "; (%d) T", ++count);
      fprintf(fw, "; (%d) KE", ++count);
      fprintf(fw, "; (%d) P", ++count);
      fprintf(fw, "; (%d) W", ++count);
      fprintf(fw, "; (%d) avg beads in the bin", ++count);
      putc('\n', fw);
      // int avg = 1;
      for (int i = 0; i < bin[0]; i++) {
        double avg_KE = 0, avg_T = 0, avg_vir = 0;
        int avg_count = 0;
        for (int j = (i-avg); j <= (i+avg); j++) {
          int k = j;
          if (k < 0) {
            k = bin[0] + k;
          } else if (k >= bin[0]) {
            k = k - bin[0];
          }
          avg_KE += sum_KE[k][0];
          avg_T += sum_Temp[k][0];
          avg_vir += sum_vir[k][0];
          avg_count += sum_beadcount[k][0];
        }
        int boxes = (2 * avg + 1) * count_used;
        fprintf(fw, " %8.3f", width[0] * (2 * i + 1) / 2);
        fprintf(fw, " %8.3f", avg_T / avg_count * real_T);
        fprintf(fw, " %8.3f", avg_KE / boxes * real_E);
        fprintf(fw, " %8.3f", (avg_T + avg_vir / 3) / volume);
        fprintf(fw, " %8.3f", avg_vir / avg_count);
        fprintf(fw, " %8.3f", (double)(avg_count) / boxes);
        putc('\n', fw);
      }
    //}}}
    // 2D case //{{{
    } else {
      double volume = SQR(2 * avg + 1) * width[0] * width[0];
      char *axis[2]; //{{{
      if (dim[1] == 0) {
        axis[0] = "x";
        if (dim[2] == 1) {
          axis[1] = "y";
          volume *= System.Box.Length.z;
        } else {
          axis[1] = "z";
          volume *= System.Box.Length.y;
        }
      } else if (dim[1] == 1) {
        axis[0] = "y";
        if (dim[2] == 0) {
          axis[1] = "x";
          volume *= System.Box.Length.z;
        } else {
          axis[1] = "z";
          volume *= System.Box.Length.x;
        }
      } else {
        axis[0] = "z";
        if (dim[2] == 2) {
          axis[1] = "x";
          volume *= System.Box.Length.y;
        } else {
          axis[1] = "y";
          volume *= System.Box.Length.x;
        }
      } //}}}
      count = 0;
      fprintf(fw, "# (%d) %s", ++count, axis[0]);
      fprintf(fw, "; (%d) %s", ++count, axis[1]);
      fprintf(fw, "; (%d) T", ++count);
      fprintf(fw, "; (%d) KE", ++count);
      fprintf(fw, "; (%d) P", ++count);
      fprintf(fw, "; (%d) W", ++count);
      fprintf(fw, "; (%d) avg beads in the bin", ++count);
      putc('\n', fw);
      // int avg = 0;
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          double avg_KE = 0, avg_T = 0, avg_vir = 0;
          int avg_count = 0;
          for (int k = (i - avg); k <= (i + avg); k++) {
            int m = k;
            if (m < 0) {
              m = bin[0] + m;
            } else if (m >= bin[0]) {
              m = m - bin[0];
            }
            for (int l = (j - avg); l <= (j + avg); l++) {
              int n = l;
              if (n < 0) {
                n = bin[1] + n;
              } else if (n >= bin[1]) {
                n = n - bin[1];
              }
              avg_KE += sum_KE[m][n];
              avg_T += sum_Temp[m][n];
              avg_vir += sum_vir[m][n];
              avg_count += sum_beadcount[m][n];
            }
          }
          int boxes = SQR(2 * avg + 1) * count_used;
          fprintf(fw, " %8.3f", width[0] * (2 * i + 1) / 2);
          fprintf(fw, " %8.3f", width[1] * (2 * j + 1) / 2);
          fprintf(fw, " %8.3f", avg_T / avg_count * real_T);
          fprintf(fw, " %8.3f", avg_KE / boxes * real_E);
          fprintf(fw, " %8.3f", (avg_T + avg_vir / 3) / volume);
          fprintf(fw, " %8.3f", avg_vir / avg_count);
          fprintf(fw, " %8.3f\n", (double)(avg_count) / boxes);
        }
        fprintf(fw, "\n");
      }
    } //}}}
    fclose(fw);
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
  for (int i = 0; i < bin[0]; i++) {
    free(bead_count[i]);
    free(KE[i]);
    free(sum_KE[i]);
    free(Temp[i]);
    free(sum_Temp[i]);
    free(vir[i]);
    free(sum_vir[i]);
    free(sum_beadcount[i]);
  }
  free(bead_count);
  free(KE);
  free(sum_KE);
  free(Temp);
  free(sum_Temp);
  free(vir);
  free(sum_vir);
  free(sum_beadcount);
  free(write);
  FreeSystem(&System);
  //}}}

  fclose(fw_data);
  return 0;
}
