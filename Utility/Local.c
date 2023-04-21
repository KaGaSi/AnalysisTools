#include "../AnalysisTools.h"

typedef struct SC {
  int n, m, CG;
  double a, c, eps;
} SC;

typedef struct constant {
  // data from lammps/src/update.cpp
  double mvv2e, // convert mv^2 to E (eV)
         boltz, // Boltzmann constant in eV/K
         nktv2p;
} CONST;

// reduced length - CG: Al 3^3, Zr 2^3
// ...not sure whether this is in any sense right
static const double red_dist[2] = {1, 1};

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "TEXT TO BE ADDED\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <axis> <with> <output> [options]\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file"
          " (vcf, vtf, lammpstrj, or xyz format)\n");
  fprintf(ptr, "   <axis>            1 or 2 axis for 1D or 2D calculation"
          " (e.g., x or xy)\n");
  fprintf(ptr, "   <width>           width of a single bin"
          " (in the 1 or both dimensions)\n");
  fprintf(ptr, "   <output>          output file with local variables\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --metal        assume lammps 'metal' units\n");
  fprintf(ptr, "      --per-step     save a new file for each step"
          "(adds '-<step>.txt' to <output>)\n");
  fprintf(ptr, "      -avg <n>       number of adjacent bins to average\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 12, all = common + 3, count = 0,
      req_arg = 4;
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
  strcpy(option[count++], "--metal");
  strcpy(option[count++], "--per-step");
  strcpy(option[count++], "--avg");
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

  // <output> - output text file
  char out_file[LINE] = "";
  snprintf(out_file, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int start= 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);
  int style = 0; // reduced units
  if (BoolOption(argc, argv, "--metal")) {
    style = 1; // lammps 'metal' units
  }
  bool per_step = BoolOption(argc, argv, "--per-step");
  int arr[1] = {0}, trash;
  if (IntegerOption(argc, argv, 1, "-avg", &trash, arr)) {
    exit(1);
  }
  int avg = arr[0]; //}}}

  CONST c; //{{{
  if (style == 1) { // lammps metal units
    c.mvv2e = 1.0364269e-4;
    c.boltz = 8.617343e-5;
    c.nktv2p = 1.6021765e6;
  } else { // reduced units (lammps lj units)
    c.mvv2e = 1;
    c.boltz = 1;
    c.nktv2p = 1;
  } //}}}

  SC SC[2];
  // Al
  SC[0].n = 6;
  SC[0].m = 8;
  SC[0].CG = 3;
  // SC[0].a = 1;      // 4.049Å
  SC[0].a = 4.049;
  SC[0].eps = 1;    // 230.5K
  SC[0].eps = 1.986e-02;
  SC[0].c = 25.398; // dimensionless
  // Zr
  SC[1].n = 6;
  SC[1].m = 8;
  SC[1].CG = 2;
  SC[1].a = 1.118795;   // 4.530Å
  SC[1].eps = 0.850759; // 196.1K
  SC[1].c = 25.398;     // dimensionless

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz);
  VECTOR centre;
  centre.x = System.Box.Length.x / 2;
  centre.y = System.Box.Length.y / 2;
  centre.z = System.Box.Length.z / 2;

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

  // print initial stuff to output vcf file (if not -per-step) //{{{
  if (!per_step) {
    FILE *out = OpenFile(out_file, "w");
    PrintByline(out, argc, argv);
    fclose(out);
  } //}}}

  // main loop //{{{
  FILE *coor = OpenFile(coor_file, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file

  // when 1D is calculated, only Temp[i][0] is used
  long **bead_count = calloc(bin[0], sizeof **bead_count); // count beads
  double **Temp = calloc(bin[0], sizeof **Temp); // temperature
  double **sum_Temp = calloc(bin[0], sizeof **sum_Temp);
  double **vir = calloc(bin[0], sizeof **vir); // virial
  double **sum_vir = calloc(bin[0], sizeof **sum_vir);
  double **sum_beadcount = calloc(bin[0], sizeof **sum_beadcount);
  for (int i = 0; i < bin[0]; i++) {
    bead_count[i] = calloc(bin[1], sizeof *bead_count);
    Temp[i] = calloc(bin[1], sizeof *Temp);
    sum_Temp[i] = calloc(bin[1], sizeof *sum_Temp);
    vir[i] = calloc(bin[1], sizeof *vir);
    sum_vir[i] = calloc(bin[1], sizeof *sum_vir);
    sum_beadcount[i] = calloc(bin[1], sizeof *sum_beadcount);
  }
  InitDouble2DArray(sum_Temp, bin[0], bin[1], 0);

  // TODO: file to save some global stuff
  FILE *t_out = OpenFile("t.txt", "w");
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
      // define width based on instantenous box size
      for (int j = 0; j < 2; j++) {
        if (dim[j+1] == 0) {
          width[j] = System.Box.Length.x / bin[j];
        } else if (dim[j+1] == 1) {
          width[j] = System.Box.Length.y / bin[j];
        } else if (dim[j+1] == 2) {
          width[j] = System.Box.Length.z / bin[j];
        }
      }
      WrapJoinCoordinates(&System, true, false); // restore pbc
      InitLong2DArray(bead_count, bin[0], bin[1], 0);
      InitDouble2DArray(Temp, bin[0], bin[1], 0);
      // calculate the local properties //{{{
      double temperature = 0, energy_kin = 0, virial = 0;
      for (int i = 0; i < System.Count.BeadCoor; i++) {
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
        if (dim[0] == 1) {
          if (dim[1] == 0) {
            centre.x = width[0] * (2 * pos[0] + 1) / 2;
          } else if (dim[1] == 1) {
            centre.y = width[0] * (2 * pos[0] + 1) / 2;
          } else if (dim[1] == 2) {
            centre.z = width[0] * (2 * pos[0] + 1) / 2;
          }
        } else { // dim[0] must be 2 (for 2D)
          if (dim[1] == 0) {
            centre.x = width[0] * (2 * pos[0] + 1) / 2;
            if (dim[2] == 1) {
              centre.y = width[1] * (2 * pos[1] + 1) / 2;
            } else {
              centre.z = width[1] * (2 * pos[1] + 1) / 2;
            }
          } else if (dim[1] == 1) {
            centre.y = width[0] * (2 * pos[0] + 1) / 2;
            if (dim[2] == 0) {
              centre.x = width[1] * (2 * pos[1] + 1) / 2;
            } else {
              centre.z = width[1] * (2 * pos[1] + 1) / 2;
            }
          } else { // dim[2] must be 2
            centre.z = width[0] * (2 * pos[0] + 1) / 2;
            if (dim[2] == 0) {
              centre.x = width[1] * (2 * pos[1] + 1) / 2;
            } else {
              centre.x = width[1] * (2 * pos[1] + 1) / 2;
            }
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
        double vel2 = SQR(b->Velocity.x / red_dist[b->Type]) +
                      SQR(b->Velocity.y / red_dist[b->Type]) +
                      SQR(b->Velocity.z / red_dist[b->Type]);
        Temp[pos[0]][pos[1]] += bt->Mass * vel2;
        vir[pos[0]][pos[1]] += b->Force.x * (b->Position.x - centre.x + 10) +
                               b->Force.y * (b->Position.y - centre.y + 11) +
                               b->Force.z * (b->Position.z - centre.z + 12);
        centre.x = System.Box.Length.x / 2;
        centre.y = System.Box.Length.y / 2;
        centre.z = System.Box.Length.z / 2;
        virial += b->Force.x * b->Position.x +
                  b->Force.y * b->Position.y +
                  b->Force.z * b->Position.z;
        temperature += bt->Mass * vel2;
        energy_kin += 0.5 * bt->Mass * vel2;
      }
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          Temp[i][j] *= c.mvv2e / (3 * c.boltz);
          sum_Temp[i][j] += Temp[i][j];
          sum_vir[i][j] += vir[i][j];
          sum_beadcount[i][j] += bead_count[i][j];
        }
      }
      temperature *= c.mvv2e / (3 * System.Count.BeadCoor * c.boltz);
      double pressure = (temperature * c.boltz * System.Count.BeadCoor +
                        virial) / (3 * System.Box.Volume) * c.nktv2p;
      energy_kin *= c.mvv2e; //}}}
      fprintf(t_out, " %8d", count_coor);
      fprintf(t_out, " %8.4f", temperature);
      fprintf(t_out, " %8.4f", energy_kin);
      fprintf(t_out, " %8.4f", virial);
      fprintf(t_out, " %8.4f", pressure);
      putc('\n', t_out);
      // int avg = 0; // use data from bins -avg to avg
      // TODO: add some averaging: sum up bins x-avg to x+avg, plot to x
      // save values (if using --pers-step) //{{{
      if (per_step) {
        char file[LINE];
        if (snprintf(file, LINE, "%s-%04d.txt", out_file, count_coor) < 0) {
          strcpy(ERROR_MSG, "something wrong with snprintf()");
          PrintErrorFile(file, "\0", "\0");
          exit(1);
        }
        FILE *fw = OpenFile(file, "w");
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
          fprintf(fw, "; (%d) W", ++count);
          fprintf(fw, "; (%d) P", ++count);
          fprintf(fw, "; (%d) avg beads in the bin", ++count);
          putc('\n', fw);
          for (int i = 0; i < bin[0]; i++) {
            double avg_T = 0, avg_vir = 0;
            int avg_count = 0;
            for (int j = (i - avg); j <= (i + avg); j++) {
              int k = j;
              if (k < 0) {
                k = bin[0] + k;
              } else if (k >= bin[0]) {
                k = k - bin[0];
              }
              avg_T += Temp[k][0];
              avg_vir += sum_vir[k][0];
              avg_count += bead_count[k][0];
            }
            int n = 2 * avg + 1;
            fprintf(fw, " %8.4f", width[0] * (2 * i + 1) / 2);
            fprintf(fw, " %8.4f", avg_T / avg_count);
            fprintf(fw, " %8.4f", avg_vir / avg_count);
            fprintf(fw, " %8.4f", (avg_T  + avg_vir / 3) / volume);
            fprintf(fw, " %8.4f", (double)(avg_count) / (count_used * n));
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
          fprintf(fw, "# (%d) %s", ++count, axis[0]);
          fprintf(fw, "; (%d) %s", ++count, axis[1]);
          fprintf(fw, "; (%d) T", ++count);
          fprintf(fw, "; (%d) W", ++count);
          fprintf(fw, "; (%d) P", ++count);
          fprintf(fw, "; (%d) avg beads in the bin", ++count);
          putc('\n', fw);
          for (int i = 0; i < bin[0]; i++) {
            for (int j = 0; j < bin[1]; j++) {
              double avg_T = 0, avg_vir = 0;
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
                  avg_T += Temp[m][n];
                  avg_vir += vir[m][n];
                  avg_count += bead_count[m][n];
                }
              }
              int n = SQR(2 * avg + 1);
              fprintf(fw, " %8.4f", width[0] * (2 * i + 1) / 2);
              fprintf(fw, " %8.4f", width[1] * (2 * j + 1) / 2);
              fprintf(fw, " %8.4f", avg_T / avg_count);
              fprintf(fw, " %8.4f", avg_vir / avg_count);
              fprintf(fw, " %8.4f", (avg_T + avg_vir / 3) / volume);
              fprintf(fw, " %8.4f\n", (double)(avg_count) / n);
              putc('\n', fw);
            }
            putc('\n', fw);
          }
        }
        //}}}
        fclose(fw);
      } //}}}
      // // calculate SC energy, etc.
      // // TODO: eventually add to the 'calculate the local properties' loop
      // long double *e_rep = calloc(System.Count.BeadCoor, sizeof *e_rep),
      //             *e_rho = calloc(System.Count.BeadCoor, sizeof *e_rho);
      // for (int i = 0; i < System.Count.BeadCoor; i++) {
      //   for (int j = 0; j < System.Count.BeadCoor; j++) {
      //     if (i != j) {
      //       int id1 = System.BeadCoor[i],
      //           id2 = System.BeadCoor[j];
      //       BEAD *b1 = &System.Bead[id1],
      //            *b2 = &System.Bead[id2];
      //       BEADTYPE *bt1 = &System.BeadType[b1->Type],
      //                *bt2 = &System.BeadType[b2->Type];
      //       VECTOR d = Distance(b1->Position, b2->Position, System.Box.Length);
      //       d.x = VectorLength(d);
      //       // TODO: check ml's code
      //     }
      //   }
      // }
      // free(e_rep);
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
  fclose(t_out);
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
    FILE *fw = OpenFile(out_file, "a");
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
      fprintf(fw, "; (%d) W", ++count);
      fprintf(fw, "; (%d) P", ++count);
      fprintf(fw, "; (%d) avg beads in the bin", ++count);
      putc('\n', fw);
      // int avg = 1;
      for (int i = 0; i < bin[0]; i++) {
        double avg_T = 0, avg_vir = 0;
        int avg_count = 0;
        for (int j = (i-avg); j <= (i+avg); j++) {
          int k = j;
          if (k < 0) {
            k = bin[0] + k;
          } else if (k >= bin[0]) {
            k = k - bin[0];
          }
          avg_T += sum_Temp[k][0];
          avg_vir += sum_vir[k][0];
          avg_count += sum_beadcount[k][0];
        }
        int n = 2 * avg + 1;
        fprintf(fw, " %8.4f", width[0] * (2 * i + 1) / 2);
        fprintf(fw, " %8.4f", avg_T / avg_count);
        fprintf(fw, " %8.4f", avg_vir / avg_count);
        fprintf(fw, " %8.4f", (avg_T + avg_vir / 3) / volume);
        fprintf(fw, " %8.4f", (double)(avg_count) / (count_used * n));
        putc('\n', fw);
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
      fprintf(fw, "# (%d) %s", ++count, axis[0]);
      fprintf(fw, "; (%d) %s", ++count, axis[1]);
      fprintf(fw, "; (%d) T", ++count);
      fprintf(fw, "; (%d) W", ++count);
      fprintf(fw, "; (%d) P", ++count);
      fprintf(fw, "; (%d) avg beads in the bin", ++count);
      putc('\n', fw);
      // int avg = 0;
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          double avg_T = 0, avg_vir = 0;
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
              avg_T += sum_Temp[m][n];
              avg_vir += sum_vir[m][n];
              avg_count += sum_beadcount[m][n];
            }
          }
          int n = SQR(2 * avg + 1);
          fprintf(fw, " %8.4f", width[0] * (2 * i + 1) / 2);
          fprintf(fw, " %8.4f", width[1] * (2 * j + 1) / 2);
          fprintf(fw, " %8.4f", avg_T / avg_count);
          fprintf(fw, " %8.4f", avg_vir / avg_count);
          fprintf(fw, " %8.4f", (avg_T + avg_vir / 3) / volume);
          fprintf(fw, " %8.4f\n", (double)(avg_count) / (count_used * n));
        }
        fprintf(fw, "\n");
      }
    } //}}}
    fclose(fw);
  } //}}}

  // free memory //{{{
  for (int i = 0; i < bin[0]; i++) {
    free(bead_count[i]);
    free(Temp[i]);
    free(sum_Temp[i]);
    free(vir[i]);
    free(sum_vir[i]);
    free(sum_beadcount[i]);
  }
  free(bead_count);
  free(Temp);
  free(sum_Temp);
  free(vir);
  free(sum_vir);
  free(sum_beadcount);
  FreeSystem(&System);
  //}}}

  return 0;
}
