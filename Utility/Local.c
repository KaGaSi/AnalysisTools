#include "../AnalysisTools.h"

typedef struct constant {
  // data from lammps/src/update.cpp
  double mvv2e, // convert mv^2 to E (eV)
         boltz, // Boltzmann constant in eV/K
         nktv2p;
} CONST;

// reduced length - CG: Al 3^3, Zr 2^3
// ...not sure whether this is in any sense right
static const double red_dist[2] = {1, 1};

void Help(char cmd[50], bool error) { //{{{
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
  int common = 12;
  char option[common][OPT_LENGTH];
  strcpy(option[ 0], "-st");
  strcpy(option[ 1], "-e");
  strcpy(option[ 2], "-sk");
  strcpy(option[ 3], "-i");
  strcpy(option[ 4], "--variable");
  strcpy(option[ 5], "-pbc");
  strcpy(option[ 6], "-ltrj");
  strcpy(option[ 7], "--detailed");
  strcpy(option[ 8], "-v");
  strcpy(option[ 9], "--silent");
  strcpy(option[10], "--help");
  strcpy(option[11], "--version");
  CommonHelp(error, common, option);
} //}}}

int main(int argc, char *argv[]) {
  // --help/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
    count++;
  }
  // possible to omit <bead name(s)> if '--reverse' is used
  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && strcmp(argv[i], "--metal") != 0 &&
        strcmp(argv[i], "--per_step") != 0 && strcmp(argv[i], "-avg") != 0 &&
        strcmp(argv[i], "-st") != 0 && strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 && strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "--variable") != 0 && strcmp(argv[i], "-pbc") != 0 &&
        strcmp(argv[i], "-ltrj") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-v") != 0 && strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "--help") != 0 && strcmp(argv[i], "--version") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

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
    Help(argv[0], true);
  } //}}}

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  double width[2];
  if (!IsPosRealNumber(argv[++count], &width[0])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  width[1] = width[0]; // TODO: for now
  //}}}

  // <output> - output text file
  char out_file[LINE] = "";
  snprintf(out_file, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var_coor;
  int start= 1, end = -1, skip = 0, pbc_xyz = -1, ltrj_start_id = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var_coor,
                &pbc_xyz, &ltrj_start_id, &start, &end, &skip);
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

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var_coor,
                                pbc_xyz, &ltrj_start_id);
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
      if (!ReadTimestep(coor_type, coor, coor_file, &System, &line_count,
                        ltrj_start_id, vtf_var_coor)) {
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
      // calculate the local properties
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
        virial += b->Force.x * (b->Position.x - centre.x + 10) +
                  b->Force.y * (b->Position.y - centre.y + 11) +
                  b->Force.z * (b->Position.z - centre.z + 12);
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
      energy_kin *= c.mvv2e;
      fprintf(t_out, " %8d", count_coor);
      fprintf(t_out, " %8.4f", temperature);
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
    strcpy(ERROR_MSG, "no coordinates written (starting timestep higher \
than the number of timestep)");
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
