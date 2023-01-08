#include "../AnalysisTools.h"

typedef struct constant {
  // data from lammps/src/update.cpp
  double mvv2e, // convert mv^2 to E (eV)
         k_B; // Boltzmann constant in eV/K
} CONST;

// TODO: -e X --last leads to leaving at step X and saving that as the last;
//       is that okay?
// TODO: make some in_struct akin to in_coor (in InputCoorStruct() function)

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
  fprintf(ptr, "   <dimensions>      1 or 2 axis for 1D or 2D calculation"
          " (e.g., x or xy)\n");
  fprintf(ptr, "   <width>           width of a single bin"
          " (in the 1 or both dimensions)\n");
  fprintf(ptr, "   <output>          output file with local variables\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -sk <skip>     leave out every 'skip' steps\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps \
(if --last option is used, save also the last timestep)\n");
  fprintf(ptr, "      --last         use only the last step \
(-st/-e options are ignored; -n option is not)\n");
  fprintf(ptr, "      --metal        assume lammps 'metal' units\n");
  fprintf(ptr, "      --per-step     save a new file for each step"
          "(adds '-<ste>.txt' to <output>)\n");
  fprintf(ptr, "      -avg <n>       number of adjacent bins to average\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {
  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
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
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "--silent") != 0 && strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 && strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 && strcmp(argv[i], "-sk") != 0 &&
        strcmp(argv[i], "-n") != 0 && strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "--last") != 0 && strcmp(argv[i], "--metal") != 0 &&
        strcmp(argv[i], "--per-step") != 0 && strcmp(argv[i], "-avg") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char in_coor[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(in_coor, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, in_coor, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // <axis> //{{{
  count++;
  bool err = false;
  int dim[3] = {-1, -1, -1}; // [0]..1/2 for 1D/2D, [0..1]..0, 1, 2 for x, y, z
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
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed);
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
    exit(1);
  }
  skip++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  // should output coordinates be joined?
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // use only the last step?
  bool last = BoolOption(argc, argv, "--last");
  int style = 0; // reduced units
  if (BoolOption(argc, argv, "--metal")) {
    style = 1; // lammps 'metal' units
  }
  bool per_step = BoolOption(argc, argv, "--per-step");
  int avg = 0;
  if (IntegerOption(argc, argv, "-avg", &avg)) {
    exit(1);
  }
  //}}}

  CONST c; //{{{
  if (style == 1) { // lammps metal units
    c.mvv2e = 1.0364269e-4;
    c.k_B = 8.617343e-5;
  } else { // reduced units (lammps lj units)
    c.mvv2e = 1;
    c.k_B = 1;
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // TODO: make reading various structure types into a function in Read.c
  SYSTEM System = ReadStructure(struct_type, struct_file, detailed);

  // read pbc - must be known to calculate number of bins
  if (coor_type == VCF_FILE) {
    System.Box = VtfReadPBC(in_coor);
  } else if (coor_type == LTRJ_FILE) {
    System.Box = LtrjReadPBC(in_coor);
  }

  // number of bins //{{{
  // TODO: variable box size? For now, assume at most twice the size
  //       nah, just make the bins variable width based on instantenous box size
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    // PrintErrorFile(); // TODO add file(s)
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

  // '-n' option - specify timestep ids //{{{
  int n_opt_save[100] = {0}, n_opt_number = -1;
  if (MultiIntegerOption(argc, argv, "-n", &n_opt_number, n_opt_save)) {
    exit(1);
  }
  SortArray(n_opt_save, n_opt_number, 0); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  // warning - unspecified mass //{{{
  for (int i = 0; i < System.Count.BeadType; i++) {
    if (System.BeadType[i].Mass == MASS) {
      strcpy(ERROR_MSG, "unspecified mass");
      PrintWarnFile(struct_file, in_coor, "\0");
      fprintf(stderr, "%s, bead type %s%s%s\n", ErrCyan(), ErrYellow(),
              System.BeadType[i].Name, ErrColourReset());
    }
  } //}}}

  // print initial stuff to output vcf file (if not -per-step) //{{{
  if (!per_step) {
    FILE *out = OpenFile(out_file, "w");
    PrintByline(out, argc, argv);
    fclose(out);
  } //}}}

  // main loop //{{{
  FILE *coor = OpenFile(in_coor, "r");
  fpos_t position1, position2; // two file pointers for finding the last step
  int n_opt_count = 0,         // count saved steps if -n option is used
      count_coor = 0,          // count steps in the vcf file
      count_used = 0,         // count steps in output file
      file_line_count = 0;     // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  // when 1D is calculated, only Temp[i][0] is used

  long **bead_count = calloc(bin[0], sizeof **bead_count); // count beads
  double **Temp = calloc(bin[0], sizeof **Temp); // temperature
  double **sum_Temp = calloc(bin[0], sizeof **sum_Temp);
  double **sum_beadcount = calloc(bin[0], sizeof **sum_beadcount);
  for (int i = 0; i < bin[0]; i++) {
    bead_count[i] = calloc(bin[1], sizeof *bead_count);
    Temp[i] = calloc(bin[1], sizeof *Temp);
    sum_Temp[i] = calloc(bin[1], sizeof *sum_Temp);
    sum_beadcount[i] = calloc(bin[1], sizeof *sum_beadcount);
  }
  InitDouble2DArray(sum_Temp, bin[0], bin[1], 0);

  // TODO: file to save some global stuff
  FILE *t_out = OpenFile("t.txt", "w");
  while (true) {
    count_coor++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (last) {
        fprintf(stdout, "\rDiscarding step: %d", count_coor);
      } else {
        if (count_coor == start) {
          fprintf(stdout, "\rStarting step: %d\n", start);
        }
        fprintf(stdout, "\rStep: %d", count_coor);
      }
      fflush(stdout);
    } //}}}
    // decide whether this timestep is to be used //{{{
    bool use = false;
    /* no -n option - use if timestep
     *    1) is between start (-st option) and end (-e option)
     *    and
     *    2) isn't skipped (-sk option); skipping starts counting with 'start'
     */
    if (n_opt_number == -1) {
      if ((count_coor >= start && (count_coor <= end || end == -1)) && // 1)
          ((count_coor - start) % skip) == 0) {                        // 2)
        use = true;
      } else {
        use = false;
      }
      // definitely not use, if --last option is used
      if (last) {
        use = false;
      }
      // -n option is used - save the timestep if it's in the list
    } else if (n_opt_count < n_opt_number &&
               n_opt_save[n_opt_count] == count_coor) {
      use = true;
      n_opt_count++;
    }
    //}}}
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (!ReadTimestep(coor_type, coor, in_coor, &System, &file_line_count,
                        stuff)) {
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
      // calculate the local properties
      InitLong2DArray(bead_count, bin[0], bin[1], 0);
      InitDouble2DArray(Temp, bin[0], bin[1], 0);
      double temperature = 0, energy_kin = 0;
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
        bead_count[pos[0]][pos[1]]++;
        double vel2 = SQR(b->Velocity.x) +
                      SQR(b->Velocity.y) +
                      SQR(b->Velocity.z);
        Temp[pos[0]][pos[1]] += bt->Mass * vel2;
        temperature += bt->Mass * vel2;
        energy_kin += 0.5 * bt->Mass * vel2;
      }
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          Temp[i][j] *= c.mvv2e / (3 * c.k_B);
          sum_Temp[i][j] += Temp[i][j];
          sum_beadcount[i][j] += bead_count[i][j];
        }
      }
      temperature *= c.mvv2e / (3 * System.Count.BeadCoor * c.k_B);
      energy_kin *= c.mvv2e;
      fprintf(t_out, "%8d %8.4f %8.4f\n", count_coor, temperature, energy_kin);
      // int avg = 0; // use data from bins -avg to avg
      // TODO: add some averaging: sum up bins x-avg to x+avg, plot to x
      // save values (if using --pers-step) //{{{
      if (per_step) {
        char file[LINE];
        snprintf(file, LINE, "%s-%04d.txt", out_file, count_coor);
        FILE *out = OpenFile(file, "w");
        // 1D case //{{{
        if (dim[0] == 1) {
          char *axis;
          if (dim[1] == 0) {
            axis = "x";
          } else if (dim[1] == 1) {
            axis = "y";
          } else {
            axis = "z";
          }
          count = 0;
          fprintf(out, "# (%d) %s coordinate", ++count, axis);
          fprintf(out, "; (%d) temperature", ++count);
          fprintf(out, "; (%d) avg beads in the bin", ++count);
          putc('\n', out);
          for (int i = 0; i < bin[0]; i++) {
            double avg_T = 0;
            int avg_count = 0;
            for (int j = (i-avg); j <= (i+avg); j++) {
              int k = j;
              if (k < 0) {
                k = bin[0] + k;
              } else if (k >= bin[0]) {
                k = k - bin[0];
              }
              avg_T += Temp[k][0];
              avg_count += bead_count[k][0];
            }
            int n = 2 * avg + 1;
            fprintf(out, "%8.4f %8.4f %8.4f\n", width[0]*(2*i+1)/2,
                    avg_T/avg_count, (double)(avg_count)/(count_used*n));
          }
        //}}}
        // 2D case //{{{
        } else {
          char *axis[2];
          if (dim[1] == 0) {
            axis[0] = "x";
          } else if (dim[1] == 1) {
            axis[0] = "y";
          } else {
            axis[0] = "z";
          }
          if (dim[2] == 0) {
            axis[1] = "x";
          } else if (dim[2] == 1) {
            axis[1] = "y";
          } else {
            axis[1] = "z";
          }
          count = 0;
          fprintf(out, "# (%d) %s coordinate", ++count, axis[0]);
          fprintf(out, "; (%d) %s coordinate", ++count, axis[1]);
          fprintf(out, "; (%d) temperature", ++count);
          fprintf(out, "; (%d) avg beads in the bin", ++count);
          putc('\n', out);
          for (int i = 0; i < bin[0]; i++) {
            for (int j = 0; j < bin[1]; j++) {
              double avg_T = 0;
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
                  avg_count += bead_count[m][n];
                }
              }
              int n = SQR(2 * avg + 1);
              fprintf(out, "%8.4f", width[0]*(2*i+1)/2);
              fprintf(out, " %8.4f", width[1]*(2*j+1)/2);
              fprintf(out, " %8.4f", avg_T/avg_count);
              fprintf(out, " %8.4f\n", (double)(avg_count)/n);
            }
            fprintf(out, "\n");
          }
        }
        //}}}
        fclose(out);
      } //}}}
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(coor_type, coor, in_coor,
                        struct_file, &file_line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // save file position (last two because of --last) //{{{
    if ((count_coor % 2) == 0) {
      fgetpos(coor, &position1);
    } else {
      fgetpos(coor, &position2);
    } //}}}
    // decide whether to exit the main loop //{{{
    /* break the loop if
     *    1) all timesteps in the -n option are saved (and --last isn't used)
     *    or
     *    2) end timestep was reached (-e option)
     */
    if ((n_opt_count == n_opt_number && !last) || // 1)
        count_coor == end) {                      // 2)
      break;
    } //}}}
  }   //}}}
  fclose(t_out);
  // if --last option is used, read & save the last timestep //{{{
  if (last) {
    // restore file pointer
    if ((count_coor % 2) == 1) {
      fsetpos(coor, &position1);
    } else {
      fsetpos(coor, &position2);
    }
    ReadTimestep(coor_type, coor, in_coor, &System, &file_line_count, stuff);
    count_used++;
    WrapJoinCoordinates(&System, true, false); // restore pbc
    // calculate the local properties
    InitLong2DArray(bead_count, bin[0], bin[1], 0);
    InitDouble2DArray(Temp, bin[0], bin[1], 0);
    double temperature = 0, energy_kin = 0;
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
      bead_count[pos[0]][pos[1]]++;
      double vel2 = SQR(b->Velocity.x) +
                    SQR(b->Velocity.y) +
                    SQR(b->Velocity.z);
      Temp[pos[0]][pos[1]] += bt->Mass * vel2;
      temperature += bt->Mass * vel2;
      energy_kin += 0.5 * bt->Mass * vel2;
    }
    for (int i = 0; i < bin[0]; i++) {
      for (int j = 0; j < bin[1]; j++) {
        Temp[i][j] *= c.mvv2e / (3 * c.k_B);
        sum_Temp[i][j] += Temp[i][j];
        sum_beadcount[i][j] += bead_count[i][j];
      }
    }
  } //}}}
  fclose(coor);
  if (count_coor == 0) { // error - input file without a valid timestep //{{{
    strcpy(ERROR_MSG, "no valid timestep found");
    PrintError();
    ErrorPrintFile(in_coor, "\0", "\0");
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
    FILE *out = OpenFile(out_file, "a");
    // 1D case //{{{
    if (dim[0] == 1) {
      char *axis;
      if (dim[1] == 0) {
        axis = "x";
      } else if (dim[1] == 1) {
        axis = "y";
      } else {
        axis = "z";
      }
      count = 0;
      fprintf(out, "# (%d) %s coordinate", ++count, axis);
      fprintf(out, "; (%d) temperature", ++count);
      fprintf(out, "; (%d) avg beads in the bin", ++count);
      putc('\n', out);
      // int avg = 1;
      for (int i = 0; i < bin[0]; i++) {
        double avg_T = 0;
        int avg_count = 0;
        for (int j = (i-avg); j <= (i+avg); j++) {
          int k = j;
          if (k < 0) {
            k = bin[0] + k;
          } else if (k >= bin[0]) {
            k = k - bin[0];
          }
          avg_T += sum_Temp[k][0];
          avg_count += sum_beadcount[k][0];
        }
        int n = 2 * avg + 1;
        fprintf(out, "%8.4f %8.4f %8.4f\n", width[0]*(2*i+1)/2,
                avg_T/avg_count, (double)(avg_count)/(count_used*n));
      }
    //}}}
    // 2D case //{{{
    } else {
      char *axis[2];
      if (dim[1] == 0) {
        axis[0] = "x";
      } else if (dim[1] == 1) {
        axis[0] = "y";
      } else {
        axis[0] = "z";
      }
      if (dim[2] == 0) {
        axis[1] = "x";
      } else if (dim[2] == 1) {
        axis[1] = "y";
      } else {
        axis[1] = "z";
      }
      count = 0;
      fprintf(out, "# (%d) %s coordinate", ++count, axis[0]);
      fprintf(out, "; (%d) %s coordinate", ++count, axis[1]);
      fprintf(out, "; (%d) temperature", ++count);
      fprintf(out, "; (%d) avg beads in the bin", ++count);
      putc('\n', out);
      // int avg = 0;
      for (int i = 0; i < bin[0]; i++) {
        for (int j = 0; j < bin[1]; j++) {
          double avg_T = 0;
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
              avg_T += sum_Temp[m][n];
              avg_count += sum_beadcount[m][n];
            }
          }
          int n = SQR(2 * avg + 1);
          fprintf(out, "%8.4f", width[0]*(2*i+1)/2);
          fprintf(out, " %8.4f", width[1]*(2*j+1)/2);
          fprintf(out, " %8.4f", avg_T/avg_count);
          fprintf(out, " %8.4f\n", (double)(avg_count)/(count_used*n));
        }
        fprintf(out, "\n");
      }
    } //}}}
    fclose(out);
  } //}}}

  // free memory //{{{
  for (int i = 0; i < bin[0]; i++) {
    free(bead_count[i]);
    free(Temp[i]);
    free(sum_beadcount[i]);
    free(sum_Temp[i]);
  }
  free(bead_count);
  free(Temp);
  free(sum_beadcount);
  free(sum_Temp);
  FreeSystem(&System);
  free(stuff);
  //}}}

  return 0;
}
