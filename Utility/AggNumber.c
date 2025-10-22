#include "../AnalysisTools.h"

// print header for an avg output file (-a option)
static void PrintFileAvgHeader(int argc, char *argv[],
                               OPT *opt, SYSTEM System);
// print header for a distr output file (-d option)
static void PrintFileDistrHeader(int argc, char *argv[],
                                 OPT *opt, SYSTEM System);
// append overall averages to file(s) (-a/-d options)
static void AppendOverallAvg(char *f, SYSTEM System, double As_sum[3][2],
                             double mass_sum[3][2], int *count_agg_per_size,
                             int timesteps, ArrNDi *molecules_sum);
// helper functions for AppendOverallAvg()
static void PrintOverallAvgHeader(FILE *f, SYSTEM System);
static void PrintOverallAvg(FILE *fw, SYSTEM System, double As_sum[3][2],
                            double mass_sum[3][2], int sum_aggs,
                            int timesteps, int *molecules_sum);

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
AggNumber calculates time evolutions and distributions of aggregate sizes or \
average composition of aggregates. The definition of aggregate size (and the \
used range of sizes) is flexible. Besides distribution of sizes, it can also \
calculate composition distribution for specified aggregate size(s), i.e., \
distribution of numbers of various molecules \
in aggregates of given size(s).\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <in.agg> <distr file> <avg file> "
          "[options]\n\n", cmd);

  fprintf(ptr, "<in.stru>           input structure file\n");
  fprintf(ptr, "<in.agg>            input agg file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -n <size> <size>  use aggregate sizes in a given range\n");
  fprintf(ptr, "  -m <name(s)>      use number of specified molecule type(s) "
          "as aggrete size\n");
  fprintf(ptr, "  -x <name(s)>      exclude aggregates containing only "
          "specified molecule(s)\n");
  fprintf(ptr, "  -only <name(s)>   use only aggregates composed of "
          "specified molecule type(s)\n");
  fprintf(ptr, "  -d <distr>        output file with distributions\n");
  fprintf(ptr, "  -a <avg>          output file with per-timestep averages\n");
  fprintf(ptr, "  -c <file> <size(s)>\n");
  fprintf(ptr, "                    write composition distributions for "
          "aggregate size(s) to two <file>s with automatic endings '-#.txt' "
          "and '_r-#.txt'\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct comp { // for composition distribution
  int size[100], // aggregate sizes
      count;     // number of sizes
  char f[LINE];  // filename
};
struct OPT {
  AGG_PICKER agg;     // -x, -only, -m, and -n arrays
  struct comp comp;   // -c
  char f_distr[LINE], // -d
       f_avg[LINE];   // -a
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 7, all = common + 7, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "--verbose", "--silent", "--help",
               "--version", "-n", "-m", "-x", "-only", "-c", "-d", "-a");

  // commad line arguments before reading the structure //{{{
  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // <input> - input structure file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, argv[++count], LINE);
  in.stru.type = StructureFileType(in.stru.name);
  // <in.agg> - input aggregate file
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);

  // -d <distr> - file with distribution of aggregation numbers
  FileOption(argc, argv, "-d", opt->f_distr);
  // -a <avg> - file with per-timestep average aggregation numbers
  FileOption(argc, argv, "-a", opt->f_avg);
  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);
  // -c option
  FileNumbersOption(argc, argv, 1, 100, "-c", opt->comp.size,
                    &opt->comp.count, opt->comp.f, 'i');
  //}}}
  // Error - at least one of -d/-a/-c must be used
  if (opt->f_distr[0] == '\0' &&
      opt->f_avg[0] == '\0' &&
      opt->comp.f[0] == '\0') {
    err_msg("at least one must be specified (or no output would be generated)");
    PrintErrorOption("-d/-a/-c");
    exit(1);
  }

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  AggPickerOptions(argc, argv, &opt->agg, System);

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // arrays for distributions //{{{
  // number distribution
  long double *ndistr = calloc(Count->Molecule, sizeof *ndistr);
  /* weight and z distributions:
   *   [][0] = mass of mols according to options
   *   [][1] = mass of whole agg
   */
  ArrNDd *wdistr = CreateArr2Dd(Count->Molecule, 2);
  ArrNDd *zdistr = CreateArr2Dd(Count->Molecule, 2);
  // TODO: huh? what's molecules_sum?
  // molecule types in aggs: [agg size][mol type][number or Square(number)]
  ArrNDi *molecules_sum = CreateArr2Di(Count->Molecule, Count->MoleculeType);
  if (!wdistr || !zdistr || !molecules_sum) {
    err_msg("ArrNDd constructor failed (wdistr/zdistr/molecules_sum)");
    PrintError();
    exit(1);
  }
  // number of aggregates throughout simulation
  int *count_agg = calloc(Count->Molecule, sizeof *count_agg);
  // arrays for composition distribution
  ArrNDli *comp_distr = NULL; // [c_size][moltype][number of mols]
  ArrNDli *ratio_distr = NULL; // [c_size][moltype1][moltype2][num1][num2]
  long int *comp_agg_count = NULL;
  int *link_c_sizes = NULL;
  if (opt->comp.f[0] != '\0') {
    // array for 1D composition distribution
    if (!(comp_distr = CreateArr3Dli(opt->comp.count, Count->MoleculeType,
                                     Count->Molecule + 1))) {
      err_msg("ArrNDli constructor failed (comp_distr)");
      PrintError();
      exit(1);
    }

    // array for 2D composition distribution
    size_t shape_ratio_distr[5];
    shape_ratio_distr[0] = opt->comp.count;
    shape_ratio_distr[1] = Count->MoleculeType;
    shape_ratio_distr[2] = Count->MoleculeType;
    // +1 as it goes from no molecules to N molecules in the agg
    shape_ratio_distr[3] = Count->Molecule + 1;
    shape_ratio_distr[4] = Count->Molecule + 1;
    if (!(ratio_distr = CreateArrNDli(5, shape_ratio_distr))) {
      err_msg("ArrNDli constructor failed (ratio_distr)");
      PrintError();
      exit(1);
    }

    link_c_sizes = malloc(Count->Molecule * sizeof *link_c_sizes);
    InitIntArray(link_c_sizes, Count->Molecule, -1);
    comp_agg_count = calloc(opt->comp.count, sizeof *comp_agg_count);
    for (int i = 0; i < opt->comp.count; i++) {
      for (int j = 0; j < Count->Molecule; j++) {
        if (j == opt->comp.size[i]) {
          link_c_sizes[j] = i;
        }
      }
    }
  }
  // zeroize arrays
  FillArrND(wdistr, 0);
  FillArrND(zdistr, 0);
  for (int i = 0; i < Count->Molecule; i++) {
    ndistr[i] = 0;
    count_agg[i] = 0;
  } //}}}

  // open <in.agg> and skip the first two lines //{{{
  FILE *fr = OpenFile(input_agg, "r");
  while (getc(fr) != '\n')
    ;
  while (getc(fr) != '\n')
    ; //}}}

  // if -a is used, print the file header
  PrintFileAvgHeader(argc, argv, opt, System);

  // main loop //{{{
  int count_step = 0,
      count_used = 0,
      agg_lines = 2; // first two lines already read (skipped)
  /*
   * Mass and aggregate size sums
   *   [0][] = simple sum, [1][] = sum of squares, [2][] = sume of cubes
   *   [][0] = mass of mols in agg from options
   *   [][1] = mass of the whole aggregate
   */
  double mass_sum[3][2] = {{0}}, As_sum[3][2] = {{0}};
  while (true) { // cycle ends with 'Last Step' line in agg file
    PrintStep(&count_step, opt->c.start, opt->c.silent);

    // decide whether this timestep is to be used
    bool use = false;
    if (UseStep(opt->c, count_step)) {
      use = true;
    }
    if (use) { //{{{
      if (ReadAggregates(fr, input_agg, &System, Aggregate, &agg_lines) < 0) {
        count_step--;
        break;
      }
      count_used++; // just to print at the end
      int aggs_step = 0; // number of eligible aggregates per step
      double avg_mass_n_step[2] = {0}, // per-step mass averages
             avg_mass_w_step[2] = {0}, //  [0] ... from options
             avg_mass_z_step[2] = {0}, //  [1] ... for whole aggregates
             avg_As_n_step = 0,      // per-step As averages
             avg_As_w_step[2] = {0}, //  [0] ... from options
             avg_As_z_step[2] = {0}, //  [1] ... for whole aggregates
             molecules_step[Count->MoleculeType];
      // zeroize per-step counts of molecule types
      InitDoubleArray(molecules_step, Count->MoleculeType, 0);
      for (int i = 0; i < Count->Aggregate; i++) { //{{{
        // skip aggregates that shouldn't be used
        int agg_size;
        double agg_mass;
        if (!UseAggregate(System, Aggregate, i, opt->agg,
                          &agg_size, &agg_mass)) {
          continue;
        }
        // number of used aggregates in the step
        aggs_step++;
        // average aggregate mass during the step
        avg_mass_n_step[0] += agg_mass;
        avg_mass_w_step[0] += Square(agg_mass);
        avg_mass_z_step[0] += Cube(agg_mass); // unused
        avg_mass_n_step[1] += Aggregate[i].Mass;
        avg_mass_w_step[1] += Square(Aggregate[i].Mass);
        avg_mass_z_step[1] += Cube(Aggregate[i].Mass);
        // average aggregation number during the step
        avg_As_n_step += agg_size;
        avg_As_w_step[0] += agg_size * agg_mass;
        avg_As_z_step[0] += agg_size * Square(agg_mass);
        avg_As_w_step[1] += agg_size * Aggregate[i].Mass;
        avg_As_z_step[1] += agg_size * Square(Aggregate[i].Mass);
        // molecule species numbers
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mtype = System.Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_step[mtype]++;
        }

        // overall number of aggregates of given size
        count_agg[agg_size-1]++;
        // distributions
        ndistr[agg_size-1]++;
        AddArr2D(wdistr, agg_size - 1, 0, agg_mass);
        AddArr2D(zdistr, agg_size - 1, 0, Square(agg_mass));
        AddArr2D(wdistr, agg_size - 1, 1, Aggregate[i].Mass);
        AddArr2D(zdistr, agg_size - 1, 1, Square(Aggregate[i].Mass));
        // overall numbers of molecules of each species in each aggregate size
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
          AddArr2D(molecules_sum, agg_size - 1, mol_type, 1);
        }

        // composition distribution (-c option) //{{{
        if (opt->comp.f[0] != '\0' && link_c_sizes[agg_size] != -1) {
          comp_agg_count[link_c_sizes[agg_size]]++;
          int comp_aux[Count->MoleculeType];
          InitIntArray(comp_aux, Count->MoleculeType, 0);
          // count molecule types in the aggregate
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mtype = System.Molecule[Aggregate[i].Molecule[j]].Type;
            comp_aux[mtype]++;
          }
          // increment the distribution
          for (int j = 0; j < Count->MoleculeType; j++) {
            int id = link_c_sizes[agg_size];
            AddArr3D(comp_distr, id, j, comp_aux[j], 1);
            for (int k = (j + 1); k < Count->MoleculeType; k++) {
              size_t index5d[5] = {id, j, k, comp_aux[j], comp_aux[k]};
              AddArrND(ratio_distr, index5d, 1);
            }
          }
        } //}}}
      } //}}}

      // sums for overall averages //{{{
      // aggregate sizes
      As_sum[0][0] += avg_As_n_step;    // <full size>
      As_sum[1][0] += avg_As_w_step[0]; // <partial size> * <partial mass>
      As_sum[2][0] += avg_As_z_step[0]; // <partial size> * <partial mass>^2
      As_sum[1][1] += avg_As_w_step[1]; // <partial size> * <total mass>
      As_sum[2][1] += avg_As_z_step[1]; // <partial size> * <total mass>^2
      // aggregate masses
      mass_sum[0][0] += avg_mass_n_step[0]; // <partial mass>
      mass_sum[1][0] += avg_mass_w_step[0]; // <partial mass>^2
      mass_sum[2][0] += avg_mass_z_step[0]; // <partial mass>^3 - unused
      mass_sum[0][1] += avg_mass_n_step[1]; // <total mass>
      mass_sum[1][1] += avg_mass_w_step[1]; // <total mass>^2
      mass_sum[2][1] += avg_mass_z_step[1]; // <total mass>^3 //}}}

      // print averages to output file //{{{
      if (opt->f_avg[0] != '\0') {
        FILE *fw = OpenFile(opt->f_avg, "a");
        fprintf(fw, "%5d", count_step); // step
        if (aggs_step > 0) {
          fprintf(fw, " %10.5f", avg_As_n_step/aggs_step); // <As>_n
          fprintf(fw, " %10.5f", avg_As_w_step[0]/avg_mass_n_step[0]); // <As>_w
          fprintf(fw, " %10.5f", avg_As_w_step[1]/avg_mass_n_step[1]); //
          fprintf(fw, " %10.5f", avg_As_z_step[0]/avg_mass_w_step[0]); // <As>_z
          fprintf(fw, " %10.5f", avg_As_z_step[1]/avg_mass_w_step[1]); //
          fprintf(fw, " %10.5f", avg_mass_n_step[1]/aggs_step); // <m>_n
          fprintf(fw, " %10.5f", avg_mass_w_step[1]/avg_mass_n_step[1]); // <m>_w
          fprintf(fw, " %10.5f", avg_mass_z_step[1]/avg_mass_w_step[1]); // <m>_z
          for (int i = 0; i < Count->MoleculeType; i++) {
            fprintf(fw, " %10.5f", molecules_step[i]/aggs_step);
          }
        } else { // zero everywhere if there are no aggregates of the specified type
          fprintf(fw, " %10.5f", 0.0); // <mass>_n
          fprintf(fw, " %10.5f", 0.0); // <mass>_w
          fprintf(fw, " %10.5f", 0.0); // <mass>_z
          fprintf(fw, " %10.5f", 0.0); // <As>_n
          fprintf(fw, " %10.5f", 0.0); // <As>_w
          fprintf(fw, " %10.5f", 0.0); // <As>_z
          for (int i = 0; i < Count->MoleculeType; i++) {
            fprintf(fw, " %10.5f", 0.0);
          }
        }
        fprintf(fw, " %10d", aggs_step); // number of aggregates in the step
        // numbers of species
        putc('\n', fw);
        fclose(fw);
      } //}}}
      ReInitAggregate(System, Aggregate);
      //}}}
    } else {
      if (!SkipAggregates(fr, input_agg, &agg_lines)) {
        count_step--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_step == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  // print last step //{{{
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d", count_step);
    fprintf(stdout, " (%d used)\n",
            count_used);
  } //}}}
  //}}}

  // print distributions to output file //{{{
  if (opt->f_distr[0] != '\0') {
    PrintFileDistrHeader(argc, argv, opt, System);
    // normalization factors
    long int ndistr_norm = 0,
             wdistr_norm[2] = {0},
             zdistr_norm[2] = {0};
    for (int i = 0; i < Count->Molecule; i++) {
      ndistr_norm += ndistr[i];
      wdistr_norm[0] += GetArr2D(wdistr, i, 0);
      zdistr_norm[0] += GetArr2D(zdistr, i, 0);
      wdistr_norm[1] += GetArr2D(wdistr, i, 1);
      zdistr_norm[1] += GetArr2D(zdistr, i, 1);
    }
    // collate data //{{{
    int ncols = Count->MoleculeType + 7; // As, 5xF(As), n_agg
    int nrows = Count->Molecule;
    ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
    if (!data) {
      err_msg("ArrNDd constructor failed (data)");
      PrintError();
      exit(1);
    }
    for (int i = 0; i < nrows; i++) {
      if (count_agg[i] > 0) {
        count = 0;
        // agg size
        SetArr2D(data, i, count++, i + 1);
        // number distribution
        SetArr2D(data, i, count++, ndistr[i] / ndistr_norm);
        // weight distribution... selected mass, then real mass
        double val = (double)(GetArr2D(wdistr, i, 0)) / wdistr_norm[0];
        SetArr2D(data, i, count++, val);
        val = (double)(GetArr2D(wdistr, i, 1)) / wdistr_norm[1];
        SetArr2D(data, i, count++, val);
        // z distribution... selected mass, then real mass
        val = (double)(GetArr2D(zdistr, i, 0)) / zdistr_norm[0];
        SetArr2D(data, i, count++, val);
        val = (double)(GetArr2D(zdistr, i, 1)) / zdistr_norm[1];
        SetArr2D(data, i, count++, val);
        // number of aggregates
        SetArr2D(data, i, count++, (double)(count_agg[i]));
        // average number of molecules in aggregates
        for (int j = 0; j < Count->MoleculeType; j++) {
          double val = (double)(GetArr2D(molecules_sum, i, j));
          val /= count_agg[i];
          SetArr2D(data, i, count++, val);
        }
      }
    } //}}}
    ComputeColumnWidths(nrows, ncols, data, 6);
    FILE *fw = OpenFile(opt->f_distr, "a");
    for (int i = 0; i < nrows; i++) {
      if (count_agg[i] > 0) {
        PrintDataRow(fw, nrows, i, ncols, data);
      }
    }
    fclose(fw);
    FreeArrND(data);
  } //}}}

  AppendOverallAvg(opt->f_distr, System, As_sum, mass_sum,
                   count_agg, count_used, molecules_sum);
  AppendOverallAvg(opt->f_avg, System, As_sum, mass_sum,
                   count_agg, count_used, molecules_sum);

  // print composition distribution(s) (-c option) //{{{
  if (opt->comp.f[0] != '\0') {
    for (int i = 0; i < opt->comp.count; i++) {
      // print the distribution //{{{
      char file[LINE];
      if (snprintf(file, LINE, "%s-%03d.txt",
                   opt->comp.f, opt->comp.size[i]) < 0) {
        ErrorSnprintf();
      }
      FILE *fw = PrintBylineOpenFile(file, argc, argv);
      // print header
      fprintf(fw, "# total number of aggregates with size %d: %ld\n",
              opt->comp.size[i], comp_agg_count[i]);
      fprintf(fw, "# (1) number of molecules of given type;");
      fprintf(fw, " fraction of aggregates with that many molecules of type:");
      for (int j = 0; j < Count->MoleculeType; j++) {
        fprintf(fw, " (%d) %s", j + 2, System.MoleculeType[j].Name);
        if (j != (Count->MoleculeType - 1)) {
          putc(',', fw);
        }
      }
      putc('\n', fw);
      // print data
      if (comp_agg_count[i] > 0) {
        int ncols = Count->MoleculeType + 1;
        int nrows = opt->comp.size[i] + 1;
        ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
        if (!data) {
          err_msg("ArrNDd constructor failed (data)");
          PrintError();
          exit(1);
        }
        for (int j = 0; j < nrows; j++) {
          count = 0;
          // data.d[idx2d(data, j, count++)] = j;
          SetArr2D(data, j, count++, (double)j);
          for (int k = 0; k < Count->MoleculeType; k++) {
            double val = GetArr3D(comp_distr, i, k, j);
            val /= comp_agg_count[i];
            SetArr2D(data, j, count++, val);
          }
        }

        ComputeColumnWidths(nrows, ncols, data, 6);
        PrintDataAll(fw, nrows, ncols, data);
        FreeArrND(data);
      } else {
        snprintf(ERROR_MSG, LINE, "no aggregates with size %s%d%s found "
                 "(may also be due to -m/-x/-only/-n options)",
                 ErrYellow(), opt->comp.size[i], ErrCyan());
        PrintWarning();
      }
      fclose(fw); //}}}
      // print the ratios //{{{
      if (snprintf(file, LINE, "%s_r-%03d.txt",
                   opt->comp.f, opt->comp.size[i]) < 0) {
        ErrorSnprintf();
      }
      fw = PrintBylineOpenFile(file, argc, argv);
      // print header
      fprintf(fw, "# total number of aggregates with size %d: %ld\n",
              opt->comp.size[i], comp_agg_count[i]);
      fprintf(fw, "# (1-2) number of molecules:");
      count = 3;
      for (int j = 0; j < Count->MoleculeType; j++) {
        for (int k = (j + 1); k < Count->MoleculeType; k++) {
          fprintf(fw, " (%d) %s-%s", count++, System.MoleculeType[j].Name,
                                              System.MoleculeType[k].Name);
          if (j != (Count->MoleculeType - 2) ||
              k != (Count->MoleculeType - 1)) {
            putc(',', fw);
          }
        }
      }
      putc('\n', fw);
      // print data
      if (comp_agg_count[i] > 0) {
        int ncols = Count->MoleculeType * (Count->MoleculeType - 1) / 2 + 2;
        int nrows = Square(Count->Molecule + 1);
        // collate data //{{{
        ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
        if (!data) {
          err_msg("ArrNDd constructor failed (data)");
          PrintError();
          exit(1);
        }
        int count_lines = 0;
        for (int j = 0; j <= Count->Molecule; j++) {
          for (int k = 0; k <= Count->Molecule; k++) {
            bool use = false;
            for (int l = 0; l < Count->MoleculeType; l++) {
              for (int m = (l + 1); m < Count->MoleculeType; m++) {
                size_t index5d[5] = {i, l, m, j, k};
                if (GetArrND(ratio_distr, index5d) > 0) {
                  use = true;
                  break;
                }
              }
              if (use) {
                break;
              }
            }
            if (!use) {
              continue;
            }
            count = 0;
            SetArr2D(data, count_lines, count++, j);
            SetArr2D(data, count_lines, count++, k);
            for (int l = 0; l < Count->MoleculeType; l++) {
              for (int m = (l + 1); m < Count->MoleculeType; m++) {
                size_t index5d[5] = {i, l, m, j, k};
                double avg = (double)(GetArrND(ratio_distr, index5d));
                avg /= comp_agg_count[i];
                SetArr2D(data, count_lines, count++, avg);
              }
            }
            count_lines++;
          }
        } //}}}
        ComputeColumnWidths(nrows, ncols, data, 6);
        for (int j = 0; j < count_lines; j++) {
          for (int col = 0; col < ncols; col++) {
            if (col < 2 || fabs(GetArr2D(data, j, col)) > 1e-5) {
              PrintDataValue(fw, nrows, j , col, data);
            } else {
              fprintf(fw, " %*s", (int)GetArr2D(data, nrows, col), "?");
            }
          }
          putc('\n', fw);
          if (j < (count_lines - 1) &&
              GetArr2D(data, j, 0) != GetArr2D(data, j + 1, 0)) {
            putc('\n', fw);
          }
        }
        // free(precisions);
        FreeArrND(data);
      }
      fclose(fw); //}}}
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  FreeArrNDi(molecules_sum);
  FreeAggPicker(&opt->agg);
  if (opt->comp.f[0] != '\0') {
    FreeArrNDli(comp_distr);
    FreeArrNDli(ratio_distr);
    free(comp_agg_count);
    free(link_c_sizes);
  }
  free(ndistr);
  FreeArrNDd(wdistr);
  FreeArrNDd(zdistr);
  free(count_agg);
  free(opt); //}}}

  return 0;
}

// print header for an avg output file (-a option) //{{{
static void PrintFileAvgHeader(int argc, char *argv[],
                               OPT *opt, SYSTEM Sys) {
  if (opt->f_avg[0] == '\0') {
    return;
  }
  FILE *fw = PrintBylineOpenFile(opt->f_avg, argc, argv);
  int count = 1;
  fprintf(fw, "# Column: ");
  fprintf(fw, "(%d) step, ", count++);
  fprintf(fw, "(%d) <As>_n, ", count++);
  fprintf(fw, "(%d) <As>_w (partial mass), ", count++);
  fprintf(fw, "(%d) <As>_w (total mass), ", count++);
  fprintf(fw, "(%d) <As>_z (partial mass), ", count++);
  fprintf(fw, "(%d) <As>_z (total mass), ", count++);
  fprintf(fw, "(%d) <M>_n, ", count++);
  fprintf(fw, "(%d) <M>_w, ", count++);
  fprintf(fw, "(%d) <M>_z, ", count++);
  for (int i = 0; i < Sys.Count.MoleculeType; i++) {
    fprintf(fw, "(%d) <%s>_n, ", count++, Sys.MoleculeType[i].Name);
  }
  fprintf(fw, "(%d) n_agg", count++);
  putc('\n', fw);
  if (!opt->agg.m_flag) {
    fprintf(fw, "# Note: The -m option was not used; therefore, 'partial mass'"
            " includes all molecules, and "
            "<As>_w/z (partial mass) = <As>_w/z (total mass)\n");
  } else {
    fprintf(fw, "# Note: 'partial mass' includes molecules specifid by -m (");
    bool first = true;
    for (int i = 0; i < Sys.Count.MoleculeType; i++) {
      if (opt->agg.m[i]) {
        if (!first) {
          fprintf(fw, ", ");
        }
        fprintf(fw, "%s", Sys.MoleculeType[i].Name);
        first = false;
      }
    }
    fprintf(fw, ")\n");
    fprintf(fw, "#       'total mass' includes "
            "all molecules present in each aggregate\n");
  }
  fclose(fw);
} //}}}
// print header for a distr output file (-d option) //{{{
static void PrintFileDistrHeader(int argc, char *argv[],
                                 OPT *opt, SYSTEM System) {
  if (opt->f_distr[0] == '\0') {
    return;
  }
  FILE *fw = PrintBylineOpenFile(opt->f_distr, argc, argv);
  int count = 1;
  fprintf(fw, "# column: ");
  fprintf(fw, "(%d) As, ", count++);
  fprintf(fw, "(%d) F_n, ", count++);
  fprintf(fw, "(%d) F_w (partial mass), ", count++);
  fprintf(fw, "(%d) F_w (total mass), ", count++);
  fprintf(fw, "(%d) F_z (partial mass), ", count++);
  fprintf(fw, "(%d) F_z (total mass), ", count++);
  fprintf(fw, "(%d) n_agg,", count++);
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    fprintf(fw, " (%d) <%s>_n", i + count, System.MoleculeType[i].Name);
    if (i != (System.Count.MoleculeType-1)) {
      putc(',', fw);
    }
  }
  putc('\n', fw);
  if (!opt->agg.m_flag) {
    fprintf(fw, "# Note: The -m option was not used; therefore, 'partial mass'"
            " includes all molecules, and "
            "<As>_w/z (partial mass) = <As>_w/z (total mass)\n");
  } else {
    fprintf(fw, "# Note: 'partial mass' includes molecules specifid by -m (");
    bool first = true;
    for (int i = 0; i < System.Count.MoleculeType; i++) {
      if (opt->agg.m[i]) {
        if (!first) {
          fprintf(fw, ", ");
        }
        fprintf(fw, "%s", System.MoleculeType[i].Name);
        first = false;
      }
    }
    fprintf(fw, ")\n");
    fprintf(fw, "#       'total mass' includes "
            "all molecules present in each aggregate\n");
  }
  fclose(fw);
} //}}}
// append overall averages to file(s) (-a/-d options) //{{{
static void AppendOverallAvg(char *f, SYSTEM System, double As_sum[3][2],
                             double mass_sum[3][2], int *count_agg_per_size,
                             int timesteps, ArrNDi *molecules_sum) {
  if (f[0] == '\0') {
    return;
  }
  int *mol_sum_per_size = calloc(System.Count.MoleculeType,
                                 sizeof *mol_sum_per_size);
  int sum_aggs = 0;
  for (int i = 0; i < System.Count.Molecule; i++) {
    sum_aggs += count_agg_per_size[i];
    for (int j = 0; j < System.Count.MoleculeType; j++) {
      mol_sum_per_size[j] += GetArr2D(molecules_sum, i, j);
    }
  }
  FILE *fa = OpenFile(f, "a");
  PrintOverallAvgHeader(fa, System);
  PrintOverallAvg(fa, System, As_sum, mass_sum, sum_aggs,
                  timesteps, mol_sum_per_size);
  fclose(fa);
  free(mol_sum_per_size);
}
// print header
static void PrintOverallAvgHeader(FILE *f, SYSTEM System) {
  // distr file
  int count = 1;
  putc('#', f);
  fprintf(f, " (%d) <As>_n,", count++);
  fprintf(f, " (%d) <As>_w (partial mass),", count++);
  fprintf(f, " (%d) <As>_w (total mass),", count++);
  fprintf(f, " (%d) <As>_z (partial mass),", count++);
  fprintf(f, " (%d) <As>_z (total mass),", count++);
  fprintf(f, " (%d) <M>_n,", count++);
  fprintf(f, " (%d) <M>_w,", count++);
  fprintf(f, " (%d) <M>_z,", count++);
  for (int j = 0; j < System.Count.MoleculeType; j++) {
    fprintf(f, " (%d) <%s>_n,", count++, System.MoleculeType[j].Name);
  }
  fprintf(f, " (%d) <n_agg>", count++);
  putc('\n', f);
}
// print data
static void PrintOverallAvg(FILE *f, SYSTEM System, double As_sum[3][2],
                            double mass_sum[3][2], int sum_aggs,
                            int timesteps, int *molecules_sum) {
  fprintf(f, "#");
  if (sum_aggs > 0) {
    fprintf(f, " %lf", As_sum[0][0]/sum_aggs); // <As>_n
    fprintf(f, " %lf", As_sum[1][0]/mass_sum[0][0]); // <As>_w
    fprintf(f, " %lf", As_sum[1][1]/mass_sum[0][1]); //
    fprintf(f, " %lf", As_sum[2][0]/mass_sum[1][0]); // <As>_z
    fprintf(f, " %lf", As_sum[2][1]/mass_sum[1][1]); //

    fprintf(f, " %lf", mass_sum[0][1]/sum_aggs); // <M>_n
    fprintf(f, " %lf", mass_sum[1][1]/mass_sum[0][1]); // <M>_w
    fprintf(f, " %lf", mass_sum[2][1]/mass_sum[1][1]); // <M>_z
    for (int j = 0; j < System.Count.MoleculeType; j++) {
      // <species>_n
      double val = (double)(molecules_sum[j]) / sum_aggs;
      fprintf(f, " %lf", val);
    }
    fprintf(f, " %lf", (double)(sum_aggs)/timesteps); // <n_agg>
  } else { // zero everywhere if no aggregates found
    fprintf(f, " 0.0  0.0  0.0  0.0  0.0  0.0  0.0");
    for (int j = 0; j < System.Count.MoleculeType; j++) {
      fprintf(f, "  0.0");
    }
  }
  putc('\n', f);
} //}}}
