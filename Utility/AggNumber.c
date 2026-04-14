#include "../src/AnalysisTools.h"

//TODO: warning if 0 mass leading to nan in output file

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AggNumber calculates time evolutions and distributions of aggregate sizes or"
  "average composition of aggregates. The definition of aggregate size (and the"
  "used range of sizes) is flexible. Besides distribution of sizes, it can also"
  "calculate composition distribution for specified aggregate size(s), i.e.,"
  "distribution of numbers of various molecules"
  "in aggregates of given size(s).",

  "Usage: AggNumber <in.stru> <in.agg> [options]",
  .args = 2, // number of mandatory arguments
  .all = 16, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_VERSION],
  {"<in.stru>", NULL, "input structure file", OPT_ARG},
  {"<in.agg>", NULL, "input agg file", OPT_ARG},
  {"-a", "<file>", "calculate per-timestep averages", OPT_EXTRA},
  {"-d", "<file>", "calculate distributions of aggregate sizes", OPT_EXTRA},
  {"-c", "<file> <size(s)>", "composition distributions for given size(s) (two <file>s with automatic endings '-#.txt' and '_r-#.txt')", OPT_EXTRA},
  {"-n", "<size> <size>", "use aggregate sizes in a given range", OPT_EXTRA},
  {"-m", "<name(s)>", "use number of specified molecule type(s) as aggrete size", OPT_EXTRA},
  {"-x", "<name(s)>", "exclude aggregates containing only specified molecule(s)", OPT_EXTRA},
  {"-only", "<name(s)>", "use only aggregates composed of specified molecule type(s)", OPT_EXTRA},
  {NULL},
}; //}}}

// helper functions //{{{
// print header for an avg output file (-a option)
static void PrintAvgHeader(int argc, char *argv[], OPT opt, SYSTEM System);
// print header for a distr output file (-d option)
static void PrintDistrHeader(int argc, char *argv[], OPT opt, SYSTEM System);
// print header for a copmposition output file (-c option)
static void PrintCompSimpleHeader(int argc, char *argv[], OPT opt,
                                  SYSTEM System, int size,
                                  long int *comp_agg_count);
// print header for a 2D copmposition output file (-c option)
static void PrintComp2DHeader(int argc, char *argv[], OPT opt, SYSTEM System,
                              int size, long int *comp_agg_count);
// print the note about the two different mass definition
static void PrintHeaderMassNote(FILE *fw, SYSTEM System, OPT opt);
// append overall averages to file(s) (-a/-d options)
static void AppendOverallAvg(char *f, SYSTEM System, double As_sum[3][2],
                             double mass_sum[3][2], int *count_agg_per_size,
                             int timesteps, ArrNDi *molecules_sum);
// helper functions for AppendOverallAvg() (-a/-d options)
static void PrintOverallAvgHeader(FILE *f, SYSTEM System);
static void PrintOverallAvg(FILE *fw, SYSTEM System, double As_sum[3][2],
                            double mass_sum[3][2], int sum_aggs,
                            int timesteps, int *molecules_sum);
// create As-based filenames (-c option)
static void FilenameCompSimple(OPT opt, int size, char filename[LINE]);
static void FilenameCode2D(OPT opt, int size, char filename[LINE]); //}}}

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
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input structure file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, argv[++count], LINE);
  in.stru.type = StructureFileType(in.stru.name);
  // <in.agg> - input aggregate file
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // -a <avg> - file with per-timestep average aggregation numbers
  FileOption(argc, argv, "-a", opt.f_avg);
  // -d <distr> - file with distribution of aggregation numbers
  FileOption(argc, argv, "-d", opt.f_distr);
  // -c option
  FileNumbersOption(argc, argv, 1, 100, "-c", opt.comp.size,
                    &opt.comp.count, opt.comp.f, 'i');
  // Error - at least one of -a/-d/-c must be used
  if (opt.f_avg[0] == '\0' &&
      opt.f_distr[0] == '\0' &&
      opt.comp.f[0] == '\0') {
    err_msg("at least one must be specified (or no output would be generated)");
    PrintErrorOption("-a/-d/-c");
    exit(1);
  }
  //}}}

  // print command to stdout
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  AggPickerOptions(argc, argv, &opt.agg, System);

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // arrays for distribution (-d option) //{{{
  // number distribution - useful whatever is calculated
  long double *ndistr = calloc(Count->Molecule, sizeof *ndistr);
  ArrNDd *wdistr = NULL; // only for -d option
  ArrNDd *zdistr = NULL; //
  // molecule types in aggs: [agg size][mol type]; needed for overall averages
  ArrNDi *molecules_sum = CreateArr2Di(Count->Molecule, Count->MoleculeType);
  // number of aggregates throughout simulation - always useful
  int *count_agg = calloc(Count->Molecule, sizeof *count_agg);
  if (!count_agg || !ndistr || !molecules_sum) {
    ErrorAlloc("ndstr/count_agg/molecules_sum");
  }
  if (opt.f_distr[0] != '\0') {
    /*
     * weight and z distributions:
     *   [][0] = mass of mols according to options
     *   [][1] = mass of whole agg
     */
    if (!(wdistr = CreateArr2Dd(Count->Molecule, 2)) ||
        !(zdistr = CreateArr2Dd(Count->Molecule, 2))) {
      ErrorAlloc("wdistr/zdistr");
    }
  }
  //}}}

  // arrays for composition distribution (-c option) //{{{
  ArrNDli *comp_distr = NULL; // [c_size][moltype][number of mols]
  ArrNDli *ratio_distr = NULL; // [c_size][moltype1][moltype2][num1][num2]
  long int *comp_agg_count = NULL;
  int *link_c_sizes = NULL;
  if (opt.comp.f[0] != '\0') {
    // array for 1D composition distribution
    if (!(comp_distr = CreateArr3Dli(opt.comp.count, Count->MoleculeType,
                                     Count->Molecule + 1))) {
      ErrorAlloc("comp_distr");
    }

    // array for 2D composition distribution
    size_t shape_ratio_distr[5];
    shape_ratio_distr[0] = opt.comp.count;
    shape_ratio_distr[1] = Count->MoleculeType;
    shape_ratio_distr[2] = Count->MoleculeType;
    // +1 as it goes from no molecules to N molecules in the agg
    shape_ratio_distr[3] = Count->Molecule + 1;
    shape_ratio_distr[4] = Count->Molecule + 1;
    if (!(ratio_distr = CreateArrNDli(5, shape_ratio_distr))) {
      ErrorAlloc("ratio_distr");
    }

    link_c_sizes = malloc(Count->Molecule * sizeof *link_c_sizes);
    comp_agg_count = calloc(opt.comp.count, sizeof *comp_agg_count);
    if (!link_c_sizes || !comp_agg_count) {
      ErrorAlloc("link_c_sizes/comp_agg_count");
    }
    InitIntArray(link_c_sizes, Count->Molecule, -1);
    for (int i = 0; i < opt.comp.count; i++) {
      for (int j = 0; j < Count->Molecule; j++) {
        if (j == opt.comp.size[i]) {
          link_c_sizes[j] = i;
        }
      }
    }
  }
  // zeroize arrays
  if (opt.f_distr[0] != '\0') {
    FillArrND(wdistr, 0);
    FillArrND(zdistr, 0);
  }
  for (int i = 0; i < Count->Molecule; i++) {
    ndistr[i] = 0;
    count_agg[i] = 0;
  } //}}}

  // open <in.agg> and skip the first two lines //{{{
  FILE *fr = OpenFile(input_agg, "r");
  SkipLine(fr);
  SkipLine(fr);
  //}}}

  // if -a and/or -d is used, print the file header(s)
  PrintAvgHeader(argc, argv, opt, System);
  PrintDistrHeader(argc, argv, opt, System);

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
    PrintStep(&count_step, commons.start, commons.silent);

    // decide whether this timestep is to be used
    bool use = false;
    if (UseStep(commons, count_step)) {
      use = true;
    }
    if (use) { //{{{
      if (ReadAggregates(fr, input_agg, &System, Aggregate, &agg_lines) < 0) {
        count_step--;
        break;
      }
      count_used++; // just to print at the end
      int aggs_step = 0; // number of eligible aggregates per step
      double avg_mass_n_step[2] = {0, 0}, // per-step mass averages
             avg_mass_w_step[2] = {0, 0}, //  [0] ... from options
             avg_mass_z_step[2] = {0, 0}, //  [1] ... for whole aggregates
             avg_As_n_step = 0,      // per-step As averages
             avg_As_w_step[2] = {0, 0}, //  [0] ... from options
             avg_As_z_step[2] = {0, 0}, //  [1] ... for whole aggregates
             molecules_step[Count->MoleculeType];
      // zeroize per-step counts of molecule types
      InitDoubleArray(molecules_step, Count->MoleculeType, 0);
      for (int i = 0; i < Count->Aggregate; i++) { //{{{
        // skip aggregates that shouldn't be used
        int agg_size;
        double agg_mass;
        if (!UseAggregate(System, Aggregate, i, opt.agg,
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
          int mtype = System.Molecule[AggGetMol(&Aggregate[i], j)].Type;
          molecules_step[mtype]++;
        }

        // overall number of aggregates of given size
        count_agg[agg_size-1]++;
        // distributions
        ndistr[agg_size-1]++;
        // overall numbers of molecules of each species in each aggregate size
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = System.Molecule[AggGetMol(&Aggregate[i], j)].Type;
          AddArr2D(molecules_sum, agg_size - 1, mol_type, 1);
        }
        if (opt.f_distr[0] != '\0') {
          AddArr2D(wdistr, agg_size - 1, 0, agg_mass);
          AddArr2D(zdistr, agg_size - 1, 0, Square(agg_mass));
          AddArr2D(wdistr, agg_size - 1, 1, Aggregate[i].Mass);
          AddArr2D(zdistr, agg_size - 1, 1, Square(Aggregate[i].Mass));
        }

        // composition distribution (-c option) //{{{
        if (opt.comp.f[0] != '\0' && link_c_sizes[agg_size] != -1) {
          comp_agg_count[link_c_sizes[agg_size]]++;
          int comp_aux[Count->MoleculeType];
          InitIntArray(comp_aux, Count->MoleculeType, 0);
          // count molecule types in the aggregate
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mtype = System.Molecule[AggGetMol(&Aggregate[i], j)].Type;
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
      if (opt.f_avg[0] != '\0') {
        FILE *fw = OpenFile(opt.f_avg, "a");
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
    if (count_step == commons.end) {
      break;
    }
  }
  fclose(fr);
  // print last step //{{{
  if (!commons.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d", count_step);
    fprintf(stdout, " (%d used)\n", count_used);
  } //}}}
  //}}}

  // print distributions to output file //{{{
  if (opt.f_distr[0] != '\0') {
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
      ErrorAlloc("data");
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
    FILE *fw = OpenFile(opt.f_distr, "a");
    for (int i = 0; i < nrows; i++) {
      if (count_agg[i] > 0) {
        PrintDataRow(fw, nrows, i, ncols, data);
      }
    }
    fclose(fw);
    FreeArrND(data);
    // append overall averages
    AppendOverallAvg(opt.f_distr, System, As_sum, mass_sum,
                     count_agg, count_used, molecules_sum);
  } //}}}

  // append overall averages (for -a option)
  if (opt.f_avg[0] != '\0') {
    AppendOverallAvg(opt.f_avg, System, As_sum, mass_sum,
                     count_agg, count_used, molecules_sum);
  }

  // print composition distribution(s) (-c option) //{{{
  if (opt.comp.f[0] != '\0') {
    for (int i = 0; i < opt.comp.count; i++) {
      // print the distribution //{{{
      PrintCompSimpleHeader(argc, argv, opt, System, i, comp_agg_count);
      char file[LINE];
      FilenameCompSimple(opt, i, file);
      FILE *fw = OpenFile(file, "a");
      // print data
      if (comp_agg_count[i] > 0) {
        int ncols = Count->MoleculeType + 1;
        int nrows = opt.comp.size[i] + 1;
        ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
        if (!data) {
          ErrorAlloc("data");
        }
        for (int j = 0; j < nrows; j++) {
          count = 0;
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
                 ErrYellow(), opt.comp.size[i], ErrCyan());
        PrintWarning();
      }
      fclose(fw); //}}}
      // print the ratios //{{{
      FilenameCode2D(opt, i, file);
      PrintComp2DHeader(argc, argv, opt, System, i, comp_agg_count);
      fw = OpenFile(file, "a");
      // print data
      if (comp_agg_count[i] > 0) {
        int ncols = Count->MoleculeType * (Count->MoleculeType - 1) / 2 + 2;
        int nrows = Square(Count->Molecule + 1);
        // collate data //{{{
        ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
        if (!data) {
          ErrorAlloc("data");
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

  // free memory //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  FreeAggPicker(&opt.agg);
  if (opt.comp.f[0] != '\0') {
    FreeArrNDli(comp_distr);
    FreeArrNDli(ratio_distr);
    free(comp_agg_count);
    free(link_c_sizes);
  }
  free(ndistr);
  if (opt.f_distr[0] != '\0') {
    FreeArrNDd(wdistr);
    FreeArrNDd(zdistr);
    FreeArrNDi(molecules_sum);
  }
  free(count_agg);
  //}}}

  return 0;
}

// helper functions
// print header for an avg output file (-a option) //{{{
static void PrintAvgHeader(int argc, char *argv[], OPT opt, SYSTEM System) {
  if (opt.f_avg[0] == '\0') {
    return;
  }
  FILE *fw = PrintBylineOpenFile(opt.f_avg, argc, argv);
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
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    fprintf(fw, "(%d) <%s>_n, ", count++, System.MoleculeType[i].Name);
  }
  fprintf(fw, "(%d) n_agg", count++);
  putc('\n', fw);
  PrintHeaderMassNote(fw, System, opt);
  fclose(fw);
} //}}}
// print header for a distr output file (-d option) //{{{
static void PrintDistrHeader(int argc, char *argv[], OPT opt, SYSTEM System) {
  if (opt.f_distr[0] == '\0') {
    return;
  }
  FILE *fw = PrintBylineOpenFile(opt.f_distr, argc, argv);
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
  PrintHeaderMassNote(fw, System, opt);
  fclose(fw);
} //}}}
// print header for a simple copmposition output file (-c option) //{{{
static void PrintCompSimpleHeader(int argc, char *argv[], OPT opt,
                                  SYSTEM System, int size,
                                  long int *comp_agg_count) {
  char file[LINE];
  FilenameCompSimple(opt, size, file);
  FILE *fw = PrintBylineOpenFile(file, argc, argv);
  // print header
  fprintf(fw, "# total number of aggregates with size %d: %ld\n",
          opt.comp.size[size], comp_agg_count[size]);
  fprintf(fw, "# (1) number of molecules of given type;");
  fprintf(fw, " fraction of aggregates with that many molecules of type:");
  for (int j = 0; j < System.Count.MoleculeType; j++) {
    fprintf(fw, " (%d) %s", j + 2, System.MoleculeType[j].Name);
    if (j != (System.Count.MoleculeType - 1)) {
      putc(',', fw);
    }
  }
  putc('\n', fw);
  fclose(fw);
} //}}}
// print header for a 2D copmposition output file (-c option) //{{{
static void PrintComp2DHeader(int argc, char *argv[], OPT opt, SYSTEM System,
                              int size, long int *comp_agg_count) {
  char file[LINE];
  FilenameCode2D(opt, size, file);
  FILE *fw = PrintBylineOpenFile(file, argc, argv);
  // print header
  fprintf(fw, "# total number of aggregates with size %d: %ld\n",
          opt.comp.size[size], comp_agg_count[size]);
  fprintf(fw, "# (1-2) number of molecules:");
  int count = 3;
  COUNT *Count = &System.Count;
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
  fclose(fw);
}
//}}}
// print the note about the two different mass definition //{{{
static void PrintHeaderMassNote(FILE *fw, SYSTEM System, OPT opt) {
  if (!opt.agg.m_flag) {
    fprintf(fw, "# Note: The -m option was not used; therefore, 'partial mass'"
            " includes all molecules, and "
            "<As>_w/z (partial mass) = <As>_w/z (total mass)\n");
  } else {
    fprintf(fw, "# Note: 'partial mass' includes molecules specifid by -m (");
    bool first = true;
    for (int i = 0; i < System.Count.MoleculeType; i++) {
      if (opt.agg.m[i]) {
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
} //}}}
// append overall averages to file(s) (-a/-d options) //{{{
static void AppendOverallAvg(char *f, SYSTEM System, double As_sum[3][2],
                             double mass_sum[3][2], int *count_agg_per_size,
                             int timesteps, ArrNDi *molecules_sum) {
  int *mol_sum_per_size = calloc(System.Count.MoleculeType,
                                 sizeof *mol_sum_per_size);
  if (!mol_sum_per_size) {
    ErrorAlloc("mol_sum_per_size");
  }
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
// As-based filenames (-c option) //{{{
static void FilenameCompSimple(OPT opt, int size, char filename[LINE]) {
  if (snprintf(filename, LINE, "%s-%03d.txt",
               opt.comp.f, opt.comp.size[size]) < 0) {
    ErrorSnprintf();
  }
}
static void FilenameCode2D(OPT opt, int size, char filename[LINE]) {
  if (snprintf(filename, LINE, "%s_r-%03d.txt",
               opt.comp.f, opt.comp.size[size]) < 0) {
    ErrorSnprintf();
  }
} //}}}
