#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
        TODO: rewrite\n\n\
DistrAgg calculates average aggregation numbers and aggregate masses during \
the simulation run (i.e., time evolution) as well as overall distributions. \
The definition of aggregate size is quite \
flexible and only a specified range can be used. Also, this utility can \
fanalyse composition of specified aggreget size(s) and write compositition \
distribution (i.e., distribution of numbers of different molecules in over \
all aggregates with given size).\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <in.agg> <distr file> <avg file> "
          "[options]\n\n", cmd);

  fprintf(ptr, "<in.stru>           input structure file\n");
  fprintf(ptr, "<in.agg>            input agg file\n");
  fprintf(ptr, "<distr>             output file with distributions\n");
  fprintf(ptr, "<avg>               output file with timestep averages\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -n <size> <size>  use aggregate sizes in a given range\n");
  fprintf(ptr, "  -m <name(s)>      use number of specified molecule type(s) "
          "as aggrete size\n");
  fprintf(ptr, "  -x <name(s)>      exclude aggregates containing only "
          "specified molecule(s)\n");
  fprintf(ptr, "  -only <name(s)>   use only aggregates composed of "
          "specified molecule type(s)\n");
  fprintf(ptr, "  -c <file> <size(s)>\n");
  fprintf(ptr, "                    write composition distributions for "
          "aggregate size(s) to two <file>s with automatic endings '-#.txt' "
          "and '-ratio_#.txt'\n");
  fprintf(ptr, "  -w <width>        bin width for the ratios in -c option "
          "(default: 0.1)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool x, only, m;                // were -x, -only, and/or -m specified?
  bool *x_mol, *only_mol, *m_mol; // arrays for which molecules to use
  int range[2],                   // -n; used with 1 & Count.Molecules if no -n
      c_s[100],                   // -c; aggregate sizes
      c_c;                        // -c; number of sizes
  char c_f[LINE];                 // -c; filename
  double w;                       // -w; bin width for -c
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 7, all = common + 6, count = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "--verbose", "--silent", "--help",
               "--version", "-n", "-m", "-x", "-only", "-c", "-w");

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
  // <distr file> - file with distribution of aggregation numbers
  char out_distr[LINE] = "";
  s_strcpy(out_distr, argv[++count], LINE);
  // <avg file> - file with per-timestep average aggregation numbers
  char out_avg[LINE] = "";
  s_strcpy(out_avg, argv[++count], LINE);
  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);
  // -c option
  FileNumbersOption(argc, argv, 1, 100, "-c", opt->c_s,
                    &opt->c_c, opt->c_f, 'i');
  if (!OneNumberOption(argc, argv, "-w", &opt->w, 'd')) {
    opt->w = 0.1;
  } else if (opt->c_c == 0) {
    err_msg("no effect when -c option is missing");
    PrintWarnOption("-w");
  } //}}}

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-n' option //{{{
  opt->range[0] = 1;
  opt->range[1] = Count->Molecule;
  TwoNumbersOption(argc, argv, "-n", opt->range, 'i');
  if (opt->range[0] > opt->range[1]) {
    SwapInt(&opt->range[0], &opt->range[1]);
  } //}}}
  // '-m' option //{{{
  opt->m_mol = calloc(Count->MoleculeType, sizeof *opt->m_mol);
  opt->m = true;
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->m_mol, System)) {
    opt->m = false;
    InitBoolArray(opt->m_mol, Count->MoleculeType, true);
  } //}}}
  // '-only' option //{{{
  opt->only_mol = calloc(Count->MoleculeType, sizeof *opt->only_mol);
  opt->only = true;
  if (!TypeOption(argc, argv, "-only", 'm', true, opt->only_mol, System)) {
    opt->only = false;
    InitBoolArray(opt->only_mol, Count->MoleculeType, true);
  } //}}}
  // '-x' option //{{{
  opt->x_mol = calloc(Count->MoleculeType, sizeof *opt->x_mol);
  opt->x = true;
  if (!TypeOption(argc, argv, "-x", 'm', true, opt->x_mol, System)) {
    opt->x = false;
  }
  // error - all molecule specified //{{{
  if (opt->x) {
    bool overlap = true; // are all molecule types specified by -x?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (!opt->x_mol[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("with all molecules listed, no aggregates would be detected");
      PrintErrorOption("-x");
      exit(1);
    }
  } //}}}
  //}}}
  // error - molecules specified by -m and -only do not overlap //{{{
  if (opt->m && opt->only) {
    bool overlap = false;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->m_mol[i] && opt->only_mol[i]) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      err_msg("for any aggregate to be used, at least one molecule "
              "must be specified in both options");
      PrintErrorOption("-m/-only");
      exit(1);
    }
  } //}}}
  // error - molecules specified by -only and -x must differ //{{{
  if (opt->only && opt->x) {
    bool overlap = true; // do the two array fully overlap?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->x_mol[i] != opt->only_mol[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("the lists of molecules must be different");
      PrintErrorOption("-x/-only");
      exit(1);
    }
  } //}}}

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // arrays for distributions //{{{
  // number distribution
  long double ndistr[Count->Molecule];
  /* weight and z distributions:
   *   [][0] = mass of mols according to options - TODO: implement?
   *   [][1] = mass of whole agg
   */
  long double (*wdistr)[2] = calloc(Count->Molecule, sizeof *wdistr);
  long double (*zdistr)[2] = calloc(Count->Molecule, sizeof *zdistr);
  // number of aggregates throughout simulation
  int *count_agg = calloc(Count->Molecule, sizeof *count_agg);
  // molecule typs in aggregates: [agg size][mol type][number or Square(number)]
  int **molecules_sum = malloc(Count->Molecule * sizeof *molecules_sum);
  for (int i = 0; i < Count->Molecule; i++) {
    molecules_sum[i] = calloc(Count->MoleculeType, sizeof *molecules_sum[i]);
  }
  // arrays for composition distribution
  long int ***comp_distr = NULL; // [c_size][moltype][number of mols]
  long int ****ratio_distr = NULL; // [c_size][moltype1][moltype2][ratio]
  long int *comp_agg_count = NULL;
  int *link_c_sizes = NULL;
  double bin_r = Count->Molecule / opt->w + 1; // +1 for N/0
  if (opt->c_c > 0) {
    link_c_sizes = malloc(Count->Molecule * sizeof *link_c_sizes);
    InitIntArray(link_c_sizes, Count->Molecule, -1);
    comp_distr = malloc(opt->c_c * sizeof *comp_distr);
    ratio_distr = malloc(opt->c_c * sizeof *ratio_distr);
    comp_agg_count = calloc(opt->c_c, sizeof *comp_agg_count);
    for (int i = 0; i < opt->c_c; i++) {
      comp_distr[i] = malloc(Count->MoleculeType * sizeof *comp_distr[i]);
      ratio_distr[i] = malloc(Count->MoleculeType * sizeof *ratio_distr[i]);
      for (int j = 0; j < Count->MoleculeType; j++) {
        // +1 as it goes from no molecules to N molecules in the agg
        comp_distr[i][j] = calloc(Count->Molecule + 1,
                                  sizeof *comp_distr[i][j]);
        ratio_distr[i][j] = malloc(Count->MoleculeType *
                                   sizeof *ratio_distr[i][j]);
        for (int k = 0; k < Count->MoleculeType; k++) {
          ratio_distr[i][j][k] = calloc(bin_r, sizeof *ratio_distr[i][j][k]);
        }
      }
      for (int j = 0; j < Count->Molecule; j++) {
        if (j == opt->c_s[i]) {
          link_c_sizes[j] = i;
        }
      }
    }
  }

  // zeroize arrays
  for (int i = 0; i < Count->Molecule; i++) {
    ndistr[i] = 0;
    wdistr[i][0] = 0;
    wdistr[i][1] = 0;
    zdistr[i][0] = 0;
    zdistr[i][1] = 0;
    count_agg[i] = 0;
  } //}}}

  // print the first two lines to output file with per-step averages //{{{
  PrintByline(out_avg, argc, argv);
  FILE *fw = OpenFile(out_avg, "a");
  count = 1;
  fprintf(fw, "# column: (%d) step, ", count++);
  fprintf(fw, "(%d) <As>_n, ", count++);
  fprintf(fw, "(%d) <As>_w, ", count++);
  fprintf(fw, "(%d) <As>_z, ", count++);
  fprintf(fw, "(%d) <M>_n, ", count++);
  fprintf(fw, "(%d) <M>_w, ", count++);
  fprintf(fw, "(%d) <M>_z, ", count++);
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(fw, "(%d) <%s>_n, ", count++, System.MoleculeType[i].Name);
  }
  fprintf(fw, "(%d) n_agg", count++);
  putc('\n', fw);
  fclose(fw); //}}}

  // open <in.agg> and skip the first two lines //{{{
  FILE *fr = OpenFile(input_agg, "r");
  while (getc(fr) != '\n')
    ;
  while (getc(fr) != '\n')
    ; //}}}

  // main loop //{{{
  int count_step = 0,
      count_used = 0,
      agg_lines = 2; // first two lines already read (skipped)
  /*
   * Mass and aggregate size sums
   *   [0][] = simple sum, [1][] = sum of squares, [2][] = sume of cubes
   *   [][0] = mass of mols in agg from options - TODO: implement?
   *   [][1] = mass of the whole aggregate
   */
  double mass_sum[3][2] = {{0}}, As_sum[3][2] = {{0}};
  while (true) { // cycle ends with 'Last Step' line in agg file
    PrintStep(&count_step, opt->c.start, opt->c.silent);

    // decide whether this timestep is to be used for averages and distributions
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
             avg_mass_w_step[2] = {0}, // [0] ... from options TODO: implement?
             avg_mass_z_step[2] = {0}, // [1] ... for whole aggregates
             avg_As_n_step[2] = {0}, // per-step As averages
             avg_As_w_step[2] = {0}, // [0] ... from options TODO: implement?
             avg_As_z_step[2] = {0}, // [1] ... for whole aggregates
             molecules_step[Count->MoleculeType];
      // zeroize per-step counts of molecule types
      InitDoubleArray(molecules_step, Count->MoleculeType, 0);
      for (int i = 0; i < Count->Aggregate; i++) { //{{{
        // decide whether to use the aggregate based on used options //{{{
        int size = 0; // -m option-adjusted aggregate size
        double agg_mass = 0; // -m option-adjusted aggregate mass
        bool only_opt = true, // acceptable composition (-only option)?
             x_opt = false; // acceptable composition (-x option)?
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          MOLECULE *mol = &System.Molecule[Aggregate[i].Molecule[j]];
          int mtype = mol->Type;
          MOLECULETYPE *mt = &System.MoleculeType[mtype];
          if (opt->m_mol[mtype]) {
            size++;
            agg_mass += mt->Mass;
          }
          // if at least one unwanted molecule is present, don't use aggregate
          if (!opt->only_mol[mtype]) {
            only_opt = false;
          }
          // if at least one molecule isn't exluded, use aggregate
          if (!opt->x_mol[mtype]) {
            x_opt = true;
          }
        }
        if (size == 0 || // -m: discarded all molecules
            size < opt->range[0] || size > opt->range[1] || // -n: not in range
            !only_opt || // -only: found molecule that weren't supposed to be in
            !x_opt) { // -x: didn't find any un-excluded molecules
          continue;
        } //}}}
        // average aggregate mass during the step
        avg_mass_n_step[0] += agg_mass;
        avg_mass_w_step[0] += Square(agg_mass);
        avg_mass_z_step[0] += Cube(agg_mass);
        avg_mass_n_step[1] += Aggregate[i].Mass;
        avg_mass_w_step[1] += Square(Aggregate[i].Mass);
        avg_mass_z_step[1] += Cube(Aggregate[i].Mass);
        // average aggregation number during the step
        avg_As_n_step[0] += size;
        avg_As_w_step[0] += size * agg_mass;
        avg_As_z_step[0] += size * Square(agg_mass);
        avg_As_n_step[1] += Aggregate[i].nMolecules;
        avg_As_w_step[1] += Aggregate[i].nMolecules * Aggregate[i].Mass;
        avg_As_z_step[1] += Aggregate[i].nMolecules * Square(Aggregate[i].Mass);
        // molecule species numbers
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mtype = System.Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_step[mtype]++;
        }

        aggs_step++;

        // use the step for averages and distributions?
        if (use) {
          count_agg[size-1]++;
          // distribution
          ndistr[size-1]++;
          wdistr[size-1][0] += agg_mass;
          wdistr[size-1][1] += Aggregate[i].Mass;
          zdistr[size-1][0] += Square(agg_mass);
          zdistr[size-1][1] += Square(Aggregate[i].Mass);
          // summed up sizes
          As_sum[0][0] += size;
          As_sum[1][0] += size * agg_mass;
          As_sum[2][0] += size * Square(agg_mass);
          As_sum[0][1] += Aggregate[i].nMolecules;
          As_sum[1][1] += Aggregate[i].nMolecules * Aggregate[i].Mass;
          As_sum[2][1] += Aggregate[i].nMolecules * Square(Aggregate[i].Mass);
          // summed up masses
          mass_sum[0][0] += agg_mass;
          mass_sum[1][0] += Square(agg_mass);
          mass_sum[2][0] += Cube(agg_mass);
          mass_sum[0][1] += Aggregate[i].Mass;
          mass_sum[1][1] += Square(Aggregate[i].Mass);
          mass_sum[2][1] += Cube(Aggregate[i].Mass);
          // molecule species numbers
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
            molecules_sum[size-1][mol_type]++;
          }
          // composition distribution (-c option)
          if (opt->c_c > 0 && link_c_sizes[size] != -1) {
            comp_agg_count[link_c_sizes[size]]++;
            double comp_aux[Count->MoleculeType];
            InitDoubleArray(comp_aux, Count->MoleculeType, 0);
            // count molecule types in the aggregate
            for (int j = 0; j < Aggregate[i].nMolecules; j++) {
              int mtype = System.Molecule[Aggregate[i].Molecule[j]].Type;
              comp_aux[mtype]++;
            }
            // increment the distribution
            for (int j = 0; j < Count->MoleculeType; j++) {
              int id = link_c_sizes[size];
              comp_distr[id][j][(int)comp_aux[j]]++;
              for (int k = (j+1); k < Count->MoleculeType; k++) {
                double a;
                if (comp_aux[k] == 0) {
                  a = size;
                } else {
                  a = comp_aux[j] / comp_aux[k];
                }
                a /= opt->w;
                ratio_distr[id][j][k][(int)a]++;
              }
            }
          }
        }
      } //}}}
      // print averages to output file //{{{
      fw = OpenFile(out_avg, "a");
      fprintf(fw, "%5d", count_step); // step
      if (aggs_step > 0) {
        fprintf(fw, " %10.5f", avg_As_n_step[0]/aggs_step); // <As>_n
        fprintf(fw, " %10.5f", avg_As_w_step[0]/avg_mass_n_step[0]); // <As>_w
        fprintf(fw, " %10.5f", avg_As_z_step[0]/avg_mass_w_step[0]); // <As>_z
        fprintf(fw, " %10.5f", avg_mass_n_step[0]/aggs_step); // <mass>_n
        fprintf(fw, " %10.5f", avg_mass_w_step[0]/avg_mass_n_step[0]); // <mass>_w
        fprintf(fw, " %10.5f", avg_mass_z_step[0]/avg_mass_w_step[0]); // <mass>_z
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
      fclose(fw); //}}}
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
    fprintf(stdout, " (%d used for distributions and overall averages)\n",
            count_used);
  } //}}}
  //}}}

  // print the first two lines to output file with distributions //{{{
  PrintByline(out_distr, argc, argv);
  fw = OpenFile(out_distr, "a");
  count = 1;
  fprintf(fw, "# column: ");
  fprintf(fw, "(%d) As, ", count++);
  fprintf(fw, "(%d) F_n(As), ", count++);
  fprintf(fw, "(%d) F_w(As), ", count++);
  fprintf(fw, "(%d) F_z(As), ", count++);
  fprintf(fw, "(%d) n_agg,", count++);
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(fw, " (%d) <%s>_n", i + count, System.MoleculeType[i].Name);
    if (i != (Count->MoleculeType-1)) {
      putc(',', fw);
    }
  }
  putc('\n', fw);
  fclose(fw); //}}}

  // print distributions to output file //{{{
  fw = OpenFile(out_distr, "a");

  if (opt->c.end == -1) {
    count_step = count_step - opt->c.start + 1;
  } else {
    count_step = count_step - (opt->c.start - 1) - (opt->c.end - 1);
  }

  // normalization factors
  long int ndistr_norm = 0, wdistr_norm[2] = {0}, zdistr_norm[2] = {0};
  for (int i = 0; i < Count->Molecule; i++) {
    ndistr_norm += ndistr[i];
    wdistr_norm[0] += wdistr[i][0];
    wdistr_norm[1] += wdistr[i][1];
    zdistr_norm[0] += zdistr[i][0];
    zdistr_norm[1] += zdistr[i][1];
  }

  // determine width of each column //{{{
  int columns = Count->MoleculeType + 5;
  int digits[columns][2];
  InitInt2DArray((int *)digits, columns, 2, 0);
  double *data[Count->Molecule]; // array for data
  for (int i = 0; i < Count->Molecule; i++) {
    data[i] = calloc(columns, sizeof data[i]);
    if (count_agg[i] > 0) {
      count = -1;
      data[i][++count] = i + 1;
      data[i][++count] = (double)(ndistr[i]) / ndistr_norm;
      data[i][++count] = (double)(wdistr[i][0]) / wdistr_norm[0];
      data[i][++count] = (double)(zdistr[i][0]) / zdistr_norm[0];
      data[i][++count] = count_agg[i];
      for (int j = 0; j < Count->MoleculeType; j++) {
        data[i][++count] = (double)(molecules_sum[i][j])/count_agg[i];
      }
    }
  }
  FillMaxDigits(columns, Count->Molecule, data, digits); //}}}
  // print data (only for aggregate sizes that actually exist)
  for (int i = 0; i < Count->Molecule; i++) {
    if (count_agg[i] > 0) {
      for (int col = 0; col < columns; col++) {
        Fprintf1(fw, data[i][col], digits[col]);
      }
      putc('\n', fw);
    }
    free(data[i]);
  }
  fclose(fw); //}}}

  // print overall averages to avg and distr output files //{{{
  // count total numbers of molecules of each type and of aggregates
  for (int i = 1; i < Count->Molecule; i++) {
    count_agg[0] += count_agg[i];
    for (int j = 0; j < Count->MoleculeType; j++) {
      molecules_sum[0][j] += molecules_sum[i][j];
    }
  }

  FILE *f[2];
  f[0] = OpenFile(out_distr, "a");
  f[1] = OpenFile(out_avg, "a");
  for (int i = 0; i < 2; i++) { // go over the two files
    // print legend (with column numbers)
    putc('#', f[i]);
    // distr file
    count = 1;
    fprintf(f[i], " (%d) <As>_n,", count++);
    fprintf(f[i], " (%d) <As>_w,", count++);
    fprintf(f[i], " (%d) <As>_z,", count++);
    fprintf(f[i], " (%d) <M>_n,", count++);
    fprintf(f[i], " (%d) <M>_w,", count++);
    fprintf(f[i], " (%d) <M>_z,", count++);
    for (int j = 0; j < Count->MoleculeType; j++) {
      fprintf(f[i], " (%d) <%s>_n,", count++, System.MoleculeType[j].Name);
    }
    fprintf(f[i], " (%d) <n_agg>", count++);
    putc('\n', f[i]);
    // print the averages
    fprintf(f[i], "#");
    if (count_agg[0] > 0) {
      fprintf(f[i], " %lf", As_sum[0][0]/count_agg[0]); // <As>_n
      fprintf(f[i], " %lf", As_sum[1][0]/mass_sum[0][0]); // <As>_w
      fprintf(f[i], " %lf", As_sum[2][0]/mass_sum[1][0]); // <As>_z

      fprintf(f[i], " %lf", mass_sum[0][0]/count_agg[0]); // <M>_n
      fprintf(f[i], " %lf", mass_sum[1][0]/mass_sum[0][0]); // <M>_w
      fprintf(f[i], " %lf", mass_sum[2][0]/mass_sum[1][0]); // <M>_z
      for (int j = 0; j < Count->MoleculeType; j++) {
        // <species>_n
        fprintf(f[i], " %lf", (double)(molecules_sum[0][j])/count_agg[0]);
      }
      fprintf(f[i], " %lf", (double)(count_agg[0])/count_step); // <n_agg>
    } else { // zero everywhere if no aggregates found
      fprintf(f[i], " 0.0  0.0  0.0  0.0  0.0  0.0  0.0");
      for (int j = 0; j < Count->MoleculeType; j++) {
        fprintf(f[i], "  0.0");
      }
    }
    putc('\n', f[i]);
  }
  fclose(f[0]);
  fclose(f[1]);
  //}}}

  // print composition distribution(s) (-c option)
  if (opt->c_c > 0) {
    for (int i = 0; i < opt->c_c; i++) {
      // print the distribution //{{{
      char file[LINE];
      if (snprintf(file, LINE, "%s-%03d.txt", opt->c_f, opt->c_s[i]) < 0) {
        ErrorSnprintf();
      }
      PrintByline(file, argc, argv);
      fw = OpenFile(file, "a");
      // print header
      fprintf(fw, "# total number of aggregates with size %d: %ld\n",
              opt->c_s[i], comp_agg_count[i]);
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
        // determine width of each column & collate data //{{{
        int columns = Count->MoleculeType + 1;
        int digits[columns][2];
        InitInt2DArray((int *)digits, columns, 3, 0);
        double *data[opt->c_s[i]+1]; // array for data
        for (int j = 0; j <= opt->c_s[i]; j++) {
          data[j] = calloc(columns, sizeof data[j]);
          count = -1;
          data[j][++count] = j;
          for (int k = 0; k < Count->MoleculeType; k++) {
            data[j][++count] = (double)(comp_distr[i][k][j]) /
                               comp_agg_count[i];
          }
        }
        FillMaxDigits(columns, opt->c_s[i], data, digits); //}}}
        for (int j = 0; j <= opt->c_s[i]; j++) {
          WriteFormatedDataLine(fw, columns, data[j], digits);
          free(data[j]);
        }
      } else {
        snprintf(ERROR_MSG, LINE, "no aggregates with size %s%d%s found",
                 ErrYellow(), opt->c_s[i], ErrCyan());
        PrintWarning();
      }
      fclose(fw); //}}}
      // print the ratios //{{{
      if (snprintf(file, LINE, "%s-ratio_%03d.txt",
                   opt->c_f, opt->c_s[i]) < 0) {
        ErrorSnprintf();
      }
      PrintByline(file, argc, argv);
      fw = OpenFile(file, "a");
      // print header
      fprintf(fw, "# total number of aggregates with size %d: %ld\n",
              opt->c_s[i], comp_agg_count[i]);
      fprintf(fw, "# (1) ratio of:");
      count = 2;
      for (int j = 0; j < (Count->MoleculeType - 1); j++) {
        for (int k = (j + 1); k < Count->MoleculeType; k++) {
          fprintf(fw, " (%d) %s/%s", count++, System.MoleculeType[j].Name,
                                              System.MoleculeType[k].Name);
          if (!(j == (Count->MoleculeType - 2) && k == (j + 1))) {
            putc(',', fw);
          }
        }
      }
      putc('\n', fw);
      // print data
      if (comp_agg_count[i] > 0) {
        // determine width of each column & collate data //{{{
        int columns = Count->MoleculeType * (Count->MoleculeType - 1) / 2 + 1;
        int digits[columns][2];
        int lines = opt->c_s[i] / opt->w + 1;
        InitInt2DArray((int *)digits, columns, 2, 0);
        double *data[lines]; // array for data
        for (int j = 0; j < lines; j++) {
          data[j] = calloc(columns, sizeof data[j]);
          count = -1;
          // data[j][++count] = j * opt->w + opt->w / 2;
          data[j][++count] = opt->w * (j + 0.5);
          for (int k = 0; k < Count->MoleculeType; k++) {
            for (int l = (k + 1); l < Count->MoleculeType; l++) {
              data[j][++count] = (double)(ratio_distr[i][k][l][j]) /
                                 comp_agg_count[i];
            }
          }
        }
        FillMaxDigits(columns, lines, data, digits); //}}}
        for (int j = 0; j < lines; j++) {
          for (int col = 0; col < columns; col++) {
            if (fabs(data[j][col]) > 1e-5) {
              Fprintf1(fw, data[j][col], digits[col]);
            } else {
              fprintf(fw, "%*s", digits[col][0] + digits[col][1] + 1, "?");
            }
          }
          putc('\n', fw);
          // WriteFormatedDataLine(fw, columns, data[j], digits);
          free(data[j]);
        }
      }
      fclose(fw); //}}}
    }
  }

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  for (int i = 0; i < Count->Molecule; i++) {
    free(molecules_sum[i]);
  }
  free(molecules_sum);
  free(opt->m_mol);
  free(opt->only_mol);
  free(opt->x_mol);
  if (opt->c_c > 0) {
    for (int i = 0; i < opt->c_c; i++) {
      for (int j = 0; j < Count->MoleculeType; j++) {
        for (int k = 0; k < Count->MoleculeType; k++) {
          free(ratio_distr[i][j][k]);
        }
        free(comp_distr[i][j]);
        free(ratio_distr[i][j]);
      }
      free(comp_distr[i]);
      free(ratio_distr[i]);
    }
    free(comp_distr);
    free(ratio_distr);
    free(comp_agg_count);
    free(link_c_sizes);
  }
  free(wdistr);
  free(zdistr);
  free(count_agg);
  free(opt); //}}}

  return 0;
}
