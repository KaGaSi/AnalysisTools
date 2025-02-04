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

  fprintf(ptr, "<input>             input structure file\n");
  fprintf(ptr, "<in.agg>            input agg file\n");
  fprintf(ptr, "<distr file>        output file with distributions\n");
  fprintf(ptr, "<avg file>          output file with timestep averages\n");
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
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 7, all = common + 5, count = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "--verbose", "--silent", "--help",
               "--version", "-n", "-m", "-x", "-only", "-c");

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
  int c_sizes[100] = {0}, c_count = 0;
  char c_file[LINE] = "";
  FileNumbersOption(argc, argv, 1, 100, "-c", c_sizes, &c_count, c_file, 'i');
  //}}}

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-n' option //{{{
  int range_As[2] = {1, Count->Molecule};
  TwoNumbersOption(argc, argv, "-n", range_As, 'i');
  if (range_As[0] > range_As[1]) {
    SwapInt(&range_As[0], &range_As[1]);
  } //}}}
  // '-m' option //{{{
  bool *mtype_As = calloc(Count->MoleculeType, sizeof *mtype_As);
  if (!TypeOption(argc, argv, "-m", 'm', true, mtype_As, System)) {
    InitBoolArray(mtype_As, Count->MoleculeType, true);
  } //}}}
  // '-only' option //{{{
  bool *mtype_only_opt = calloc(Count->MoleculeType, sizeof *mtype_only_opt);
  if (!TypeOption(argc, argv, "-only", 'm', true, mtype_only_opt, System)) {
    InitBoolArray(mtype_only_opt, Count->MoleculeType, true);
  } //}}}
  // error - molecules specified by -m and -only do not overlap //{{{
  bool overlap = false;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (mtype_As[i] && mtype_only_opt[i]) {
      overlap = true;
      break;
    }
  }
  if (!overlap) {
    err_msg("for any aggregate to be used, at least one molecule "
            "must be specified in both options");
    PrintErrorOption("-m/-only");
    exit(1);
  } //}}}
  // '-x' option //{{{
  bool *mtype_x_opt = calloc(Count->MoleculeType, sizeof *mtype_x_opt);
  TypeOption(argc, argv, "-x", 'm', true, mtype_x_opt, System); //}}}
  // error - molecules specified by -only and -x must differ //{{{
  overlap = true; // do the two array fully overlap?
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (mtype_x_opt[i] != mtype_only_opt[i]) {
      overlap = false;
      break;
    }
  }
  if (overlap) {
    err_msg("the lists of molecules must be different");
    PrintErrorOption("-x/-only");
    exit(1);
  } //}}}
  // error - -x specifies all molecules in the system //{{{
  overlap = true; // are all molecule types specified by -x?
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!mtype_x_opt[i]) {
      overlap = false;
      break;
    }
  }
  if (overlap) {
    err_msg("with all molecules listed, no aggregates would be detected");
    PrintErrorOption("-x");
    exit(1);
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
  double width_r = 0.1;
  double bin_r = Count->Molecule / width_r + 1; // +1 for N/0
  if (c_count > 0) {
    link_c_sizes = malloc(Count->Molecule * sizeof *link_c_sizes);
    InitIntArray(link_c_sizes, Count->Molecule, -1);
    comp_distr = malloc(c_count * sizeof *comp_distr);
    ratio_distr = malloc(c_count * sizeof *ratio_distr);
    comp_agg_count = calloc(c_count, sizeof *comp_agg_count);
    for (int i = 0; i < c_count; i++) {
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
        if (j == c_sizes[i]) {
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
  fprintf(fw, "(%d) <M>_n, ", count++);
  fprintf(fw, "(%d) <M>_w, ", count++);
  fprintf(fw, "(%d) <M>_z, ", count++);
  fprintf(fw, "(%d) <As>_n, ", count++);
  fprintf(fw, "(%d) <As>_w, ", count++);
  fprintf(fw, "(%d) <As>_z, ", count++);
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
   *   [][0] = mass of mols in agg from options - TODO: implemen?
   *   [][1] = mass of the whole aggregate
   */
  double mass_sum[3][2] = {{0}}, As_sum[3][2] = {{0}};
  while (true) { // cycle ends with 'Last Step' line in agg file
    PrintStep(&count_step, opt->c.start, opt->c.silent);
    if (ReadAggregates(fr, input_agg, &System, Aggregate, &agg_lines) < 0) {
      count_step--;
      break;
    }

    // decide whether this timestep is to be used for averages and distributions
    bool use = false;
    if (UseStep(opt->c, count_step)) {
      use = true;
    }
    if (use) { //{{{
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
        bool only_opt = true, // acceptable composition (-only option)
             x_opt = true; // acceptable composition (-x option)
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mtype = System.Molecule[Aggregate[i].Molecule[j]].Type;
          if (mtype_As[mtype]) {
            size++;
            agg_mass += System.MoleculeType[mtype].Mass;
          }
          if (!mtype_only_opt[mtype]) {
            only_opt = false;
          }
          if (mtype_x_opt[mtype]) {
            x_opt = false;
          }
        }
        if (size == 0 || size < range_As[0] || size > range_As[1] ||
            !only_opt || !x_opt) {
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
          if (c_count > 0 && link_c_sizes[size] != -1) {
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
                a /= width_r;
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
        fprintf(fw, " %10.5f", avg_mass_n_step[0]/aggs_step); // <mass>_n
        fprintf(fw, " %10.5f", avg_mass_w_step[0]/avg_mass_n_step[0]); // <mass>_w
        fprintf(fw, " %10.5f", avg_mass_z_step[0]/avg_mass_w_step[0]); // <mass>_z
        fprintf(fw, " %10.5f", avg_As_n_step[0]/aggs_step); // <As>_n
        fprintf(fw, " %10.5f", avg_As_w_step[0]/avg_mass_n_step[0]); // <As>_w
        fprintf(fw, " %10.5f", avg_As_z_step[0]/avg_mass_w_step[0]); // <As>_z
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
    }
      //}}}
    // } else {
    //   if (!SkipTimestep(in, fr, &agg_lines)) {
    //     count_step--;
    //     break;
    //   }
    // }
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
    fprintf(fw, " (%d) <%s>_n", i+count, System.MoleculeType[i].Name);
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

  for (int i = 0; i < Count->Molecule; i++) {
    if (count_agg[i] > 0) {
      fprintf(fw, "%4d", i+1); // As
      fprintf(fw, " %lf", (double)(ndistr[i]) / ndistr_norm);
      fprintf(fw, " %lf", (double)(wdistr[i][0]) / wdistr_norm[0]);
      fprintf(fw, " %lf", (double)(zdistr[i][0]) / zdistr_norm[0]);
      fprintf(fw, " %10d", count_agg[i]); // number of aggregates
      // print average number of molecule types in aggregates
      for (int j = 0; j < Count->MoleculeType; j++) {
        fprintf(fw, " %10.5f", (double)(molecules_sum[i][j])/count_agg[i]);
      }
      putc('\n', fw);
    }
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

  // print composition distribution(s) (-c option) //{{{
  if (c_count > 0) {
    for (int i = 0; i < c_count; i++) {
      char file[LINE];
      if (snprintf(file, LINE, "%s-%03d.txt", c_file, c_sizes[i]) < 0) {
        ErrorSnprintf();
      }
      PrintByline(file, argc, argv);
      fw = OpenFile(file, "a");
      fprintf(fw, "# total number of aggregates with size %d: %ld\n",
              c_sizes[i], comp_agg_count[i]);
      fprintf(fw, "# (1) number of molecules of given type;");
      fprintf(fw, " fraction of aggregates with that many molecules of type:");
      for (int j = 0; j < Count->MoleculeType; j++) {
        fprintf(fw, " (%d) %s", j + 2, System.MoleculeType[j].Name);
        if (j != (Count->MoleculeType - 1)) {
          putc(',', fw);
        }
      }
      putc('\n', fw);
      if (comp_agg_count[i] > 0) {
        // print the distribution
        for (int j = 0; j <= c_sizes[i]; j++) {
          fprintf(fw, "%3d", j);
          for (int k = 0; k < Count->MoleculeType; k++) {
            fprintf(fw, " %lf",
                    (double)(comp_distr[i][k][j]) / comp_agg_count[i]);
          }
          putc('\n', fw);
        }
      } else {
        snprintf(ERROR_MSG, LINE, "no aggregates with size %s%d%s found",
                 ErrYellow(), c_sizes[i], ErrCyan());
        PrintWarning();
      }
      fclose(fw);
      if (snprintf(file, LINE, "%s-ratio_%03d.txt", c_file, c_sizes[i]) < 0) {
        ErrorSnprintf();
      }
      PrintByline(file, argc, argv);
      fw = OpenFile(file, "a");
      fprintf(fw, "# total number of aggregates with size %d: %ld\n",
              c_sizes[i], comp_agg_count[i]);
      fprintf(fw, "# (1) ratio of:");
      count = 2;
      for (int j = 0; j < (Count->MoleculeType-1); j++) {
        for (int k = (j+1); k < Count->MoleculeType; k++) {
          fprintf(fw, " (%d) %s/%s", count++, System.MoleculeType[j].Name,
                                              System.MoleculeType[k].Name);
          if (!(j == (Count->MoleculeType - 2) && k == (j + 1))) {
            putc(',', fw);
          }
        }
      }
      putc('\n', fw);
      if (comp_agg_count[i] > 0) {
        // print the distribution
        for (int j = 0; j <= (c_sizes[i] / width_r); j++) {
          fprintf(fw, "%lf", j * width_r);
          for (int k = 0; k < Count->MoleculeType; k++) {
            for (int l = (k + 1); l < Count->MoleculeType; l++) {
              fprintf(fw, " %lf",
                      (double)(ratio_distr[i][k][l][j])/comp_agg_count[i]);
            }
          }
          putc('\n', fw);
        }
      }
      fclose(fw);
    }
  }
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  for (int i = 0; i < Count->Molecule; i++) {
    free(molecules_sum[i]);
  }
  free(molecules_sum);
  free(mtype_As);
  free(mtype_only_opt);
  free(mtype_x_opt);
  if (c_count > 0) {
    for (int i = 0; i < c_count; i++) {
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
