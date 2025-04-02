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
GyrationAggregates calculates gyration tensor for aggregates and determines \
shape descriptors like radius of gyration, acylindricity, asphericity, or \
relative shape anisotropy. By default, it calculates per-timestep averages, \
but per-size averages can also be determined. Overall averages are appended \
to the output file. The definition of aggregate size is quite flexible and \
the calculation can also be made only for aggregate sizes in a given \
range.\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <in.agg> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<in.agg>            input agg file\n");
  fprintf(ptr, "<output>            output file with per-timestep data\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --joined          specify that <input> contains joined "
          "coordinates\n");
  fprintf(ptr, "  -bt               specify bead types used for calculation "
          "(default is all)\n");
  fprintf(ptr, "  -m <name(s)>      agg size means number of <name(s)> "
          "molecules in an aggregate\n");
  fprintf(ptr, "  -only <name(s)>   use only aggregates composed "
          "of specified molecule(s)\n");
  fprintf(ptr, "  -ps <file>        save per-size averages to a <file>\n");
  fprintf(ptr, "  -n <int> <int>    calculate for aggregate sizes in "
          "given range\n");
  fprintf(ptr, "  -st <int>         starting timestep for calculation "
          "(only affects per-size and overall averages)\n");
  fprintf(ptr, "  -e <end>          ending timestep for calculation "
          "(only affects per-size and overall averages)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool joined,        // --joined
       *mt_m,         // -m
       *mt_only,      // -only
       only,          // -only (is the option present?)
       *bt;           // -bt (number of types; list of the types)
  int range_As[2];    // -n
  char ps_file[LINE]; // -ps
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 6, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version", "--joined", "-bt", "-m", "-only", "-ps", "-n");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);
  // test if <in.agg> ends with '.agg'
  int ext = 1;
  char extension[2][EXTENSION];
  s_strcpy(extension[0], ".agg", EXTENSION);
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  } //}}}

  // <output> - filename with data during simulation run
  char output[LINE];
  s_strcpy(output, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt->joined = false; // joined coordinates supplied, so no need to join
  } else {
    opt->joined = true; // molecules need to be joined
  } //}}}
  if (!FileOption(argc, argv, "-ps", opt->ps_file)) {
    opt->ps_file[0] = '\0';
  }
  //}}}

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-n' option - range of aggregation numbers //{{{
  opt->range_As[0] = 1;
  opt->range_As[1] = Count->Molecule;
  TwoNumbersOption(argc, argv, "-n", opt->range_As, 'i');
  if (opt->range_As[0] > opt->range_As[1]) {
    SwapInt(&opt->range_As[0], &opt->range_As[1]);
  } //}}}
  // '-m' option //{{{
  opt->mt_m = calloc(Count->MoleculeType, sizeof *opt->mt_m);
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->mt_m, System)) {
    InitBoolArray(opt->mt_m, Count->MoleculeType, true);
  } //}}}
  // '-only' option //{{{
  opt->mt_only = calloc(Count->MoleculeType, sizeof *opt->mt_only);
  opt->only = TypeOption(argc, argv, "-only", 'm', true, opt->mt_only, System);
  if (!opt->only) {
    InitBoolArray(opt->mt_only, Count->MoleculeType, true);
  } //}}}
  // -bt option //{{{
  opt->bt = calloc(Count->BeadType, sizeof *opt->bt);
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt->bt, System)) {
    InitBoolArray(opt->bt, Count->BeadType, true);
  } //}}}

  // // TODO those ridiculous flags are everywhere! //{{{
  // // copy Use flag to Write (for '-x' option)
  // for (int i = 0; i < Count->MoleculeType; i++) {
  //   MoleculeType[i].Write = MoleculeType[i].Flag;
  // }
  // // count total number of chains in excluded aggs
  // long int exclude_count_chains = 0;
  // // count total number of excluded aggs
  // long int exclude_count_agg = 0; //}}}

  // write initial stuff to output file //{{{
  PrintByline(output, argc, argv);
  FILE *out = OpenFile(output, "a");
  // print legend line to output file
  count = 1;
  fprintf(out, "# column: ");
  fprintf(out, "(%d) timestep", count++);
  fprintf(out, ", (%d) <Rg>_n", count++);
  fprintf(out, ", (%d) <Rg>_w", count++);
  fprintf(out, ", (%d) <Rg>_z", count++);
  fprintf(out, ", (%d) <Rg^2>_n", count++);
  fprintf(out, ", (%d) <Rg^2>_w", count++);
  fprintf(out, ", (%d) <Rg^2>_z", count++);
  fprintf(out, ", (%d) <Anis>_n", count++);
  fprintf(out, ", (%d) <Acyl>_n", count++);
  fprintf(out, ", (%d) <Aspher>_n", count++);
  fprintf(out, ", (%d) <eigen[0]>_n", count++);
  fprintf(out, ", (%d) <eigen[1]>_n", count++);
  fprintf(out, ", (%d) <eigen[2]>_n", count++);
  putc('\n', out);
  fclose(out); //}}}

  // open input aggregate file and skip the first lines (Aggregate command & blank line) //{{{
  double distance = 1; // TODO: read from agg file
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO go for while(fgets()); treatment
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // TODO memory allocation... Aaargh!
  // allocate memory for sum of various things //{{{
  // numbers of aggregates of all possibe sizes (maximum size is Count.Molecule)
  int *agg_counts_sum = calloc(Count->Molecule, sizeof *agg_counts_sum);
  // total radius of gyration: [size][0] normal sum, [size][1] sum of Rg*mass, [size][2] Rg*mass^2
  double (*Rg_sum)[3] = calloc(Count->Molecule, sizeof *Rg_sum);
  // total square of radius of gyration: [size][0] normal sum, [size][1] sum of Rg^2*mass, [size][2] Rg^2*mass^2
  double (*sqrRg_sum)[3] = calloc(Count->Molecule, sizeof *sqrRg_sum);
  // relative shape anisotropy: only normal sum
  double *Anis_sum = calloc(Count->Molecule, sizeof *Anis_sum);
  // acylindricity: only normal sum
  double *Acyl_sum = calloc(Count->Molecule, sizeof *Acyl_sum);
  // asphericity: only normal sum
  double *Aspher_sum = calloc(Count->Molecule, sizeof *Aspher_sum);
  // gyration tensor eigenvalues
  double (*eigen_sum)[3] = calloc(Count->Molecule, sizeof *eigen_sum);
  // total mass of aggregates: [size][0] normal sum, [size][1] sum of squares
  long int (*mass_sum)[2] = calloc(Count->Molecule, sizeof *mass_sum);
  // number of molecule types in aggregates: [size][mol type] only normal sum
  int **molecules_sum = malloc(Count->Molecule*sizeof(int *));
  for (int i = 0; i < Count->Molecule; i++) {
    molecules_sum[i] = calloc(Count->MoleculeType,sizeof(int));
  } //}}}

  // main loop //{{{
  FILE *coor = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0, // count lines in the vcf file
      line_count_agg = 0; // count lines in the agg file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);

    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(in, coor, &System, &line_count) ||
          ReadAggregates(agg, input_agg, &System,
                          Aggregate, &line_count_agg) < 0) {
        count_coor--;
        break;
      }
      count_used++;
      if (!opt->joined) {
        WrapJoinCoordinates(&System, false, true);
        RemovePBCAggregates(distance, Aggregate, &System);
      }

      // TODO: allocation... Aaargh!
      // allocate arrays for the timestep //{{{
      int *agg_counts_step = calloc(Count->Molecule, sizeof *agg_counts_step);
      double (*Rg_step)[3] = calloc(Count->Molecule, sizeof *Rg_step);
      double (*sqrRg_step)[3] = calloc(Count->Molecule, sizeof *sqrRg_step);
      double *Anis_step = calloc(Count->Molecule, sizeof *Anis_step);
      double *Acyl_step = calloc(Count->Molecule, sizeof *Acyl_step);
      double *Aspher_step = calloc(Count->Molecule,sizeof *Aspher_step);
      double (*eigen_step)[3] = calloc(Count->Molecule, sizeof *eigen_step); //}}}

      // calculate shape descriptors //{{{
      double mass_step[2] = {0}; // total mass of aggregates in a step: [0] normal, [1] sum of squares
      for (int i = 0; i < Count->Aggregate; i++) {
        // test if aggregate 'i' should be used //{{{
        int size = 0;
        // agg size = number of molecules of type 'specific_moltype_for_size'
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
          if (opt->mt_m[mol_type]) {
            size++;
          }
        }
        int correct_size = size - 1; // all agg sizes are used - relic of old times
        // make calculations only if agg size is well defined and within given range
        if (size == 0 || size < opt->range_As[0] || size > opt->range_As[1]) {
          continue;
        } //}}}
        // if '-only' is used, use only aggregates composed the specified molecule(s) //{{{
        bool test = true;
        if (opt->only) {
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int id = Aggregate[i].Molecule[j];
            if (opt->mt_only[System.Molecule[id].Type] == 0) {
              test = false; // a molecule is not of the required type
              break;
            }
          }
          if (!test) { // should the rest of the for loop agg i be skipped?
            continue;
          }
        } //}}}
        // // if '-x' option is used, discount aggregates with only specified molecule type(s) //{{{
        // test = false;
        // for (int j = 0; j < size; j++) {
        //   int moltype = System.Molecule[Aggregate[i].Molecule[j]].Type;
        //   if (opt->mt_x[moltype]) {
        //     test = true; // a molecule that shouldn't be in agg 'i' is there
        //     break;
        //   }
        // }
        // if (!test) { // should the rest of the for loop agg i be skipped?
        //   exclude_count_chains += Aggregate[i].nMolecules;
        //   exclude_count_agg++;
        //   continue;
        // } //}}}

        if (correct_size != -1) {
          // copy bead ids to a separate array //{{{
          int *list = malloc(Aggregate[i].nBeads * sizeof *list);
          int n = 0;
          double agg_mass = 0;
          for (int j = 0; j < Aggregate[i].nBeads; j++) {
            int id = Aggregate[i].Bead[j];
            int btype = System.Bead[id].Type;
            if (opt->bt[btype]) {
              list[n] = id;
              n++;
              agg_mass += System.BeadType[System.Bead[id].Type].Mass;
            }
          } //}}}

          double eigen[3];
          Gyration(n, list, &System, eigen);
          free(list); // free array of bead ids for gyration calculation

          double Rgi = sqrt(eigen[0] + eigen[1] + eigen[2]);
          // agg masses
          mass_step[0] += agg_mass; // for this timestep
          mass_step[1] += Square(agg_mass); // for this timestep
          // radius of gyration
          Rg_step[correct_size][0] += Rgi; // for number avg
          Rg_step[correct_size][1] += Rgi * agg_mass; // for weight average
          Rg_step[correct_size][2] += Rgi * Square(agg_mass); // for z-average
          // squared radius of gyration
          sqrRg_step[correct_size][0] += Square(Rgi); // for number avg
          sqrRg_step[correct_size][1] += Square(Rgi) * agg_mass; // for weight average
          sqrRg_step[correct_size][2] += Square(Rgi) * Square(agg_mass); // for z-average
          // relative shape anisotropy
          double temp[2];
          temp[0] = Square(eigen[0]) + Square(eigen[1]) + Square(eigen[2]);
          temp[1] = Square(eigen[0] + eigen[1] + eigen[2]);
          Anis_step[correct_size] += 1.5 * temp[0] / temp[1] - 0.5;
          // acylindricity
          Acyl_step[correct_size] += eigen[1] - eigen[0];
          // asphericity
          Aspher_step[correct_size] += eigen[2] - 0.5 * (eigen[0] + eigen[1]);
          // gyration vector eigenvalues
          eigen_step[correct_size][0] += eigen[0];
          eigen_step[correct_size][1] += eigen[1];
          eigen_step[correct_size][2] += eigen[2];
          // aggregate count
          agg_counts_step[correct_size]++;

          // sum molecules and aggregates (if step between start and end)
          if (count >= opt->c.start && (opt->c.end == -1 || count <= opt->c.end)) {
            agg_counts_sum[correct_size]++;
            // aggregate mass
            mass_sum[correct_size][0] += agg_mass;
            mass_sum[correct_size][1] += Square(agg_mass);
            for (int j = 0; j < Aggregate[i].nMolecules; j++) {
              int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
              molecules_sum[correct_size][mol_type]++;
            }
          }
        }
      } //}}}

      // add values to sums //{{{
      if (count >= opt->c.start && (opt->c.end == -1 || count <= opt->c.end)) {
        for (int i = 0; i < Count->Molecule; i++) {
          Rg_sum[i][0] += Rg_step[i][0];
          Rg_sum[i][1] += Rg_step[i][1];
          Rg_sum[i][2] += Rg_step[i][2];
          sqrRg_sum[i][0] += sqrRg_step[i][0];
          sqrRg_sum[i][1] += sqrRg_step[i][1];
          sqrRg_sum[i][2] += sqrRg_step[i][2];
          Anis_sum[i] += Anis_step[i];
          Acyl_sum[i] += Acyl_step[i];
          Aspher_sum[i] += Aspher_step[i];
          eigen_sum[i][0] += eigen_step[i][0];
          eigen_sum[i][1] += eigen_step[i][1];
          eigen_sum[i][2] += eigen_step[i][2];
        }
      } //}}}

      // print data to output file //{{{
      // sum up contributions from all aggregate sizes
      for (int i = 1; i < Count->Molecule; i++) {
        Rg_step[0][0] += Rg_step[i][0];
        Rg_step[0][1] += Rg_step[i][1];
        Rg_step[0][2] += Rg_step[i][2];
        sqrRg_step[0][0] += sqrRg_step[i][0];
        sqrRg_step[0][1] += sqrRg_step[i][1];
        sqrRg_step[0][2] += sqrRg_step[i][2];
        Anis_step[0] += Anis_step[i];
        Acyl_step[0] += Acyl_step[i];
        Aspher_step[0] += Aspher_step[i];
        eigen_step[0][0] += eigen_step[i][0];
        eigen_step[0][1] += eigen_step[i][1];
        eigen_step[0][2] += eigen_step[i][2];

        agg_counts_step[0] += agg_counts_step[i];
      }
      if (agg_counts_step[0] > 0) {
        out = OpenFile(output, "a");
        fprintf(out, "%5d", count_used); // timestep
        // <R_G>
        fprintf(out, " %lf %lf %lf", Rg_step[0][0]/agg_counts_step[0],
                                     Rg_step[0][1]/mass_step[0],
                                     Rg_step[0][2]/mass_step[1]);
        // <R_G^2>
        fprintf(out, " %lf %lf %lf", sqrRg_step[0][0]/agg_counts_step[0],
                                     sqrRg_step[0][1]/mass_step[0],
                                     sqrRg_step[0][2]/mass_step[1]);
        // relative shape anisotropy
        fprintf(out, " %lf", Anis_step[0]/agg_counts_step[0]);
        // acylindricity
        fprintf(out, " %lf", Acyl_step[0]/agg_counts_step[0]);
        // asphericity
        fprintf(out, " %lf", Aspher_step[0]/agg_counts_step[0]);
        // eigenvalues
        fprintf(out, " %lf %lf %lf", eigen_step[0][0]/agg_counts_step[0],
                                     eigen_step[0][1]/agg_counts_step[0],
                                     eigen_step[0][2]/agg_counts_step[0]);
        putc('\n', out);
        fclose(out);
      } //}}}

      // free memory //{{{
      free(agg_counts_step);
      free(Rg_step);
      free(sqrRg_step);
      free(Anis_step);
      free(Acyl_step);
      free(Aspher_step);
      free(eigen_step); //}}}
    } else { // skip the timestep - both coordinate and agg file //{{{
      if (!SkipTimestep(in, coor, &line_count) ||
          !SkipAggregates(agg, input_agg, &line_count_agg)) {
        count_coor--;
        break;
      }
    } //}}}

    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(coor);
  fclose(agg);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // calculate per-size averages? //{{{
  if (opt->ps_file[0] != '\0') {
    out = OpenFile(opt->ps_file, "w");
    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);

    fprintf(out, "# column: (1) agg size");
    for (int i = 0; i < Count->MoleculeType; i++) {
      fprintf(out, " (%d) <%s>", i+2, System.MoleculeType[i].Name);
    }
    fprintf(out, " (%d) <Rg>, ", Count->MoleculeType+2);
    fprintf(out, " (%d) <Rg^2>, ", Count->MoleculeType+3);
    fprintf(out, " (%d) <Anis>, ", Count->MoleculeType+4);
    fprintf(out, " (%d) <Acyl>, ", Count->MoleculeType+5);
    fprintf(out, " (%d) <Aspher>, ", Count->MoleculeType+6);
    fprintf(out, " (%d) <eigen[0]>, ", Count->MoleculeType+7);
    fprintf(out, " (%d) <eigen[1]>, ", Count->MoleculeType+8);
    fprintf(out, " (%d) <eigen[2]>, ", Count->MoleculeType+9);
    fprintf(out, " (%d) number of aggs", Count->MoleculeType+10);
    putc('\n', out);
    // determine width of each column & collate data //{{{
    int columns = Count->MoleculeType + 10;
    int digits[columns][2];
    InitInt2DArray((int *)digits, columns, 3, 0);
    // double data[columns][Count->Molecule];
    double *data[Count->Molecule];
    for (int i = 0; i < Count->Molecule; i++) {
      data[i] = calloc(columns, sizeof data[i]);
      if (agg_counts_sum[i] > 0) {
        count = -1;
        data[i][++count] = i + 1;
        for (int j = 0; j < Count->MoleculeType; j++) {
          data[i][++count] = (double)(molecules_sum[i][j]) / agg_counts_sum[i];
        }
        data[i][++count] = Rg_sum[i][0]/agg_counts_sum[i];
        data[i][++count] = sqrRg_sum[i][0]/agg_counts_sum[i];
        data[i][++count] = Anis_sum[i]/agg_counts_sum[i];
        data[i][++count] = Acyl_sum[i]/agg_counts_sum[i];
        data[i][++count] = Aspher_sum[i]/agg_counts_sum[i];
        data[i][++count] = eigen_sum[i][0]/agg_counts_sum[i];
        data[i][++count] = eigen_sum[i][1]/agg_counts_sum[i];
        data[i][++count] = eigen_sum[i][2]/agg_counts_sum[i];
        data[i][++count] = (double)(agg_counts_sum[i]);
      }
    }
    FillMaxDigits(columns, Count->Molecule, data, digits); //}}}

    // go over all possible sizes
    for (int i = 0; i < Count->Molecule; i++) {
      // is that size in the data?
      if (agg_counts_sum[i] > 0) {
        WriteFormatedDataLine(out, columns, data[i], digits);
      }
      free(data[i]);
    }

    fclose(out);
  } //}}}

  // total averages //{{{
  for (int i = 1; i < Count->Molecule; i++) {
    Rg_sum[0][0] += Rg_sum[i][0];
    Rg_sum[0][1] += Rg_sum[i][1];
    Rg_sum[0][2] += Rg_sum[i][2];
    sqrRg_sum[0][0] += sqrRg_sum[i][0];
    sqrRg_sum[0][1] += sqrRg_sum[i][1];
    sqrRg_sum[0][2] += sqrRg_sum[i][2];
    Anis_sum[0] += Anis_sum[i];
    Acyl_sum[0] += Acyl_sum[i];
    Aspher_sum[0] += Aspher_sum[i];
    eigen_sum[0][0] += eigen_sum[i][0];
    eigen_sum[0][1] += eigen_sum[i][1];
    eigen_sum[0][2] += eigen_sum[i][2];

    agg_counts_sum[0] += agg_counts_sum[i];

    mass_sum[0][0] += mass_sum[i][0];
    mass_sum[0][1] += mass_sum[i][1];

    for (int j = 0; j < Count->MoleculeType; j++) {
      molecules_sum[0][j] += molecules_sum[i][j];
    }
  }

  // print to output file
  out = OpenFile(output, "a");

  fprintf(out, "# (1) <M>_n, (2) <M>_w ");
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(out, "(%d) <%s>, ", i+3, System.MoleculeType[i].Name);
  }
  fprintf(out, "(%d) <Rg>_n, ", Count->MoleculeType+3);
  fprintf(out, "(%d) <Rg>_w, ", Count->MoleculeType+4);
  fprintf(out, "(%d) <Rg>_z, ", Count->MoleculeType+5);
  fprintf(out, "(%d) <Rg^2>_n, ", Count->MoleculeType+6);
  fprintf(out, "(%d) <Rg^2>_w, ", Count->MoleculeType+7);
  fprintf(out, "(%d) <Rg^2>_z, ", Count->MoleculeType+8);
  fprintf(out, "(%d) <Anis>, ", Count->MoleculeType+9);
  fprintf(out, "(%d) <Acyl>, ", Count->MoleculeType+10);
  fprintf(out, "(%d) <Aspher>, ", Count->MoleculeType+11);
  fprintf(out, "(%d) <eigen[0]>, ", Count->MoleculeType+12);
  fprintf(out, "(%d) <eigen[1]>, ", Count->MoleculeType+13);
  fprintf(out, "(%d) <eigen[2]>, ", Count->MoleculeType+14);
  putc('\n', out);
  fprintf(out, "# %lf", (double)(mass_sum[0][0])/agg_counts_sum[0]); //<M_As>_n
  fprintf(out, " %lf", (double)(mass_sum[0][1])/mass_sum[0][0]); //<M_As>_w
  // molecule types
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(out, " %lf", (double)(molecules_sum)[0][i]/agg_counts_sum[0]);
  }
  fprintf(out, " %lf", Rg_sum[0][0]/agg_counts_sum[0]); // <Rg>_n
  fprintf(out, " %lf", Rg_sum[0][1]/mass_sum[0][0]); // <Rg>_w
  fprintf(out, " %lf", Rg_sum[0][2]/mass_sum[0][1]); // <Rg>_z
  fprintf(out, " %lf", sqrRg_sum[0][0]/agg_counts_sum[0]); // <Rg^2>_n
  fprintf(out, " %lf", sqrRg_sum[0][1]/mass_sum[0][0]); // <Rg^2>_w
  fprintf(out, " %lf", sqrRg_sum[0][2]/mass_sum[0][1]); // <Rg^2>_z
  fprintf(out, " %lf", Anis_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %lf", Acyl_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %lf", Aspher_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0][0]/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0][1]/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0][2]/agg_counts_sum[0]);
  putc('\n', out);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  for (int i = 0; i < Count->Molecule; i++) {
    free(molecules_sum[i]);
  }
  free(molecules_sum);
  free(mass_sum);
  free(agg_counts_sum);
  free(Rg_sum);
  free(sqrRg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(eigen_sum);
  free(opt->mt_m);
  free(opt->mt_only);
  free(opt->bt);
  free(opt); //}}}

  return 0;
}
