#include "../src/AnalysisTools.h"
#include <stdbool.h>
// TODO: --joined function when RemovePBCAggregates() works
// TODO: --joined --> --join; make default expectation of joined coordinates
// TODO: two masses - -bt defined + always total (for contributions of given
//       subset of beads to the total gyration tensor)
// TODO: output printing
// TODO: arrays

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "ASSUMES JOINED COORDINATES!\n"
  "GyrationAggregates calculates gyration tensor for aggregates and determines "
  "shape descriptors like radius of gyration, acylindricity, asphericity, or "
  "relative shape anisotropy. By default, it calculates per-timestep averages, "
  "but per-size averages can also be determined. Overall averages are appended "
  "to the output file. The definition of aggregate size is quite flexible and "
  "the calculation can also be made only for aggregate sizes in a given "
  "range.",

  "Usage: GyrationAggregates <input> <in.agg> <output> [options]",
  .args = 3,
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", NULL, "input coordinate file"},
  {"<in.agg>", NULL, "input agg file"},
  {"<output>", NULL, "output file with per-timestep data"},
  // {"--joined", NULL, "<input> contains joined coordinates"},
  {"-bt", NULL, "bead types used for calculation (default: all)"},
  {"-m", "<name(s)>", "agg size defined as number of <name(s)> molecules in an aggregate"},
  {"-only", "<name(s)>", "use only aggregates composed of specified molecule(s)"},
  {"-n", "<int> <int>", "calculate for aggregate sizes in given range"},
  {"-ps", "<file>", "save per-size averages to a <file>"},
  {NULL}
}; //}}}

// Help() //{{{
void Help_old(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  AGG_PICKER agg;     // -x, -only, -m, and -n arrays
  bool join,          // --joined
       *bt;           // -bt (number of types; list of the types)
  char ps_file[LINE]; // -ps
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
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
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <output> - filename with data during simulation run
  char output[LINE];
  s_strcpy(output, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt.join = false; // joined coordinates supplied, so no need to join
  } else {
    opt.join = true; // molecules need to be joined
  } //}}}
  if (!FileOption(argc, argv, "-ps", opt.ps_file)) {
    opt.ps_file[0] = '\0';
  }
  //}}}

  // print command to stdout
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  AggPickerOptions(argc, argv, &opt.agg, System);

  // -bt option //{{{
  opt.bt = calloc(Count->BeadType, sizeof *opt.bt);
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System)) {
    InitBoolArray(opt.bt, Count->BeadType, true);
  } //}}}

  // // TODO: those ridiculous flags are everywhere! //{{{
  // // copy Use flag to Write (for '-x' option)
  // for (int i = 0; i < Count->MoleculeType; i++) {
  //   MoleculeType[i].Write = MoleculeType[i].Flag;
  // }
  // // count total number of chains in excluded aggs
  // long int exclude_count_chains = 0;
  // // count total number of excluded aggs
  // long int exclude_count_agg = 0; //}}}

  // write initial stuff to output file //{{{
  FILE *out = PrintBylineOpenFile(output, argc, argv);
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
  // TODO: make into function (other *Aggregates utils need it)
  if (!ReadAndSplitLine(agg, SPL_STR, " \t\n")) {
    if (snprintf(ERROR_MSG, LINE, "empty %s%s%s line",
                 ErrYellow(), input_agg, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  }
  // bead types for connecting aggregates
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], System);
    if (type != -1) { // TODO: don't use Flag (RemovePBCAggregates function)
      System.BeadType[type].Flag = true;
    }
  }
  // redefine distance if -d option is present
  for (int i = 5; i < words; i++) {
    if (strcmp(split[i], "-d") == 0 && (i + 1) < words) {
      if (!IsPosRealNumber(split[i+1], &distance)) {
        distance = 1;
      }
      break;
    }
  }
  //}}}

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
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
  int count_step = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0, // count lines in the vcf file
      line_count_agg = 0; // count lines in the agg file
  while (true) {
    PrintStep(&count_step, commons.start, commons.silent);

    bool use = false;
    if (UseStep(commons, count_step)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, coor, &System, &line_count) ||
          ReadAggregates(agg, input_agg, &System,
                         Aggregate, &line_count_agg) < 0) {
        count_step--;
        break;
      }
      count_used++;
      // TODO: return when it works
      // if (opt.join) {
      //   RemovePBCAggregates(distance, Aggregate, &System);
      // }

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
        // skip aggregates that shouldn't be used
        int agg_size;
        double rubbish; // unused
        if (!UseAggregate(System, Aggregate, i, opt.agg, &agg_size, &rubbish) ||
            Aggregate[i].nBeads == 1) { // cannot do gyration for 1 point
          continue;
        }

        // copy bead ids to a separate array //{{{
        int *list = malloc(Aggregate[i].nBeads * sizeof *list);
        int n = 0;
        double agg_mass = 0;
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          int btype = System.Bead[id].Type;
          if (opt.bt[btype]) {
            list[n] = id;
            n++;
            agg_mass += System.BeadType[System.Bead[id].Type].Mass;
          }
        } //}}}

        vec3d eigen = Gyration(n, list, &System);
        free(list); // free array of bead ids for gyration calculation
        // skip case of no size (particles fully collapsed onto each other)
        if (eigen.x == 0 && eigen.y == 0 && eigen.z == 0) {
          continue;
        }

        double Rgi = sqrt(eigen.x + eigen.y + eigen.z);
        // agg masses
        mass_step[0] += agg_mass; // for this timestep
        mass_step[1] += Square(agg_mass); // for this timestep
        // radius of gyration
        Rg_step[agg_size][0] += Rgi; // for number avg
        Rg_step[agg_size][1] += Rgi * agg_mass; // for weight average
        Rg_step[agg_size][2] += Rgi * Square(agg_mass); // for z-average
        // squared radius of gyration
        sqrRg_step[agg_size][0] += Square(Rgi); // for number avg
        sqrRg_step[agg_size][1] += Square(Rgi) * agg_mass; // for weight average
        sqrRg_step[agg_size][2] += Square(Rgi) * Square(agg_mass); // for z-average
        // relative shape anisotropy
        double temp[2];
        temp[0] = SqVectLength(eigen);
        temp[1] = Square(eigen.x + eigen.y + eigen.z);
        Anis_step[agg_size] += 1.5 * temp[0] / temp[1] - 0.5;
        // acylindricity
        Acyl_step[agg_size] += eigen.y - eigen.x;
        // asphericity
        Aspher_step[agg_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
        // gyration vector eigenvalues
        for (int dd = 0; dd < 3; dd++) {
          eigen_step[agg_size][dd] += eigen.v[dd];
        }
        // aggregate count
        agg_counts_step[agg_size]++;

        // count molecules and aggregates
        agg_counts_sum[agg_size]++;
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_sum[agg_size][mol_type]++;
        }
        // sum aggregate mass
        mass_sum[agg_size][0] += agg_mass;
        mass_sum[agg_size][1] += Square(agg_mass);
      } //}}}

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
        fprintf(out, "%5d", count_step); // timestep
        // <R_G>
        fprintf(out, " %14f %14f %14f", Rg_step[0][0]/agg_counts_step[0],
                                        Rg_step[0][1]/mass_step[0],
                                        Rg_step[0][2]/mass_step[1]);
        // <R_G^2>
        fprintf(out, " %14f %14f %14f", sqrRg_step[0][0]/agg_counts_step[0],
                                        sqrRg_step[0][1]/mass_step[0],
                                        sqrRg_step[0][2]/mass_step[1]);
        // relative shape anisotropy
        fprintf(out, " %14f", Anis_step[0]/agg_counts_step[0]);
        // acylindricity
        fprintf(out, " %14f", Acyl_step[0]/agg_counts_step[0]);
        // asphericity
        fprintf(out, " %14f", Aspher_step[0]/agg_counts_step[0]);
        // eigenvalues
        fprintf(out, " %14f %14f %14f", eigen_step[0][0]/agg_counts_step[0],
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
    //}}}
    } else {
      if (!SkipTimestep(in, coor, &line_count) ||
          !SkipAggregates(agg, input_agg, &line_count_agg)) {
        count_step--;
        break;
      }
    }

    // exit the main loop if reached user-specied end timestep
    if (count_step == commons.end) {
      break;
    }
  }
  fclose(coor);
  fclose(agg);
  PrintLastStep(count_step, count_used, commons.silent); //}}}

  // calculate per-size averages? //{{{
  if (opt.ps_file[0] != '\0') {
    out = OpenFile(opt.ps_file, "w");
    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);

    count = 1;
    fprintf(out, "# column: (%d) agg size", count++);
    fprintf(out, " (%d) <Rg>, ", count++);
    fprintf(out, " (%d) <Rg^2>, ", count++);
    fprintf(out, " (%d) <Anis>, ", count++);
    fprintf(out, " (%d) <Acyl>, ", count++);
    fprintf(out, " (%d) <Aspher>, ", count++);
    fprintf(out, " (%d) <eigen[0]>, ", count++);
    fprintf(out, " (%d) <eigen[1]>, ", count++);
    fprintf(out, " (%d) <eigen[2]>, ", count++);
    fprintf(out, " (%d) number of aggs", count++);
    for (int i = 0; i < Count->MoleculeType; i++) {
      fprintf(out, " (%d) <%s>_n", count++, System.MoleculeType[i].Name);
    }
    putc('\n', out);
    // determine width of each column & collate data //{{{
    int columns = Count->MoleculeType + 10;
    int digits[columns][2];
    InitInt2DArray((int *)digits, columns, 2, 0);
    // double data[columns][Count->Molecule];
    double *data[Count->Molecule];
    for (int i = 0; i < Count->Molecule; i++) {
      data[i] = calloc(columns, sizeof data[i]);
      if (agg_counts_sum[i] > 0) {
        count = -1;
        data[i][++count] = i + 1;
        data[i][++count] = Rg_sum[i][0]/agg_counts_sum[i];
        data[i][++count] = sqrRg_sum[i][0]/agg_counts_sum[i];
        data[i][++count] = Anis_sum[i]/agg_counts_sum[i];
        data[i][++count] = Acyl_sum[i]/agg_counts_sum[i];
        data[i][++count] = Aspher_sum[i]/agg_counts_sum[i];
        data[i][++count] = eigen_sum[i][0]/agg_counts_sum[i];
        data[i][++count] = eigen_sum[i][1]/agg_counts_sum[i];
        data[i][++count] = eigen_sum[i][2]/agg_counts_sum[i];
        data[i][++count] = (double)(agg_counts_sum[i]);
        for (int j = 0; j < Count->MoleculeType; j++) {
          data[i][++count] = (double)(molecules_sum[i][j]) / agg_counts_sum[i];
        }
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

  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(out, "(%d) <%s>, ", count++, System.MoleculeType[i].Name);
  }
  fprintf(out, "(%d) <Rg>_n, ", count++);
  fprintf(out, "(%d) <Rg>_w, ", count++);
  fprintf(out, "(%d) <Rg>_z, ", count++);
  fprintf(out, "(%d) <Rg^2>_n, ", count++);
  fprintf(out, "(%d) <Rg^2>_w, ", count++);
  fprintf(out, "(%d) <Rg^2>_z, ", count++);
  fprintf(out, "(%d) <Anis>, ", count++);
  fprintf(out, "(%d) <Acyl>, ", count++);
  fprintf(out, "(%d) <Aspher>, ", count++);
  fprintf(out, "(%d) <eigen.x>, ", count++);
  fprintf(out, "(%d) <eigen.y>, ", count++);
  fprintf(out, "(%d) <eigen.z>, ", count++);
  putc('\n', out);

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
  free(opt.bt);
  FreeAggPicker(&opt.agg);
  //}}}

  return 0;
}
