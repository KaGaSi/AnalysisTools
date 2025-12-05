#include "../src/AnalysisTools.h"
#include <stdbool.h>
// TODO: --joined function when RemovePBCAggregates() works
// TODO: --joined --> --join; make default expectation of joined coordinates
// TODO: two masses - -bt defined + always total (for contributions of given
//       subset of beads to the total gyration tensor)
// TODO: output printing

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
  if (!(opt.bt = calloc(Count->BeadType, sizeof *opt.bt))) {
    ErrorAlloc("op.bt");
  }
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
  fprintf(out, ", (%d) <eigen.x>_n", count++);
  fprintf(out, ", (%d) <eigen.y>_n", count++);
  fprintf(out, ", (%d) <eigen.z>_n", count++);
  putc('\n', out);
  fclose(out); //}}}

  // open input aggregate file and skip the first lines (Aggregate command & blank line) //{{{
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
  double distance;
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

  // allocate memory for sum of various things //{{{
  // numbers of aggregates of all possibe sizes (maximum size is Count.Molecule)
  int *agg_counts_sum = calloc(Count->Molecule, sizeof *agg_counts_sum);
  // total radius of gyration: [size][0] normal sum, [size][1] sum of Rg*mass, [size][2] Rg*mass^2
  ArrNDd *Rg_sum = CreateArr2Dd(Count->Molecule, 3);
  // total square of radius of gyration
  // ...[size][0] normal sum, [size][1] sum of Rg^2*mass, [size][2] Rg^2*mass^2
  ArrNDd *sqrRg_sum = CreateArr2Dd(Count->Molecule, 3);
  // relative shape anisotropy: only normal sum
  double *Anis_sum = calloc(Count->Molecule, sizeof *Anis_sum);
  // acylindricity: only normal sum
  double *Acyl_sum = calloc(Count->Molecule, sizeof *Acyl_sum);
  // asphericity: only normal sum
  double *Aspher_sum = calloc(Count->Molecule, sizeof *Aspher_sum);
  // gyration tensor eigenvalues
  ArrNDd *eigen_sum = CreateArr2Dd(Count->Molecule, 3);
  // total mass of aggregates: [size][0] normal sum, [size][1] sum of squares
  ArrNDli *mass_sum = CreateArr2Dli(Count->Molecule, 2);
  // number of molecule types in aggregates: [size][mol type] only normal sum
  ArrNDi *molecules_sum = CreateArr2Di(Count->Molecule, Count->MoleculeType);
  if (!agg_counts_sum || !Rg_sum || !sqrRg_sum || !Anis_sum || !Acyl_sum ||
      !Aspher_sum || !eigen_sum || !mass_sum || !molecules_sum) {
    ErrorAlloc("molecules_sum");
  }
  //}}}

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

      // allocate arrays for the timestep //{{{
      int *agg_counts_step = calloc(Count->Molecule, sizeof *agg_counts_step);
      ArrNDd *Rg_step = CreateArr2Dd(Count->Molecule, 3);
      ArrNDd *sqrRg_step = CreateArr2Dd(Count->Molecule, 3);
      double *Anis_step = calloc(Count->Molecule, sizeof *Anis_step);
      double *Acyl_step = calloc(Count->Molecule, sizeof *Acyl_step);
      double *Aspher_step = calloc(Count->Molecule,sizeof *Aspher_step);
      ArrNDd *eigen_step = CreateArr2Dd(Count->Molecule, 3);
      if (!agg_counts_step || !Rg_step || !sqrRg_step || !Anis_step ||
          !Acyl_step || !Aspher_step || !eigen_step) {
        ErrorAlloc("step arrays");
      } //}}}

      // calculate shape descriptors //{{{
      double mass_step[2] = {0}; // [0] normal agg mass, [1] sum of squares
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
        if (!list) {
          ErrorAlloc("list");
        }
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
        AddArr2D(Rg_step, agg_size, 0, Rgi);
        AddArr2D(Rg_step, agg_size, 1, Rgi * agg_mass);
        AddArr2D(Rg_step, agg_size, 2, Rgi * Square(agg_mass));
        // squared radius of gyration
        AddArr2D(sqrRg_step, agg_size, 0, Square(Rgi));
        AddArr2D(sqrRg_step, agg_size, 1, Square(Rgi) * agg_mass);
        AddArr2D(sqrRg_step, agg_size, 2, Square(Rgi) * Square(agg_mass));
        // relative shape anisotropy
        Anis_step[agg_size] += 1.5 * SqVectLength(eigen) /
                               Square(eigen.x + eigen.y + eigen.z) - 0.5;
        // acylindricity
        Acyl_step[agg_size] += eigen.y - eigen.x;
        // asphericity
        Aspher_step[agg_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
        // gyration vector eigenvalues
        for (int dd = 0; dd < 3; dd++) {
          AddArr2D(eigen_step, agg_size, dd, eigen.v[dd]);
        }
        // aggregate count
        agg_counts_step[agg_size]++;

        // count molecules and aggregates
        agg_counts_sum[agg_size]++;
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
          AddArr2D(molecules_sum, agg_size, mol_type, 1);
        }
        // sum aggregate mass
        AddArr2D(mass_sum, agg_size, 0, agg_mass);
        AddArr2D(mass_sum, agg_size, 1, Square(agg_mass));
      } //}}}

      for (int i = 0; i < Count->Molecule; i++) {
        for (int dd = 0; dd < 3; dd++) {
          AddArr2D(Rg_sum, i, dd, GetArr2D(Rg_step, i, dd));
          AddArr2D(sqrRg_sum, i, dd, GetArr2D(sqrRg_step, i, dd));
          AddArr2D(eigen_sum, i, dd, GetArr2D(eigen_step, i, dd));
        }
        Anis_sum[i] += Anis_step[i];
        Acyl_sum[i] += Acyl_step[i];
        Aspher_sum[i] += Aspher_step[i];
      }

      // print data to output file //{{{
      // sum up contributions from all aggregate sizes
      for (int i = 1; i < Count->Molecule; i++) {
        for (int dd = 0; dd < 3; dd++) {
          AddArr2D(Rg_step, 0, dd, GetArr2D(Rg_step, i, dd));
          AddArr2D(sqrRg_step, 0, dd, GetArr2D(sqrRg_step, i, dd));
          AddArr2D(eigen_step, 0, dd, GetArr2D(eigen_step, i, dd));
        }
        Anis_step[0] += Anis_step[i];
        Acyl_step[0] += Acyl_step[i];
        Aspher_step[0] += Aspher_step[i];

        agg_counts_step[0] += agg_counts_step[i];
      }
      if (agg_counts_step[0] > 0) {
        out = OpenFile(output, "a");
        fprintf(out, "%d", count_step); // timestep
        // <R_G>
        vec3d val;
        val.v[0] = GetArr2D(Rg_step, 0, 0) / agg_counts_step[0];
        val.v[1] = GetArr2D(Rg_step, 0, 1) / mass_step[0];
        val.v[2] = GetArr2D(Rg_step, 0, 2) / mass_step[1];
        fprintf(out, " %lf %lf %lf", val.v[0], val.v[1], val.v[2]);
        // <R_G^2>
        val.v[0] = GetArr2D(sqrRg_step, 0, 0) / agg_counts_step[0];
        val.v[1] = GetArr2D(sqrRg_step, 0, 1) / mass_step[0];
        val.v[2] = GetArr2D(sqrRg_step, 0, 2) / mass_step[1];
        fprintf(out, " %lf %lf %lf", val.v[0], val.v[1], val.v[2]);
        // relative shape anisotropy
        fprintf(out, " %lf", Anis_step[0]/agg_counts_step[0]);
        // acylindricity
        fprintf(out, " %lf", Acyl_step[0]/agg_counts_step[0]);
        // asphericity
        fprintf(out, " %lf", Aspher_step[0]/agg_counts_step[0]);
        // eigenvalues
        val.v[0] = GetArr2D(eigen_step, 0, 0) / agg_counts_step[0];
        val.v[1] = GetArr2D(eigen_step, 0, 1) / agg_counts_step[0];
        val.v[2] = GetArr2D(eigen_step, 0, 2) / agg_counts_step[0];
        fprintf(out, " %lf %lf %lf", val.v[0], val.v[1], val.v[2]);
        putc('\n', out);
        fclose(out);
      } //}}}

      FreeArrND(Rg_step);
      FreeArrND(sqrRg_step);
      FreeArrND(eigen_step);
      free(agg_counts_step);
      free(Anis_step);
      free(Acyl_step);
      free(Aspher_step);
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

    // collate data //{{{
    int ncols = count - 1; // header print ends with count++, therefore - 1
    int nrows = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      if (agg_counts_sum[i] > 0) {
        nrows++;
      }
    }
    ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
    if (!data) {
      ErrorAlloc("data");
    }
    int row_count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      if (agg_counts_sum[i] > 0) {
        count = 0;
        SetArr2D(data, row_count, count++, i + 1);
        double val = GetArr2D(Rg_sum, i, 0) / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        val = GetArr2D(sqrRg_sum, i, 0) / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        val = Anis_sum[i] / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        val = Acyl_sum[i] / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        val = Aspher_sum[i] / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        for (int dd = 0; dd < 3; dd++) {
          val = GetArr2D(eigen_sum, i, dd) / agg_counts_sum[i];
          SetArr2D(data, row_count, count++, val);
        }
        SetArr2D(data, row_count, count++, agg_counts_sum[i]);
        for (int j = 0; j < Count->MoleculeType; j++) {
          val = (double)GetArr2D(molecules_sum, i, j) / agg_counts_sum[i];
          SetArr2D(data, row_count, count++, val);
        }
        row_count++;
      }
    } //}}}

    ComputeColumnWidths(nrows, ncols, data, 6);
    PrintDataAll(out, nrows, ncols, data);
    FreeArrND(data);

    fclose(out);
  } //}}}

  // total averages //{{{
  for (int i = 1; i < Count->Molecule; i++) {
    for (int dd = 0; dd < 3; dd++) {
      AddArr2D(Rg_sum, 0, dd, GetArr2D(Rg_sum, i, dd));
      AddArr2D(sqrRg_sum, 0, dd, GetArr2D(sqrRg_sum, i, dd));
      AddArr2D(eigen_sum, 0, dd, GetArr2D(eigen_sum, i, dd));
    }
    Anis_sum[0] += Anis_sum[i];
    Acyl_sum[0] += Acyl_sum[i];
    Aspher_sum[0] += Aspher_sum[i];

    agg_counts_sum[0] += agg_counts_sum[i];

    AddArr2D(mass_sum, 0, 0, GetArr2D(mass_sum, i, 0));
    AddArr2D(mass_sum, 0, 1, GetArr2D(mass_sum, i, 1));

    for (int j = 0; j < Count->MoleculeType; j++) {
      AddArr2D(molecules_sum, 0, j,  GetArr2D(molecules_sum, i, j));
    }
  }

  // print to output file
  out = OpenFile(output, "a");

  count = 1;
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
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(out, "(%d) <%s>", count++, System.MoleculeType[i].Name);
    if (i != (Count->MoleculeType - 1)) {
      fprintf(out, ", ");
    }
  }
  putc('\n', out);

  fprintf(out, " %lf", GetArr2D(Rg_sum, 0, 0) / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(Rg_sum, 0, 1) / GetArr2D(mass_sum, 0, 0));
  fprintf(out, " %lf", GetArr2D(Rg_sum, 0, 2) / GetArr2D(mass_sum, 0, 1));
  fprintf(out, " %lf", GetArr2D(sqrRg_sum, 0, 0) / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(sqrRg_sum, 0, 1) / GetArr2D(mass_sum, 0, 0));
  fprintf(out, " %lf", GetArr2D(sqrRg_sum, 0, 2) / GetArr2D(mass_sum, 0, 1));
  fprintf(out, " %lf", Anis_sum[0] / agg_counts_sum[0]);
  fprintf(out, " %lf", Acyl_sum[0] / agg_counts_sum[0]);
  fprintf(out, " %lf", Aspher_sum[0] / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(eigen_sum, 0, 0) / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(eigen_sum, 0, 1) / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(eigen_sum, 0, 2) / agg_counts_sum[0]);
  // molecule types
  for (int i = 0; i < Count->MoleculeType; i++) {
    fprintf(out, " %lf", (double)GetArr2D(molecules_sum, 0, i) /
                         agg_counts_sum[0]);
  }
  putc('\n', out);
  fclose(out); //}}}

  // free memory //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  FreeArrND(molecules_sum);
  FreeArrND(mass_sum);
  free(agg_counts_sum);
  FreeArrND(Rg_sum);
  FreeArrND(sqrRg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  FreeArrND(eigen_sum);
  free(opt.bt);
  FreeAggPicker(&opt.agg);
  //}}}

  return 0;
}
