#include "../src/AnalysisTools.h"
// TODO: two masses - -bt defined + always total (for contributions of given
//       subset of beads to the total gyration tensor)
// TODO: output printing

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "GyrationAggregates calculates gyration tensor for aggregates and determines "
  "shape descriptors like radius of gyration, acylindricity, asphericity, or "
  "relative shape anisotropy. By default, it calculates per-timestep averages, "
  "but per-size averages can also be determined. Overall averages are appended "
  "to the output file. The definition of aggregate size is quite flexible and "
  "the calculation can also be made only for aggregate sizes in a given "
  "range.",

  "Usage: GyrationAggregates <input> <in.agg> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 18, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<in.agg>", nullptr, "input agg file", OPT_ARG},
  {"<output>", nullptr, "output file with per-timestep data", OPT_ARG},
  {"--joined", nullptr, "<input> contains joined coordinates", OPT_EXTRA},
  {"-bt", nullptr, "bead types used for calculation (default: all)", OPT_EXTRA},
  {"-m", "<name(s)>", "aggregate size defined as number of specified molecules "
    "in an aggregate", OPT_EXTRA},
  {"-only", "<name(s)>", "use only aggregates composed of specified molecules",
    OPT_EXTRA},
  {"-n", "<int> <int>", "calculate for aggregate sizes in the given range",
    OPT_EXTRA},
  {"-ps", "<file>", "save per-size averages to a <file>", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  AGG_PICKER agg;     // -x, -only, -m, and -n arrays
  bool join,          // --joined
       *bt;           // -bt (number of types; list of the types)
  char ps_file[LINE]; // -ps
}; //}}}

// column indices for the shape-descriptor array (single 'normal sum' each)
enum { ANIS, ACYL, ASPH, NUM };

// all state shared between main() and the per-timestep callback //{{{
struct user_data {
  OPT opt;
  const char *output;
  double distance;      // -d from the agg file's Aggregates command
  bool *join_bt;        // -bt from the agg file's Aggregates command
  AGGREGATE *Aggregate;
  // overall sums
  int *agg_counts_sum;
  ArrNDd *Rg_sum, *sqrRg_sum, *shape_sum, *eigen_sum;
  ArrNDli *mass_sum;
  ArrNDi *molecules_sum;
  // per-timestep buffers (allocated once, zeroed each call)
  int *agg_counts_step;
  ArrNDd *Rg_step, *sqrRg_step, *shape_step, *eigen_step;
}; //}}}

// per-timestep calculation and output //{{{
static void Calculation(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *ud = userdata;
  const OPT *opt = &ud->opt;
  AGGREGATE *Aggregate = ud->Aggregate;
  COUNT *Count = &System->Count;

  if (opt->join) {
    RemovePBCAggregates(ud->distance, Aggregate, System, ud->join_bt);
  }

  // zero per-step buffers (reuse pre-allocated memory) //{{{
  int *agg_counts_step = ud->agg_counts_step;
  ArrNDd *Rg_step = ud->Rg_step;
  ArrNDd *sqrRg_step = ud->sqrRg_step;
  ArrNDd *shape_step = ud->shape_step;
  ArrNDd *eigen_step = ud->eigen_step;
  memset(agg_counts_step, 0, Count->Molecule * sizeof *agg_counts_step);
  FillArrND(Rg_step, 0.0);
  FillArrND(sqrRg_step, 0.0);
  FillArrND(shape_step, 0.0);
  FillArrND(eigen_step, 0.0); //}}}

  // calculate shape descriptors //{{{
  double mass_step[2] = {0}; // [0] normal agg mass, [1] sum of squares
  for (int i = 0; i < Count->Aggregate; i++) {
    // skip aggregates that shouldn't be used
    int agg_size;
    double rubbish; // unused
    if (!UseAggregate(*System, Aggregate, i, opt->agg, &agg_size, &rubbish) ||
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
      int btype = System->Bead[id].Type;
      if (opt->bt[btype]) {
        list[n] = id;
        n++;
        agg_mass += System->BeadType[System->Bead[id].Type].Mass;
      }
    } //}}}

    vec3d eigen = Gyration(n, list, System);
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
    double val = 1.5 * SqVectLength(eigen) /
                 Square(eigen.x + eigen.y + eigen.z) - 0.5;
    AddArr2D(shape_step, agg_size, ANIS, val);
    // acylindricity
    val = eigen.y - eigen.x;
    AddArr2D(shape_step, agg_size, ACYL, val);
    // asphericity
    val = eigen.z - 0.5 * (eigen.x + eigen.y);
    AddArr2D(shape_step, agg_size, ASPH, val);
    // gyration vector eigenvalues
    for (int dd = 0; dd < 3; dd++) {
      AddArr2D(eigen_step, agg_size, dd, eigen.v[dd]);
    }
    // aggregate count
    agg_counts_step[agg_size]++;

    // count molecules and aggregates
    ud->agg_counts_sum[agg_size]++;
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol_type = System->Molecule[AggGetMol(&Aggregate[i], j)].Type;
      AddArr2D(ud->molecules_sum, agg_size, mol_type, 1);
    }
    // sum aggregate mass
    AddArr2D(ud->mass_sum, agg_size, 0, agg_mass);
    AddArr2D(ud->mass_sum, agg_size, 1, Square(agg_mass));
  } //}}}

  for (int i = 0; i < Count->Molecule; i++) {
    for (int dd = 0; dd < 3; dd++) {
      AddArr2D(ud->Rg_sum, i, dd, GetArr2D(Rg_step, i, dd));
      AddArr2D(ud->sqrRg_sum, i, dd, GetArr2D(sqrRg_step, i, dd));
      AddArr2D(ud->eigen_sum, i, dd, GetArr2D(eigen_step, i, dd));
    }
    for (int dd = 0; dd < NUM; dd++) {
      AddArr2D(ud->shape_sum, i, dd, GetArr2D(shape_step, i, dd));
    }
  }

  // print data to output file //{{{
  // sum up contributions from all aggregate sizes
  for (int i = 1; i < Count->Molecule; i++) {
    for (int dd = 0; dd < 3; dd++) {
      AddArr2D(Rg_step, 0, dd, GetArr2D(Rg_step, i, dd));
      AddArr2D(sqrRg_step, 0, dd, GetArr2D(sqrRg_step, i, dd));
      AddArr2D(eigen_step, 0, dd, GetArr2D(eigen_step, i, dd));
    }
    for (int dd = 0; dd < NUM; dd++) {
      AddArr2D(shape_step, 0, dd, GetArr2D(shape_step, i, dd));
    }
    agg_counts_step[0] += agg_counts_step[i];
  }
  if (agg_counts_step[0] > 0) {
    FILE *out = OpenFile(ud->output, "a");
    fprintf(out, "%d", step->coor); // timestep
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
    fprintf(out, " %lf", GetArr2D(shape_step, 0, ANIS) / agg_counts_step[0]);
    // acylindricity
    fprintf(out, " %lf", GetArr2D(shape_step, 0, ACYL) / agg_counts_step[0]);
    // asphericity
    fprintf(out, " %lf", GetArr2D(shape_step, 0, ASPH) / agg_counts_step[0]);
    // eigenvalues
    val.v[0] = GetArr2D(eigen_step, 0, 0) / agg_counts_step[0];
    val.v[1] = GetArr2D(eigen_step, 0, 1) / agg_counts_step[0];
    val.v[2] = GetArr2D(eigen_step, 0, 2) / agg_counts_step[0];
    fprintf(out, " %lf %lf %lf", val.v[0], val.v[1], val.v[2]);
    putc('\n', out);
    fclose(out);
  } //}}}
} //}}}

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
  // --joined option (opt.join == true -> needs joining)
  opt.join = !BoolOption(argc, argv, "--joined");
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
  bool *join_bt = calloc(Count->BeadType, sizeof *join_bt);
  if (!join_bt) {
    ErrorAlloc("join_bt");
  }
  for (int i = 5; i < words; i++) {
    if (strcmp(split[i], "-bt") == 0) {
      for (i++; i < words && split[i][0] != '-'; i++) {
        int type = FindBeadType(split[i], System);
        if (type != -1) {
          join_bt[type] = true;
        }
      }
      break;
    }
  }
  // redefine distance if -d option is present
  double distance = 1;
  for (int i = 5; i < words; i++) {
    if (strcmp(split[i], "-d") == 0 && (i + 1) < words) {
      if (!IsPosRealNumber(split[i+1], &distance)) {
        distance = 1;
      }
      break;
    }
  }
  fclose(agg); // the main loop reopens the file itself
  //}}}

  AGGREGATE *Aggregate = nullptr;
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
  // shape descriptors [size][SH_ANIS|SH_ACYL|SH_ASPHER]: only normal sum
  ArrNDd *shape_sum = CreateArr2Dd(Count->Molecule, NUM);
  // gyration tensor eigenvalues
  ArrNDd *eigen_sum = CreateArr2Dd(Count->Molecule, 3);
  // total mass of aggregates: [size][0] normal sum, [size][1] sum of squares
  ArrNDli *mass_sum = CreateArr2Dli(Count->Molecule, 2);
  // number of molecule types in aggregates: [size][mol type] only normal sum
  ArrNDi *molecules_sum = CreateArr2Di(Count->Molecule, Count->MoleculeType);
  // per-timestep buffers (allocated once, zeroed each call in Calculation)
  int *agg_counts_step = calloc(Count->Molecule, sizeof *agg_counts_step);
  ArrNDd *Rg_step = CreateArr2Dd(Count->Molecule, 3);
  ArrNDd *sqrRg_step = CreateArr2Dd(Count->Molecule, 3);
  ArrNDd *shape_step = CreateArr2Dd(Count->Molecule, NUM);
  ArrNDd *eigen_step = CreateArr2Dd(Count->Molecule, 3);
  if (!agg_counts_sum || !Rg_sum || !sqrRg_sum || !shape_sum || !eigen_sum ||
      !mass_sum || !molecules_sum || !agg_counts_step || !Rg_step ||
      !sqrRg_step || !shape_step || !eigen_step) {
    ErrorAlloc("sum/step arrays");
  }
  //}}}

  // main loop //{{{
  struct user_data ud = {
    .opt = opt,
    .output = output,
    .distance = distance,
    .join_bt = join_bt,
    .Aggregate = Aggregate,
    .agg_counts_sum = agg_counts_sum,
    .Rg_sum = Rg_sum,
    .sqrRg_sum = sqrRg_sum,
    .shape_sum = shape_sum,
    .eigen_sum = eigen_sum,
    .mass_sum = mass_sum,
    .molecules_sum = molecules_sum,
    .agg_counts_step = agg_counts_step,
    .Rg_step = Rg_step,
    .sqrRg_step = sqrRg_step,
    .shape_step = shape_step,
    .eigen_step = eigen_step,
  };
  STEP step = InitStep;
  MainLoopCoorAgg(&System, in, input_agg, commons, &step, Aggregate,
                  Calculation, &ud); //}}}

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
        val = GetArr2D(shape_sum, i, ANIS) / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        val = GetArr2D(shape_sum, i, ACYL) / agg_counts_sum[i];
        SetArr2D(data, row_count, count++, val);
        val = GetArr2D(shape_sum, i, ASPH) / agg_counts_sum[i];
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
    for (int dd = 0; dd < NUM; dd++) {
      AddArr2D(shape_sum, 0, dd, GetArr2D(shape_sum, i, dd));
    }
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
  fprintf(out, " %lf", GetArr2D(shape_sum, 0, ANIS) / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(shape_sum, 0, ACYL) / agg_counts_sum[0]);
  fprintf(out, " %lf", GetArr2D(shape_sum, 0, ASPH) / agg_counts_sum[0]);
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
  FreeArrND(shape_sum);
  FreeArrND(eigen_sum);
  free(agg_counts_step);
  FreeArrND(Rg_step);
  FreeArrND(sqrRg_step);
  FreeArrND(shape_step);
  FreeArrND(eigen_step);
  free(opt.bt);
  free(join_bt);
  FreeAggPicker(&opt.agg);
  //}}}

  return 0;
}
