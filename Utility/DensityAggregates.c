#include "../src/AnalysisTools.h"
// TODO: -m_id option?

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "DensityAggregates utility calculates radial bead density for aggregates "
  "of given size(s) from their centre of mass. For beads in molecules, "
  "it takes into account only beads from the current aggregate, not from "
  "any other aggregate. The utility does not check what molecule type a "
  "given bead belongs to; therefore, if the same bead type appears in "
  "more molecule types, the resulting density for that bead type will be "
  "averaged without regard for the various types of molecule it comes from "
  "(i.e., if 'mol1' and 'mol2' both both contain bead 'A', there will be "
  "only one column for 'A' bead type).",

  "Usage: DensityAggregates <input> <in.agg> <width> <output> <size(s)> "
  "[options]",
  .args = 5, // number of mandatory arguments
  .all = 19, // number of valid lines OptSpec (not counting last {NULL})
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
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<in.agg>", NULL, "input agg file", OPT_ARG},
  {"<width>", NULL, "width of a single bin", OPT_ARG},
  {"<output>", NULL, "output density file (one per size; automatic '<size>.rho' ending)", OPT_ARG},
  {"<size(s)>", NULL, "aggregate sizes for density calculation", OPT_ARG},
  {"--joined", NULL, "specify that <input> contains joined coordinates", OPT_EXTRA},
  {"-m", "<name(s)>", "use number of specified molecule type(s) as aggrete size", OPT_EXTRA},
  {"-x", "<name(s)>", "exclude aggregates containing only specified molecule(s)", OPT_EXTRA},
  {"-only", "<name(s)>", "use only aggregates composed of specified molecule type(s)", OPT_EXTRA},
  {"-n", "<int> <int>", "calculate for aggregate sizes in given range", OPT_EXTRA},
  // {"-m_id", "<int>", "calculate only for aggregate containing the <int> molecule (by resid numbering in vsf)", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  AGG_PICKER agg; // -x, -only, -m, and -n arrays
  bool join;      // --joined
  FILE_TYPE fout; // -o
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
  // <width> - width of a single bin //{{{
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <output> - filename with bead densities
  char fout[LINE];
  s_strcpy(fout, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt.join = false; // joined coordinates supplied, so no need to join
  } else {
    opt.join = true; // molecules need to be joined
  } //}}}
  //}}}

  // print command to stdout
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  vec3d box = System.Box.OrthoLength;

  AggPickerOptions(argc, argv, &opt.agg, System);

  // <agg sizes> - aggregate sizes for calculation //{{{
  ArrNDi *agg_sizes = CreateArr2Di(Count->Molecule, 2);
  ArrNDi *agg_mols = CreateArr2Di(Count->Molecule, Count->MoleculeType);
  if (!agg_sizes || !agg_mols) {
    ErrorAlloc("agg_sizes/agg_mols");
  }
  int aggs = 0;
  while (++count < argc && argv[count][0] != '-') {
    // Error - non-numeric argument //{{{
    long val;
    if (!IsNaturalNumber(argv[count], &val)) {
      ErrorNaN("<agg size(s)>");
      Help(true, HelpDesc, opts);
      exit(1);
    } //}}}
    SetArr2D(agg_sizes, aggs, 0, atoi(argv[count]));
    aggs++; // number of aggregate sizes
    // ensure output string isn't too long for attaching <size>.vcf
    if (GetArr2D(agg_sizes, aggs, 0) < 10) {
      fout[LINE-1-4] = '\0';
    } else if (GetArr2D(agg_sizes, aggs, 0) < 100) {
      fout[LINE-2-4] = '\0';
    } else if (GetArr2D(agg_sizes, aggs, 0) < 1000) {
      fout[LINE-3-4] = '\0';
    } else if (GetArr2D(agg_sizes, aggs, 0) < 10000) {
      fout[LINE-4-4] = '\0';
    } else {
      fout[LINE-100] = '\0';
    }
  } //}}}

  double distance = 1; // <distance> parameter from Aggregate command

  // open input aggregate file and skip the first two lines
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO: go for while(fgetc()!='\n'); treatment
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
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], System);
    if (type != -1) {
      join_bt[type] = true;
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

  // number of bins
  double max_dist = 0.5 * Max3(box.x, box.y, box.z);
  int bins = ceil(max_dist / width);

  // allocate memory for density arrays
  ArrNDd *rho = CreateArr3Dd(Count->BeadType, aggs, bins);
  ArrNDd *rho_2 = CreateArr3Dd(Count->BeadType, aggs, bins);
  if (!rho || !rho_2) {
    ErrorAlloc("rho/rho_2");
  }

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_step = 0,
      count_used = 0,
      line_count = 0,
      line_count_agg = 0; // count lines in the agg file
  while (true) {
    PrintStep(&count_step, commons.start, commons.silent);

    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, count_step)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count) ||
          ReadAggregates(agg, input_agg, &System,
                          Aggregate, &line_count_agg) < 0) {
        count_step--;
        break;
      }
      count_used++;
      if (opt.join) {
        RemovePBCAggregates(distance, Aggregate, &System, join_bt);
      }
      ArrNDd *rho_temp = CreateArr3Dd(Count->BeadType, aggs, bins);
      if (!rho_temp) {
        ErrorAlloc("rho_temp");
      }

      // calculate densities //{{{
      for (int i = 0; i < Count->Aggregate; i++) {
        // // TODO: -m_id option? //{{{
        // // if '-m_id' is used, check if specified resid from vsf is in aggregate
        // if (m_id != -1) {
        //   for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        //     int id = Aggregate[i].Molecule[j];
        //     if (m_id == (id+1)) { // resname in vsf start from 1
        //       correct_size = 0;
        //       break;
        //     }
        //   }
        // } //}}}


        int agg_size;
        double agg_mass;
        if (!UseAggregate(System, Aggregate, i, opt.agg,
                          &agg_size, &agg_mass)) {
          continue;
        }
        // is agg_size in provided list?
        int correct_size = -1;
        for (int j = 0; j < aggs; j++) {
          if (GetArr2D(agg_sizes, j, 0) == agg_size) {
            correct_size = j;
          }
        }
        if (correct_size == -1) {
          continue;
        }

        vec3d com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead,
                                 System);

        FillArrND(rho_temp, 0);

        // aggregate beads //{{{
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          vec3d dist = DistancePBC(System.Bead[id].Position, com, &System.Box);
          dist.v[0] = VectLength(dist);

          if (dist.v[0] < max_dist) {
            int k = dist.v[0] / width;
            AddArr3D(rho_temp, System.Bead[id].Type, correct_size, k, 1);
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Count->Unbonded; j++) {
          int id = System.Unbonded[j];
          vec3d dist = Distance(System.Bead[id].Position, com, box);
          dist.v[0] = VectLength(dist);

          if (dist.v[0] < max_dist) {
            int k = dist.v[0] / width;
            // temp_rho[System.Bead[id].Type][correct_size][k]++;
            AddArr3D(rho_temp, System.Bead[id].Type, correct_size, k, 1);
          }
        } //}}}

        AddArr2D(agg_sizes, correct_size, 1, 1);

        // add from temporary density array to global density arrays //{{{
        for (int j = 0; j < Count->BeadType; j++) {
          for (int k = 0; k < bins; k++) {
            int val = GetArr3D(rho_temp, j, correct_size, k);
            AddArr3D(rho, j, correct_size, k, val);
            AddArr3D(rho_2, j, correct_size, k, Square(val));
          }
        } //}}}

        // count mol types in the aggregate
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol = AggGetMol(&Aggregate[i], j);
          // agg_mols[correct_size][System.Molecule[mol].Type]++;
          AddArr2D(agg_mols, correct_size, System.Molecule[mol].Type, 1);
        }
      } //}}}

      FreeArrND(rho_temp);
    //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count) ||
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
  fclose(fr);
  fclose(agg);
  PrintLastStep(count_step, count_used, commons.silent); //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < aggs; i++) {
    FILE *out;
    char str[LINE];
    // assemble correct name
    if (snprintf(str, LINE, "%s%d.rho", fout, GetArr2D(agg_sizes, i, 0)) < 0) {
      ErrorSnprintf();
    }
    // write initial stuff to the file //{{{
    out = PrintBylineOpenFile(str, argc, argv);
    // print bead type names to output file
    count = 1;
    fprintf(out, "# 2 radial profiles (columns) per bead type: ");
    fprintf(out, "(%d) number profile", count++);
    fprintf(out, ", (%d) density profile", count++);
    putc('\n', out);
    count = 1;
    fprintf(out, "# Columns: (%d) distance;", count++);
    for (int j = 0; j < Count->BeadType; j++) {
      fprintf(out, " (%d) %s", count, System.BeadType[j].Name);
      // 2 columns per bead type
      count += 2;
      if (j != (Count->BeadType-1)) {
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    int ncols = 2 * Count->BeadType + 1;
    int nrows = bins;
    ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
    if (!data) {
      err_msg("ArrNDd constructor failed (data)");
      PrintError();
      exit(1);
    }

    // calculate rdf
    for (int j = 0; j < nrows; j++) {
      double shell = 4 * PI * Cube(width) * (Cube(j + 1) - Cube(j)) / 3;
      // fprintf(out, "%.2f", width * (2 * j + 1) / 2);
      count = 0;
      SetArr2D(data, j, count++, width * (2 * j + 1) / 2);

      for (int k = 0; k < Count->BeadType; k++) {
        // radial number profile
        double val = GetArr3D(rho, k, i, j) / GetArr2D(agg_sizes, i, 1);
        SetArr2D(data, j, count++, val);
        // radial density profile
        val /= shell;
        SetArr2D(data, j, count++, val);
      }
    }
    ComputeColumnWidths(nrows, ncols, data, 6);
    PrintDataAll(out, nrows, ncols, data);
    FreeArrND(data);

    fprintf(out, "# (1) Number of aggregates; Average numbers of molecules:");
    count = 2;
    for (int j = 0; j < Count->MoleculeType; j++) {
      fprintf(out," (%d) %s", count++, System.MoleculeType[j].Name);
      if (i < (Count->MoleculeType - 1)) {
        putc(';', out);
      }
    }
    putc('\n', out);
    putc('#', out);
    fprintf(out," %d", GetArr2D(agg_sizes, i, 1));
    for (int j = 0; j < Count->MoleculeType; j++) {
      double val = (double)(GetArr2D(agg_mols, i, j)) /
                   GetArr2D(agg_sizes, i, 1);
      fprintf(out," %lf", val);
    }
    putc('\n', out);

    fclose(out);
  } //}}}

  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  FreeAggPicker(&opt.agg);
  FreeArrND(agg_sizes);
  FreeArrND(agg_mols);
  FreeArrND(rho);
  FreeArrND(rho_2);
  free(join_bt);
  return 0;
}
