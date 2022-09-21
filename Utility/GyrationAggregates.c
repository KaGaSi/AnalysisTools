#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
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

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <in.agg> <output> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <in.agg>   input agg file\n");
  fprintf(ptr, "   <output>   output file with per-timestep data\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined          specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -bt               specify bead types to be used for calculation (default is all)\n");
  fprintf(ptr, "      -m <name(s)>      agg size means number of <name(s)> \
molecules in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>      exclude aggregates containing only specified molecule(s)\n");
  fprintf(ptr, "      -only <name(s)>   use only aggregates composed \
of specified molecule(s)\n");
  fprintf(ptr, "      -ps <file>        save per-size averages to a <file>\n");
  fprintf(ptr, "      -n <int> <int>    calculate for aggregate sizes in \
given range\n");
  fprintf(ptr, "      -st <int>         starting timestep for calculation \
(only affects per-size and overall averages)\n");
  fprintf(ptr, "      -e <end>          ending timestep for calculation \
(only affects per-size and overall averages)\n");
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
  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }
  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-only") != 0 &&
        strcmp(argv[i], "-ps") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE] = "", input_vsf[LINE] = "";
  snprintf(input_coor, LINE, "%s", argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor_old(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  snprintf(input_agg, LINE, "%s", argv[++count]);
  // test if <in.agg> ends with '.agg'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - filename with data during simulation run
  char output[LINE];
  snprintf(output, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");
  // write per-agg averages to a file?
  char per_size_file[LINE] = "";
  if (FileOption(argc, argv, "-ps", per_size_file, LINE)) {
    exit(0);
  }
  int start, end;
  StartEndTime(argc, argv, &start, &end); //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // '-n' option - range of aggregation numbers //{{{
  int range_As[100], // TODO: 100 only until IntegerOption is changed
      test = 2;                   //
  range_As[0] = 1;                // by default, use all aggregation numbers
  range_As[1] = Counts.Molecules; //
  if (MultiIntegerOption(argc, argv, "-n", &test, range_As)) {
    exit(1);
  }
  if (test != 2) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-n");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " option requires two numeric arguments\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  // make sure first number is larger
  if (range_As[0] > range_As[1]) {
    SwapInt(&range_As[0], &range_As[1]);
  } //}}}

  // '-m' option //{{{
  int specific_moltype_for_size[Counts.TypesOfMolecules];
  // set all to be used when '-m' option is missing
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", specific_moltype_for_size,
                          Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption_old(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }
  // TODO those ridiculous flags are everywhere!
  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  }
  // count total number of chains in excluded aggs
  long int exclude_count_chains = 0;
  // count total number of excluded aggs
  long int exclude_count_agg = 0; //}}}

  // '-only' option //{{{
  int only_specific_moltype_aggregates[Counts.TypesOfMolecules];
  if (MoleculeTypeOption2(argc, argv, "-only", only_specific_moltype_aggregates,
                          Counts, &MoleculeType)) {
    exit(1);
  }
  // is '-only' used?
  bool only = false;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (only_specific_moltype_aggregates[i] == 1) {
      only = true;
      break;
    }
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption_old_old(argc, argv, "-bt", true, Counts, &BeadType)) {
    exit(0);
  }

  bool types_for_gyr[Counts.TypesOfBeads];
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    types_for_gyr[i] = BeadType[i].Use;
    BeadType[i].Use = false;
  } //}}}

  // write initial stuff to output file //{{{
  FILE *out = OpenFile(output, "w");
  PrintByline(out, argc, argv);
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

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // open input aggregate file and skip the first lines (Aggregate command & blank line) //{{{
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO go for while(fgets()); treatment
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}


  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules, sizeof (AGGREGATE));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded, sizeof *Aggregate[i].Monomer);
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded, sizeof *Aggregate[i].Bead);
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules, sizeof *Aggregate[i].Molecule);
  } //}}}

  if (verbose) { //{{{
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // TODO memory allocation... Aaargh!
  // allocate memory for sum of various things //{{{
  // numbers of aggregates of all possibe sizes
  int *agg_counts_sum = calloc(Counts.Molecules, sizeof *agg_counts_sum);
  // total radius of gyration: [size][0] normal sum, [size][1] sum of Rg*mass, [size][2] Rg*mass^2
  double (*Rg_sum)[3] = calloc(Counts.Molecules, sizeof *Rg_sum);
  // total square of radius of gyration: [size][0] normal sum, [size][1] sum of Rg^2*mass, [size][2] Rg^2*mass^2
  double (*sqrRg_sum)[3] = calloc(Counts.Molecules, sizeof *sqrRg_sum);
  // relative shape anisotropy: only normal sum
  double *Anis_sum = calloc(Counts.Molecules, sizeof *Anis_sum);
  // acylindricity: only normal sum
  double *Acyl_sum = calloc(Counts.Molecules, sizeof *Acyl_sum);
  // asphericity: only normal sum
  double *Aspher_sum = calloc(Counts.Molecules, sizeof *Aspher_sum);
  // gyration tensor eigenvalues
  VECTOR *eigen_sum = calloc(Counts.Molecules,sizeof (VECTOR));
  // total mass of aggregates: [size][0] normal sum, [size][1] sum of squares
  long int (*mass_sum)[2] = calloc(Counts.Molecules, sizeof *mass_sum);
  // number of molecule types in aggregates: [size][mol type] only normal sum
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  } //}}}

  // open input coordinate file
  FILE *vcf= OpenFile(input_coor, "r");
  SkipVtfStructure(vcf, struct_lines);

  // main loop //{{{
  count = 0; // count timesteps
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    // TODO: will change (probably)
    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead,
                   MoleculeType, &Molecule, Index);
    // TODO: fractionals - should the calculatio be in fractionals? I guess so.
    //       why though? The aggregates are joined and that's that...
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    if (!joined) {
      RemovePBCMolecules2(Counts, Box, BeadType, &Bead,
                         MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, Box.Length,
                          BeadType, &Bead, MoleculeType, Molecule);
    }

    // TODO: allocation... Aaargh!
    // allocate arrays for the timestep //{{{
    int *agg_counts_step = calloc(Counts.Molecules, sizeof *agg_counts_step);
    double (*Rg_step)[3] = calloc(Counts.Molecules, sizeof *Rg_step);
    double (*sqrRg_step)[3] = calloc(Counts.Molecules, sizeof *sqrRg_step);
    double *Anis_step = calloc(Counts.Molecules, sizeof *Anis_step);
    double *Acyl_step = calloc(Counts.Molecules, sizeof *Acyl_step);
    double *Aspher_step = calloc(Counts.Molecules,sizeof *Aspher_step);
    VECTOR *eigen_step = calloc(Counts.Molecules, sizeof (VECTOR)); //}}}

  // TODO check
    // calculate shape descriptors //{{{
    double mass_step[2] = {0}; // total mass of aggregates in a step: [0] normal, [1] sum of squares
    for (int i = 0; i < Counts.Aggregates; i++) {

      // test if aggregate 'i' should be used //{{{
      int size = 0;
      // agg size = number of molecules of type 'specific_moltype_for_size'
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
        }
      }
      int correct_size = size - 1; // all agg sizes are used - relic of old times
      // make calculations only if agg size is well defined and within given range
      if (size == 0 || size < range_As[0] || size > range_As[1]) {
        continue;
      } //}}}

      // if '-only' is used, use only aggregates composed the specified molecule(s) //{{{
      test = true;
      if (only) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (only_specific_moltype_aggregates[Molecule[id].Type] == 0) {
            test = false; // a molecule is not of the required type
            break;
          }
        }
        if (!test) { // should the rest of the for loop agg i be skipped?
          continue;
        }
      } //}}}

      // if '-x' option is used, discount aggregates with only specified molecule type(s) //{{{
      test = false;
      for (int j = 0; j < size; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Write) {
          test = true; // a molecule that shouldn't be in agg 'i' is there
          break;
        }
      }
      if (!test) { // should the rest of the for loop agg i be skipped?
        exclude_count_chains += Aggregate[i].nMolecules;
        exclude_count_agg++;
        continue;
      } //}}}

      if (correct_size != -1) {
        // copy bead ids to a separate array //{{{
        int *list = malloc(Aggregate[i].nBeads*sizeof(int));
        int n = 0;
        double agg_mass = 0;
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          int btype = Bead[id].Type;
          if (types_for_gyr[btype]) {
            list[n] = id;
            n++;
            agg_mass += BeadType[Bead[id].Type].Mass;
          }
        } //}}}

      // TODO fractionals
        VECTOR eigen = Gyration(n, list, Counts, BeadType, &Bead);

        free(list); // free array of bead ids for gyration calculation

        double Rgi = sqrt(eigen.x + eigen.y + eigen.z);

        if (eigen.x < 0 || eigen.y < 0 || eigen.z < 0) {
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "Warning: negative eigenvalues (");
          ColourChange(STDERR_FILENO, CYAN);
          fprintf(stderr, "%lf", eigen.x);
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, ", ");
          ColourChange(STDERR_FILENO, CYAN);
          fprintf(stderr, "%lf", eigen.y);
          ColourChange(STDERR_FILENO, CYAN);
          fprintf(stderr, ", ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%lf", eigen.z);
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, ")\n\n");
          ColourReset(STDERR_FILENO);
        }
        // agg masses
        mass_step[0] += agg_mass; // for this timestep
        mass_step[1] += SQR(agg_mass); // for this timestep
        // radius of gyration
        Rg_step[correct_size][0] += Rgi; // for number avg
        Rg_step[correct_size][1] += Rgi * agg_mass; // for weight average
        Rg_step[correct_size][2] += Rgi * SQR(agg_mass); // for z-average
        // squared radius of gyration
        sqrRg_step[correct_size][0] += SQR(Rgi); // for number avg
        sqrRg_step[correct_size][1] += SQR(Rgi) * agg_mass; // for weight average
        sqrRg_step[correct_size][2] += SQR(Rgi) * SQR(agg_mass); // for z-average
        // relative shape anisotropy
        Anis_step[correct_size] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
        // acylindricity
        Acyl_step[correct_size] += eigen.y - eigen.x;
        // asphericity
        Aspher_step[correct_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
        // gyration vector eigenvalues
        eigen_step[correct_size].x += eigen.x;
        eigen_step[correct_size].y += eigen.y;
        eigen_step[correct_size].z += eigen.z;
        // aggregate count
        agg_counts_step[correct_size]++;

        // sum molecules and aggregates (if step between start and end)
        if (count >= start && (end == -1 || count <= end)) {
          agg_counts_sum[correct_size]++;
          // aggregate mass
          mass_sum[correct_size][0] += agg_mass;
          mass_sum[correct_size][1] += SQR(agg_mass);
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
            molecules_sum[correct_size][mol_type]++;
          }
        }
      }
    } //}}}

    // add values to sums //{{{
    if (count >= start && (end == -1 || count <= end)) {
      for (int i = 0; i < Counts.Molecules; i++) {
        Rg_sum[i][0] += Rg_step[i][0];
        Rg_sum[i][1] += Rg_step[i][1];
        Rg_sum[i][2] += Rg_step[i][2];
        sqrRg_sum[i][0] += sqrRg_step[i][0];
        sqrRg_sum[i][1] += sqrRg_step[i][1];
        sqrRg_sum[i][2] += sqrRg_step[i][2];
        Anis_sum[i] += Anis_step[i];
        Acyl_sum[i] += Acyl_step[i];
        Aspher_sum[i] += Aspher_step[i];
        eigen_sum[i].x += eigen_step[i].x;
        eigen_sum[i].y += eigen_step[i].y;
        eigen_sum[i].z += eigen_step[i].z;
      }
    } //}}}

    // print data to output file
    out = OpenFile(output, "a");

    fprintf(out, "%5d", count); // timestep
    // sum up contributions from all aggregate sizes
    for (int i = 1; i < Counts.Molecules; i++) {
      Rg_step[0][0] += Rg_step[i][0];
      Rg_step[0][1] += Rg_step[i][1];
      Rg_step[0][2] += Rg_step[i][2];
      sqrRg_step[0][0] += sqrRg_step[i][0];
      sqrRg_step[0][1] += sqrRg_step[i][1];
      sqrRg_step[0][2] += sqrRg_step[i][2];
      Anis_step[0] += Anis_step[i];
      Acyl_step[0] += Acyl_step[i];
      Aspher_step[0] += Aspher_step[i];
      eigen_step[0].x += eigen_step[i].x;
      eigen_step[0].y += eigen_step[i].y;
      eigen_step[0].z += eigen_step[i].z;

      agg_counts_step[0] += agg_counts_step[i];
    }
    // <R_G>
    fprintf(out, " %lf %lf %lf", Rg_step[0][0]/agg_counts_step[0], Rg_step[0][1]/mass_step[0], Rg_step[0][2]/mass_step[1]);
    // <R_G^2>
    fprintf(out, " %lf %lf %lf", sqrRg_step[0][0]/agg_counts_step[0], sqrRg_step[0][1]/mass_step[0], sqrRg_step[0][2]/mass_step[1]);
    // relative shape anisotropy
    fprintf(out, " %lf", Anis_step[0]/agg_counts_step[0]);
    // acylindricity
    fprintf(out, " %lf", Acyl_step[0]/agg_counts_step[0]);
    // asphericity
    fprintf(out, " %lf", Aspher_step[0]/agg_counts_step[0]);
    // eigenvalues
    fprintf(out, " %lf %lf %lf", eigen_step[0].x/agg_counts_step[0], eigen_step[0].y/agg_counts_step[0], eigen_step[0].z/agg_counts_step[0]);
    putc('\n', out);

    fclose(out); //}}}

    // free memory //{{{
    free(agg_counts_step);
    free(Rg_step);
    free(sqrRg_step);
    free(Anis_step);
    free(Acyl_step);
    free(Aspher_step);
    free(eigen_step); //}}}

    // if there's no additional timestep, exit the while loop
    if (LastStep(vcf, NULL)) {
      break;
    }
  }
  fclose(vcf);
  fclose(agg);
  // print last step count?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count);
  } //}}}

// TODO check
  // calculate per-size averages? //{{{
  if (per_size_file[0] != '\0') {
    out = OpenFile(per_size_file, "w");
    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);

    fprintf(out, "# column: (1) agg size");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, " (%d) <%s>", i+2, MoleculeType[i].Name);
    }
    fprintf(out, " (%d) <Rg>, ", Counts.TypesOfMolecules+2);
    fprintf(out, " (%d) <Rg^2>, ", Counts.TypesOfMolecules+3);
    fprintf(out, " (%d) <Anis>, ", Counts.TypesOfMolecules+4);
    fprintf(out, " (%d) <Acyl>, ", Counts.TypesOfMolecules+5);
    fprintf(out, " (%d) <Aspher>, ", Counts.TypesOfMolecules+6);
    fprintf(out, " (%d) <eigen.x>, ", Counts.TypesOfMolecules+7);
    fprintf(out, " (%d) <eigen.y>, ", Counts.TypesOfMolecules+8);
    fprintf(out, " (%d) <eigen.z>, ", Counts.TypesOfMolecules+9);
    fprintf(out, " (%d) number of aggs", Counts.TypesOfMolecules+10);
    putc('\n', out);
    for (int i = 0; i < Counts.Molecules; i++) {
      if (agg_counts_sum[i] > 0) {
        fprintf(out, "%4d", i+1);
        for (int j = 0; j < Counts.TypesOfMolecules; j++) {
          fprintf(out, " %7.3f", (double)(molecules_sum[i][j])/agg_counts_sum[i]);
        }
        fprintf(out, " %7.3f", Rg_sum[i][0]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", sqrRg_sum[i][0]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", Anis_sum[i]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", Acyl_sum[i]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", Aspher_sum[i]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", eigen_sum[i].x/agg_counts_sum[i]);
        fprintf(out, " %7.3f", eigen_sum[i].y/agg_counts_sum[i]);
        fprintf(out, " %7.3f", eigen_sum[i].z/agg_counts_sum[i]);
        fprintf(out, " %d", agg_counts_sum[i]);
        putc('\n', out);
      }
    }

    fclose(out);
  } //}}}

  // total averages //{{{
  for (int i = 1; i < Counts.Molecules; i++) {
    Rg_sum[0][0] += Rg_sum[i][0];
    Rg_sum[0][1] += Rg_sum[i][1];
    Rg_sum[0][2] += Rg_sum[i][2];
    sqrRg_sum[0][0] += sqrRg_sum[i][0];
    sqrRg_sum[0][1] += sqrRg_sum[i][1];
    sqrRg_sum[0][2] += sqrRg_sum[i][2];
    Anis_sum[0] += Anis_sum[i];
    Acyl_sum[0] += Acyl_sum[i];
    Aspher_sum[0] += Aspher_sum[i];
    eigen_sum[0].x += eigen_sum[i].x;
    eigen_sum[0].y += eigen_sum[i].y;
    eigen_sum[0].z += eigen_sum[i].z;

    agg_counts_sum[0] += agg_counts_sum[i];

    mass_sum[0][0] += mass_sum[i][0];
    mass_sum[0][1] += mass_sum[i][1];

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules_sum[0][j] += molecules_sum[i][j];
    }
  }

  // print to output file
  out = OpenFile(output, "a");

  fprintf(out, "# (1) <M>_n, (2) <M>_w ");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, "(%d) <%s>, ", i+3, MoleculeType[i].Name);
  }
  fprintf(out, "(%d) <Rg>_n, ", Counts.TypesOfMolecules+3);
  fprintf(out, "(%d) <Rg>_w, ", Counts.TypesOfMolecules+4);
  fprintf(out, "(%d) <Rg>_z, ", Counts.TypesOfMolecules+5);
  fprintf(out, "(%d) <Rg^2>_n, ", Counts.TypesOfMolecules+6);
  fprintf(out, "(%d) <Rg^2>_w, ", Counts.TypesOfMolecules+7);
  fprintf(out, "(%d) <Rg^2>_z, ", Counts.TypesOfMolecules+8);
  fprintf(out, "(%d) <Anis>, ", Counts.TypesOfMolecules+9);
  fprintf(out, "(%d) <Acyl>, ", Counts.TypesOfMolecules+10);
  fprintf(out, "(%d) <Aspher>, ", Counts.TypesOfMolecules+11);
  fprintf(out, "(%d) <eigen.x>, ", Counts.TypesOfMolecules+12);
  fprintf(out, "(%d) <eigen.y>, ", Counts.TypesOfMolecules+13);
  fprintf(out, "(%d) <eigen.z>, ", Counts.TypesOfMolecules+14);
  putc('\n', out);
  fprintf(out, "# %lf", (double)(mass_sum[0][0])/agg_counts_sum[0]); //<M_As>_n
  fprintf(out, " %lf", (double)(mass_sum[0][1])/mass_sum[0][0]); //<M_As>_w
  // molecule types
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
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
  fprintf(out, " %lf", eigen_sum[0].x/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0].y/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0].z/agg_counts_sum[0]);
  putc('\n', out);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(Counts, &Aggregate);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  for (int i = 0; i < Counts.Molecules; i++) {
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
  free(stuff); //}}}

  return 0;
}
