#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
        TODO: -m_id option!!!\n\n\
DistrAgg calculates average aggregation numbers and aggregate masses during \
the simulation run (i.e., time evolution) as well as overall distributions \
and volume fractions of aggregates. The definition of aggregate size is quite \
flexible and only a specified range can be used. Also, this utility can \
fanalyse composition of specified aggreget size(s) and write compositition \
distribution (i.e., distribution of numbers of different molecules in over \
all aggregates with given size).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <in.agg> <distr file> <avg file> [options>\n\n", cmd);

  fprintf(ptr, "   <in.agg>       input agg file\n");
  fprintf(ptr, "   <distr file>   output file with distribution \
of aggregate sizes\n");
  fprintf(ptr, "   <avg file>     output file with per-timestep averages\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -st <int>            starting timestep \
for calculation\n");
  fprintf(ptr, "      -e <end>             ending timestep for calculation\n");
  fprintf(ptr, "      -n <int> <int>       use aggregate sizes in \
a given range\n");
  fprintf(ptr, "      -m <name(s)>         agg size means number of <name(s)> \
molecules in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>         exclude aggregates containing \
only specified molecule(s)\n");
  fprintf(ptr, "      --only <name(s)>     use only aggregates composed \
of specified molecule(s)\n");
  fprintf(ptr, "      -c <name> <int(s)>   write composition distribution \
of aggregate size(s) to <name>\n");
  fprintf(ptr, "      -nc <name> <name>    two molecule names \
for composition calculation\n");
  fprintf(ptr, "      -m_id <int>          calculate only for aggregates \
containing the molecule with vsf index <int>\n");
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
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "--only") != 0 &&
        strcmp(argv[i], "-nc") != 0 &&
        strcmp(argv[i], "-c") != 0 &&
        strcmp(argv[i], "-m_id") != 0 ) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  snprintf(input_agg, LINE, "%s", argv[++count]);
  // test if <in.agg> ends with '.agg'
  int ext = 1;
  char extension[2][EXTENSION];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

// TODO this will change with the change to agg file
// TODO why open the file here? Later should be more than enough
  // open <in.agg> and skip the first two lines //{{{
  FILE *agg = OpenFile(input_agg, "r");
  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // <distr file> - file with distribution of aggregation numbers
  char out_distr[LINE] = "";
  snprintf(out_distr, LINE, "%s", argv[++count]);

  // <avg file> - file with per-timestep average aggregation numbers
  char out_avg[LINE] = "";
  snprintf(out_distr, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  char input_vsf[LINE] = "";
  snprintf(input_vsf, LINE, "traject.vsf");
  CommonOptions_old(argc, argv, input_vsf, &verbose, &silent, LINE);
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
  FullVtfRead(input_vsf, "\0", false, false, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // '-m_id' option (TODO: consider what it should override) //{{{
  int m_id = -1; // no -m_id option
  if (IntegerOption(argc, argv, "-m_id", &m_id)) {
    exit(1);
  } //}}}

  // '-n' option - range of aggregation numbers //{{{
  int range_As[2], test = 2;
  range_As[0] = 1;
  range_As[1] = Counts.Molecules;
  if (MultiIntegerOption(argc, argv, "-n", &test, range_As)) {
    exit(1);
  }
  // TODO colours
  if (test != 2) {
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\033[1;33m-n\033[1;31m option requires two numeric arguments\n\n");
    fprintf(stderr, "\033[0m");
    exit(1);
  }

  // make sure first number is smaller
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
  }

  // TODO why? just use specific_moltype_for_size[] everywhere...
  // is '-m' used?
  bool m_option = false;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (specific_moltype_for_size[i] == 0) {
      m_option = true;
      break;
    }
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption_old(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  }

  // count total number of chains in excluded aggs
  long int exclude_count_chains = 0;
  // count total number of excluded aggs
  long int exclude_count_aggs = 0; //}}}

  // '-c' option - aggregate sizes (may be complemented by '-nc' option) //{{{
  int multiply = 1000;
  int composition[100] = {0}, // list of agg sizes
      comp_number_of_sizes = 0, // number of agg sizes
      types[2][2] = {{-1},{-1}}; // [x][0]: mol type; [x][1]: number of mols
  char output_comp[LINE];
  if (FileIntsOption(argc, argv, "-c", composition, &comp_number_of_sizes, output_comp)) {
    exit(1);
  }
  if (comp_number_of_sizes > 0) {
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (specific_moltype_for_size[i] && types[0][0] == -1) {
        types[0][0] = i;
      } else if (specific_moltype_for_size[i] && types[0][0] != -1 && types[1][0] == -1) {
        types[1][0] = i;
      } else if (specific_moltype_for_size[i]) {
        // TODO colours
        ErrorPrintError_old();
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\033[1;33m-c\033[1;31m option - more than two molecule types for composition distribution\n\n");
        fprintf(stderr, "\033[0m");
        exit(1);
      }
    }
    if (types[0][0] == -1 || types[1][0] == -1) {
      // TODO colours
      ErrorPrintError_old();
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\033[1;33m-c\033[1;31m option - less than two molecule types for composition distribution\n");
      fprintf(stderr, "\033[0m");
      exit(1);
    }
    for (int i = 0; i < comp_number_of_sizes; i++) {
      composition[i]+=0;
    }
  } //}}}

  // '-nc' option must be complemented by '-c' option //{{{
  int composition_mol_names[Counts.TypesOfMolecules];
  // TODO zeroize composition_mol_names[]?
  if (MoleculeTypeOption2(argc, argv, "-nc", composition_mol_names,
                          Counts, &MoleculeType)) {
    exit(1);
  }
  // check number of provided names
  count = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (composition_mol_names[i] == 1) {
      count++;
    }
  }
  // error if wrong number of names
  // TODO colours
  if (count != 0 && count != 2) {
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\033[1;33m-nc\033[1;31m option - exactly two molecule names are required\n\n");
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  // assign provided molecule names
  if (count == 2) {
    // TODO colours
    printf("%d %d\n", composition_mol_names[0], composition_mol_names[1]);
    printf("%s %s\n", MoleculeType[types[0][0]].Name, MoleculeType[types[1][0]].Name);
    types[0][0] = composition_mol_names[0];
    types[1][0] = composition_mol_names[1];
    printf("%s %s\n", MoleculeType[types[0][0]].Name, MoleculeType[types[1][0]].Name);
  }
  // warning if -nc is used, but not -c
  if (count == 2 && comp_number_of_sizes == 0 && !silent) {
    // TODO colours
    fprintf(stderr, "\033[1;33m");
    fprintf(stdout, "\nWarning: '-nc' option has no effect if '-c' option is not present\n\n");
    fprintf(stderr, "\033[0m");
  }
  //}}}

  // '--only' option //{{{
  int only_specific_moltype_aggregates[Counts.TypesOfMolecules];
  // TODO zeroize the array
  if (MoleculeTypeOption2(argc, argv, "--only",
                          only_specific_moltype_aggregates,
                          Counts, &MoleculeType)) {
    exit(1);
  }
  // is '--only' used?
  bool only = false;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (only_specific_moltype_aggregates[i] == 1) {
      only = true;
      break;
    }
  }

  // count total number chains in specified aggs
  long int only_count_chains = 0;
  // count total number of specified aggs
  long int only_count_aggs = 0; //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules, sizeof (AGGREGATE));
  for (int i = 0; i < Counts.Molecules; i++) {
    // at most, all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,
                                   sizeof *Aggregate[i].Molecule);
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded, sizeof *Aggregate[i].Bead);
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded,
                                  sizeof *Aggregate[i].Monomer);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    Box.Length.x = -1;
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

// TODO memory allocation... Aaargh!
  // arrays for distribution //{{{
  // number distribution
  long double ndistr[Counts.Molecules];
  // weight and z distributions - [][0] = mass of mols according to options; [][1] = mass of whole agg
  long double wdistr[Counts.Molecules][2]; //TODO: = {{0},{0}}; ?
  long double zdistr[Counts.Molecules][2];
  // volume distribution - probably works, but not really needed, so not sure
  // [][0] = volume of mols according to options; [][1] = volume of whole agg
  long double voldistr[Counts.Molecules][2];
  // number distribution of agg composition
  long int *ndistr_comp[Counts.Molecules]; // TODO: should have size of the number of provided sizes
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr_comp[i] = malloc((multiply*Counts.Molecules)*sizeof(long int));
  }
  // number of aggregates throughout simulation
  int count_agg[Counts.Molecules];
  // molecule typs in aggregates: [agg size][mol type][number or SQR(number)]
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  }

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr[i] = 0;
    wdistr[i][0] = 0;
    wdistr[i][1] = 0;
    zdistr[i][0] = 0;
    zdistr[i][1] = 0;
    voldistr[i][0] = 0;
    voldistr[i][1] = 0;
    count_agg[i] = 0;
  } //}}}

  // print the first two lines to output file with per-step averages //{{{
  FILE *out = OpenFile(out_avg, "w");
  PrintByline(out, argc, argv);
  count = 1;
  fprintf(out, "# column: (%d) step, ", count++);
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out, "(%d) <M>_n, ", count++);
    count++;
    fprintf(out, "(%d) <M>_w, ", count++);
    count++;
    fprintf(out, "(%d) <M>_z, ", count++);
    count++;
    fprintf(out, "(%d) <As>_n, ", count++);
    count++;
    fprintf(out, "(%d) <As>_w, ", count++);
    count++;
    fprintf(out, "(%d) <As>_z, ", count++);
    count++;
  } else {
    fprintf(out, "(%d) <M>_n, ", count++);
    fprintf(out, "(%d) <M>_w, ", count++);
    fprintf(out, "(%d) <M>_z, ", count++);
    fprintf(out, "(%d) <As>_n, ", count++);
    fprintf(out, "(%d) <As>_w, ", count++);
    fprintf(out, "(%d) <As>_z, ", count++);
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, "(%d) <%s>_n, ", count++, MoleculeType[i].Name);
  }
  fprintf(out, "(%d) number of aggregates", count++);
  putc('\n', out);
  fclose(out); //}}}

  // main loop //{{{
  int count_step = 0;
  // [0][] = simple sum, [1][] = sum of squares, [2][] = sume of cubes
  // [][0] = mass of mols in agg from options, [][1] = mass of the whole aggregate
  double mass_sum[3][2] = {{0}}, As_sum[3][2] = {{0}};
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    count_step++;
    // print step? //{{{
    if (count_step < start) {
      if (!silent && isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\rDiscarding Step: %d", count_step);
      }
    } else if (count_step == start) {
      if (!silent) {
        if (isatty(STDOUT_FILENO)) {
          fflush(stdout);
          fprintf(stdout, "\r                          \r");
        }
        fprintf(stdout, "Starting step: %d\n", start);
      }
    } else {
      if (!silent && isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\rStep: %d", count_step);
      }
    } //}}}

    // TODO: will change (probably)
    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead,
                   MoleculeType, &Molecule, Index);

    // TODO check
    // go through all aggregates
    int aggs_step = 0; // number of eligible aggregates per step
    double avg_mass_n_step[2] = {0}, avg_mass_w_step[2] = {0}, avg_mass_z_step[2] = {0}, // per-step mass averages
           avg_As_n_step[2] = {0}, avg_As_w_step[2] = {0},  avg_As_z_step[2] = {0}, // per-step As averages
           molecules_step[Counts.TypesOfMolecules];
    // zeroize per-step counts of molecule types
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      molecules_step[i] = 0;
    }
    for (int i = 0; i < Counts.Aggregates; i++) {

      /* determine aggregate size and mass: //{{{
       * if `-m` is used, the size is the number of 'specific_moltype_for_size'
       * mols and the mass is the mass only of those molecules */
      int size = 0;
      double agg_mass = 0, agg_vol = 0; // mass and volume according to options
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
          agg_mass += MoleculeType[mol_type].Mass;
          agg_vol += MoleculeType[mol_type].nBeads; // all DPD beads have same volume
        }
      }
      // make calculations only if agg size is well defined and within given range
      if (size == 0 || size < range_As[0] || size > range_As[1]) {
        continue;
      } //}}}

      // if '--only' is used, use only aggregates composed the specified molecule(s) //{{{
      bool use_agg = true;
      if (only) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (only_specific_moltype_aggregates[Molecule[id].Type] == 0) {
            use_agg = false; // a molecule is not of the required type
            break;
          }
        }
        if (!use_agg) { // should the rest of the for loop agg i be skipped?
          continue;
        }
      } //}}}

      // if '-x' option is used, discount aggregates with only specified molecule type(s) //{{{
      use_agg = false;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Write) {
          use_agg = true; // a molecule that shouldn't be in agg 'i' is there
          break;
        }
      }
      if (!use_agg) { // should the rest of the for loop agg i be skipped?
        exclude_count_chains += Aggregate[i].nMolecules;
        exclude_count_aggs++;
        continue;
      } //}}}

      // if '-m_id' is used, check if specified resid from vsf is in aggregate
      if (m_id != -1) {
        test = false;
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (m_id == (id+1)) { // resname in vsf start from 1
            test = true;
            break;
          }
        }
        if (!test) {
          continue;
        }
      }

      // if '-c' is used, calculate composition distribution //{{{
      if (comp_number_of_sizes > 0) {
        // use the size?
        int use_size_comp = -1;
        for (int j = 0; j < comp_number_of_sizes; j++) {
          if (size == composition[j]) {
            use_size_comp = j;
            break;
          }
        }

        // TODO: ndistr_comp array size (should be number of provided sizes)
        if (use_size_comp != -1) {
          // count numbers of the two mol types in aggregate 'i'
          types[0][1] = 0;
          types[1][1] = 0;
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            // molecule type
            int id = Molecule[Aggregate[i].Molecule[j]].Type;
            if (types[0][0] == id) {
              types[0][1]++;
            } else if (types[1][0] == id) {
              types[1][1]++;
            }
          }
          double ratio;
          if (types[1][1] == 0) {
            ratio = 1;
          } else {
            ratio = (double)(types[0][1]) / types[1][1];
          }
          ndistr_comp[size-1][(int)(ratio*multiply)]++;
          printf("%d / %d = %lf; ndistr_comp[%d][%d] = %ld\n", types[0][1], types[1][1], ratio,
                                                               size-1, (int)(ratio*multiply),
                                                               ndistr_comp[size-1][(int)(ratio*multiply)]);
        }
      } //}}}

      // for average aggregate mass during the step
      avg_mass_n_step[0] += agg_mass;
      avg_mass_w_step[0] += SQR(agg_mass);
      avg_mass_z_step[0] += CUBE(agg_mass);
      avg_mass_n_step[1] += Aggregate[i].Mass;
      avg_mass_w_step[1] += SQR(Aggregate[i].Mass);
      avg_mass_z_step[1] += CUBE(Aggregate[i].Mass);

      // for average aggregation number during the step
      avg_As_n_step[0] += size;
      avg_As_w_step[0] += size * agg_mass;
      avg_As_z_step[0] += size * SQR(agg_mass);
      avg_As_n_step[1] += Aggregate[i].nMolecules;
      avg_As_w_step[1] += Aggregate[i].nMolecules * Aggregate[i].Mass;
      avg_As_z_step[1] += Aggregate[i].nMolecules * SQR(Aggregate[i].Mass);
      // molecule species numbers
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mtype = Molecule[Aggregate[i].Molecule[j]].Type;
        molecules_step[mtype]++;
      }

      aggs_step++;

      // start calculation of averages and distributions from specified 'start' timestep
      if (count_step >= start && (end == -1 || count_step <= end)) {
        count_agg[size-1]++;
        // distribution
        ndistr[size-1]++;
        wdistr[size-1][0] += agg_mass;
        wdistr[size-1][1] += Aggregate[i].Mass;
        zdistr[size-1][0] += SQR(agg_mass);
        zdistr[size-1][1] += SQR(Aggregate[i].Mass);
        voldistr[size-1][0] += agg_vol;
        voldistr[size-1][1] += Aggregate[i].nBeads;
        // summed up sizes
        As_sum[0][0] += size;
        As_sum[1][0] += size * agg_mass;
        As_sum[2][0] += size * SQR(agg_mass);
        As_sum[0][1] += Aggregate[i].nMolecules;
        As_sum[1][1] += Aggregate[i].nMolecules * Aggregate[i].Mass;
        As_sum[2][1] += Aggregate[i].nMolecules * SQR(Aggregate[i].Mass);
        // summed up masses
        mass_sum[0][0] += agg_mass;
        mass_sum[1][0] += SQR(agg_mass);
        mass_sum[2][0] += CUBE(agg_mass);
        mass_sum[0][1] += Aggregate[i].Mass;
        mass_sum[1][1] += SQR(Aggregate[i].Mass);
        mass_sum[2][1] += CUBE(Aggregate[i].Mass);
        // molecule species numbers
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_sum[size-1][mol_type]++;
        }

        // count chains if `--only` is used //{{{
        if (only) {
          only_count_chains += Aggregate[i].nMolecules; // superfluous
          only_count_aggs++;
        } //}}}
      }
    }

    // print averages to output file //{{{
    out = OpenFile(out_avg, "a");
    fprintf(out, "%5d", count_step); // step
    if (aggs_step > 0) {
      fprintf(out, " %10.5f", avg_mass_n_step[0]/aggs_step); // <mass>_n
      if (m_option) {
        fprintf(out, " %10.5f", avg_mass_n_step[1]/aggs_step); // <mass>_n (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_mass_w_step[0]/avg_mass_n_step[0]); // <mass>_w
      if (m_option) {
        fprintf(out, " %10.5f", avg_mass_w_step[1]/avg_mass_n_step[1]); // <mass>_w (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_mass_z_step[0]/avg_mass_w_step[0]); // <mass>_z
      if (m_option) {
        fprintf(out, " %10.5f", avg_mass_z_step[1]/avg_mass_w_step[1]); // <mass>_z (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_As_n_step[0]/aggs_step); // <As>_n
      if (m_option) {
        fprintf(out, " %10.5f", avg_As_n_step[1]/aggs_step); // <As>_n (whole agg)
      }
      fprintf(out, " %10.5f", avg_As_w_step[0]/avg_mass_n_step[0]); // <As>_w
      if (m_option) {
        fprintf(out, " %10.5f", avg_As_w_step[1]/avg_mass_n_step[1]); // <As>_w (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_As_z_step[0]/avg_mass_w_step[0]); // <As>_z
      if (m_option) {
        fprintf(out, " %10.5f", avg_As_z_step[1]/avg_mass_w_step[1]); // <As>_z (whole agg mass)
      }
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        fprintf(out, " %10.5f", molecules_step[i]/aggs_step);
      }
    } else { // zero everywhere if there are no aggregates of the specified type
      fprintf(out, " %10.5f", 0.0); // <mass>_n
      if (m_option) {
        fprintf(out, " %10.5f", 0.0); // <mass>_n (whole agg mass)
      }
      fprintf(out, " %10.5f", 0.0); // <mass>_w
      if (m_option) {
        fprintf(out, " %10.5f", 0.0); // <mass>_w (whole agg mass)
      }
      fprintf(out, " %10.5f", 0.0); // <mass>_z
      if (m_option) {
        fprintf(out, " %10.5f", 0.0); // <mass>_z (whole agg mass)
      }
      fprintf(out, " %10.5f", 0.0); // <As>_n
      if (m_option) {
        fprintf(out, " %10.5f", 0.0); // <As>_n (whole agg)
      }
      fprintf(out, " %10.5f", 0.0); // <As>_w
      if (m_option) {
        fprintf(out, " %10.5f", 0.0); // <As>_w (whole agg mass)
      }
      fprintf(out, " %10.5f", 0.0); // <As>_z
      if (m_option) {
        fprintf(out, " %10.5f", 0.0); // <As>_z (whole agg mass)
      }
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        fprintf(out, " %10.5f", 0.0);
      }
    }
    fprintf(out, " %5d", aggs_step); // number of aggregates in the step
    // numbers of species
    putc('\n', out);
    fclose(out); //}}}
  }
  fclose(agg);
  // print last step //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_step);
  } //}}}
  //}}}

  // print the first two lines to output file with distributions //{{{
  out = OpenFile(out_distr, "w");
  PrintByline(out, argc, argv);
  count = 1;
  fprintf(out, "# column: ");
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out, "(%d) As, ", count++);
    fprintf(out, "(%d) F_n(As), ", count++);
    fprintf(out, "(%d) F_w(As), ", count++);
    count++;
    fprintf(out, "(%d) F_z(As), ", count++);
    count++;
    fprintf(out, "(%d) <volume distribution>, ", count++);
    count++;
  } else {
    fprintf(out, "(%d) As, ", count++);
    fprintf(out, "(%d) F_n(As), ", count++);
    fprintf(out, "(%d) F_w(As), ", count++);
    fprintf(out, "(%d) F_z(As), ", count++);
    fprintf(out, "(%d) <volume distribution>, ", count++);
  }
  fprintf(out, "(%d) number of aggs,", count++);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " (%d) <%s>_n", i+count, MoleculeType[i].Name);
    if (i != (Counts.TypesOfMolecules-1)) {
      putc(',', out);
    }
  }
  putc('\n', out);
  fclose(out); //}}}

  // print distributions to output file //{{{
  out = OpenFile(out_distr, "a");

  if (end == -1) {
    count_step = count_step - start + 1;
  } else {
    count_step = count_step - (start - 1) - (end - 1);
  }

  // normalization factors
  long int ndistr_norm = 0, wdistr_norm[2] = {0}, zdistr_norm[2] = {0}, voldistr_norm[2] = {0};
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr_norm += ndistr[i];

    wdistr_norm[0] += wdistr[i][0];
    wdistr_norm[1] += wdistr[i][1];

    zdistr_norm[0] += zdistr[i][0];
    zdistr_norm[1] += zdistr[i][1];

    voldistr_norm[0] += voldistr[i][0];
    voldistr_norm[1] += voldistr[i][1];
  }

  for (int i = 0; i < Counts.Molecules; i++) {
    if (count_agg[i] > 0) {
      fprintf(out, "%4d", i+1); // As
      fprintf(out, " %10.5f", (double)(ndistr[i])/ndistr_norm); // number As distr
      fprintf(out, " %10.5f", (double)(wdistr[i][0])/wdistr_norm[0]); // weight distr
      if (m_option) {
        fprintf(out, " %10.5f", (double)(wdistr[i][1])/wdistr_norm[1]); // weight distr (whole agg mass)
      }
      fprintf(out, " %10.5f", (double)(zdistr[i][0])/zdistr_norm[0]); // z distr
      if (m_option) {
        fprintf(out, " %10.5f", (double)(zdistr[i][1])/zdistr_norm[1]); // z distr (whole agg mass)
      }
      fprintf(out, " %10.5f", (double)(voldistr[i][0])/voldistr_norm[0]); // volume distribution
      if (m_option) {
        fprintf(out, " %10.5f", (double)(voldistr[i][1])/voldistr_norm[1]); // volume distribution (whole agg mass)
      }
      fprintf(out, " %6d", count_agg[i]); // number of aggregates
      // print average number of molecule types in aggregates
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        if (count_agg[i] == 0) {
          fprintf(out, "%8s", "?");
        } else {
          fprintf(out, " %10.5f", (double)(molecules_sum[i][j])/count_agg[i]);
        }
      }
      putc('\n', out);
    }
  }
  fclose(out); //}}}

  // print overall averages to avg and distr output files //{{{
  for (int i = 1; i < Counts.Molecules; i++) {
    count_agg[0] += count_agg[i]; // total number of aggregates
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules_sum[0][j] += molecules_sum[i][j]; // total number molecules of each type
    }
  }

  // open output files
  // distr file
  out = OpenFile(out_distr, "a");
  // avg file
  FILE *out2 = OpenFile(out_avg, "a");

  // print legend (with column numbers) //{{{
  putc('#', out);
  putc('#', out2);
  // distr file
  count = 1;
  fprintf(out, " (%d) <As>_n,", count++);
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out, " (%d) <As>_w,", count++);
    count++;
    fprintf(out, " (%d) <As>_z,", count++);
    count++;
    fprintf(out, " (%d) <M>_n,", count++);
    count++;
    fprintf(out, " (%d) <M>_w,", count++);
    count++;
    fprintf(out, " (%d) <M>_z,", count++);
    count++;
  } else {
    fprintf(out, " (%d) <As>_w,", count++);
    fprintf(out, " (%d) <As>_z,", count++);
    fprintf(out, " (%d) <M>_n,", count++);
    fprintf(out, " (%d) <M>_w,", count++);
    fprintf(out, " (%d) <M>_z,", count++);
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " (%d) <%s>_n,", ++count, MoleculeType[i].Name);
  }
  fprintf(out, " (%d) <n_agg>", count++);
  // avg file
  count = 1;
  fprintf(out2, " (%d) <As>_n,", count++);
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out2, " (%d) <As>_w,", count++);
    count++;
    fprintf(out2, " (%d) <As>_z,", count++);
    count++;
    fprintf(out2, " (%d) <M>_n,", count++);
    count++;
    fprintf(out2, " (%d) <M>_w,", count++);
    count++;
    fprintf(out2, " (%d) <M>_z,", count++);
    count++;
  } else {
    fprintf(out2, " (%d) <As>_w,", count++);
    fprintf(out2, " (%d) <As>_z,", count++);
    fprintf(out2, " (%d) <M>_n,", count++);
    fprintf(out2, " (%d) <M>_w,", count++);
    fprintf(out2, " (%d) <M>_z,", count++);
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out2, " (%d) <%s>_n,", count++, MoleculeType[i].Name);
  }
  fprintf(out2, " (%d) <n_agg>", count++);
  putc('\n', out);
  putc('\n', out2); //}}}

  // print the averages //{{{
  // distr file
  fprintf(out, "#");
  fprintf(out2, "#");

  // distr file
  if (count_agg[0] > 0) {
    fprintf(out, " %lf", As_sum[0][0]/count_agg[0]); // <As>_n
    fprintf(out, " %lf", As_sum[1][0]/mass_sum[0][0]); // <As>_w
    if (m_option) {
      fprintf(out, " %lf", As_sum[1][1]/mass_sum[0][1]); // <As>_w (whole agg)
    }
    fprintf(out, " %lf", As_sum[2][0]/mass_sum[1][0]); // <As>_z
    if (m_option) {
      fprintf(out, " %lf", As_sum[2][1]/mass_sum[1][1]); // <As>_z (whole agg)
    }

    fprintf(out, " %lf", mass_sum[0][0]/count_agg[0]); // <M>_n
    if (m_option) {
      fprintf(out, " %lf", mass_sum[0][1]/count_agg[0]); // <M>_n (whole agg)
    }
    fprintf(out, " %lf", mass_sum[1][0]/mass_sum[0][0]); // <M>_w
    if (m_option) {
      fprintf(out, " %lf", mass_sum[1][1]/mass_sum[0][1]); // <M>_w (whole agg)
    }
    fprintf(out, " %lf", mass_sum[2][0]/mass_sum[1][0]); // <M>_z
    if (m_option) {
      fprintf(out, " %lf", mass_sum[2][1]/mass_sum[1][1]); // <M>_z (whole agg)
    }
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, " %lf", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
    }
    fprintf(out, " %lf", (double)(count_agg[0])/count_step); // number of aggregates
  } else { // zero everywhere if no aggregates found
    fprintf(out, " 0.0  0.0  0.0  0.0  0.0  0.0  0.0");
    if (m_option) {
      fprintf(out, "  0.0  0.0  0.0  0.0  0.0  0.0");
    }
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, "  0.0");
    }
  }
  // avg file
  if (count_agg[0] > 0) {
    fprintf(out2, " %lf", As_sum[0][0]/count_agg[0]); // <As>_n
    fprintf(out2, " %lf", As_sum[1][0]/mass_sum[0][0]); // <As>_w
    if (m_option) {
      fprintf(out2, " %lf", As_sum[1][1]/mass_sum[0][1]); // <As>_w (whole agg)
    }
    fprintf(out2, " %lf", As_sum[2][0]/mass_sum[1][0]); // <As>_z
    if (m_option) {
      fprintf(out2, " %lf", As_sum[2][1]/mass_sum[1][1]); // <As>_z (whole agg)
    }

    fprintf(out2, " %lf", mass_sum[0][0]/count_agg[0]); // <M>_n
    if (m_option) {
      fprintf(out2, " %lf", mass_sum[0][1]/count_agg[0]); // <M>_n (whole agg)
    }
    fprintf(out2, " %lf", mass_sum[1][0]/mass_sum[0][0]); // <M>_w
    if (m_option) {
      fprintf(out2, " %lf", mass_sum[1][1]/mass_sum[0][1]); // <M>_w (whole agg)
    }
    fprintf(out2, " %lf", mass_sum[2][0]/mass_sum[1][0]); // <M>_z
    if (m_option) {
      fprintf(out2, " %lf", mass_sum[2][1]/mass_sum[1][1]); // <M>_z (whole agg)
    }
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out2, " %lf", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
    }
    fprintf(out2, " %lf", (double)(count_agg[0])/count_step); // number of aggregates
  } else { // zero everywhere if no aggregates found
    fprintf(out2, " 0.0  0.0  0.0  0.0  0.0  0.0  0.0");
    if (m_option) {
      fprintf(out2, "  0.0  0.0  0.0  0.0  0.0  0.0");
    }
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out2, "  0.0");
    }
  }

  putc('\n', out);
  putc('\n', out2); //}}}

  fclose(out);
  fclose(out2);
  //}}}

  // print composition distribution //{{{
  if (comp_number_of_sizes > 0) {
    out = OpenFile(output_comp, "w");

    // print command
    putc('#', out);
    PrintCommand(out, argc, argv);

    // print header line //{{{
    fprintf(out, "# column: (1) ");
    count = 0;
    // molecule names
    fprintf(out, "%s/%s; agg sizes - ", MoleculeType[types[0][0]].Name, MoleculeType[types[1][0]].Name);
    // aggregate sizes
    for (int i = 0; i < comp_number_of_sizes; i++) {
      fprintf(out, " (%d) %d;", i+2, composition[i]);
    }
    putc('\n', out); //}}}

    // normalization factor (simple sum)
    long int *sum = calloc(Counts.Molecules,sizeof(long int));
    int low = multiply * Counts.Molecules;
    bool *include = calloc((multiply*Counts.Molecules),sizeof(bool));
    for (int i = 0; i < (multiply*Counts.Molecules); i++) {
      for (int j = 0; j < comp_number_of_sizes; j++) {
        if (ndistr_comp[composition[j]-1][i] > 0 && low == (multiply*Counts.Molecules)) {
          low = i;
        }
        if (ndistr_comp[composition[j]-1][i] > 0) {
          include[i] = true;
        }
        sum[composition[j]-1] += ndistr_comp[composition[j]-1][i];
      }
    }

    // print data
    for (int i = low; i < (multiply*Counts.Molecules); i++) {
      if (include[i]) {
        fprintf(out, "%10.5f", (double)(i)/multiply);
        for (int j = 0; j < comp_number_of_sizes; j++) {
          if (ndistr_comp[composition[j]-1][i] != 0) {
            printf("%ld %ld\n", ndistr_comp[composition[j]-1][i], sum[composition[j]-1]);
            fprintf(out, " %7.5f", (double)(ndistr_comp[composition[j]-1][i])/sum[composition[j]-1]);
          } else {
            fprintf(out, "       ?");
          }
        }
        putc('\n', out);
      }
    }
    fclose(out);

    free(sum);
    free(include);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(Counts, &Aggregate);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(molecules_sum[i]);
    free(ndistr_comp[i]);
  }
  free(molecules_sum); //}}}

  return 0;
}
