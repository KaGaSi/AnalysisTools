#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
        REWRITE!!!\n\n\
DensityAggregates utility calculates radial bead density for aggregates \
of given size(s) from their centre of mass. For beads in molecules, \
it takes into account only beads from the current aggregate, not from \
any other aggregate. The utility does not check what molecule type a \
given bead belongs to; therefore, if the same bead type appears in \
more molecule types, the resulting density for that bead type will be \
averaged without regard for the various types of molecule it comes from \
(i.e., if 'mol1' and 'mol2' both both contain bead 'A', there will be \
only one column for 'A' bead type).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <in.agg> <width> <output> <size(s)> \
[options]\n\n", cmd);

  fprintf(ptr, "   <input>     input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <in.agg>    input agg file\n");
  fprintf(ptr, "   <width>     width of a single bin\n");
  fprintf(ptr, "   <output>    output orientation file \
(automatic '<size>.txt' ending)\n");
  fprintf(ptr, "   <size(s)>   aggregate sizes to calculate orientation for\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -n <ints>      bead indices (multiple of 2 <ints>; default: 1 2)\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> \
molecules in an aggregate\n");
  CommonHelp(error);
} //}}}

// function to test if provided indices are too high for a given molecule type //{{{
int test_ids(int id1, int id2, int mbeads) {
  if (id1 >= mbeads &&
      id2 >= mbeads) { // both higher than number beads - skip
    return 0;
  } else if (id1 >= mbeads &&
             id2 < mbeads) { // 'id1' too high, 'id2' okay
    return 1;
  } else if (id1 < mbeads &&
             id2 >= mbeads) { // 'id1' okay, 'id2' too high
    return 2;
  } else { // both are okay
    return 3;
  }
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
  int req_args = 5; //}}}

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
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-m") != 0) {

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
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
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

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with orientational parameter distribution
  char output[LINE];
  snprintf(output, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");
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

  // '-n' option - specify bead ids //{{{
  int bead[100] = {0}, // specified bead indices
      number_of_beads = 0, // number of parameters to -n
      beads_per_set = 2; // the numbers must come in pairs
  bead[0] = 1; // default ids
  bead[1] = 2;
  if (MultiIntegerOption(argc, argv, "-n", &number_of_beads, bead)) {
    exit(1);
  }
  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero
  }
  // errors & warnings //{{{
  // Error: wrong number of integers
  if (number_of_beads == 0 || (number_of_beads%beads_per_set) != 0) {
    ErrorPrintError_old();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "-n");
    RedText(STDERR_FILENO);
    fprintf(stderr, " - bead ids must be in pairs\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  // check all provided ids
  for (int i = 0; i < number_of_beads; i += beads_per_set) {
    // Error - two same beads
    if (bead[i] == bead[i+1]) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "-n");
      RedText(STDERR_FILENO);
      fprintf(stderr, " - two different bead ids are required (invalid pair: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d %d", bead[i], bead[i+1]);
      RedText(STDERR_FILENO);
      fprintf(stderr, ")\n");
      ResetColour(STDERR_FILENO);
      Help(argv[0], true);
      exit(1);
    }
    // Warning - id too high
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && (bead[i] >= MoleculeType[j].nBeads ||
                                  bead[i+1] >= MoleculeType[j].nBeads)) {
        YellowText(STDERR_FILENO);
        fprintf(stderr, "\nWarning: ");
        CyanText(STDERR_FILENO);
        fprintf(stderr, "-n");
        YellowText(STDERR_FILENO);
        fprintf(stderr, " - index in ");
        CyanText(STDERR_FILENO);
        fprintf(stderr, "%d %d", bead[i], bead[i+1]);
        YellowText(STDERR_FILENO);
        fprintf(stderr, " pair is larger than the number of beads in ");
        CyanText(STDERR_FILENO);
        fprintf(stderr, "%s", MoleculeType[j].Name);
        YellowText(STDERR_FILENO);
        fprintf(stderr, " molecule;");
    // TODO: does it really use the highest molecule's id
        fprintf(stderr, " using %d instead\n", MoleculeType[j].nBeads);
        ResetColour(STDERR_FILENO);
      }
    }
  } //}}}
  int number_of_pairs = number_of_beads / beads_per_set;
  //}}}

// TODO: allocation... aaargh
  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(sizeof(int *) * Counts.Molecules);
  int **agg_mols = malloc(sizeof(int *) * Counts.Molecules);
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2, sizeof *agg_sizes[i]);
    agg_mols[i] = calloc(Counts.TypesOfMolecules, sizeof *agg_mols[i]);
  }
  int aggs = 0; // number of aggregate sizes
  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (!IsInteger(argv[count]) || atoi(argv[count]) == 0) {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}
    agg_sizes[aggs][0] = atoi(argv[count]);
    // ensure output string isn't too long for attaching <size>.vcf
    if (agg_sizes[aggs][0] < 10) {
      output[LINE-1-4] = '\0';
    } else if (agg_sizes[aggs][0] < 100) {
      output[LINE-2-4] = '\0';
    } else if (agg_sizes[aggs][0] < 1000) {
      output[LINE-3-4] = '\0';
    } else if (agg_sizes[aggs][0] < 10000) {
      output[LINE-4-4] = '\0';
    } else {
      output[LINE-100] = '\0';
    }
    aggs++; // number of aggregate sizes
  } //}}}

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // open input aggregate file and skip the first lines
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  char line[LINE];
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  // number of bins
  int bins = 1.5 / width;

  // array for distributions //{{{
  double *orientation[Counts.TypesOfMolecules][number_of_pairs][aggs],
         avg_orientation[Counts.TypesOfMolecules][number_of_pairs][aggs];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      for (int k = 0; k < aggs; k++) {
        orientation[i][j][k] = calloc(bins, sizeof *orientation[i][j][k]);
        avg_orientation[i][j][k] = 0;
      }
    }
  } //}}}

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
    VerboseOutput(input_coor, Counts, Box, BeadType, Bead,
                  MoleculeType, Molecule);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

  count = SkipCoorAggSteps(vcf, input_coor, agg,
                           input_agg, Counts, start, silent);

  // use z axis as a normal TODO: option for axis //{{{
  VECTOR normal;
  normal.x = 0;
  normal.y = 0;
  normal.z = 1; //}}}

  // main loop //{{{
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // TODO: will change (probably)
    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead,
                   MoleculeType, &Molecule, Index);
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    if (!joined) {
      // transform coordinates into fractional ones for non-orthogonal box
      ToFractionalCoor(Counts.Beads, &Bead, Box);
      RemovePBCMolecules(Counts, Box, BeadType, &Bead,
                         MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, Box.Length,
                          BeadType, &Bead, MoleculeType, Molecule);
      // transform back to 'normal' coordinates for non-orthogonal box
      FromFractionalCoor(Counts.Beads, &Bead, Box);
    }

    // calculate orientation //{{{
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
      // is 'size' in provided list?
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j][0] == size) {
          correct_size = j;
        }
      } //}}}
      if (correct_size != -1) {
        agg_sizes[correct_size][1]++;

        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol = Aggregate[i].Molecule[j];
          int mtype = Molecule[mol].Type;

          for (int k = 0; k < number_of_beads; k += 2) {
            int mbeads = MoleculeType[mtype].nBeads;
            int id1 = -2, id2 = -1;
            // check provided ids for given molecule type //{{{
            switch(test_ids(bead[k], bead[k+1], mbeads)) {
              case 0: // bead[k] & bead[k+1] too high
                continue;
              case 1: // bead[k] too high & bead[k+1] okay
                id1 = Molecule[mol].Bead[mbeads-1];
                id2 = Molecule[mol].Bead[bead[k+1]];
                break;
              case 2: // bead[k] okay & bead[k+1] too high
                id1 = Molecule[mol].Bead[bead[k]];
                id2 = Molecule[mol].Bead[mbeads-1];
                break;
              case 3: // bead[k] & bead[k+1] okay
                id1 = Molecule[mol].Bead[bead[k]];
                id2 = Molecule[mol].Bead[bead[k+1]];
                break;
            } //}}}
            // if 'k' or 'k+1' is too high, it can lead to both being the highest id in molecule
            if (id1 == id2) {
              continue;
            }
            double length = SQR(Bead[id1].Position.x - Bead[id2].Position.x) +
                            SQR(Bead[id1].Position.y - Bead[id2].Position.y) +
                            SQR(Bead[id1].Position.z - Bead[id2].Position.z);
            length = sqrt(length);
            VECTOR unit;
            unit.x = (Bead[id1].Position.x - Bead[id2].Position.x) / length;
            unit.y = (Bead[id1].Position.y - Bead[id2].Position.y) / length;
            unit.z = (Bead[id1].Position.z - Bead[id2].Position.z) / length;

            double value = unit.x * normal.x + unit.y * normal.y + unit.z * normal.z; // cos(angle)
            value = 0.5 * (3 * SQR(value) - 1);

            // add to overall averages
            avg_orientation[mtype][k/2][correct_size] += value;

            // place into appropriate distribution bin
            int bin = (value + 0.5) / width; // +0.5, because its range is <-0.5,1>
            if (k < bins) {
              orientation[mtype][k/2][correct_size][bin]++;
            } else if (k == bins) { // in case of value=1 exactly
              orientation[mtype][k/2][correct_size][bins-1]++;
            }
          }
        }

        // count mol types in the aggregate
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol = Aggregate[i].Molecule[j];
          int mtype = Molecule[mol].Type;
          agg_mols[correct_size][mtype]++;
        }
      }
    } //}}}

    // exit the while loop if there's no more coordinates or -e step was reached
    if (LastStep(vcf, NULL) || end == count_vcf) {
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
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // write distribution of orientation parameters //{{{
  for (int temp = 0; temp < aggs; temp++) {
    FILE *out;
    char str[LINE];
    // assemble correct name & open the file
    snprintf(str, LINE, "%s%d.txt", output, agg_sizes[temp][0]);
    if ((out = fopen(str, "a")) == NULL) {
      ErrorFileOpen(str, 'a');
      exit(1);
    }
    // write initial stuff to the file //{{{
    PrintByline(out, argc, argv);
    // print first lines of output file - molecule names and beadtype pairs
    fprintf(out, "# columns: (1) Orientation order parameter");
    int count_print = 1;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, "\n#          %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += 2) {
        // are ids right?
        switch(test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads)) {
          case 0: // bead[k] & bead[k+1] too high
            continue;
          case 1: // bead[k] too high & bead[k+1] okay
            if ((bead[j+1]+1) != MoleculeType[i].nBeads) {
              fprintf(out, " (%d) %d-%d,", ++count_print, MoleculeType[i].nBeads, bead[j+1]+1);
            }
            break;
          case 2: // bead[k] okay & bead[k+1] too high
            if ((bead[j]+1) != MoleculeType[i].nBeads) {
              fprintf(out, " (%d) %d-%d,", ++count_print, bead[j]+1, MoleculeType[i].nBeads);
            }
            break;
          case 3: // bead[k] & bead[k+1] okay
            fprintf(out, " (%d) %d-%d,", ++count_print, bead[j]+1, bead[j+1]+1);
            break;
        }
      }
    }
    putc('\n', out);

    fclose(out); //}}}

    // calculate normalisation factor
    long int normalisation[Counts.TypesOfMolecules][number_of_beads/2];
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      for (int j = 0; j < number_of_beads; j+=2) {
        normalisation[i][j/2] = 0;
        for (int k = 0; k < bins; k++) {
          normalisation[i][j/2] += orientation[i][j/2][temp][k];
        }
      }
    }

    // write distribution to output file //{{{
    for (int i = 0; i < bins; i++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2-0.5); // orientation order parameter: <-0.5,1>
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        for (int k = 0; k < number_of_beads; k += 2) {
          switch(test_ids(bead[k], bead[k+1], MoleculeType[j].nBeads)) {
            case 0: // bead[k] & bead[k+1] too high
              continue;
            case 1: // bead[k] too high & bead[k+1] okay
              if ((bead[k+1]+1) == MoleculeType[j].nBeads) {
                continue;
              }
            case 2: // bead[k] okay & bead[k+1] too high
              if ((bead[k]+1) == MoleculeType[j].nBeads) {
                continue;
              }
            case 3: // bead[k] & bead[k+1] okay
              fprintf(out, "%10f", orientation[j][k/2][temp][i]/normalisation[j][k/2]);
              break;
          }
        }
      }
      putc('\n', out);
    } //}}}
    // write overall averages //{{{
    // legend line //{{{
    fprintf(out, "# Average orientation order parameter");
    count = 0;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, "\n#    %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += 2) {
        if (bead[j] < MoleculeType[i].nBeads && bead[j+1] < MoleculeType[i].nBeads) {
            fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
        } else if (bead[j] >= MoleculeType[i].nBeads && bead[j+1] < MoleculeType[i].nBeads) {
            fprintf(out, " (%d) %d-%d,", ++count, MoleculeType[i].nBeads, bead[j+1]+1);
        } else if (bead[j+1] >= MoleculeType[i].nBeads && bead[j] < MoleculeType[i].nBeads) {
            fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, MoleculeType[i].nBeads);
        }
      }
    }
    putc('\n', out); //}}}
    // data line
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      for (int j = 0; j < number_of_beads; j += 2) {
        int mbeads = MoleculeType[i].nBeads;
        // check provided ids for given molecule type
        if (bead[j] >= mbeads &&
            bead[j+1] >= mbeads) { // both higher than number beads - skip
          continue;
        }
        fprintf(out, " %10f", avg_orientation[i][j/2][temp]/agg_mols[temp][i]);
      }
    }
    putc('\n', out); //}}}
    fclose(out);
  } //}}}
printf("agg_sizes[0] = (%d, %d)\n", agg_sizes[0][0], agg_sizes[0][1]);

  // free memory - to make valgrind happy //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      for (int k = 0; k < aggs; k++) {
        free(orientation[i][j][k]);
      }
    }
  }
  FreeAggregate(Counts, &Aggregate);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
    free(agg_mols[i]);
  }
  free(agg_sizes);
  free(agg_mols); //}}}

  return 0;
}
