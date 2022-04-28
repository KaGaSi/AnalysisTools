#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Orientation utility calculates orientational order parameter for specified \
molecules and specified bead pairs in the molecule. Besides average value, \
it calculates distribution of the parameter.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>    width of a single distribution bin\n");
  fprintf(ptr, "   <output>   output file with distribution of the orientation \
parameter\n");
  fprintf(ptr, "   <mol(s)>   molecule name(s) to calculate orientational \
parameters for (optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --all       all molecules (overwrites <mol(s)>)\n");
  fprintf(ptr, "      --joined    specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -st <int>   starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>    ending timestep for calculation\n");
  fprintf(ptr, "      -a x/y/z    axis to use as the normal vector \
(default: z)\n");
  fprintf(ptr, "      -n <ints>   bead indices (multiple of 2 <ints>; \
default: 1 2)\n");
  CommonHelp(error);
} //}}}

// function to test if indices are too high for a given molecule type //{{{
int test_ids(int id1, int id2, int mol_beads) {
  if (id1 >= mol_beads &&
      id2 >= mol_beads) { // both higher than number beads - skip
    return 0;
  } else if (id1 >= mol_beads &&
             id2 < mol_beads) { // 'id1' too high, 'id2' okay
    return 1;
  } else if (id1 < mol_beads &&
             id2 >= mol_beads) { // 'id1' okay, 'id2' too high
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
  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }

  // use all molecules? ...do now to check correct number of arguments
  bool all = BoolOption(argc, argv, "--all");

  if (count < (req_args-1) || (count == (req_args-1) && !all)) {
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
        strcmp(argv[i], "--all") != 0 &&
        strcmp(argv[i], "--joined") &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-a") != 0 &&
        strcmp(argv[i], "-n") != 0) {

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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal_old(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with orientational parameter distribution //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");

  // axis to use as a normal //{{{
  int test = 'z'; // default value
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-a") == 0) {
      // Error - missing argument //{{{
      if ((i+1) >= argc) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "-a");
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing an x/y/z argument\n\n");
        ColourReset(STDERR_FILENO);
        Help(argv[0], true);
        exit(1);
      } //}}}
      test = argv[i+1][0];
    }
  }
  VECTOR normal;
  switch(test) {
    case 'x': // x
      // TODO: normal = {1,0,0} possible?
      normal.x = 1;
      normal.y = 0;
      normal.z = 0;
      break;
    case 'y': // y
      normal.x = 0;
      normal.y = 1;
      normal.z = 0;
      break;
    case 'z': // z
      normal.x = 0;
      normal.y = 0;
      normal.z = 1;
      break;
    default:
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "-a");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - argument must be 'x', 'y', or 'z'\n\n");
      ColourReset(STDERR_FILENO);
      Help(argv[0], true);
      exit(1);
  } //}}}
  //}}}

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

  // <molecule names> - types of molecules for calculation //{{{
  if (!all) { // --all option not used
    while (++count < argc && argv[count][0] != '-') {
      int mtype = FindMoleculeType(argv[count], Counts, MoleculeType);
      // error - nonexistent molecule  //{{{
      if (mtype == -1) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", input_coor);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - non-existent molecule type ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", argv[count]);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, "\n");
        ColourReset(STDERR_FILENO);
        ErrorMoleculeType(Counts, MoleculeType);
        exit(1);
      } //}}}
      MoleculeType[mtype].Use = true;
    }
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      MoleculeType[i].Use = true;
    }
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
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-n");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - bead ids must be in pairs\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  // check all provided ids
  for (int i = 0; i < number_of_beads; i += beads_per_set) {
    // Error - two same beads
    if (bead[i] == bead[i+1]) {
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, "\nWarning: ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "-n");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - two different bead ids are required (invalid pair: ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%d %d", bead[i], bead[i+1]);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, ")\n");
      ColourReset(STDERR_FILENO);
      Help(argv[0], true);
      exit(1);
    }
    // Warning - id too high
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && (bead[i] >= MoleculeType[j].nBeads ||
                                  bead[i+1] >= MoleculeType[j].nBeads)) {
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "\nWarning: ");
        ColourChange(STDERR_FILENO, CYAN);
        fprintf(stderr, "-n");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, " - index in ");
        ColourChange(STDERR_FILENO, CYAN);
        fprintf(stderr, "%d %d", bead[i], bead[i+1]);
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, " pair is larger than the number of beads in ");
        ColourChange(STDERR_FILENO, CYAN);
        fprintf(stderr, "%s", MoleculeType[j].Name);
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, " molecule;");
        fprintf(stderr, " using %d instead\n", MoleculeType[j].nBeads);
        ColourReset(STDERR_FILENO);
      }
    }
  } //}}}
  int number_of_pairs = number_of_beads / beads_per_set;
  //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // number of bins
  int bins = 1.5 / width; // orientation order is <-0.5,1>

  // array for distributions //{{{
  double *orientation[Counts.TypesOfMolecules][number_of_pairs],
         avg_orientation[Counts.TypesOfMolecules][number_of_pairs];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      orientation[i][j] = calloc(bins, sizeof *orientation[i][j]);
      avg_orientation[i][j] = 0;
    }
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  SkipVtfStructure(vcf, struct_lines);

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

  // main loop //{{{
  int count_step = 0; // count calculated timesteps
  int count_vcf = start - 1; // count timesteps from the beginning
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count_step++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}
    // read & join molecules //{{{
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // join molecules if un-joined coordinates provided
    if (!joined) {
      // transform coordinates into fractional ones for non-orthogonal box
      ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      RemovePBCMolecules(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

  // TODO: check
    // calculate orientation parameter //{{{
    // go through all molecules
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      if (MoleculeType[mtype].Use) { // use only specified molecule types
        for (int j = 0; j < number_of_beads; j += beads_per_set) {
          int pair = j / 2,
              mol_beads = MoleculeType[mtype].nBeads,
              id1 = -2, id2 = -1;
          // TODO: switch was fine; is the if else?
          // check provided ids for given molecule type
          if (test_ids(bead[j], bead[j+1], mol_beads) > 0 &&
              bead[j] != bead[j+1]) {
            id1 = Molecule[i].Bead[bead[j]];
            id2 = Molecule[i].Bead[bead[j+1]];
        //switch(test_ids(bead[j], bead[j+1], mol_beads)) { //{{{
        //  case 0: // bead[k] & bead[k+1] too high
        //    continue;
        //  case 1: // bead[k] too high & bead[k+1] okay
        //    if ((bead[j+1]+1) == mol_beads) {
        //      continue;
        //    } else {
        //      id1 = Molecule[i].Bead[mol_beads-1];
        //      id2 = Molecule[i].Bead[bead[j+1]];
        //    }
        //    break;
        //  case 2: // bead[k] okay & bead[k+1] too high
        //    if ((bead[j]+1) == mol_beads) {
        //      continue;
        //    } else {
        //      id1 = Molecule[i].Bead[bead[j]];
        //      id2 = Molecule[i].Bead[mol_beads-1];
        //    }
        //    break;
        //  case 3: // bead[k] & bead[k+1] okay
        //    id1 = Molecule[i].Bead[bead[j]];
        //    id2 = Molecule[i].Bead[bead[j+1]];
        //    break;
        //} //}}}
            // TODO: are fractionals working?
            VECTOR unit;
            unit.x = Bead[id1].Position.x - Bead[id2].Position.x;
            unit.y = Bead[id1].Position.y - Bead[id2].Position.y;
            unit.z = Bead[id1].Position.z - Bead[id2].Position.z;
            double length = Length(FromFractional(unit, Box));
            unit.x /= length;
            unit.y /= length;
            unit.z /= length;
            // cos(angle)
            double value = unit.x * normal.x +
                           unit.y * normal.y +
                           unit.z * normal.z;
            value = 0.5 * (3 * SQR(value) - 1);
            // add to overall averages
            avg_orientation[mtype][pair] += value;
            // place into appropriate distribution bin
            int k = (value + 0.5) / width; // +0.5, because its range is <-0.5,1>
            if (k < bins) {
              orientation[mtype][pair][k]++;
            } else if (k == bins) { // in case of value=1 exactly
              orientation[mtype][pair][bins-1]++;
            }
          }
        }
      }
    } //}}}

    // exit the while loop if there's no more coordinates or -e step was reached
    if (LastStep(vcf, NULL) || end == count_vcf) {
      break;
    }
  }
  fclose(vcf);
  // print last step?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // write distribution of orientation parameters //{{{
  // open output file for writing
  FILE *out = OpenFile(output, "w");
  PrintByline(out, argc, argv);
  // print first lines of output file - molecule names and beadtype pairs //{{{
  fprintf(out, "# columns: (1) Orientation order parameter");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "\n#          %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += 2) {
        // TODO: switch was fine; is the if else?
        // are ids right?
//      switch(test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads)) { //{{{
//        case 0: // bead[k] & bead[k+1] too high, skip
//          continue;
//        case 1: // bead[k] too high & bead[k+1] okay
//        case 2: // bead[k] okay & bead[k+1] too high
//        case 3: // bead[k] & bead[k+1] okay
//          if ((bead[j]+1) != bead[j]) {
//            fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
//          }
//          break;
//      } //}}}
        if (test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads) > 0 &&
            bead[j] != bead[j+1]) {
          fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
        }
      }
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.4f", width*(2*i+1)/2-0.5); // orientation order parameter: <-0.5,1>
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {
        for (int k = 0; k < number_of_beads; k += beads_per_set) {
        // TODO: switch was fine; is the if else?
          if (test_ids(bead[k], bead[k+1], MoleculeType[j].nBeads) > 0 &&
              bead[k+1] != bead[k]) {
            int pair = k / beads_per_set;
            double value = orientation[j][pair][i] /
                           (count_step * MoleculeType[j].Number);
            fprintf(out, "%10f", value);
        //switch(test_ids(bead[k], bead[k+1], MoleculeType[j].nBeads)) { //{{{
        //  case 0: // bead[k] & bead[k+1] too high
        //    continue;
        //  case 1: // bead[k] too high & bead[k+1] okay
        //    if ((bead[k+1]+1) == MoleculeType[j].nBeads) {
        //      continue;
        //    }
        //  case 2: // bead[k] okay & bead[k+1] too high
        //    if ((bead[k]+1) == MoleculeType[j].nBeads) {
        //      continue;
        //    }
        //  case 3: // bead[k] & bead[k+1] okay
        //    fprintf(out, "%10f", orientation[j][k/2][i]/(count_step*MoleculeType[j].Number));
        //    break;
        //} //}}}
          }
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write averages //{{{
  // legend line
  fprintf(out, "# Average orientation order parameter");
  count = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "\n#    %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += beads_per_set) {
        // TODO: the commented if else works; does this one?
        if (test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads) > 0 &&
            bead[j+1] != bead[j]) {
              fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
        }
//      if (bead[j] < MoleculeType[i].nBeads && bead[j+1] < MoleculeType[i].nBeads) {
//          fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
//      } else if (bead[j] >= MoleculeType[i].nBeads && bead[j+1] < MoleculeType[i].nBeads) {
//          fprintf(out, " (%d) %d-%d,", ++count, MoleculeType[i].nBeads, bead[j+1]+1);
//      } else if (bead[j+1] >= MoleculeType[i].nBeads && bead[j] < MoleculeType[i].nBeads) {
//          fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, MoleculeType[i].nBeads);
//      }
      }
    }
  }
  putc('\n', out);
  // data line
  putc('#', out);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      for (int j = 0; j < number_of_beads; j += beads_per_set) {
        // check provided ids for given molecule type
        // TODO: again, commented-out if else work, does the other one?
        if (test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads) > 0 &&
            bead[j+1] != bead[j]) {
          int pair = j / beads_per_set;
          double value = avg_orientation[i][pair] /
                         (count_step * MoleculeType[i].Number);
          fprintf(out, "%10f", value);
        }
      //int mol_beads = MoleculeType[i].nBeads;
      //if (bead[j] >= mol_beads &&
      //    bead[j+1] >= mol_beads) { // both higher than number beads - skip
      //  continue;
      //}
      //fprintf(out, "%10f", avg_orientation[i][j/2]/(count_step*MoleculeType[i].Number));
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      free(orientation[i][j]);
    }
  } //}}}

  return 0;
}
