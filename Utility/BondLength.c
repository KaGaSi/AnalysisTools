#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
BondLength utility calculates distribution of bond lengths for all bead type \
pairs in specified molecule type(s). It can also calculate distribution of \
distances between any two beads in those molecule type(s).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>    width of a single distribution bin\n");
  fprintf(ptr, "   <output>   output file with the distribution of \
bond lengths\n");
  fprintf(ptr, "   <mol(s)>   molecule name(s) to calculate bond lengths for \
(optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --all             use all molecules \
(overwrites <mol(s)>)\n");
  fprintf(ptr, "      --joined          specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -st <int>         starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>          ending timestep for calculation\n");
  fprintf(ptr, "      -d <out> [ints]   write distribution of distances \
between specified beads to file <out> (if no [ints] are \
provided, the molecule's first and last bead are used)\n");
  fprintf(ptr, "      -w <float>        warn if the length exceeds <float> \n");
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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-d") != 0 &&
        strcmp(argv[i], "-w") != 0) {

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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal_old(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with bond length distribution //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions_old(argc, argv, input_vsf, &verbose, &silent, LINE);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

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

  // <mol(s)> - names of molecule types to use //{{{
  if (!all) { // --all option not used
    while (++count < argc && argv[count][0] != '-') {
      int type = FindMoleculeType_old(argv[count], Counts, MoleculeType);
      // error - nonexistent molecule  //{{{
      if (type == -1) {
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
        ErrorMoleculeType_old(Counts, MoleculeType);
        exit(1);
      } //}}}
      MoleculeType[type].Use = true;
    }
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      MoleculeType[i].Use = true;
    }
  } //}}}

  // '-d' option - specify bead ids to calculate distance between //{{{
  int bead[100] = {0},
      number_of_beads, // number of parameters to -d
      beads_per_set = 2; // the numbers must come in pairs
  char output_d[LINE] = "";
  if (FileIntsOption(argc, argv, "-d", bead, &number_of_beads, output_d)) {
    exit(1);
  }
  // if '-d' is present, but without numbers - use first and last for each molecule
  if (output_d[0] != '\0' && number_of_beads == 0) {
    number_of_beads = 2;
    bead[0] = 1;
    bead[1] = 1000000;
  }
  // Error: wrong number of integers //{{{
  if (output_d[0] != '\0' && (number_of_beads%beads_per_set) != 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-d");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - number of bead ids must be even\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < number_of_beads; i += 2) {
    if (bead[i] == bead[i+1] || bead[i] == 0 || bead[i+1] == 0) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "-d");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - the bead indices must be non-zero and different\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    }
  } //}}}
  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero (or is -1 if none specified)
  } //}}}

  int number_of_pairs = number_of_beads / beads_per_set;

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // '-w' option - bond length warning //{{{
  double warn = -1;
  if (DoubleOption(argc, argv, "-w", &warn)) {
    exit(1);
  } //}}}

  // number of bins
  int bins = Max3(Box.Length.x, Box.Length.y, Box.Length.z) / (2 * width);

  // arrays for distributions //{{{
  double *length[Counts.TypesOfMolecules]
                [Counts.TypesOfBeads][Counts.TypesOfBeads];
  double min_max[Counts.TypesOfMolecules]
                [Counts.TypesOfBeads][Counts.TypesOfBeads][2];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        length[i][j][k] = calloc(bins, sizeof *length[i][j][k]);
        min_max[i][j][k][0] = 10 * Max3(Box.Length.x,
                                        Box.Length.y, Box.Length.z);
        min_max[i][j][k][1] = 0;
      }
    }
  }
  // extra arrays for -d option
  double *distance[Counts.TypesOfMolecules][number_of_pairs];
  double min_max_d_option[Counts.TypesOfMolecules][number_of_pairs][2];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      distance[i][j] = calloc(bins, sizeof *distance[i][j]);
      min_max_d_option[i][j][0] = 10 * Max3(Box.Length.x,
                                            Box.Length.y, Box.Length.z);
      min_max_d_option[i][j][1] = 0;
    }
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  SkipVtfStructure(vcf, struct_lines);

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

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
    // read & join molecules //{{{
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // join molecules if un-joined coordinates provided
    if (!joined) {
      // transform coordinates into fractional ones for non-orthogonal box
      ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      RemovePBCMolecules2(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
      // transform back to 'normal' coordinates for non-orthogonal box
      FromFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    } //}}}
    // calculate bond lengths //{{{
    // go through all molecules
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      if (MoleculeType[mtype].Use) { // use only specified molecule types
        for (int j = 0; j < MoleculeType[mtype].nBonds; j++) {

          // bead ids in the bond
          int id1 = Molecule[i].Bead[MoleculeType[mtype].Bond[j][0]];
          int id2 = Molecule[i].Bead[MoleculeType[mtype].Bond[j][1]];
          int btype1 = Bead[id1].Type;
          int btype2 = Bead[id2].Type;

          // bond length //{{{
          VECTOR bond;
          bond.x = Bead[id1].Position.x - Bead[id2].Position.x;
          bond.y = Bead[id1].Position.y - Bead[id2].Position.y;
          bond.z = Bead[id1].Position.z - Bead[id2].Position.z;
          bond.x = Length(bond); //}}}

          // warn if bond is too long //{{{
          if (warn > -1 && bond.x > warn) {
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "\nWarning (-w option): long bond (");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%lf", bond.x);
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, " )\n   Step: ");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%d", count_vcf);
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "\n   Beads: ");
            // first bead
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%6d", Bead[id1].Indexxx);
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, " (");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%s", BeadType[btype1].Name);
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "): ");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%lf %lf %lf\n", Bead[id1].Position.x,
                                             Bead[id1].Position.y,
                                             Bead[id1].Position.z);
            // second bead
            fprintf(stderr, "          %6d", Bead[id1].Indexxx);
            fprintf(stderr, " (");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%s", BeadType[btype1].Name);
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "): ");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stderr, "%lf %lf %lf\n", Bead[id1].Position.x,
                                             Bead[id1].Position.y,
                                             Bead[id1].Position.z);
            ColourReset(STDERR_FILENO);
          } //}}}

          // btype1 must be lower then btype2
          if (btype1 > btype2) {
            SwapInt(&btype1, &btype2);
          }

          // mins & maxes //{{{
          if (bond.x < min_max[mtype][btype1][btype2][0]) {
            min_max[mtype][btype1][btype2][0] = bond.x;
          } else if (bond.x > min_max[mtype][btype1][btype2][1]) {
            min_max[mtype][btype1][btype2][1] = bond.x;
          } //}}}

          int k = bond.x / width;
          if (k < bins) {
            length[mtype][btype1][btype2][k]++;
          }
        }
      }
    } //}}}

    // calculate distance (-d option) //{{{
    if (output_d[0] != '\0') {
      for (int i = 0; i < Counts.Molecules; i++) {
        int mtype = Molecule[i].Type;
        if (MoleculeType[mtype].Use) { // use only specified molecule types
          for (int j = 0; j < number_of_beads; j += beads_per_set) {
            // bead ids to use //{{{
            int id1, id2;
            // use first molecule bead if bead index too high or -1
            if (bead[j] == -1 || bead[j] >= MoleculeType[mtype].nBeads) {
              id1 = Molecule[i].Bead[0];
            } else{ // use specified index otherwise
              id1 = Molecule[i].Bead[bead[j]];
            }
            // use last molecule bead if bead index too high or -1
            if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[mtype].nBeads) {
              id2 = Molecule[i].Bead[MoleculeType[mtype].nBeads-1];
            } else{ // use specified index otherwise
              id2 = Molecule[i].Bead[bead[j+1]];
            } //}}}
            // distance calculation //{{{
            VECTOR dist;
            dist.x = Bead[id1].Position.x - Bead[id2].Position.x;
            dist.y = Bead[id1].Position.y - Bead[id2].Position.y;
            dist.z = Bead[id1].Position.z - Bead[id2].Position.z;
            dist.x = Length(dist); //}}}
            // distance mins & maxes //{{{
            int pair = j / 2;
            if (dist.x < min_max_d_option[mtype][pair][0]) {
              min_max_d_option[mtype][pair][0] = dist.x;
            } else if (dist.x > min_max_d_option[mtype][pair][1]) {
              min_max_d_option[mtype][pair][1] = dist.x;
            } //}}}
            int k = dist.x / width; // distribution 'bin'
            if (k < bins) {
              distance[mtype][pair][k]++;
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

  // count number of calculated steps //{{{
  int steps;
  if (end != -1) {
    steps = end - start + 1;
  } else {
    steps = count_vcf - start + 1;
  } //}}}

  // count total number of bonds in molecules //{{{
  int bonds[Counts.TypesOfMolecules][Counts.TypesOfBeads][Counts.TypesOfBeads];
  // zeroize the array
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] = 0;
        }
      }
    }
  }
  // count the bonds
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] += length[i][j][k][l];
        }
      }
    }
  } //}}}

  // write distribution of bond lengths //{{{
  FILE *out = OpenFile(output, "w");
  PrintByline(out, argc, argv);
  // print first line of output file - molecule names and beadtype pairs //{{{
  fprintf(out, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, " %s molecule:", MoleculeType[i].Name);
      for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < MoleculeType[i].nBTypes; k++) {
          int btype1 = MoleculeType[i].BType[j];
          int btype2 = MoleculeType[i].BType[k];
          if (bonds[i][btype1][btype2] > 0) {
            fprintf(out, " (%d) %s-%s", count, BeadType[btype1].Name,
                                               BeadType[btype2].Name);
            count++;
            // add semicolon for molecule's last pair; add comma otherwise
            if (k == (MoleculeType[i].nBTypes-1)) {
              putc(';', out);
            } else {
              putc(',', out);
            }
          }
        }
      }
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.4f", width*(2*i+1)/2);
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {
        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < MoleculeType[j].nBTypes; k++) {
          for (int l = k; l < MoleculeType[j].nBTypes; l++) {
            int btype1 = MoleculeType[j].BType[k];
            int btype2 = MoleculeType[j].BType[l];
            // btype1 must be lower than btype2
            if (btype1 > btype2) {
              SwapInt(&btype1, &btype2);
            }
            if (bonds[j][btype1][btype2] > 0) {
              double value = length[j][btype1][btype2][i] /
                             bonds[j][btype1][btype2];
              fprintf(out, "%10f", value);
            }
          }
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write mins and maxes //{{{
  // legend line
  fprintf(out, "# mins(odd columns)/maxes(even columns) -");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, " %s molecule:", MoleculeType[i].Name);
      for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < MoleculeType[i].nBTypes; k++) {
          int btype1 = MoleculeType[i].BType[j];
          int btype2 = MoleculeType[i].BType[k];
          if (bonds[i][btype1][btype2] > 0) {
            fprintf(out, " (%d) %s-%s", count, BeadType[btype1].Name,
                                               BeadType[btype2].Name);
            count += 2;
            if (k == (MoleculeType[i].nBTypes-1)) {
              if (i != (Counts.TypesOfMolecules-1)) {
                putc(';', out); // molecule's last bead type pair
              }
            } else {
              putc(',', out); // not the last pair
            }
          }
        }
      }
    }
  }
  putc('\n', out);
  // data line
  putc('#', out);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        if (min_max[i][j][k][1] > 0) {
          fprintf(out, " %lf %lf", min_max[i][j][k][0], min_max[i][j][k][1]);
        }
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  // write distribution of distances from '-d' option //{{{
  if (output_d[0] != '\0') {
    out = OpenFile(output_d, "w");
    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);
    // print the first line - molecule names with bead order //{{{
    fprintf(out, "# bead order in molecule(s) -");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s:", MoleculeType[i].Name);
        for (int j = 0; j < MoleculeType[i].nBeads; j++) {
          int btype = MoleculeType[i].Bead[j];
          fprintf(out, " %s", BeadType[btype].Name);
        }
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    // print the second line of output file - molecule names and indices with column numbers //{{{
    fprintf(out, "# (1) distance;");
    count = 1;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s:", MoleculeType[i].Name);

        for (int j = 0; j < number_of_beads; j += beads_per_set) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= MoleculeType[i].nBeads) {
            id1 = 1;
          } else{ // use specified index otherwise
            id1 = bead[j]+1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[i].nBeads) {
            id2 = MoleculeType[i].nBeads;
          } else{ // use specified index otherwise
            id2 = bead[j+1]+1;
          } //}}}

          fprintf(out, " (%d) %d-%d", ++count, id1, id2);
          // add semicolon for the molecule's last pair, add comma otherwise
          if (j == (number_of_beads-beads_per_set)) {
            putc(';', out);
          } else {
            putc(',', out);
          }
        }
      }
    }
    putc('\n', out); //}}}

    // write the distribution to output file //{{{
    for (int i = 0; i < bins; i++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2);
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        if (MoleculeType[j].Use) {
          for (int k = 0; k < number_of_beads; k += beads_per_set) {
            double value = (double)(distance[j][k/2][i]) /
                           (steps * MoleculeType[j].Number);
            fprintf(out, " %10f", value);
          }
        }
      }
      putc('\n', out);
    } //}}}

    // write mins and maxes //{{{
    // legend line
    fprintf(out, "# mins(odd columns)/maxes(even columns) -");
    count = 1;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s:", MoleculeType[i].Name);
        for (int j = 0; j < number_of_beads; j += 2) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= MoleculeType[i].nBeads) {
            id1 = 1;
          } else{ // use specified index otherwise
            id1 = bead[j]+1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[i].nBeads) {
            id2 = MoleculeType[i].nBeads;
          } else{ // use specified index otherwise
            id2 = bead[j+1]+1;
          } //}}}

          fprintf(out, " (%d) %d-%d", count, id1, id2);
          count += 2;
          // add semicolon if this is the last pair for this molecule, add comma otherwise
          if (j == (number_of_beads-2)) {
            putc(';', out);
          } else {
            putc(',', out);
          }
        }
      }
    }
    putc('\n', out);

    // data line
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        for (int j = 0; j < number_of_beads; j += beads_per_set) {
          int pair = j / 2;
          fprintf(out, " %lf %lf", min_max_d_option[i][pair][0],
                                   min_max_d_option[i][pair][1]);
        }
      }
    }
    putc('\n', out); //}}}

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        free(length[i][j][k]);
      }
    }
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      free(distance[i][j]);
    }
  } //}}}

  return 0;
}
