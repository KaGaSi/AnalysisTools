#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
AngleMolecules utility calculates angles for specified beads (three beads per \
angle) for molecules of specified type(s). The specified beads do not have \
to be connected by bonds.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol(s)> [options]\n\n", cmd);
  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>    width of a single bin in degrees\n");
  fprintf(ptr, "   <output>   output file with distribution of angles\n");
  fprintf(ptr, "   <mol(s)>   molecule name(s) to calculate angles for \
(optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --all       use all molecules (overwrites <mol(s)>)\n");
  fprintf(ptr, "      --joined    specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -n <ints>   bead indices (multiple of 3 <ints>) \
for angle calculation (default: 1 2 3)\n");
  fprintf(ptr, "      -a <out>    write angles of all molecules in \
all timesteps to <out>\n");
  fprintf(ptr, "      -st <int>   starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>    ending timestep for calculation\n");
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
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-a") != 0 &&
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
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]);
  // number of bins between 0 and 180 deg
  int bins = ceil(180 / width); //}}}

  // <output> - file name with angle distribution //{{{
  char output_distr[LINE] = "";
  snprintf(output_distr, LINE, "%s", argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent, verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");
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
  BOX Box = InitBox; // box dimensions
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // <mol(s)> - names of molecule types to use //{{{
  if (!all) { // --all option not used
    while (++count < argc && argv[count][0] != '-') {
      int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);
      if (mol_type == -1) {
        ErrorPrintError_old();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - non-existent molecule");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", argv[count]);
        RedText(STDERR_FILENO);
        fprintf(stderr, "\n");
        ResetColour(STDERR_FILENO);
        ErrorMoleculeType(Counts, MoleculeType);
        exit(1);
      } else {
        MoleculeType[mol_type].Use = true;
      }
    }
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      MoleculeType[i].Use = true;
    }
  } //}}}

  // '-n' option - specify bead ids //{{{
  int bead[100] = {0}, beads_per_set = 3, number_of_beads = 0;
  bead[0] = 1; // default ids for angle
  bead[1] = 2;
  bead[2] = 3;
  if (MultiIntegerOption(argc, argv, "-n", &number_of_beads, bead)) {
    exit(1);
  }
  if (number_of_beads == 0) { // -n is missing
    number_of_beads = 3;
  }
  // Error: wrong number of integers //{{{
  if ((number_of_beads%beads_per_set) != 0) {
    ErrorPrintError_old();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "-n");
    RedText(STDERR_FILENO);
    fprintf(stderr, " - number of bead ids must be divisible by three\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero
    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] >= MoleculeType[j].nBeads) {
        ErrorPrintError_old();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "-n");
        RedText(STDERR_FILENO);
        fprintf(stderr, " - index ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%d", bead[i]+1);
        RedText(STDERR_FILENO);
        fprintf(stderr, " is larger than the number of beads in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n\n", MoleculeType[j].Name);
        ResetColour(STDERR_FILENO);
        Help(argv[0], true);
        exit(1);
      }
    } //}}}
  } //}}}

  int number_of_angles = number_of_beads / beads_per_set;

  // '-a' option - write angles for all molecules //{{{
  char output[LINE] = "";
  if (FileOption(argc, argv, "-a", output, LINE)) {
    exit(1);
  }
  // write initial stuff to output if '-a' is used
  if (output[0] != '\0') {
    // open output //{{{
    FILE *out;
    if ((out = fopen(output, "w")) == NULL) {
      ErrorFileOpen(output, 'w');
      exit(1);
    } //}}}
    PrintByline(out, argc, argv);
    // print molecule names & bead ids //{{{
    fprintf(out, "# angles between beads:");
    for (int j = 0; j < number_of_beads; j += beads_per_set) {
      fprintf(out, " (%d) %d-%d-%d;", (j/beads_per_set)+1, bead[j]+1,
                                                             bead[j+1]+1,
                                                             bead[j+2]+1);
    }
    putc('\n', out);
    fprintf(out, "# columns: (1) step;");
    int j = 2;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        int angles_per_mtype = MoleculeType[i].Number * number_of_angles;
        if (angles_per_mtype == 1) {
          fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
        } else {
          int high = j + angles_per_mtype - 1;
          fprintf(out, " (%d) to (%d) %s molecules;", j, high,
                                                      MoleculeType[i].Name);
        }
        j += angles_per_mtype;
      }
    }
    putc('\n', out); //}}}
    fclose(out);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
    fprintf(stdout, "\nAngle-specifying bead trios:\n");
    for (int i = 0; i < number_of_beads; i += beads_per_set) {
      fprintf(stdout, "  %d-%d-%d\n", bead[i], bead[i+1], bead[i+2]);
    }
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

  // allocate array for distribution of angles //{{{
  double avg_angle[Counts.TypesOfMolecules][number_of_angles];
  double *distr[Counts.TypesOfMolecules][number_of_angles];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_angles; j++) {
      avg_angle[i][j] = 0;
      distr[i][j] = calloc(bins, sizeof *distr[i][j]);
    }
  } //}}}

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
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    // join molecules if un-joined coordinates provided
    if (!joined) {
      RemovePBCMolecules(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
    }
    // transform back to 'normal' coordinates for non-orthogonal box
    FromFractionalCoor(Counts.BeadsCoor, &Bead, Box); //}}}

    // calculate angles //{{{
    double angle[Counts.Molecules][number_of_angles];
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type = Molecule[i].Type;
      if (MoleculeType[mol_type].Use) {
        for (int j = 0; j < number_of_beads; j += beads_per_set) {
          int count_angle = j / beads_per_set,
              id1 = Molecule[i].Bead[bead[j+0]],
              id2 = Molecule[i].Bead[bead[j+1]],
              id3 = Molecule[i].Bead[bead[j+2]];
          VECTOR u, v;
          // first vector
          u.x = Bead[id1].Position.x - Bead[id2].Position.x;
          u.y = Bead[id1].Position.y - Bead[id2].Position.y;
          u.z = Bead[id1].Position.z - Bead[id2].Position.z;
          // second vector
          v.x = Bead[id3].Position.x - Bead[id2].Position.x;
          v.y = Bead[id3].Position.y - Bead[id2].Position.y;
          v.z = Bead[id3].Position.z - Bead[id2].Position.z;
          // calculate angle between the two vectors
          double scalar = u.x * v.x + u.y * v.y + u.z * v.z;
          angle[i][count_angle] = acos(scalar / (Length(u) * Length(v)));
          angle[i][count_angle] *= 180 / PI; // from radians to degrees

          // add to average
          avg_angle[mol_type][count_angle] += angle[i][count_angle];
          int k = angle[i][count_angle] / width; // distribution bin
          if (k < bins) {
            distr[mol_type][count_angle][k]++;
          } else {
            YellowText(STDERR_FILENO);
            fprintf(stdout, "\nWarning - weird angle: ");
            CyanText(STDERR_FILENO);
            fprintf(stdout, "%lf degrees\n", angle[i][count_angle]);
            ResetColour(STDERR_FILENO);
          }
        }
      }
    } //}}}

    // write all angles to output if '-a' is used //{{{
    if (output[0] != '\0') {
      FILE *out;
      if ((out = fopen(output, "a")) == NULL) {
        ErrorFileOpen(output, 'a');
        exit(1);
      }
      fprintf(out, "%6d", count_vcf);
      for (int i = 0; i < Counts.Molecules; i++) {
        int mol_type = Molecule[i].Type;
        if (MoleculeType[mol_type].Use) {
          for (int j = 0; j < number_of_beads; j += beads_per_set){
            int count_angle = j / beads_per_set;
            fprintf(out, " %10.6f", angle[i][count_angle]);
          }
        }
      }
      putc('\n', out);
      fclose(out);
    } //}}}

    // exit the while loop if there's no more coordinates or -e step was reached
    if (LastStep(vcf, NULL) || end == count_vcf) {
      break;
    }
  }
  fclose(vcf);
  // print last step? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}
  //}}}

  // write distribution of angles //{{{
  // open output file for appending //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    ErrorFileOpen(output_distr, 'w');
    exit(1);
  } //}}}
  PrintByline(out, argc, argv);
  // print molecule names & bead ids //{{{
  fprintf(out, "# angles between beads:");
  for (int j = 0; j < number_of_beads; j += beads_per_set) {
    int count_angle = j / beads_per_set;
    fprintf(out, " (%d) %d-%d-%d;", count_angle+1, bead[j]+1,
                                    bead[j+1]+1, bead[j+2]+1);
  }
  putc('\n', out);
  fprintf(out, "# columns: (1) angle [deg];");
  int j = 2;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if (number_of_angles == 1) {
        fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_angles-1,
                                                    MoleculeType[i].Name);
      }
      j += number_of_angles;
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%5.1f", width*(2*i+1)/2);
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {
        for (int k = 0; k < (number_of_beads/beads_per_set); k++) {
          double value = (double)(distr[j][k][i]) /
                         (count_step * MoleculeType[j].Number);
          fprintf(out, "%10f", value);
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write to output average angles //{{{
  fprintf(out, "# simple averages:");
  j = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_set) == 1) {
        fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_angles-1,
                                                    MoleculeType[i].Name);
      }
      j += number_of_angles;
    }
  }
  fprintf(out, "\n#");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      for (int j = 0; j < (number_of_beads/beads_per_set); j++) {
        double value = avg_angle[i][j] / (count_step * MoleculeType[i].Number);
        fprintf(out, " %7.3f", value);
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out);
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_set); j++) {
      free(distr[i][j]);
    }
  } //}}}

  return 0;
}
