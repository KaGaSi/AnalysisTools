#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
DihedralMolecules utility calculates angle between planes specified by beads \
in each molecule of specified molecule type(s). Each plane can be specified \
by an arbitrary trio of beads in the molecules.\n\n");
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
  fprintf(ptr, "      -n <ints>   bead indices (multiple of 6 <ints>) \
for angle calculation (default: 1 2 3 2 3 4)\n");
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
  double width = atof(argv[count]);
  // number of bins between 0 and 180 deg
  int bins = ceil(180 / width); //}}}

  // <output> - file name with angle distribution //{{{
  char output_distr[LINE] = "";
  snprintf(output_distr, LINE, "%s", argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent, verbose;
  CommonOptions_old(argc, argv, input_vsf, &verbose, &silent, LINE);
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
      int mol_type = FindMoleculeType_old(argv[count], Counts, MoleculeType);
      if (mol_type == -1) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", input_coor);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - non-existent molecule");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", argv[count]);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, "\n");
        ColourReset(STDERR_FILENO);
        ErrorMoleculeType_old(Counts, MoleculeType);
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
  int bead[100] = {0}, beads_per_angle = 6, number_of_beads = 0;
  bead[0] = 1; // default planes: 1-2-3 & 2 3 4
  bead[1] = 2;
  bead[2] = 3;
  bead[3] = 2;
  bead[4] = 3;
  bead[5] = 4;
  if (MultiIntegerOption(argc, argv, "-n", &number_of_beads, bead)) {
    exit(1);
  }
  if (number_of_beads == 0) { // -n is missing
    number_of_beads = 6;
  }
  // Error: wrong number of integers //{{{
  if ((number_of_beads%beads_per_angle) != 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-n");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - number of bead ids must be divisible by six.\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero
    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] >= MoleculeType[j].nBeads) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "-n");
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - index ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%d", bead[i]+1);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " is larger than the number of beads in molecule ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s\n\n", MoleculeType[j].Name);
        ColourReset(STDERR_FILENO);
        Help(argv[0], true);
        exit(1);
      }
    } //}}}
  }
  // Error - same bead id used in specifying a plane
  for (int i = 0; i < number_of_beads; i += (beads_per_angle/2)) {
    if (bead[i] == bead[i+1] ||
        bead[i] == bead[i+2] ||
        bead[i+1] == bead[i+2]) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "-n");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - a plane must be specified by \
three different beads (wrong trio: ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%d %d %d\n\n", bead[i], bead[i+1], bead[i+2]);
      ColourReset(STDERR_FILENO);
      exit(1);
    }
  } //}}}

  int number_of_angles = number_of_beads / beads_per_angle;

  // '-a' option - write angles for all molecules //{{{
  char output[LINE] = "";
  if (FileOption(argc, argv, "-a", output, LINE)) {
    exit(1);
  }
  // write initial stuff to output if '-a' is used
  if (output[0] != '\0') {
    // open output
    FILE *out = OpenFile(output, "w");
    PrintByline(out, argc, argv);
    // print molecule names & bead ids //{{{
    fprintf(out, "# dihedral angles between planes specifief by:");
    for (int j = 0; j < number_of_beads; j += beads_per_angle) {
      fprintf(out, " (%d) %d-%d-%d & %d-%d-%d;", j/beads_per_angle+1, bead[j]+1,
                                                 bead[j+1]+1, bead[j+2]+1,
                                                 bead[j+3]+1, bead[j+4]+1,
                                                 bead[j+5]+1);
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
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
    fprintf(stdout, "\nPlanes for which to calculate angles:\n");
    for (int i = 0; i < number_of_beads; i += beads_per_angle) {
      fprintf(stdout, "  %d-%d-%d &", bead[i], bead[i+1], bead[i+2]);
      fprintf(stdout, " %d-%d-%d\n", bead[i+3], bead[i+4], bead[i+5]);
    }
    putchar('\n');
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  SkipVtfStructure(vcf, struct_lines);

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
      RemovePBCMolecules2(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
    }
    // transform back to 'normal' coordinates for non-orthogonal box
    FromFractionalCoor(Counts.BeadsCoor, &Bead, Box); //}}}

    // calculate angles //{{{
    double angle[Counts.Molecules][number_of_angles];
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type = Molecule[i].Type;
      if (MoleculeType[mol_type].Use) {

        // calculate normal vectors to specified planes
        // first plane is given by 0 1 2 and second plane by 1 2 3
        for (int j = 0; j < number_of_beads; j += beads_per_angle) {
          int count_angle = j / beads_per_angle,
              id1 = Molecule[i].Bead[bead[j+0]],
              id2 = Molecule[i].Bead[bead[j+1]],
              id3 = Molecule[i].Bead[bead[j+2]],
              id4 = Molecule[i].Bead[bead[j+3]],
              id5 = Molecule[i].Bead[bead[j+4]],
              id6 = Molecule[i].Bead[bead[j+5]];
          // plane normals //{{{
          VECTOR u[2], v[2], n[2];
          // vectors in first plane (points 0 1 2): u=1-2; v=1-2
          u[0].x = Bead[id2].Position.x - Bead[id3].Position.x;
          u[0].y = Bead[id2].Position.y - Bead[id3].Position.y;
          u[0].z = Bead[id2].Position.z - Bead[id3].Position.z;
          v[0].x = Bead[id1].Position.x - Bead[id3].Position.x;
          v[0].y = Bead[id1].Position.y - Bead[id3].Position.y;
          v[0].z = Bead[id1].Position.z - Bead[id3].Position.z;
          // normal in first plane
          n[0].x = u[0].y * v[0].z - u[0].z * v[0].y;
          n[0].y = u[0].z * v[0].x - u[0].x * v[0].z;
          n[0].z = u[0].x * v[0].y - u[0].y * v[0].x;
          // vectors in second plane (points 3 4 5): u=4-3; v=5-3
          u[1].x = Bead[id5].Position.x - Bead[id4].Position.x;
          u[1].y = Bead[id5].Position.y - Bead[id4].Position.y;
          u[1].z = Bead[id5].Position.z - Bead[id4].Position.z;
          v[1].x = Bead[id6].Position.x - Bead[id4].Position.x;
          v[1].y = Bead[id6].Position.y - Bead[id4].Position.y;
          v[1].z = Bead[id6].Position.z - Bead[id4].Position.z;
          // normal in second plane
          n[1].x = u[1].y * v[1].z - u[1].z * v[1].y;
          n[1].y = u[1].z * v[1].x - u[1].x * v[1].z;
          n[1].z = u[1].x * v[1].y - u[1].y * v[1].x; //}}}
          // calculate angle between the two normals //{{{
          double scalar = n[0].x * n[1].x + n[0].y * n[1].y + n[0].z * n[1].z,
                 cosine = scalar / (VectorLength(n[0]) * VectorLength(n[1]));
          // too close to 180 or 0 degrees to call
          if (scalar < 0 && fabs(cosine+1) <= 0.00000001) {
            angle[i][count_angle] = PI - 0.00000001;
          } else if (scalar > 0 && fabs(cosine-1) <= 0.00000001) {
            angle[i][count_angle] = 0;
          } else { // 'normal' value
            angle[i][count_angle] = acos(cosine);
          }
          angle[i][count_angle] *= 180 / PI; // in degrees //}}}
          // add to average
          avg_angle[mol_type][count_angle] += angle[i][count_angle];
          int k = angle[i][count_angle] / width; // distribution bin
          if (k < bins) {
            distr[mol_type][count_angle][k]++;
          } else {
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stdout, "\nWarning - weird angle: ");
            ColourChange(STDERR_FILENO, CYAN);
            fprintf(stdout, "%lf degrees\n", angle[i][count_angle]);
            ColourReset(STDERR_FILENO);
          }
        }
      }
    } //}}}

    // write all angles to output if '-a' is used //{{{
    if (output[0] != '\0') {
      FILE *out = OpenFile(output, "a");
      fprintf(out, "%6d", count_vcf);
      for (int i = 0; i < Counts.Molecules; i++) {
        int mol_type = Molecule[i].Type;
        if (MoleculeType[mol_type].Use) {
          for (int j = 0; j < number_of_beads; j += beads_per_angle){
            int count_angle = j / beads_per_angle;
            // write angle between planes, not normals
            fprintf(out, " %10.6f", 180-angle[i][count_angle]);
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
  FILE *out = OpenFile(output, "w");
  PrintByline(out, argc, argv);
  // print molecule names & bead ids //{{{
  fprintf(out, "# angles between planes specifief by:");
  for (int j = 0; j < number_of_beads; j += beads_per_angle) {
    int count_angle = j / beads_per_angle;
    fprintf(out, " (%d) %d-%d-%d & %d-%d-%d;", count_angle+1, bead[j]+1,
                                               bead[j+1]+1, bead[j+2]+1,
                                               bead[j+3]+1, bead[j+4]+1,
                                               bead[j+5]+1);
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
        for (int k = 0; k < (number_of_beads/beads_per_angle); k++) {
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
      if ((number_of_beads/beads_per_angle) == 1) {
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
      for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
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
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      free(distr[i][j]);
    }
  } //}}}

  return 0;
}
