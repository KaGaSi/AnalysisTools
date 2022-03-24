#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
DensityBox utility calculates number \
density for all bead types in the direction of specified axis (x, y, or z).\
The utility works properly only for orthogonal boxes that do not change \
size \n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <axis> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>    width of a single bin\n");
  fprintf(ptr, "   <output>   output density file (automatic ending \
'<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>     calculate along x, y, or z axis\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
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
        strcmp(argv[i], "-x") != 0) {

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
  double width = atof(argv[count]); //}}}

  // <output> - filename with bead densities //{{{
  char output_rho[LINE] = "";
  snprintf(output_rho, LINE, "%s", argv[++count]); //}}}

  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {
    axis = argv[count][0];
    // Error - not x/y/z
    if (axis != 'x' && axis != 'y' && axis != 'z') {
      ErrorPrintError_old();
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "<axis>");
      ColourText(STDERR_FILENO, RED);
      fprintf(stderr, " - use 'x', 'y', or 'z'\n\n");
      ColourReset(STDERR_FILENO);
      Help(argv[0], true);
      exit(1);
    }
    // add <axis>.rho to the output filename
    char str[LINE];
    output_rho[LINE-5] = '\0'; // ensure output_rho isn't too long
    snprintf(str, LINE, "%s%c.rho", output_rho, axis);
    strcpy(output_rho, str);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
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
              &MoleculeType, &Molecule);
  // warning - works properly only for orthogonal box
  if (Box.alpha != 90 ||
      Box.beta != 90 ||
      Box.gamma != 90) {
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "\nWarning: non-orthogonal box; angles - ");
    ColourText(STDERR_FILENO, CYAN);
    fprintf(stderr, "%lf %lf %lf", Box.alpha, Box.beta, Box.gamma);
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, ".\n         %s works properly only for \
orthogonal box.\n", argv[0]);
    ColourReset(STDERR_FILENO);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // write initial stuff to output density file //{{{
  FILE *out;
  if ((out = fopen(output_rho, "w")) == NULL) {
    ErrorFileOpen(output_rho, 'w');
    exit(1);
  }
  PrintByline(out, argc, argv);
  // print bead type names to output file
  fprintf(out, "# columns: (1) distance;");
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(out, " (%d) %s", i+2, BeadType[i].Name);
    if (i != (Counts.TypesOfBeads-1)) {
      putc(';', out);
    }
  }
  putc('\n', out);
  fclose(out); //}}}

  // number of bins //{{{
  double max_dist;
  switch(axis) {
    case 'x':
      max_dist = Box.Length.x;
      break;
    case 'y':
      max_dist = Box.Length.y;
      break;
    case 'z':
      max_dist = Box.Length.z;
      break;
  }
  int bins = ceil(max_dist / width); //}}}

// TODO: sizeof ...argh!
  // allocate memory for density arrays //{{{
//long int **rho = malloc(Counts.TypesOfBeads*sizeof(long int *));
//long int **rho_2 = malloc(Counts.TypesOfBeads*sizeof(long int *));
  long int **rho = malloc(Counts.TypesOfBeads * sizeof **rho);
  long int **rho_2 = malloc(Counts.TypesOfBeads * sizeof **rho_2);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  rho[i] = calloc(bins, sizeof(long int));
//  rho_2[i] = calloc(bins, sizeof(long int));
    rho[i] = calloc(bins, sizeof *rho[i]);
    rho_2[i] = calloc(bins, sizeof *rho_2[i]);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

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
    // read coordinates and wrap box - assumes orthogonal box
    BOX test = InitBox; // check the box didn't change
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // warning - box size/angles //{{{
    if (Box.Length.x != test.Length.x ||
        Box.Length.y != test.Length.y ||
        Box.Length.z != test.Length.z) {
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: box side lengths changed from ");
      ColourText(STDERR_FILENO, CYAN);
      fprintf(stderr, "%lf %lf %lf", Box.Length.x, Box.Length.y, Box.Length.z);
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, " to ");
      ColourText(STDERR_FILENO, CYAN);
      fprintf(stderr, "%lf %lf %lf", test.Length.x, test.Length.y,
                                     test.Length.z);
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, ".\n         %s works properly only when the box \
size is constant.\n", argv[0]);
      ColourReset(STDERR_FILENO);
    }
    if (Box.alpha != 90 ||
        Box.beta != 90 ||
        Box.gamma != 90) {
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: non-orthogonal box; angles - ");
      ColourText(STDERR_FILENO, CYAN);
      fprintf(stderr, "%lf %lf %lf", Box.alpha, Box.beta, Box.gamma);
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, ".\n         %s works properly only for \
orthogonal box.\n", argv[0]);
      ColourReset(STDERR_FILENO);
    } //}}}
    RestorePBC(Counts.BeadsCoor, Box, &Bead);

  // TODO: sizeof ...argh!
    // allocate memory for temporary density arrays //{{{
//  int **temp_rho = malloc(Counts.TypesOfBeads*sizeof(int *));
    int **temp_rho = malloc(Counts.TypesOfBeads * sizeof **temp_rho);
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
//    temp_rho[i] = calloc(bins,sizeof(int));
      temp_rho[i] = calloc(bins, sizeof temp_rho[i]);
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Counts.BeadsCoor; i++) {
      bool use = true;
      int mol = Bead[i].Molecule;
      if (mol != -1) { // do not use excluded molecules (-x option)
        int mtype = Molecule[mol].Type;
        use = MoleculeType[mtype].Use;
      }
      if (use) {
        if (axis == 'x') {
          int j = Bead[i].Position.x / width;
          temp_rho[Bead[i].Type][j]++;
        } else if (axis == 'y') {
          int j = Bead[i].Position.y / width;
          temp_rho[Bead[i].Type][j]++;
        } else {
          int j = Bead[i].Position.z / width;
          temp_rho[Bead[i].Type][j]++;
        }
      }
    } //}}}

    // add from temporary density array to global density arrays //{{{
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < bins; k++) {
        rho[j][k] += temp_rho[j][k];
        rho_2[j][k] += SQR(temp_rho[j][k]);
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      free(temp_rho[i]);
    }
    free(temp_rho); //}}}

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

  // write densities to output file(s) //{{{
  if ((out = fopen(output_rho, "a")) == NULL) {
    ErrorFileOpen(output_rho, 'a');
    exit(1);
  }

  // calculate rdf
  double volume = width;
  if (axis == 'x') {
    volume *= Box.Length.y * Box.Length.z;
  } else if (axis == 'y') {
    volume *= Box.Length.x * Box.Length.z;
  } else {
    volume *= Box.Length.x * Box.Length.y;
  }
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.3f", width*(2*i+1)/2);
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      double temp_rho = rho[j][i] / (volume * count);
      // print average value to output file
      fprintf(out, " %10f", temp_rho);
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
