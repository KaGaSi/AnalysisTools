#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Surface utility determines the first bead (either from the box's centre or \
edges) in each square prism defined by the given <width> parameter, thus \
defining a surface of, e.g., polymer brush or lipid bilayer. The <width> \
slices the box into square prisms along the chosen axis (i.e., if z is the \
chosen axis, the xy plane is chopped into squares, creating \
<width>*<width>*<box length in z> prisms). In each such prism, two beads \
are found corresponding to the two surfaces (e.g., polymer brush on both box \
edges or the two surfaces of a lipid bilayer inside the box).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <axis> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>    width of a single bin\n");
  fprintf(ptr, "   <output>   output file\n");
  fprintf(ptr, "   <axis>     calculate along x, y, or z axis\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -in             start from the box's centre \
instead of edges\n");
  fprintf(ptr, "      -m <mol(s)>     molecule type(s) to use\n");
  fprintf(ptr, "      -bt <name(s)>   bead type(s) to use\n");
  fprintf(ptr, "      -st <int>       starting timestep for calculation\n");
  fprintf(ptr, "      -e <int>        number of timestep to end with\n");
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
        strcmp(argv[i], "-in") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
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
  double width = atof(argv[count]); //}}}

  // <output> - output filename //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]); //}}}

  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {
    axis = argv[count][0];
    // Error - not x/y/z
    if (axis != 'x' && axis != 'y' && axis != 'z') {
      ErrorPrintError_old();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "<axis>");
      RedText(STDERR_FILENO);
      fprintf(stderr, " - use 'x', 'y', or 'z'\n\n");
      ResetColour(STDERR_FILENO);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // go from box edges to the box centre instead of from centre to edges //{{{
  bool in = BoolOption(argc, argv, "-in"); //}}}
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
    YellowText(STDERR_FILENO);
    fprintf(stderr, "\nWarning: non-orthogonal box; angles - ");
    CyanText(STDERR_FILENO);
    fprintf(stderr, "%lf %lf %lf", Box.alpha, Box.beta, Box.gamma);
    YellowText(STDERR_FILENO);
    fprintf(stderr, ".\n         %s works properly only for \
orthogonal box.\n", argv[0]);
    ResetColour(STDERR_FILENO);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", true, Counts, &BeadType)) {
    exit(1);
  } //}}}

  // -m <name(s)> - specify what molecule types to use //{{{
//int use = malloc(Counts.TypesOfMolecules * sizeof(int));
  int use[Counts.TypesOfMolecules];
  // if -m not present, use all
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    use[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", use, Counts, &MoleculeType)) {
    exit(1);
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Use = use[i];
  } //}}}

  // set maximum/minimum as half a box length in the given direction //{{{
  double range[2];
  switch(axis) {
    case 'x':
      range[0] = 0;
      range[1] = Box.Length.x;
      break;
    case 'y':
      range[0] = 0;
      range[1] = Box.Length.y;
      break;
    case 'z':
      range[0] = 0;
      range[1] = Box.Length.z;
      break;
  } //}}}

  // number of bins //{{{
  int bins[2];
  switch(axis) {
    case 'x':
      bins[0] = Box.Length.y / width + 1;
      bins[1] = Box.Length.z / width + 1;
      break;
    case 'y':
      bins[0] = Box.Length.x / width + 1;
      bins[1] = Box.Length.z / width + 1;
      break;
    case 'z':
      bins[0] = Box.Length.x / width + 1;
      bins[1] = Box.Length.y / width + 1;
      break;
  }
  //}}}

// TODO: sizeof ...argh!
  // allocate memory for density arrays //{{{
//double ***surf = calloc(bins[0], sizeof(double **)); // sum
  double ***surf = calloc(bins[0], sizeof ***surf); // sum
//int ***values = calloc(bins[0], sizeof(int **)); // number of values
  int ***values = calloc(bins[0], sizeof ***values); // number of values
  for (int i = 0; i < bins[0]; i++) {
//  surf[i] = calloc(bins[1], sizeof(double *));
    surf[i] = calloc(bins[1], sizeof **surf[i]);
//  values[i] = calloc(bins[0], sizeof(int *));
    values[i] = calloc(bins[0], sizeof **values[i]);
    for (int j = 0; j < bins[1]; j++) {
//    surf[i][j] = calloc(2, sizeof(double));
      surf[i][j] = calloc(2, sizeof *surf[i][j]);
//    values[i][j] = calloc(2, sizeof(int));
      values[i][j] = calloc(2, sizeof *values[i][j]);
    }
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
    // read coordinates and wrap box - assumes orthogonal box
    BOX test = InitBox; // check the box didn't change
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // warning - box size/angles //{{{
    if (Box.Length.x != test.Length.x ||
        Box.Length.y != test.Length.y ||
        Box.Length.z != test.Length.z) {
      YellowText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: box side lengths changed from ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%lf %lf %lf", Box.Length.x, Box.Length.y, Box.Length.z);
      YellowText(STDERR_FILENO);
      fprintf(stderr, " to ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%lf %lf %lf", test.Length.x, test.Length.y,
                                     test.Length.z);
      YellowText(STDERR_FILENO);
      fprintf(stderr, ".\n         %s works properly only when the box \
size is constant.\n", argv[0]);
      ResetColour(STDERR_FILENO);
    }
    if (Box.alpha != 90 ||
        Box.beta != 90 ||
        Box.gamma != 90) {
      YellowText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: non-orthogonal box; angles - ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%lf %lf %lf", Box.alpha, Box.beta, Box.gamma);
      YellowText(STDERR_FILENO);
      fprintf(stderr, ".\n         %s works properly only for \
orthogonal box.\n", argv[0]);
      ResetColour(STDERR_FILENO);
    } //}}}
    RestorePBC(Counts.BeadsCoor, Box, &Bead);

  // TODO: sizeof ...argh!
    // allocate memory for temporary arrays //{{{
//  double ***temp = calloc(bins[0], sizeof(double **));
//  TODO: change to (**temp)[2]
    double ***temp = calloc(bins[0], sizeof ***temp);
    for (int i = 0; i < bins[0]; i++) {
//    temp[i] = calloc(bins[1], sizeof(double *));
      temp[i] = calloc(bins[1], sizeof **temp[i]);
      for (int j = 0; j < bins[1]; j++) {
//      temp[i][j] = calloc(2, sizeof(double));
        temp[i][j] = calloc(2, sizeof *temp[i][j]);
        if (!in) {
          temp[i][j][0] = 0;
          switch(axis) {
            case 'x':
              temp[i][j][1] = Box.Length.x;
              break;
            case 'y':
              temp[i][j][1] = Box.Length.y;
              break;
            case 'z':
              temp[i][j][1] = Box.Length.z;
              break;
          }
        } else {
          switch(axis) {
            case 'x':
              temp[i][j][0] = Box.Length.x;
              break;
            case 'y':
              temp[i][j][0] = Box.Length.y;
              break;
            case 'z':
              temp[i][j][0] = Box.Length.z;
              break;
          }
          temp[i][j][1] = 0;
        }
      }
    } //}}}

  // TODO: check
    // calculate surface //{{{
    for (int i = 0; i < Counts.BeadsCoor; i++) {
      int btype = Bead[i].Type,
          mol = Bead[i].Molecule;
      if (mol == -1) { // consider only beads in molecules
        continue;
      }
      int mtype = Molecule[mol].Type;
      if (Bead[i].Molecule != -1 && MoleculeType[mtype].Use && BeadType[btype].Use) {
        double coor[3];
        switch(axis) {
          case 'x':
            coor[0] = Bead[i].Position.y;
            coor[1] = Bead[i].Position.z;
            coor[2] = Bead[i].Position.x;
            break;
          case 'y':
            coor[0] = Bead[i].Position.x;
            coor[1] = Bead[i].Position.z;
            coor[2] = Bead[i].Position.y;
            break;
          case 'z':
            coor[0] = Bead[i].Position.x;
            coor[1] = Bead[i].Position.y;
            coor[2] = Bead[i].Position.z;
            break;
        }
        int bin[2];
        bin[0] = coor[0] / width;
        bin[1] = coor[1] / width;
        if (!in) { // go from box centre to edges
          if (coor[2] >= temp[bin[0]][bin[1]][0] && coor[2] >= range[0] && coor[2] <= ((range[0]+range[1])/2)) {
            temp[bin[0]][bin[1]][0] = coor[2];
          }
          if (coor[2] <= temp[bin[0]][bin[1]][1] && coor[2] >= ((range[0]+range[1])/2) && coor[2] <= range[1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
          }
        } else if (coor[2] >= range[0] && coor[2] <= range[1]) { // go from box edges to centre
          if (coor[2] <= temp[bin[0]][bin[1]][0]) {
            temp[bin[0]][bin[1]][0] = coor[2];
          }
          if (coor[2] >= temp[bin[0]][bin[1]][1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
          }
        }
      }
    } //}}}

    // add to sums //{{{
    for (int i = 0; i < bins[0]; i++) {
      for (int j = 0; j < bins[1]; j++) {
        if (!in) {
          if (temp[i][j][0] > 0) {
            surf[i][j][0] += temp[i][j][0];
            values[i][j][0]++;
          }
          switch(axis) {
            case 'x':
              if (temp[i][j][1] < Box.Length.x) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'y':
              if (temp[i][j][1] < Box.Length.y) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'z':
              if (temp[i][j][1] < Box.Length.z) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
          }
        } else {
          switch(axis) {
            case 'x':
              if (temp[i][j][0] != (Box.Length.x/2)) {
                surf[i][j][0] += temp[i][j][0];
                values[i][j][0]++;
              }
              if (temp[i][j][1] != (Box.Length.x/2)) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'y':
              if (temp[i][j][0] != (Box.Length.y/2)) {
                surf[i][j][0] += temp[i][j][0];
                values[i][j][0]++;
              }
              if (temp[i][j][1] != (Box.Length.y/2)) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'z':
              if (temp[i][j][0] < Box.Length.z) {
                surf[i][j][0] += temp[i][j][0];
                values[i][j][0]++;
              }
              if (temp[i][j][1] > 0) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
          }
        }
      }
    } //}}}

    // free the temporary array //{{{
    for (int i = 0; i < bins[0]; i++) {
      for (int j = 0; j < bins[1]; j++) {
        free(temp[i][j]);
      }
      free(temp[i]);
    }
    free(temp); //}}}

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

  // write surface to output file //{{{
  // open output file //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  } //}}}
  PrintByline(out, argc, argv);
  // print legend
  switch(axis) {
    case 'x':
      fprintf(out, "# (1) y coordinate; (2) z coordinate;");
      break;
    case 'y':
      fprintf(out, "# (1) x coordinate; (2) z coordinate;");
      break;
    case 'z':
      fprintf(out, "# (1) x coordinate; (2) y coordinate;");
      break;
  }
  fprintf(out, " (3) surface 1; (4) surface 2\n");

  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2);
      fprintf(out, " %7.4f", width*(2*j+1)/2);
      if (values[i][j][0] > 0) {
        fprintf(out, " %7.4f", surf[i][j][0]/values[i][j][0]);
      } else {
        fprintf(out, "       ?");
      }
      if (values[i][j][1] > 0) {
        fprintf(out, " %7.4f", surf[i][j][1]/values[i][j][1]);
      } else {
        fprintf(out, "       ?");
      }
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      free(surf[i][j]);
      free(values[i][j]);
    }
    free(surf[i]);
    free(values[i]);
  }
  free(surf);
  free(values); //}}}

  return 0;
}
