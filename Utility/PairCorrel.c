#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
PairCorrel utility calculates pair correlation function for specified \
bead types. All pairs of bead types (including same type pairs) are \
calculated - given A and B types, pcf between A-A, A-B and B-B are \
calculated.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <bead(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>     input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>     width of a single bin\n");
  fprintf(ptr, "   <output>    output file with pair correlation \
function(s)\n");
  fprintf(ptr, "   <bead(s)>   bead name(s) for calculation \
(optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --all       use all bead types (overwrites <bead(s)>)\n");
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
  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }

  // use all bead types? ...do now to check correct number of arguments
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

  // <output> - filename with pcf(s) //{{{
  char output_pcf[LINE] = "";
  snprintf(output_pcf, LINE, "%s", argv[++count]); //}}}

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
              &MoleculeType, &Molecule); //}}}

  // <bead(s)> - names of bead types to use //{{{
  if (!all) { // --all option not used
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], Counts, BeadType);
      // error - nonexistent bead  //{{{
      if (type == -1) {
        ErrorPrintError_old();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - non-existent bead name ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n", argv[count]);
        ResetColour(STDERR_FILENO);
        ErrorBeadType(Counts, BeadType);
        exit(1);
      } //}}}
      BeadType[type].Use = true;
    } //}}}
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      BeadType[i].Use = true;
    }
  }

  // write initial stuff to output pcf file //{{{
  FILE *out;
  if ((out = fopen(output_pcf, "w")) == NULL) {
    ErrorFileOpen(output_pcf, 'w');
    exit(1);
  }
  PrintByline(out, argc, argv);
  // print bead type names to output file //{{{
  fprintf(out, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = i; j < Counts.TypesOfBeads; j++) {
      if (BeadType[i].Use && BeadType[j].Use) {
        count++;
        fprintf(out, " (%d) %s-%s", count, BeadType[i].Name, BeadType[j].Name);
        if (i != (Counts.TypesOfBeads-1) || j != (Counts.TypesOfBeads-1)) {
          putc(';', out);
        }
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  // number of bins //{{{
  double max_dist = 0.5 * Min3(Box.Length.x, Box.Length.y, Box.Length.z);
  int bins = ceil(max_dist / width); //}}}

// TODO: sizeof ...argh!
  // allocate memory //{{{
  // array counting number of pairs
  int *counter = calloc(Counts.TypesOfBeads, sizeof *counter);
  // pair correlation function
  // TODO: possibly change to long int?
//double ***pcf = malloc(Counts.TypesOfBeads * sizeof(double **));
  double ***pcf = malloc(Counts.TypesOfBeads * sizeof **pcf);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  pcf[i] = malloc(Counts.TypesOfBeads * sizeof(double *));
    pcf[i] = malloc(Counts.TypesOfBeads * sizeof *pcf[i]);
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
//    pcf[i][j] = calloc(bins, sizeof(double));
      pcf[i][j] = calloc(bins, sizeof *pcf[i][j]);
    }
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
    // read coordinates & wrap box
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    RestorePBC(Counts.BeadsCoor, Box, &Bead);
    // TODO: check + fractionals?
    // calculate pair correlation function //{{{
    for (int j = 0; j < (Counts.Bonded+Counts.Unbonded); j++) {
      if (BeadType[Bead[j].Type].Use) {

        for (int k = (j+1); k < (Counts.Bonded+Counts.Unbonded); k++) {
          if (BeadType[Bead[k].Type].Use) {

            int bead1 = j;
            int bead2 = k;

            int type1 = Bead[bead1].Type;
            int type2 = Bead[bead2].Type;

            // type1 shouldn't be larger then type2 //{{{
            if (type1 > type2) {
              int temp = type1;
              type1 = type2;
              type2 = temp;

              temp = bead1;
              bead1 = bead2;
              bead2 = temp;
            } //}}}

            counter[type2]++;

            // distance between bead1 and bead2
            // TODO: fractional coordinates?
            VECTOR rij = Distance(Bead[bead1].Position, Bead[bead2].Position, Box.Length);
            rij.x = Length(rij);

            // count only distances up to half of the shortest box length
            if (rij.x < max_dist) {
              int l = rij.x / width;
              pcf[type1][type2][l]++;
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

// TODO: check
  // write data to output file(s) //{{{
  if ((out = fopen(output_pcf, "a")) == NULL) {
    ErrorFileOpen(output_pcf, 'a');
    exit(1);
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    counter[0] = 0;
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Use) {
      counter[0] += BeadType[i].Number;
    }
  }

  // calculate pcf
  for (int j = 0; j < bins; j++) {

    // calculate volume of every shell that will be averaged
    double shell;
    shell = 4.0 / 3 * PI * CUBE(width) * (CUBE(j+1) - CUBE(j));
    fprintf(out, "%8.5f", width*(2*j+1)/2);

    // TODO: volume of triclinic?
    double volume = Box.Length.x * Box.Length.y * Box.Length.z;
    for (int k = 0; k < Counts.TypesOfBeads; k++) {
      for (int l = k; l < Counts.TypesOfBeads; l++) {
        if (BeadType[k].Use && BeadType[l].Use) {

          // sum up pcfs from all shells to be averaged
//        for (int m = 0; m < 1; m++) {
            double pairs;
            if (k == l) {
              pairs = ((SQR(BeadType[k].Number) - BeadType[k].Number)) / 2;
            } else {
              pairs = BeadType[k].Number * BeadType[l].Number;
            }
            // for normalisation
            double pair_den = volume / pairs;
            double norm_factor = pair_den / shell / count_step;
            double temp = pcf[k][l][j] * norm_factor;
//        }

          // print average value to output file
          fprintf(out, " %10f", temp);
        }
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      free(pcf[i][j]);
    }
    free(pcf[i]);
  }
  free(pcf);
  free(counter); //}}}

  return 0;
}
