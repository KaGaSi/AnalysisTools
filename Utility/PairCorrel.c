#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, " \
PairCorrel utility calculates pair correlation function for specified \
bead types. All pairs of bead types (including same type pairs) are \
calculated - given A and B types, pcf between A-A, A-B and B-B are \
calculated.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <bead type(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>          input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <width>          width of a single bin\n");
  fprintf(ptr, "   <output>         output file with pair correlation function(s)\n");
  fprintf(ptr, "   <bead type(s)>   bead type name(s) for pcf calculation\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -n <int>      number of bins to average\n");
  fprintf(ptr, "      -st <int>     starting timestep for calculation (default: 1)\n");
  fprintf(ptr, "      -e <end>      ending timestep for calculation (default: none)\n");
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
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
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
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
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
  char input_coor[LINE];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  char *input_vsf = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (!IsPosDouble(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double bin_width = atof(argv[count]); //}}}

  // <output> - filename with pcf(s) //{{{
  char output_pcf[LINE];
  strcpy(output_pcf, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // ending timestep //{{{
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  } //}}}

  // error if ending step is lower than starging step //{{{
  if (end != -1 && start > end) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - non-existent bead name \033[1;33m%s\033[1;31m\n", input_coor, argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // write initial stuff to output pcf file //{{{
  FILE *out;
  if ((out = fopen(output_pcf, "w")) == NULL) {
    ErrorFileOpen(output_pcf, 'w');
    exit(1);
  }

  // print command to output file
  putc('#', out);
  PrintCommand(out, argc, argv);

  // print bead type names to output file //{{{
  fprintf(out, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = i; j < Counts.TypesOfBeads; j++) {
      if (BeadType[i].Use && BeadType[j].Use) {
        fprintf(out, " (%d) %s-%s", ++count, BeadType[i].Name, BeadType[j].Name);
        if (i != (Counts.TypesOfBeads-1) || j != (Counts.TypesOfBeads-1)) {
          putc(';', out);
        }
      }
    }
  }
  putc('\n', out); //}}}

  fclose(out); //}}}

  // number of bins - maximum distance is taken as half of the shortes BoxLength //{{{
  double max_dist = Min3(BoxLength.x, BoxLength.y, BoxLength.z) / 2;
  int bins = ceil(max_dist / bin_width); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // allocate memory for pcf arrays //{{{
  int *counter = calloc(Counts.TypesOfBeads,sizeof(int *)); // to count number of pairs
  double ***pcf = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    pcf[i] = malloc(Counts.TypesOfBeads*sizeof(double *));
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      pcf[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
      exit(1);
    }
  }
  if (test == EOF) {
    fflush(stdout);
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file - %d steps to discard, but %d steps in all\n\n", input_coor, start, count);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  // print number of starting step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rStarting step: %d\n", start);
    }
  } //}}}
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  }
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int count_vcf = start - 1;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      ErrorCoorRead(input_coor, test, count_vcf, stuff);
      exit(1);
    } //}}}

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
            VECTOR rij = Distance(Bead[bead1].Position, Bead[bead2].Position, BoxLength);
            rij.x = Length(rij);

            // count only distances up to half of the shortest box length
            if (rij.x < max_dist) {
              int l = rij.x / bin_width;
              pcf[type1][type2][l]++;
            }
          }
        }
      }
    } //}}}

    if (end == count_vcf)
      break;
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  } //}}}

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
  for (int j = 0; j < (bins-avg); j++) {

    // calculate volume of every shell that will be averaged
    double shell[avg];
    for (int k = 0; k < avg; k++) {
      shell[k] = 4.0 / 3 * PI * CUBE(bin_width) * (CUBE(j+k+1) - CUBE(j+k));
    }

    fprintf(out, "%8.5f", bin_width*(j+0.5*avg));

    double volume = BoxLength.x * BoxLength.y * BoxLength.z;
    for (int k = 0; k < Counts.TypesOfBeads; k++) {
      for (int l = k; l < Counts.TypesOfBeads; l++) {
        if (BeadType[k].Use && BeadType[l].Use) {

          double temp = 0; // for normalisation

          // sum up pcfs from all shells to be averaged
          for (int m = 0; m < avg; m++) {
            double pairs;
            if (k == l) {
              pairs = ((SQR(BeadType[k].Number) - BeadType[k].Number)) / 2;
            } else {
              pairs = BeadType[k].Number * BeadType[l].Number;
            }
            // for normalisation
            double pair_den = volume / pairs;
            double norm_factor = pair_den / shell[m] / count;
            temp += pcf[k][l][j+m] * norm_factor;
          }

          // print average value to output file
          fprintf(out, " %10f", temp/avg);
        }
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
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
