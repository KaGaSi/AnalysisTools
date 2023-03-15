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

  // <output> - filename with pcf(s) //{{{
  char output_pcf[LINE] = "";
  snprintf(output_pcf, LINE, "%s", argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions_old(argc, argv, input_vsf, &verbose, &silent, LINE);
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
  int *Index_mol; // same as Index, but for molecules
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  VtfReadStruct_old(input_vsf, false, &Counts, &BeadType, &Bead, &Index,
                &MoleculeType, &Molecule, &Index_mol);
  InFile = calloc(Counts.BeadsTotal, sizeof *InFile); //}}}

  // <bead(s)> - names of bead types to use //{{{
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Use = false;
  }
  if (!all) { // --all option not used
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType_old(argv[count], Counts, BeadType);
      // error - nonexistent bead  //{{{
      if (type == -1) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", input_coor);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - non-existent bead name ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s\n", argv[count]);
        ColourReset(STDERR_FILENO);
        ErrorBeadType_old(Counts, BeadType);
        exit(1);
      } //}}}
      BeadType[type].Use = true;
    } //}}}
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      BeadType[i].Use = true;
    }
  }
  for (int i = 0; i < Counts.BeadsTotal; i++) {
    int type = Bead[i].Type;
    if (BeadType[type].Use) {
      Bead[i].Use = true;
    } else {
      Bead[i].Use = false;
    }
  }
//for (int i = 0; i < Counts.BeadsTotal; i++) {
//  printf("%s %d\n", BeadType[Bead[i].Type].Name, Bead[i].Use);
//}

  // write initial stuff to output pcf file //{{{
  FILE *out = OpenFile(output_pcf, "w");
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

  VtfReadPBC(input_coor, &Box);
  if (!TriclinicCellData(&Box)) {
//  ErrorPrintFull(vcf_file, *file_line_count, split, words);
    fprintf(stderr, "ERROR WITH INITIAL BOX STUFF\n");
    exit(1);
  }
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
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");

  // main loop //{{{
  int count_vcf = 0, // count timesteps from the beginning
      count_used = 0, // count steps used for calculation
      file_line_count = 0; // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count_vcf++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      if (count_vcf == start) {
        fprintf(stdout, "\rStarting step: %d\n", start);
      } else {
        fprintf(stdout, "\rStep: %d", count_vcf);
      }
    } //}}}
    // decide whether this timestep is to be used //{{{
    bool use = false;
    /* use if timestep
     *    1) is between start (-st option) and end (-e option)
     *    and
     *    TODO 2) isn't skipped (-sk option); skipping starts counting with 'start'
     */
    if (count_vcf >= start && (count_vcf <= end || end == -1)) {
//  if ((count_vcf >= start && (count_vcf <= end || end == -1)) && // 1)
//      ((count_vcf-start)%skip) == 0) { // 2)
      use = true;
    } else {
      use = false;
    } //}}}
    // work with the timestep, if it's to be used //{{{
    if (use) {
      if (!VtfReadTimestep_old(vcf, input_coor, &Box, &Counts, BeadType, &Bead,
                           Index, MoleculeType, Molecule,
                           &file_line_count, count_vcf)) {
        count_vcf--;
        break;
      }
      count_used++;

      for (int i = 0; i < (Counts.BeadsCoor-1); i++) {
        int id1 = InFile[i];
        if (Bead[id1].InTimestep && Bead[id1].Use) {
          for (int j = (i+1); j < Counts.BeadsCoor; j++) {
            int id2 = InFile[j];
            if (Bead[id2].InTimestep && Bead[id2].Use) {
              // bead types
              int type1 = Bead[id1].Type;
              int type2 = Bead[id2].Type;
              // type1 shouldn't be larger then type2
              // TODO why?
              if (type1 > type2) {
                SwapInt(&type1, &type2);
              }
              counter[type2]++;
              VECTOR rij = Distance(Bead[id1].Position, Bead[id2].Position,
                                    Box.Length);
              rij = FromFractional(rij, Box);
              rij.x = VectorLength(rij);
              // count only distances up to half of the shortest box length
              if (rij.x < max_dist) {
                int l = rij.x / width;
                pcf[type1][type2][l]++;
              }
            }
          }
        }
      }
    //}}}
    // skip the timestep, if it shouldn't be saved //{{{
    } else {
      if (!VtfSkipTimestep(vcf, input_coor, &file_line_count, count_vcf)) {
        count_vcf--;
        break;
      }
    } //}}}
    // decide whether to exit the main loop //{{{
    // break the loop if end timeste was reached (-e option)
    if (count_vcf == end) {
      break;
    } //}}}
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
  out = OpenFile(output_pcf, "a");
  printf("%s\n", output_pcf);
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

    // TODO: volume of triclinic? ...volume should be calculated per-step
    for (int k = 0; k < Counts.TypesOfBeads; k++) {
      for (int l = k; l < Counts.TypesOfBeads; l++) {
        if (BeadType[k].Use && BeadType[l].Use) {
          double pairs;
          if (k == l) {
            pairs = ((SQR(BeadType[k].Number) - BeadType[k].Number)) / 2;
          } else {
            pairs = BeadType[k].Number * BeadType[l].Number;
          }
          // for normalisation
          double pair_den = Box.Volume / pairs;
          double norm_factor = pair_den / shell / count_used;
          double temp = pcf[k][l][j] * norm_factor;
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
