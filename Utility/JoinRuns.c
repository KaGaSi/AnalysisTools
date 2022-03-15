#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
JoinRuns joins two simulation runs with different .vsf files. The first .vsf is \
assumed to be traject.vsf (if not, use '-i' option) and the FIELD file has to \
be the same for both simulation runs. Bead types in both .vcf files must be the \
same, but only selected bead types are saved to output.vcf file.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <1st input.vcf> <2nd input.vsf> <2nd input.vcf> ", cmd);
  fprintf(ptr, "<output.vcf> <type names> <options>\n\n");

  fprintf(ptr, "   <1st input.vcf>   input coordinate file of 1st run (vcf or vtf format)\n");
  fprintf(ptr, "   <2nd input.vsf>   input coordinate file of 2nd run (vsf or vtf format)\n");
  fprintf(ptr, "   <2nd input.vcf>   input coordinate file of 2nd run (vcf or vtf format)\n");
  fprintf(ptr, "   <output.vcf>      output filename (vcf format)\n");
  fprintf(ptr, "   <type names>      names of bead types to save (optional if '-r' used)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -r             reverse <type name(s)>, i.e., exclude the specified bead types\n");
  fprintf(ptr, "      --join         join molecules (remove pbc)\n");
  fprintf(ptr, "      -st1 <int>     starting timestep from 1st run\n");
  fprintf(ptr, "      -st2 <int>     starting timestep from 2nd run\n");
  fprintf(ptr, "      -e1 <end>      ending timestep from the 1st run\n");
  fprintf(ptr, "      -e2 <end>      ending timestep from the 2st runh\n");
  fprintf(ptr, "      -sk1 <int>     leave out every <int> steps from 1st run\n");
  fprintf(ptr, "      -sk2 <int>     leave out every <int> steps from 2st run\n");
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
  int req_args = 5; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "-r");

  // possible to omit <type name(s)> if '-r' is used
  if (count < (req_args-1) || (count == (req_args-1) && !reverse)) {
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
        strcmp(argv[i], "-r") != 0 &&
        strcmp(argv[i], "--join") != 0 &&
        strcmp(argv[i], "-st1") != 0 &&
        strcmp(argv[i], "-st2") != 0 &&
        strcmp(argv[i], "-e1") != 0 &&
        strcmp(argv[i], "-e2") != 0 &&
        strcmp(argv[i], "-sk1") != 0 &&
        strcmp(argv[i], "-sk2") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <1st input> - first input coordinate file //{{{
  char input_coor_1[LINE];
  strcpy(input_coor_1, argv[++count]);

  // test if <1st input> ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor_1, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  char *input_vsf_1 = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor_1, '.'),".vtf") == 0) {
    strcpy(input_vsf_1, input_coor_1);
  } else {
    strcpy(input_vsf_1, "traject.vsf");
  } //}}}

  // <2nd input.vsf> - second structure file //{{{
  char *input_vsf_2 = calloc(LINE,sizeof(char));
  strcpy(input_vsf_2, argv[++count]);

  // test if <2nd input.vsf> ends with '.vsf' (required by VMD)
  ext = 2;
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_vsf_2, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <2nd input> - second input coordinate file //{{{
  char input_coor_2[LINE];
  strcpy(input_coor_2, argv[++count]);

  // test if <2nd input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor_2, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[LINE];
  strcpy(output_vcf, argv[++count]);

  // test if output coordinate file ends with '.vcf' (required by vmd)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, &input_vsf_1, &verbose, &silent);

  // save coordinates of joined aggregates //{{{
  char joined_vcf[LINE];
  bool error = JoinCoorOption(argc, argv, joined_vcf);
  if (error) {
    exit(1);
  } //}}}

  // should output coordinates be joined? //{{{
  bool join = BoolOption(argc, argv, "--join"); //}}}

  // starting timesteps for both simulations //{{{
  int start_1 = 1, start_2 = 1;
  if (IntegerOption(argc, argv, "-st1", &start_1)) {
    exit(1);
  }
  if (IntegerOption(argc, argv, "-st2", &start_2)) {
    exit(1);
  } //}}}

  // ending timestep //{{{
  int end_1 = -1, end_2 = -1;
  if (IntegerOption(argc, argv, "-e1", &end_1)) {
    exit(1);
  }
  if (IntegerOption(argc, argv, "-e2", &end_2)) {
    exit(1);
  } //}}}

  // error if ending step is lower than starging step //{{{
  // TODO colours
  if (end_1 != -1 && start_1 > end_1) {
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "For the first run, starting step (%d) is higher than ending step (%d)\n", start_1, end_1);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  // TODO colours
  if (end_2 != -1 && start_2 > end_2) {
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "For the first run, starting step (%d) is higher than ending step (%d)\n", start_2, end_2);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // skipped timesteps per used step //{{{
  int skip_1 = 0, skip_2 = 0;
  if (IntegerOption(argc, argv, "-sk1", &skip_1)) {
    exit(1);
  }
  if (IntegerOption(argc, argv, "-sk2", &skip_2)) {
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}
  //}}}

  // variables - structures //{{{
  // data from 1st run
  BEADTYPE *BeadType1; // structure with info about all bead types
  MOLECULETYPE *MoleculeType1; // structure with info about all molecule types
  BEAD *Bead1; // structure with info about every bead
  int *Index1; // link between indices (i.e., Index[Bead[i].Index]=i)
  MOLECULE *Molecule1; // structure with info about every molecule
  // data from 2nd run
  BEADTYPE *BeadType2; // structure with info about all bead types
  MOLECULETYPE *MoleculeType2; // structure with info about all molecule types
  BEAD *Bead2; // structure with info about every bead
  int *Index2; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule2; // structure with info about every molecule
  // Counts is the same for both runs
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // read system information //{{{
  bool indexed = ReadStructure(input_vsf_1, input_coor_1, &Counts, &BeadType1, &Bead1, &Index1, &MoleculeType1, &Molecule1);
  ReadStructure(input_vsf_2, input_coor_2, &Counts, &BeadType2, &Bead2, &Index2, &MoleculeType2, &Molecule2);

  // vsf files are not needed anymore
  free(input_vsf_1);
  free(input_vsf_2);

  // set all molecule types to write to output.vcf
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType1[i].Write = true;
    MoleculeType2[i].Write = true;
  } //}}}

  // <type names> - names of bead types to save //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType1);
    // TODO colours
    if (type == -1) {
      ErrorPrintError_old();
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "non-existent bead name \033[1;33m%s\033[1;31m\n", argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorBeadType(Counts, BeadType1);
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "Note that all bead types must be in both coordinate files.\n\n");
      fprintf(stderr, "\033[0m");
      exit(1);
    }
    BeadType1[type].Write = true;
  } //}}}

  // if '-r' is used, switch Write bools for all bead types //{{{
  if (reverse) {
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      if (BeadType1[i].Write) {
        BeadType1[i].Write = false;
      } else {
        BeadType1[i].Write = true;
      }
    }
  } //}}}

  // print selected bead type names to output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  }

  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType1[i].Write) {
      fprintf(out, "# %s\n", BeadType1[i].Name);
    }
  }

  fclose(out); //}}}

  // open input coordinate files //{{{
  FILE *vcf_1, *vcf_2;
  if ((vcf_1 = fopen(input_coor_1, "r")) == NULL) {
    ErrorFileOpen(input_coor_1, 'r');
    exit(1);
  }
  if ((vcf_2 = fopen(input_coor_2, "r")) == NULL) {
    ErrorFileOpen(input_coor_2, 'r');
    exit(1);
  } //}}}

  // get to the pbc line in both input files //{{{
  VECTOR BoxLength = GetPBC(vcf_1, input_coor_1);
  VECTOR BoxLength_2 = GetPBC(vcf_2, input_coor_2);
  // check that the box sizes are the same
  // TODO colours
  if (BoxLength.x != BoxLength_2.x ||
      BoxLength.y != BoxLength_2.y ||
      BoxLength.z != BoxLength_2.z) {
    ErrorPrintError_old();
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "- different box sizes in provided coordinate files\n");
    fprintf(stderr, "        %s: %lf %lf %lf\n", input_coor_1, BoxLength.x, BoxLength.y, BoxLength.z);
    fprintf(stderr, "        %s: %lf %lf %lf\n\n", input_coor_1, BoxLength_2.x, BoxLength_2.y, BoxLength_2.z);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
    VerboseOutput(input_coor_1, Counts, BoxLength, BeadType1, Bead1, MoleculeType1, Molecule1);
    fprintf(stdout, "\n   Starting from %d. (%d.) timestep\n", start_1, start_2);
    fprintf(stdout, "   Every %d. (%d.) timestep used\n", skip_1+1, skip_2+1);
  } //}}}

  // print pbc to output .vcf file //{{{
  if ((out = fopen(output_vcf, "a")) == NULL) {
    ErrorFileOpen(output_vcf, 'a');
    exit(1);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  fclose(out); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // start first run with start-th step //{{{
  int test;
  // first run
  count = 0;
  for (int i = 1; i < start_1; i++) {
    count++;

    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded from 1st run: %d", count);
    }

    SkipVcfCoor(vcf_1, input_coor_1, Counts, &stuff);
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step for 1st run: %d\n", start_1);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                                        ");
      fprintf(stdout, "\rStarting step for 1st run: %d\n", start_1);
    }
  }
  // is the vcf file continuing?
  if (ErrorDiscard(start_1, count, input_coor_1, vcf_1)) {
    exit(1);
  }
  //}}}

  // main loop - 1st run //{{{
  int count_vcf = start_1 - 1;
  while ((test = getc(vcf_1)) != EOF) {
    ungetc(test, vcf_1);

    ReadCoordinates(indexed, input_coor_1, vcf_1, Counts, Index1, &Bead1, &stuff);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep from 1st run: %d", count_vcf);
    } //}}}

    // join molecules? //{{{
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType1, &Bead1, MoleculeType1, Molecule1);
    } //}}}

    // open output .vcf file for appending //{{{
    if ((out = fopen(output_vcf, "a")) == NULL) {
      ErrorFileOpen(output_vcf, 'a');
      exit(1);
    } //}}}

    WriteCoorIndexed(out, Counts, BeadType1, Bead1, MoleculeType1, Molecule1, stuff);

    fclose(out);

    // skip every 'skip' steps //{{{
    for (int i = 0; i < skip_1; i++) {
      // test whether at vcf's eof //{{{
      if ((test = getc(vcf_1)) == EOF) {
        break;
      }
      ungetc(test, vcf_1); //}}}

      count++;
      count_vcf++;

      ReadCoordinates(indexed, input_coor_1, vcf_1, Counts, Index1, &Bead1, &stuff);
    } //}}}

    if (end_1 == count_vcf)
      break;
  }

  // print last step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Last step from 1st run: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                                       ");
      fprintf(stdout, "\rLast step from 1st run: %d\n", count_vcf);
    }
  } //}}}

  fclose(vcf_1); //}}}

  // connect Bead2 with Bead1 via Index arrays //{{{
  if (!silent && !script) {
    fprintf(stderr, "\033[1;33m");
    fprintf(stdout, "\nWarning: depending on the number of beads in the system, connecting bead indices from the two runs may take a long time\n\n");
    fprintf(stderr, "\033[0m");
  }
  bool used[Counts.BeadsCoor];
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    used[i] = false;
  }

  for (int i = 0; i < Counts.BeadsCoor; i++) {
    if (Bead1[i].Molecule == -1) { // monomer bead
      for (int j = 0; j < Counts.BeadsCoor; j++) {
        if (Bead2[i].Molecule == -1 &&
            strcmp(BeadType2[Bead2[j].Type].Name, BeadType1[Bead1[i].Type].Name) == 0 &&
            !used[j]) {
          Index1[Bead1[i].Index] = Index2[Bead2[j].Index];

          used[j] = true;

          break;
        }
      }
    } else { // molecular bead
      for (int j = 0; j < Counts.BeadsCoor; j++) {
        if (Bead2[j].Molecule != -1) {
          // mol from run 1
          int mol_id_1 = Bead1[i].Molecule;
          int mol_type_1 = Molecule1[mol_id_1].Type;
          // mol from run 2
          int mol_id_2 = Bead2[j].Molecule;
          int mol_type_2 = Molecule2[mol_id_2].Type;

          if (!used[j] && strcmp(MoleculeType1[mol_type_1].Name, MoleculeType2[mol_type_2].Name) == 0) {
            test = -1;
            // find where exactly in the molecule the bead is
            for (int k = 0; k <= MoleculeType2[mol_type_2].nBeads; k++) {
              if (j == Molecule2[mol_id_2].Bead[k]) {
                test = k;
                break;
              }
            }
            test = Bead1[Molecule1[mol_id_1].Bead[test]].Index;
            Index1[test] = Index2[Bead2[j].Index];
            used[j] = true;
            break;
          }
        }
      }
    }
    if (!silent && !script) {
      fflush(stdout);
      printf("\r%d. bead done", i);
    }
  } //}}}

  // start second run with start-th step //{{{
  count = 0;
  for (int i = 1; i < start_2; i++) {
    count++;

    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded from 2nd run: %d", count);
    }

    SkipVcfCoor(vcf_2, input_coor_2, Counts, &stuff);
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step for 2nd run: %d\n", start_2);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                                        ");
      fprintf(stdout, "\rStarting step for 2nd run: %d\n", start_2);
    }
  }
  // is the vcf file continuing?
  if (ErrorDiscard(start_2, count, input_coor_2, vcf_2)) {
    exit(1);
  }
  //}}}

  // main loop - 2nd run //{{{
  count_vcf = start_2 - 1;
  while ((test = getc(vcf_2)) != EOF) {
    ungetc(test, vcf_2);

    ReadCoordinates(indexed, input_coor_2, vcf_2, Counts, Index2, &Bead2, &stuff);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep from 2nd run: %d", count_vcf);
    } //}}}

    // join molecules? //{{{
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType2, &Bead2, MoleculeType2, Molecule2);
    } //}}}

    // open output .vcf file for appending //{{{
    if ((out = fopen(output_vcf, "a")) == NULL) {
      ErrorFileOpen(output_vcf, 'a');
      exit(1);
    } //}}}

    for (int i = 0; i < Counts.BeadsCoor; i++) {
        Bead1[i].Position.x = Bead2[Index1[Bead1[i].Index]].Position.x;
        Bead1[i].Position.y = Bead2[Index1[Bead1[i].Index]].Position.y;
        Bead1[i].Position.z = Bead2[Index1[Bead1[i].Index]].Position.z;
    }

    // Molecule(Type)1 and Bead(Type)1 are used, because *2 were copied to them
    WriteCoorIndexed(out, Counts, BeadType1, Bead1, MoleculeType1, Molecule1, stuff);

    fclose(out);

    // skip every 'skip' steps //{{{
    for (int i = 0; i < skip_2; i++) {
      // test whether at vcf's eof //{{{
      if ((test = getc(vcf_2)) == EOF) {
        break;
      }
      ungetc(test, vcf_2); //}}}

      count++;
      count_vcf++;

      ReadCoordinates(indexed, input_coor_2, vcf_2, Counts, Index2, &Bead2, &stuff);
    } //}}}

    if (end_2 == count_vcf)
      break;
  }
  fclose(vcf_2);

  // print last step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Last step from 2nd run: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                                       ");
      fprintf(stdout, "\rLast step from 2nd run: %d\n", count_vcf);
    }
  } //}}}
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType1, &Molecule1, &BeadType1, &Bead1, &Index1);
  FreeSystemInfo(Counts, &MoleculeType2, &Molecule2, &BeadType2, &Bead2, &Index2);
  free(stuff);
  //}}}

  return 0;
}
