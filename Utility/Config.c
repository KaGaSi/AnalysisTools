#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Config utility generates CONFIG file from given step of a vcf file. If the \
given timestep is larger than the number of steps the coordinate file, the \
last step is used. Coordinate file needs to contain all beads in the \
simulation for it to work with dl_meso (no such check is made).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <step>     timestep for creating CONFIG (default: last)\n");
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
  int req_args = 1; //}}}

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
        strcmp(argv[i], "-st") != 0) {

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

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // timestep to create CONFIG file from //{{{
  int timestep = -1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // warn if not all beads //{{{
  if (Counts.Beads != Counts.BeadsInVsf) {
    fprintf(stderr, "\033[1;33m");
    fprintf(stdout, "\nWarning: '%s' does not contain all beads from '%s'\n\n", input_coor, input_vsf);
    fprintf(stderr, "\033[0m");
  } //}}}

  // vsf file is not needed anymore
  free(input_vsf);

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // main loop //{{{
  fpos_t pos, pos_old; // for saving pointer position in vcf file
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF && count != timestep) {
    ungetc(test, vcf);

    count++;
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    }

    // save pointer position in file
    pos_old = pos;
    fgetpos(vcf, &pos);

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
      pos = pos_old;
      count--;
    }
  }

  // restore pointer position in FIELD file
  fsetpos(vcf, &pos);

  // read coordinates //{{{
  if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
    // print newline to stdout if Step... doesn't end with one
    ErrorCoorRead(input_coor, test, count, stuff);
    exit(1);
  } //}}}

  if (!silent) {
    if (script) {
      fprintf(stdout, "CONFIG Step: %d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                         ");
      fprintf(stdout, "\rCCONFIG Step: %d\n", count);
    }
  }

  fclose(vcf); //}}}

  // create CONFIG //{{{
  // open output CONFIG file for writing //{{{
  FILE *out;
  if ((out = fopen("CONFIG", "w")) == NULL) {
    ErrorFileOpen("CONFIG", 'w');
    exit(1);
  } //}}}

  // print CONFIG file initial stuff //{{{
  fprintf(out, "NAME\n       0       1\n");
  fprintf(out, "%lf 0.000000 0.000000\n", BoxLength.x);
  fprintf(out, "0.000000 %lf 0.000000\n", BoxLength.y);
  fprintf(out, "0.000000 0.000000 %lf\n", BoxLength.z); //}}}

  // bead coordinates //{{{
  // unbonded beads must be first (dl_meso requirement)
  count = 0;
  for (int i = 0; i < Counts.Beads; i++) {
    if (Bead[i].Molecule == -1) {
      fprintf(out, "%s %d\n", BeadType[Bead[i].Type].Name, ++count);
      fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x-BoxLength.x/2,
                                    Bead[i].Position.y-BoxLength.y/2,
                                    Bead[i].Position.z-BoxLength.z/2);
    }
  }
  // bonded beads follow
  for (int i = 0; i < Counts.Beads; i++) {
    if (Bead[i].Molecule != -1) {
      fprintf(out, "%s %d\n", BeadType[Bead[i].Type].Name, ++count);
      fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x-BoxLength.x/2,
                                    Bead[i].Position.y-BoxLength.y/2,
                                    Bead[i].Position.z-BoxLength.z/2);
    }
  } //}}}

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  //}}}

  return 0;
}
