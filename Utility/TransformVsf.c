#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

// TODO: --change?
// TODO: add --detailed?

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
TransformVsf reads information from FIELD and traject.vsf files and creates \
.vsf structure file used for visualisation of trajectory (.vcf files) via VMD \
visualisation tool.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <output.vsf> [options]\n\n", cmd);
  fprintf(ptr, "   <output.vsf>   output structure file (*.vsf)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -i <name>   use other structure file than traject.vsf\n");
  fprintf(ptr, "      -v          verbose output\n");
  fprintf(ptr, "      -h          print this help and exit\n");
  fprintf(ptr, "      --change    transform molecules according to the first one of each type\n");
  fprintf(ptr, "      --version   print version number and exit\n");
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
        strcmp(argv[i], "--change") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // <output.vsf> - output structure file (must end with .vsf) //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(output, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data
  // use .vsf file other than traject.vsf?
  char input_vsf[LINE] = "";
  if (FileOption(argc, argv, "-i", input_vsf, LINE)) {
    exit(1);
  }
  if (input_vsf[0] == '\0') {
    strcpy(input_vsf, "traject.vsf");
  }

  // test if structure file ends with '.vsf' //{{{
  ext = 2;
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print new vsf so that all molecules always contain the same beads as the
  // first one of their type in the vsf
  // TODO: what exactly is that?
  bool change = BoolOption(argc, argv, "--change");
  // output verbosity
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output

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
  // TODO: --detailed option
  bool detailed = false;
  FullVtfRead(input_vsf, "\0", detailed, false, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // print information - verbose option
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // create & fill output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule, change);

  // free memory - to make valgrind happy
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);

  return 0;
}
