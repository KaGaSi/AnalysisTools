#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info simply prints information about a system in the provided structure file. \
The verbose option prints detailed information about every molecule as \
well.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> [options]\n\n", cmd);
  fprintf(ptr, "   <input>   vtf input structure file\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --detailed   differentiate bead types not just \
by names\n");
  fprintf(ptr, "      -v           verbose output\n");
  fprintf(ptr, "      -h           print this help and exit\n");
  fprintf(ptr, "      --version    print version number and exit\n");
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
        strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // <input> - input structure file (must end with .vsf or .vtf) //{{{
  char input_vsf[LINE] = "";
  snprintf(input_vsf, LINE, "%s", argv[++count]);
  // test if <input> filename ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data //{{{
  bool verbose = BoolOption(argc, argv, "-v");
  bool detailed = BoolOption(argc, argv, "--detailed");
  //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  int *Index_mol; // same as Index, but for molecules
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VtfReadStruct(input_vsf, detailed, &Counts, &BeadType, &Bead, &Index,
                &MoleculeType, &Molecule, &Index_mol);
  InFile = calloc(Counts.BeadsTotal, sizeof *InFile); //}}}

  // print information
  VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  // TODO: if beads in vsf != beads in vcf, write out vcf
  if (verbose) { //{{{
    fprintf(stdout, "\nInformation about every bead:\n");
    PrintBead2(Counts.BeadsCoor, Index, BeadType, Bead);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(Counts.Molecules, MoleculeType, Molecule, BeadType, Bead);
  } //}}}

  // free memory - to make valgrind happy
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);

  return 0;
}
