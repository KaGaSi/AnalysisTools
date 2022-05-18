#include "../AnalysisTools.h"
int *InFile;

// TODO triclinic box

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
lmp_vsf takes lammps data file and transforms it into vsf and vcf files. \
Different bead types are named bead1 to beadN (where N is the number of bead \
types) or as specified in comment in the Masses section of the data file. \
Molecule types are named mol1 to molN (where N is the number of \
molecule types). Molecule types are determined according to bead order as \
well as according to molecule connectivity. Angles are disregarded.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out.vsf> <out.vcf> [options]\n\n", cmd);

  fprintf(ptr, "   <input>     input lammps data file\n");
  fprintf(ptr, "   <out.vsf>   output vsf structure file\n");
  fprintf(ptr, "   <out.vcf>   output vcf coordinate file\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -v          verbose output\n");
  fprintf(ptr, "      -h          print this help and exit\n");
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
  int req_args = 3; //}}}

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
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file
  char input[LINE] = "";
  snprintf(input, LINE, "%s", argv[++count]);

  // <out.vsf> - output structure file //{{{
  char output_vsf[LINE] = "";
  snprintf(output_vsf, LINE, "%s", argv[++count]);
  // test if <out.vsf> ends with '.vsf'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output vcf file //{{{
  char output_vcf[LINE] = "";
  snprintf(output_vcf, LINE, "%s", argv[++count]);
  // test if <out.vcf> ends with '.vsf'
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  // TODO add --silent option
  bool silent = false;
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  // options before reading system data
  bool verbose = BoolOption(argc, argv, "-v");

  // some variables //{{{
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BEAD *Bead;
  BEADTYPE *BeadType;
  MOLECULE *Molecule;
  MOLECULETYPE *MoleculeType;
  BOX Box;
  VECTOR box_lo; // {x,y,z}lo from data file to place beads in (0, BoxLength>
  int bonds = 0; // total number of bonds
  int angles = 0; // total number of angles
  PARAMS *bond_type; // information about bond types
  PARAMS *angle_type; // information about angle types
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i) //}}}

  // TODO add dihedrals
  ReadLmpData(input, &bonds, &bond_type, &angles, &angle_type,
              &Box.Length, &box_lo, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
    // TODO bond & angle & dihedral types into VerboseOutput
    PrintBondTypes2(Counts.TypesOfBonds, bond_type);
    PrintAngleTypes2(Counts.TypesOfAngles, angle_type);
  } //}}}

  // write out.vcf file //{{{
  FILE *fw = OpenFile(output_vcf, "w");
  PrintByline(fw, argc, argv);
  WriteCoorIndexed(fw, Counts, BeadType, Bead,
                   MoleculeType, Molecule, "\0", Box);
  fclose(fw); //}}}

  // create & fill output vsf file
  WriteVsf_old(output_vsf, Counts, BeadType, Bead, MoleculeType, Molecule, false);

  // free memory (to make valgrind happy) //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(angle_type);
  free(bond_type); //}}}
}
