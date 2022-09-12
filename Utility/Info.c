#include "../AnalysisTools.h"
void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info prints information about the provided system and can also create a new \
vtf structure file. The new structure file can differ from the original one if \
--detailed switch is used. A FIELD-like file or lammps data file \
can be used to gather additional \
information (bond types, angles, dihedrals, improper dihedral, and, \
optionally, exchange the beads from the vtf structure file); \
applies to molecules with the same name in both files.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s [options]\n\n", cmd);
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "    input files:\n");
  fprintf(ptr, "      -vs_in[!] <file.vsf> vtf structure file\n");
  fprintf(ptr, "      --detailed           differentiate bead types not just \
by names (works only with -vs_in)\n");
  fprintf(ptr, "      -vc_in <file.vcf>    vtf coordinate file\n");
  fprintf(ptr, "      -x_in <file.xyz>     xyz coordinate file\n");
  fprintf(ptr, "      -st <int>            what timestep to use; default: 1; \
for last, use '0' (works with -x_in or -vc_in)\n");
  fprintf(ptr, "      -f_in[!] <file>      input FIELD-like file\n");
  fprintf(ptr, "      -l_in[!] <file>      input lammps data file\n");
  fprintf(ptr, "    output files:\n");
  fprintf(ptr, "      -vs_out <file.vsf>   output vtf structure file\n");
  fprintf(ptr, "      -def <bead name>     default bead type \
(works with -vs_out)\n");
  fprintf(ptr, "      -f_out <file>        output FIELD-like file\n");
  fprintf(ptr, "      -l_out <file>        output lammps data file\n");
  fprintf(ptr, "      --mass               define lammps atom types by mass, \
while printing per-atom charges in Atoms section (works with -l_out)\n");
  fprintf(ptr, "      -vc_out <file.vcf>   output vtf coordinate file\n");
  fprintf(ptr, "      -x_out <file.xyz>    output xyz coordinate file\n");
  fprintf(ptr, "    general options:\n");
  fprintf(ptr, "      -v                   more verbose output\n");
  fprintf(ptr, "      -h                   print this help and exit\n");
  fprintf(ptr, "      --version            print version number and exit\n");
} //}}}

// TODO: implement to choose box - vcf/lmp
// TODO: implement to choose coordinates - vcf/lmp/xyz
// TODO: ! options not working - produce bead-less system

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
  int req_args = 0; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
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
        strcmp(argv[i], "-vs_in") != 0 &&
        strcmp(argv[i], "-vs_in!") != 0 &&
        strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-vc_in") != 0 &&
        strcmp(argv[i], "-x_in") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-f_in") != 0 &&
        strcmp(argv[i], "-f_in!") != 0 &&
        strcmp(argv[i], "-l_in") != 0 &&
        strcmp(argv[i], "-l_in!") != 0 &&
        strcmp(argv[i], "-vs_out") != 0 &&
        strcmp(argv[i], "-def") != 0 &&
        strcmp(argv[i], "-f_out") != 0 &&
        strcmp(argv[i], "-l_out") != 0 &&
        strcmp(argv[i], "--mass") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // input file names & other options //{{{
  int ext;
  char extension[2][5];
  // -vs_in[!] option //{{{
  char input_vsf[LINE] = "\0";
  bool change_beads_vsf = false;
  if (FileOption(argc, argv, "-vs_in", input_vsf, LINE)) {
    exit(1);
  }
  if (input_vsf[0] == '\0') {
    if (FileOption(argc, argv, "-vs_in!", input_vsf, LINE)) {
      exit(1);
    }
    if (input_vsf[0] != '\0') {
      change_beads_vsf = true;
    }
  }
  if (input_vsf[0] != '\0') {
    ext = 2;
    strcpy(extension[0], ".vsf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(input_vsf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  bool detailed = BoolOption(argc, argv, "--detailed");
  // -vc_in option //{{{
  char input_vcf[LINE] = "\0";
  if (FileOption(argc, argv, "-vc_in", input_vcf, LINE)) {
    exit(1);
  }
  if (input_vcf[0] != '\0') {
    ext = 2;
    strcpy(extension[0], ".vcf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(input_vcf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  // -x_in option //{{{
  char input_xyz[LINE] = "\0";
  if (FileOption(argc, argv, "-x_in", input_xyz, LINE)) {
    exit(1);
  }
  if (input_xyz[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".xyz");
    if (ErrorExtension(input_xyz, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  // timestep to use coordinates from //{{{
  int timestep = 1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  } //}}}
  // -f_in[!] option //{{{
  char input_field[LINE] = "\0";
  bool change_beads_field = false;
  if (FileOption(argc, argv, "-f_in", input_field, LINE)) {
    exit(1);
  }
  if (input_field[0] == '\0') {
    if (FileOption(argc, argv, "-f_in!", input_field, LINE)) {
      exit(1);
    }
    if (input_field[0] != '\0') {
      change_beads_field = true;
    }
  } //}}}
  // -l_in[!] option //{{{
  char input_lmp[LINE] = "\0";
  bool change_beads_lmp = false;
  if (FileOption(argc, argv, "-l_in", input_lmp, LINE)) {
    exit(1);
  }
  if (input_lmp[0] == '\0') {
    if (FileOption(argc, argv, "-l_in!", input_lmp, LINE)) {
      exit(1);
    }
    if (input_lmp[0] != '\0') {
      change_beads_lmp = true;
    }
  }
  //}}}
  // error - no input file with structure information //{{{
  if (input_vsf[0] == '\0' && input_field[0] == '\0' && input_lmp[0] == '\0') {
    strcpy(ERROR_MSG, "input vtf structure file, FIELD-like file, and/or \
lammps data file must be specified");
    PrintError();
    Help(argv[0], true);
    exit(1);
  } //}}}
  //}}}

  // options before reading system data
  bool verbose = BoolOption(argc, argv, "-v");

  // read information from input file(s) //{{{
  SYSTEM vsf, field, lmp;
  SYSTEM *System; // pointer to one of the above SYSTEMs
  // find the first structure file and read it //{{{
  int vs_in = 1e2, f_in = 1e2, l_in = 1e2;
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-vs_in") == 0 ||
        strcmp(argv[i], "-vs_in!") == 0) {
      vs_in = i;
    }
    if (strcmp(argv[i], "-f_in") == 0 ||
        strcmp(argv[i], "-f_in!") == 0) {
      f_in = i;
    }
    if (strcmp(argv[i], "-l_in") == 0 ||
        strcmp(argv[i], "-l_in!") == 0) {
      l_in = i;
    }
  }
  int primary = Min3(vs_in, f_in, l_in);
  if (primary == vs_in) {
    System = &vsf;
    *System = VtfReadStruct(input_vsf, detailed);
    if (verbose) {
      printf("System in %s:\n", input_vsf);
      VerboseOutput(*System);
    }
  } else if (primary == f_in) {
    System = &field;
    *System = FieldRead(input_field);
    if (verbose) {
      printf("System in %s:\n", input_field);
      VerboseOutput(*System);
    }
  } else {
    System = &lmp;
    *System = LmpDataRead(input_lmp);
    if (verbose) {
      printf("System in %s:\n", input_lmp);
      VerboseOutput(*System);
    }
  } //}}}
  for (int i = 0; i < System->Count.Bead; i++) {
    System->Bead[i].InTimestep = true;
  }
  // vsf input (if present) //{{{
  if (input_vsf[0] != '\0' && primary != vs_in) {
    vsf = VtfReadStruct(input_vsf, detailed);
    ChangeMolecules(System, vsf, change_beads_vsf, false);
    CheckSystem(*System, input_vsf);
    if (verbose) {
      printf("System in %s:\n", input_vsf);
      VerboseOutput(vsf);
    }
  } //}}}
  // FIELD input (if present) //{{{
  if (input_field[0] != '\0' && primary != f_in) {
    field = FieldRead(input_field);
    ChangeMolecules(System, field, change_beads_field, false);
    CheckSystem(*System, input_field);
    if (verbose) {
      printf("System in %s:\n", input_field);
      VerboseOutput(field);
    }
  } //}}}
  // lammps input (if present) //{{{
  if (input_lmp[0] != '\0' && primary != l_in) {
    lmp = LmpDataRead(input_lmp);
    ChangeMolecules(System, lmp, change_beads_lmp, false);
    CheckSystem(*System, input_lmp);
    if (verbose) {
      printf("System in %s:\n", input_lmp);
      VerboseOutput(lmp);
    }
  } //}}}
//// vcf coordinates (if present) //{{{
//if (input_vcf[0] != '\0') {
//  FILE *coor = OpenFile(input_vcf, "r");
//  char stuff[LINE];
//  int step_count = 0, file_line_count = 0;
//  VtfReadTimestep(coor, input_vcf, "\0", &vsf, &file_line_count,
//                  step_count, stuff);
//  fclose(coor);
//  PrintBox(vsf.Box);
//} else {
//  for (int i = 0; i < vsf.Count.Bead; i++) {
//    vsf.Bead[i].InTimestep = true;
//  }
//} //}}}
  // TODO: picking Box between vcf and lmp
  // use Box from lmp if Box is unspecified in System //{{{
  if (System->Box.Volume == -1) {
    if (primary != l_in && input_lmp[0] != '\0' && lmp.Box.Volume != -1) {
      vsf.Box = lmp.Box;
    }
  } //}}}
  // TODO: picking coordinates between vcf, xyz, and lmp
  // use coordinates from lmp if all coordinates in System are 0 //{{{
  bool coor = false;
  for (int i = 0; i < System->Count.Bead; i++) {
    BEAD *b_i = &System->Bead[i];
    if (b_i->Position.x != 0 || b_i->Position.y != 0 || b_i->Position.z != 0) {
      coor = true;
      break;
    }
  }
  if (!coor) {
    if (primary != l_in && input_lmp[0] != '\0') {
      for (int i = 0; i < System->Count.Bead && i < lmp.Count.Bead; i++) {
        System->Bead[i].Position = lmp.Bead[i].Position;
      }
    }
  } //}}}
  // TODO: specify correctly structur input files
  WarnChargedSystem(*System, input_vsf, input_field);
  //}}}

  // output file names & other options //{{{
  // -vs_out option //{{{
  char output_vsf[LINE] = "\0";
  if (FileOption(argc, argv, "-vs_out", output_vsf, LINE)) {
    exit(1);
  }
  if (output_vsf[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".vsf");
    if (ErrorExtension(output_vsf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  // -def option //{{{
  bool *def_type = calloc(System->Count.BeadType, sizeof *def_type);
  if (BeadTypeOption(argc, argv, "-def", false, def_type, System)) {
    exit(1);
  }
  int default_type = -1;
  for (int i = 0; i < System->Count.BeadType; i++) {
    if (def_type[i]) {
      default_type = i;
      break;
    }
  }
  free(def_type); //}}}
  // -f_out option //{{{
  char output_field[LINE] = "\0";
  if (FileOption(argc, argv, "-f_out", output_field, LINE)) {
    exit(1);
  } //}}}
  // -l_out option //{{{
  char output_lmp[LINE] = "\0";
  if (FileOption(argc, argv, "-l_out", output_lmp, LINE)) {
    exit(1);
  } //}}}
  // use mass only for atom type definition (for -l_out)
  bool mass = BoolOption(argc, argv, "--mass");
  // -vc_out option //{{{
  char output_vcf[LINE] = "\0";
  if (FileOption(argc, argv, "-vc_out", output_vcf, LINE)) {
    exit(1);
  }
  if (output_vcf[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".vcf");
    if (ErrorExtension(output_vcf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  // -x_out option //{{{
  char output_xyz[LINE] = "\0";
  if (FileOption(argc, argv, "-x_out", output_xyz, LINE)) {
    exit(1);
  }
  if (output_xyz[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".xyz");
    if (ErrorExtension(output_xyz, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  //}}}

  PruneSystem(System);

  // print information //{{{
  printf("Final system composition:\n");
  VerboseOutput(*System);
  if (verbose) { // -v option
    fprintf(stdout, "\nInformation about every bead:\n");
    PrintBead(*System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(*System);
  } //}}}

  // write output file(s)? //{{{
  if (output_vsf[0] != '\0') {
    VtfWriteStruct(output_vsf, *System, default_type);
  }
  if (output_field[0] != '\0') {
    WriteField(*System, output_field);
  }
  if (output_lmp[0] != '\0') {
    WriteLmpData(*System, output_lmp, false, mass);
  } //}}}

  // free memory //{{{
  if (input_vsf[0] != '\0') {
    FreeSystem(&vsf);
  }
  if (input_field[0] != '\0') {
    FreeSystem(&field);
  }
  if (input_lmp[0] != '\0') {
    FreeSystem(&lmp);
  } //}}}

  return 0;
}
