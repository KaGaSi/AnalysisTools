#include "../AnalysisTools.h"
void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info gathers data about the system from provided file(s), printing it \
to standard output and creating new file(s) if prompted to. It can read \
the data from vtf file(s), lammps data file, and/or dl_meso FIELD file. \
When multiple files are provided, the following rules apply:\n\
  1) The underlying system (i.e., numbers of beads and molecules) is taken \
from the first structure file in the Info command.\n\
  2) The second file is used to supply missing data to molecules, i.e., \
connectivity, angles, etc; the extra information is added only when \
original molecule has none. The molecules from the two files must share \
only the number of beads for the extra information to be added. \
In case of '!' option, the bead types in the molecules are also \
switched for those from the second file.\n\
  3) coordinates?\n\
  4) box?\n\
\n\
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
  fprintf(ptr, "      -c_out <filename>    output CONFIG file \
for dl_software\n");
  fprintf(ptr, "      -v                   more verbose output\n");
  fprintf(ptr, "      -h                   print this help and exit\n");
  fprintf(ptr, "      --version            print version number and exit\n");
} //}}}

// TODO: implement to choose coordinates - vcf/lmp/xyz

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
    if (argv[i][0] == '-' && strcmp(argv[i], "-x_in") != 0 &&
        strcmp(argv[i], "-vs_in") != 0 && strcmp(argv[i], "-vs_in!") != 0 &&
        strcmp(argv[i], "-vc_in") != 0 && strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-f_in") != 0 && strcmp(argv[i], "-f_in!") != 0 &&
        strcmp(argv[i], "-l_in") != 0 && strcmp(argv[i], "-l_in!") != 0 &&
        strcmp(argv[i], "-vs_out") != 0 && strcmp(argv[i], "-def") != 0 &&
        strcmp(argv[i], "-f_out") != 0 && strcmp(argv[i], "-l_out") != 0 &&
        strcmp(argv[i], "--mass") != 0 && strcmp(argv[i], "-vc_out") != 0 &&
        strcmp(argv[i], "-x_out") != 0 && strcmp(argv[i], "-c_out") != 0 &&
        strcmp(argv[i], "-v") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-h") != 0 && strcmp(argv[i], "--version") != 0) {

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
  if (input_vsf[0] == '\0' && input_field[0] == '\0' &&
      input_lmp[0] == '\0' && input_xyz[0] == '\0') {
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
  SYSTEM vsf, field, lmp, xyz;
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
  } //}}}
  int primary = Min3(vs_in, f_in, l_in);
  // read vsf input if present //{{{
  if (input_vsf[0] != '\0') {
    vsf = VtfReadStruct(input_vsf, detailed);
    if (verbose) {
      printf("System in %s:\n", input_vsf);
      VerboseOutput(vsf);
    }
  } //}}}
  // read FIELD input if present //{{{
  if (input_field[0] != '\0') {
    field = FieldRead(input_field);
    if (verbose) {
      printf("System in %s:\n", input_field);
      VerboseOutput(field);
    }
  } //}}}
  // read lammps input if present //{{{
  if (input_lmp[0] != '\0') {
    lmp = LmpDataRead(input_lmp);
    if (verbose) {
      printf("System in %s:\n", input_lmp);
      VerboseOutput(lmp);
    }
  } //}}}
  // read xyz input if present //{{{
  if (input_xyz[0] != '\0') {
    xyz = XYZReadStruct(input_xyz);
  } //}}}
  // assign primary system //{{{
  char *struct_in;
  if (primary == 100) { // xyz if no 'real' structure file
    System = &xyz;
  } else if (primary == vs_in) { // vsf structure file
    System = &vsf;
    struct_in = input_vsf;
  } else if (primary == f_in) { // field structure file
    System = &field;
    struct_in = input_field;
  } else { // lammps data file
    System = &lmp;
    struct_in = input_lmp;
  } //}}}
  // vsf input (if present) //{{{
  if (primary == vs_in) {
    if (input_field[0] != '\0') {
      ChangeMolecules(System, field, change_beads_field, false);
      CheckSystem(*System, input_field);
    }
    if (input_lmp[0] != '\0') {
      ChangeMolecules(System, lmp, change_beads_lmp, false);
      CheckSystem(*System, input_lmp);
    }
  } else if (primary == l_in) {
    if (input_field[0] != '\0') {
      ChangeMolecules(System, field, change_beads_field, false);
      CheckSystem(*System, input_field);
    }
    if (input_vsf[0] != '\0') {
      ChangeMolecules(System, vsf, change_beads_vsf, false);
      CheckSystem(*System, input_vsf);
    }
  } else if (primary == f_in) {
    if (input_lmp[0] != '\0') {
      ChangeMolecules(System, lmp, change_beads_lmp, false);
      CheckSystem(*System, input_lmp);
    }
    if (input_vsf[0] != '\0') {
      ChangeMolecules(System, vsf, change_beads_vsf, false);
      CheckSystem(*System, input_vsf);
    }
  }
  if (verbose) {
    printf("System in %s:\n", struct_in);
    VerboseOutput(*System);
  } //}}}
  // first, all beads are in the timestep; revised if coordinates supplied //{{{
  for (int i = 0; i < System->Count.Bead; i++) {
    System->Bead[i].InTimestep = true;
  } //}}}
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
  // use vcf coordinate if provided //{{{
  char stuff[LINE];
  if (input_vcf[0] != '\0') {
    SYSTEM Sys_new = CopySystem(*System);
    VtfReadPBC(input_vcf, &Sys_new.Box);
    TriclinicCellData(&Sys_new.Box, 0);
    int l_count = 0;
    FILE *fr = OpenFile(input_vcf, "r");
    if (!VtfReadTimestep(fr, input_vcf, &Sys_new, &l_count, stuff)) {
      strcpy(ERROR_MSG, "not all coordinates from vcf file could be read; \
not using vcf coordinates");
      PrintWarning();
    } else {
      for (int i = 0; i < System->Count.Bead; i++) {
        System->Bead[i].InTimestep = false;
      }
      for (int i = 0; i < Sys_new.Count.BeadCoor; i++) {
        int id = Sys_new.BeadCoor[i];
        System->Bead[id].Position = Sys_new.Bead[id].Position;
        System->Bead[id].InTimestep = true;
        System->BeadCoor[i] = id;
      }
    }
    System->Box = Sys_new.Box;
    fclose(fr);
    FreeSystem(&Sys_new);
  } //}}}
  // TODO: picking Box between vtf and lmp
  // use Box from lmp if Box is unspecified in System //{{{
  if (System->Box.Volume == -1) {
    if (input_lmp[0] != '\0' && lmp.Box.Volume != -1) {
      System->Box = lmp.Box;
    }
    if (input_vsf[0] != '\0' && vsf.Box.Volume != -1) {
      System->Box = vsf.Box;
    }
  } //}}}
  // check electroneutrality //{{{
  char *second, *third;
  if (primary == vs_in) {
    second = input_lmp;
    third = input_field;
  } else if (primary == l_in) {
    second = input_vsf;
    third = input_field;
  } else {
    second = input_vsf;
    third = input_lmp;
  }
  WarnChargedSystem(*System, struct_in, second, third); //}}}
  //}}}

  // output file names & other options //{{{
  // -vs_out option //{{{
  char out_vsf[LINE] = "\0";
  if (FileOption(argc, argv, "-vs_out", out_vsf, LINE)) {
    exit(1);
  }
  if (out_vsf[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".vsf");
    if (ErrorExtension(out_vsf, ext, extension) == -1) {
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
  char out_field[LINE] = "\0";
  if (FileOption(argc, argv, "-f_out", out_field, LINE)) {
    exit(1);
  } //}}}
  // -l_out option //{{{
  char out_lmp[LINE] = "\0";
  if (FileOption(argc, argv, "-l_out", out_lmp, LINE)) {
    exit(1);
  } //}}}
  // use mass only for atom type definition (for -l_out)
  bool mass = BoolOption(argc, argv, "--mass");
  // -vc_out option //{{{
  char out_vcf[LINE] = "\0";
  if (FileOption(argc, argv, "-vc_out", out_vcf, LINE)) {
    exit(1);
  }
  if (out_vcf[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".vcf");
    if (ErrorExtension(out_vcf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  // -x_out option //{{{
  char out_xyz[LINE] = "\0";
  if (FileOption(argc, argv, "-x_out", out_xyz, LINE)) {
    exit(1);
  }
  if (out_xyz[0] != '\0') {
    ext = 1;
    strcpy(extension[0], ".xyz");
    if (ErrorExtension(out_xyz, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  // -c_out option //{{{
  char out_config[LINE] = "\0";
  if (FileOption(argc, argv, "-c_out", out_config, LINE)) {
    exit(1);
  } //}}}
  //}}}

  PruneSystem(System);

  // print information //{{{
  printf("Final system composition:\n");
  VerboseOutput(*System);
  if (verbose) { // -v option
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(*System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(*System);
  } //}}}

  // write output file(s)? //{{{
  strcpy(stuff, "# Created via Info utility from AnalysisTools");
  if (out_vsf[0] != '\0') {
    VtfWriteStruct(out_vsf, *System, default_type);
  }
  if (out_field[0] != '\0') {
    WriteField(*System, out_field);
  }
  if (out_lmp[0] != '\0') {
    WriteLmpData(*System, out_lmp, false, mass);
  }
  bool *write = malloc(sizeof *write * System->Count.Bead);
  for (int i = 0; i < System->Count.Bead; i++) {
    write[i] = true;
  }
  if (out_vcf[0] != '\0') {
    FILE *vcf = OpenFile(out_vcf, "w");
    VtfWriteCoorIndexed(vcf, stuff, write, *System);
    fclose(vcf);
  }
  if (out_xyz[0] != '\0') {
    FILE *xyz = OpenFile(out_xyz, "w");
    XyzWriteCoor(xyz, write, stuff, *System);
    fclose(xyz);
  }
  if (out_config[0] != '\0') {
    WriteConfig(*System, out_config);
  }
  free(write); //}}}

  // free memory //{{{
  if (input_vsf[0] != '\0') {
    FreeSystem(&vsf);
  }
  if (input_field[0] != '\0') {
    FreeSystem(&field);
  }
  if (input_lmp[0] != '\0') {
    FreeSystem(&lmp);
  }
  if (input_xyz[0] != '\0') {
    FreeSystem(&xyz);
  } //}}}

  return 0;
}
