#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
lmp_data utility generates lammps data file from the FIELD and coordinate \
files. It assumes molecules have bonds and can also have angles and \
dihedrals. For all cases, the a harmonic potential is assumed (the \
dihedrals are written to the lmp data file as impropers).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out.data> [options]\n\n", cmd);

  fprintf(ptr, "   <input>      input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <out.data>   output lammps data file\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -f[!] <name> FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      --srp        add one more bead type for srp\n");
  fprintf(ptr, "      --mass       define atom types by mass (printing their\
different charges in the Atoms section)\n");
  fprintf(ptr, "      -st <step>   timestep for creating the output file \
(default: last)\n");
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
  int req_args = 2; //}}}

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
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "--srp") != 0 &&
        strcmp(argv[i], "--mass") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-f!") != 0 &&
        strcmp(argv[i], "-st") != 0) {

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

  // <out.data> - output lammps data file
  char file_out_lmp[LINE] = "\0";
  snprintf(file_out_lmp, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent, verbose, detailed = false;
  CommonOptions(argc, argv, input_vsf, LINE, &verbose, &silent, &detailed);
  // -f[!] option //{{{
  char input_field[LINE] = "\0";
  bool change_beads = false;
  if (FileOption(argc, argv, "-f", input_field, LINE)) {
    exit(1);
  }
  if (input_field[0] == '\0') {
    if (FileOption(argc, argv, "-f!", input_field, LINE)) {
      exit(1);
    }
    if (input_field[0] != '\0') {
      change_beads = true;
    }
  } //}}}
  // add srp bead type
  bool srp = BoolOption(argc, argv, "--srp");
  // use mass only for atom type definition
  bool mass = BoolOption(argc, argv, "--mass");
  // timestep to create CONFIG file from
  int timestep = 1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  SYSTEM System = VtfReadStruct(input_vsf, detailed);

  // read coordinates //{{{
  FILE *vcf = OpenFile(input_coor, "r");
  int count_vcf = 0, // count steps in the vcf file
      file_line_count = 0; // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count_vcf++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (count_vcf == timestep) {
        fprintf(stdout, "\r                         ");
        fprintf(stdout, "\rUsing step: %d\n", timestep);
      } else {
        fprintf(stdout, "\rDiscarding step: %d", timestep);
      }
      fflush(stdout);
    } //}}}
    if (!VtfReadTimestep(vcf, input_coor, input_vsf, &System,
                         &file_line_count, count_vcf, stuff)) {
      count_vcf--;
      break;
    }
    if (count_vcf == timestep) {
      break;
    }
  }
  fclose(vcf); //}}}

  SYSTEM field = FieldReadFull(input_field);
  ChangeMolecules(&System, field, change_beads);
  CheckSystem(System, input_field);

  if (verbose) {
    VerboseOutput(System);
  }

  WriteLmpData(System, file_out_lmp, srp, mass);

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  FreeSystem(&field);
  free(stuff);
  //}}}

  return 0;
}
