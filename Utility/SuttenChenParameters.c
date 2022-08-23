#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
SuttenChenParameters calculates parameters for simulations using Sutten-Chen \
potentials (https://doi.org/10.1080/09500839008206493).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> [options]\n\n", cmd);

  fprintf(ptr, "      <input>       input coordinate file \
(vcf or vtf format)\n");
  fprintf(ptr, "   [options]\n");
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
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }
  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "--reverse");
  // possible to omit <bead name(s)> if '--reverse' is used
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
        strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0) {
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

  // options before reading system data
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, input_vsf, LINE, &verbose, &silent, &detailed);

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  SYSTEM System = VtfReadStruct(input_vsf, detailed);
  VtfReadPBC(input_coor, input_vsf, &System.Box);
  if (!TriclinicCellData(&System.Box)) {
    strcpy(ERROR_MSG, "wrong pbc data");
    PrintError();
    exit(1);
  } //}}}

  WarnChargedSystem(System, input_vsf, "\0");

  // read coordinates //{{{
  FILE *vcf = OpenFile(input_coor, "r");
  int count_vcf = 0, // count steps in the vcf file
      file_line_count = 0; // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  count_vcf++;
  if (!VtfReadTimestep(vcf, input_coor, input_vsf, &System,
                       &file_line_count, count_vcf, stuff)) {
    count_vcf--;
  }
  fclose(vcf); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  COUNT *Count = &System.Count;

  // Nickel
  long double Mw = 58.693e-3; // kg/mol
  long double rho = 8938.6e3; // kg/m^3
  long double E_coh = 428.4e3; // cohesive (sublimation) energy in J/mol
  long double B_mod = 187.5e9; // bulk modulus in Pa
  int n_unit = 4; // number of atoms in repeating unit

  // common constants
  long double N_A = 6.0221367e23; // Avogadro
  long double k_B = 1.380658e-23; // Boltzman

  // system quantities
  // TODO: check ml's stuff & simplify
  long double U = Count->BeadCoor / N_A * E_coh; // energy in J
  long double vol = Count->BeadCoor * Mw / (N_A * rho); // volume in m^3
  long double lattice_const = pow(n_unit*vol/Count->BeadCoor, 1.0/3); // in m^3

  int n = 6, m; // Sutton-Chen coefficients
  long double sum_n = 0,
              sum_m = 0;
  int id1 = System.BeadCoor[0];
  VECTOR *first = &System.Bead[id1].Position;
  for (int i = 1; i < Count->BeadCoor; i++) {
    int id2 = System.BeadCoor[i];
    if (id1 != id2) {
      VECTOR dist = Distance(*first, System.Bead[i].Position, System.Box.Length);
      dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
      sum_n += 1 / pow(dist.x, n);
      sum_m += 1 / pow(dist.x, m);
    }
  }
  sum_n *= pow(lattice_const, n);
  sum_m *= pow(lattice_const, m);

  // free memory
  FreeSystem(&System);
  free(stuff);

  return 0;
}
