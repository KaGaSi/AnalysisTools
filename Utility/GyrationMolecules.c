#include "../AnalysisTools.h"
char ERROR_MSG[LINE];

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
GyrationMolecules calculates the gyration tensor for molecules and \
determines their shape descriptors like the radius of gyration, \
acylindricity, asphericity, or relative shape anisotropy. It writes \
per-timestep averages to the output file and appends overall averages to \
that file.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <output> <mol(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <output>   output file with shape descriptors \
(automatic ending '<molecule_name>.rho' added)\n");
  fprintf(ptr, "   <mol(s)>   molecule types to calculate shape \
descriptors for (optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --all       all molecules (overwrites <mol(s)>)\n");
  fprintf(ptr, "      --joined    specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -bt         specify bead types to be used for \
calculation (default is all)\n");
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
  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }

  // use all molecules? ...do now to check correct number of arguments
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
        strcmp(argv[i], "--joined") &&
        strcmp(argv[i], "-bt") != 0 &&
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
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - filename with shape descriptors //{{{
  char output[LINE-MOL_NAME-4] = ""; // enough space for attaching <name>.txt
  snprintf(output, LINE-MOL_NAME-4, "%s", argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  int start, end;
  StartEndTime(argc, argv, &start, &end);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

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
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // <mol(s)> - names of molecule types to use //{{{
  if (!all) { // --all option not used
    while (++count < argc && argv[count][0] != '-') {
      int mtype = FindMoleculeType(argv[count], Counts, MoleculeType);
      // error - nonexistent molecule  //{{{
      if (mtype == -1) {
        ErrorPrintError_old();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - non-existent molecule type ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", argv[count]);
        RedText(STDERR_FILENO);
        fprintf(stderr, "\n");
        ResetColour(STDERR_FILENO);
        ErrorMoleculeType(Counts, MoleculeType);
        exit(1);
      } //}}}
      MoleculeType[mtype].Use = true;
    }
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      MoleculeType[i].Use = true;
    }
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", true, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // write initial stuff to output file //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      FILE *out;
      char str[LINE];
      sprintf(str, "%s%s.txt", output, MoleculeType[i].Name);
      if ((out = fopen(str, "w")) == NULL) {
        ErrorFileOpen(str, 'w');
        exit(1);
      }
      PrintByline(out, argc, argv);
      fprintf(out, "# (1) dt, (2) <Rg>, (3) <Rg^2>, ");
      fprintf(out, "(4) <Anis>, (5) <Acyl>, (6) <Aspher>");
      fprintf(out, "(7) <eigen.x>, (8) <eigen.y>, (9) <eigen.z>");
      putc('\n', out);

      fclose(out);
    }
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // allocate memory for sums of shape descriptors //{{{
  double *Rg_sum = calloc(Counts.TypesOfMolecules, sizeof *Rg_sum);
  double *sqrRg_sum = calloc(Counts.TypesOfMolecules, sizeof *sqrRg_sum);
  double *Anis_sum = calloc(Counts.TypesOfMolecules, sizeof *Anis_sum);
  double *Acyl_sum = calloc(Counts.TypesOfMolecules, sizeof *Acyl_sum);
  double *Aspher_sum = calloc(Counts.TypesOfMolecules, sizeof *Aspher_sum);
  VECTOR *eigen_sum = calloc(Counts.TypesOfMolecules, sizeof *eigen_sum); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

  // main loop //{{{
  int count_vcf = 0; // count all timesteps
  count = 0; // timesteps from -st to -e
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}
    // read & join molecules //{{{
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // join molecules if un-joined coordinates provided
    if (!joined) {
      // transform coordinates into fractional ones for non-orthogonal box
      ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      RemovePBCMolecules(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}


    // allocate arrays for the timestep //{{{
    double *Rg_step = calloc(Counts.TypesOfMolecules, sizeof *Rg_step);
    double *sqrRg_step = calloc(Counts.TypesOfMolecules, sizeof *sqrRg_step);
    double *Anis_step = calloc(Counts.TypesOfMolecules, sizeof *Anis_step);
    double *Acyl_step = calloc(Counts.TypesOfMolecules, sizeof *Acyl_step);
    double *Aspher_step = calloc(Counts.TypesOfMolecules, sizeof *Aspher_step);
    VECTOR *eigen_step = calloc(Counts.TypesOfMolecules, sizeof *eigen_step);
    //}}}

    // calculate shape descriptors //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;

      if (MoleculeType[mtype].Use) {
        // copy bead ids to a separate array //{{{
        int list[MoleculeType[mtype].nBeads],
            n = 0;
        for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
          int bead_id = Molecule[i].Bead[j];
          if (BeadType[Bead[bead_id].Type].Use) {
            list[n] = Molecule[i].Bead[j];
            n++;
          }
        } //}}}
        // TODO: fractionals working?
        VECTOR eigen = Gyration(n, list, Counts, BeadType, &Bead);
        eigen = FromFractional(eigen, Box);
        // warning - negative eigenvalues //{{{
        if (eigen.x < 0 || eigen.y < 0 || eigen.z < 0) {
          YellowText(STDERR_FILENO);
          fprintf(stderr, "\nWarning: negative eigenvalues (");
          CyanText(STDERR_FILENO);
          fprintf(stderr, "%lf", eigen.x);
          YellowText(STDERR_FILENO);
          fprintf(stderr, ", ");
          CyanText(STDERR_FILENO);
          fprintf(stderr, "%lf", eigen.y);
          YellowText(STDERR_FILENO);
          fprintf(stderr, ", ");
          CyanText(STDERR_FILENO);
          fprintf(stderr, "%lf", eigen.z);
          YellowText(STDERR_FILENO);
          fprintf(stderr, ")\n");
          ResetColour(STDERR_FILENO);
        } //}}}
        double Rgi = sqrt(eigen.x + eigen.y + eigen.z);
        // radius of gyration
        Rg_step[mtype] += Rgi; // for number avg
        // squared radius of gyration
        sqrRg_step[mtype] += SQR(Rgi); // for number avg
        // relative shape anisotropy
        Anis_step[mtype] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) /
                            SQR(eigen.x + eigen.y + eigen.z) - 0.5;
        // acylindricity
        Acyl_step[mtype] += eigen.y - eigen.x;
        // asphericity
        Aspher_step[mtype] += eigen.z - 0.5 * (eigen.x + eigen.y);
        // eigenvalues
        eigen_step[mtype].x += eigen.x;
        eigen_step[mtype].y += eigen.y;
        eigen_step[mtype].z += eigen.z;
      }
    } //}}}

    // add values to sums //{{{
    if (count_vcf >= start && (end == -1 || count_vcf <= end)) {
      count++;
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        Rg_sum[i] += Rg_step[i];
        sqrRg_sum[i] += sqrRg_step[i];
        Anis_sum[i] += Anis_step[i];
        Acyl_sum[i] += Acyl_step[i];
        Aspher_sum[i] += Aspher_step[i];
        eigen_sum[i].x += eigen_step[i].x;
        eigen_sum[i].y += eigen_step[i].y;
        eigen_sum[i].z += eigen_step[i].z;
      }
    } //}}}

    // print shape descriptors to output file(s) //{{{
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        char str[LINE];
        sprintf(str, "%s%s.txt", output, MoleculeType[i].Name);
        FILE *out;
        if ((out = fopen(str, "a")) == NULL) {
          ErrorFileOpen(str, 'a');
          exit(1);
        }
        fprintf(out, "%5d", count_vcf);
        fprintf(out, " %8.5f", Rg_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", sqrRg_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Anis_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Acyl_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Aspher_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", eigen_step[i].x/MoleculeType[i].Number);
        fprintf(out, " %8.5f", eigen_step[i].y/MoleculeType[i].Number);
        fprintf(out, " %8.5f", eigen_step[i].z/MoleculeType[i].Number);
        putc('\n', out);
        fclose(out);
      }
    } //}}}

    // free memory //{{{
    free(Rg_step);
    free(sqrRg_step);
    free(Anis_step);
    free(Acyl_step);
    free(Aspher_step);
    free(eigen_step); //}}}

    // exit the while loop if there's no more coordinates or -e step was reached
    if (LastStep(vcf, NULL) || end == count_vcf) {
      break;
    }
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

  // write simple averages to <output> //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      FILE *out;
      char str[LINE]; // file name <output><mol_name>.txt
      sprintf(str, "%s%s.txt", output, MoleculeType[i].Name);
      if ((out = fopen(str, "a")) == NULL) {
        ErrorFileOpen(str, 'a');
        exit(1);
      }
      // write legend line to output file //{{{
      fprintf(out, "# %s\n", MoleculeType[i].Name);
      fprintf(out, "# simple averages (steps %d to %d): ", start, count_vcf);
      int counter = 1;
      fprintf(out, "(%d) <Rg>, ", counter++);
      fprintf(out, "(%d) <Anis>, ", counter++);
      fprintf(out, "(%d) <Acyl>, ", counter++);
      fprintf(out, "(%d) <Aspher>", counter++);
      fprintf(out, "(%d) <eigen.x>, ", counter++);
      fprintf(out, "(%d) <eigen.y>, ", counter++);
      fprintf(out, "(%d) <eigen.z>, ", counter++);
      putc('\n', out); //}}}

      // write averages to output file //{{{
      putc('#', out);
      fprintf(out, " %lf", Rg_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", Anis_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", Acyl_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", Aspher_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", eigen_sum[i].x/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", eigen_sum[i].y/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", eigen_sum[i].z/(count*MoleculeType[i].Number));
      putc('\n', out); //}}}
      fclose(out);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(Rg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(stuff); //}}}

  return 0;
}
