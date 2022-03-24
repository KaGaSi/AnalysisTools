#include "../AnalysisTools.h"
int *InFile;

// TODO: change double rho arrays into (long) int?

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
DensityMolecules utility calculates number density for all bead types from \
the centre of mass (or any bead in a molecule) of specified molecules \
similarly to how DensityAggregates calculates the density for \
aggregates.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol(s)> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <width>    width of a single distribution bin\n");
  fprintf(ptr, "   <output>   output density file \
(automatic ending '<molecule_name>.rho' added)\n");
  fprintf(ptr, "   <mol(s)>   molecule name(s) to calculate density for \
(optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --all            use all molecules \
(overwrites <mol(s)>)\n");
  fprintf(ptr, "      --joined         specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -st <int>        starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>         ending timestep for calculation\n");
  fprintf(ptr, "      -c <out> <int>   use <int>-th molecule bead instead \
of centre of mass\n");
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
  int req_args = 4; //}}}

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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-c") != 0) {

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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - filename with bead densities //{{{
  // leave space for <mol_name>.rho ending within LINE length
  char output_rho[LINE-MOL_NAME-4] = "";
  snprintf(output_rho, LINE-MOL_NAME-4, "%s", argv[++count]); //}}}

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
      int type = FindMoleculeType(argv[count], Counts, MoleculeType);
      // error - nonexistent molecule  //{{{
      if (type == -1) {
        ErrorPrintError_old();
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", input_coor);
        ColourText(STDERR_FILENO, RED);
        fprintf(stderr, " - non-existent molecule type ");
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s\n", argv[count]);
        ColourReset(STDERR_FILENO);
        ErrorMoleculeType(Counts, MoleculeType);
        exit(1);
      } //}}}
      MoleculeType[type].Use = true;
    }
  } else { // --all option is used
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      MoleculeType[i].Use = true;
    }
  } //}}}

  // -c option - specify which bead to use as a molecule centre //{{{
  // array for considering whether to use COM or specified bead number
  int centre[Counts.TypesOfMolecules];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    centre[i] = -1; // helper value
  }
  // the -c option can be used multiple times
  for (int i = 1; i < argc; i++) {
    int mtype = -1, id = -1;
    if (MoleculeTypeIntOption(argc, i, argv, "-c", &mtype, &id,
          Counts, MoleculeType)) {
      exit(1);
    }
    if (mtype > -1) {
      if (!MoleculeType[mtype].Use) {
        // warning - this molecule not specified in <mol(s)> //{{{
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, "\nWarning: ");
        ColourText(STDERR_FILENO, CYAN);
        fprintf(stderr, "-c");
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, " - molecule ");
        ColourText(STDERR_FILENO, CYAN);
        fprintf(stderr, "%s", MoleculeType[mtype].Name);
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, " not specified in <mol(s)>;");
        fprintf(stderr, " this option will be ignored\n");
        ColourReset(STDERR_FILENO); //}}}
      } else if (id >= MoleculeType[mtype].nBeads) {
        // warning - too high an id; using the last bead //{{{
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, "\nWarning: ");
        ColourText(STDERR_FILENO, CYAN);
        fprintf(stderr, "-c");
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, " - index ");
        ColourText(STDERR_FILENO, CYAN);
        fprintf(stderr, "%d", id);
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, " is larger than the number of beads in molecule ");
        ColourText(STDERR_FILENO, CYAN);
        fprintf(stderr, "%s", MoleculeType[mtype].Name);
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, "; last bead will be used insted (i.e., index");
        ColourText(STDERR_FILENO, CYAN);
        fprintf(stderr, "%d", MoleculeType[mtype].nBeads-1);
        ColourText(STDERR_FILENO, YELLOW);
        fprintf(stderr, ")\n");
        centre[mtype] = MoleculeType[i].nBeads - 1; //}}}
      } else {
        // 'proper' index used
        centre[mtype] = id - 1; // the index should start from zero
      }
    }
  } //}}}

  // number of bins //{{{
  double max_dist = 0.5 * Min3(Box.Length.x, Box.Length.y, Box.Length.z);
  int bins = ceil(max_dist / width); //}}}

// TODO: sizeof ...argh!
  // allocate memory for density arrays //{{{
//double ***rho = malloc(Counts.TypesOfMolecules * sizeof(double **));
//double ***rho_2 = malloc(Counts.TypesOfMolecules * sizeof(double **));
  double ***rho = malloc(Counts.TypesOfBeads * sizeof ***rho);
  double ***rho_2 = malloc(Counts.TypesOfBeads * sizeof ***rho_2);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  rho[i] = malloc(Counts.TypesOfMolecules * sizeof(double *));
//  rho_2[i] = malloc(Counts.TypesOfMolecules * sizeof(double *));
    rho[i] = malloc(Counts.TypesOfMolecules * sizeof **rho[i]);
    rho_2[i] = malloc(Counts.TypesOfMolecules * sizeof **rho_2[i]);
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
//    rho[i][j] = calloc(bins,sizeof(double));
//    rho_2[i][j] = calloc(bins,sizeof(double));
      rho[i][j] = calloc(bins, sizeof *rho[i][j]);
      rho_2[i][j] = calloc(bins, sizeof *rho_2[i][j]);
    }
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

  // main loop //{{{
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count++;
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

  // TODO: sizeof ...argh!
    // allocate memory for temporary density arrays //{{{
//  double ***temp_rho = malloc(Counts.TypesOfBeads*sizeof(double **));
    double ***temp_rho = malloc(Counts.TypesOfBeads * sizeof ***temp_rho);
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
//    temp_rho[i] = malloc(Counts.TypesOfMolecules*sizeof(double *));
      temp_rho[i] = malloc(Counts.TypesOfMolecules * sizeof **temp_rho[i]);
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
//      temp_rho[i][j] = calloc(bins,sizeof(double));
        temp_rho[i][j] = calloc(bins,sizeof *temp_rho[i][j]);
      }
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      if (MoleculeType[mtype].Use) {
        // determine centre to calculate densities from //{{{
        VECTOR com;
        if (centre[mtype] == -1 ) { // use molecule's centre of mass
          com = CentreOfMass(MoleculeType[mtype].nBeads, Molecule[i].Bead,
                             Bead, BeadType);
        } else { // use centre[mtype]-th molecule's bead as com
          int id = Molecule[i].Bead[centre[mtype]];
          com.x = Bead[id].Position.x;
          com.y = Bead[id].Position.y;
          com.z = Bead[id].Position.z;
        } //}}}

        // zeroize temporary density array //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            for (int l = 0; l < bins; l++) {
              temp_rho[j][k][l] = 0;
            }
          }
        } //}}}

        // molecule beads //{{{
        for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
          int id = Molecule[i].Bead[j];

          VECTOR dist = Distance(Bead[id].Position, com, Box.Length);
          dist.x = Length(FromFractional(dist, Box));
          if (dist.x < max_dist) {
            int k = dist.x / width;
            temp_rho[Bead[id].Type][mtype][k]++;
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Counts.Unbonded; j++) {
          VECTOR dist = Distance(Bead[j].Position, com, Box.Length);
          dist.x = Length(FromFractional(dist, Box));
          if (dist.x < max_dist) {
            int k = dist.x / width;
            temp_rho[Bead[j].Type][mtype][k]++;
          }
        } //}}}

        // add from temporary density array to global density arrays //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < bins; k++) {
            rho[j][mtype][k] += temp_rho[j][mtype][k];
            rho_2[j][mtype][k] += SQR(temp_rho[j][mtype][k]);
          }
        } //}}}
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        free(temp_rho[i][j]);
      }
      free(temp_rho[i]);
    }
    free(temp_rho); //}}}

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

  // write initial stuff to output density file(s) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      FILE *out;
      char str[LINE];
      sprintf(str, "%s%s.rho", output_rho, MoleculeType[i].Name);
      if ((out = fopen(str, "w")) == NULL) {
        ErrorFileOpen(str, 'w');
        exit(1);
      }
      fclose(out);
    }
  } //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      FILE *out;
      char str[LINE]; // file name <output_rho><mol_name>.rho
      sprintf(str, "%s%s.rho", output_rho, MoleculeType[i].Name);
      if ((out = fopen(str, "a")) == NULL) {
        ErrorFileOpen(str, 'a');
        exit(1);
      }
      // initial stuff to the file //{{{
      PrintByline(out, argc, argv);
      fprintf(out, "# for each bead type: ");
      fprintf(out, "(1) rdp; ");
      fprintf(out, "(2) stderr; ");
      fprintf(out, "(3) rnp; ");
      fprintf(out, "(4) stderr\n");
      fprintf(out, "# columns: (1) distance;");
      for (int i = 0; i < Counts.TypesOfBeads; i++) {
        fprintf(out, " (%d) %s", 4*i+2, BeadType[i].Name);
        if (i != (Counts.TypesOfBeads-1)) {
          putc(';', out);
        }
      }
      putc('\n', out); //}}}

      // calculate rdf
      for (int j = 0; j < bins; j++) {
        // calculate volume of every shell that will be averaged
        double shell;
        shell = 4 * PI * CUBE(width) *(CUBE(j+1) - CUBE(j)) / 3;

        fprintf(out, "%.3f", width*(2*j+1)/2);

        for (int k = 0; k < Counts.TypesOfBeads; k++) {
          // print rnp/rdp to file only if there's something the given bin
          if (rho[k][i][j] > 0) {
            double temp_rdp = 0, temp_rnp = 0,
                   temp_rdp_err = 0, temp_rnp_err = 0;
            // radial number profile
            temp_rnp = rho[k][i][j] / (MoleculeType[i].Number * count);
            temp_rnp_err = rho_2[k][i][j] / (MoleculeType[i].Number * count);
            // radial density profile
            temp_rdp = temp_rnp / shell;
            temp_rdp_err = temp_rnp_err / shell;
            // errors
            temp_rdp_err = sqrt(temp_rdp_err - temp_rdp);
            temp_rnp_err = sqrt(temp_rnp_err - temp_rnp);
            // print to file
            fprintf(out, " %10f %10f", temp_rdp, temp_rdp_err);
            fprintf(out, " %10f %10f", temp_rnp, temp_rnp_err);
          } else { // otherwise print question marks (for gnuplot)
            fprintf(out, "          ?          ?");
            fprintf(out, "          ?          ?");
          }
        }
        putc('\n',out);
      }
      fclose(out);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
