#include "../AnalysisTools.h"
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
        TODO: -m_id option!!!\n\n\
DensityAggregates utility calculates radial bead density for aggregates \
of given size(s) from their centre of mass. For beads in molecules, \
it takes into account only beads from the current aggregate, not from \
any other aggregate. The utility does not check what molecule type a \
given bead belongs to; therefore, if the same bead type appears in \
more molecule types, the resulting density for that bead type will be \
averaged without regard for the various types of molecule it comes from \
(i.e., if 'mol1' and 'mol2' both both contain bead 'A', there will be \
only one column for 'A' bead type).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <in.agg> <width> <output> <size(s)> \
[options]\n\n", cmd);

  fprintf(ptr, "   <input>     input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <in.agg>    input agg file\n");
  fprintf(ptr, "   <width>     width of a single bin\n");
  fprintf(ptr, "   <output>    output density file - one per agg size \
(automatic '<size>.rho' ending)\n");
  fprintf(ptr, "   <size(s)>   aggregate sizes to calculate density for\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> \
molecules in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      -m_id <int>    calculate only for aggregate containing \
the <int> molecule (by resid numbering in vsf)\n");
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
  int req_args = 5; //}}}

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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-m_id") != 0 ) {

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

  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  snprintf(input_agg, LINE, "%s", argv[++count]);
  // test if <in.agg> ends with '.agg'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (!IsPosReal(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - filename with bead densities
  char output_rho[LINE];
  snprintf(output_rho, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");
  int start, end;
  StartEndTime(argc, argv, &start, &end); //}}}

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

  // '-m_id' option (TODO: consider what it should override) //{{{
  int m_id = -1; // no -m_id option
  if (IntegerOption(argc, argv, "-m_id", &m_id)) {
    exit(1);
  } //}}}

  // '-m' option //{{{
  int specific_moltype_for_size[Counts.TypesOfMolecules];
  // set all to be used when '-m' option is missing
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", specific_moltype_for_size,
                          Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

// TODO: allocation... aaargh
  // <agg sizes> - aggregate sizes for calculation //{{{
//int **agg_sizes = malloc(Counts.Molecules * sizeof(int *));
  int (*agg_sizes)[2] = malloc(Counts.Molecules * sizeof *agg_sizes);
//int **agg_mols = malloc(Counts.Molecules * sizeof(int *));
  int **agg_mols = malloc(Counts.Molecules * sizeof **agg_mols);
  for (int i = 0; i < Counts.Molecules; i++) {
//  agg_sizes[i] = calloc(2,sizeof(int));
    agg_sizes[i][0] = 0;
    agg_sizes[i][1] = 0;
//  agg_mols[i] = calloc(Counts.TypesOfMolecules, sizeof(int));
    agg_mols[i] = calloc(Counts.TypesOfMolecules, sizeof *agg_mols[i]);
  }
  int aggs = 0;
  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (!IsInteger(argv[count]) || atoi(argv[count]) == 0) {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}
    agg_sizes[aggs][0] = atoi(argv[count]);
    // ensure output string isn't too long for attaching <size>.vcf
    if (agg_sizes[aggs][0] < 10) {
      output_rho[LINE-1-4] = '\0';
    } else if (agg_sizes[aggs][0] < 100) {
      output_rho[LINE-2-4] = '\0';
    } else if (agg_sizes[aggs][0] < 1000) {
      output_rho[LINE-3-4] = '\0';
    } else if (agg_sizes[aggs][0] < 10000) {
      output_rho[LINE-4-4] = '\0';
    } else {
      output_rho[LINE-100] = '\0';
    }
    aggs++; // number of aggregate sizes
  } //}}}

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // TODO: will change when the agg format changes (at least when Byline is
  //       added to Aggregates*)
  // open input aggregate file and skip the first two lines
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO go for while(fgetc()!='\n'); treatment
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  // number of bins
  double max_dist = 0.5 * Max3(Box.Length.x, Box.Length.y, Box.Length.z);
  int bins = ceil(max_dist / width);

  // TODO: allocation... Aaargh!
  // allocate memory for density arrays //{{{
  double ***rho = malloc(Counts.TypesOfBeads*sizeof(double **));
  double ***rho_2 = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = malloc(aggs*sizeof(double *));
    rho_2[i] = malloc(aggs*sizeof(double *));
    for (int j = 0; j < aggs; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
      rho_2[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules, sizeof (AGGREGATE));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded, sizeof *Aggregate[i].Monomer);
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded, sizeof *Aggregate[i].Bead);
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules, sizeof *Aggregate[i].Molecule);
  } //}}}

  if (verbose) { //{{{
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  SkipVtfStructure(vcf, struct_lines);

  count = SkipCoorAggSteps(vcf, input_coor, agg,
                           input_agg, Counts, start, silent);

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

    // TODO: will change (probably)
    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead,
                   MoleculeType, &Molecule, Index);
    // TODO: fractionals - should the calculatio be in fractionals? I guess so.
    //       why though? The aggregates are joined and that's that; monomeric
    //       beads?
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    if (!joined) {
      RemovePBCMolecules(Counts, Box, BeadType, &Bead,
                         MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, Box.Length,
                          BeadType, &Bead, MoleculeType, Molecule);
    }

    // TODO: allocation... Aaargh!
    // allocate memory for temporary density arrays //{{{
    double ***temp_rho = malloc(Counts.TypesOfBeads*sizeof(double **));
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      temp_rho[i] = malloc(aggs*sizeof(double *));
      for (int j = 0; j < aggs; j++) {
        temp_rho[i][j] = calloc(bins,sizeof(double));
      }
    } //}}}

  // TODO: check
    // calculate densities //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {
      // test if aggregate 'i' should be used //{{{
      int size = 0;
      // agg size = number of molecules of type 'specific_moltype_for_size'
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
        }
      }
      // is 'size' in provided list?
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j][0] == size) {
          correct_size = j;
        }
      } //}}}

      // if '-m_id' is used, check if specified resid from vsf is in aggregate
      if (m_id != -1) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (m_id == (id+1)) { // resname in vsf start from 1
            correct_size = 0;
            break;
          }
        }
      }

      if (correct_size != -1) {
        VECTOR com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead, BeadType);

        // free temporary density array //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < aggs; k++) {
            for (int l = 0; l < bins; l++) {
              temp_rho[j][k][l] = 0;
            }
          }
        } //}}}

        // aggregate beads //{{{
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int bead = Aggregate[i].Bead[j];
          int mol = Bead[bead].Molecule;
          int moltype = Molecule[mol].Type;

          if (MoleculeType[moltype].Use) {
            // TODO: fractionals!
            VECTOR dist = Distance(Bead[Aggregate[i].Bead[j]].Position,
                                   com, Box.Length);
            dist.x = Length(dist);

            if (dist.x < max_dist) {
              int k = dist.x / width;

              temp_rho[Bead[Aggregate[i].Bead[j]].Type][correct_size][k]++;
            }
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Counts.Unbonded; j++) {
          // fractionals!
          VECTOR dist = Distance(Bead[j].Position, com, Box.Length);
          dist.x = Length(dist);

          if (dist.x < max_dist) {
            int k = dist.x / width;

            temp_rho[Bead[j].Type][correct_size][k]++;
          }
        } //}}}

        agg_sizes[correct_size][1]++;

        // add from temporary density array to global density arrays //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < bins; k++) {
            rho[j][correct_size][k] += temp_rho[j][correct_size][k];
            rho_2[j][correct_size][k] += SQR(temp_rho[j][correct_size][k]);
          }
        } //}}}

        // count mol types in the aggregate
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          agg_mols[correct_size][Molecule[Aggregate[i].Molecule[j]].Type]++;
        }
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      for (int j = 0; j < aggs; j++) {
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
  fclose(agg);
  // print last step count?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < aggs; i++) {
    FILE *out;
    char str[LINE];
    // assemble correct name
    snprintf(str, LINE, "%s%d.rho", output_rho, agg_sizes[i][0]);
    out = OpenFile(str, "a");
    // write initial stuff to the file //{{{
    PrintByline(out, argc, argv);
    // print bead type names to output file
    fprintf(out, "# for each bead type: (1) rdp; (2) stderr; (3) rnp; (4) stderr\n");
    fprintf(out, "# columns: (1) distance;");
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      fprintf(out, " (%d) %s", 4*i+2, BeadType[j].Name);
      if (j != (Counts.TypesOfBeads-1)) {
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    // calculate rdf
    for (int j = 0; j < bins; j++) {
      double shell = 4 * PI * CUBE(width) * (CUBE(j+1) - CUBE(j)) / 3;
      fprintf(out, "%.2f", width*(2*j+1)/2);

      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        // print rnp/rdp to file only if there's something the given bin
        if (rho[k][i][j] > 0) {
          // radial number profile
          double temp_rnp = rho[k][i][j] / agg_sizes[i][1];
          double temp_rnp_err = rho_2[k][i][j] / agg_sizes[i][1];
          // radial density profile
          double temp_rdp = temp_rnp / shell;
          double temp_rdp_err = temp_rnp_err / shell;
          // errors
          temp_rdp_err = sqrt(temp_rdp_err - temp_rdp);
          temp_rnp_err = sqrt(temp_rnp_err - temp_rnp);
          // print average value to output file
          fprintf(out, " %10f %10f", temp_rdp, temp_rdp_err);
          fprintf(out, " %10f %10f", temp_rnp, temp_rnp_err);
        } else { // otherwise print question marks (for gnuplot)
          fprintf(out, "          ?          ?");
          fprintf(out, "          ?          ?");
        }
      }
      putc('\n',out);
    }

    fprintf(out, "# (1) Number of aggregates; Average numbers of molecules:");
    count = 2;
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      fprintf(out," (%d) %s", count++, MoleculeType[j].Name);
      if (i < (Counts.TypesOfMolecules-1)) {
        putc(';', out);
      }
    }
    putc('\n', out);
    putc('#', out);
    fprintf(out," %d", agg_sizes[i][1]);
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      fprintf(out," %lf", (double)(agg_mols[i][j])/agg_sizes[i][1]);
    }
    putc('\n', out);

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(Counts, &Aggregate);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
    free(agg_mols[i]);
  }
  free(agg_sizes);
  free(agg_mols);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < aggs; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho);
  free(rho_2); //}}}

  return 0;
}
