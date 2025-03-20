#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
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
//   fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> "
//           "molecules in an aggregate\n");
//   fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
//   fprintf(ptr, "      -m_id <int>    calculate only for aggregate containing "
//           "the <int> molecule (by resid numbering in vsf)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join;      // --joined
  FILE_TYPE fout; // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 1, count = 0,
      req_arg = 5;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent",
              "--help", "--version", "--joined");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);
  // test if <in.agg> ends with '.agg'
  int ext = 1;
  char extension[2][EXTENSION];
  s_strcpy(extension[0], ".agg", EXTENSION);
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output> - filename with bead densities
  char output_rho[LINE];
  s_strcpy(output_rho, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt->join = false; // joined coordinates supplied, so no need to join
  } else {
    opt->join = true; // molecules need to be joined
  } //}}}
  //}}}

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  double *box = System.Box.Length;

  // <agg sizes> - aggregate sizes for calculation //{{{
  // TODO: proper allocation...
  int (*agg_sizes)[2] = malloc(Count->Molecule * sizeof *agg_sizes);
  int **agg_mols = malloc(Count->Molecule * sizeof *agg_mols);
  for (int i = 0; i < Count->Molecule; i++) {
    agg_sizes[i][0] = 0;
    agg_sizes[i][1] = 0;
    agg_mols[i] = calloc(Count->MoleculeType, sizeof *agg_mols[i]);
  }
  int aggs = 0;
  while (++count < argc && argv[count][0] != '-') {
    // Error - non-numeric argument //{{{
    long val;
    if (!IsNaturalNumber(argv[count], &val)) {
      ErrorNaN("<agg size(s)>");
      Help(argv[0], true, common, option);
      exit(1);
    } //}}}
    agg_sizes[aggs][0] = atoi(argv[count]);
    aggs++; // number of aggregate sizes
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
  } //}}}

  // // TODO: read this from agg file and use?
  // double distance = 1; // <distance> parameter from Aggregate command
  // int contacts = 1; // <contacts> parameter from Aggregate command - not used here

  // open input aggregate file and skip the first two lines
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO: go for while(fgetc()!='\n'); treatment
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  // number of bins
  double max_dist = 0.5 * Max3(box[0], box[1], box[2]);
  int bins = ceil(max_dist / width);

  // TODO: allocation... Aaargh!
  // allocate memory for density arrays //{{{
  double ***rho = malloc(Count->BeadType*sizeof(double **));
  double ***rho_2 = malloc(Count->BeadType*sizeof(double **));
  for (int i = 0; i < Count->BeadType; i++) {
    rho[i] = malloc(aggs*sizeof(double *));
    rho_2[i] = malloc(aggs*sizeof(double *));
    for (int j = 0; j < aggs; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
      rho_2[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0, // count lines in the vcf file
      line_count_agg = 0; // count lines in the agg file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);

    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      // TODO: join aggregates
      ReadAggregates(agg, input_agg, &System, Aggregate, &line_count_agg);

      // TODO: proper allocation...
      // allocate memory for temporary density arrays //{{{
      double ***temp_rho = malloc(Count->BeadType*sizeof(double **));
      for (int i = 0; i < Count->BeadType; i++) {
        temp_rho[i] = malloc(aggs*sizeof(double *));
        for (int j = 0; j < aggs; j++) {
          temp_rho[i][j] = calloc(bins,sizeof(double));
        }
      } //}}}

      // calculate densities //{{{
      for (int i = 0; i < Count->Aggregate; i++) {
        // test if aggregate 'i' should be used //{{{
        // // TODO: -m option (I think) //{{{
        // // agg size = number of molecules of type 'specific_moltype_for_size'
        // int size = 0;
        // for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        //   int mol_type = System.Molecule[Aggregate[i].Molecule[j]].Type;
        //   if (specific_moltype_for_size[mol_type]) {
        //     size++;
        //   }
        // } //}}}
        // is 'size' in provided list?
        int size = Aggregate[i].nMolecules;
        int correct_size = -1;
        for (int j = 0; j < aggs; j++) {
          if (agg_sizes[j][0] == size) {
            correct_size = j;
          }
        } //}}}

        // // TODO: -m_id option //{{{
        // // if '-m_id' is used, check if specified resid from vsf is in aggregate
        // if (m_id != -1) {
        //   for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        //     int id = Aggregate[i].Molecule[j];
        //     if (m_id == (id+1)) { // resname in vsf start from 1
        //       correct_size = 0;
        //       break;
        //     }
        //   }
        // } //}}}

        if (correct_size != -1) {
          double com[3];
          CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, System, com);

          // free temporary density array //{{{
          for (int j = 0; j < Count->BeadType; j++) {
            for (int k = 0; k < aggs; k++) {
              for (int l = 0; l < bins; l++) {
                temp_rho[j][k][l] = 0;
              }
            }
          } //}}}

          // aggregate beads //{{{
          for (int j = 0; j < Aggregate[i].nBeads; j++) {
            int bead = Aggregate[i].Bead[j];
            // // TODO: all the jazz about which moltypes to use
            // int mol = System.Bead[bead].Molecule;
            // int moltype = System.Molecule[mol].Type;
            // if (System.MoleculeType[moltype].Flag) {
              double dist[3];
              Distance(System.Bead[bead].Position, com, box, dist);
              dist[0] = VectLength(dist);

              if (dist[0] < max_dist) {
                int k = dist[0] / width;

                temp_rho[System.Bead[Aggregate[i].Bead[j]].Type][correct_size][k]++;
              }
            // }
          } //}}}

          // monomeric beads //{{{
          for (int j = 0; j < Count->Unbonded; j++) {
            int id = System.Unbonded[j];
            double dist[3];
            Distance(System.Bead[id].Position, com, box, dist);
            dist[0] = VectLength(dist);

            if (dist[0] < max_dist) {
              int k = dist[0] / width;
              temp_rho[System.Bead[id].Type][correct_size][k]++;
            }
          } //}}}

          agg_sizes[correct_size][1]++;

          // add from temporary density array to global density arrays //{{{
          for (int j = 0; j < Count->BeadType; j++) {
            for (int k = 0; k < bins; k++) {
              rho[j][correct_size][k] += temp_rho[j][correct_size][k];
              rho_2[j][correct_size][k] += Square(temp_rho[j][correct_size][k]);
            }
          } //}}}

          // count mol types in the aggregate
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            agg_mols[correct_size][System.Molecule[Aggregate[i].Molecule[j]].Type]++;
          }
        }
      } //}}}

      // free temporary density array //{{{
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = 0; j < aggs; j++) {
          free(temp_rho[i][j]);
        }
        free(temp_rho[i]);
      }
      free(temp_rho); //}}}
    //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
      // TODO: skip agg file step
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  fclose(agg);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < aggs; i++) {
    FILE *out;
    char str[LINE];
    // assemble correct name
    if (snprintf(str, LINE, "%s%d.rho", output_rho, agg_sizes[i][0]) < 0) {
      ErrorSnprintf();
    }
    // write initial stuff to the file //{{{
    PrintByline(str, argc, argv);
    // print bead type names to output file
    out = OpenFile(str, "a");
    fprintf(out, "# for each bead type: (1) rdp; (2) rnp\n");
    fprintf(out, "# columns: (1) distance;");
    for (int j = 0; j < Count->BeadType; j++) {
      fprintf(out, " (%d) %s", 4*i+2, System.BeadType[j].Name);
      if (j != (Count->BeadType-1)) {
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    // calculate rdf
    for (int j = 0; j < bins; j++) {
      double shell = 4 * PI * Cube(width) * (Cube(j+1) - Cube(j)) / 3;
      fprintf(out, "%.2f", width*(2*j+1)/2);

      for (int k = 0; k < Count->BeadType; k++) {
        // print rnp/rdp to file only if there's something the given bin
        if (rho[k][i][j] > 0) {
          // radial number profile
          double temp_rnp = rho[k][i][j] / agg_sizes[i][1];
          // radial density profile
          double temp_rdp = temp_rnp / shell;
          // print average value to output file
          fprintf(out, " %10f", temp_rdp);
          fprintf(out, " %10f", temp_rnp);
        } else { // otherwise print question marks (for gnuplot)
          fprintf(out, "          ?");
          fprintf(out, "          ?");
        }
      }
      putc('\n',out);
    }

    fprintf(out, "# (1) Number of aggregates; Average numbers of molecules:");
    count = 2;
    for (int j = 0; j < Count->MoleculeType; j++) {
      fprintf(out," (%d) %s", count++, System.MoleculeType[j].Name);
      if (i < (Count->MoleculeType-1)) {
        putc(';', out);
      }
    }
    putc('\n', out);
    putc('#', out);
    fprintf(out," %d", agg_sizes[i][1]);
    for (int j = 0; j < Count->MoleculeType; j++) {
      fprintf(out," %lf", (double)(agg_mols[i][j])/agg_sizes[i][1]);
    }
    putc('\n', out);

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  for (int i = 0; i < Count->Molecule; i++) {
    free(agg_mols[i]);
  }
  free(agg_sizes);
  free(agg_mols);
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = 0; j < aggs; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho);
  free(rho_2);
  free(opt); //}}}
  return 0;
}
