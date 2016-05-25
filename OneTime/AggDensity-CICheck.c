#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>         input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <width>             width of a single bin\n");
  fprintf(stderr, "   <output.rho>        output density file (automatic ending 'agg#.rho' added)\n");
  fprintf(stderr, "   <agg sizes>         aggregate sizes to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j               specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -n <average>     number of bins to average\n");
  fprintf(stderr, "      -i <name>        use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>        file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v               verbose output\n");
  fprintf(stderr, "      -V               verbose output with more information\n");
  fprintf(stderr, "      -h               print this help and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("AggDensity utility calculates bead density for aggregates of given size(s)  \n");
      printf("from their center of mass. Beside unbonded beads it takes into account only \n");
      printf("beads from the current aggregate, not from any other aggregate.             \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>\n\n", argv[0]);

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <input.agg>         input filename with information about aggregates (agg format)\n");
      printf("   <width>             width of a single bin\n");
      printf("   <output.rho>        output density file (automatic ending 'agg#.rho' added)\n");
      printf("   <agg sizes>         aggregate sizes to calculate density for\n");
      printf("   <options>\n");
      printf("      -j               specify that aggregates with joined coordinates are used\n");
      printf("      -n <average>     number of bins to average\n");
      printf("      -i <name>        use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>        file containing bond alternatives to FIELD\n");
      printf("      -v               verbose output\n");
      printf("      -V               verbose output with more information\n");
      printf("      -h               print this help and exit\n");
      exit(0);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  printf("\n\n"); //}}}

  // check if correct number of arguments //{{{
  if (argc < 6) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -j option - coordinates are joined //{{{
  bool joined = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      joined = true;
    }
  } //}}}

  // -n option - number of bins to average //{{{
  int avg = 1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {

      // Error - non-numeric argument
      if (argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      avg = atoi(argv[i+1]);
    }
  } //}}}

  // -i <name> option - filename of input structure file //{{{
  char vsf_file[32];
  vsf_file[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      // check if .vsf ending is present
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(vsf_file, argv[i+1]);
    }
  }

  // -i option is not used
  if (vsf_file[0] == '\0') {
    strcpy(vsf_file, "dl_meso.vsf");
  } //}}}

  // -b <name> option - filename of input bond file //{{{
  char bonds_file[32];
  bonds_file[0] = '\0'; // check if -b option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-b") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-b' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(bonds_file, argv[i+1]);
    }
  } //}}}

  // -v option - verbose output //{{{
  bool verbose = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      verbose = true;

      break;
    }
  } //}}}

  // -V option - verbose output with comments from input .vcf file //{{{
  bool verbose2 = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-V") == 0) {
      verbose = true;
      verbose2 = true;

      break;
    }
  } //}}}

  int count = 0; // count mandatory arguments

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input.vcf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(input_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<input.vcf> '%s' does not have .vcf ending!\n", input_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <input.agg> - filename of input file with aggregate information //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]); //}}}

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <width>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[16];
  strcpy(output_rho, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
  }

  int aggs = 0;

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (argv[count][0] < '1' || argv[count][0] > '9') {
      fprintf(stderr, "Non-numeric option in <agg sizes>!\n");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[32];

    sprintf(str, "%s%d.rho", output_rho, agg_sizes[aggs][0]);
    if ((out = fopen(str, "w")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", str);
      exit(1);
    }

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print bead type names to output file //{{{
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      fprintf(out, " %s", BeadType[i].Name);
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    aggs++; // number of aggregate sizes
  } //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
  double distance;
  fscanf(agg, "%*s %*s %lf", &distance);

  // skip <contacts> and <output.agg> in Aggregates command
  fscanf(agg, "%*s %*s");

  // read <type names> from Aggregates command //{{{
  int test;
  while ((test = getc(agg)) != '-') {
    ungetc(test, agg);

    char name[10];
    fscanf(agg, "%s ", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", name, input_vcf);
      exit(1);
    }

    BeadType[type].Use = true;
  } //}}}

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf);
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read string from '%s' file!\n", input_vcf);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;
  // skip blank line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    printf("   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // number of bins //{{{
  double max_dist = 0.5 * Min3(BoxLength.x, BoxLength.y, BoxLength.z);
  int bins = ceil(max_dist / width); //}}}

  // allocate memory for density arrays //{{{
  double ***rho = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = malloc(aggs*sizeof(double *));
    for (int j = 0; j < aggs; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    printf("Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      printf(" %d", agg_sizes[i][0]);
    }
    putchar('\n');
  } //}}}

  // counterion arrays
  int *CounterP = calloc(aggs,sizeof(int));
  int *CounterM = calloc(aggs,sizeof(int));
  // number of unimers
  int Uni[2];
  Uni[0] = 0;
  Uni[1] = 0;

  // main loop //{{{
  count = 0; // count timesteps
  int other = 0;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    fflush(stdout);
    printf("\rStep: %6d", ++count);

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {

      if (Aggregate[i].nMolecules == 1) {
        if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[0]].Type].Name, "Diblock") == 0) {
          Uni[0]++;
        } else if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[0]].Type].Name, "Homo") == 0) {
          Uni[1]++;
        } else {
          fprintf(stderr, "Error - uni neither 'Diblock' nor 'Homo'!\n");
          exit(1);
        }
      }

      // test if aggregate size should be used //{{{
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j][0] == Aggregate[i].nMolecules) {
          correct_size = j;
        }
      } //}}}

      int Diblock = 0, Homo = 0;

      // go through all molecules in an aggregate to determine number of diblocks
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name,"Diblock") == 0) {
          Diblock++;
        } else if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name,"Homo") == 0) {
          Homo++;
        }
      }

      if (Homo == 20) {
        agg_sizes[0][1]++;

        other += Diblock;

        Vector com = CenterOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead, BeadType);

        // aggregate beads
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          Vector dist = Distance(Bead[Aggregate[i].Bead[j]].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[Aggregate[i].Bead[j]].Type][0][k]++;
          }
        }

        // monomeric beads
        for (int j = 0; j < Counts.Unbonded; j++) {
          Vector dist = Distance(Bead[j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[j].Type][0][k]++;
          }
        }
      }

      if (correct_size != -1) {
        // calculate number of near counterions
        for (int j = 0; j < Counts.Unbonded; j++) {
          if (BeadType[Bead[j].Type].Charge != 0) {
            bool near = false;

            for (int k = 0; k < Aggregate[i].nBeads && !near; k++) {
              if (BeadType[Bead[Aggregate[i].Bead[k]].Type].Charge != 0) {
                Vector dist = Distance(Bead[j].Position, Bead[Aggregate[i].Bead[k]].Position, BoxLength);
                dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

                if (dist.x < 1) {
                  near = true;

                  if (BeadType[Bead[j].Type].Charge > 0) {
                    CounterP[correct_size]++;
                  } else {
                    CounterM[correct_size]++;
                  }
                }
              }
            }
          }
        }

//      agg_sizes[correct_size][1]++;
      }
    } //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }
  fclose(vcf);
  fclose(agg);

  fflush(stdout);
  printf("\rLast Step: %6d\n", count); //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < aggs && i < 1; i++) {
    FILE *out;
    char str[32];

    sprintf(str, "%s%d.rho", output_rho, agg_sizes[i][0]);
    if ((out = fopen(str, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", str);
      exit(1);
    }

    // calculate rdf
    for (int j = 0; j < bins; j++) {

      // calculate volume of every shell that will be averaged
      double shell[avg];
      for (int k = 0; k < avg; k++) {
        shell[k] = 4 * PI * CUBE(width) *(CUBE(j+k+1) - CUBE(j+k)) / 3;
      }

      fprintf(out, "%.2f", width*(j+0.5*avg));

      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        double temp = 0;

        // sum rdfs from all shells to be averaged
        for (int l = 0; l < avg; l++) {
          temp += rho[k][i][j+l] / (shell[l] * agg_sizes[i][1]);
        }

        // print average value to output file
        fprintf(out, " %10f", temp/avg);
      }
      putc('\n',out);
    }

    fclose(out);

    // print average number of counterions
    int CIp = 0, CIm = 0;
    for (int j = 0; j < bins && (width*j) < 5; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        if (strcmp(BeadType[k].Name, "CounterM") == 0) {
          CIm += rho[k][i][j];
        } else if (strcmp(BeadType[k].Name, "CounterP") == 0) {
          CIp += rho[k][i][j];
        }
      }
    }

    printf("A_S = %d (%d)\n", agg_sizes[i][0], agg_sizes[i][1]);
    printf("The other: %lf\n", (double)(other)/agg_sizes[i][1]);
    printf("   CI+ : %lf CI- : %lf\n", (double)(CounterP[i])/agg_sizes[i][1], (double)(CounterM[i])/agg_sizes[i][1]);
  } //}}}

  printf("Diblock unimers: %lf per step\n", (double)(Uni[0])/count);
  printf("Homopol unimers: %lf per step\n", (double)(Uni[1])/count);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
  }
  free(agg_sizes);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < aggs; j++) {
      free(rho[i][j]);
    }
    free(rho[i]);
  }
  free(rho);
  free(CounterP);
  free(CounterM); //}}}

  return 0;
}
