#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
PotentialAggregates calculates electrostatic potential for aggregates of \
specified size from their centre of mass. It uses Coulomb potential beyond \
short-ranged cut-off and a potential for exponentially smeared charged \
clouds for shorter distances. Parameters for the potential are hard-coded \
in the source code for now. The utility takes into account periodic images \
of the simulation box. Note that the calculation is extremely slow, \
because the utility does not use any Ewald sum-based method.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <input.agg> <width> <output.txt> <agg sizes> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <input.agg>       input agg file\n");
  fprintf(ptr, "   <width>           width of a single bin\n");
  fprintf(ptr, "   <output>          output file with electrostatic potential\n");
  fprintf(ptr, "   <agg size(s)>     aggregate size(s) to calculate density for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 5; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
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
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-m") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // elstat parameters
  double bjerrum = 1.1,
         lambda = 0.2,
         r_c = 3.0,
         beta = (5 * r_c) / (8 * lambda),
         images = 5;
  int point_density = 2;
  int max_points = 100;

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  char *input_vsf = calloc(LINE,sizeof(char));
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // ending timestep //{{{
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  } //}}}

  // error if ending step is lower than starging step //{{{
  if (end != -1 && start > end) {
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <input.agg> - filename of input file with aggregate information //{{{
  char input_agg[LINE];
  strcpy(input_agg, argv[++count]);

  // test if <input.agg> ends with '.agg'
  ext = 1;
  strcpy(extension[0], ".agg");
  if (!ErrorExtension(input_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - filename with bead densities //{{{
  char output_elstat[LINE];
  strcpy(output_elstat, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // '-m' option //{{{
  int *specific_moltype_for_size;
  specific_moltype_for_size = malloc(Counts.TypesOfMolecules*sizeof(int *));
  // all are to be used without '-m' option
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", &specific_moltype_for_size, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
  }

  int aggs = 0;

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (argv[count][0] < '1' || argv[count][0] > '9') {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    aggs++; // number of aggregate sizes
  } //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  // open for the first time to read the distance and type names
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
  double distance;
  fscanf(agg, "%*s %*s %lf", &distance);

  // skip <contacts> and <output.agg> in Aggregates command
  fscanf(agg, "%*s %*s");

  // read <type names> from Aggregates command //{{{
  int test;
  // reading ends if next argument (beginning with '-') or the following empty line is read
  while ((test = getc(agg)) != '-' && test != '\n') {
    ungetc(test, agg);

    char name[LINE];
    fscanf(agg, "%s", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\nError: bead type '%s' is not in %s file\n\n", name, input_coor);
      exit(1);
    }

    BeadType[type].Use = true;

    while ((test = getc(agg)) == ' ')
      ;
    ungetc(test, agg);
  } //}}}
  fclose(agg);

  // open again for production run - to ensure the pointer position in file is correct (at first 'Step')
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  // skip line with command to produce the agg file
  while (getc(agg) != '\n')
    ;
  // skip empty line
  while (getc(agg) != '\n')
    ; //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[LINE];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s\n\n", input_coor);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // number of bins
  int bins = ceil(Min3(BoxLength.x, BoxLength.y, BoxLength.z) / (3 * width));

  double max_dist = (0.5 + images) * Min3(BoxLength.x, BoxLength.y, BoxLength.z);

  // allocate memory for density arrays //{{{
  double **elstat_potential = malloc(aggs*sizeof(double *));
  double **elstat_potential_sqr = malloc(aggs*sizeof(double *));
  for (int i = 0; i < aggs; i++) {
    elstat_potential[i] = calloc(bins,sizeof(double));
    elstat_potential_sqr[i] = calloc(bins,sizeof(double));
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

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
    VerboseOutput(input_coor, Counts, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      fprintf(stdout, " %d", agg_sizes[i][0]);
    }
    putchar('\n');
  } //}}}

  // skip first start-1 steps //{{{
  count = 0;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %6d", count);
    } //}}}

    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      fprintf(stderr, "\nError: premature end of %s file (after %d. step - '%s')\n\n", input_agg, count, stuff);
      test = '\0';
      exit(1);
    }
    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file (%d. step - '%s')\n\n", input_coor, --count, stuff);
      exit(1);
    }
  }

  // print number of starting step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %6d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rStarting step: %6d   \n", start);
    }
  } //}}}
  //}}}

  // create array of bead indices //{{{
  int *index = calloc(Counts.BeadsInVsf,sizeof(int));
  for (int i = 0; i < Counts.Beads; i++) {
    index[Bead[i].Index] = i;
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int count_vcf = start - 1;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %6d", count_vcf);
    } //}}}

    // read aggregates //{{{
    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      count_vcf--;
      fprintf(stderr, "\nError: premature end of %s file (after %d. step - '%s')\n", input_agg, count_vcf, stuff);
      break;
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      ErrorCoorRead(input_coor, test, count_vcf, stuff, input_vsf);
      exit(1);
    } //}}}

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate potential //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {

      // allocate memory for temporary elstat array
      double *temp_elstat = calloc(bins,sizeof(double));

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

      if (correct_size != -1) {
        // zeroize temporary array //{{{
        for (int j = 0; j < bins; j++) {
          temp_elstat[j] = 0;
        } //}}}

        Vector com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead, BeadType);

        // move beads so that com is in the box's centre
        for (int j = 0; j < Counts.Beads; j++) {
          Bead[j].Position.x -= com.x;
          Bead[j].Position.y -= com.y;
          Bead[j].Position.z -= com.z;
        }

        for (int j = 1; j < bins; j++) {
          double r = width * j;

          // start calculations from distance 0.1
          if (r > 0.1) {
            int N = point_density * 4 * PI * SQR(r);
            if (N > max_points) {
              N = max_points;
            }
//          N = 1000;

            int N_count = 0; // number of points

            double a = 4 * PI / N;
            double d = sqrt(a);
            int M1 = round(PI / d);
            double d1 = PI / M1;
            double d2 = a / d1;

            for (int mm = 0; mm < M1; mm++) {
              double angle1 = PI * (mm + 0.5) / M1;
              int M2 = round(2 * PI * sin(angle1) / d2);
              for (int nn = 0; nn < M2; nn++) {
                double angle2 = 2 * PI * nn / M2;
                Vector point;
                point.x = r * sin(angle1) * cos(angle2);
                point.y = r * sin(angle1) * sin(angle2);
                point.z = r * cos(angle1);
                N_count++;

                // TODO: the calculation goes over all beads; how much quicker
                // would it be, if an array of only charged beads was used? Not
                // much, as the slow part is the calculation itself? Then
                // again, this loop is inside other loops.
                for (int l = 0; l < Counts.Beads; l++) {
                  if (BeadType[Bead[l].Type].Charge != 0) {
                    Vector dist;
                    dist = Distance(Bead[l].Position, point, BoxLength);
                    double rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                    double coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                    if (rij < r_c) { // short ranged part
                      temp_elstat[j] += coulomb * (1 - (1 + beta) * exp(-2 * beta * rij));
                    } else { // long ranged part
                      temp_elstat[j] += coulomb;
                    }
                    if (max_dist < 0);

                    // periodic images //{{{
                    Vector dist_orig;
                    dist_orig.x = dist.x;
                    dist_orig.y = dist.y;
                    dist_orig.z = dist.z;
                    for (int m = 1; m <= images; m++) {
                      // +mz
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +0y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -1y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -1y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -1y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0z
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -1y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -1y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -1y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mz
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +0y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                    } //}}}
                  }
                }
              }
            }

            // add from temporary elstat array to global elstat array
            elstat_potential[correct_size][j] += temp_elstat[j] / N_count;
            elstat_potential_sqr[correct_size][j] += SQR(temp_elstat[j]) / N_count; // for error calculation
          }
        }

        agg_sizes[correct_size][1]++;
      }

      // free temporary elstat array
      free(temp_elstat);
    } //}}}

    if (end == count_vcf)
      break;
  }
  fclose(vcf);
  fclose(agg);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count_vcf);
    }
  } //}}}

  // write elstat to output file(s) //{{{
  FILE *out;
  if ((out = fopen(output_elstat, "w")) == NULL) {
    ErrorFileOpen(output_elstat, 'w');
    exit(1);
  }

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print aggregate sizes //{{{
  fprintf(out, "# (1) distance;");
  for (int i = 0; i < aggs; i++) {
    fprintf(out, " (%d) %d", i+2, agg_sizes[i][0]);
    if (i != (aggs-1)) {
      putc(';', out);
    }
  } //}}}

  // print electrostatic potential
  for (int j = 0; j < bins; j++) {
    if ((width*j) > 0.1) {
      // print distance from aggregate com
      fprintf(out, "%.2f", width*(j+0.5));

      for (int i = 0; i < aggs; i++) {
        // print average value and standard deviation to output file
        fprintf(out, " %lf", elstat_potential[i][j]/agg_sizes[i][1]);
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
  }
  free(agg_sizes);
  for (int i = 0; i < aggs; i++) {
    free(elstat_potential[i]);
    free(elstat_potential_sqr[i]);
  }
  free(elstat_potential_sqr);
  free(specific_moltype_for_size);
  free(index); //}}}

  return 0;
}
