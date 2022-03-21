#include "../AnalysisTools.h"
int *InFile;

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
  fprintf(ptr, "   %s <input> <in.agg> <width> <output> <size(s)> \
[options]\n\n", cmd);

  fprintf(ptr, "   <input>     input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <in.agg>    input agg file\n");
  fprintf(ptr, "   <width>     width of a single bin\n");
  fprintf(ptr, "   <output>    output file with electrostatic potential\n");
  fprintf(ptr, "   <size(s)>   aggregate size(s) to calculate density for\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> \
molecules in an aggregate\n");
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
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-m") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

// TODO: all parameters as options
  // elstat parameters
  double bjerrum = 1.1,
         lambda = 0.2,
         r_c = 3.0,
         beta = (5 * r_c) / (8 * lambda),
         images = 5;
  int point_density = 2;
  int max_points = 100;

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
  // test if <input.agg> ends with '.agg'
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

  // <output> - filename with electrostatic potential
  char output_elstat[LINE];
  snprintf(output_elstat, LINE, "%s", argv[++count]);

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

// TODO: allocation... aaargh
  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
  }
  int aggs = 0;
  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (!IsInteger(argv[count]) || atoi(argv[count]) == 0) {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}
    agg_sizes[aggs][0] = atoi(argv[count]);

    aggs++; // number of aggregate sizes
  } //}}}

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // TODO: will change when the agg format changes (at least when Byline is
  //       added to Aggregates*)
  // open input aggregate file and skip the first two lines
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  char line[LINE];
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  // number of bins
  int bins = ceil(Min3(Box.Length.x, Box.Length.y, Box.Length.z) / (3 * width));
  double max_dist = (0.5 + images) * Min3(Box.Length.x, Box.Length.y,
                                          Box.Length.z);

  // TODO: allocation... Aaargh!
  // allocate memory for density arrays //{{{
  double **elstat_potential = malloc(aggs*sizeof(double *));
  double **elstat_potential_sqr = malloc(aggs*sizeof(double *));
  for (int i = 0; i < aggs; i++) {
    elstat_potential[i] = calloc(bins,sizeof(double));
    elstat_potential_sqr[i] = calloc(bins,sizeof(double));
  } //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
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
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

  count = SkipCoorAggSteps(vcf, input_coor, agg,
                           input_agg, Counts, start, silent);

  // create array of bead indices //{{{
  int *index = calloc(Counts.BeadsTotal,sizeof(int));
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    index[Bead[i].Index] = i;
  } //}}}

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
    // TODO: fractionals - should the calculatio be in fractionals? I guess so.
    //       why though? The aggregates are joined and that's that
    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead,
                   MoleculeType, &Molecule, Index);
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

    // TODO: check + fractionals
    // calculate potential //{{{
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
      if (correct_size != -1) {
        // allocate memory for temporary elstat array
        double *temp_elstat = calloc(bins, sizeof *temp_elstat);
        // TODO: ...com should probably be fine in fractional coordinates, eh?
        VECTOR com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead,
                                  BeadType);
        // move beads so that com is in the box's centre
        for (int j = 0; j < Counts.BeadsCoor; j++) {
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
                VECTOR point;
                point.x = r * sin(angle1) * cos(angle2);
                point.y = r * sin(angle1) * sin(angle2);
                point.z = r * cos(angle1);
                N_count++;

                // TODO: the calculation goes over all beads; how much quicker
                // would it be, if an array of only charged beads was used? Not
                // much, as the slow part is the calculation itself? Then
                // again, this loop is inside other loops.
                for (int l = 0; l < Counts.BeadsCoor; l++) {
                  if (BeadType[Bead[l].Type].Charge != 0) {
                    VECTOR dist;
                    dist = Distance(Bead[l].Position, point, Box.Length);
                    double rij = Length(dist);
                    double coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                    if (rij < r_c) { // short ranged part
                      temp_elstat[j] += coulomb * (1 - (1 + beta) * exp(-2 * beta * rij));
                    } else { // long ranged part
                      temp_elstat[j] += coulomb;
                    }
                    if (max_dist < 0);

                    // periodic images //{{{
                    VECTOR dist_orig;
                    dist_orig.x = dist.x;
                    dist_orig.y = dist.y;
                    dist_orig.z = dist.z;
                    for (int m = 1; m <= images; m++) {
                      // +mz
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +0y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -1y //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -1y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -1y //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z + m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0z
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -1y //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -1y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -1y //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mz
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y + m * Box.Length.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +0y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -my //{{{
                      dist.x = dist_orig.x - m * Box.Length.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -my //{{{
                      dist.x = dist_orig.x + m * Box.Length.x;
                      dist.y = dist_orig.y - m * Box.Length.y;
                      dist.z = dist_orig.z - m * Box.Length.z;
                      rij = Length(dist);
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
        // free temporary elstat array
        free(temp_elstat);
      }
    } //}}}

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

  // write elstat to output file //{{{
  FILE *out;
  if ((out = fopen(output_elstat, "w")) == NULL) {
    ErrorFileOpen(output_elstat, 'w');
    exit(1);
  }
  PrintByline(out, argc, argv);
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
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  FreeAggregate(Counts, &Aggregate);
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
  free(index); //}}}

  return 0;
}
