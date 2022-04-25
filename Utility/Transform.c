#include "../AnalysisTools.h"
int nanosleep(const struct timespec *req, struct timespec *rem);
int *InFile;
// TODO: option define reflecting plane by three series of bead indices (change
//       <mode>)?

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
    fprintf(ptr, "<mode>:\n");
  } else {
    ptr = stdout;
//bead indices defining the normal to the reflecting plane \n
    fprintf(ptr, "\
Transform utility transforms individual molecules by reflecting, inverting, or \
rotating them. The different operations are defined by the <mode> argument \
that contains indices corresponding to molecule's internal bead indices (as \
given in the input structure file); these indices are, therefore, always \
between 1 and the number of beads in target molecule. For each target \
molecule, the coordinates of these beads then determine the transforming \
elements; whenever more than one index is supplied to define a point, \
an average coordinate is used.\n");
  }
  fprintf(ptr, "\
i) <mode> for reflection: 'ref A <int(s)> B <int(s)> C <int(s)>'\n\
  A <int(s)> & B <int(s)> .. bead indices for two points defining the normal \
to the reflecting plane\
  C <int(s)> .. one or more indices defining a point on that plane\n\
ii) <mode> for inversion: 'inv <int(s)>'\n\
   <int(s)> .. bead indices defining the centre of inversion (if more than \
than one index is given, an averaged coordinate is used)\n\
iii) <mode> for rotation: 'rot A <int(s)> B <int(s)> [N <int(s)>] \
C <int(s)> P <angle>'\n\
    A <int(s)> & B <int(s)> .. bead indices for two points defining rotation \
axis orientation\n\
    [N <int(s)>] .. optional argument for the third point; if it is present, \
the normal to the plane formed by points A, B, and N is used as \
the rotation axis\n\
    C <int(s)> .. one or more indices defining a point lying \
on the rotation axis\n\
    P <angle> .. the rotation angle in degrees\n\
\n\n");

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <output.vcf> <mode> [options]\n\n", cmd);

  fprintf(ptr, "   <input>        input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <output.vcf>   output coordinate file (vcf format)\n");
  fprintf(ptr, "   <mode>         how to transform molecules (see the help \
text or manual for details)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined       are provided coordinates joined?\n");
  fprintf(ptr, "      -n <int(s)>    use specified molecules (vsf indices)\n");
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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-n") != 0) {
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

  // <output.vcf> - output vcf file //{{{
  char output_vcf[LINE] = "";
  snprintf(output_vcf, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vcf'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <mode> - how to transform molecules //{{{
  count++;
  int max_ids = 100,
      mode = -1; // -1..error, 0..ref, 1..inv, 2..rot
  double angle = 0;
  long A[max_ids+1], B[max_ids+1], C[max_ids+1], N[max_ids+1];
  for (int i = 0; i < max_ids; i++) {
    A[i] = -1;
    B[i] = -1;
    C[i] = -1;
    N[i] = -1;
  }
  A[max_ids] = 0;
  B[max_ids] = 0;
  C[max_ids] = 0;
  N[max_ids] = 0;
  // find what mode to use //{{{
  if (strncmp(argv[count],"reflection", 3) == 0) {
    mode = 0;
  } else if (strncmp(argv[count],"inversion", 3) == 0) {
    mode = 1;
  } else if (strncmp(argv[count],"rotation", 3) == 0) {
    mode = 2;
  } //}}}
  int i = 0;
  switch (mode) {
    // <mode> 'ref[lection]' //{{{
    /* the reflecting plane is specified by a normal vector constructed from two
     * bead indices and a point taken as either the third provided bead index or
     * average position from indices 3+ (when more than three are provided)
     * TODO: what about a normal perpendicular to a molecule or some such?
     */
    case 0:
      // first point (normal to the reflecting plane): 'A <int(s)>' //{{{
      if (argv[++count][0] != 'A') {
        strcpy(ERROR_MSG, "<mode> 'reflection': missing 'A' point\n");
        PrintError();
        Help(argv[0], true);
        exit(1);
      }
      for (int i = 0;
           ++count < argc && argv[count][0] != 'B' && i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &A[i]) || A[i] < 1) {
          strcpy(ERROR_MSG, "<mode>: bead index must be positive integer \
(or 'B' is missing)");
          PrintError();
          fprintf(stderr, "%swrong argument for %sA%s point:%s %s%s\n",
                  Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
          Help(argv[0], true);
          exit(1);
        }
        A[i]--;
        A[max_ids]++; // counts the number of indices for A point
      } //}}}
      // second point (normal to the reflecting plane): 'B <int(s)>' //{{{
      for (int i = 0;
           ++count < argc && argv[count][0] != 'C' && i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &B[i]) || B[i] < 1) {
          strcpy(ERROR_MSG, "<mode>: bead index must be positive integer \
(or 'C' is missing)");
          PrintError();
          fprintf(stderr, "%swrong argument for %sB%s point:%s %s%s\n",
                  Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
          Help(argv[0], true);
          exit(1);
        }
        B[i]--;
        B[max_ids]++; // counts the number of indices for B point
      } //}}}
      // third point (on reflecting plane): 'C <int(s)>' //{{{
      for (int i = 0;
           ++count < argc && argv[count][0] != '-' && i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &C[i]) || C[i] < 1) {
          strcpy(ERROR_MSG, "<mode>: bead index must be positive integer");
          PrintError();
          fprintf(stderr, "%swrong argument for %sC%s point:%s %s%s\n",
                  Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
          Help(argv[0], true);
          exit(1);
        }
        C[i]--;
        C[max_ids]++; // counts the number of indices for C point
      } //}}}
      if (A[max_ids] == 0 || B[max_ids] == 0 || C[max_ids] == 0) {
        strcpy(ERROR_MSG, "<mode>: at least one index must be provided for \
each of A, B, and C");
        PrintError();
        Help(argv[0], true);
        exit(1);
      }
      break; //}}}
    // <mode> 'inv[ersion]' //{{{
    /* the centre of inversion is defined as the averaged coordinates of all the
     * given bead indices
     */
    case 1:
      for (i = 0; ++count < argc && argv[count][0] != '-' && i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &A[i]) || A[i] < 0) {
          strcpy(ERROR_MSG, "<mode> 'inversion': at least one \
positive number is required as an argument");
          PrintError();
          Help(argv[0], true);
          exit(1);
        }
        A[i]--;
        A[max_ids]++;
      }
      break; //}}}
    // <mode> 'rot[ation]' //{{{
    /* the centre of inversion is defined as the averaged coordinates of all the
     * given bead indices
     */
    case 2:
      // first point: 'A <int(s)>' //{{{
      if (argv[++count][0] != 'A') {
        strcpy(ERROR_MSG, "<mode> 'rotation': missing 'A' point\n");
        PrintError();
        Help(argv[0], true);
        exit(1);
      }
      for (int i = 0;
           ++count < argc && argv[count][0] != 'B' && i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &A[i]) || A[i] < 1) {
          strcpy(ERROR_MSG, "<mode>: bead index must be positive integer \
(or 'B' is missing)");
          PrintError();
          fprintf(stderr, "%swrong argument for %sA%s point:%s %s%s\n",
                  Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
          Help(argv[0], true);
          exit(1);
        }
        A[i]--;
        A[max_ids]++; // counts the number of indices for A point
      } //}}}
      // second point: 'B <int(s)>' //{{{
      for (int i = 0; ++count < argc &&
           argv[count][0] != 'N' && // optional 'N <int(s)>' or
           argv[count][0] != 'C' && // mandatory 'C <int(s)>' follows
           i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &B[i]) || B[i] < 1) {
          strcpy(ERROR_MSG, "<mode>: bead index must be positive integer \
(or 'P' is missing)");
          PrintError();
          fprintf(stderr, "%swrong argument for %sB%s point:%s %s%s\n",
                  Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
          Help(argv[0], true);
          exit(1);
        }
        B[i]--;
        B[max_ids]++; // counts the number of indices for B point
      } //}}}
      // third point (optional for axis normal to plane): 'N <int(s)>' //{{{
      if (argv[count][0] == 'N') {
        for (int i = 0;
             ++count < argc && argv[count][0] != 'C' && i < max_ids; i++) {
          // error - not a non-negative number
          if (!IsInteger2(argv[count], &N[i]) || N[i] < 1) {
            strcpy(ERROR_MSG, "<mode>: bead index must be positive integer");
            PrintError();
            fprintf(stderr, "%swrong argument for %sN%s point:%s %s%s\n",
                    Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
            Help(argv[0], true);
            exit(1);
          }
          N[i]--;
          N[max_ids]++; // counts the number of indices for N point
        }
        if (N[max_ids] == 0) {
          strcpy(ERROR_MSG, "<mode>: at least one index must be provided for \
the optional N");
          PrintError();
          Help(argv[0], true);
          exit(1);
        }
      } //}}}
      // point on the axis: 'C <int(s)>' //{{{
      for (int i = 0;
           ++count < argc && argv[count][0] != 'P' && i < max_ids; i++) {
        // error - not a non-negative number
        if (!IsInteger2(argv[count], &C[i]) || C[i] < 1) {
          strcpy(ERROR_MSG, "<mode>: bead index must be positive integer \
(or 'P' is missing)");
          PrintError();
          fprintf(stderr, "%swrong argument for %sC%s point:%s %s%s\n",
                  Red(), Yellow(), Red(), Yellow(), argv[count], ColourReset());
          Help(argv[0], true);
          exit(1);
        }
        C[i]--;
        C[max_ids]++; // counts the number of indices for C point
      } //}}}
      // angle: 'P <double>' //{{{
      if (++count >= argc) {
        strcpy(ERROR_MSG, "<mode>: missing rotation angle");
        PrintError();
        Help(argv[0], true);
        exit(1);
      }
      if (!IsReal2(argv[count], &angle)) {
        strcpy(ERROR_MSG, "<mode>: angle must be a real number");
        PrintError();
        fprintf(stderr, "%swrong argument:%s %s%s\n", Red(), Yellow(),
                                                      argv[count],
                                                      ColourReset());
        Help(argv[0], true);
        exit(1);
      } //}}}
      if (A[max_ids] == 0 || B[max_ids] == 0 || C[max_ids] == 0) {
        strcpy(ERROR_MSG, "<mode>: at least one index must be provided for \
each of A, B, and C");
        PrintError();
        Help(argv[0], true);
        exit(1);
      }
      angle *= PI / 180; // transform degrees to radians
      break; //}}}
    // error - incorrect <mode> //{{{
    case -1:
      strcpy(ERROR_MSG, "<mode> must be 'ref[lection]', 'inv[ersion]', \
or 'rot[ation]'");
      PrintError();
      Help(argv[0], true);
      exit(1); //}}}
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  // are joined coordinates provided?
  bool joined = BoolOption(argc, argv, "--joined");
  // use the last step?
  bool last = BoolOption(argc, argv, "--last");
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  int *Index_mol; // same as Index, but for molecules
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  VtfReadStruct(input_vsf, false, &Counts, &BeadType, &Bead, &Index,
                &MoleculeType, &Molecule, &Index_mol);
  InFile = calloc(Counts.BeadsTotal, sizeof *InFile); //}}}

  // '-n' option - specify molecule ids //{{{
  int n_opt_id[100] = {0}, n_opt_number = -1;
  if (MultiIntegerOption(argc, argv, "-n", &n_opt_number, n_opt_id)) {
    exit(1);
  }
  SortArray(n_opt_id, n_opt_number, 0);
  // test the provided indices are good
  for (int i = 0; i < n_opt_number; i++) {
    int id = n_opt_id[i];
    if (id > Counts.HighestResid || Index_mol[id] == -1) {
      strcpy(ERROR_MSG, "non-existent molecule index");
      PrintErrorOption("-n");
      ErrorPrintFile(input_vsf);
      fprintf(stderr, "%s, molecule %s%d%s\n", ErrRed(), ErrYellow(), id,
                                               ErrColourReset());
      exit(1);
    }
  } //}}}

  // save all beads
  for (int i = 0; i < Counts.BeadsTotal; i++) {
    Bead[i].Use = true;
  }

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // print initial stuff to output vcf file
  FILE *out = OpenFile(output_vcf, "w");
  PrintByline(out, argc, argv);
  fclose(out);

  // main loop //{{{
  int count_vcf = 0, // count steps in the vcf file
      file_line_count = 0; // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  while (true) {
    count_vcf++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (last) {
        fprintf(stdout, "\rDiscarding step: %d", count_vcf);
      } else {
        fprintf(stdout, "\rStep: %d", count_vcf);
        fflush(stdout);
      }
    } //}}}
    // decide whether this timestep is to be saved
    bool use = true; // TODO maybe later there'll be something
    // read timestep and calculate stuff, if the timestep should be used //{{{
    if (use) {
      if (!VtfReadTimestep(vcf, input_coor, &Box, &Counts, BeadType, &Bead,
                           Index, MoleculeType, Molecule,
                           &file_line_count, count_vcf)) {
        count_vcf--;
        break;
      }
      // remove pbc if provided molecules aren't joined
      if (!joined) {
        ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
        RemovePBCMolecules_new(Counts, Box, BeadType, &Bead,
                           MoleculeType, Molecule);
        FromFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      }

      for (int i = 0; i < Counts.Molecules; i++) {
        use = false;
        for (int j = 0; j < n_opt_number; j++) {
          if (Molecule[i].Index == n_opt_id[j]) {
            use = true;
            break;
          }
        }
        if (use) {
          int mtype = Molecule[i].Type;
          VECTOR point[3] = {{0, 0, 0},{0, 0, 0},{0, 0, 0}}, normal;
          switch (mode) {
            // mirror the molecule //{{{
            case 0:
              // normal to mirroring plane
              // TODO error when index[] is higher than i's number of beads
              // point A (on the normal) - point[0] //{{{
              for (int j = 0; j < A[max_ids]; j++) {
                int id = Molecule[i].Bead[A[j]];
                point[0].x += Bead[id].Position.x;
                point[0].y += Bead[id].Position.y;
                point[0].z += Bead[id].Position.z;
              }
              point[0].x /= A[max_ids];
              point[0].y /= A[max_ids];
              point[0].z /= A[max_ids]; //}}}
              // point B (on the normal) - point[1] //{{{
              for (int j = 0; j < B[max_ids]; j++) {
                int id = Molecule[i].Bead[B[j]];
                point[1].x += Bead[id].Position.x;
                point[1].y += Bead[id].Position.y;
                point[1].z += Bead[id].Position.z;
              }
              point[1].x /= B[max_ids];
              point[1].y /= B[max_ids];
              point[1].z /= B[max_ids]; //}}}
              // point C (on the plane) - point[2] //{{{
              for (int j = 0; j < C[max_ids]; j++) {
                int id = Molecule[i].Bead[C[j]];
                point[2].x += Bead[id].Position.x;
                point[2].y += Bead[id].Position.y;
                point[2].z += Bead[id].Position.z;
              }
              point[2].x /= C[max_ids];
              point[2].y /= C[max_ids];
              point[2].z /= C[max_ids]; //}}}
              // equation of the reflecting plane //{{{
              // normal vector
              normal.x = point[0].x - point[1].x;
              normal.y = point[0].y - point[1].y;
              normal.z = point[0].z - point[1].z;
              // d from ax+by+cz+d=0
              double d_plane = -(normal.x * point[2].x +
                                 normal.y * point[2].y +
                                 normal.z * point[2].z); //}}}
              // reflect the molecule
              for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
                int id = Molecule[i].Bead[j];
                // intersection of the line containing bead id with the plane
                double t = -(d_plane + normal.x * Bead[id].Position.x +
                                       normal.y * Bead[id].Position.y +
                                       normal.z * Bead[id].Position.z);
                t = t / (SQR(normal.x) + SQR(normal.y) + SQR(normal.z));
                VECTOR intersect;
                intersect.x = Bead[id].Position.x + normal.x * t;
                intersect.y = Bead[id].Position.y + normal.y * t;
                intersect.z = Bead[id].Position.z + normal.z * t;
                Bead[id].Position.x = 2 * intersect.x - Bead[id].Position.x;
                Bead[id].Position.y = 2 * intersect.y - Bead[id].Position.y;
                Bead[id].Position.z = 2 * intersect.z - Bead[id].Position.z;
              }
              break; //}}}
            // invert the molcule  //{{{
            case 1:
              // centre of inversion //{{{
              putchar('\n');
              for (int j = 0; j < A[max_ids]; j++) {
              // TODO error when index[] is higher than i's number of beads
                int id = Molecule[i].Bead[A[j]];
                point[0].x += Bead[id].Position.x;
                point[0].y += Bead[id].Position.y;
                point[0].z += Bead[id].Position.z;
              }
              point[0].x /= A[max_ids];
              point[0].y /= A[max_ids];
              point[0].z /= A[max_ids]; //}}}
              // invert the molecule
              for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
                int id = Molecule[i].Bead[j];
                Bead[id].Position.x += 2 * (point[0].x - Bead[id].Position.x);
                Bead[id].Position.y += 2 * (point[0].y - Bead[id].Position.y);
                Bead[id].Position.z += 2 * (point[0].z - Bead[id].Position.z);
              }
              break; //}}}
            // rotate the molcule  //{{{
            case 2:
              // rotation axis //{{{
              // first point
              for (int j = 0; j < A[max_ids]; j++) {
                int id = Molecule[i].Bead[A[j]];
                point[0].x += Bead[id].Position.x;
                point[0].y += Bead[id].Position.y;
                point[0].z += Bead[id].Position.z;
              }
              point[0].x /= A[max_ids];
              point[0].y /= A[max_ids];
              point[0].z /= A[max_ids];
              // secont point
              for (int j = 0; j < B[max_ids]; j++) {
                int id = Molecule[i].Bead[B[j]];
                point[1].x += Bead[id].Position.x;
                point[1].y += Bead[id].Position.y;
                point[1].z += Bead[id].Position.z;
              }
              point[1].x /= B[max_ids];
              point[1].y /= B[max_ids];
              point[1].z /= B[max_ids];
              if (N[max_ids] == 0) { // N not provided, use A and B as the axis
                normal.x = point[1].x - point[0].x;
                normal.y = point[1].y - point[0].y;
                normal.z = point[1].z - point[0].z;
              } else { // N provided, use normal to ABN plane as the axis
                // third point
                for (int j = 0; j < N[max_ids]; j++) {
                  int id = Molecule[i].Bead[N[j]];
                  point[2].x += Bead[id].Position.x;
                  point[2].y += Bead[id].Position.y;
                  point[2].z += Bead[id].Position.z;
                }
                point[2].x /= N[max_ids];
                point[2].y /= N[max_ids];
                point[2].z /= N[max_ids];
                VECTOR u, v;
                u.x = point[1].x - point[0].x;
                u.y = point[1].y - point[0].y;
                u.z = point[1].z - point[0].z;
                v.x = point[2].x - point[0].x;
                v.y = point[2].y - point[0].y;
                v.z = point[2].z - point[0].z;
                // normal to the ABN plane (u x v)
                normal.x = u.y * v.z - u.z * v.y;
                normal.y = u.z * v.x - u.x * v.z;
                normal.z = u.x * v.y - u.y * v.x;
              }
              // create unit rotation vector
              double dist = Length(normal);
              normal.x /= dist;
              normal.y /= dist;
              normal.z /= dist; //}}}
              // create rotation matrix //{{{
              struct Tensor {
                VECTOR x, y, z;
              } rot;
              double c = 1 - cos(angle);
              rot.x.x = cos(angle) + SQR(normal.x) * c;
              rot.x.y = normal.x * normal.y * c - normal.z * sin(angle);
              rot.x.z = normal.x * normal.z * c + normal.y * sin(angle);

              rot.y.x = normal.x * normal.y * c + normal.z * sin(angle);
              rot.y.y = cos(angle) + SQR(normal.y) * c;
              rot.y.z = normal.y * normal.z * c - normal.x * sin(angle);

              rot.z.x = normal.x * normal.z * c - normal.y * sin(angle);
              rot.z.y = normal.y * normal.z * c + normal.x * sin(angle);
              rot.z.z = cos(angle) + SQR(normal.z) * c; //}}}
              // rotate the molecule
              // i) move the beginning to C
              point[0].x = 0;
              point[0].y = 0;
              point[0].z = 0;
              for (int j = 0; j < C[max_ids]; j++) {
                int id = Molecule[i].Bead[C[j]];
                point[0].x += Bead[id].Position.x;
                point[0].y += Bead[id].Position.y;
                point[0].z += Bead[id].Position.z;
              }
              point[0].x /= C[max_ids];
              point[0].y /= C[max_ids];
              point[0].z /= C[max_ids];
              for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
                int id = Molecule[i].Bead[j];
                Bead[id].Position.x -= point[0].x;
                Bead[id].Position.y -= point[0].y;
                Bead[id].Position.z -= point[0].z;
              }
              // ii) actual rotation
              for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
                int id = Molecule[i].Bead[j];
                VECTOR new;
                new.x = rot.x.x * Bead[id].Position.x +
                        rot.x.y * Bead[id].Position.y +
                        rot.x.z * Bead[id].Position.z;
                new.y = rot.y.x * Bead[id].Position.x +
                        rot.y.y * Bead[id].Position.y +
                        rot.y.z * Bead[id].Position.z;
                new.z = rot.z.x * Bead[id].Position.x +
                        rot.z.y * Bead[id].Position.y +
                        rot.z.z * Bead[id].Position.z;
                Bead[id].Position = new;
              }
              // iii) move the beginning back to where it was
              for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
                int id = Molecule[i].Bead[j];
                Bead[id].Position.x += point[0].x;
                Bead[id].Position.y += point[0].y;
                Bead[id].Position.z += point[0].z;
              }
              break; //}}}
          }
        }
      }

      // write to output .vcf file
      out = OpenFile(output_vcf, "a");
      VtfWriteCoorIndexed(out, stuff, Counts, Bead, Box);
      fclose(out);
      //}}}
    // skip the timestep, if it shouldn't be saved //{{{
    } else {
      if (!VtfSkipTimestep(vcf, input_coor, &file_line_count, count_vcf)) {
        count_vcf--;
        break;
      }
    } //}}}
    // definitely break the loop - later, there'll be option for specific step
    break;
  }
  fclose(vcf); //}}}

  // print last step count? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // free memory
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(Index_mol);
  free(stuff);
  free(InFile);

  return 0;
}
