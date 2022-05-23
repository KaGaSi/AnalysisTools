#include "../AnalysisTools.h"
// TODO: split into two utilities - adding existing (vtf) configuration and
//       generating addition from FIELD?
// TODO: --random for -vtf option (i.e., place added system's components
//       randomly in a new box)

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
AddToSystem either creates a system from scratch or adds unbonded beads \
and/or molecules to an existing system. The new components are defined either \
by a FIELD-like file (an input file for DL_MESO simulation program) or by \
vsf/vcf files (-vtf option). In the first case, new components are placed \
randomly (with several possible constraints), while in the second case, the \
provided coordinates are used as is.\n\n");
  }
  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.vcf> <out.vsf> <out.vcf> [options]\n\n", cmd);

  fprintf(ptr, "      <input>/--     input coordinate file (vcf or vtf format) \
or '--' to create a system from scratch\n");
  fprintf(ptr, "      <out.vsf>      output structure file (vsf format)\n");
  fprintf(ptr, "      <out.vcf>      output coordinate file (vcf format)\n");
  fprintf(ptr, "   [general options]\n");
  fprintf(ptr, "      -st <int>            timestep to use (default: 1)\n");
  fprintf(ptr, "      -xyz <name>          save coordinates to an xyz too\n");
  fprintf(ptr, "      -xb <bead name(s)>   replace original beads instead of \
increasing the total number of beads\n");
  fprintf(ptr, "      -b <x> <y> <z>       size of the new cuboid box\n");
  fprintf(ptr, "      -ba <a> <b> <g>      angles of the new triclinic box\n");
  fprintf(ptr, "      -f <name>            FIELD-like file with molecules \
to add (default: FIELD)\n");
  fprintf(ptr, "      -ld <float>          specify lowest distance from \
chosen bead types (default: none)\n");
  fprintf(ptr, "      -hd <float>          specify highest distance from \
chosen bead types (default: none)\n");
  fprintf(ptr, "      -bt <name(s)>        specify bead types new beads \
should be far from/near to (default: none)\n");
  fprintf(ptr, "      -cx <num> <num2>     constrain x coordinate of \
randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -cy <num> <num2>     constrain y coordinate of \
randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -cz <num> <num2>     constrain z coordinate of \
randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -gc                  use molecule's geometric centre \
for the distance check instead of its first bead\n");
  fprintf(ptr, "      -sd <int>            seed for the random number \
generator (default: clock-based seed)\n");
  fprintf(ptr, "      --no-rotate          do not randomly rotate added \
molecules\n");
  fprintf(ptr, "      -vtf <vsf> <vcf>     use vtf files instead of \
FIELD (divided to vsf and vcf files)\n");
  fprintf(ptr, "      -offset <x> <y> <z>  offset the added system \
by given amount in all directions\n");
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
  while ((count+1) < argc &&
         (argv[count+1][0] != '-' || strcmp(argv[count+1], "--") == 0)) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = req_args; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-vtf") != 0 &&
        strcmp(argv[i], "-offset") != 0 &&
        (argv[i][1] < '0' || argv[i][1] > '9') &&
        strcmp(argv[i], "-sd") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-xyz") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-ld") != 0 &&
        strcmp(argv[i], "-hd") != 0 &&
        strcmp(argv[i], "-cx") != 0 &&
        strcmp(argv[i], "-cy") != 0 &&
        strcmp(argv[i], "-cz") != 0 &&
        strcmp(argv[i], "-gc") != 0 &&
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-ba") != 0 &&
        strcmp(argv[i], "--no-rotate") != 0 &&
        strcmp(argv[i], "-xb") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE] = "", // unchanged => new system
       input_vsf[LINE] = "";
  bool vtf = false;
  if (strcmp(argv[++count], "--") != 0) {
    // test that <input> filename ends with '.vcf' or '.vtf'
    snprintf(input_coor, LINE, "%s", argv[count]);
    if (!InputCoor(&vtf, input_coor, input_vsf)) {
      Help(argv[0], true);
      exit(1);
    }
  } else {
    strcpy(input_vsf, "in.vsf"); // won't be used
  } //}}}

  // <out.vsf> - output vsf file //{{{
  char output_vsf[LINE] = "";
  snprintf(output_vsf, LINE, "%s", argv[++count]);
  // test that <out.vsf> filename ends with '.vsf'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output vcf file //{{{
  char output_vcf[LINE] = "";
  snprintf(output_vcf, LINE, "%s", argv[++count]);
  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, input_vsf, LINE, &verbose, &silent, &detailed);

  // -f <add> - FIELD-like file with molecules to add //{{{
  char input_add[LINE] = "";
  if (FileOption(argc, argv, "-f", input_add, LINE)) {
    exit(1);
  }
  if (input_add[0] == '\0') {
    strcpy(input_add, "FIELD");
  } //}}}

  // -vtf <vsf> <vcf> - vtf file(s) to use instead of FIELD //{{{
  char add_vsf[LINE] = "", input_coor_add[LINE] = "";
  // 1) vsf file
  if (FileOption(argc, argv, "-vtf", add_vsf, LINE)) {
    exit(1);
  }
  // 2) if vsf file exists, look for vcf
//bool vtf_add = true; // if -vtf is present present, is the file a vtf format?
  if (strlen(add_vsf) > 0) {
    ext = 2;
    strcpy(extension[0], ".vsf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(add_vsf, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
    if (add_vsf[strlen(add_vsf)-2] == 't') { // if *.vtf file, use it as vcf too
      snprintf(input_coor_add, LINE, "%s", add_vsf);
    } else { // if *.vsf file, read the vcf file
      for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-vtf") == 0) {
          if (argc > (i+2)) { // enough arguments for a vcf file?
            char temp[LINE];
            // save vsf filename
            snprintf(temp, LINE, "%s", argv[i+1]);
            // copy vcf filename to (i+1)th place - required by FileOption()
            snprintf(argv[i+1], LINE, "%s", argv[i+2]);
            // read vcf file name
            if (FileOption(argc, argv, "-vtf", input_coor_add, LINE)) {
              exit(1);
            }
            // restore vsf filename so the command in unchanged
            strcpy(argv[i+1], temp);
            // coordinate file must be vcf, because there's already vsf
            ext = 1;
            strcpy(extension[0], ".vcf");
            ext = ErrorExtension(input_coor_add, ext, extension);
            if (ext == -1) {
              Help(argv[0], true);
              exit(1);
//          } else {
//            vtf_add = false;
            }
            break;
          } else { // missing vcf file name
            ErrorPrintError_old();
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "-vtf");
            ColourChange(STDERR_FILENO, RED);
            fprintf(stderr, " - missing second file (vcf format;");
            fprintf(stderr, " cannot be full vtf, because vsf is used)\n\n");
            ColourReset(STDERR_FILENO);
            exit(1);
          }
        }
      }
    }
  } //}}}

  // -offset <x> <y> <z> define offset for -vtf file //{{{
  double offset[100] = {1000000};
  if (MultiDoubleOption(argc, argv, "-offset", &count, offset)) {
    exit(1);
  }
  if (count != 3) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-offset");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - three numbers required\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // Warning - missing -vtf option
  if (strlen(add_vsf) == 0 && offset[0] != 1000000) {
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "\nWarning: ");
    ColourChange(STDERR_FILENO, CYAN);
    fprintf(stderr, "-offset");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, " option has no effect if -vtf is not used\n");
    ColourReset(STDERR_FILENO);
  } //}}}

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // save into xyz file? //{{{
  char output_xyz[LINE] = "";
  if (FileOption(argc, argv, "-xyz", output_xyz, LINE)) {
    exit(1);
  } //}}}

  // lowest and/or highest distance from beads of type specified by '-bt' //{{{
  double lowest_dist = -1, highest_dist = -1;
  if (DoubleOption(argc, argv, "-ld", &lowest_dist)) {
    exit(1);
  }
  if (DoubleOption(argc, argv, "-hd", &highest_dist)) {
    exit(1);
  }
  // errors: 1) if a new system is generated, it cannot be used
  //         2) if '-ld' and/or '-hd' are present, '-bt' must be too
  if (highest_dist != -1 || lowest_dist != -1) {
    // 1)
    if (strcmp(argv[1],"--") == 0) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, "if new system is generated,");
      fprintf(stderr, "there cannot be -ld/-hd/-bt options present\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    }
    // 2)
    bool bt = false;
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-bt") == 0) {
        bt = true;
      }
    }
    if (!bt) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, "if '-ld' and/or '-hd' is used,");
      fprintf(stderr, "'-bt' must be specified as well\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    }
  } //}}}

  // coordinate constraints //{{{
  // x direction //{{{
  int test = 2;
  double range[2] = {0, 0};
  if (MultiDoubleOption(argc, argv, "-cx", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-cx");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - two non-negative numbers required");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  VECTOR constraint[2];
  constraint[0].x = range[0];
  constraint[1].x = range[1]; //}}}
  // y direction //{{{
  test = 2;
  range[0] = range[1] = 0;
  if (MultiDoubleOption(argc, argv, "-cy", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-cy");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - two non-negative numbers required");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  constraint[0].y = range[0];
  constraint[1].y = range[1]; //}}}
  // z direction //{{{
  test = 2;
  range[0] = range[1] = 0;
  if (MultiDoubleOption(argc, argv, "-cz", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-cz");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - two non-negative numbers required");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  constraint[0].z = range[0];
  constraint[1].z = range[1]; //}}}
  //}}}

  // use centre of mass for distance check of new molecules
  bool com = BoolOption(argc, argv, "-gc");

  // rotate added molecules?
  bool no_rot = BoolOption(argc, argv, "--no-rotate");

  // define new box size //{{{
  double box_option[100] = {-1};
  if (MultiDoubleOption(argc, argv, "-b", &count, box_option)) {
    exit(1);
  }
  if (count != 3) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-b");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - three non-negative numbers required\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  double box_angle_option[100] = {-1};
  if (MultiDoubleOption(argc, argv, "-ba", &count, box_angle_option)) {
    exit(1);
  }
  if (count != 3) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-ba");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - three non-negative numbers required\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  //}}}

  // seed for the random number generator //{{{
  int seed = -1; // not present
  if (IntegerOption(argc, argv, "-sd", &seed)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from input vtf file(s) if present //{{{
  SYSTEM S_orig;
  if (strlen(input_coor) > 0) { // is there an input coordinate file?
    S_orig = VtfReadStruct(input_vsf, detailed);
  } //}}}

  // -xb <name(s)> - specify what bead types to exchange //{{{
  bool sw = BoolOption(argc, argv, "-xb"); // is -xb present?
  // error - if -xb is used, 
  if (sw && strlen(input_coor) == 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-xb");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - <input> file must be present\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // which beads to exchange?
  if (BeadTypeOption(argc, argv, "-xb", false, &S_orig)) {
    exit(0);
  }
  // use Write flag to decide which bead types to use
  bool all_false = true; // no '-xb' option
  for (int i = 0; i < S_orig.Count.nBeadTypes; i++) {
    S_orig.BeadType[i].Write = S_orig.BeadType[i].Use;
    S_orig.BeadType[i].Use = false; // this flag may be used later
    if (S_orig.BeadType[i].Write) {
      all_false = false; // '-xb' option is present
    }
  }
  if (all_false) {
    for (int i = 0; i < S_orig.Count.nBeadTypes; i++) {
      if (S_orig.BeadType[i].Charge == 0) {
        S_orig.BeadType[i].Write = true;
      }
    }
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", false, &S_orig)) {
    exit(0);
  } //}}}

  // seed random number generator //{{{
  if (seed > -1) {
    srand(seed);
  } else {
    srand(time(0));
  } //}}}

  // print original system (if present) //{{{
  if (verbose && strlen(input_coor) > 0) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(S_orig);
    if (start > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", start);
    }
  } //}}}

  // array for the timestep preamble
  char *stuff = calloc(LINE, sizeof *stuff);

  // open input coordinate file //{{{
  FILE *vcf;
  if (strlen(input_coor) > 0) {
    vcf = OpenFile(input_coor, "r");
    count= 0;
    if (!silent) {
      fprintf(stdout, "Using step %6d\n", ++count);
    }
    int file_line_count = 0, count_vcf = 0;
    VtfReadTimestep(vcf, input_coor, &S_orig, &file_line_count,
                    count_vcf, stuff);
    fclose(vcf);
//  PrintCounts(Counts_orig);
  } //}}}

  SYSTEM S_add;

  // TODO FIELD must be completely redone
  if (strlen(add_vsf) == 0) { // read stuff to be added from FIELD //{{{
//  ReadField(input_add, '\0', &Counts_add, &S_add.BeadType, &S_add.Bead,
//            &Index_add, &S_add.MoleculeType, &S_add.Molecule,
//            &bond_type, &angle_type, &dihedral_type);
//  S_add.Box.Length = S_orig.Box.Length; //}}}
  } else { // read stuff to add from vtf file(s) ('-vtf' option) //{{{
    S_add = VtfReadStruct(add_vsf, false);
    // read coordinates
    vcf = OpenFile(input_coor_add, "r");
    int file_line_count = 0, count_vcf = 0;
    VtfReadTimestep(vcf, input_coor_add, &S_add, &file_line_count, count_vcf, stuff);
    fclose(vcf);
    // TODO: !no_rot? ...shouldn't -vtf be this by default?
    VECTOR rotated[S_add.Count.nBeadsCoor];
    if (!no_rot) { //{{{
      // random rotation axis
      VECTOR random = {0};
      random.x = (double)rand() / ((double)RAND_MAX) * 2 - 1; // a number <-1,1>
      random.y = (double)rand() / ((double)RAND_MAX) * 2 - 1;
      random.z = (double)rand() / ((double)RAND_MAX) * 2 - 1;
      double dist = Length(random);
      random.x /= dist;
      random.y /= dist;
      random.z /= dist;
      // random rotation angle
      double angle = (double)rand() / ((double)RAND_MAX) * PI;
      // create rotation matrix
      struct Tensor {
        VECTOR x, y, z;
      } rot;
      rot.x.x = cos(angle) + SQR(random.x) * (1 - cos(angle));
      rot.x.y = random.x * random.y * (1 - cos(angle)) - random.z * sin(angle);
      rot.x.z = random.x * random.z * (1 - cos(angle)) + random.y * sin(angle);

      rot.y.x = random.x * random.y * (1 - cos(angle)) + random.z * sin(angle);
      rot.y.y = cos(angle) + SQR(random.y) * (1 - cos(angle));
      rot.y.z = random.y * random.z * (1 - cos(angle)) - random.x * sin(angle);

      rot.z.x = random.x * random.z * (1 - cos(angle)) - random.y * sin(angle);
      rot.z.y = random.y * random.z * (1 - cos(angle)) + random.x * sin(angle);
      rot.z.z = cos(angle) + SQR(random.z) * (1 - cos(angle));
      printf("xxXxx\n");
      printf("%lf %lf %lf\n", random.x, random.y, random.z);
      printf("%lf\n", (double)rand());
      printf("XxXxX\n");
      printf("%lf %lf %lf\n", rot.x.x, rot.x.y, rot.x.z);
      printf("%lf %lf %lf\n", rot.y.x, rot.y.y, rot.y.z);
      printf("%lf %lf %lf\n", rot.z.x, rot.z.y, rot.z.z);
      // transform the prototype molecule (rotation matrix * coordinates)
      for (int i = 0; i < S_add.Count.nBeadsCoor; i++) {
        rotated[i].x = rot.x.x * (S_add.Bead[i].Position.x - S_add.Box.Length.x / 2)
                     + rot.x.y * (S_add.Bead[i].Position.y - S_add.Box.Length.y / 2)
                     + rot.x.z * (S_add.Bead[i].Position.z - S_add.Box.Length.z / 2);
        rotated[i].y = rot.y.x * (S_add.Bead[i].Position.x - S_add.Box.Length.x / 2)
                     + rot.y.y * (S_add.Bead[i].Position.y - S_add.Box.Length.y / 2)
                     + rot.y.z * (S_add.Bead[i].Position.z - S_add.Box.Length.z / 2);
        rotated[i].z = rot.z.x * (S_add.Bead[i].Position.x - S_add.Box.Length.x / 2)
                     + rot.z.y * (S_add.Bead[i].Position.y - S_add.Box.Length.y / 2)
                     + rot.z.z * (S_add.Bead[i].Position.z - S_add.Box.Length.z / 2);
      }
      for (int i = 0; i < S_add.Count.nBeadsCoor; i++) {
        S_add.Bead[i].Position.x = rotated[i].x + offset[0] + S_add.Box.Length.x / 2;
        S_add.Bead[i].Position.y = rotated[i].y + offset[1] + S_add.Box.Length.y / 2;
        S_add.Bead[i].Position.z = rotated[i].z + offset[2] + S_add.Box.Length.z / 2;
      }
     //}}}
    } else { // don't rotate //{{{
      ToFractionalCoor(S_add.Count.nBeadsTotal, &S_add.Bead, S_add.Box);
      for (int i = 0; i < S_add.Count.nBeadsCoor; i++) {
        int id = S_add.BeadsCoor[i];
        S_add.Bead[id].Position.x += offset[0];
        S_add.Bead[id].Position.y += offset[1];
        S_add.Bead[id].Position.z += offset[2];
      }
      FromFractionalCoor(S_add.Count.nBeadsTotal, &S_add.Bead, S_add.Box);
    } //}}}
    // allocate memory only to free it later
//  bond_type = calloc(1, sizeof (PARAMS));
//  angle_type = calloc(1, sizeof (PARAMS));
//  dihedral_type = calloc(1, sizeof (PARAMS));
  } //}}}

  SYSTEM S_new;

  // set final box size //{{{
  /*
   * i -b isn't used, set box size as the larger of the dimensions from
   * original and to-be-added systems
   */
  S_new.Box = InitBox;
  S_new.Box.Length.x = 0;
  S_new.Box.Length.y = 0;
  S_new.Box.Length.z = 0;
  if (box_option[0] == -1) {
    if (S_add.Box.Length.x > S_orig.Box.Length.x) {
      S_new.Box.Length.x = S_add.Box.Length.x;
    } else {
      S_new.Box.Length.x = S_orig.Box.Length.x;
    }
    if (S_add.Box.Length.y > S_orig.Box.Length.y) {
      S_new.Box.Length.y = S_add.Box.Length.y;
    } else {
      S_new.Box.Length.y = S_orig.Box.Length.y;
    }
    if (S_add.Box.Length.z > S_orig.Box.Length.z) {
      S_new.Box.Length.z = S_add.Box.Length.z;
    } else {
      S_new.Box.Length.z = S_orig.Box.Length.z;
    }
  } else {
    S_new.Box = S_orig.Box; // TODO just for now
    S_new.Box.Length.x = box_option[0];
    S_new.Box.Length.y = box_option[1];
    S_new.Box.Length.z = box_option[2];
  }
  if (box_angle_option[0] != -1) {
    S_new.Box.alpha = box_angle_option[0];
    S_new.Box.beta = box_angle_option[1];
    S_new.Box.gamma = box_angle_option[2];
  }
  //}}}

  // define 'box' for additions using constraints (-c{x,y,z} options)//{{{
  BOX constraint_box = InitBox;
  if (constraint[1].x != 0) {
    constraint_box.Length.x = constraint[1].x - constraint[0].x;
  } else {
    constraint_box.Length.x = S_new.Box.Length.x;
  }
  if (constraint[1].y != 0) {
    constraint_box.Length.y = constraint[1].y - constraint[0].y;
  } else {
    constraint_box.Length.y = S_new.Box.Length.y;
  }
  if (constraint[1].z != 0) {
    constraint_box.Length.z = constraint[1].z - constraint[0].z;
  } else {
    constraint_box.Length.z = S_new.Box.Length.z;
  } //}}}

  // error - no box size (should never trigger) //{{{
  if (S_new.Box.Length.x == 0 ||
      S_new.Box.Length.y == 0 ||
      S_new.Box.Length.z == 0) {
    strcpy(ERROR_MSG, "zero box size for the new system");
    PrintError();
    exit(1);
  } //}}}

  // check number of exchangeable beads //{{{
  int can_be_exchanged = 0;
  for (int i = 0; i < S_orig.Count.nBeadsTotal; i++) {
    int btype = S_orig.Bead[i].Type;
    if (S_orig.Bead[i].Molecule == -1 && S_orig.BeadType[btype].Write) {
      can_be_exchanged++;
    }
  }
  // count beads to be added
  if (sw && S_add.Count.nBeadsCoor > can_be_exchanged) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, "insufficient beads to exchange for new ones\n");
    fprintf(stderr, "     Exchangeable beads in the original system: ");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d\n", can_be_exchanged);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, "     Beads to be added: ");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d\n\n", S_add.Count.nBeadsCoor);
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}

  // if '-gc' is used, put prototypes' geometric centres to (0,0,0) //{{{
  if (com) {
    for (int i = 0; i < S_add.Count.nMolecules; i++) {
      int mtype = S_add.Molecule[i].Type;
      VECTOR geom_centre;
      geom_centre.x = 0;
      geom_centre.y = 0;
      geom_centre.z = 0;
      for (int j = 0; j < S_add.MoleculeType[mtype].nBeads; j++) {
        int id = S_add.Molecule[i].Bead[j];
        geom_centre.x += S_add.Bead[id].Position.x;
        geom_centre.y += S_add.Bead[id].Position.y;
        geom_centre.z += S_add.Bead[id].Position.z;
      }
      geom_centre.x /= S_add.MoleculeType[mtype].nBeads;
      geom_centre.y /= S_add.MoleculeType[mtype].nBeads;
      geom_centre.z /= S_add.MoleculeType[mtype].nBeads;
      for (int j = 0; j < S_add.MoleculeType[mtype].nBeads; j++) {
        int id = S_add.Molecule[i].Bead[j];
        S_add.Bead[id].Position.x -= geom_centre.x;
        S_add.Bead[id].Position.y -= geom_centre.y;
        S_add.Bead[id].Position.z -= geom_centre.z;
      }
    }
  } //}}}

  // print what is to be added //{{{
  if (verbose) {
    fprintf(stdout, "\nBEADS AND MOLECULES TO ADD\n");
    VerboseOutput(S_add);
  } //}}}

  /* decide which beads to exchange //{{{
   * i.e., give them Bead[].Flag = true); has effect only if --switch is used
   */
  // zeroize Bead[].Flag //{{{
  for (int i = 0; i < S_orig.Count.nBeadsCoor; i++) {
    S_orig.Bead[i].Use = false;
  } //}}}
  count = 0; // counts bead in the original Bead[] struct
  for (int i = 0; i < S_add.Count.nBeadsCoor; i++) {
    for (; count < S_orig.Count.nBeadsCoor; count++) {
      int type = S_orig.Bead[count].Type;
      if (S_orig.BeadType[type].Write && S_orig.Bead[count].Molecule == -1) {
        S_orig.Bead[count].Use = true; // exchange bead 'count'
        break;
      }
    }
    count++; // loop didn't update count because of the break
  } //}}}

  // join original and added systems (depending on '--switch' mode)
  if (sw) { // switch old beads for new ones? //{{{
    S_new.Count.nBeadsCoor = S_orig.Count.nBeadsCoor;
    S_new.Count.nBeadsTotal = S_orig.Count.nBeadsTotal;
    S_new.Count.nBonded = S_orig.Count.nBonded + S_add.Count.nBonded;
    S_new.Count.nUnbonded = S_orig.Count.nBeadsCoor - S_new.Count.nBonded;
    S_new.Count.nBondTypes = S_add.Count.nBondTypes;
    S_new.Count.nAngleTypes = S_add.Count.nAngleTypes;
    S_new.Count.nMolecules = S_orig.Count.nMolecules + S_add.Count.nMolecules;
    // fill BeadType struct for the new system
    S_new.Count.nBeadTypes = S_orig.Count.nBeadTypes;
    // 1) copy original BeadType
    CopyBeadType(S_new.Count.nBeadTypes, &S_new.BeadType, S_orig.BeadType, 3);
    // 2) add new bead types - the check is based only on Name //{{{
    for (int i = 0; i < S_add.Count.nBeadTypes; i++) {
      bool new = true;
      for (int j = 0; j < S_orig.Count.nBeadTypes; j++) {
        if (strcmp(S_add.BeadType[i].Name, S_orig.BeadType[j].Name) == 0) {
          new = false;
          S_new.BeadType[j].Number += S_add.BeadType[i].Number;
          break;
        }
      }
      if (new) {
        int type = S_new.Count.nBeadTypes;
        S_new.BeadType = realloc(S_new.BeadType, sizeof (BEADTYPE) * (type + 1));
        S_new.BeadType[type] = S_add.BeadType[i];
        S_new.Count.nBeadTypes++;
      }
    } //}}}
    // fill MoleculeType struct for the new system
    S_new.Count.nMoleculeTypes = S_orig.Count.nMoleculeTypes;
    S_new.MoleculeType = realloc(S_new.MoleculeType, sizeof (MOLECULETYPE) * S_new.Count.nMoleculeTypes);
    // copy original MoleculeType to _new //{{{
    for (int i = 0; i < S_new.Count.nMoleculeTypes; i++) {
      S_new.MoleculeType[i] = S_orig.MoleculeType[i];
      S_new.MoleculeType[i].Bead = malloc(sizeof *S_new.MoleculeType[i].Bead * S_new.MoleculeType[i].nBeads);
      for (int j = 0; j < S_new.MoleculeType[i].nBeads; j++) {
        S_new.MoleculeType[i].Bead[j] = S_orig.MoleculeType[i].Bead[j];
      }
      S_new.MoleculeType[i].Bond = malloc(sizeof *S_new.MoleculeType[i].Bond * S_new.MoleculeType[i].nBonds);
      for (int j = 0; j < S_new.MoleculeType[i].nBonds; j++) {
        S_new.MoleculeType[i].Bond[j][0] = S_orig.MoleculeType[i].Bond[j][0];
        S_new.MoleculeType[i].Bond[j][1] = S_orig.MoleculeType[i].Bond[j][1];
        S_new.MoleculeType[i].Bond[j][2] = S_orig.MoleculeType[i].Bond[j][2];
      }
    } //}}}
    // add new molecule types - check if their the same based only on Name //{{{
    for (int i = 0; i < S_add.Count.nMoleculeTypes; i++) {
      bool new = true;
      for (int j = 0; j < S_orig.Count.nMoleculeTypes; j++) {
        if (strcmp(S_add.MoleculeType[i].Name, S_orig.MoleculeType[j].Name) == 0) {
          new = false;
          S_new.MoleculeType[j].Number += S_add.MoleculeType[i].Number;
          break;
        }
      }
      if (new) {
        int type = S_new.Count.nMoleculeTypes;
        S_new.MoleculeType = realloc(S_new.MoleculeType, sizeof (MOLECULETYPE) * (type + 1));
        S_new.MoleculeType[type] = S_add.MoleculeType[i];
        S_new.MoleculeType[type].Bead = malloc(sizeof *S_new.MoleculeType[type].Bead *
                                   S_new.MoleculeType[type].nBeads);
        for (int j = 0; j < S_new.MoleculeType[type].nBeads; j++) {
          int old_type = S_add.MoleculeType[i].Bead[j];
          int btype = FindBeadType(S_add.BeadType[old_type].Name, S_new);
          S_new.MoleculeType[type].Bead[j] = btype;
        }
        S_new.MoleculeType[type].Bond = malloc(sizeof *S_new.MoleculeType[i].Bond *
                                   S_new.MoleculeType[type].nBonds);
        for (int j = 0; j < S_new.MoleculeType[type].nBonds; j++) {
          S_new.MoleculeType[type].Bond[j][0] = S_add.MoleculeType[i].Bond[j][0];
          S_new.MoleculeType[type].Bond[j][1] = S_add.MoleculeType[i].Bond[j][1];
          S_new.MoleculeType[type].Bond[j][2] = S_add.MoleculeType[i].Bond[j][2];
        }
        if (S_new.MoleculeType[i].nAngles > 0) {
          S_new.MoleculeType[type].Angle = malloc(sizeof *S_new.MoleculeType[type].Angle *
                                      S_new.MoleculeType[type].nAngles);
          for (int j = 0; j < S_new.MoleculeType[type].nAngles; j++) {
            S_new.MoleculeType[type].Angle[j][0] = S_add.MoleculeType[i].Angle[j][0];
            S_new.MoleculeType[type].Angle[j][1] = S_add.MoleculeType[i].Angle[j][1];
            S_new.MoleculeType[type].Angle[j][2] = S_add.MoleculeType[i].Angle[j][2];
            S_new.MoleculeType[type].Angle[j][3] = S_add.MoleculeType[i].Angle[j][3];
          }
        }
        if (S_new.MoleculeType[i].nDihedrals > 0) {
          S_new.MoleculeType[type].Dihedral = malloc(sizeof *S_new.MoleculeType[type].Dihedral *
                                         S_new.MoleculeType[type].nDihedrals);
          for (int j = 0; j < S_new.MoleculeType[type].nDihedrals; j++) {
            S_new.MoleculeType[type].Dihedral[j][0] = S_add.MoleculeType[i].Dihedral[j][0];
            S_new.MoleculeType[type].Dihedral[j][1] = S_add.MoleculeType[i].Dihedral[j][1];
            S_new.MoleculeType[type].Dihedral[j][2] = S_add.MoleculeType[i].Dihedral[j][2];
            S_new.MoleculeType[type].Dihedral[j][3] = S_add.MoleculeType[i].Dihedral[j][3];
            S_new.MoleculeType[type].Dihedral[j][4] = S_add.MoleculeType[i].Dihedral[j][4];
          }
        }
        S_new.Count.nMoleculeTypes++;
      }
    } //}}}
    // fill Bead struct for the new system
    S_new.Bead = realloc(S_new.Bead, sizeof (BEAD) * S_new.Count.nBeadsCoor);
    // copy unbonded beads not to be exchanged to the start of S_new.Bead //{{{
    // TODO: assumes unbonded beads are before bonded beads
    count = 0; // counts copied beads
    for (int i = 0; i < S_orig.Count.nUnbonded; i++) {
      // first, copy only beads of the type that's not to be exchange
      if (!S_orig.Bead[i].Use) {
        S_new.Bead[count] = S_orig.Bead[i];
        S_new.Bead[count].Molecule = -1;
        S_new.Bead[count].Use = false; // do not rewrite, obviously
        count++;
      }
    }
    // count ended at <number of unbonded original beads> - <added beads> //}}}
    // put unbonded beads to be added beyond the unchanged unbonded beads //{{{
    count = S_orig.Count.nUnbonded - S_add.Count.nBeadsCoor; // just to be sure
    for (int i = 0; i < S_add.Count.nUnbonded; i++) {
      int type = S_add.Bead[i].Type;
      int new_type = FindBeadType(S_add.BeadType[type].Name, S_new);
      S_new.Bead[count] = S_orig.Bead[i];
      S_new.Bead[count].Type = new_type;
      S_new.Bead[count].Molecule = -1;
      S_new.Bead[count].Use = true; // coordinates to be rewritten
      count++; // use count to make it consistent & easy to read
    } //}}}
    // copy the original bonded beads //{{{
    count = S_new.Count.nUnbonded;
    for (int i = S_orig.Count.nUnbonded; i < S_orig.Count.nBeadsCoor; i++) {
      S_new.Bead[count] = S_orig.Bead[i];
      S_new.Bead[count].Use = false; // coordinates to be rewritten
      count++;
    } //}}}
    // put bonded beads to be added at the very end //{{{
    count = S_new.Count.nUnbonded + S_orig.Count.nBonded;
    for (int i = S_add.Count.nUnbonded; i < S_add.Count.nBeadsCoor; i++) {
      int type = S_add.Bead[i].Type;
      int new_type = FindBeadType(S_add.BeadType[type].Name, S_new);
      S_new.Bead[count] = S_add.Bead[i];
      S_new.Bead[count].Type = new_type;
      S_new.Bead[count].Molecule = S_add.Bead[i].Molecule + S_orig.Count.nMolecules;
      S_new.Bead[count].Use = true; // coordinates to be rewritten
      count++; // use count to make it consistent & easy to read
    } //}}}
    // alocate new molecule struct
    S_new.Molecule = realloc(S_new.Molecule, sizeof (MOLECULE) * S_new.Count.nMolecules);
    // copy original molecules to _new struct //{{{
    for (int i = 0; i < S_orig.Count.nMolecules; i++) {
      int type = S_orig.Molecule[i].Type;
      S_new.Molecule[i].Type = type;
      S_new.Molecule[i].Bead = malloc(sizeof *S_new.Molecule[i].Bead * S_new.MoleculeType[type].nBeads);
      for (int j = 0; j < S_new.MoleculeType[type].nBeads; j++) {
        S_new.Molecule[i].Bead[j] = S_orig.Molecule[i].Bead[j] - S_add.Count.nBonded;
      }
    } //}}}
    // put _add molecules into _new struct //{{{
    count = S_new.Count.nBeadsCoor - S_add.Count.nBonded;
    for (int i = 0; i < S_add.Count.nMolecules; i++) {
      int add_type = S_add.Molecule[i].Type;
      int new_type = FindMoleculeType(S_add.MoleculeType[add_type].Name, S_new);
      int new_i = S_orig.Count.nMolecules + i;
      S_new.Molecule[new_i].Type = new_type;
      S_new.Molecule[new_i].Bead = malloc(sizeof *S_new.Molecule[new_i].Bead *
                                   S_new.MoleculeType[new_type].nBeads);
      for (int j = 0; j < S_new.MoleculeType[new_type].nBeads; j++) {
        S_new.Molecule[new_i].Bead[j] = count;
        count++;
      }
    } //}}}
    //}}}
  } else { // or add beads to the system? //{{{
    BOX Box_new = S_new.Box;
    S_new = CopySystem(S_orig);
    ConcatenateSystems(&S_new, S_add, Box_new);
//  FreeSystem(&S_orig);
//  PruneSystem(&S_new);
  } //}}}
//FillMolBTypes(S_new.TypesOfMolecules, &S_new.MoleculeType);
//FillMolMassCharge(S_new.nMoleculeTypes, &S_new.MoleculeType, S_new.BeadType);

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(S_new);
//  PrintBondTypes2(S_new.TypesOfBonds, bond_type);
  } //}}}

  // add beads randomly if FIELD-like file is used //{{{
  double dist;
  if (strlen(add_vsf) == 0) {
    count = 0;
    // add monomeric beads //{{{
    for (int i = 0; i < S_add.Count.nUnbonded; i++) {
      VECTOR random;
      if (lowest_dist != -1 || highest_dist != -1) {
        double min_dist;
        int tries = 0;
        do {
          tries++;
          if (tries == 1000000) {
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "\nWarning: million attempts");
            fprintf(stderr, " to place a bead failed. Are the constraints");
            fprintf(stderr, " (-cx/-cy/-cz options) correct?\n");
            ColourReset(STDERR_FILENO);
          }
          double number = (double)rand() / ((double)RAND_MAX + 1);
          random.x = number * constraint_box.Length.x + constraint[0].x;
          number = (double)rand() / ((double)RAND_MAX + 1);
          random.y = number * constraint_box.Length.y + constraint[0].y;
          number = (double)rand() / ((double)RAND_MAX + 1);
          random.z = number * constraint_box.Length.z + constraint[0].z;

          min_dist = SQR(S_orig.Box.Length.x * 100);
          for (int j = 0; j < S_orig.Count.nBeadsCoor; j++) {
            int btype = S_orig.Bead[j].Type;
            /*
             * j can be added monomeric bead, so it's type can be higher than
             * the number of types
             */
            if (btype < S_orig.Count.nBeadTypes && S_orig.BeadType[btype].Use) {
              VECTOR dist;
              dist = Distance(S_orig.Bead[j].Position, random, S_new.Box.Length);
              dist.x = SQR(dist.x) + SQR(dist.y) + SQR(dist.z);
              if (dist.x < min_dist) {
                min_dist = dist.x;
              }
            }
          }
        } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
                 (highest_dist != -1 && highest_dist <= min_dist));
      } else {
        double number = (double)rand() / ((double)RAND_MAX + 1);
        random.x = number * constraint_box.Length.x + constraint[0].x;
        number = (double)rand() / ((double)RAND_MAX + 1);
        random.y = number * constraint_box.Length.y + constraint[0].y;
        number = (double)rand() / ((double)RAND_MAX + 1);
        random.z = number * constraint_box.Length.z + constraint[0].z;
      }

      // determine index of the added bead
      int id = -1;
      if (!sw) { // added beads (no --switch option)
        id = S_orig.Count.nUnbonded + i;
      } else { // switched beds (--switch option)
        for (int j = count; j < S_new.Count.nUnbonded; j++) {
          if (S_new.Bead[j].Use) { // is this an original bead to be exchanged?
            id = j;
            S_new.Bead[j].Use = false; // just exchanged (only pro forma)
            count = j + 1;
            break;
          }
        }
      }
      if (id == -1) {
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, "!!!SOME ERROR!!!");
        fprintf(stderr, "...very useful.");
        ColourReset(STDERR_FILENO);
        exit(1);
      }

      // add the new coordinate
      S_new.Bead[id].Position.x = random.x;
      S_new.Bead[id].Position.y = random.y;
      S_new.Bead[id].Position.z = random.z;

      // print number of placed beads? //{{{
      if (!silent && isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\rMonomers placed: %d", i+1);
      } //}}}
    } //}}}
    // print total number of placed beads? //{{{
    if (!silent) {
      if (isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\r                           \r");
      }
      fprintf(stdout, "\rMonomer placed: %d\n", S_add.Count.nUnbonded);
    } //}}}
    // add molecules //{{{
    // doesn't depend on --switch option as it's determined by the S_new.Molecule
    // array established earlier
    count = 0;
    for (int i = S_orig.Count.nMolecules; i < S_new.Count.nMolecules; i++) {
      int mtype = S_new.Molecule[i].Type;

      VECTOR rotated[S_new.MoleculeType[mtype].nBeads];
      VECTOR random = {0};

      // rotate the molecule randomly if desired //{{{
      if (!no_rot) {
        // random rotation axis
        random.x = (double)rand() / ((double)RAND_MAX) * 2 - 1; // a number <-1,1>
        random.y = (double)rand() / ((double)RAND_MAX) * 2 - 1;
        random.z = (double)rand() / ((double)RAND_MAX) * 2 - 1;
        dist = Length(random);
        random.x /= dist;
        random.y /= dist;
        random.z /= dist;
        // random rotation angle
        double angle = (double)rand() / ((double)RAND_MAX) * PI;
        // create rotation matrix
        struct Tensor {
          VECTOR x, y, z;
        } rot;
        double c = 1 - cos(angle);
        rot.x.x = cos(angle) + SQR(random.x) * c;
        rot.x.y = random.x * random.y * c - random.z * sin(angle);
        rot.x.z = random.x * random.z * c + random.y * sin(angle);

        rot.y.x = random.x * random.y * c + random.z * sin(angle);
        rot.y.y = cos(angle) + SQR(random.y) * c;
        rot.y.z = random.y * random.z * c - random.x * sin(angle);

        rot.z.x = random.x * random.z * c - random.y * sin(angle);
        rot.z.y = random.y * random.z * c + random.x * sin(angle);
        rot.z.z = cos(angle) + SQR(random.z) * c;
        // transform the prototype molecule (rotation matrix * coordinates)
        for (int j = 0; j < S_new.MoleculeType[mtype].nBeads; j++) {
          int id = S_new.Molecule[i].Bead[j];
          rotated[j].x = rot.x.x * S_new.Bead[id].Position.x
                       + rot.x.y * S_new.Bead[id].Position.y
                       + rot.x.z * S_new.Bead[id].Position.z;
          rotated[j].y = rot.y.x * S_new.Bead[id].Position.x
                       + rot.y.y * S_new.Bead[id].Position.y
                       + rot.y.z * S_new.Bead[id].Position.z;
          rotated[j].z = rot.z.x * S_new.Bead[id].Position.x
                       + rot.z.y * S_new.Bead[id].Position.y
                       + rot.z.z * S_new.Bead[id].Position.z;
        }
      } else { // don't rotate
        for (int j = 0; j < S_new.MoleculeType[mtype].nBeads; j++) {
          int id = S_new.Molecule[i].Bead[j];
          rotated[j].x = S_new.Bead[id].Position.x;
          rotated[j].y = S_new.Bead[id].Position.y;
          rotated[j].z = S_new.Bead[id].Position.z;
        }
      } //}}}

      // first bead's distance from specified bead typtes is checked //{{{
      // first bead can have coordinates [0,0,0] or such that the molecule's geometric centre is [0,0,0] (if -gc is used)
      if (lowest_dist != -1 || highest_dist != -1) {
        int tries = 0;
        double min_dist;
        do {
          tries++;
          if (tries == 1000000) {
            ColourChange(STDERR_FILENO, YELLOW);
            fprintf(stderr, "\nWarning: million attempts");
            fprintf(stderr, " to place a bead failed. Are the constraints");
            fprintf(stderr, " (-cx/-cy/-cz options) correct?\n");
            ColourReset(STDERR_FILENO);
          }
          double number = (double)rand() / ((double)RAND_MAX + 1);
          random.x = number * constraint_box.Length.x + constraint[0].x;
          number = (double)rand() / ((double)RAND_MAX + 1);
          random.y = number * constraint_box.Length.y + constraint[0].y;
          number = (double)rand() / ((double)RAND_MAX + 1);
          random.z = number * constraint_box.Length.z + constraint[0].z;

          min_dist = SQR(S_new.Box.Length.x) +
                     SQR(S_new.Box.Length.y) +
                     SQR(S_new.Box.Length.z);
          for (int j = 0; j < S_orig.Count.nBeadsCoor; j++) {
            int btype_j = S_orig.Bead[j].Type;
            /*
             * j can be added monomeric bead, so it's type can be higher than
             * the number of types
             */
            if (btype_j < S_orig.Count.nBeadTypes && S_orig.BeadType[btype_j].Use) {
              dist = Length(Distance(S_orig.Bead[j].Position,
                                     random, S_orig.Box.Length));
              if (dist < min_dist) {
                min_dist = dist;
              }
            }
          }
        } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
                 (highest_dist != -1 && highest_dist <= min_dist));
      } else { // no '-ld' or '-hd' options
        double number = (double)rand() / ((double)RAND_MAX + 1);
        random.x = number * constraint_box.Length.x + constraint[0].x;
        number = (double)rand() / ((double)RAND_MAX + 1);
        random.y = number * constraint_box.Length.y + constraint[0].y;
        number = (double)rand() / ((double)RAND_MAX + 1);
        random.z = number * constraint_box.Length.z + constraint[0].z;
      } //}}}

      // place the rest of the molecule //{{{
      for (int j = 0; j < S_new.MoleculeType[mtype].nBeads; j++) {
        int id = S_new.Molecule[i].Bead[j];
        S_new.Bead[id].Position.x = random.x + rotated[j].x;
        S_new.Bead[id].Position.y = random.y + rotated[j].y;
        S_new.Bead[id].Position.z = random.z + rotated[j].z;
      } //}}}

      // print number of placed molecules? //{{{
      if (!silent && isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\rMolecules placed: %d", i-S_orig.Count.nMolecules+1);
      } //}}}
    } //}}}
    // print total number of placed molecules? //{{{
    if (!silent) {
      if (isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\r                                             \r");
      }
      fprintf(stdout, "Molecules placed: %3d\n", S_add.Count.nMolecules);
    } //}}}
  } //}}}

  // write data to output files //{{{
  // vsf file
  VtfWriteStruct(output_vsf, S_new);
  // .vcf file
  FILE *out = OpenFile(output_vcf, "w");
  PrintByline(out, argc, argv);
  for (int i = 0; i < S_new.Count.nBeadsTotal; i++) {
    S_new.Bead[i].Use = true; // TODO change somewhere (use different flag for sw)
  }
  VtfWriteCoorIndexed(out, stuff, S_new);
  fclose(out);
  // xyz file (if -xyz option is present)
  if (strlen(output_xyz) > 0) {
    FILE *xyz = OpenFile(output_xyz, "w");
    XyzWriteCoor(xyz, S_new);
    fclose(xyz);
  }
  //}}}

  count=0;
  SYSTEM sys_test = CopySystem(S_orig);
//VerboseOutput(S_new);
//PrintMoleculeType(S_orig);
//PrintBox(S_orig.Box);
//PrintBeadType(S_orig);
//PrintCounts(S_orig);
//PrintCounts(sys_test);
//PrintBeadType(sys_test);
//PrintBox(sys_test.Box);
  ConcatenateSystems(&sys_test, S_add, S_add.Box);
//VerboseOutput(sys_test);
//PrintCounts(sys_test);
//PrintBeadType(sys_test);

//printf("%sORIGINAL\n", Green());
//PrintCounts(S_orig);
//PrintBeadType(S_orig);
//PrintBead(S_orig);
//printf("%sCONCATENATED\n", Magenta());
//PrintCounts(sys_test);
//PrintBeadType(sys_test);
//PrintBead(sys_test);
//puts(ColourReset());
  // TODO double free error appears when sys_test is freed
  //      maybe I shouldn't use S_out = S_in in ConcatenateSystems()?
  FreeSystem(&sys_test);

  // free memory - to make valgrind happy //{{{
  FreeSystem(&S_orig);
  FreeSystem(&S_add);
  FreeSystem(&S_new);
  free(stuff);
  //}}}

  return 0;
}
