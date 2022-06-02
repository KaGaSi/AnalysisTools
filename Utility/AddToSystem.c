#include "../AnalysisTools.h"
// TODO: split into two utilities - adding existing (vtf) configuration and
//       generating addition from FIELD?
// TODO: --random for -vtf option (i.e., place added system's components
//       randomly in a new box)
// TODO: remove --detailed option

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
        strcmp(argv[i], "--detailed") != 0 &&
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
  char file_coor[LINE] = "", // unchanged => new system
       file_struct[LINE] = "";
  bool vtf = false;
  if (strcmp(argv[++count], "--") != 0) {
    // test that <input> filename ends with '.vcf' or '.vtf'
    snprintf(file_coor, LINE, "%s", argv[count]);
    if (!InputCoor(&vtf, file_coor, file_struct)) {
      Help(argv[0], true);
      exit(1);
    }
  } else {
    strcpy(file_struct, "in.vsf"); // won't be used
  } //}}}

  // <out.vsf> - output vsf file //{{{
  char file_out_struct[LINE] = "";
  snprintf(file_out_struct, LINE, "%s", argv[++count]);
  // test that <out.vsf> filename ends with '.vsf'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(file_out_struct, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output vcf file //{{{
  char file_out_coor[LINE] = "";
  snprintf(file_out_coor, LINE, "%s", argv[++count]);
  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(file_out_coor, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions(argc, argv, file_struct, LINE, &verbose, &silent, &detailed);

  // -f <add> - FIELD-like file with molecules to add //{{{
  char file_add_field[LINE] = "";
  if (FileOption(argc, argv, "-f", file_add_field, LINE)) {
    exit(1);
  }
  if (file_add_field[0] == '\0') {
    strcpy(file_add_field, "FIELD");
  } //}}}

  // -vtf <vsf> <vcf> - vtf file(s) to use instead of FIELD //{{{
  char file_add_struct[LINE] = "", file_add_coor[LINE] = "";
  // 1) vsf file
  if (FileOption(argc, argv, "-vtf", file_add_struct, LINE)) {
    exit(1);
  }
  // 2) if vsf file exists, look for vcf
//bool vtf_add = true; // if -vtf is present present, is the file a vtf format?
  if (strlen(file_add_struct) > 0) {
    ext = 2;
    strcpy(extension[0], ".vsf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(file_add_struct, ext, extension) == -1) {
      Help(argv[0], true);
      exit(1);
    }
    // TODO: change == to != and remove else
    if (file_add_struct[strlen(file_add_struct)-2] == 't') {
    // if *.vtf file, use it as vcf too
      snprintf(file_add_coor, LINE, "%s", file_add_struct);
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
            if (FileOption(argc, argv, "-vtf", file_add_coor, LINE)) {
              exit(1);
            }
            // restore vsf filename so the command in unchanged
            strcpy(argv[i+1], temp);
            // coordinate file must be vcf, because there's already vsf
            ext = 1;
            strcpy(extension[0], ".vcf");
            ext = ErrorExtension(file_add_coor, ext, extension);
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
  if (strlen(file_add_struct) == 0 && offset[0] != 1000000) {
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
  char file_out_xyz[LINE] = "";
  if (FileOption(argc, argv, "-xyz", file_out_xyz, LINE)) {
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
  if (strlen(file_coor) > 0) { // is there an input coordinate file?
    S_orig = VtfReadStruct_old(file_struct, detailed);
  } else {
    InitSystem(&S_orig);
  } //}}}

  // -xb <name(s)> - specify what bead types to exchange //{{{
  bool sw = BoolOption(argc, argv, "-xb"); // is -xb present?
  // error - if -xb is used, input system must be present
  if (sw && strlen(file_coor) == 0) {
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
  for (int i = 0; i < S_orig.Count.BeadType; i++) {
    S_orig.BeadType[i].Write = S_orig.BeadType[i].Use;
    S_orig.BeadType[i].Use = false; // this flag may be used later
    if (S_orig.BeadType[i].Write) {
      all_false = false; // '-xb' option is present
    }
  }
  if (all_false) {
    for (int i = 0; i < S_orig.Count.BeadType; i++) {
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

  // array for the timestep preamble
  char *stuff = calloc(LINE, sizeof *stuff);

  // open input coordinate file //{{{
  FILE *vcf;
  if (strlen(file_coor) > 0) {
    vcf = OpenFile(file_coor, "r");
    count= 0;
    if (!silent) {
      fprintf(stdout, "Using step %6d\n", ++count);
    }
    int file_line_count = 0, count_vcf = 0;
    VtfReadTimestep(vcf, file_coor, file_struct, &S_orig, &file_line_count,
                    count_vcf, stuff);
    if (S_orig.Count.BeadCoor == -1) {
      // TODO shouldn't this be handled better in the Read.c?
      S_orig.Count.BeadCoor = 0;
      strcpy(ERROR_MSG, "no coordinate for input system found");
      PrintWarning();
      WarnPrintFile(file_coor, file_struct);
      putc('\n', stderr);
    }
    fclose(vcf);
  } //}}}

  // print original system (if present) //{{{
  if (verbose && strlen(file_coor) > 0) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(S_orig);
    if (start > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", start);
    }
  } //}}}

  SYSTEM S_add;

  // TODO FIELD must be completely redone
  if (strlen(file_add_struct) == 0) { // read stuff to be added from FIELD //{{{
    InitSystem(&S_add);
//  ReadField(input_add, '\0', &Counts_add, &S_add.BeadType, &S_add.Bead,
//            &Index_add, &S_add.MoleculeType, &S_add.Molecule,
//            &bond_type, &angle_type, &dihedral_type);
//  S_add.Box.Length = S_orig.Box.Length; //}}}
  } else { // read stuff to add from vtf file(s) ('-vtf' option) //{{{
    S_add = VtfReadStruct_old(file_add_struct, false);
    // read coordinates
    vcf = OpenFile(file_add_coor, "r");
    int file_line_count = 0, count_vcf = 0;
    VtfReadTimestep(vcf, file_add_coor, file_add_struct, &S_add, &file_line_count, count_vcf, stuff);
    fclose(vcf);
    VECTOR rotated[S_add.Count.BeadCoor];
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
      // transform the added system (rotation matrix * coordinates)
      ToFractionalCoor(S_add.Count.Bead, &S_add.Bead, S_add.Box);
      for (int i = 0; i < S_add.Count.BeadCoor; i++) {
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
      for (int i = 0; i < S_add.Count.BeadCoor; i++) {
        S_add.Bead[i].Position.x = rotated[i].x + offset[0] + S_add.Box.Length.x / 2;
        S_add.Bead[i].Position.y = rotated[i].y + offset[1] + S_add.Box.Length.y / 2;
        S_add.Bead[i].Position.z = rotated[i].z + offset[2] + S_add.Box.Length.z / 2;
      }
      FromFractionalCoor(S_add.Count.Bead, &S_add.Bead, S_add.Box);
      //}}}
    } else { // don't rotate //{{{
      ToFractionalCoor(S_add.Count.Bead, &S_add.Bead, S_add.Box);
      for (int i = 0; i < S_add.Count.BeadCoor; i++) {
        int id = S_add.BeadCoor[i];
        S_add.Bead[id].Position.x += offset[0];
        S_add.Bead[id].Position.y += offset[1];
        S_add.Bead[id].Position.z += offset[2];
      }
      FromFractionalCoor(S_add.Count.Bead, &S_add.Bead, S_add.Box);
    } //}}}
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
  for (int i = 0; i < S_orig.Count.Bead; i++) {
    int btype = S_orig.Bead[i].Type;
    if (S_orig.Bead[i].Molecule == -1 && S_orig.BeadType[btype].Write) {
      can_be_exchanged++;
    }
  }
  // count beads to be added
  if (sw && S_add.Count.BeadCoor > can_be_exchanged) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, "insufficient beads to exchange for new ones\n");
    fprintf(stderr, "     Exchangeable beads in the original system: ");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d\n", can_be_exchanged);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, "     Beads to be added: ");
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d\n\n", S_add.Count.BeadCoor);
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}

  // if '-gc' is used, put prototypes' geometric centres to (0,0,0) //{{{
  if (com) {
    for (int i = 0; i < S_add.Count.Molecule; i++) {
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
  for (int i = 0; i < S_orig.Count.BeadCoor; i++) {
    S_orig.Bead[i].Use = false;
  } //}}}
  count = 0; // counts bead in the original Bead[] struct
  for (int i = 0; i < S_add.Count.BeadCoor; i++) {
    for (; count < S_orig.Count.BeadCoor; count++) {
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
    S_new.Count.BeadCoor = S_orig.Count.BeadCoor;
    S_new.Count.Bead = S_orig.Count.Bead;
    S_new.Count.Bonded = S_orig.Count.Bonded + S_add.Count.Bonded;
    S_new.Count.Unbonded = S_orig.Count.BeadCoor - S_new.Count.Bonded;
    S_new.Count.BondType = S_add.Count.BondType;
    S_new.Count.AngleType = S_add.Count.AngleType;
    S_new.Count.Molecule = S_orig.Count.Molecule + S_add.Count.Molecule;
    // fill BeadType struct for the new system
    S_new.Count.BeadType = S_orig.Count.BeadType;
    // 1) copy original BeadType
    CopyBeadType(S_new.Count.BeadType, &S_new.BeadType, S_orig.BeadType, 3);
    // 2) add new bead types - the check is based only on Name //{{{
    for (int i = 0; i < S_add.Count.BeadType; i++) {
      bool new = true;
      for (int j = 0; j < S_orig.Count.BeadType; j++) {
        if (strcmp(S_add.BeadType[i].Name, S_orig.BeadType[j].Name) == 0) {
          new = false;
          S_new.BeadType[j].Number += S_add.BeadType[i].Number;
          break;
        }
      }
      if (new) {
        int type = S_new.Count.BeadType;
        S_new.BeadType = realloc(S_new.BeadType, sizeof (BEADTYPE) * (type + 1));
        S_new.BeadType[type] = S_add.BeadType[i];
        S_new.Count.BeadType++;
      }
    } //}}}
    // fill MoleculeType struct for the new system
    S_new.Count.MoleculeType = S_orig.Count.MoleculeType;
    S_new.MoleculeType = realloc(S_new.MoleculeType, sizeof (MOLECULETYPE) * S_new.Count.MoleculeType);
    // copy original MoleculeType to _new //{{{
    for (int i = 0; i < S_new.Count.MoleculeType; i++) {
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
    for (int i = 0; i < S_add.Count.MoleculeType; i++) {
      bool new = true;
      for (int j = 0; j < S_orig.Count.MoleculeType; j++) {
        if (strcmp(S_add.MoleculeType[i].Name, S_orig.MoleculeType[j].Name) == 0) {
          new = false;
          S_new.MoleculeType[j].Number += S_add.MoleculeType[i].Number;
          break;
        }
      }
      if (new) {
        int type = S_new.Count.MoleculeType;
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
        S_new.Count.MoleculeType++;
      }
    } //}}}
    // fill Bead struct for the new system
    S_new.Bead = realloc(S_new.Bead, sizeof (BEAD) * S_new.Count.BeadCoor);
    // copy unbonded beads not to be exchanged to the start of S_new.Bead //{{{
    // TODO: assumes unbonded beads are before bonded beads
    count = 0; // counts copied beads
    for (int i = 0; i < S_orig.Count.Unbonded; i++) {
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
    count = S_orig.Count.Unbonded - S_add.Count.BeadCoor; // just to be sure
    for (int i = 0; i < S_add.Count.Unbonded; i++) {
      int type = S_add.Bead[i].Type;
      int new_type = FindBeadType(S_add.BeadType[type].Name, S_new);
      S_new.Bead[count] = S_orig.Bead[i];
      S_new.Bead[count].Type = new_type;
      S_new.Bead[count].Molecule = -1;
      S_new.Bead[count].Use = true; // coordinates to be rewritten
      count++; // use count to make it consistent & easy to read
    } //}}}
    // copy the original bonded beads //{{{
    count = S_new.Count.Unbonded;
    for (int i = S_orig.Count.Unbonded; i < S_orig.Count.BeadCoor; i++) {
      S_new.Bead[count] = S_orig.Bead[i];
      S_new.Bead[count].Use = false; // coordinates to be rewritten
      count++;
    } //}}}
    // put bonded beads to be added at the very end //{{{
    count = S_new.Count.Unbonded + S_orig.Count.Bonded;
    for (int i = S_add.Count.Unbonded; i < S_add.Count.BeadCoor; i++) {
      int type = S_add.Bead[i].Type;
      int new_type = FindBeadType(S_add.BeadType[type].Name, S_new);
      S_new.Bead[count] = S_add.Bead[i];
      S_new.Bead[count].Type = new_type;
      S_new.Bead[count].Molecule = S_add.Bead[i].Molecule + S_orig.Count.Molecule;
      S_new.Bead[count].Use = true; // coordinates to be rewritten
      count++; // use count to make it consistent & easy to read
    } //}}}
    // alocate new molecule struct
    S_new.Molecule = realloc(S_new.Molecule, sizeof (MOLECULE) * S_new.Count.Molecule);
    // copy original molecules to _new struct //{{{
    for (int i = 0; i < S_orig.Count.Molecule; i++) {
      int type = S_orig.Molecule[i].Type;
      S_new.Molecule[i].Type = type;
      S_new.Molecule[i].Bead = malloc(sizeof *S_new.Molecule[i].Bead * S_new.MoleculeType[type].nBeads);
      for (int j = 0; j < S_new.MoleculeType[type].nBeads; j++) {
        S_new.Molecule[i].Bead[j] = S_orig.Molecule[i].Bead[j] - S_add.Count.Bonded;
      }
    } //}}}
    // put _add molecules into _new struct //{{{
    count = S_new.Count.BeadCoor - S_add.Count.Bonded;
    for (int i = 0; i < S_add.Count.Molecule; i++) {
      int add_type = S_add.Molecule[i].Type;
      int new_type = FindMoleculeType(S_add.MoleculeType[add_type].Name, S_new);
      int new_i = S_orig.Count.Molecule + i;
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
    BOX Box_new = S_new.Box; // TODO box set up previously; maybe make it a different BOX
    S_new = CopySystem(S_orig);
    fflush(stdout);
    ConcatenateSystems(&S_new, S_add, Box_new);
    PruneSystem(&S_new);
//  VerboseOutput(S_new);
//  PrintMolecule(S_new);
  } //}}}

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(S_new);
//  PrintBondTypes2(S_new.TypesOfBonds, bond_type);
  } //}}}

  // add beads randomly if FIELD-like file is used //{{{
  double dist;
  if (strlen(file_add_struct) == 0) {
    count = 0;
    // add monomeric beads //{{{
    for (int i = 0; i < S_add.Count.Unbonded; i++) {
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
          for (int j = 0; j < S_orig.Count.BeadCoor; j++) {
            int btype = S_orig.Bead[j].Type;
            /*
             * j can be added monomeric bead, so it's type can be higher than
             * the number of types
             */
            if (btype < S_orig.Count.BeadType && S_orig.BeadType[btype].Use) {
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
        id = S_orig.Count.Unbonded + i;
      } else { // switched beds (--switch option)
        for (int j = count; j < S_new.Count.Unbonded; j++) {
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
      fprintf(stdout, "\rMonomer placed: %d\n", S_add.Count.Unbonded);
    } //}}}
    // add molecules //{{{
    // doesn't depend on --switch option as it's determined by the S_new.Molecule
    // array established earlier
    count = 0;
    for (int i = S_orig.Count.Molecule; i < S_new.Count.Molecule; i++) {
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
          for (int j = 0; j < S_orig.Count.BeadCoor; j++) {
            int btype_j = S_orig.Bead[j].Type;
            /*
             * j can be added monomeric bead, so it's type can be higher than
             * the number of types
             */
            if (btype_j < S_orig.Count.BeadType && S_orig.BeadType[btype_j].Use) {
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
        fprintf(stdout, "\rMolecules placed: %d", i-S_orig.Count.Molecule+1);
      } //}}}
    } //}}}
    // print total number of placed molecules? //{{{
    if (!silent) {
      if (isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\r                                             \r");
      }
      fprintf(stdout, "Molecules placed: %3d\n", S_add.Count.Molecule);
    } //}}}
  } //}}}

  // write data to output files //{{{
  // vsf file
  VtfWriteStruct(file_out_struct, S_new);
  // .vcf file
  FILE *out = OpenFile(file_out_coor, "w");
  PrintByline(out, argc, argv);
  for (int i = 0; i < S_new.Count.Bead; i++) {
    S_new.Bead[i].Use = true; // TODO change somewhere (use different flag for sw)
  }
  VtfWriteCoorIndexed(out, stuff, S_new);
  fclose(out);
  // xyz file (if -xyz option is present)
  if (strlen(file_out_xyz) > 0) {
    FILE *xyz = OpenFile(file_out_xyz, "w");
    XyzWriteCoor(xyz, S_new);
    fclose(xyz);
  }
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&S_orig);
  FreeSystem(&S_add);
  FreeSystem(&S_new);
  free(stuff);
  //}}}

  return 0;
}
