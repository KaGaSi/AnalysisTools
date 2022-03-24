#include "../AnalysisTools.h"
int *InFile;
// TODO: split into two utilities - adding existing (vtf) configuration and
//       generating addition from FIELD?
// TODO: -offset for -vtf option
// TODO: --random for -vtf option (i.e., place added system's components
//       randomly in a new box)
// TODO: triclinic box

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
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);

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
  bool vtf_add = true; // if -vtf is present present, is the file a vtf format?
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
            } else {
              vtf_add = false;
            }
            break;
          } else { // missing vcf file name
            ErrorPrintError_old();
            ColourText(STDERR_FILENO, YELLOW);
            fprintf(stderr, "-vtf");
            ColourText(STDERR_FILENO, RED);
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
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-offset");
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, " - three numbers required\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // Warning - missing -vtf option
  if (strlen(add_vsf) == 0 && offset[0] != 1000000) {
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "\nWarning: ");
    ColourText(STDERR_FILENO, CYAN);
    fprintf(stderr, "-offset");
    ColourText(STDERR_FILENO, YELLOW);
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
      ColourText(STDERR_FILENO, RED);
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
      ColourText(STDERR_FILENO, RED);
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
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-cx");
    ColourText(STDERR_FILENO, RED);
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
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-cy");
    ColourText(STDERR_FILENO, RED);
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
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-cz");
    ColourText(STDERR_FILENO, RED);
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

  // use centre of mass for distance check of new molecules //{{{
  bool com = BoolOption(argc, argv, "-gc"); //}}}

  // rotate added molecules?
  bool no_rot = BoolOption(argc, argv, "--no-rotate");

  // define new box size //{{{
  double box_option[100] = {-1};
  if (MultiDoubleOption(argc, argv, "-b", &count, box_option)) {
    exit(1);
  }
  if (count != 3) {
    ErrorPrintError_old();
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-b");
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, " - three non-negative numbers required\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  } //}}}

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
  BEADTYPE *bt_orig;
  MOLECULETYPE *mt_orig;
  BEAD *bead_orig;
  int *Index_orig;
  MOLECULE *mol_orig;
  COUNTS Counts_orig = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box_orig = InitBox; // triclinic box dimensions and angles
  Box_orig.Length.x = 0;
  Box_orig.Length.y = 0; // TODO: why?
  Box_orig.Length.z = 0;
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  if (strlen(input_coor) != 0) { // is there an input coordinate file?
    FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
                &Box_orig, &Counts_orig, &bt_orig, &bead_orig, &Index_orig,
                &mt_orig, &mol_orig);
  } else { // if there's no input coordinate file, just allocate some memory
    bt_orig = calloc(1, sizeof (BEADTYPE));
    mt_orig = calloc(1, sizeof (MOLECULETYPE));
    bead_orig = calloc(1, sizeof (BEAD));
    Index_orig = calloc(1, sizeof *Index_orig);
    mol_orig = calloc(1, sizeof (MOLECULE));
  } //}}}

  // -xb <name(s)> - specify what bead types to exchange //{{{
  bool sw = BoolOption(argc, argv, "-xb"); // is -xb present?
  // error - if -xb is used, 
  if (sw && strlen(input_coor) == 0) {
    ErrorPrintError_old();
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-xb");
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, " - <input> file must be present\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], true);
    exit(1);
  }
  // which beads to exchange?
  if (BeadTypeOption(argc, argv, "-xb", false, Counts_orig, &bt_orig)) {
    exit(0);
  }
  // use Write flag to decide which bead types to use
  bool all_false = true; // no '-xb' option
  for (int i = 0; i < Counts_orig.TypesOfBeads; i++) {
    bt_orig[i].Write = bt_orig[i].Use;
    bt_orig[i].Use = false; // this flag may be used later
    if (bt_orig[i].Write) {
      all_false = false; // '-xb' option is present
    }
  }
  if (all_false) {
    for (int i = 0; i < Counts_orig.TypesOfBeads; i++) {
      if (bt_orig[i].Charge == 0) {
        bt_orig[i].Write = true;
      }
    }
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", false, Counts_orig, &bt_orig)) {
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
    VerboseOutput(Counts_orig, bt_orig, bead_orig, mt_orig, mol_orig);
    if (start > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", start);
    }
  } //}}}

  // array for the timestep preamble
  char *stuff = calloc(LINE, sizeof *stuff);

  // TODO: BoxLength jsut to add to reading coors
  BOX Box;
  // open input coordinate file //{{{
  FILE *vcf;
  if (strlen(input_coor) > 0) {
    if ((vcf = fopen(input_coor, "r")) == NULL) {
      ErrorFileOpen(input_coor, 'r');
      exit(1);
    }
    SkipVtfStructure(vcf, struct_lines);
    count = SkipCoorSteps(vcf, input_coor, Counts_orig, start, silent);
    if (!silent) {
      fprintf(stdout, "Using step %6d\n", ++count);
    }
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box, Counts_orig,
                       Index_orig, &bead_orig, &stuff);
    fclose(vcf);
  } //}}}

  // create structures for added stuff //{{{
  COUNTS Counts_add = InitCounts;
  MOLECULE *mol_add;
  MOLECULETYPE *mt_add;
  BEADTYPE *bt_add;
  BEAD *bead_add;
  int *Index_add;
  PARAMS *bond_type;
  PARAMS *angle_type;
  PARAMS *dihedral_type;
  BOX Box_add;
  //}}}

  if (strlen(add_vsf) == 0) { // read stuff to be added from FIELD //{{{
    ReadField(input_add, '\0', &Counts_add, &bt_add, &bead_add,
              &Index_add, &mt_add, &mol_add,
              &bond_type, &angle_type, &dihedral_type);
    Box_add.Length = Box_orig.Length; //}}}
  } else { // read stuff to add from vtf file(s) ('-vtf' option) //{{{
    bool indexed_add;
    int struct_lines_add;
    FullVtfRead(add_vsf, input_coor_add, false, vtf_add, &indexed_add,
                &struct_lines_add, &Box_add, &Counts_add,
                &bt_add, &bead_add, &Index_add, &mt_add, &mol_add);
    // read coordinates
    if ((vcf = fopen(input_coor_add, "r")) == NULL) {
      ErrorFileOpen(input_coor_add, 'r');
      exit(1);
    }
    SkipVtfStructure(vcf, struct_lines_add);
    ReadVcfCoordinates(indexed_add, input_coor_add, vcf, &Box, Counts_add,
                       Index_add, &bead_add, &stuff);
    fclose(vcf);
    // TODO: !no_rot? ...shouldn't -vtf be this by default?
    VECTOR rotated[Counts_add.BeadsCoor];
    if (!no_rot) {
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
      for (int i = 0; i < Counts_add.BeadsCoor; i++) {
        rotated[i].x = rot.x.x * (bead_add[i].Position.x - Box_add.Length.x / 2)
                     + rot.x.y * (bead_add[i].Position.y - Box_add.Length.y / 2)
                     + rot.x.z * (bead_add[i].Position.z - Box_add.Length.z / 2);
        rotated[i].y = rot.y.x * (bead_add[i].Position.x - Box_add.Length.x / 2)
                     + rot.y.y * (bead_add[i].Position.y - Box_add.Length.y / 2)
                     + rot.y.z * (bead_add[i].Position.z - Box_add.Length.z / 2);
        rotated[i].z = rot.z.x * (bead_add[i].Position.x - Box_add.Length.x / 2)
                     + rot.z.y * (bead_add[i].Position.y - Box_add.Length.y / 2)
                     + rot.z.z * (bead_add[i].Position.z - Box_add.Length.z / 2);
      }
      for (int i = 0; i < Counts_add.BeadsCoor; i++) {
        bead_add[i].Position.x = rotated[i].x + offset[0] + Box_add.Length.x / 2;
        bead_add[i].Position.y = rotated[i].y + offset[1] + Box_add.Length.y / 2;
        bead_add[i].Position.z = rotated[i].z + offset[2] + Box_add.Length.z / 2;
      }
    } else { // don't rotate
      for (int i = 0; i < Counts_add.BeadsCoor; i++) {
        bead_add[i].Position.x += offset[0];
        bead_add[i].Position.y += offset[1];
        bead_add[i].Position.z += offset[2];
      }
    }
    // allocate memory only to free it later
    bond_type = calloc(1, sizeof (PARAMS));
    angle_type = calloc(1, sizeof (PARAMS));
    dihedral_type = calloc(1, sizeof (PARAMS));
  } //}}}

  // set final box size //{{{
  /*
   * if -b isn't used, set box size as the larger of the dimensions from
   * original and to-be-added systems
   */
  BOX Box_new = InitBox;
  Box_new.Length.x = 0;
  Box_new.Length.y = 0;
  Box_new.Length.z = 0;
  if (box_option[0] == -1) {
    if (Box_add.Length.x > Box_orig.Length.x) {
      Box_new.Length.x = Box_add.Length.x;
    } else {
      Box_new.Length.x = Box_orig.Length.x;
    }
    if (Box_add.Length.y > Box_orig.Length.y) {
      Box_new.Length.y = Box_add.Length.y;
    } else {
      Box_new.Length.y = Box_orig.Length.y;
    }
    if (Box_add.Length.z > Box_orig.Length.z) {
      Box_new.Length.z = Box_add.Length.z;
    } else {
      Box_new.Length.z = Box_orig.Length.z;
    }
  } else {
    Box_new.Length.x = box_option[0];
    Box_new.Length.y = box_option[1];
    Box_new.Length.z = box_option[2];
  } //}}}

  // define 'box' for additions using constraints (-c{x,y,z} options)//{{{
  BOX constraint_box = InitBox;
  if (constraint[1].x != 0) {
    constraint_box.Length.x = constraint[1].x - constraint[0].x;
  } else {
    constraint_box.Length.x = Box_new.Length.x;
  }
  if (constraint[1].y != 0) {
    constraint_box.Length.y = constraint[1].y - constraint[0].y;
  } else {
    constraint_box.Length.y = Box_new.Length.y;
  }
  if (constraint[1].z != 0) {
    constraint_box.Length.z = constraint[1].z - constraint[0].z;
  } else {
    constraint_box.Length.z = Box_new.Length.z;
  } //}}}

  // error - no box size //{{{
  if (Box_new.Length.x == 0 || Box_new.Length.y == 0 || Box_new.Length.z == 0) {
    ErrorPrintError_old();
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, "zero box size for the new system\n\n");
    ColourReset(STDERR_FILENO);
    Help(argv[0], 1);
    exit(1);
  } //}}}

  // check number of exchangeable beads //{{{
  int can_be_exchanged = 0;
  for (int i = 0; i < Counts_orig.BeadsTotal; i++) {
    int btype = bead_orig[i].Type;
    if (bead_orig[i].Molecule == -1 && bt_orig[btype].Write) {
      can_be_exchanged++;
    }
  }
  // count beads to be added
  if (sw && Counts_add.BeadsCoor > can_be_exchanged) {
    ErrorPrintError_old();
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, "insufficient beads to exchange for new ones\n");
    fprintf(stderr, "     Exchangeable beads in the original system: ");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d\n", can_be_exchanged);
    ColourText(STDERR_FILENO, RED);
    fprintf(stderr, "     Beads to be added: ");
    ColourText(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%d\n\n", Counts_add.BeadsCoor);
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}

  // if '-gc' is used, put prototypes' geometric centres to (0,0,0) //{{{
  if (com) {
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int mtype = mol_add[i].Type;
      VECTOR geom_centre;
      geom_centre.x = 0;
      geom_centre.y = 0;
      geom_centre.z = 0;
      for (int j = 0; j < mt_add[mtype].nBeads; j++) {
        int id = mol_add[i].Bead[j];
        geom_centre.x += bead_add[id].Position.x;
        geom_centre.y += bead_add[id].Position.y;
        geom_centre.z += bead_add[id].Position.z;
      }
      geom_centre.x /= mt_add[mtype].nBeads;
      geom_centre.y /= mt_add[mtype].nBeads;
      geom_centre.z /= mt_add[mtype].nBeads;
      for (int j = 0; j < mt_add[mtype].nBeads; j++) {
        int id = mol_add[i].Bead[j];
        bead_add[id].Position.x -= geom_centre.x;
        bead_add[id].Position.y -= geom_centre.y;
        bead_add[id].Position.z -= geom_centre.z;
      }
    }
  } //}}}

  // print what is to be added //{{{
  if (verbose) {
    fprintf(stdout, "\nBEADS AND MOLECULES TO ADD\n");
    VerboseOutput(Counts_add, bt_add, bead_add, mt_add, mol_add);
  } //}}}

  /* decide which beads to exchange //{{{
   * i.e., give them Bead[].Flag = true); has effect only if --switch is used
   */
  // zeroize Bead[].Flag //{{{
  for (int i = 0; i < Counts_orig.BeadsCoor; i++) {
    bead_orig[i].Use = false;
  } //}}}
  count = 0; // counts bead in the original Bead[] struct
  for (int i = 0; i < Counts_add.BeadsCoor; i++) {
    for (; count < Counts_orig.BeadsCoor; count++) {
      int type = bead_orig[count].Type;
      if (bt_orig[type].Write && bead_orig[count].Molecule == -1) {
        bead_orig[count].Use = true; // exchange bead 'count'
        break;
      }
    }
    count++; // loop didn't update count because of the break
  } //}}}

  // create structures for new system //{{{
  COUNTS Counts_new = InitCounts;
  BEADTYPE *bt_new = calloc(1, sizeof (BEADTYPE));
  MOLECULETYPE *mt_new = calloc(1, sizeof (MOLECULETYPE));
  BEAD *bead_new = calloc(1, sizeof (BEAD));
  MOLECULE *mol_new = calloc(1, sizeof (MOLECULE));
  int *Index_new = calloc(1, sizeof *Index_new); //}}}

  // join original and added systems (depending on '--switch' mode)
  if (sw) { // switch old beads for new ones? //{{{
    Counts_new.BeadsCoor = Counts_orig.BeadsCoor;
    Counts_new.BeadsTotal = Counts_orig.BeadsTotal;
    Counts_new.Bonded = Counts_orig.Bonded + Counts_add.Bonded;
    Counts_new.Unbonded = Counts_orig.BeadsCoor - Counts_new.Bonded;
    Counts_new.TypesOfBonds = Counts_add.TypesOfBonds;
    Counts_new.TypesOfAngles = Counts_add.TypesOfAngles;
    Counts_new.Molecules = Counts_orig.Molecules + Counts_add.Molecules;
    // fill BeadType struct for the new system
    Counts_new.TypesOfBeads = Counts_orig.TypesOfBeads;
    // 1) copy original BeadType
    CopyBeadType(Counts_new.TypesOfBeads, &bt_new, bt_orig, 3);
    // 2) add new bead types - the check is based only on Name //{{{
    for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
      bool new = true;
      for (int j = 0; j < Counts_orig.TypesOfBeads; j++) {
        if (strcmp(bt_add[i].Name, bt_orig[j].Name) == 0) {
          new = false;
          bt_new[j].Number += bt_add[i].Number;
          break;
        }
      }
      if (new) {
        int type = Counts_new.TypesOfBeads;
        bt_new = realloc(bt_new, sizeof (BEADTYPE) * (type + 1));
        bt_new[type] = bt_add[i];
        Counts_new.TypesOfBeads++;
      }
    } //}}}
    // fill MoleculeType struct for the new system
    Counts_new.TypesOfMolecules = Counts_orig.TypesOfMolecules;
    mt_new = realloc(mt_new, sizeof (MOLECULETYPE) * Counts_new.TypesOfMolecules);
    // copy original MoleculeType to _new //{{{
    for (int i = 0; i < Counts_new.TypesOfMolecules; i++) {
      mt_new[i] = mt_orig[i];
      mt_new[i].Bead = malloc(sizeof *mt_new[i].Bead * mt_new[i].nBeads);
      for (int j = 0; j < mt_new[i].nBeads; j++) {
        mt_new[i].Bead[j] = mt_orig[i].Bead[j];
      }
      mt_new[i].Bond = malloc(sizeof *mt_new[i].Bond * mt_new[i].nBonds);
      for (int j = 0; j < mt_new[i].nBonds; j++) {
        mt_new[i].Bond[j][0] = mt_orig[i].Bond[j][0];
        mt_new[i].Bond[j][1] = mt_orig[i].Bond[j][1];
        mt_new[i].Bond[j][2] = mt_orig[i].Bond[j][2];
      }
    } //}}}
    // add new molecule types - check if their the same based only on Name //{{{
    for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
      bool new = true;
      for (int j = 0; j < Counts_orig.TypesOfMolecules; j++) {
        if (strcmp(mt_add[i].Name, mt_orig[j].Name) == 0) {
          new = false;
          mt_new[j].Number += mt_add[i].Number;
          break;
        }
      }
      if (new) {
        int type = Counts_new.TypesOfMolecules;
        mt_new = realloc(mt_new, sizeof (MOLECULETYPE) * (type + 1));
        mt_new[type] = mt_add[i];
        mt_new[type].Bead = malloc(sizeof *mt_new[type].Bead *
                                   mt_new[type].nBeads);
        for (int j = 0; j < mt_new[type].nBeads; j++) {
          int old_type = mt_add[i].Bead[j];
          int btype = FindBeadType(bt_add[old_type].Name,
                                   Counts_new, bt_new);
          mt_new[type].Bead[j] = btype;
        }
        mt_new[type].Bond = malloc(sizeof *mt_new[i].Bond *
                                   mt_new[type].nBonds);
        for (int j = 0; j < mt_new[type].nBonds; j++) {
          mt_new[type].Bond[j][0] = mt_add[i].Bond[j][0];
          mt_new[type].Bond[j][1] = mt_add[i].Bond[j][1];
          mt_new[type].Bond[j][2] = mt_add[i].Bond[j][2];
        }
        if (mt_new[i].nAngles > 0) {
          mt_new[type].Angle = malloc(sizeof *mt_new[type].Angle *
                                      mt_new[type].nAngles);
          for (int j = 0; j < mt_new[type].nAngles; j++) {
            mt_new[type].Angle[j][0] = mt_add[i].Angle[j][0];
            mt_new[type].Angle[j][1] = mt_add[i].Angle[j][1];
            mt_new[type].Angle[j][2] = mt_add[i].Angle[j][2];
            mt_new[type].Angle[j][3] = mt_add[i].Angle[j][3];
          }
        }
        if (mt_new[i].nDihedrals > 0) {
          mt_new[type].Dihedral = malloc(sizeof *mt_new[type].Dihedral *
                                         mt_new[type].nDihedrals);
          for (int j = 0; j < mt_new[type].nDihedrals; j++) {
            mt_new[type].Dihedral[j][0] = mt_add[i].Dihedral[j][0];
            mt_new[type].Dihedral[j][1] = mt_add[i].Dihedral[j][1];
            mt_new[type].Dihedral[j][2] = mt_add[i].Dihedral[j][2];
            mt_new[type].Dihedral[j][3] = mt_add[i].Dihedral[j][3];
            mt_new[type].Dihedral[j][4] = mt_add[i].Dihedral[j][4];
          }
        }
        Counts_new.TypesOfMolecules++;
      }
    } //}}}
    // fill Bead struct for the new system
    bead_new = realloc(bead_new, sizeof (BEAD) * Counts_new.BeadsCoor);
    Index_new = realloc(Index_new, sizeof *Index_new * Counts_new.BeadsCoor);
    // copy unbonded beads not to be exchanged to the start of bead_new //{{{
    // TODO: assumes unbonded beads are before bonded beads
    count = 0; // counts copied beads
    for (int i = 0; i < Counts_orig.Unbonded; i++) {
      // first, copy only beads of the type that's not to be exchange
      if (!bead_orig[i].Use) {
        bead_new[count] = bead_orig[i];
        bead_new[count].Molecule = -1;
        bead_new[count].Index = count;
        bead_new[count].Use = false; // do not rewrite, obviously
        bead_new[count].Aggregate = malloc(sizeof *bead_new[count].Aggregate *
                                           1); // just to free later
        Index_new[count] = count;
        count++;
      }
    }
    // count ended at <number of unbonded original beads> - <added beads> //}}}
    // put unbonded beads to be added beyond the unchanged unbonded beads //{{{
    count = Counts_orig.Unbonded - Counts_add.BeadsCoor; // just to be sure
    for (int i = 0; i < Counts_add.Unbonded; i++) {
      int type = bead_add[i].Type;
      int new_type = FindBeadType(bt_add[type].Name, Counts_new, bt_new);
      bead_new[count] = bead_orig[i];
      bead_new[count].Type = new_type;
      bead_new[count].Molecule = -1;
      bead_new[count].Index = count;
      bead_new[count].Use = true; // coordinates to be rewritten
      bead_new[count].Aggregate = malloc(sizeof *bead_new[count].Aggregate *
                                         1); // just to free later
      Index_new[count] = count;
      count++; // use count to make it consistent & easy to read
    } //}}}
    // copy the original bonded beads //{{{
    count = Counts_new.Unbonded;
    for (int i = Counts_orig.Unbonded; i < Counts_orig.BeadsCoor; i++) {
      bead_new[count] = bead_orig[i];
      bead_new[count].Index = count;
      bead_new[count].Use = false; // coordinates to be rewritten
      bead_new[count].Aggregate = malloc(sizeof *bead_new[count].Aggregate *
                                         1); // just to free later
      Index_new[count] = count;
      count++;
    } //}}}
    // put bonded beads to be added at the very end //{{{
    count = Counts_new.Unbonded + Counts_orig.Bonded;
    for (int i = Counts_add.Unbonded; i < Counts_add.BeadsCoor; i++) {
      int type = bead_add[i].Type;
      int new_type = FindBeadType(bt_add[type].Name, Counts_new, bt_new);
      bead_new[count] = bead_add[i];
      bead_new[count].Type = new_type;
      bead_new[count].Molecule = bead_add[i].Molecule + Counts_orig.Molecules;
      bead_new[count].Index = count;
      bead_new[count].Use = true; // coordinates to be rewritten
      bead_new[count].Aggregate = malloc(sizeof *bead_new[count].Aggregate *
                                         1); // just to free later
      Index_new[count] = count;
      count++; // use count to make it consistent & easy to read
    } //}}}
    // alocate new molecule struct
    mol_new = realloc(mol_new, sizeof (MOLECULE) * Counts_new.Molecules);
    // copy original molecules to _new struct //{{{
    for (int i = 0; i < Counts_orig.Molecules; i++) {
      int type = mol_orig[i].Type;
      mol_new[i].Type = type;
      mol_new[i].Bead = malloc(sizeof *mol_new[i].Bead * mt_new[type].nBeads);
      for (int j = 0; j < mt_new[type].nBeads; j++) {
        mol_new[i].Bead[j] = mol_orig[i].Bead[j] - Counts_add.Bonded;
      }
    } //}}}
    // put _add molecules into _new struct //{{{
    count = Counts_new.BeadsCoor - Counts_add.Bonded;
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int add_type = mol_add[i].Type;
      int new_type = FindMoleculeType(mt_add[add_type].Name,
                                      Counts_new, mt_new);
      int new_i = Counts_orig.Molecules + i;
      mol_new[new_i].Type = new_type;
      mol_new[new_i].Bead = malloc(sizeof *mol_new[new_i].Bead *
                                   mt_new[new_type].nBeads);
      for (int j = 0; j < mt_new[new_type].nBeads; j++) {
        mol_new[new_i].Bead[j] = count;
        count++;
      }
    } //}}}
    //}}}
  } else { // or add beads to the system? //{{{
    Counts_new.BeadsCoor = Counts_orig.BeadsCoor + Counts_add.BeadsCoor;
    Counts_new.BeadsTotal = Counts_orig.BeadsTotal + Counts_add.BeadsTotal;
    Counts_new.Bonded = Counts_orig.Bonded + Counts_add.Bonded;
    Counts_new.Unbonded = Counts_new.BeadsCoor - Counts_new.Bonded;
    Counts_new.TypesOfBonds = Counts_add.TypesOfBonds;
    Counts_new.TypesOfAngles = Counts_add.TypesOfAngles;
    Counts_new.Molecules = Counts_orig.Molecules + Counts_add.Molecules;
    // fill BeadType struct for the new system //{{{
    Counts_new.TypesOfBeads = Counts_orig.TypesOfBeads;
    // 1) copy original BeadType (if there is an input system)
    if (Counts_new.TypesOfBeads > 0) {
      CopyBeadType(Counts_new.TypesOfBeads, &bt_new, bt_orig, 3);
    }
    // 2) add new bead types - check is based only on Name //{{{
    for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
      bool new = true;
      for (int j = 0; j < Counts_orig.TypesOfBeads; j++) {
        if (strcmp(bt_add[i].Name, bt_orig[j].Name) == 0) {
          new = false;
          // increase old type's number of beads
          bt_new[j].Number += bt_add[i].Number;
          break;
        }
      }
      if (new) { // create new type
        int type = Counts_new.TypesOfBeads;
        bt_new = realloc(bt_new, sizeof (BEADTYPE) * (type + 1));
        bt_new[type] = bt_add[i];
        Counts_new.TypesOfBeads++;
      }
    } //}}}
    //}}}
    // fill MoleculeType struct for the new system //{{{
    Counts_new.TypesOfMolecules = Counts_orig.TypesOfMolecules;
    if (Counts_new.TypesOfMolecules > 0) { // are there input coordinates?
      mt_new = realloc(mt_new, sizeof (MOLECULETYPE) * Counts_new.TypesOfMolecules);
      // copy original MoleculeType to _new //{{{
      for (int i = 0; i < Counts_new.TypesOfMolecules; i++) {
        mt_new[i] = mt_orig[i]; // copy simple variables
        // copy Bead array
        mt_new[i].Bead = malloc(sizeof *mt_new[i].Bead * mt_new[i].nBeads);
        for (int j = 0; j < mt_new[i].nBeads; j++) {
          mt_new[i].Bead[j] = mt_orig[i].Bead[j];
        }
        // copy Bond array
        mt_new[i].Bond = malloc(sizeof *mt_new[i].Bond * mt_new[i].nBonds);
        for (int j = 0; j < mt_new[i].nBonds; j++) {
          mt_new[i].Bond[j][0] = mt_orig[i].Bond[j][0];
          mt_new[i].Bond[j][1] = mt_orig[i].Bond[j][1];
          mt_new[i].Bond[j][2] = mt_orig[i].Bond[j][2];
        }
      } //}}}
    }
    // add new molecule types - check if their the same based only on Name //{{{
    for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
      bool new = true;
      for (int j = 0; j < Counts_orig.TypesOfMolecules; j++) {
        if (strcmp(mt_add[i].Name, mt_orig[j].Name) == 0) {
          new = false;
          mt_new[j].Number += mt_add[i].Number;
          break;
        }
      }
      if (new) {
        int type = Counts_new.TypesOfMolecules;
        mt_new = realloc(mt_new, sizeof (MOLECULETYPE) * (type + 1));
        mt_new[type] = mt_add[i];
        mt_new[type].Bead = malloc(sizeof *mt_new[type].Bead *
                                   mt_new[type].nBeads);
        for (int j = 0; j < mt_new[type].nBeads; j++) {
          int old_type = mt_add[i].Bead[j];
          int btype = FindBeadType(bt_add[old_type].Name,
                                   Counts_new, bt_new);
          mt_new[type].Bead[j] = btype;
        }
        if (mt_new[type].nBonds > 0) {
          mt_new[type].Bond = malloc(sizeof *mt_new[type].Bond *
                                     mt_new[type].nBonds);
          for (int j = 0; j < mt_new[type].nBonds; j++) {
            mt_new[type].Bond[j][0] = mt_add[i].Bond[j][0];
            mt_new[type].Bond[j][1] = mt_add[i].Bond[j][1];
            mt_new[type].Bond[j][2] = mt_add[i].Bond[j][2];
          }
        }
        if (mt_new[type].nAngles > 0) {
          mt_new[type].Angle = malloc(sizeof *mt_new[type].Angle *
                                      mt_new[type].nAngles);
          for (int j = 0; j < mt_new[type].nAngles; j++) {
            mt_new[type].Angle[j][0] = mt_add[i].Angle[j][0];
            mt_new[type].Angle[j][1] = mt_add[i].Angle[j][1];
            mt_new[type].Angle[j][2] = mt_add[i].Angle[j][2];
            mt_new[type].Angle[j][3] = mt_add[i].Angle[j][3];
          }
        }
        if (mt_new[type].nDihedrals > 0) {
          mt_new[type].Dihedral = malloc(sizeof *mt_new[type].Dihedral *
                                      mt_new[type].nDihedrals);
          for (int j = 0; j < mt_new[type].nDihedrals; j++) {
            mt_new[type].Dihedral[j][0] = mt_add[i].Dihedral[j][0];
            mt_new[type].Dihedral[j][1] = mt_add[i].Dihedral[j][1];
            mt_new[type].Dihedral[j][2] = mt_add[i].Dihedral[j][2];
            mt_new[type].Dihedral[j][3] = mt_add[i].Dihedral[j][3];
            mt_new[type].Dihedral[j][4] = mt_add[i].Dihedral[j][4];
          }
        }
        Counts_new.TypesOfMolecules++;
      }
    } //}}}
    //}}}
    // fill Bead struct for the new system
    bead_new = realloc(bead_new, sizeof (BEAD) * Counts_new.BeadsCoor);
    Index_new = realloc(Index_new, sizeof *Index_new * Counts_new.BeadsCoor);
    // copy original unbonded beads to the start of bead_new //{{{
    // TODO: assumes unbonded beads are before bonded beads
    for (int i = 0; i < Counts_orig.Unbonded; i++) {
      bead_new[i] = bead_orig[i];
      bead_new[i].Molecule = -1;
      bead_new[i].Index = i;
      bead_new[i].Use = false; // do not rewrite, obviously
      bead_new[i].Aggregate = malloc(sizeof *bead_new[i].Aggregate * 1);
      Index_new[i] = i;
    } //}}}
    // put unbonded beads to be added beyond the original unbonded beads //{{{
    for (int i = Counts_orig.Unbonded; i < Counts_new.Unbonded; i++) {
      int id_add = i - Counts_orig.Unbonded;
      int type = bead_add[id_add].Type;
      int new_type = FindBeadType(bt_add[type].Name, Counts_new, bt_new);
      bead_new[i] = bead_add[id_add];
      bead_new[i].Type = new_type;
      bead_new[i].Molecule = -1;
      bead_new[i].Index = i;
      bead_new[i].Aggregate = malloc(sizeof *bead_new[i].Aggregate * 1);
      Index_new[i] = i;
    } //}}}
    // copy the original bonded beads //{{{
    for (int i = Counts_orig.Unbonded; i < Counts_orig.BeadsCoor; i++) {
      // id goes from Counts_new.Unbonded to (Counts_new.Unbonded+Counts_orig.Bonded)
      int id = Counts_new.Unbonded + i - Counts_orig.Unbonded;
      bead_new[id] = bead_orig[i];
      bead_new[id].Index = id;
      bead_new[id].Aggregate = malloc(sizeof *bead_new[id].Aggregate * 1);
      Index_new[id] = id;
    } //}}}
    // put bonded beads to be added at the very end //{{{
    for (int i = Counts_add.Unbonded; i < Counts_add.BeadsCoor; i++) {
      int type = bead_add[i].Type;
      int new_type = FindBeadType(bt_add[type].Name, Counts_new, bt_new);
      int id = Counts_new.BeadsCoor - Counts_add.BeadsCoor + i;
      bead_new[id] = bead_add[i];
      bead_new[id].Type = new_type;
      bead_new[id].Molecule = bead_add[i].Molecule + Counts_orig.Molecules;
      bead_new[id].Index = id;
      bead_new[id].Aggregate = malloc(sizeof *bead_new[id].Aggregate * 1);
      Index_new[id] = id;
    } //}}}
    // alocate new molecule struct
    mol_new = realloc(mol_new, sizeof (MOLECULE) * Counts_new.Molecules);
    // copy original molecules to _new struct //{{{
    for (int i = 0; i < Counts_orig.Molecules; i++) {
      int type = mol_orig[i].Type;
      mol_new[i].Type = type;
      mol_new[i].Bead = malloc(sizeof *mol_new[i].Bead * mt_new[type].nBeads);
      for (int j = 0; j < mt_new[type].nBeads; j++) {
        mol_new[i].Bead[j] = mol_orig[i].Bead[j] + Counts_add.Unbonded;
      }
    } //}}}
    // put _add molecules into _new struct //{{{
    count = Counts_new.BeadsCoor - Counts_add.Bonded;
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int add_type = mol_add[i].Type;
      int new_type = FindMoleculeType(mt_add[add_type].Name,
                                      Counts_new, mt_new);
      int new_i = Counts_orig.Molecules + i;
      mol_new[new_i].Type = new_type;
      mol_new[new_i].Bead = malloc(sizeof *mol_new[new_i].Bead *
                                   mt_new[new_type].nBeads);
      for (int j = 0; j < mt_new[new_type].nBeads; j++) {
        mol_new[new_i].Bead[j] = count;
        count++;
      }
    } //}}}
  } //}}}
  FillMolBTypes(Counts_new.TypesOfMolecules, &mt_new);
  FillMolMassCharge(Counts_new.TypesOfMolecules, &mt_new, bt_new);

//PrintCounts(Counts_new);
//PrintMoleculeType2(Counts_new.TypesOfMolecules, bt_new, mt_new);
//PrintMolecule(Counts_new.Molecules, mt_new, mol_new, bt_new, bead_new);

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(Counts_new, bt_new, bead_new, mt_new, mol_new);
    PrintBondTypes2(Counts_new.TypesOfBonds, bond_type);
  } //}}}

  // add beads randomly if FIELD-like file is used //{{{
  double dist;
  if (strlen(add_vsf) == 0) {
    count = 0;
    // add monomeric beads //{{{
    for (int i = 0; i < Counts_add.Unbonded; i++) {
      VECTOR random;
      if (lowest_dist != -1 || highest_dist != -1) {
        double min_dist;
        int tries = 0;
        do {
          tries++;
          if (tries == 1000000) {
            ColourText(STDERR_FILENO, YELLOW);
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

          min_dist = SQR(Box_orig.Length.x * 100);
          for (int j = 0; j < Counts_orig.BeadsCoor; j++) {
            int btype = bead_orig[j].Type;
            /*
             * j can be added monomeric bead, so it's type can be higher than
             * the number of types
             */
            if (btype < Counts_orig.TypesOfBeads && bt_orig[btype].Use) {
              VECTOR dist;
              dist = Distance(bead_orig[j].Position, random, Box_new.Length);
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
        id = Counts_orig.Unbonded + i;
      } else { // switched beds (--switch option)
        for (int j = count; j < Counts_new.Unbonded; j++) {
          if (bead_new[j].Use) { // is this an original bead to be exchanged?
            id = j;
            bead_new[j].Use = false; // just exchanged (only pro forma)
            count = j + 1;
            break;
          }
        }
      }
      if (id == -1) {
        ColourText(STDERR_FILENO, RED);
        fprintf(stderr, "!!!SOME ERROR!!!");
        fprintf(stderr, "...very useful.");
        ColourReset(STDERR_FILENO);
        exit(1);
      }

      // add the new coordinate
      bead_new[id].Position.x = random.x;
      bead_new[id].Position.y = random.y;
      bead_new[id].Position.z = random.z;

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
      fprintf(stdout, "\rMonomer placed: %d\n", Counts_add.Unbonded);
    } //}}}
    // add molecules //{{{
    // doesn't depend on --switch option as it's determined by the mol_new
    // array established earlier
    count = 0;
    for (int i = Counts_orig.Molecules; i < Counts_new.Molecules; i++) {
      int mtype = mol_new[i].Type;

      VECTOR rotated[mt_new[mtype].nBeads];
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
        for (int j = 0; j < mt_new[mtype].nBeads; j++) {
          int id = mol_new[i].Bead[j];
          rotated[j].x = rot.x.x * bead_new[id].Position.x
                       + rot.x.y * bead_new[id].Position.y
                       + rot.x.z * bead_new[id].Position.z;
          rotated[j].y = rot.y.x * bead_new[id].Position.x
                       + rot.y.y * bead_new[id].Position.y
                       + rot.y.z * bead_new[id].Position.z;
          rotated[j].z = rot.z.x * bead_new[id].Position.x
                       + rot.z.y * bead_new[id].Position.y
                       + rot.z.z * bead_new[id].Position.z;
        }
      } else { // don't rotate
        for (int j = 0; j < mt_new[mtype].nBeads; j++) {
          int id = mol_new[i].Bead[j];
          rotated[j].x = bead_new[id].Position.x;
          rotated[j].y = bead_new[id].Position.y;
          rotated[j].z = bead_new[id].Position.z;
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
            ColourText(STDERR_FILENO, YELLOW);
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

          min_dist = SQR(Box_new.Length.x) +
                     SQR(Box_new.Length.y) +
                     SQR(Box_new.Length.z);
          for (int j = 0; j < Counts_orig.BeadsCoor; j++) {
            int btype_j = bead_orig[j].Type;
            /*
             * j can be added monomeric bead, so it's type can be higher than
             * the number of types
             */
            if (btype_j < Counts_orig.TypesOfBeads && bt_orig[btype_j].Use) {
              dist = Length(Distance(bead_orig[j].Position,
                                     random, Box_orig.Length));
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
      for (int j = 0; j < mt_new[mtype].nBeads; j++) {
        int id = mol_new[i].Bead[j];
        bead_new[id].Position.x = random.x + rotated[j].x;
        bead_new[id].Position.y = random.y + rotated[j].y;
        bead_new[id].Position.z = random.z + rotated[j].z;
      } //}}}

      // print number of placed molecules? //{{{
      if (!silent && isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\rMolecules placed: %d", i-Counts_orig.Molecules+1);
      } //}}}
    } //}}}
    // print total number of placed molecules? //{{{
    if (!silent) {
      if (isatty(STDOUT_FILENO)) {
        fflush(stdout);
        fprintf(stdout, "\r                                             \r");
      }
      fprintf(stdout, "Molecules placed: %3d\n", Counts_add.Molecules);
    } //}}}
  } //}}}

  // open output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  // print command, bead type names & box size to output .vcf file //{{{
  fprintf(out, "# Generated by:");
  PrintCommand(out, argc, argv);
  fprintf(out, "# AnalysisTools version %s;", VERSION);
  fprintf(out, " https://github.com/KaGaSi/AnalysisTools/releases\n");

  fprintf(out, "\npbc %lf %lf %lf\n", Box_new.Length.x,
                                      Box_new.Length.y,
                                      Box_new.Length.z); //}}}

  // print coordinates to output .vcf file //{{{
  // write all beads (Write flag was used with '-xb' option)
  for (int i = 0; i < Counts_new.TypesOfBeads; i++) {
    bt_new[i].Write = true;
  }
  // write all molecules (basically just to make sure)
  for (int i = 0; i < Counts_new.TypesOfMolecules; i++) {
    mt_new[i].Write = true;
  }
  WriteCoorIndexed(out, Counts_new, bt_new, bead_new,
                   mt_new, mol_new, stuff, Box_new); //}}}

  // print coordinates to xyz file (if -xyz option is present) //{{{
  FILE *xyz = NULL;
  if (strlen(output_xyz) > 0) {
    // open output .xyz file for reading
    if ((xyz = fopen(output_xyz, "w")) == NULL) {
      ErrorFileOpen(output_xyz, 'w');
      exit(1);
    }
    WriteCoorXYZ(xyz, Counts_new, bt_new, bead_new);
  } //}}}

  // close output files //{{{
  fclose(out);
  if (output_xyz[0] != '\0') {
    fclose(xyz);
  } //}}}

  // create output vsf file
  WriteVsf(output_vsf, Counts_new, bt_new, bead_new, mt_new, mol_new, false);

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts_orig, &mt_orig, &mol_orig,
                 &bt_orig, &bead_orig, &Index_orig);
  FreeSystemInfo(Counts_add, &mt_add, &mol_add,
                 &bt_add, &bead_add, &Index_add);
  FreeSystemInfo(Counts_new, &mt_new, &mol_new,
                 &bt_new, &bead_new, &Index_new);
  free(stuff);
  free(bond_type);
  free(angle_type);
  free(dihedral_type);
  //}}}

  return 0;
}
