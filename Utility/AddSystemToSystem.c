#include "../AnalysisTools.h"
// TODO: 'option' to create vtf instead of vsf+vcf
// TODO: 'option' to read vtf instead of vsf+vcf (to-be-added system)
// TODO: add -rot <angle> option
// TODO: don't forget - new system molecules' ids from the to-be-added system
//       are S_add.Molecule[].Index + S_orig.Count.Molecule (I think)

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
AddSystemToSystem adds new beads with predefined coordinates into the original \
system. Therefore the new system must be fully defined (to randomly add \
monmeric beads and/or molecules, use AddRandomToSystem utility). \
The added system can be offset or rotated with regards to the original one, \
and the dimension of the new system can be defined separately.\n\n");
  }
  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.vcf> <in.vsf> <in.vcf> \
<out.vsf> <out.vcf> [options]\n\n", cmd);

  fprintf(ptr, "      <input>        input original coordinate file \
(vcf or vtf format)\n");
  fprintf(ptr, "      <in.vsf>       input structure file for added system \
(vsf or vtf format)\n");
  fprintf(ptr, "      <in.vcf>       input coordinate file for added system \
(vcf or vtf format)\n");
  fprintf(ptr, "      <out.vsf>      output structure file (vsf format)\n");
  fprintf(ptr, "      <out.vcf>      output coordinate file (vcf format)\n");
  fprintf(ptr, "   [general options]\n");
  fprintf(ptr, "      -st_orig <int>       timestep to use from the original \
system (default: 1)\n");
  fprintf(ptr, "      -st_add <int>        timestep to use from the \
to-be-added system (default: 1)\n");
  fprintf(ptr, "      -b <x> <y> <z>       size of the new cuboid box\n");
  fprintf(ptr, "      -ba <a> <b> <g>      angles of the new triclinic box\n");
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
  int req_args = 5; //}}}

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
        strcmp(argv[i], "-offset") != 0 &&
        (argv[i][1] < '0' || argv[i][1] > '9') &&
        strcmp(argv[i], "-st_orig") != 0 &&
        strcmp(argv[i], "-st_add") != 0 &&
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-ba") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char file_in_coor[LINE] = "",
       file_in_struct[LINE] = "";
  snprintf(file_in_coor, LINE, "%s", argv[++count]);
  bool vtf = false;
  if (!InputCoor_old(&vtf, file_in_coor, file_in_struct)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <in.vsf> & <in.vcf> - system to be added //{{{
  char file_add_struct[LINE] = "",
       file_add_coor[LINE] = "";
  // <in.vsf>
  snprintf(file_add_struct, LINE, "%s", argv[++count]);
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(file_add_struct, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  }
  // <in.vcf>
  strcpy(extension[0], ".vcf");
  snprintf(file_add_coor, LINE, "%s", argv[++count]);
  if (ErrorExtension(file_add_coor, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> & <out.vsf> - new system //{{{
  char file_out_struct[LINE] = "",
       file_out_coor[LINE] = "";
  // <out.vsf>
  snprintf(file_out_struct, LINE, "%s", argv[++count]);
  ext = 1;
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(file_out_struct, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  }
  // <out.vcf>
  snprintf(file_out_coor, LINE, "%s", argv[++count]);
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(file_out_coor, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed;
  CommonOptions_old(argc, argv, file_in_struct, LINE, &verbose, &silent, &detailed);

  // -offset <x> <y> <z> define offset for -vtf file //{{{
  double offset[100] = {0};
  if (MultiDoubleOption(argc, argv, "-offset", &count, offset)) {
    exit(1);
  }
  if (count != 0 && count != 3) {
    strcpy(ERROR_MSG, "three numbers required");
    PrintErrorOption("-offset");
    Help(argv[0], true);
    exit(1);
  } //}}}

  // starting timestep(s) //{{{
  int st_orig = 1, st_add = 1;
  if (IntegerOption(argc, argv, "-st_orig", &st_orig)) {
    exit(1);
  }
  if (IntegerOption(argc, argv, "-st_add", &st_add)) {
    exit(1);
  } //}}}

  // define new box size //{{{
  double box_option[100] = {-1};
  if (MultiDoubleOption(argc, argv, "-b", &count, box_option)) {
    exit(1);
  }
  if (count != 0 && count != 3) {
    strcpy(ERROR_MSG, "three non-negative numbers required");
    PrintErrorOption("-b");
    Help(argv[0], true);
    exit(1);
  }
  double box_angle_option[100] = {-1};
  if (MultiDoubleOption(argc, argv, "-ba", &count, box_angle_option)) {
    exit(1);
  }
  if (count != 0 && count != 3) {
    strcpy(ERROR_MSG, "three non-negative numbers required");
    PrintErrorOption("-ba");
    Help(argv[0], true);
    exit(1);
  }
  //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // array for the timestep preamble
  char *stuff = calloc(LINE, sizeof *stuff);

  // read information about the original system //{{{
  SYSTEM S_orig = VtfReadStruct(file_in_struct, detailed);
  FILE *vcf;
  vcf = OpenFile(file_in_coor, "r");
  int file_line_count = 0, count_vcf = 0;
  while (count_vcf < st_orig) {
    count_vcf++;
    if (!VtfReadTimestep(vcf, file_in_coor, file_in_struct, &S_orig,
                         &file_line_count, stuff)) {
      count_vcf--;
      break;
    }
  }
  fclose(vcf);
  if (count_vcf < st_orig) {
    strcpy(ERROR_MSG, "starting step for the original system is too high");
    PrintWarning();
    WarnPrintFile(file_in_coor, file_in_struct, "\0");
    fprintf(stderr, "%s contains %s%d%s steps (using the last one)%s\n",
            ErrCyan(), ErrYellow(), count_vcf, ErrCyan(), ErrColourReset());
  } //}}}

  // read information about the to-be-added system //{{{
  SYSTEM S_add = VtfReadStruct(file_add_struct, false);
  vcf = OpenFile(file_add_coor, "r");
  file_line_count = 0;
  count_vcf = 0;
  while (count_vcf < st_add) {
    count_vcf++;
    if (!VtfReadTimestep(vcf, file_add_coor, file_add_struct, &S_add,
                         &file_line_count, stuff)) {
      count_vcf--;
      break;
    }
  }
  fclose(vcf);
  if (count_vcf < st_add) {
    strcpy(ERROR_MSG, "starting step for the to-be-added system is too high");
    PrintWarning();
    WarnPrintFile(file_add_coor, file_add_struct, "\0");
    fprintf(stderr, "%s contains %s%d%s steps (using the last one)%s\n",
            ErrCyan(), ErrYellow(), count_vcf, ErrCyan(), ErrColourReset());
  }
  ToFractionalCoor(&S_add);
  for (int i = 0; i < S_add.Count.BeadCoor; i++) {
    int id = S_add.BeadCoor[i];
    S_add.Bead[id].Position.x += offset[0];
    S_add.Bead[id].Position.y += offset[1];
    S_add.Bead[id].Position.z += offset[2];
  }
  FromFractionalCoor(&S_add); //}}}

#if 0 // there will be rotation one day //{{{
  VECTOR rotated[S_add.Count.BeadCoor];
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
      // transform the added system (rotation matrix * coordinates)
      ToFractionalCoor(&S_add);
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
      FromFractionalCoor(&S_add);
    }
#endif //}}}

  // set new box size as either the original one or via -b & -ba options //{{{
  BOX Box_new = S_orig.Box;
  if (box_option[0] != -1) {
    Box_new.Length.x = box_option[0];
    Box_new.Length.y = box_option[1];
    Box_new.Length.z = box_option[2];
  }
  if (box_angle_option[0] != -1) {
    Box_new.alpha = box_angle_option[0];
    Box_new.beta = box_angle_option[1];
    Box_new.gamma = box_angle_option[2];
  }
  // error - no box size (should never trigger) //{{{
  if (Box_new.Length.x == 0 || Box_new.Length.y == 0 || Box_new.Length.z == 0) {
    strcpy(ERROR_MSG, "zero box size for the new system");
    PrintError();
    exit(1);
  } //}}}
  //}}}

  // join original and added systems
  SYSTEM S_new = CopySystem(S_orig);
  ConcatenateSystems(&S_new, S_add, Box_new);
//PruneSystem(&S_new);

  // verbose output //{{{
  if (verbose) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(S_orig);
    if (st_orig > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", st_orig);
    }
    fprintf(stdout, "\nSYSTEM TO ADD\n");
    VerboseOutput(S_add);
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(S_new);
  } //}}}

  // write data to output files //{{{
  // vsf file
  VtfWriteStruct(file_out_struct, S_new, -1);
  // .vcf file
  FILE *out = OpenFile(file_out_coor, "w");
  PrintByline(out, argc, argv);
  bool *write = malloc(sizeof *write * S_new.Count.Bead);
  for (int i = 0; i < S_new.Count.Bead; i++) {
    write[i] = true;
  }
  VtfWriteCoorIndexed(out, stuff, write, S_new);
  fclose(out);
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&S_orig);
  FreeSystem(&S_add);
  FreeSystem(&S_new);
  free(stuff);
  free(write);
  //}}}

  return 0;
}
