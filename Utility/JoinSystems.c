#include "../AnalysisTools.h"

// TODO: final box size is somewhat hard

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
JoinSystems connects two coordinate files, creating a new systems consisting \
of both systems.\n\n");
  }

  fprintf(ptr, "Usage: %s <input1> <input2> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input1>            first input coordinate file\n");
  fprintf(ptr, "<input2>            second input coordinate file\n");
  fprintf(ptr, "<output>            output coordinate file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -i1 <filename>    structure file for <input1>\n");
  fprintf(ptr, "  -i2 <filename>    structure file for <input2>\n");
  fprintf(ptr, "  -o <filename>     output structure file\n");
  fprintf(ptr, "  -offset 3x<float> how much offset the second system"
               " against to the first\n");
  fprintf(ptr, "  --centre          centre the second system on the first"
               "(overrides -offset option)\n");
  fprintf(ptr, "  -box 3x<float>    output box dimensions (orthogonal)\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 4, all = common + 6, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-i1");
  strcpy(option[count++], "-i2");
  strcpy(option[count++], "-o");
  strcpy(option[count++], "-offset");
  strcpy(option[count++], "--centre");
  strcpy(option[count++], "-box");
  if (count != all) {
    strcpy(ERROR_MSG, "Debug - numbers of arguments");
    PrintError();
  }
  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // <input1> - first input coordinate file //{{{
  char coor_file_1[LINE] = "", struct_file_1[LINE] = "";
  int coor_type_1, struct_type_1 = 0;
  snprintf(coor_file_1, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file_1, &coor_type_1,
                       struct_file_1, &struct_type_1)) {
    exit(1);
  } //}}}
  // <input2> - first input coordinate file //{{{
  char coor_file_2[LINE] = "", struct_file_2[LINE] = "";
  int coor_type_2, struct_type_2 = 0;
  snprintf(coor_file_2, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file_2, &coor_type_2,
                       struct_file_2, &struct_type_2)) {
    exit(1);
  } //}}}
  // <output> - output coordinate file //{{{
  char coor_out_file[LINE] = "";
  snprintf(coor_out_file, LINE, "%s", argv[++count]);
  int coor_out_type = CoordinateFileType(coor_out_file, 1); //}}}

  // options before reading system data //{{{
  // output structure file (-o option) //{{{
  char struct_file_out[LINE] = "";
  int struct_type_out = -1;
  if (FileOption(argc, argv, "-o", struct_file_out)) {
    exit(1);
  }
  if (struct_file_out[0] != '\0') {
    struct_type_out = StructureFileType(struct_file_out, 1);
  } //}}}
  // input structure files (-i1/-i2 options) //{{{
  if (FileOption(argc, argv, "-i1", struct_file_1)) {
    exit(1);
  }
  if (struct_file_1[0] != '\0') {
    struct_type_1 = StructureFileType(struct_file_1, 0);
  }
  if (FileOption(argc, argv, "-i2", struct_file_2)) {
    exit(1);
  }
  if (struct_file_2[0] != '\0') {
    struct_type_2 = StructureFileType(struct_file_2, 0);
  } //}}}
  // TODO: --detailed1 & --detailed2 options
  //       --pbc_xyz1 & --pbc_xyz2 options
  //       -st1 & -st2 options
  bool silent, verbose, detailed;
  int start = 1,
      trash, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &start, &trash, &trash);
  // // -offset option //{{{
  // double temp[3] = {0};
  // VECTOR offset = {0, 0, 0};
  // // bool off = false;
  // if (DoubleOption(argc, argv, 3, "-offset", &count, temp)) {
  //   if (count != 3) {
  //     strcpy(ERROR_MSG, "three numbers required");
  //     PrintErrorOption("-offset");
  //     exit(1);
  //   }
  //   // off = true;
  //   offset.x = temp[0];
  //   offset.y = temp[1];
  //   offset.z = temp[2];
  // } //}}}
  // -offset option //{{{
  VECTOR offset = {-11111, -11111, -11111};
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-offset") == 0) {
      if (argc < (i + 3) ||
          (argv[i+1][0] != 'c' && !IsRealNumber(argv[i+1], &offset.x)) ||
          (argv[i+2][0] != 'c' && !IsRealNumber(argv[i+2], &offset.y)) ||
          (argv[i+3][0] != 'c' && !IsRealNumber(argv[i+3], &offset.z))) {
        strcpy(ERROR_MSG, "wrong/missing arguments (either number or 'c')");
        PrintErrorOption("-offset");
        exit(1);
      }
      break;
    }
  } //}}}
  printf("%lf %lf %lf\n", offset.x, offset.y, offset.z);
  // output box dimensions //{{{
  VECTOR box_opt = {0, 0, 0};
  double temp[3] = {0};
  if (DoubleOption(argc, argv, 3, "-box", &count, temp)) {
    if (count != 3) {
      strcpy(ERROR_MSG, "three numbers required");
      PrintErrorOption("-box");
      exit(1);
    }
    box_opt.x = temp[0];
    box_opt.y = temp[1];
    box_opt.z = temp[2];
  } //}}}
  bool centre = BoolOption(argc, argv, "--centre");
  //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM Sys_1 = ReadStructure(struct_type_1, struct_file_1,
                               coor_type_1, coor_file_1, detailed, pbc_xyz);
  SYSTEM Sys_2 = ReadStructure(struct_type_2, struct_file_2,
                               coor_type_2, coor_file_2, detailed, pbc_xyz);

  // verbose output describing the two systems to be joined //{{{
  if (verbose) {
    printf("\n==================================================");
    printf("\nFirst sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys_1);
    printf("\n==================================================");
    printf("\nSecond sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys_2);
  } //}}}

  // read first coordinate file //{{{
  FILE *fr = OpenFile(coor_file_1, "r");
  int line_count = 0;     // count lines in the vcf file
  if (!ReadTimestep(coor_type_1, fr, coor_file_1, &Sys_1, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_1, "\0", "\0");
    exit(1);
  }
  fclose(fr); //}}}
  // read second coordinate file //{{{
  fr = OpenFile(coor_file_2, "r");
  if (!ReadTimestep(coor_type_2, fr, coor_file_2, &Sys_2, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_2, "\0", "\0");
    exit(1);
  }
  fclose(fr); //}}}

  // make proper offset if 'c' is used in -offset option //{{{
  if (offset.x == -11111) {
    double c1 = Sys_1.Box.Low.x + 0.5 * Sys_1.Box.Length.x,
           c2 = Sys_2.Box.Low.x + 0.5 * Sys_2.Box.Length.x;
    offset.x = c1 - c2;
  }
  if (offset.y == -11111) {
    double c1 = Sys_1.Box.Low.y + 0.5 * Sys_1.Box.Length.y,
           c2 = Sys_2.Box.Low.y + 0.5 * Sys_2.Box.Length.y;
    offset.y = c1 - c2;
  }
  if (offset.z == -11111) {
    double c1 = Sys_1.Box.Low.z + 0.5 * Sys_1.Box.Length.z,
           c2 = Sys_2.Box.Low.z + 0.5 * Sys_2.Box.Length.z;
    offset.z = c1 - c2;
  } //}}}

  // move the beads of the second system //{{{
  for (int i = 0; i < Sys_2.Count.Bead; i++) {
    int id = Sys_2.BeadCoor[i];
    Sys_2.Bead[id].Position.x += offset.x;
    Sys_2.Bead[id].Position.y += offset.y;
    Sys_2.Bead[id].Position.z += offset.z;
  }
  Sys_2.Box.Length.x += offset.x;
  Sys_2.Box.Length.y += offset.y;
  Sys_2.Box.Length.z += offset.z;
  CalculateBoxData(&Sys_2.Box, 0); //}}}

  // create output system
  SYSTEM Sys_out = CopySystem(Sys_1);
  // pick box size as the larger dimensions from the initial systems //{{{
  BOX box_out = InitBox,
      *box1 = &Sys_1.Box,
      *box2 = &Sys_2.Box;
  // ...assumes orthogonal box
  // set output Length to maxima of the two systems
  if (box1->Length.x > box2->Length.x) {
    box_out.Length.x = box1->Length.x;
  } else {
    box_out.Length.x = box2->Length.x;
  }
  if (box1->Length.y > box2->Length.y) {
    box_out.Length.y = box1->Length.y;
  } else {
    box_out.Length.y = box2->Length.y;
  }
  if (box1->Length.z > box2->Length.z) {
    box_out.Length.z = box1->Length.z;
  } else {
    box_out.Length.z = box2->Length.z;
  }
  // set out Low to the minima of the two systems
  if (box1->Low.x < box2->Low.x) {
    box_out.Low.x = box1->Low.x;
  } else {
    box_out.Low.x = box2->Low.x;
  }
  if (box1->Low.y < box2->Low.y) {
    box_out.Low.y = box1->Low.y;
  } else {
    box_out.Low.y = box2->Low.y;
  }
  if (box1->Low.z < box2->Low.z) {
    box_out.Low.z = box1->Low.z;
  } else {
    box_out.Low.z = box2->Low.z;
  }
  // assume orthogonal box
  box_out.alpha = 90;
  box_out.beta = 90;
  box_out.gamma = 90;
  CalculateBoxData(&box_out, 0); //}}}
  ConcatenateSystems(&Sys_out, Sys_2, box_out);
  // if -box option is present, use it as box size //{{{
  if (box_opt.x != 0) {
    // put the system into the box centre (-centre option)
    if (centre) {
      VECTOR c1;
      c1.x = 0.5 * box_opt.x;
      c1.y = 0.5 * box_opt.y;
      c1.z = 0.5 * box_opt.z;
      VECTOR c2;
      c2.x = Sys_out.Box.Low.x + 0.5 * Sys_out.Box.Length.x;
      c2.y = Sys_out.Box.Low.y + 0.5 * Sys_out.Box.Length.y;
      c2.z = Sys_out.Box.Low.z + 0.5 * Sys_out.Box.Length.z;
      VECTOR delta;
      delta.x = c1.x - c2.x;
      delta.y = c1.y - c2.y;
      delta.z = c1.z - c2.z;
      for (int i = 0; i < Sys_out.Count.Bead; i++) {
        int id = Sys_out.BeadCoor[i];
        Sys_out.Bead[id].Position.x += delta.x;
        Sys_out.Bead[id].Position.y += delta.y;
        Sys_out.Bead[id].Position.z += delta.z;
      }
    }
    Sys_out.Box.Length = box_opt;
    Sys_out.Box.Low.x = 0;
    Sys_out.Box.Low.y = 0;
    Sys_out.Box.Low.z = 0;
    CalculateBoxData(&Sys_out.Box, 0);
  } //}}}

  // verbose output describing the output system //{{{
  if (verbose) {
    printf("\n==================================================");
    printf("\nNew sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys_out);
  } //}}}

  // print initial stuff to output coordinate file //{{{
  if (coor_out_type == VCF_FILE) {
    PrintByline(coor_out_file, argc, argv);
  } else if (coor_out_type == VTF_FILE) {
    PrintByline(coor_out_file, argc, argv);
    WriteStructure(VSF_FILE, coor_out_file, Sys_out, -1, false);
    coor_out_type = VCF_FILE;
  } else {
    FILE *out = OpenFile(coor_out_file, "w");
    fclose(out);
  } //}}}

  SubtractLow(&Sys_out); // make coordinates from 0 to Box.Length
  bool *write = malloc(Sys_out.Count.Bead * sizeof *write);
  InitBoolArray(write, Sys_out.Count.Bead, true);
  WriteTimestep(coor_out_type, coor_out_file, Sys_out, 1, write);
  if (struct_file_out[0] != '\0') {
    WriteStructure(struct_type_out, struct_file_out, Sys_out, -1, false);
  }

  // free memory
  FreeSystem(&Sys_1);
  FreeSystem(&Sys_2);
  FreeSystem(&Sys_out);
  free(write);

  return 0;
}
