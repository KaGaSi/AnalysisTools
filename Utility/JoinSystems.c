#include "../AnalysisTools.h"

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
  fprintf(ptr, "  -offset 3x<float>|c\n"
          "                    offset of the second system against "
          "the first ('c' to place it in the centre of the first system)\n");
  fprintf(ptr, "  --centre          centre the second system on the first "
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
  // -offset option //{{{
  VECTOR offset = {0, 0, 0};
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
      if (argv[i+1][0] == 'c') {
        offset.x = -11111;
      }
      if (argv[i+2][0] == 'c') {
        offset.y = -11111;
      }
      if (argv[i+3][0] == 'c') {
        offset.z = -11111;
      }
      break;
    }
  } //}}}
  printf("%lf %lf %lf\n", offset.x, offset.y, offset.z);
  // output box dimensions //{{{
  VECTOR box_opt = {0, 0, 0};
  double temp[3] = {0};
  if (DoubleOption(argc, argv, 3, "-box", &count, temp)) {
    box_opt.x = temp[0];
    box_opt.y = temp[1];
    box_opt.z = temp[2];
    if (count != 3 || box_opt.x <= 0 || box_opt.y <= 0 || box_opt.z <= 0) {
      strcpy(ERROR_MSG, "three positive numbers required");
      PrintErrorOption("-box");
      Help(argv[0], true, common, option);
      exit(1);
    }
  } //}}}
  bool centre = BoolOption(argc, argv, "--centre");
  if (centre && box_opt.x == 0) {
    strcpy(ERROR_MSG, "no effect when -box option is not used");
    PrintWarnOption("--centre");
  }
  //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM Sys_1 = ReadStructure(struct_type_1, struct_file_1,
                               coor_type_1, coor_file_1, detailed, pbc_xyz);
  SYSTEM Sys_2 = ReadStructure(struct_type_2, struct_file_2,
                               coor_type_2, coor_file_2, detailed, pbc_xyz);
  BOX *box1 = &Sys_1.Box,
      *box2 = &Sys_2.Box;

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
  fclose(fr);
  AddLow(&Sys_1); //}}}
  // read second coordinate file //{{{
  fr = OpenFile(coor_file_2, "r");
  if (!ReadTimestep(coor_type_2, fr, coor_file_2, &Sys_2, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_2, "\0", "\0");
    exit(1);
  }
  fclose(fr);
  AddLow(&Sys_2); //}}}

  // make proper offset vector //{{{
  if (offset.x == -11111) { // a)
    offset.x = (box1->Low.x + 0.5 * box1->Length.x) -
               (box2->Low.x + 0.5 * box2->Length.x);
  }
  if (offset.y == -11111) { // a)
    offset.y = (box1->Low.y + 0.5 * box1->Length.y) -
               (box2->Low.y + 0.5 * box2->Length.y);
  }
  if (offset.z == -11111) { // a)
    offset.z = (box1->Low.z + 0.5 * box1->Length.z) -
               (box2->Low.z + 0.5 * box2->Length.z);
  } //}}}

  printf("offset: %lf %lf %lf\n", offset.x, offset.y, offset.z);
  // move the beads of the second system //{{{
  for (int i = 0; i < Sys_2.Count.Bead; i++) {
    int id = Sys_2.BeadCoor[i];
    Sys_2.Bead[id].Position.x += offset.x;
    Sys_2.Bead[id].Position.y += offset.y;
    Sys_2.Bead[id].Position.z += offset.z;
  } //}}}

  // create output system
  SYSTEM Sys_out = CopySystem(Sys_1);
  // pick box size as the larger dimensions from the initial systems //{{{
  BOX box_out = InitBox;
  printf("BOX 1: \n");
  PrintBox(*box1);
  printf("BOX 2: \n");
  PrintBox(*box2);
  // ...assumes orthogonal box
  double Low1[3] = {box1->Low.x, box1->Low.y, box1->Low.z},
         Low2[3] = {box2->Low.x, box2->Low.y, box2->Low.z},
         Low3[3] = {0, 0, 0}, // output box lower bound
         Length1[3] = {box1->Length.x, box1->Length.y, box1->Length.z},
         Length2[3] = {box2->Length.x, box2->Length.y, box2->Length.z},
         Length3[3] = {0, 0, 0}, // output box sidelengths
         off[3] = {offset.x, offset.y, offset.z};
         // off[3] = {-20, offset.y, offset.z};
  for (int i = 0; i < 3; i++) {
    if (Low1[i] < (Low2[i] + off[i])) {
      Low3[i] = Low1[i];
    } else {
      Low3[i] = Low2[i] + off[i];
    }
    if ((Low1[i] + Length1[i]) > (Low2[i] + Length2[i] + off[i])) {
      Length3[i] = Low1[i] + Length1[i];
    } else {
      Length3[i] = Low2[i] + Length2[i] + off[i];
    }
    Length3[i] -= Low3[i];
  }
  // fill output box Low & Length
  box_out.Length.x = Length3[0];
  box_out.Length.y = Length3[1];
  box_out.Length.z = Length3[2];
  box_out.Low.x = Low3[0];
  box_out.Low.y = Low3[1];
  box_out.Low.z = Low3[2];
  // assume orthogonal box
  box_out.alpha = 90;
  box_out.beta = 90;
  box_out.gamma = 90;
  CalculateBoxData(&box_out, 0); //}}}
  printf("BOX OUT: \n");
  PrintBox(box_out);
  ConcatenateSystems(&Sys_out, Sys_2, box_out);

  // if -box option is present, use it as box size //{{{
  if (box_opt.x != 0) {
    printf("-box & --centre OPTS\n");
    PrintBox(Sys_out.Box);
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
    // Sys_out.Box.Low.x = 0;
    // Sys_out.Box.Low.y = 0;
    // Sys_out.Box.Low.z = 0;
    CalculateBoxData(&Sys_out.Box, 0);
    PrintBox(Sys_out.Box);
    printf("-box & --centre OPTS\n");
    AddLow(&Sys_out);
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
