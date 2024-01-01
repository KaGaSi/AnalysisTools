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
  fprintf(ptr, "<output>            output structure/coordinate file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -o <filename>     output extra structure file\n");
  fprintf(ptr, "  -off 3×<float>|c  offset of the second system against "
          "the first ('c' to place it in the centre of the first system)\n");
  fprintf(ptr, "  -b 3×<float>      output box dimensions (orthogonal)\n");
  fprintf(ptr, "  -i1/-i2 <file>    structure file for <input1>/<input2>\n");
  fprintf(ptr, "  --detailed1/--detailed2\n"
               "                    detailed bead type recognistion "
               "(vtf input file)\n");
  fprintf(ptr, "  -pbc1/-pbc2 <file>\n"
               "                    position of pbc in input xyz file's "
               "comment line\n");
  fprintf(ptr, "  -st1/-st2 <int>   starting timestep for the input files\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 4, all = common + 9, count = 0, req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-o");
  strcpy(option[count++], "-off");
  strcpy(option[count++], "-b");
  strcpy(option[count++], "-i1");
  strcpy(option[count++], "-i2");
  strcpy(option[count++], "--detailed1");
  strcpy(option[count++], "--detailed2");
  strcpy(option[count++], "-st1");
  strcpy(option[count++], "-st2");
  OptionCheck(argc, argv, count, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // <input1> - first input coordinate file //{{{
  char coor_file_1[LINE] = "", struct_file_1[LINE] = "";
  int coor_type_1, struct_type_1 = 0;
  snprintf(coor_file_1, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file_1, &coor_type_1, struct_file_1,
                       &struct_type_1)) {
    exit(1);
  } //}}}
  // <input2> - first input coordinate file //{{{
  char coor_file_2[LINE] = "", struct_file_2[LINE] = "";
  int coor_type_2, struct_type_2 = 0;
  snprintf(coor_file_2, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file_2, &coor_type_2, struct_file_2,
                       &struct_type_2)) {
    exit(1);
  } //}}}
  // <output> - output coordinate file //{{{
  char out_file[LINE] = "";
  snprintf(out_file, LINE, "%s", argv[++count]);
  // int out_type = FullFileType(out_file, 1);
  int out_type = CoordinateFileType(out_file); //}}}

  // options before reading system data //{{{
  // output extra file (-o option) //{{{
  char out2_file[LINE] = "";
  int out2_type = -1;
  if (FileOption(argc, argv, "-o", out2_file)) {
    exit(1);
  }
  if (out2_file[0] != '\0') {
    out2_type = FileType(out2_file);
  } //}}}
  // input structure files (-i1/-i2 options) //{{{
  if (FileOption(argc, argv, "-i1", struct_file_1)) {
    exit(1);
  }
  if (struct_file_1[0] != '\0') {
    struct_type_1 = StructureFileType(struct_file_1);
  }
  if (FileOption(argc, argv, "-i2", struct_file_2)) {
    exit(1);
  }
  if (struct_file_2[0] != '\0') {
    struct_type_2 = StructureFileType(struct_file_2);
  } //}}}
  bool silent, verbose, rubbish;
  int trash;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &rubbish,
                &trash, &trash, &trash);
  // --detailed option for both input systems; copied from CommonOptions()
  bool detailed1 = BoolOption(argc, argv, "--detailed1"),
       detailed2 = BoolOption(argc, argv, "--detailed2");
  // -st option for both input systems; copied from CommonOptions()
  int start1 = 1, start2 = 1;
  if (IntegerOption(argc, argv, 1, "-st1", &trash, &start1)) {
    fprintf(stderr, "%sCommand: %s", ErrRed(), ErrColourReset());
    PrintCommand(stderr, argc, argv);
    exit(1);
  }
  if (IntegerOption(argc, argv, 1, "-st2", &trash, &start2)) {
    fprintf(stderr, "%sCommand: %s", ErrRed(), ErrColourReset());
    PrintCommand(stderr, argc, argv);
    exit(1);
  }
  // -off option //{{{
  double offset[3] = {0, 0, 0};
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-off") == 0) {
      if (argc < (i + 3) ||
          (argv[i+1][0] != 'c' && !IsRealNumber(argv[i + 1], &offset[0])) ||
          (argv[i+2][0] != 'c' && !IsRealNumber(argv[i + 2], &offset[1])) ||
          (argv[i+3][0] != 'c' && !IsRealNumber(argv[i + 3], &offset[2]))) {
        strcpy(ERROR_MSG, "wrong/missing arguments (either number or 'c')");
        PrintErrorOption("-off");
        exit(1);
      }
      for (int dd = 0; dd < 3; dd++) {
        if (argv[i+dd+1][0] == 'c') {
          offset[dd] = -11111;
        }
      }
      break;
    }
  } //}}}
  // output box dimensions //{{{
  double box_opt[3] = {0, 0, 0};
  if (DoubleOption3(argc, argv, "-b", box_opt)) {
    if (box_opt[0] <= 0 || box_opt[1] <= 0 || box_opt[2] <= 0) {
      strcpy(ERROR_MSG, "three positive numbers required");
      PrintErrorOption("-b");
      Help(argv[0], true, common, option);
      exit(1);
    }
  } //}}}
  //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM Sys_1 = ReadStructure(struct_type_1, struct_file_1, coor_type_1,
                               coor_file_1, detailed1);
  SYSTEM Sys_2 = ReadStructure(struct_type_2, struct_file_2, coor_type_2,
                               coor_file_2, detailed2);
  BOX *box1 = &Sys_1.Box, *box2 = &Sys_2.Box;

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
  int line_count = 0; // count lines in the vcf file
  for (int i = 0; i < (start1 - 1); i++) {
    if (!SkipTimestep(coor_type_1, fr, coor_file_1, struct_file_1,
                      &line_count)) {
      break;
    }
  }
  if (!ReadTimestep(coor_type_1, fr, coor_file_1, &Sys_1, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_1, "\0", "\0");
    exit(1);
  }
  fclose(fr);
  AddLow(&Sys_1); //}}}
  // read second coordinate file //{{{
  fr = OpenFile(coor_file_2, "r");
  for (int i = 0; i < (start2 - 1); i++) {
    if (!SkipTimestep(coor_type_1, fr, coor_file_1, struct_file_1,
                      &line_count)) {
      break;
    }
  }
  if (!ReadTimestep(coor_type_2, fr, coor_file_2, &Sys_2, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_2, "\0", "\0");
    exit(1);
  }
  fclose(fr);
  AddLow(&Sys_2); //}}}

  // make proper offset vector
  for (int dd = 0; dd < 3; dd++) {
    if (offset[dd] == -11111) { // a)
      offset[dd] = (box1->Low[dd] + 0.5 * box1->Length[dd]) -
                   (box2->Low[dd] + 0.5 * box2->Length[dd]);
    }
  }

  // move the beads of the second system //{{{
  for (int i = 0; i < Sys_2.Count.Bead; i++) {
    int id = Sys_2.BeadCoor[i];
    Sys_2.Bead[id].Position[0] += offset[0];
    Sys_2.Bead[id].Position[1] += offset[1];
    Sys_2.Bead[id].Position[2] += offset[2];
  } //}}}

  // create output system(s) //{{{
  SYSTEM S_out1 = CopySystem(Sys_1);
  // pick box size as the larger dimensions from the initial systems //{{{
  BOX box_out = InitBox;
  // ...assumes orthogonal box
  double Low1[3] = {box1->Low[0], box1->Low[1], box1->Low[2]},
         Low2[3] = {box2->Low[0], box2->Low[1], box2->Low[2]},
         Low3[3] = {0, 0, 0}, // output box lower bound
      Length1[3] = {box1->Length[0], box1->Length[1], box1->Length[2]},
         Length2[3] = {box2->Length[0], box2->Length[1], box2->Length[2]},
         Length3[3] = {0, 0, 0}; // output box sidelengths
  for (int dd = 0; dd < 3; dd++) {
    if (Low1[dd] < (Low2[dd] + offset[dd])) {
      Low3[dd] = Low1[dd];
    } else {
      Low3[dd] = Low2[dd] + offset[dd];
    }
    if ((Low1[dd] + Length1[dd]) > (Low2[dd] + Length2[dd] + offset[dd])) {
      Length3[dd] = Low1[dd] + Length1[dd];
    } else {
      Length3[dd] = Low2[dd] + Length2[dd] + offset[dd];
    }
    Length3[dd] -= Low3[dd];
  }
  // fill output box Low & Length
  for (int dd = 0; dd < 3; dd++) {
    box_out.Length[dd] = Length3[dd];
    box_out.Low[dd] = Low3[dd];
  }
  // assume orthogonal box
  box_out.alpha = 90;
  box_out.beta = 90;
  box_out.gamma = 90;
  CalculateBoxData(&box_out, 0); //}}}
  SYSTEM S_in;
  SYSTEM S_out2;
  if (out2_file[0] != '\0') {
    S_out2 = CopySystem(Sys_1);
    S_in = CopySystem(Sys_2);
    if (out2_type == VCF_FILE ||
        out2_type == VSF_FILE ||
        out2_type == VTF_FILE) {
      VtfSystem(&S_out2);
      VtfSystem(&S_in);
    }
    ConcatenateSystems(&S_out2, S_in, box_out);
  }
  if (out2_type == VCF_FILE ||
      out2_type == VSF_FILE ||
      out2_type == VTF_FILE) {
    VtfSystem(&S_out1);
    VtfSystem(&Sys_2);
  }
  ConcatenateSystems(&S_out1, Sys_2, box_out); //}}}

  // if -b option is present, use it as box size //{{{
  if (box_opt[0] != 0) {
    // align the centre of box_opt with the centre of the original output box
    for (int dd = 0; dd < 3; dd++) {
      S_out1.Box.Low[dd] += 0.5 * (S_out1.Box.Length[dd] - box_opt[dd]);
      S_out1.Box.Length[dd] = box_opt[dd];
    }
    CalculateBoxData(&S_out1.Box, 0);
  } //}}}

  // verbose output describing the output system //{{{
  if (verbose) {
    printf("\n==================================================");
    printf("\nNew sytem");
    printf("\n==================================================\n");
    VerboseOutput(S_out1);
  } //}}}

  // write data to output file(s) //{{{
  // make coordinates from 0 to Box.Length (Lows are added if needed)
  SubtractLow(&S_out1);
  // save all beads
  bool *write = malloc(S_out1.Count.Bead * sizeof *write);
  InitBoolArray(write, S_out1.Count.Bead, true);
  // create vsf file if output file is vcf format
  if (out_type == VCF_FILE) {
    PrintByline(out_file, argc, argv); // byline to vcf file
    out_file[strlen(out_file)-2] = 's';
    WriteStructure(VSF_FILE, out_file, S_out1, -1, false, argc, argv);
    out_file[strlen(out_file)-2] = 'c';
  } else if (out_type == VTF_FILE) {
    WriteStructure(out_type, out_file, S_out1, -1, false, argc, argv);
  } else { // some formats 'append' coordinates, not 'write' them
    FILE *out = OpenFile(out_file, "w");
    fclose(out);
  }
  WriteTimestep(out_type, out_file, S_out1, 0, write);
  if (out2_file[0] != '\0') {
    if (out2_type == VTF_FILE ||
        out2_type == VSF_FILE ||
        out2_type == FIELD_FILE) {
      WriteStructure(out2_type, out2_file, S_out2, -1, false, argc, argv);
    }
    if (out2_type == VTF_FILE ||
        out2_type == VCF_FILE ||
        out2_type == LTRJ_FILE ||
        out2_type == LDATA_FILE ||
        out2_type == CONFIG_FILE) {
      WriteTimestep(out2_type, out2_file, S_out2, 0, write);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&Sys_1);
  FreeSystem(&Sys_2);
  FreeSystem(&S_out1);
  if (out2_file[0] != '\0') {
    FreeSystem(&S_out2);
    FreeSystem(&S_in);
  }
  free(write); //}}}

  return 0;
}
