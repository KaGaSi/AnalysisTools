#include "../AnalysisTools.h"

// TODO: --real switch for the -b and -off options

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
  fprintf(ptr, "  -st1/-st2 <int>   starting timestep for the input files\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool detailed1, detailed2; // --detailed1 --detailed2
  int start1, start2;        // -st1 -st2
  double off[3], box[3];     // -off -b
  FILE_TYPE fout;            // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
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
  OPT *opt = opt_create();
  // input/output files //{{{
  // <input1> - first input coordinate file
  SYS_FILES in1 = InitSysFiles;
  snprintf(in1.coor.name, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, &in1)) {
    exit(1);
  }
  // <input2> - first input coordinate file
  SYS_FILES in2 = InitSysFiles;
  snprintf(in2.coor.name, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, &in2)) {
    exit(1);
  }
  // <output> - output coordinate file
  FILE_TYPE fout = InitFile;
  snprintf(fout.name, LINE, "%s", argv[++count]);
  fout.type = CoordinateFileType(fout.name); //}}}

  // options before reading system data //{{{
  // output extra file (-o option) //{{{
  opt->fout = InitFile;
  FileOption(argc, argv, "-o", opt->fout.name);
  if (opt->fout.name[0] != '\0') {
    opt->fout.type = FileType(opt->fout.name);
  } //}}}
  // input structure files (-i1/-i2 options) //{{{
  if (FileOption(argc, argv, "-i1", in1.stru.name)) {
    in1.stru.type = StructureFileType(in1.stru.name);
  }
  if (FileOption(argc, argv, "-i2", in2.stru.name)) {
    in2.stru.type = StructureFileType(in2.stru.name);
  } //}}}
  opt->c = CommonOptions(argc, argv, LINE);
  // --detailed option for both input systems; copied from CommonOptions()
  opt->detailed1 = BoolOption(argc, argv, "--detailed1"),
  opt->detailed2 = BoolOption(argc, argv, "--detailed2");
  // -st option for both input systems; copied from CommonOptions()
  opt->start1 = 1, opt->start2 = 1;
  IntegerOption1(argc, argv, "-st1", &opt->start1);
  IntegerOption1(argc, argv, "-st2", &opt->start2);
  // -off option //{{{
  InitDoubleArray(opt->off, 3, 0);
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-off") == 0) {
      if (argc < (i + 3) ||
          (argv[i+1][0] != 'c' && !IsRealNumber(argv[i + 1], &opt->off[0])) ||
          (argv[i+2][0] != 'c' && !IsRealNumber(argv[i + 2], &opt->off[1])) ||
          (argv[i+3][0] != 'c' && !IsRealNumber(argv[i + 3], &opt->off[2]))) {
        strcpy(ERROR_MSG, "wrong/missing arguments (either number or 'c')");
        PrintErrorOption("-off");
        exit(1);
      }
      for (int dd = 0; dd < 3; dd++) {
        if (argv[i+dd+1][0] == 'c') {
          opt->off[dd] = -11111;
        }
      }
      break;
    }
  } //}}}
  // output box dimensions //{{{
  InitDoubleArray(opt->box, 3, 0);
  if (DoubleOption3(argc, argv, "-b", opt->box)) {
    if (opt->box[0] <= 0 || opt->box[1] <= 0 || opt->box[2] <= 0) {
      strcpy(ERROR_MSG, "three positive numbers required");
      PrintErrorOption("-b");
      Help(argv[0], true, common, option);
      exit(1);
    }
  } //}}}
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM Sys_1 = ReadStructure(in1, opt->detailed1);
  SYSTEM Sys_2 = ReadStructure(in2, opt->detailed2);
  BOX *box1 = &Sys_1.Box, *box2 = &Sys_2.Box;

  // verbose output describing the two systems to be joined //{{{
  if (opt->c.verbose) {
    printf("\n==================================================");
    printf("\nFirst sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys_1);
    printf("\n==================================================");
    printf("\nSecond sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys_2);
  } //}}}

  // read coordinate file //{{{
  // first coordinate file
  FILE *fr = OpenFile(in1.coor.name, "r");
  int line_count = 0; // count lines in the vcf file
  for (int i = 0; i < (opt->start1 - 1); i++) {
    if (!SkipTimestep(in1, fr, &line_count)) {
      break;
    }
  }
  if (!ReadTimestep(in1, fr, &Sys_1, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(in1.coor.name, "\0", "\0");
    exit(1);
  }
  fclose(fr);
  // make the two systems be correctly oriented with respect to each other
  AddLow(&Sys_1);
  // second coordinate file
  fr = OpenFile(in2.coor.name, "r");
  for (int i = 0; i < (opt->start2 - 1); i++) {
    if (!SkipTimestep(in2, fr, &line_count)) {
      break;
    }
  }
  if (!ReadTimestep(in2, fr, &Sys_2, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(in2.coor.name, "\0", "\0");
    exit(1);
  }
  fclose(fr);
  // make the two systems be correctly oriented with respect to each other
  AddLow(&Sys_2); //}}}

  // make proper offset vector //{{{
  for (int dd = 0; dd < 3; dd++) {
    if (opt->off[dd] == -11111) { // a)
      opt->off[dd] = (box1->Low[dd] + 0.5 * box1->Length[dd]) -
                        (box2->Low[dd] + 0.5 * box2->Length[dd]);
    }
  } //}}}

  // move the beads of the second system //{{{
  for (int i = 0; i < Sys_2.Count.Bead; i++) {
    int id = Sys_2.BeadCoor[i];
    Sys_2.Bead[id].Position[0] += opt->off[0];
    Sys_2.Bead[id].Position[1] += opt->off[1];
    Sys_2.Bead[id].Position[2] += opt->off[2];
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
    if (Low1[dd] < (Low2[dd] + opt->off[dd])) {
      Low3[dd] = Low1[dd];
    } else {
      Low3[dd] = Low2[dd] + opt->off[dd];
    }
    if ((Low1[dd] + Length1[dd]) > (Low2[dd] + Length2[dd] + opt->off[dd])) {
      Length3[dd] = Low1[dd] + Length1[dd];
    } else {
      Length3[dd] = Low2[dd] + Length2[dd] + opt->off[dd];
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
  if (opt->fout.name[0] != '\0') {
    S_out2 = CopySystem(Sys_1);
    S_in = CopySystem(Sys_2);
    if (opt->fout.type == VCF_FILE ||
        opt->fout.type == VSF_FILE ||
        opt->fout.type == VTF_FILE) {
      VtfSystem(&S_out2);
      VtfSystem(&S_in);
    }
    ConcatenateSystems(&S_out2, S_in, box_out);
  }
  if (opt->fout.type == VCF_FILE ||
      opt->fout.type == VSF_FILE ||
      opt->fout.type == VTF_FILE) {
    VtfSystem(&S_out1);
    VtfSystem(&Sys_2);
  }
  ConcatenateSystems(&S_out1, Sys_2, box_out); //}}}

  // if -b option is present, use it as box size //{{{
  if (opt->box[0] != 0) {
    // align the centre of box_opt with the centre of the original output box
    for (int dd = 0; dd < 3; dd++) {
      S_out1.Box.Low[dd] += 0.5 * (S_out1.Box.Length[dd] - opt->box[dd]);
      S_out1.Box.Length[dd] = opt->box[dd];
    }
    CalculateBoxData(&S_out1.Box, 0);
  } //}}}

  // verbose output describing the output system //{{{
  if (opt->c.verbose) {
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
  if (fout.type == VCF_FILE) {
    PrintByline(fout.name, argc, argv); // byline to vcf file
    fout.name[strlen(fout.name)-2] = 's';
    fout.type = VSF_FILE;
    WriteStructure(fout, S_out1, -1, false, argc, argv);
    fout.name[strlen(fout.name)-2] = 'c';
    fout.type = VCF_FILE;
  } else if (fout.type == VTF_FILE) {
    WriteStructure(fout, S_out1, -1, false, argc, argv);
  } else { // some formats 'append' coordinates, not 'write' them
    FILE *out = OpenFile(fout.name, "w");
    fclose(out);
  }
  WriteTimestep(fout, S_out1, 1, write);
  if (opt->fout.name[0] != '\0') {
    if (opt->fout.type == VTF_FILE ||
        opt->fout.type == VSF_FILE ||
        opt->fout.type == FIELD_FILE) {
      WriteStructure(opt->fout, S_out2, -1, false, argc, argv);
    }
    if (opt->fout.type == VTF_FILE ||
        opt->fout.type == VCF_FILE ||
        opt->fout.type == LTRJ_FILE ||
        opt->fout.type == LDATA_FILE ||
        opt->fout.type == CONFIG_FILE) {
      WriteTimestep(opt->fout, S_out2, 0, write);
    }
  } //}}}

  FreeSystem(&Sys_1);
  FreeSystem(&Sys_2);
  FreeSystem(&S_out1);
  if (opt->fout.name[0] != '\0') {
    FreeSystem(&S_out2);
    FreeSystem(&S_in);
  }
  free(write);
  free(opt);

  return 0;
}
