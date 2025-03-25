#include "../AnalysisTools.h"

// TODO: --real switch for the -b and -off options

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
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
  fprintf(ptr, "  -st1/-st2 <int>   starting timestep for the input files\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int start[2];              // -st1 -st2
  double off[3], box[3];     // -off -b
  FILE_TYPE fout;            // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 4, all = common + 7, count = 0, req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "--verbose", "--silent", "--help", "--version",
               "-o", "-off", "-b", "-i1", "-i2", "-st1", "-st2");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // input/output files //{{{
  // <input1> <input2> - input coordinate (and structure) files
  SYS_FILES in[2];
  for (int s = 0; s < 2; s++) {
    in[s] = InitSysFiles;
    s_strcpy(in[s].coor.name, argv[++count], LINE);
    if (!InputCoorStruct(argc, argv, &in[s])) {
      exit(1);
    }
  }

  // <output> - output coordinate file
  FILE_TYPE fout = InitFile;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name); //}}}

  // options before reading system data //{{{
  // output extra file (-o option) //{{{
  opt->fout = InitFile;
  FileOption(argc, argv, "-o", opt->fout.name);
  if (opt->fout.name[0] != '\0') {
    opt->fout.type = FileType(opt->fout.name);
  } //}}}
  // input structure files (-i1/-i2 options) //{{{
  char tmp[LINE] = "\0";
  if (FileOption(argc, argv, "-i1", tmp)) {
    s_strcpy(in[0].stru.name, tmp, LINE);
    in[0].stru.type = StructureFileType(in[0].stru.name);
  }
  if (FileOption(argc, argv, "-i2", tmp)) {
    s_strcpy(in[1].stru.name, tmp, LINE);
    in[1].stru.type = StructureFileType(in[1].stru.name);
  } //}}}
  opt->c = CommonOptions(argc, argv, in[0]);
  // -st option for both input systems; copied from CommonOptions()
  opt->start[0] = 1, opt->start[1] = 1;
  OneNumberOption(argc, argv, "-st1", &opt->start[0], 'i');
  OneNumberOption(argc, argv, "-st2", &opt->start[1], 'i');
  // -off option //{{{
  InitDoubleArray(opt->off, 3, 0);
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-off") == 0) {
      if (argc < (i + 3) ||
          (argv[i+1][0] != 'c' && !IsRealNumber(argv[i + 1], &opt->off[0])) ||
          (argv[i+2][0] != 'c' && !IsRealNumber(argv[i + 2], &opt->off[1])) ||
          (argv[i+3][0] != 'c' && !IsRealNumber(argv[i + 3], &opt->off[2]))) {
        err_msg("wrong/missing arguments (either number or 'c')");
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
  // if (DoubleOption3(argc, argv, "-b", opt->box)) {
  if (ThreeNumbersOption(argc, argv, "-b", opt->box, 'd')) {
    if (opt->box[0] <= 0 || opt->box[1] <= 0 || opt->box[2] <= 0) {
      err_msg("three positive numbers required");
      PrintErrorOption("-b");
      Help(StripPath(argv[0]), true, common, option);
      exit(1);
    }
  } //}}}
  // only one timestep is in lammps data file
  for (int s = 0; s < 2; s++) {
    if (in[s].coor.type == LDATA_FILE) {
      opt->start[s] = 1;
    }
  }
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM Sys[2];
  BOX *box[2] = {NULL, NULL};
  for (int s = 0; s < 2; s++) {
    Sys[s] = ReadStructure(in[s], false);
    box[s] = &(Sys[s].Box);
  }

  // verbose output describing the two systems to be joined //{{{
  if (opt->c.verbose) {
    printf("\n==================================================");
    printf("\nFirst sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys[0]);
    printf("\n==================================================");
    printf("\nSecond sytem");
    printf("\n==================================================\n");
    VerboseOutput(Sys[1]);
  } //}}}

  // read coordinate files //{{{
  for (int s = 0; s < 2; s++) {
    FILE *fr = OpenFile(in[s].coor.name, "r");
    int line_count = 0; // count lines in the vcf file
    for (int i = 1; i < opt->start[s]; i++) {
      if (!SkipTimestep(in[s], fr, &line_count)) {
        break;
      }
    }
    if (!ReadTimestep(in[s], fr, &Sys[s], &line_count)) {
      err_msg("no valid timestep; maybe starting step is too high?");
      PrintErrorFile(in[s].coor.name, "\0", "\0");
      exit(1);
    }
    fclose(fr);
    // make the two systems be correctly oriented with respect to each other
    ChangeBoxByLow(&Sys[s], +1);
  } //}}}

  // make proper offset vector //{{{
  for (int dd = 0; dd < 3; dd++) {
    if (opt->off[dd] == -11111) { // a)
      opt->off[dd] = (box[0]->Low[dd] + 0.5 * box[0]->Length[dd]) -
                     (box[1]->Low[dd] + 0.5 * box[1]->Length[dd]);
    }
  } //}}}

  // move the beads of the second system //{{{
  for (int i = 0; i < Sys[1].Count.Bead; i++) {
    int id = Sys[1].BeadCoor[i];
    Sys[1].Bead[id].Position[0] += opt->off[0];
    Sys[1].Bead[id].Position[1] += opt->off[1];
    Sys[1].Bead[id].Position[2] += opt->off[2];
  } //}}}

  // create output system(s) //{{{
  // pick box size as the larger dimensions from the initial systems //{{{
  BOX box_out = InitBox;
  // ...assumes orthogonal box
  double Low1[3] = {box[0]->Low[0], box[0]->Low[1], box[0]->Low[2]},
         Low2[3] = {box[1]->Low[0], box[1]->Low[1], box[1]->Low[2]},
         Low3[3] = {0, 0, 0}, // output box lower bound
         Length1[3] = {box[0]->Length[0], box[0]->Length[1], box[0]->Length[2]},
         Length2[3] = {box[1]->Length[0], box[1]->Length[1], box[1]->Length[2]},
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
  // main output file
  // SYSTEM S_in = CopySystem(Sys[1]);
  SYSTEM S_out = CopySystem(Sys[0]);
  ConcatenateSystems(&S_out, Sys[1], box_out, false);
  // TODO: the whole VtfSystem() stuff - what's it for, anyway?
  if (fout.type == VCF_FILE ||
      fout.type == VSF_FILE ||
      fout.type == VTF_FILE) {
    VtfSystem(&S_out);
  }
  PruneSystem(&S_out);
  // optional output file
  SYSTEM S_out_opt;
  if (opt->fout.name[0] != '\0') {
    S_out_opt = CopySystem(Sys[0]);
    ConcatenateSystems(&S_out_opt, Sys[1], box_out, false);
    if (opt->fout.type == VCF_FILE ||
        opt->fout.type == VSF_FILE ||
        opt->fout.type == VTF_FILE) {
      VtfSystem(&S_out_opt);
    }
    PruneSystem(&S_out_opt);
  }
  //}}}

  // if -b option is present, use it as box size //{{{
  if (opt->box[0] != 0) {
    // align the centre of box_opt with the centre of the original output box
    for (int dd = 0; dd < 3; dd++) {
      S_out.Box.Low[dd] += 0.5 * (S_out.Box.Length[dd] - opt->box[dd]);
      S_out.Box.Length[dd] = opt->box[dd];
    }
    CalculateBoxData(&S_out.Box, 0);
  } //}}}

  // verbose output describing the output system //{{{
  if (opt->c.verbose) {
    printf("\n==================================================");
    printf("\nNew sytem");
    printf("\n==================================================\n");
    VerboseOutput(S_out);
  } //}}}

  // write data to output file(s) //{{{
  // make coordinates from 0 to Box.Length
  ChangeBoxByLow(&S_out, -1);
  // save all beads
  bool *write = malloc(S_out.Count.Bead * sizeof *write);
  InitBoolArray(write, S_out.Count.Bead, true);
  // write to main output file
  WriteOutput(S_out, write, fout, false, -1, argc, argv);
  if (opt->fout.name[0] != '\0') {
    // make coordinates from 0 to Box.Length
    ChangeBoxByLow(&S_out_opt, -1);
    WriteOutput(S_out_opt, write, opt->fout, false, -1, argc, argv);
  } //}}}

  // free memory //{{{
  FreeSystem(&Sys[0]);
  FreeSystem(&Sys[1]);
  FreeSystem(&S_out);
  if (opt->fout.name[0] != '\0') {
    FreeSystem(&S_out_opt);
  }
  free(write);
  free(opt); //}}}

  return 0;
}
