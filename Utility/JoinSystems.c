#include "../src/AnalysisTools.h"
// TODO: --real switch for the -b and -off options

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "JoinSystems connects two coordinate files, creating a new systems "
  "consisting of both systems.",

  "Usage: JoinSystems <input1> <input2> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 16, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input1>", nullptr, "first input coordinate file", OPT_ARG},
  {"<input2>", nullptr, "second input coordinate file", OPT_ARG},
  {"<output>", nullptr, "output structure/coordinate file", OPT_ARG},
  {"-o", "<filename>", "output extra structure file", OPT_EXTRA},
  {"-off", "3x<float>|c", "offset of the second system against the first "
    "('c' to place it in the centre of the first system)", OPT_EXTRA},
  {"-b", "3×<float>", "output box dimensions (orthogonal)", OPT_EXTRA},
  {"--real", nullptr, "use real coordinates for -b and -off "
    "instead of fraction of first input system's box size", OPT_EXTRA},
  {"-i1", "<file>", "structure file for <input1>", OPT_EXTRA},
  {"-i2", "<file>", "structure file for <input2>", OPT_EXTRA},
  {"-st1", "<int>", "starting timestep <input1>", OPT_EXTRA},
  {"-st2", "<int>", "starting timestep <input2>", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  int start[2];          // -st1 -st2
  double off[3], box[3]; // -off -b
  bool real;             // --real
  FILE_TYPE fout;        // -o
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
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
  // options before reading system data
  // output extra file (-o option) //{{{
  opt.fout = InitFile;
  FileOption(argc, argv, "-o", opt.fout.name);
  if (opt.fout.name[0] != '\0') {
    opt.fout.type = FileType(opt.fout.name);
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
  COMMON_OPT commons = CommonOptions(argc, argv, in[0]);
  opt.real = BoolOption(argc, argv, "--real");
  // -st option for both input systems; copied from CommonOptions()
  opt.start[0] = 1, opt.start[1] = 1;
  OneNumberOption(argc, argv, "-st1", &opt.start[0], 'i');
  OneNumberOption(argc, argv, "-st2", &opt.start[1], 'i');
  // -off option //{{{
  InitDoubleArray(opt.off, 3, 0);
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-off") == 0) {
      if (argc < (i + 3) ||
          (argv[i+1][0] != 'c' && !IsRealNumber(argv[i + 1], &opt.off[0])) ||
          (argv[i+2][0] != 'c' && !IsRealNumber(argv[i + 2], &opt.off[1])) ||
          (argv[i+3][0] != 'c' && !IsRealNumber(argv[i + 3], &opt.off[2]))) {
        err_msg("wrong/missing arguments (either number or 'c')");
        PrintErrorOption("-off");
        exit(1);
      }
      for (int dd = 0; dd < 3; dd++) {
        if (argv[i+dd+1][0] == 'c') {
          opt.off[dd] = -11111;
        }
      }
      break;
    }
  } //}}}
  // output box dimensions //{{{
  InitDoubleArray(opt.box, 3, 0);
  // if (DoubleOption3(argc, argv, "-b", opt.box)) {
  if (ThreeNumbersOption(argc, argv, "-b", opt.box, 'd')) {
    if (opt.box[0] <= 0 || opt.box[1] <= 0 || opt.box[2] <= 0) {
      err_msg("three positive numbers required");
      PrintErrorOption("-b");
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}
  // only one timestep is in lammps data file
  for (int s = 0; s < 2; s++) {
    if (in[s].coor.type == LDATA_FILE) {
      opt.start[s] = 1;
    }
  }
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM Sys[2];
  BOX *box[2] = {nullptr, nullptr};
  for (int s = 0; s < 2; s++) {
    Sys[s] = ReadStructure(in[s], false);
    box[s] = &(Sys[s].Box);
  }

  // verbose output describing the two systems to be joined //{{{
  if (commons.verbose) {
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
    for (int i = 1; i < opt.start[s]; i++) {
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
    if (opt.off[dd] == -11111) { // a) put the two centres on top of each other
      opt.off[dd] = (box[0]->Low.v[dd] + 0.5 * box[0]->Length.v[dd]) -
                     (box[1]->Low.v[dd] + 0.5 * box[1]->Length.v[dd]);
    } else if (!opt.real) {
      opt.off[dd] *= box[0]->Length.v[dd];
    }
  } //}}}

  // move the beads of the second system //{{{
  for (int i = 0; i < Sys[1].Count.Bead; i++) {
    int id = Sys[1].BeadCoor[i];
    Sys[1].Bead[id].Position.v[0] += opt.off[0];
    Sys[1].Bead[id].Position.v[1] += opt.off[1];
    Sys[1].Bead[id].Position.v[2] += opt.off[2];
  } //}}}

  // create output system(s) //{{{
  // pick box size as the larger dimensions from the initial systems //{{{
  BOX box_out = InitBox;
  // ...assumes orthogonal box
  double Low1[3] = {box[0]->Low.x, box[0]->Low.y, box[0]->Low.z},
         Low2[3] = {box[1]->Low.x, box[1]->Low.y, box[1]->Low.z},
         Low3[3] = {0, 0, 0}, // output box lower bound
         Length1[3] = {box[0]->Length.x, box[0]->Length.y, box[0]->Length.z},
         Length2[3] = {box[1]->Length.x, box[1]->Length.y, box[1]->Length.z},
         Length3[3] = {0, 0, 0}; // output box sidelengths
  for (int dd = 0; dd < 3; dd++) {
    if (Low1[dd] < (Low2[dd] + opt.off[dd])) {
      Low3[dd] = Low1[dd];
    } else {
      Low3[dd] = Low2[dd] + opt.off[dd];
    }
    if ((Low1[dd] + Length1[dd]) > (Low2[dd] + Length2[dd] + opt.off[dd])) {
      Length3[dd] = Low1[dd] + Length1[dd];
    } else {
      Length3[dd] = Low2[dd] + Length2[dd] + opt.off[dd];
    }
    Length3[dd] -= Low3[dd];
  }
  // fill output box Low & Length
  for (int dd = 0; dd < 3; dd++) {
    box_out.Length.v[dd] = Length3[dd];
    box_out.Low.v[dd] = Low3[dd];
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
  PruneSystem(&S_out, nullptr);
  // optional output file
  SYSTEM S_out_opt;
  if (opt.fout.name[0] != '\0') {
    S_out_opt = CopySystem(Sys[0]);
    ConcatenateSystems(&S_out_opt, Sys[1], box_out, false);
    if (opt.fout.type == VCF_FILE ||
        opt.fout.type == VSF_FILE ||
        opt.fout.type == VTF_FILE) {
      VtfSystem(&S_out_opt);
    }
    PruneSystem(&S_out_opt, nullptr);
  }
  //}}}

  // if -b option is present, use it as box size //{{{
  if (opt.box[0] != 0) {
    // align the centre of box_opt with the centre of the original output box
    for (int dd = 0; dd < 3; dd++) {
      S_out.Box.Low.v[dd] += 0.5 * (S_out.Box.Length.v[dd] - opt.box[dd]);
      S_out.Box.Length.v[dd] = opt.box[dd];
    }
    CalculateBoxData(&S_out.Box, 0);
  } //}}}

  // verbose output describing the output system //{{{
  if (commons.verbose) {
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
  if (!write) {
    ErrorAlloc("write");
  }
  InitBoolArray(write, S_out.Count.Bead, true);
  // write to main output file
  WriteOutput(S_out, write, fout, false, -1, argc, argv);
  if (opt.fout.name[0] != '\0') {
    // make coordinates from 0 to Box.Length
    ChangeBoxByLow(&S_out_opt, -1);
    WriteOutput(S_out_opt, write, opt.fout, false, -1, argc, argv);
  } //}}}

  // free memory //{{{
  FreeSystem(&Sys[0]);
  FreeSystem(&Sys[1]);
  FreeSystem(&S_out);
  if (opt.fout.name[0] != '\0') {
    FreeSystem(&S_out_opt);
  }
  free(write); //}}}

  return 0;
}
