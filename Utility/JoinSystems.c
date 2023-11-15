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
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 4, all = common + 2, count = 0,
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
    struct_type_out = StructureFileType(struct_file_1, 0);
  }
  if (FileOption(argc, argv, "-i2", struct_file_2)) {
    exit(1);
  }
  if (struct_file_2[0] != '\0') {
    struct_type_out = StructureFileType(struct_file_2, 0);
  } //}}}
  // TODO: --detailed1 & --detailed2 options
  //       --pbc_xyz1 & --pbc_xyz2 options
  //       -st1 & -st2 options
  bool silent, verbose, detailed;
  int start = 1,
      trash, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &start, &trash, &trash); //}}}

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

  // read first coordinate file
  FILE *fr = OpenFile(coor_file_1, "r");
  int line_count = 0;     // count lines in the vcf file
  if (!ReadTimestep(coor_type_1, fr, coor_file_1, &Sys_1, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_1, "\0", "\0");
    exit(1);
  }
  fclose(fr);
  // read second coordinate file
  fr = OpenFile(coor_file_2, "r");
  if (!ReadTimestep(coor_type_2, fr, coor_file_2, &Sys_2, &line_count)) {
    strcpy(ERROR_MSG, "no valid timestep");
    PrintErrorFile(coor_file_2, "\0", "\0");
    exit(1);
  }
  fclose(fr);

  SYSTEM Sys_out = CopySystem(Sys_1);
  BOX box_out = Sys_out.Box;
  ConcatenateSystems(&Sys_out, Sys_2, box_out);

  // print initial stuff to output coordinate file //{{{
  if (coor_out_type == VCF_FILE) {
    FILE *out = OpenFile(coor_out_file, "w");
    PrintByline(out, argc, argv);
    fclose(out);
  } else if (coor_out_type == VTF_FILE) {
    WriteStructure(VSF_FILE, coor_out_file, Sys_out, -1, false);
    coor_out_type = VCF_FILE;
  } else {
    FILE *out = OpenFile(coor_out_file, "w");
    fclose(out);
  } //}}}

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
