#include "../AnalysisTools.h"
#include <string.h>

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info analyzes the provided input structure file, printing system composition \
to standard output and, optionally, producing an output structure file \
of specified format (-o option). If some information required in the output \
file is missing, '???\' is printed instead. The system from the input file \
can be modified using a second structure file (-i option) and/or \
a coordinate file (-c option); see Examples/Info folder for details.\
\n\n");
  }

  fprintf(ptr, "Usage: %s <input> [options]\n\n", cmd);
  fprintf(ptr, "<input>             input structure file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -i <file>         secondary structure file\n");
  fprintf(ptr, "  -c <file>         input coordinate file\n");
  fprintf(ptr, "  --detailed        use name as well as charge, mass, and "
          "radius to identfy bead types (input vtf structure files only)\n");
  fprintf(ptr, "  -o <file>         output structure file\n");
  fprintf(ptr, "  --unique          make all bead/molecule names unique\n");
  fprintf(ptr, "  -def <bead name>  default bead type "
          "(output vtf structure file only)\n");
  fprintf(ptr, "  --mass            define lammps atom types by mass, but "
          "print per-atom charges in Atoms section "
          "(output lammps data file only)\n");
  fprintf(ptr, "  -ebt <int>        number of extra bead types "
          "(output lammps data file only)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int vsf_def, ebt;          // -def -ebt
  bool lmp_mass, detailed;   // --mass --detailed
  FILE_TYPE fout;            // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 5, all = common + 8, count = 0, req_arg = 1;
  char option[all][OPT_LENGTH];
  // common options
  s_strcpy(option[count++], "-st", OPT_LENGTH);
  s_strcpy(option[count++], "--verbose", OPT_LENGTH);
  s_strcpy(option[count++], "--silent", OPT_LENGTH);
  s_strcpy(option[count++], "--help", OPT_LENGTH);
  s_strcpy(option[count++], "--version", OPT_LENGTH);
  // extra options
  s_strcpy(option[count++], "-i", OPT_LENGTH);
  s_strcpy(option[count++], "-c", OPT_LENGTH);
  s_strcpy(option[count++], "-o", OPT_LENGTH);
  s_strcpy(option[count++], "--unique", OPT_LENGTH);
  s_strcpy(option[count++], "--detailed", OPT_LENGTH);
  s_strcpy(option[count++], "-def", OPT_LENGTH);
  s_strcpy(option[count++], "--mass", OPT_LENGTH);
  s_strcpy(option[count++], "-ebt", OPT_LENGTH);
  OptionCheck(argc, argv, count, req_arg, common, all, option, true); //}}}

  count = 0; // count arguments
  OPT *opt = opt_create();
  SYS_FILES in = InitSysFiles;
  // snprintf(in.stru.name, LINE, "%s", argv[++count]);
  s_strcpy(in.stru.name, argv[++count], LINE);
  in.stru.type = StructureFileType(in.stru.name);

  PrintCommand(stdout, argc, argv);

  // -i option //{{{
  SYS_FILES extra = InitSysFiles;
  if (FileOption(argc, argv, "-i", extra.stru.name)) {
    extra.stru.type = StructureFileType(extra.stru.name);
  } //}}}
  // input coordinate file (-c option) //{{{
  FileOption(argc, argv, "-c", in.coor.name);
  if (in.coor.name[0] != '\0') {
    in.coor.type = CoordinateFileType(in.coor.name);
    extra.coor = in.coor;
  } //}}}
  // output file (-o option) //{{{
  opt->fout = InitFile;
  if (FileOption(argc, argv, "-o", opt->fout.name)) {
    opt->fout.type = FileType(opt->fout.name);
  } //}}}
  opt->c = CommonOptions(argc, argv, LINE, in);
  // extra bead types for data output (-ebt option)
  opt->ebt = 0;
  IntegerOption1(argc, argv, "-ebt", &opt->ebt);
  // use mass only for atom type definition for data output file
  opt->lmp_mass = BoolOption(argc, argv, "--mass");
  // base bead types on name, charge, mass, and radius (vtf input file)
  opt->detailed = BoolOption(argc, argv, "--detailed");

  // read information from input file(s) //{{{
  SYSTEM System = ReadStructure(in, opt->detailed);
  COUNT *Count = &System.Count;
  // use coordinate from a separate file (-c option)
  if (in.coor.type != -1) {
    int line_count = 0;
    FILE *fr = OpenFile(in.coor.name, "r");
    for (int i = 1; i < opt->c.start; i++) { // from 1 as timestep=1 is the first
      SkipTimestep(in, fr, &line_count);
    }
    ReadTimestep(in, fr, &System, &line_count);
    fclose(fr);
  } else {
    // all beads are in the timestep
    for (int i = 0; i < Count->Bead; i++) {
      System.Bead[i].InTimestep = true;
    }
  }
  // print initial system information only if extra file(s) are present
  if (opt->c.verbose && (extra.stru.type != -1 || in.coor.type != -1)) {
    fprintf(stdout, "\n==================================================");
    printf("\nSystem in %s", in.stru.name);
    if (in.coor.type != -1) {
      printf(" (coordinates: %s)", in.coor.name);
    }
    fprintf(stdout, "\n==================================================\n");
    VerboseOutput(System);
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecules(System);
  }
  SYSTEM Sys_extra;
  if (extra.stru.name[0] != '\0') {
    Sys_extra = ReadStructure(extra, opt->detailed);
    if (opt->c.verbose) {
      fprintf(stdout, "\n==================================================");
      printf("\nSystem in extra file (%s)", extra.stru.name);
      fprintf(stdout, "\n==================================================\n");
      VerboseOutput(Sys_extra);
      fprintf(stdout, "Information about every bead:\n");
      PrintBead(Sys_extra);
      fprintf(stdout, "\nInformation about every molecule:\n");
      PrintMolecules(Sys_extra);
    }
    // add charge, mass, and radius to bead types if possible
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt = &System.BeadType[i];
      int type_extra = FindBeadType(bt->Name, Sys_extra);
      if (type_extra != -1) {
        BEADTYPE *bt_extra = &Sys_extra.BeadType[type_extra];
        if (bt->Charge == CHARGE) {
          bt->Charge = bt_extra->Charge;
        }
        if (bt->Mass == MASS) {
          bt->Mass = bt_extra->Mass;
        }
        if (bt->Radius == RADIUS) {
          bt->Radius = bt_extra->Radius;
        }
      }
    }
    ChangeMolecules(&System, Sys_extra, false);
    CheckSystem(System, extra.stru.name);
  } //}}}

  // -def option (for vsf output file) //{{{
  bool *def_type = calloc(Count->BeadType, sizeof *def_type);
  BeadTypeOption(argc, argv, "-def", false, def_type, System);
  opt->vsf_def = -1;
  for (int i = 0; i < Count->BeadType; i++) {
    if (def_type[i]) {
      opt->vsf_def = i;
      break;
    }
  }
  free(def_type); //}}}

  PruneSystem(&System);

  // print the system information //{{{
  fprintf(stdout, "\n==================================================");
  if (extra.stru.type != -1 || in.coor.type != -1) {
    printf("\nFinal system composition");
  } else {
    printf("\nSystem composition");
  }
  fprintf(stdout, "\n==================================================\n");
  VerboseOutput(System);
  if (opt->c.verbose) { // -v option
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecules(System);
  } //}}}

  // write the output file if required (-o option) //{{{
  if (opt->fout.name[0] != '\0') {
    if (opt->fout.type == LDATA_FILE && opt->ebt > 0) {
      NewBeadType(&System.BeadType, &Count->BeadType, "extra", 0, 1, 1);
    }
    bool *write = malloc(sizeof *write * Count->Bead);
    InitBoolArray(write, Count->Bead, true);
    WriteOutput(System, write, opt->fout, opt->lmp_mass, opt->vsf_def,
                argc, argv);
    free(write);
  } //}}}

  FreeSystem(&System);
  if (extra.stru.name[0] != '\0') {
    FreeSystem(&Sys_extra);
  }
  free(opt);

  return 0;
}
