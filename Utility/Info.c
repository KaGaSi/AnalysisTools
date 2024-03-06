#include "../AnalysisTools.h"

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
can be modified using a second structure file (-i[!] option) and/or \
a coordinate file (-c option); see Examples/Info folder for details.\
\n\n");
  }

  fprintf(ptr, "Usage: %s <input> [options]\n\n", cmd);
  fprintf(ptr, "<input>             input structure file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -i[!] <file>      secondary structure file\n");
  fprintf(ptr, "  -c <file>         input coordinate file\n");
  fprintf(ptr, "  -o <file>         output structure file\n");
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
  bool lmp_mass;             // --mass
  FILE_TYPE fout;            // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 6, all = common + 7, count = 0, req_arg = 1;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-i");
  strcpy(option[count++], "-i!");
  strcpy(option[count++], "-c");
  strcpy(option[count++], "-o");
  strcpy(option[count++], "-def");
  strcpy(option[count++], "--mass");
  strcpy(option[count++], "-ebt");
  OptionCheck(argc, argv, count, req_arg, common, all, option); //}}}

  count = 0; // count arguments
  OPT *opt = opt_create();
  SYS_FILES in = InitSysFiles;
  snprintf(in.stru.name, LINE, "%s", argv[++count]);
  in.stru.type = StructureFileType(in.stru.name);

  PrintCommand(stdout, argc, argv);

  // -i[!] option //{{{
  SYS_FILES extra = InitSysFiles;
  bool change_beads = false;
  if (!FileOption(argc, argv, "-i", extra.stru.name)) { // -i present?
    if (!FileOption(argc, argv, "-i!", extra.stru.name)) { // -i! present?
      change_beads = true;
    }
  }
  if (extra.stru.name[0] != '\0') {
    extra.stru.type = StructureFileType(extra.stru.name);
  } //}}}
  // input coordinate file (-c option) //{{{
  FileOption(argc, argv, "-c", in.coor.name);
  if (in.coor.name[0] != '\0') {
    in.coor.type = CoordinateFileType(in.coor.name);
    extra.coor = in.coor;
  } //}}}
  // output structure file (-o option) //{{{
  opt->fout = InitFile;
  FileOption(argc, argv, "-o", opt->fout.name);
  if (opt->fout.name[0] != '\0') {
    opt->fout.type = StructureFileType(opt->fout.name);
  } //}}}
  opt->c = CommonOptions(argc, argv, LINE, in);
  // extra bead types for data output (-ebt option)
  opt->ebt = 0;
  IntegerOption1(argc, argv, "-ebt", &opt->ebt);

  // read information from input file(s) //{{{
  SYSTEM System = ReadStructure(in, opt->c.detailed);
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
    for (int i = 0; i < System.Count.Bead; i++) {
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
    PrintMolecule(System);
  }
  SYSTEM Sys_extra;
  if (extra.stru.name[0] != '\0') {
    Sys_extra = ReadStructure(extra, opt->c.detailed);
    if (opt->c.verbose) {
      fprintf(stdout, "\n==================================================");
      printf("\nSystem in extra file (%s)", extra.stru.name);
      fprintf(stdout, "\n==================================================\n");
      VerboseOutput(Sys_extra);
      fprintf(stdout, "Information about every bead:\n");
      PrintBead(Sys_extra);
      fprintf(stdout, "\nInformation about every molecule:\n");
      PrintMolecule(Sys_extra);
    }
    // add charge, mass, and radius to bead types if possible
    for (int i = 0; i < System.Count.BeadType; i++) {
      BEADTYPE *bt = &System.BeadType[i];
      int type_extra = FindBeadType(bt->Name, Sys_extra);
      if (type_extra != -1) {
        if (bt->Charge == CHARGE) {
          bt->Charge = Sys_extra.BeadType[type_extra].Charge;
        }
        if (bt->Mass == MASS) {
          bt->Mass = Sys_extra.BeadType[type_extra].Mass;
        }
        if (bt->Radius == RADIUS) {
          bt->Radius = Sys_extra.BeadType[type_extra].Radius;
        }
      }
    }
    ChangeMolecules(&System, Sys_extra, change_beads, false);
    CheckSystem(System, extra.stru.name);
  } //}}}

  // -def option (for vsf output file) //{{{
  bool *def_type = calloc(System.Count.BeadType, sizeof *def_type);
  BeadTypeOption(argc, argv, "-def", false, def_type, System);
  opt->vsf_def = -1;
  for (int i = 0; i < System.Count.BeadType; i++) {
    if (def_type[i]) {
      opt->vsf_def = i;
      break;
    }
  }
  free(def_type); //}}}
  // use mass only for atom type definition for data output file
  opt->lmp_mass = BoolOption(argc, argv, "--mass");

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
    PrintMolecule(System);
  } //}}}

  // write the output file if required (-o option) //{{{
  if (opt->fout.name[0] != '\0') {
    if (opt->fout.type == LDATA_FILE && opt->ebt != 0) {
      for (int i = 0; i < opt->ebt; i++) {
        NewBeadType(&System.BeadType, &System.Count.BeadType, "extra", 0, 1, 1);
      }
    }
    // test if coordinates are present for structure files with coordinates
    if (opt->fout.type == LDATA_FILE ||
        opt->fout.type == LTRJ_FILE ||
        opt->fout.type == VTF_FILE) {
      bool coor = false;
      for (int i = 0; i < System.Count.BeadCoor; i++) {
        int id = System.BeadCoor[i];
        if (fabs(System.Bead[id].Position[0]) > 0.00001 ||
            fabs(System.Bead[id].Position[1]) > 0.00001 ||
            fabs(System.Bead[id].Position[2]) > 0.00001) {
          coor = true;
          break;
        }
      }
      if (!coor) {
        strcpy(ERROR_MSG, "no coordinates loaded for lammps/vtf output");
        PrintWarnFile(opt->fout.name, in.coor.name, "\0");
      }
    }
    WriteStructure(opt->fout, System, opt->vsf_def, opt->lmp_mass, argc, argv);
    if (opt->fout.type == VTF_FILE) {
      bool *write = malloc(sizeof *write * System.Count.Bead);
      InitBoolArray(write, System.Count.Bead, true);
      WriteTimestep(opt->fout, System, 1, write);
      free(write);
    }
  } //}}}

  FreeSystem(&System);
  if (extra.stru.name[0] != '\0') {
    FreeSystem(&Sys_extra);
  }
  free(opt);

  return 0;
}
