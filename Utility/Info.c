#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info analyzes the provided input structure file, \
printing system composition to standard output and, optionally, producing \
an output structure file of specified format (-o option). If some information \
required in the output file is missing, \
'???\' is printed instead. The system from the input file can \
be modified using a second structure file (-i[!] option) and/or \
a coordinate file (-c option); see Examples/Info folder for \
details.\
\n\n");
  }

  fprintf(ptr, "Usage: %s <input> [options]\n\n", cmd);
  fprintf(ptr, "<input>             input structure file\n");
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

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 7, all = common + 7, count = 0, req_arg = 1;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  // strcpy(option[count++], "--variable"); // TODO: makes no sense, I think
  strcpy(option[count++], "-pbc");
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

  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count arguments

  char struct_file[LINE] = "";
  snprintf(struct_file, LINE, "%s", argv[++count]);
  int struct_type = StructureFileType(struct_file, 0);

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // -i[!] option //{{{
  char struct_file_extra[LINE] = "";
  int struct_type_extra = -1;
  bool change_beads = false;
  if (FileOption(argc, argv, "-i", struct_file_extra)) {
    exit(1);
  }
  if (struct_file_extra[0] == '\0') {
    if (FileOption(argc, argv, "-i!", struct_file_extra)) {
      exit(1);
    }
    if (struct_file_extra[0] != '\0') {
      change_beads = true;
    }
  }
  if (struct_file_extra[0] != '\0') {
    struct_type_extra = StructureFileType(struct_file_extra, 0);
  } //}}}
  // input coordinate file (-c option) //{{{
  char coor_file[LINE] = "";
  int coor_type = -1;
  if (FileOption(argc, argv, "-c", coor_file)) {
    exit(1);
  }
  if (coor_file[0] != '\0') {
    coor_type = CoordinateFileType(coor_file, 0);
  } //}}}
  // output structure file (-o option) //{{{
  char struct_file_out[LINE] = "";
  int struct_type_out = -1;
  if (FileOption(argc, argv, "-o", struct_file_out)) {
    exit(1);
  }
  if (struct_file_out[0] != '\0') {
    struct_type_out = StructureFileType(struct_file_out, 1);
    if (struct_type_out != VSF_FILE &&
        struct_type_out != VTF_FILE &&
        struct_type_out != LDATA_FILE &&
        struct_type_out != FIELD_FILE) {
      strcpy(ERROR_MSG, "accepted output structure file are "
             "vsf, lammps data, or FIELD file");
      if (snprintf(ERROR_MSG, LINE, "output file %s%s%s is not in accepted "
                   "format (vsf, lammps data, or FIELD file)",
                   ErrYellow(), struct_file_out, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
  } //}}}
  bool silent, verbose, detailed;
  int timestep = 1, pbc_xyz = -1,
      trash[1]; // some stuff for unused things in options
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &timestep, trash, trash);
  // extra bead types for data output (-ebt option)
  int extra_types = 0;
  if (IntegerOption(argc, argv, 1, "-ebt", trash, &extra_types)) {
    exit(1);
  }

  // read information from input file(s) //{{{
  SYSTEM System = ReadStructure(struct_type, struct_file,
                                coor_type, coor_file, detailed, pbc_xyz);
  // use coordinate from a separate file (-c option)
  if (coor_type != -1) {
    int line_count = 0;
    FILE *fr = OpenFile(coor_file, "r");
    for (int i = 1; i < timestep; i++) { // from 1 as timestep=1 is the first
      SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count);
    }
    ReadTimestep(coor_type, fr, coor_file, &System, &line_count);
    fclose(fr);
  } else {
    // all beads are in the timestep
    for (int i = 0; i < System.Count.Bead; i++) {
      System.Bead[i].InTimestep = true;
    }
  }
  if (verbose) {
    char coor[LINE] = "\0";
    if (coor_type != -1) {
      if (snprintf(coor, LINE, " (%s)", coor_file) < 0) {
        ErrorSnprintf();
      }
    }
    printf("System in %s%s:\n", struct_file, coor);
    VerboseOutput(System);
    PrintBead(System);
  }
  SYSTEM Sys_extra;
  if (struct_file_extra[0] != '\0') {
    Sys_extra = ReadStructure(struct_type_extra, struct_file_extra,
                              coor_type, coor_file, detailed, pbc_xyz);
    if (verbose) {
      printf("System in %s:\n", struct_file_extra);
      VerboseOutput(Sys_extra);
    }
  }
  // add extra info to original system
  if (struct_file_extra[0] != '\0') {
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
    CheckSystem(System, struct_file_extra);
    WarnChargedSystem(System, struct_file, struct_file_extra, "\0");
  }
 //}}}

  // -def option (for vsf output file) //{{{
  bool *def_type = calloc(System.Count.BeadType, sizeof *def_type);
  BeadTypeOption(argc, argv, "-def", false, def_type, System);
  int vsf_def_type = -1;
  for (int i = 0; i < System.Count.BeadType; i++) {
    if (def_type[i]) {
      vsf_def_type = i;
      break;
    }
  }
  free(def_type); //}}}
  // use mass only for atom type definition for data output file
  bool lmp_mass = BoolOption(argc, argv, "--mass");

  PruneSystem(&System);

  printf("Final system composition:\n");
  VerboseOutput(System);
  if (verbose) { // -v option
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(System);
  }

  if (struct_file_out[0] != '\0') {
    if (struct_type_out == LDATA_FILE && extra_types != 0) {
      for (int i = 0; i < extra_types; i++) {
        NewBeadType(&System.BeadType, &System.Count.BeadType, "extra", 0, 1, 1);
      }
    }
    WriteStructure(struct_type_out, struct_file_out, System,
                   vsf_def_type, lmp_mass);
    if (struct_type_out == VTF_FILE) {
      bool *write = malloc(sizeof *write * System.Count.Bead);
      InitBoolArray(write, System.Count.Bead, true);
      WriteTimestep(struct_type_out, struct_file_out, System, 1, write);
      free(write);
    }
  }

  FreeSystem(&System);
  if (struct_file_extra[0] != '\0') {
    FreeSystem(&Sys_extra);
  }

  return 0;
}
