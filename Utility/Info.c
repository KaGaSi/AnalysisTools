#include "../AnalysisTools.h"
// TODO: pbc_xyz in ReadTimestep

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info analyzes the provided input structure file, \
printing system composition to standard output. It can also \
modify the system using a second structure file (-i[!] option). \
Both bead types and molecule types can be modified: when a bead \
type has unspecified mass, charge, or radius, \
these values are be taken from the bead type of the same name from \
the second file; when a molecule type misses bond types, angles, \
dihedrals, impropers, or their types, these are taken from the molecule type \
of the same name and number of beads from the second file. If '!' is used, \
the bead types in the original molecules are changed for bead types from \
molecules in the second system (for molecules that share the name and \
number of beads). See the manual or Examples/Info for details and examples. \
Info can also print the resulting system into an output file of the specified \
type, including (if possible) coordinates from a different coordinate file. \
If some information required for the given output file type is missing, \
'???\' is printed instead.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> [options]\n\n", cmd);
  fprintf(ptr, "   <input>               input structure file\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "    input files:\n");
  fprintf(ptr, "      -i[!] <file>       secondary structure file\n");
  fprintf(ptr, "      --detailed         differentiate bead types"
          " not just by names (works only with vtf structure file)\n");
  fprintf(ptr, "      -c <file>          input coordinate file\n");
  fprintf(ptr, "      -st <int>          what timestep to use;"
          " default: 1; for last, use '0' (works with -x_in or -vc_in)\n");
  fprintf(ptr, "    output files:\n");
  fprintf(ptr, "      -o <file>          output structure file\n");
  fprintf(ptr, "      -def <bead name>   default bead type"
          " (works with vsf output file)\n");
  fprintf(ptr, "      --mass             define lammps atom types by mass, but"
          " print per-atom charges in Atoms section (works with data file)\n");
  fprintf(ptr, "      -pbc <int>         position of pbc in xyz file's comment"
          " line (of the first number)\n");
  fprintf(ptr, "      -v                 more verbose output\n");
  fprintf(ptr, "      -h                 print this help and exit\n");
  fprintf(ptr, "      --version          print version and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && strcmp(argv[i], "--mass") != 0 &&
        strcmp(argv[i], "-i") != 0 && strcmp(argv[i], "-i!") != 0 &&
        strcmp(argv[i], "-o") != 0 && strcmp(argv[i], "-c") != 0 &&
        strcmp(argv[i], "-st") != 0 && strcmp(argv[i], "-def") != 0 &&
        strcmp(argv[i], "--mass") != 0 && strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-h") != 0 && strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-pbc") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "--variable")) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // <input> //{{{
  char struct_file[LINE] = "";
  int struct_type;
  snprintf(struct_file, LINE, "%s", argv[++count]);
  if (strcasecmp(struct_file, "FIELD") == 0) {
    struct_type = FIELD_FILE;
  } else {
    int ext = 6;
    char extension[ext][EXTENSION];
    strcpy(extension[0], ".vsf");
    strcpy(extension[1], ".vtf");
    strcpy(extension[2], ".xyz");
    strcpy(extension[3], ".lammpstrj");
    strcpy(extension[4], ".data");
    strcpy(extension[5], ".field");
    ext = ErrorExtension(struct_file, ext, extension);
    switch (ext) {
      case 0:
        struct_type = VSF_FILE;
        break;
      case 1:
        struct_type = VSF_FILE;
        break;
      case 2:
        struct_type = XYZ_FILE;
        break;
      case 3:
        struct_type = LTRJ_FILE;
        break;
      case 4:
        struct_type = LDATA_FILE;
        break;
      case 5:
        struct_type = FIELD_FILE;
        break;
      default: // wrong extension
        exit(1);
    }
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // -i[!] option //{{{
  char struct_file_extra[LINE] = "";
  int struct_type_extra = -1;
  bool change_beads = false;
  if (FileOption(argc, argv, "-i", struct_file_extra, LINE)) {
    exit(1);
  }
  if (struct_file_extra[0] == '\0') {
    if (FileOption(argc, argv, "-i!", struct_file_extra, LINE)) {
      exit(1);
    }
    if (struct_file_extra[0] != '\0') {
      change_beads = true;
    }
  }
  if (struct_file_extra[0] != '\0') {
    if (strcasecmp(struct_file_extra, "FIELD") == 0) {
      struct_type_extra = FIELD_FILE;
    } else {
      int ext = 6;
      char extension[ext][EXTENSION];
      strcpy(extension[0], ".vsf");
      strcpy(extension[1], ".vtf");
      strcpy(extension[2], ".xyz");
      strcpy(extension[3], ".lammpstrj");
      strcpy(extension[4], ".data");
      strcpy(extension[5], ".field");
      ext = ErrorExtension(struct_file_extra, ext, extension);
      switch (ext) {
        case 0:
          struct_type_extra = VSF_FILE;
          break;
        case 1:
          struct_type_extra = VSF_FILE;
          break;
        case 2:
          struct_type_extra = XYZ_FILE;
          break;
        case 3:
          struct_type_extra = LTRJ_FILE;
          break;
        case 4:
          struct_type_extra = LDATA_FILE;
          break;
        case 5:
          struct_type_extra = FIELD_FILE;
          break;
        default: // wrong extension
          exit(1);
      }
    }
  } //}}}
  // input coordinate file (-c option) //{{{
  char coor_file[LINE] = "";
  int coor_type = -1;
  if (FileOption(argc, argv, "-c", coor_file, LINE)) {
    exit(1);
  }
  if (coor_file[0] != '\0') {
    int ext = 4;
    char extension[ext][EXTENSION];
    strcpy(extension[0], ".vcf");
    strcpy(extension[1], ".vtf");
    strcpy(extension[2], ".xyz");
    strcpy(extension[3], ".lammpstrj");
    ext = ErrorExtension(coor_file, ext, extension);
    // define coordinate type and possibly vtf structure file
    switch (ext) {
      case 0:
        coor_type = VCF_FILE;
        break;
      case 1:
        coor_type = VCF_FILE;
        break;
      case 2:
        coor_type = XYZ_FILE;
        break;
      case 3: // lammpstrj
        coor_type = LTRJ_FILE;
        break;
      default: // wrong extenstion
        exit(1);
    }
  } //}}}
  // output structure file (-o option) //{{{
  char struct_file_out[LINE] = "";
  int struct_type_out = -1;
  if (FileOption(argc, argv, "-o", struct_file_out, LINE)) {
    exit(1);
  }
  if (struct_file_out[0] != '\0') {
    if (strcasecmp(struct_file_out, "FIELD") == 0) {
      struct_type_out = FIELD_FILE;
    } else if (strcasecmp(struct_file_out, "CONFIG") == 0) {
      struct_type_out = CONFIG_FILE;
    } else {
      int ext = 3;
      char extension[ext][EXTENSION];
      strcpy(extension[0], ".vsf");
      strcpy(extension[1], ".data");
      strcpy(extension[2], ".field");
      ext = ErrorExtension(struct_file_out, ext, extension);
      switch (ext) {
        case 0:
          struct_type_out = VSF_FILE;
          break;
        case 1:
          struct_type_out = LDATA_FILE;
          break;
        case 2:
          struct_type_out = FIELD_FILE;
          break;
        default: // wrong extension
          exit(1);
      }
    }
  } //}}}
  // timestep to use coordinates from (-st option) //{{{
  int timestep = 1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  }
  timestep--; //}}}
  bool detailed = BoolOption(argc, argv, "--detailed");
  // vtf timesteps with variable number of beads
  bool vtf_var_coor = BoolOption(argc, argv, "--variable");
  bool verbose = BoolOption(argc, argv, "-v");
  // position of the first number of pbc in xyz file //{{{
  int pbc_xyz = -1;
  if (IntegerOption(argc, argv, "-pbc", &pbc_xyz)) {
    exit(1);
  }
  if (pbc_xyz == 0) {
    strcpy(ERROR_MSG, "position must be a positive number");
    PrintErrorOption("-pbc");
  } //}}}

  // read information from input file(s) //{{{
  int ltrj_start_id = -1;
  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var_coor,
                                pbc_xyz, &ltrj_start_id);
  if (verbose) {
    printf("System in %s:\n", struct_file);
    VerboseOutput(System);
  }
  SYSTEM Sys_extra;
  if (struct_file_extra[0] != '\0') {
    Sys_extra = ReadStructure(struct_type_extra, struct_file_extra,
                              coor_type, coor_file, detailed, vtf_var_coor,
                              pbc_xyz, &ltrj_start_id);
    if (verbose) {
      printf("System in %s:\n", struct_file_extra);
      VerboseOutput(Sys_extra);
    }
  }
  // all beads are in the timestep; revised if coordinates supplied
  for (int i = 0; i < System.Count.Bead; i++) {
    System.Bead[i].InTimestep = true;
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
    ChangeMolecules(&System, Sys_extra, change_beads, true);
    CheckSystem(System, struct_file_extra);
  }
  // use coordinate from a separate file (-c option)
  bool vtf_coor_var = false;
  int start_id = -1; // for lammpstrj file
  if (coor_file[0] != '\0') {
    // if (coor_type == LTRJ_FILE) {
    //   start_id = LtrjLowIndex(coor_file);
    // }
    int file_line_count = 0;
    FILE *fr = OpenFile(coor_file, "r");
    for (int i = 0; i < timestep; i++) {
      SkipTimestep(coor_type, fr, coor_file, struct_file, &file_line_count);
    }
    ReadTimestep(coor_type, fr, coor_file, &System, &file_line_count,
                 start_id, vtf_coor_var);
    fclose(fr);
  } //}}}

  // -def option (for vsf output file) //{{{
  bool *def_type = calloc(System.Count.BeadType, sizeof *def_type);
  if (BeadTypeOption(argc, argv, "-def", false, def_type, &System)) {
    exit(1);
  }
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
    WriteStructure(struct_type_out, struct_file_out, System,
                   vsf_def_type, lmp_mass);
  }

  FreeSystem(&System);
  if (struct_file_extra[0] != '\0') {
    FreeSystem(&Sys_extra);
  }

  return 0;
}
