#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Info analyzes the provided input structure file, printing system "
  "composition to standard output and, optionally, producing an output "
  "structure file of specified format (-o option). If some information "
  "required in the output file is missing, '???\' is printed instead. The "
  "system from the input file can be modified using a second structure file "
  "(-i option) and/or a coordinate file (-c option)",

  "Usage: Info <input> [options]",
  .args = 1, // number of mandatory arguments
  .all = 24, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  {"-ft", "<type>", "structure file type: vtf/vsf/xyz/data/ltrj/field/itp/pdb", OPT_COMMON},
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", nullptr, "input structure file", OPT_ARG},
  {"-i", "<file> [type]", "secondary input structure file", OPT_EXTRA},
  {"-c", "<file>", "input coordinate file", OPT_EXTRA},
  {"--detailed", nullptr, "use name, charge, mass, and radius "
    "to identfy bead types", OPT_EXTRA},
  {"-o", "<file>", "output structure file", OPT_EXTRA},
  {"--unique", nullptr, "make all bead/molecule names unique", OPT_EXTRA},
  {"-def", "<bead name>", "default bead type (output vtf structure file only)",
    OPT_EXTRA},
  {"--mol", nullptr, "make unbonded beads into molecules", OPT_EXTRA},
  {"--mass", nullptr, "define lammps atom types by mass, but print per-atom "
    "charges in Atoms section (output lammps data file only)", OPT_EXTRA},
  {"-ebt", "<int>", "number of extra bead types (output lammps data file only)",
    OPT_EXTRA},
  {"--chbt", nullptr, "change bead types using -i-provided file; "
    "molecules matched by name and bead count", OPT_EXTRA},
  {"--frag", nullptr, "split disconnected molecules into fragments; "
    "single-bead fragments become unbonded beads", OPT_EXTRA},
  {"-lib", "<dir>", "library directory: rename bead types and print DPD "
    "interactions; appends interactions block to FIELD output", OPT_EXTRA},
  {"-sys", "<file>", "system_info file: assign names to molecule types before "
    "library renaming (required when types are unnamed)", OPT_EXTRA},
  {"--check-library", nullptr, "check the library for bad names, dangling counterions, unknown type IDs, unbalanced charges and drifted row numbering, then exit (requires -lib; takes no input file)", OPT_EXTRA},
  {"--mol-info", "<mol>", "print one library molecule's role, mass, bead count, charge, z-extent and counterions, then exit (requires -lib; takes no input file)", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  int vsf_def,   // -def
      b_mol,     // --mole
      ebt;       // -ebt
  bool lmp_mass, // --mass
       detailed, // --detailed
       chbt,     // --chbt
       frag;     // --frag
  FILE_TYPE fout;          // -o
}; //}}}

// find root with path halving //{{{
static int uf_find(int *parent, int x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
} //}}}

/*
 * Answer a question about the library rather than about a system: bead counts,
 * mass, role, charge, z-extent and counterions of one molecule, as
 * 'key value' lines in the library's own grammar.
 *
 * This is what prep_simulation.sh used to work out by reading the library
 * itself, with hard-coded line numbers and a second copy of the format's rules.
 */
static void PrintMoleculeInfo(const char *lib_dir, const char *name) { //{{{
  LIBRARY lib = ReadLibrary(lib_dir);
  LIB_MOL mol;
  if (!ReadLibraryMolFile(lib_dir, name, &mol)) {
    if (snprintf(ERROR_MSG, LINE, "no library file for molecule %s%s%s",
                 ErrYellow(), name, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintErrorOption("--mol-info");
    exit(1);
  }
  // charge, ion count (sum of q^2) and z-extent of one molecule
  double charge = 0, ions = 0, zmin = 0, zmax = 0;
  for (int i = 0; i < mol.n_beads; i++) {
    int bt = FindBeadType(mol.bead_name[i], lib.System);
    if (bt == -1) {
      if (snprintf(ERROR_MSG, LINE, "bead type '%s%s%s' not in library",
                   ErrYellow(), mol.bead_name[i], ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    double q = lib.System.BeadType[bt].Charge;
    charge += q;
    ions += q * q;
    double z = mol.bead_pos[i].v[2];
    if (i == 0 || z < zmin) {
      zmin = z;
    }
    if (i == 0 || z > zmax) {
      zmax = z;
    }
  }
  const char *role = "soluble";
  if (mol.bilayer) {
    role = "bilayer";
  }
  printf("name           %s\n", mol.name);
  printf("role           %s\n", role);
  if (mol.has_M_w) {
    printf("M_w            %.4f\n", mol.M_w);
  }
  printf("n_beads        %d\n", mol.n_beads);
  printf("charge         %.4f\n", charge);
  printf("ions           %.4f\n", ions);
  printf("z_min          %.4f\n", zmin);
  printf("z_max          %.4f\n", zmax);
  // one line per counterion: name  count  n_beads  charge  ions
  for (int i = 0; i < mol.n_cion; i++) {
    LIB_MOL cion;
    if (!ReadLibraryMolFile(lib_dir, mol.cion[i].name, &cion)) {
      if (snprintf(ERROR_MSG, LINE, "no library file for counterion %s%s%s",
                   ErrYellow(), mol.cion[i].name, ErrRed()) < 0) {
        ErrorSnprintf();
      }
      PrintError();
      exit(1);
    }
    double cq = 0, cions = 0;
    for (int b = 0; b < cion.n_beads; b++) {
      int bt = FindBeadType(cion.bead_name[b], lib.System);
      if (bt == -1) {
        continue;
      }
      double q = lib.System.BeadType[bt].Charge;
      cq += q;
      cions += q * q;
    }
    printf("counterion     %s %d %d %.4f %.4f\n", mol.cion[i].name,
           mol.cion[i].count, cion.n_beads, cq, cions);
    FreeLibMol(&cion);
  }
  FreeLibMol(&mol);
  FreeLibrary(&lib);
} //}}}

// sum of bead charges over one molecule; ok=false if a bead type is unknown //{{{
static double MolCharge(const LIB_MOL *mol, const LIBRARY *lib, bool *ok) {
  double q = 0;
  for (int i = 0; i < mol->n_beads; i++) {
    int bt = FindBeadType(mol->bead_name[i], lib->System);
    if (bt == -1) {
      *ok = false;
      continue;
    }
    q += lib->System.BeadType[bt].Charge;
  }
  return q;
} //}}}
/*
 * The leading index on every row is decoration - position in the block is the
 * index and no parser reads the printed one - so nothing else would notice it
 * drifting out of step. Checking it here is what lets it stay a convenience
 * instead of becoming another thing to maintain by hand.
 */
static int CheckIndices(const char *path, const char *label) { //{{{
  FILE *fr = fopen(path, "r");
  if (!fr) {
    return 0;
  }
  char block[32] = "";
  int expect = 0, bad = 0;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (words == 0 || split[0][0] == '#') {
      continue;
    }
    if (words == 1) {
      s_strcpy(block, split[0], sizeof block);
      expect = 1;
      continue;
    }
    if (block[0] == '\0') {
      continue; // a scalar, before any block
    }
    long idx;
    if (!IsWholeNumber(split[0], &idx) || idx != expect) {
      printf("  %-24s %s row %d is numbered '%s'\n", label, block, expect,
             split[0]);
      bad++;
    }
    expect++;
  }
  fclose(fr);
  return bad;
} //}}}
/*
 * Check the library for the things loading it does not: names that disagree
 * with their filenames, counterions that point nowhere, IDs with no entry in
 * the type tables, molecules that do not balance, and drifted row numbering.
 *
 * A file that cannot be parsed at all is still fatal, here as anywhere else -
 * this reports what a library that loads can still be wrong about.
 */
static int CheckLibrary(const char *lib_dir) { //{{{
  LIBRARY lib = ReadLibrary(lib_dir);
  int problems = 0;
  int n_mols = 0, n_bare = 0;
  char (*names)[MOL_NAME] = nullptr;
  char (*files)[MOL_NAME] = nullptr;

  DIR *dir = opendir(lib_dir);
  if (!dir) {
    if (snprintf(ERROR_MSG, LINE, "cannot open library directory %s%s%s",
                 ErrYellow(), lib_dir, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  }
  printf("Checking %s\n\n", lib_dir);
  /*
   * Which species are somebody's counterion. A bare ion is charged by
   * definition and is neutralised by whatever names it, so it is the one thing
   * that must not be required to balance on its own - and being referenced is
   * what says so. Having no M_w says the same thing more weakly, and not every
   * library spells it that way: an example library may well give Cl its real
   * 35.45 while still only ever using it as a counterion.
   */
  char (*cion_of)[MOL_NAME] = nullptr;
  int n_cion_names = 0;
  struct dirent *ent;
  while ((ent = readdir(dir)) != nullptr) {
    size_t len = strlen(ent->d_name);
    if (len < 5 || strcmp(ent->d_name + len - 4, ".txt") != 0) {
      continue;
    }
    if (strncmp(ent->d_name, "list_", 5) == 0) {
      continue;
    }
    char stem[MOL_NAME];
    s_strcpy(stem, ent->d_name, MOL_NAME);
    stem[len - 4] = '\0';
    LIB_MOL mol;
    if (!ReadLibraryMolFile(lib_dir, stem, &mol)) {
      continue;
    }
    for (int c = 0; c < mol.n_cion; c++) {
      cion_of = s_realloc(cion_of, (n_cion_names + 1) * sizeof *cion_of);
      s_strcpy(cion_of[n_cion_names], mol.cion[c].name, MOL_NAME);
      n_cion_names++;
    }
    FreeLibMol(&mol);
  }
  rewinddir(dir);
  while ((ent = readdir(dir)) != nullptr) {
    size_t len = strlen(ent->d_name);
    if (len < 5 || strcmp(ent->d_name + len - 4, ".txt") != 0) {
      continue;
    }
    if (strncmp(ent->d_name, "list_", 5) == 0) {
      continue;
    }
    char stem[MOL_NAME];
    s_strcpy(stem, ent->d_name, MOL_NAME);
    stem[len - 4] = '\0';
    LIB_MOL mol;
    if (!ReadLibraryMolFile(lib_dir, stem, &mol)) {
      continue;
    }
    files = s_realloc(files, (n_mols + 1) * sizeof *files);
    names = s_realloc(names, (n_mols + 1) * sizeof *names);
    s_strcpy(files[n_mols], stem, MOL_NAME);
    s_strcpy(names[n_mols], mol.name, MOL_NAME);
    n_mols++;

    // the name inside the file is what the loader indexes by
    if (strcmp(mol.name, stem) != 0) {
      printf("  %-24s declares name '%s'\n", ent->d_name, mol.name);
      problems++;
    }
    // bead types must exist
    for (int i = 0; i < mol.n_beads; i++) {
      if (FindBeadType(mol.bead_name[i], lib.System) == -1) {
        printf("  %-24s bead %d is type '%s', not in list_parameters.txt\n",
               stem, i + 1, mol.bead_name[i]);
        problems++;
      }
    }
    // bond and angle IDs must exist
    for (int i = 0; i < mol.n_bonds; i++) {
      bool found = false;
      for (int k = 0; k < lib.n_bond_ids; k++) {
        if (strcmp(lib.bond_id[k].id, mol.bond_id[i]) == 0) {
          found = true;
          break;
        }
      }
      if (!found) {
        printf("  %-24s bond %d is type '%s', not in list_bonds.txt\n",
               stem, i + 1, mol.bond_id[i]);
        problems++;
      }
    }
    for (int i = 0; i < mol.n_angles; i++) {
      bool found = false;
      for (int k = 0; k < lib.n_angle_ids; k++) {
        if (strcmp(lib.angle_id[k].id, mol.angle_id[i]) == 0) {
          found = true;
          break;
        }
      }
      if (!found) {
        printf("  %-24s angle %d is type '%s', not in list_angles.txt\n",
               stem, i + 1, mol.angle_id[i]);
        problems++;
      }
    }
    /*
     * Counterions must resolve, and a species must balance with them. Only a
     * species that can be listed in input.txt has to: a bare ion like Cl is
     * charged by definition and is neutralised by whatever names it. Having an
     * M_w is exactly that distinction - it is the mass you would weigh out.
     */
    bool is_cion = false;
    for (int c = 0; c < n_cion_names; c++) {
      if (strcmp(cion_of[c], mol.name) == 0) {
        is_cion = true;
        break;
      }
    }
    if (is_cion || !mol.has_M_w) {
      n_bare++;
    }
    bool q_ok = !is_cion && mol.has_M_w;
    double q = MolCharge(&mol, &lib, &q_ok);
    for (int c = 0; c < mol.n_cion; c++) {
      LIB_MOL cion;
      if (!ReadLibraryMolFile(lib_dir, mol.cion[c].name, &cion)) {
        printf("  %-24s counterion '%s' has no library file\n",
               stem, mol.cion[c].name);
        problems++;
        q_ok = false;
        continue;
      }
      bool cion_ok = true;
      q += mol.cion[c].count * MolCharge(&cion, &lib, &cion_ok);
      if (!cion_ok) {
        q_ok = false;
      }
      FreeLibMol(&cion);
    }
    if (q_ok && fabs(q) > 1e-6) {
      printf("  %-24s net charge %+.2f with its counterions\n", stem, q);
      problems++;
    }
    char path[LINE];
    snprintf(path, LINE, "%s/%s", lib_dir, ent->d_name);
    problems += CheckIndices(path, stem);
    FreeLibMol(&mol);
  }
  closedir(dir);

  // names must be unique across the directory
  for (int i = 0; i < n_mols; i++) {
    for (int j = i + 1; j < n_mols; j++) {
      if (strcmp(names[i], names[j]) == 0) {
        printf("  %s.txt and %s.txt both declare the name '%s'\n",
               files[i], files[j], names[i]);
        problems++;
      }
    }
  }

  // interactions must name known bead types
  const char *tables[] = {"list_parameters.txt", "list_bonds.txt",
                          "list_angles.txt", "list_cross_interactions.txt"};
  for (int t = 0; t < 4; t++) {
    char path[LINE];
    snprintf(path, LINE, "%s/%s", lib_dir, tables[t]);
    problems += CheckIndices(path, tables[t]);
  }
  for (int i = 0; i < lib.n_inter; i++) {
    if (strcmp(lib.inter[i].name1, lib.inter[i].name2) == 0) {
      continue; // self-interactions come from the bead type table itself
    }
    if (FindBeadType(lib.inter[i].name1, lib.System) == -1 ||
        FindBeadType(lib.inter[i].name2, lib.System) == -1) {
      printf("  list_cross_interactions.txt  '%s'-'%s' names an unknown bead "
             "type\n", lib.inter[i].name1, lib.inter[i].name2);
      problems++;
    }
  }

  printf("\n%d molecules (%d of them only ever another species' counterion), "
         "%d bead types,\n           %d bond types, %d angle types, "
         "%d pair interactions\n",
         n_mols, n_bare, lib.System.Count.BeadType, lib.System.Count.BondType,
         lib.System.Count.AngleType, lib.n_inter);
  if (problems == 0) {
    printf("no problems found\n");
  } else {
    printf("%d problem(s)\n", problems);
  }
  free(names);
  free(files);
  free(cion_of);
  FreeLibrary(&lib);
  return problems;
} //}}}

int main(int argc, char *argv[]) {

  /*
   * --mol-info asks about the library, not about a system, so it needs no
   * input structure file and runs before OptionCheck() enforces one.
   */
  { //{{{
    const char *lib_q = nullptr, *mol_q = nullptr;
    bool check_q = false;
    for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "--check-library") == 0) {
        check_q = true;
      }
      if ((i + 1) >= argc) {
        continue;
      }
      if (strcmp(argv[i], "-lib") == 0) {
        lib_q = argv[i + 1];
      }
      if (strcmp(argv[i], "--mol-info") == 0) {
        mol_q = argv[i + 1];
      }
    }
    if (mol_q != nullptr || check_q) {
      if (lib_q == nullptr) {
        s_strcpy(ERROR_MSG, "--mol-info and --check-library require -lib", LINE);
        PrintError();
        exit(1);
      }
      if (check_q) {
        if (CheckLibrary(lib_q) > 0) {
          return 1;
        }
        return 0;
      }
      PrintMoleculeInfo(lib_q, mol_q);
      return 0;
    }
  } //}}}


  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, argv[++count], LINE);
  // -ft option: override structure file type without extension-based detection
  char ft_str[LINE];
  if (FileOption(argc, argv, COMMON_OPTS[C_FT].opt, ft_str)) {
    in.stru.type = FileTypeFromString(ft_str);
    if (in.stru.type != VTF_FILE && in.stru.type != VSF_FILE &&
        in.stru.type != FIELD_FILE && in.stru.type != LDATA_FILE &&
        in.stru.type != LTRJ_FILE && in.stru.type != XYZ_FILE &&
        in.stru.type != ITP_FILE && in.stru.type != PDB_FILE) {
      if (snprintf(ERROR_MSG, LINE, "not a structure file type '%s'",
                   ft_str) < 0) {
        ErrorSnprintf();
      }
      PrintErrorOption(COMMON_OPTS[C_FT].opt);
      exit(1);
    }
  } else {
    in.stru.type = StructureFileType(in.stru.name);
  }
  // -i option with optional type string as 2nd argument //{{{
  SYS_FILES extra = InitSysFiles;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG,
                 "missing file name (or file name begins with a dash)", LINE);
        PrintErrorOption("-i");
        exit(1);
      }
      s_strcpy(extra.stru.name, argv[i+1], LINE);
      if ((i+2) < argc && argv[i+2][0] != '-') {
        int ft = FileTypeFromString(argv[i+2]);
        if (ft == VTF_FILE || ft == VSF_FILE || ft == FIELD_FILE ||
            ft == LDATA_FILE || ft == LTRJ_FILE || ft == XYZ_FILE ||
            ft == ITP_FILE || ft == PDB_FILE) {
          extra.stru.type = ft;
        } else if (ft != -1) {
          if (snprintf(ERROR_MSG, LINE, "not a structure file type '%s'",
                       argv[i+2]) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption("-i");
          exit(1);
        } else {
          if (snprintf(ERROR_MSG, LINE, "unknown structure file type: '%s%s%s'",
                       ErrYellow(), argv[i+2], ErrRed()) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption("-i");
          exit(1);
        }
      } else {
        extra.stru.type = StructureFileType(extra.stru.name);
      }
      break;
    }
  } //}}}
  // input coordinate file (-c option) //{{{
  FileOption(argc, argv, "-c", in.coor.name);
  if (in.coor.name[0] != '\0') {
    in.coor.type = CoordinateFileType(in.coor.name);
    extra.coor = in.coor;
  } //}}}
  // output file (-o option) //{{{
  opt.fout = InitFile;
  if (FileOption(argc, argv, "-o", opt.fout.name)) {
    opt.fout.type = FileType(opt.fout.name);
    if (opt.fout.type == ITP_FILE || opt.fout.type == PDB_FILE) {
      const char *fmt_name = "pdb";
      if (opt.fout.type == ITP_FILE) {
        fmt_name = "itp";
      }
      snprintf(ERROR_MSG, LINE, "writing to %s%s%s format is not supported",
               ErrRed(), fmt_name, ErrYellow());
      PrintError();
      exit(1);
    }
  } //}}}
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // extra bead types for data output (-ebt option)
  opt.ebt = 0;
  OneNumberOption(argc, argv, "-ebt", &opt.ebt, 'i');
  // use mass only for atom type definition for data output file
  opt.lmp_mass = BoolOption(argc, argv, "--mass");
  // make unbonded beads into molecules (vtf output only)
  opt.b_mol = BoolOption(argc, argv, "--mol");
  // base bead types on name, charge, mass, and radius (vtf input file)
  opt.detailed = BoolOption(argc, argv, "--detailed");
  // change bead types using secondary structure file
  opt.chbt = BoolOption(argc, argv, "--chbt");
  opt.frag = BoolOption(argc, argv, "--frag"); //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // read information from input file(s) //{{{
  SYSTEM System = ReadStructure(in, opt.detailed);
  // if (in.stru.type == PDB_FILE) {
  //   FreeSystem(&System);
  //   free(opt);
  //   return 0;
  // }
  COUNT *Count = &System.Count;
  // apply -sys and -lib to rename bead types //{{{
  char sys_in[LINE] = "", lib_dir[LINE] = "";
  LIBRARY lib = {0};
  bool lib_used = false;
  FileOption(argc, argv, "-sys", sys_in);
  FileOption(argc, argv, "-lib", lib_dir);
  if (sys_in[0] != '\0') {
    ReadSysInfo(sys_in, &System);
  }
  if (lib_dir[0] != '\0') {
    if (sys_in[0] == '\0') {
      bool any_named = false;
      for (int i = 0; i < System.Count.MoleculeType; i++) {
        if (System.MoleculeType[i].Named) {
          any_named = true;
          break;
        }
      }
      if (!any_named) {
        s_strcpy(ERROR_MSG, "-lib without -sys: molecule types are unnamed; "
                 "bead type renaming will be skipped", LINE);
        PrintWarning();
      }
    }
    lib = ReadLibrary(lib_dir);
    lib_used = true;
    RenameBeadTypesFromLibrary(&System, &lib, lib_dir);
  } //}}}
  // use coordinate from a separate file (-c option)
  if (in.coor.type != -1) {
    int line_count = 0;
    FILE *fr = OpenFile(in.coor.name, "r");
    for (int i = 1; i < commons.start; i++) { // from 1 as timestep=1 is the first
      SkipTimestep(in, fr, &line_count);
    }
    ReadTimestep(in, fr, &System, &line_count);
    fclose(fr);
  } else {
    // all beads are in the timestep
    Count->BeadCoor = Count->Bead;
    for (int i = 0; i < Count->Bead; i++) {
      System.Bead[i].InTimestep = true;
      System.BeadCoor[i] = i;
    }
  }
  // print initial system information only if extra file(s) are present
  if (commons.verbose && (extra.stru.type != -1 || in.coor.type != -1)) {
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
    Sys_extra = ReadStructure(extra, opt.detailed);
    if (commons.verbose) {
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
    // exchanged bead types using extra file (--chbt option)
    if (opt.chbt) { //{{{
      for (int i = 0; i < Count->MoleculeType; i++) {
        MOLECULETYPE *mt = &System.MoleculeType[i];
        if (mt->Number == 0) {
          continue;
        }
        int type_e = FindMoleculeName(mt->Name, Sys_extra);
        if (type_e == -1) {
          continue;
        }
        MOLECULETYPE *mt_e = &Sys_extra.MoleculeType[type_e];
        if (mt->nBeads != mt_e->nBeads) {
          if (snprintf(ERROR_MSG, LINE,
                       "bead count mismatch for molecule %s%s%s "
                       "(%s%d%s vs %s%d%s); ignoring this molecule type",
                       ErrYellow(), mt->Name, ErrCyan(),
                       ErrYellow(), mt->nBeads, ErrCyan(),
                       ErrYellow(), mt_e->nBeads, ErrCyan()) < 0) {
            ErrorSnprintf();
          }
          PrintWarnOption("--chbt");
          continue;
        }
        for (int j = 0; j < mt->nBeads; j++) {
          char *extra_name = Sys_extra.BeadType[mt_e->Bead[j]].Name;
          int new_type = FindBeadType(extra_name, System);
          if (new_type == -1) {
            BEADTYPE *bt_e = &Sys_extra.BeadType[mt_e->Bead[j]];
            new_type = Count->BeadType;
            NewBeadType(&System.BeadType, &Count->BeadType,
                        extra_name, bt_e->Charge, bt_e->Mass, bt_e->Radius);
          }
          mt->Bead[j] = new_type;
          for (int k = 0; k < mt->Number; k++) {
            int mol_id = mt->Index[k];
            System.Bead[System.Molecule[mol_id].Bead[j]].Type = new_type;
          }
        }
      }
    } //}}}
    ChangeMolecules(&System, Sys_extra, false);
    CheckSystem(System, extra.stru.name);
  } //}}}

  // -def option (for vsf output file) //{{{
  bool *def_type = calloc(Count->BeadType, sizeof *def_type);
  if (!def_type) {
    ErrorAlloc("def_type");
  }
  TypeOption(argc, argv, "-def", 'b', false, def_type, System);
  opt.vsf_def = -1;
  for (int i = 0; i < Count->BeadType; i++) {
    if (def_type[i]) {
      opt.vsf_def = i;
      break;
    }
  }
  free(def_type); //}}}

  if (Count->Bead > 0) {
    PruneSystem(&System, nullptr);
  }

  // split disconnected molecules into fragments (--frag option) //{{{
  if (opt.frag) {
    int orig_nmol = Count->Molecule;
    bool *removed = calloc(orig_nmol, sizeof *removed);
    if (!removed) {
      ErrorAlloc("removed");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      if (mt->Number == 0 || mt->nBeads <= 1) {
        continue;
      }
      // build union-find over bead positions 0..nBeads-1
      int *parent = malloc(mt->nBeads * sizeof *parent);
      if (!parent) {
        ErrorAlloc("parent");
      }
      for (int j = 0; j < mt->nBeads; j++) {
        parent[j] = j;
      }
      for (int b = 0; b < mt->nBonds; b++) {
        int ra = uf_find(parent, mt->Bond[b][0]);
        int rb = uf_find(parent, mt->Bond[b][1]);
        if (ra != rb) {
          parent[ra] = rb;
        }
      }
      // assign component IDs
      int *comp = malloc(mt->nBeads * sizeof *comp);
      int *root_to_comp = malloc(mt->nBeads * sizeof *root_to_comp);
      if (!comp || !root_to_comp) {
        ErrorAlloc("comp/root_to_comp");
      }
      InitIntArray(root_to_comp, mt->nBeads, -1);
      int n_comp = 0;
      for (int j = 0; j < mt->nBeads; j++) {
        int r = uf_find(parent, j);
        if (root_to_comp[r] == -1) {
          root_to_comp[r] = n_comp++;
        }
        comp[j] = root_to_comp[r];
      }
      free(root_to_comp);
      free(parent);
      if (n_comp == 1) {
        free(comp);
        continue;
      }
      // count beads/bonds/angles/dihedrals/impropers per component
      int *comp_nbeads = calloc(n_comp, sizeof *comp_nbeads);
      int *comp_nbonds = calloc(n_comp, sizeof *comp_nbonds);
      int *comp_nangles = calloc(n_comp, sizeof *comp_nangles);
      int *comp_ndihed = calloc(n_comp, sizeof *comp_ndihed);
      int *comp_nimpro = calloc(n_comp, sizeof *comp_nimpro);
      if (!comp_nbeads || !comp_nbonds || !comp_nangles ||
          !comp_ndihed || !comp_nimpro) {
        ErrorAlloc("comp_n*");
      }
      for (int j = 0; j < mt->nBeads; j++) {
        comp_nbeads[comp[j]]++;
      }
      for (int b = 0; b < mt->nBonds; b++) {
        comp_nbonds[comp[mt->Bond[b][0]]]++;
      }
      for (int a = 0; a < mt->nAngles; a++) {
        comp_nangles[comp[mt->Angle[a][0]]]++;
      }
      for (int d = 0; d < mt->nDihedrals; d++) {
        comp_ndihed[comp[mt->Dihedral[d][0]]]++;
      }
      for (int d = 0; d < mt->nImpropers; d++) {
        comp_nimpro[comp[mt->Improper[d][0]]]++;
      }
      // old id -> new id
      int *old_to_new = malloc(mt->nBeads * sizeof *old_to_new);
      int *comp_cursor = calloc(n_comp, sizeof *comp_cursor);
      if (!old_to_new || !comp_cursor) {
        ErrorAlloc("old_to_new/comp_cursor");
      }
      for (int j = 0; j < mt->nBeads; j++) {
        old_to_new[j] = comp_cursor[comp[j]]++;
      }
      free(comp_cursor);
      // create new molecule types for each fragment
      int *comp_type = malloc(n_comp * sizeof *comp_type);
      if (!comp_type) {
        ErrorAlloc("comp_type");
      }
      int frag_num = 1;
      for (int c = 0; c < n_comp; c++) {
        if (comp_nbeads[c] <= 1) {
          comp_type[c] = -1;
          continue;
        }
        char new_name[MOL_NAME];
        if (snprintf(new_name, MOL_NAME, "%s_%d", mt->Name, frag_num) < 0) {
          ErrorSnprintf();
        }
        frag_num++;
        comp_type[c] = Count->MoleculeType;
        NewMolType(&System.MoleculeType, &Count->MoleculeType, new_name,
                   comp_nbeads[c], comp_nbonds[c], comp_nangles[c],
                   comp_ndihed[c], comp_nimpro[c]);
        mt = &System.MoleculeType[i]; // re-derive: MoleculeType was realloc'd
        MOLECULETYPE *mt_new = &System.MoleculeType[comp_type[c]];
        mt_new->Number = 0;
        mt_new->Index = nullptr;
        mt_new->InVcf = false;
        mt_new->Named = true;
        // fill Bead[]
        int pos = 0;
        for (int j = 0; j < mt->nBeads; j++) {
          if (comp[j] == c) {
            mt_new->Bead[pos] = mt->Bead[j];
            pos++;
          }
        }
        // fill Bond[] (renumber local indices)
        pos = 0;
        for (int b = 0; b < mt->nBonds; b++) {
          if (comp[mt->Bond[b][0]] == c) {
            mt_new->Bond[pos][0] = old_to_new[mt->Bond[b][0]];
            mt_new->Bond[pos][1] = old_to_new[mt->Bond[b][1]];
            mt_new->Bond[pos][2] = mt->Bond[b][2];
            pos++;
          }
        }
        // fill Angle[]
        pos = 0;
        for (int a = 0; a < mt->nAngles; a++) {
          if (comp[mt->Angle[a][0]] == c) {
            mt_new->Angle[pos][0] = old_to_new[mt->Angle[a][0]];
            mt_new->Angle[pos][1] = old_to_new[mt->Angle[a][1]];
            mt_new->Angle[pos][2] = old_to_new[mt->Angle[a][2]];
            mt_new->Angle[pos][3] = mt->Angle[a][3];
            pos++;
          }
        }
        // fill Dihedral[]
        pos = 0;
        for (int d = 0; d < mt->nDihedrals; d++) {
          if (comp[mt->Dihedral[d][0]] == c) {
            mt_new->Dihedral[pos][0] = old_to_new[mt->Dihedral[d][0]];
            mt_new->Dihedral[pos][1] = old_to_new[mt->Dihedral[d][1]];
            mt_new->Dihedral[pos][2] = old_to_new[mt->Dihedral[d][2]];
            mt_new->Dihedral[pos][3] = old_to_new[mt->Dihedral[d][3]];
            mt_new->Dihedral[pos][4] = mt->Dihedral[d][4];
            pos++;
          }
        }
        // fill Improper[]
        pos = 0;
        for (int d = 0; d < mt->nImpropers; d++) {
          if (comp[mt->Improper[d][0]] == c) {
            mt_new->Improper[pos][0] = old_to_new[mt->Improper[d][0]];
            mt_new->Improper[pos][1] = old_to_new[mt->Improper[d][1]];
            mt_new->Improper[pos][2] = old_to_new[mt->Improper[d][2]];
            mt_new->Improper[pos][3] = old_to_new[mt->Improper[d][3]];
            mt_new->Improper[pos][4] = mt->Improper[d][4];
            pos++;
          }
        }
      }
      // split each molecule instance into component molecules
      int n_mol_orig = mt->Number;
      for (int k = 0; k < n_mol_orig; k++) {
        int mol_id = mt->Index[k];
        removed[mol_id] = true;
        for (int c = 0; c < n_comp; c++) {
          if (comp_type[c] == -1) { // single bead -> unbonded (no molecule)
            for (int j = 0; j < mt->nBeads; j++) {
              if (comp[j] == c) {
                System.Bead[System.Molecule[mol_id].Bead[j]].Molecule = -1;
                break;
              }
            }
          } else { // multiple beads -> molecule
            int new_mol_id = Count->Molecule++;
            System.Molecule = s_realloc(System.Molecule,
                                        Count->Molecule * sizeof *System.Molecule);
            mt = &System.MoleculeType[i]; // re-derive after potential Molecule move
            MOLECULE *mol_new = &System.Molecule[new_mol_id];
            mol_new->Type = comp_type[c];
            mol_new->Bead = malloc(comp_nbeads[c] * sizeof *mol_new->Bead);
            if (!mol_new->Bead) {
              ErrorAlloc("mol_new->Bead");
            }
            mol_new->Index = Count->HighestResid + 1;
            Count->HighestResid++;
            mol_new->InTimestep = System.Molecule[mol_id].InTimestep;
            mol_new->Aggregate = -1;
            int pos = 0;
            for (int j = 0; j < mt->nBeads; j++) {
              if (comp[j] == c) {
                int bead_id = System.Molecule[mol_id].Bead[j];
                mol_new->Bead[pos++] = bead_id;
                System.Bead[bead_id].Molecule = new_mol_id;
              }
            }
            System.MoleculeType[comp_type[c]].Number++;
          }
        }
      }
      mt->Number = 0;
      free(comp);
      free(comp_nbeads);
      free(comp_nbonds);
      free(comp_nangles);
      free(comp_ndihed);
      free(comp_nimpro);
      free(old_to_new);
      free(comp_type);
    }
    // compact Molecule[] array: remove old fragmented instances
    int *remap_mol = malloc(Count->Molecule * sizeof *remap_mol);
    if (!remap_mol) {
      ErrorAlloc("remap_mol");
    }
    int new_mol_count = 0;
    for (int m = 0; m < Count->Molecule; m++) {
      if (m < orig_nmol && removed[m]) {
        free(System.Molecule[m].Bead);
        remap_mol[m] = -1;
      } else {
        if (new_mol_count != m) {
          System.Molecule[new_mol_count] = System.Molecule[m];
        }
        remap_mol[m] = new_mol_count++;
      }
    }
    Count->Molecule = new_mol_count;
    for (int b = 0; b < Count->Bead; b++) {
      if (System.Bead[b].Molecule != -1) {
        System.Bead[b].Molecule = remap_mol[System.Bead[b].Molecule];
      }
    }
    free(removed);
    free(remap_mol);
    if (Count->Molecule > 0) {
      System.Molecule = s_realloc(System.Molecule,
                                  Count->Molecule * sizeof *System.Molecule);
    }
    // compact MoleculeType[] array: remove zero-count (fragmented) types
    int *remap_type = malloc(Count->MoleculeType * sizeof *remap_type);
    if (!remap_type) {
      ErrorAlloc("remap_type");
    }
    int new_type_count = 0;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      if (mt->Number == 0) {
        free(mt->Index);
        FreeMoleculeTypeEssentials(mt);
        if (mt->nBTypes > 0) {
          free(mt->BType);
        }
        remap_type[i] = -1;
      } else {
        if (new_type_count != i) {
          System.MoleculeType[new_type_count] = System.MoleculeType[i];
        }
        remap_type[i] = new_type_count++;
      }
    }
    Count->MoleculeType = new_type_count;
    for (int m = 0; m < Count->Molecule; m++) {
      System.Molecule[m].Type = remap_type[System.Molecule[m].Type];
    }
    free(remap_type);
    if (Count->MoleculeType > 0) {
      System.MoleculeType = s_realloc(System.MoleculeType,
                                      Count->MoleculeType * sizeof *System.MoleculeType);
    }
    // update bonded/unbonded counts then refill arrays
    Count->Bonded = 0;
    Count->Unbonded = 0;
    for (int b = 0; b < Count->Bead; b++) {
      if (System.Bead[b].Molecule != -1) {
        Count->Bonded++;
      } else {
        Count->Unbonded++;
      }
    }
    Count->BondedCoor = Count->Bonded;
    Count->UnbondedCoor = Count->Unbonded;
    FillBondedUnbonded(&System);
    for (int i = 0; i < Count->MoleculeType; i++) {
      FillMoleculeTypeBType(&System.MoleculeType[i]);
      FillMoleculeTypeChargeMass(&System.MoleculeType[i], System.BeadType);
      SortAll(&System.MoleculeType[i]);
    }
    ReFillMoleculeTypeIndex(&System);
    CountBondAngleDihedralImproper(&System);
    if (Count->Molecule > 0) {
      System.MoleculeCoor = s_realloc(System.MoleculeCoor,
                                      Count->Molecule * sizeof *System.MoleculeCoor);
    }
    Count->MoleculeCoor = 0;
    for (int m = 0; m < Count->Molecule; m++) {
      if (System.Molecule[m].InTimestep) {
        System.MoleculeCoor[Count->MoleculeCoor++] = m;
      }
    }
  } //}}}

  // make unbonded beads into molecules //{{{
  if (opt.b_mol) {
    // first new molid is the one higher than the last one
    int mol_id = Count->Molecule;
    // first new resid is one higher than the old one, so +1
    int resid = Count->HighestResid + 1;
    Count->HighestResid += Count->Unbonded;
    Count->Molecule += Count->Unbonded;
    System.Molecule = s_realloc(System.Molecule,
                                Count->Molecule * sizeof *System.Molecule);
    for (int i = 0; i < Count->BeadType; i++) {
      count = System.BeadType[i].Number;
      bool first = false;
      for (int j = 0; j < count; j++) {
        int id = System.BeadType[i].Index[j];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        if (b->Molecule == -1) {
          int n_mt = Count->MoleculeType;
          MOLECULE *mol = &System.Molecule[mol_id];
          mol->Bead = calloc(1, sizeof *mol->Bead);
          if (!mol->Bead) {
            ErrorAlloc("mol->Bead");
          }
          if (!first) {
            NewMolType(&System.MoleculeType, &Count->MoleculeType,
                       bt->Name, 1, 0, 0, 0, 0);
            MOLECULETYPE *mt = &System.MoleculeType[n_mt];
            mt->Number = 1;
            mt->Mass = bt->Mass;
            mt->Charge = bt->Charge;
            mt->Bead[0] = b->Type;
            mt->nBTypes = 1;
            mt->BType = malloc(sizeof *mt->BType);
            if (!mt->BType) {
              ErrorAlloc("mt->BType");
            }
            mt->BType[0] = b->Type;
            mol->Type = n_mt;
            mt->Index = malloc(sizeof *mt->Index);
            if (!mt->Index) {
              ErrorAlloc("mt->Index");
            }
            first = true;
          } else {
            MOLECULETYPE *mt = &System.MoleculeType[n_mt-1];
            mt->Number++;
            mol->Type = n_mt - 1;
          }
          b->Molecule = mol_id;
          mol->Index = resid;
          mol->Bead[0] = id;
          mol->Aggregate = -1;
          mol->InTimestep = true;
          resid++;
          mol_id++;
        }
      }
    }
    Count->Bonded += Count->Unbonded;
    Count->Unbonded = 0;
    FillBondedUnbonded(&System);
    // fill all molecule type indices
    ReFillMoleculeTypeIndex(&System);
  } //}}}

  // print the system information //{{{
  fprintf(stdout, "\n==================================================");
  if (extra.stru.type != -1 || in.coor.type != -1) {
    printf("\nFinal system composition");
  } else {
    printf("\nSystem composition");
  }
  fprintf(stdout, "\n==================================================\n");
  VerboseOutput(System);
  if (commons.verbose) { // -v option
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecules(System);
  }
  if (lib_used && Count->BeadType > 0) {
    ArrNDd *pot = CreateArr3Dd(Count->BeadType, Count->BeadType, 3);
    FillPotFromLibrary(&lib, &System, pot);
    fprintf(stdout, "\nDPD interactions:\n");
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = i; j < Count->BeadType; j++) {
        fprintf(stdout, "  %10s %10s dpd %g %g %g\n",
                System.BeadType[i].Name, System.BeadType[j].Name,
                GetArr3D(pot, i, j, 0), GetArr3D(pot, i, j, 1),
                GetArr3D(pot, i, j, 2));
      }
    }
    FreeArrND(pot);
  } //}}}

  // write the output file if required (-o option) //{{{
  if (opt.fout.name[0] != '\0') {
    if (opt.fout.type == LDATA_FILE && opt.ebt > 0) {
      NewBeadType(&System.BeadType, &Count->BeadType, "extra", 0, 1, 1);
    }
    bool *write = malloc(sizeof *write * Count->Bead);
    if (!write) {
      ErrorAlloc("write");
    }
    InitBoolArray(write, Count->Bead, true);
    WriteOutput(System, write, opt.fout, opt.lmp_mass, opt.vsf_def,
                argc, argv);
    if (lib_used && opt.fout.type == FIELD_FILE) {
      AppendFieldInteractions(opt.fout.name, &System, &lib);
    }
    free(write);
  } //}}}

  FreeSystem(&System);
  if (extra.stru.name[0] != '\0') {
    FreeSystem(&Sys_extra);
  }
  if (lib_used) {
    FreeLibrary(&lib);
  }

  return 0;
}
