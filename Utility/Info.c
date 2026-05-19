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
  .all = 22, // number of valid lines OptSpec (not counting last {NULL})
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
  {"<input>", NULL, "input structure file", OPT_ARG},
  {"-i", "<file> [type]", "secondary structure file (type: vtf/vsf/xyz/data/ltrj/field/itp/pdb)", OPT_EXTRA},
  {"-c", "<file>", "input coordinate file", OPT_EXTRA},
  {"--detailed", NULL, "use name, charge, mass, and radius to identfy bead types", OPT_EXTRA},
  {"-o", "<file>", "output structure file", OPT_EXTRA},
  {"--unique", NULL, "make all bead/molecule names unique", OPT_EXTRA},
  {"-def", "<bead name>", "default bead type (output vtf structure file only)", OPT_EXTRA},
  {"--mol", NULL, "make unbonded beads into molecules", OPT_EXTRA},
  {"--mass", NULL, "define lammps atom types by mass, but print per-atom charges in Atoms section (output lammps data file only)", OPT_EXTRA},
  {"-ebt", "<int>", "number of extra bead types (output lammps data file only)", OPT_EXTRA},
  {"--chbt", NULL, "change bead types using -i-provided file; molecules matched by name and bead count", OPT_EXTRA},
  {"--frag", NULL, "split disconnected molecules into fragments; single-bead fragments become unbonded beads", OPT_EXTRA},
  {"-lib", "<dir>", "library directory: rename bead types and print DPD interactions; appends interactions block to FIELD output", OPT_EXTRA},
  {"-sys", "<file>", "system_info file: assign names to molecule types before library renaming (required when types are unnamed)", OPT_EXTRA},
  {NULL}
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

int main(int argc, char *argv[]) {

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
  if (sys_in[0] != '\0')
    ReadSysInfo(sys_in, &System);
  if (lib_dir[0] != '\0') {
    if (sys_in[0] == '\0') {
      bool any_named = false;
      for (int i = 0; i < System.Count.MoleculeType; i++) {
        if (System.MoleculeType[i].Named) { any_named = true; break; }
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
    PruneSystem(&System, NULL);
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
        mt_new->Index = NULL;
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
    if (lib_used && opt.fout.type == FIELD_FILE && Count->BeadType > 0) {
      ArrNDd *pot = CreateArr3Dd(Count->BeadType, Count->BeadType, 3);
      FillPotFromLibrary(&lib, &System, pot);
      FILE *fw = OpenFile(opt.fout.name, "a");
      int n = Count->BeadType * (Count->BeadType - 1) / 2 + Count->BeadType;
      fprintf(fw, "interactions %d <a_ij> <r_c> <gamma>\n", n);
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = i; j < Count->BeadType; j++) {
          fprintf(fw, "%10s %10s dpd %lf %lf %lf\n",
                  System.BeadType[i].Name, System.BeadType[j].Name,
                  GetArr3D(pot, i, j, 0), GetArr3D(pot, i, j, 1),
                  GetArr3D(pot, i, j, 2));
        }
      }
      fclose(fw);
      FreeArrND(pot);
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
