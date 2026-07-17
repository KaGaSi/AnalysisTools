#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "BondLength utility calculates distribution of bond lengths in specified "
  "molecule type(s) for all bonds, printing results per bond type or adding "
  "all molecules' bonds as well (--all option). "
  "Note that input structure file with defined bonds must be used. "
  "The utility can also calculate distribution of "
  "distances between any two beads in those molecule types (-n option).",

  "Usage: BondLength <input> <width> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 18, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<width>", NULL, "width of a distribution bin", OPT_ARG},
  {"<output>", NULL, "output file with the distribution", OPT_ARG},
  {"-mt", "<name(s)>", "molecule types to use (default: all)", OPT_EXTRA},
  {"--joined", NULL, "<input> contains joined coordinates", OPT_EXTRA},
  {"--all", NULL, "calculate distribution for each bond in the molecule types",
    OPT_EXTRA},
  {"-n", "<file> [ints]", "distribution of distances between specified "
    "bead pair(s) (default [ints]: first and last bead)", OPT_EXTRA},
  {"-w", "<float>", "warn if the length exceeds <float>", OPT_EXTRA},
  {"-t", "<file>", "save per-timestep data to <file>", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  bool join,         // --joined
       *mt,          // -mt
       all;          // --all
  int n_list[100],   // -n (list of bead id pairs)
      n_number;      // -n (total number of beads in the pairs)
  double warn;       // -w option
  char n_file[LINE], // -n (output file)
       t_file[LINE]; // -t (output file)
}; //}}}

// helper functions //{{{
// write mins, maxes, and averages
static void WriteMinsMaxesAvgs(FILE *fw, SYSTEM System, OPT opt,
                               ArrNDd *bond_bt_mma, ArrNDi *bond_bt_norm,
                               ArrNDd *bond_all_mma, ArrNDi *bond_all_norm);
// are the provided bead types bonded?
static bool BondedPair(MOLECULETYPE mt, int btype[2]); //}}}

// per-timestep calculation and output //{{{
void Calculation(SYSTEM *System, STEP step, OPT opt,
                 ArrNDd *bond_bt, ArrNDd *bond_bt_mma,
                 ArrNDd *bond_all, ArrNDd *bond_all_mma,
                 ArrNDd *bond_n, ArrNDd *bond_n_mma,
                 double width, int bins, int n_per_set) {
  COUNT *Count = &System->Count;
  WrapJoinCoordinates(System, true, opt.join);
  // per-timestep arrays //{{{
  ArrNDd *per_step = CreateArr3Dd(Count->MoleculeType,
                                  Count->BeadType, Count->BeadType);
  ArrNDi *per_step_c = CreateArr3Di(Count->MoleculeType,
                                    Count->BeadType, Count->BeadType);
  if (!per_step || !per_step_c) {
    ErrorAlloc("step/step_c");
  }
  ArrNDd *step_t = NULL;
  ArrNDi *step_t_c = NULL;
  int n_pair_num = opt.n_number / n_per_set;
  if (opt.t_file[0] != '\0' || opt.n_file[0] != '\0') {
    if (!(step_t = CreateArr2Dd(Count->MoleculeType, n_pair_num)) ||
        !(step_t_c = CreateArr2Di(Count->MoleculeType, n_pair_num))) {
      ErrorAlloc("step_t/step_t_c");
    }
  }
  ArrNDd *step_bond_all = NULL;
  ArrNDi *step_bond_all_c = NULL;
  if (opt.t_file[0] != '\0' && opt.all) {
    int max_b = 0;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (System->MoleculeType[i].nBonds > max_b) {
        max_b = System->MoleculeType[i].nBonds;
      }
    }
    if (!(step_bond_all = CreateArr2Dd(Count->MoleculeType, max_b)) ||
        !(step_bond_all_c = CreateArr2Di(Count->MoleculeType, max_b))) {
      ErrorAlloc("step_bond_all/step_bond_all_c");
    }
  } //}}}
  // calculate bond lengths //{{{
  // go through all molecules
  // TODO: make into for (mtype); for (mtype.index)
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol_i = &System->Molecule[i];
    MOLECULETYPE *mt_i = &System->MoleculeType[mol_i->Type];
    if (opt.mt[mol_i->Type]) { // use only specified molecule types
      for (int j = 0; j < mt_i->nBonds; j++) {
        // bead ids in the bond
        int id1 = mol_i->Bead[mt_i->Bond[j][0]],
            id2 = mol_i->Bead[mt_i->Bond[j][1]];
        BEAD *b_1 = &System->Bead[id1],
             *b_2 = &System->Bead[id2];
        // bond length
        vec3d bond = Vector(b_1->Position, b_2->Position);
        bond.v[0] = VectLength(bond);
        // warn if bond is too long //{{{
        if (opt.warn != HIGHNUM && bond.v[0] > opt.warn) {
          snprintf(ERROR_MSG, LINE, "bond between "
                   "beads %s%d%s and %s%d%s is too long (%s%lf%s)",
                   ErrYellow(), id1, ErrCyan(), ErrYellow(), id2, ErrCyan(),
                   ErrYellow(), bond.v[0], ErrCyan());
          PrintWarnOption("-w");
        } //}}}
        // btype1 must be lower than btype2
        int *id_lo, *id_hi;
        if (b_1->Type < b_2->Type) {
          id_lo = &b_1->Type;
          id_hi = &b_2->Type;
        } else {
          id_lo = &b_2->Type;
          id_hi = &b_1->Type;
        }
        AddArr3D(per_step, mol_i->Type, *id_lo, *id_hi, bond.v[0]);
        AddArr3D(per_step_c, mol_i->Type, *id_lo, *id_hi, 1);
        // mins & maxes & averages //{{{
        size_t id_0[4] = {mol_i->Type, *id_lo, *id_hi, 0};
        size_t id_1[4] = {mol_i->Type, *id_lo, *id_hi, 1};
        size_t id_2[4] = {mol_i->Type, *id_lo, *id_hi, 2};
        if (bond.v[0] < GetArrND(bond_bt_mma, id_0)) {
          SetArrND(bond_bt_mma, id_0, bond.v[0]);
        }
        if (bond.v[0] > GetArrND(bond_bt_mma, id_1)) {
          SetArrND(bond_bt_mma, id_1, bond.v[0]);
        }
        AddArrND(bond_bt_mma, id_2, bond.v[0]);
        if (opt.all) {
          if (bond.v[0] < GetArr3D(bond_all_mma, mol_i->Type, j, 0 )) {
            SetArr3D(bond_all_mma, mol_i->Type, j, 0, bond.v[0]);
          }
          if (bond.v[0] > GetArr3D(bond_all_mma, mol_i->Type, j, 1)) {
            SetArr3D(bond_all_mma, mol_i->Type, j, 1, bond.v[0]);
          }
          AddArr3D(bond_all_mma, mol_i->Type, j, 2, bond.v[0]);
          if (opt.t_file[0] != '\0') {
            AddArr2D(step_bond_all, mol_i->Type, j, bond.v[0]);
            AddArr2D(step_bond_all_c, mol_i->Type, j, 1);
          }
        }
        //}}}
        int k = bond.v[0] / width;
        if (k < bins) {
          size_t id[4] = {mol_i->Type, *id_lo, *id_hi, k};
          AddArrND(bond_bt, id, 1);
          if (opt.all) {
            AddArr3D(bond_all, mol_i->Type, j, k, 1);
          }
        }
      }
    }
  } //}}}
  // calculate distance (-n option) //{{{
  if (opt.n_file[0] != '\0') {
    for (int i = 0; i < Count->Molecule; i++) {
      MOLECULE *mol_i = &System->Molecule[i];
      MOLECULETYPE *mt_i = &System->MoleculeType[mol_i->Type];
      if (opt.mt[mol_i->Type]) { // use only specified molecule types
        for (int j = 0; j < opt.n_number; j += n_per_set) {
          // bead ids to use //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt.n_list[j] >= mt_i->nBeads) {
            id1 = mol_i->Bead[mt_i->nBeads-1];
          } else { // use specified index otherwise
            id1 = mol_i->Bead[opt.n_list[j]-1];
          }
          // use last molecule bead if bead index too high or -1
          if (opt.n_list[j+1] >= mt_i->nBeads) {
            id2 = mol_i->Bead[mt_i->nBeads-1];
          } else { // use specified index otherwise
            id2 = mol_i->Bead[opt.n_list[j+1]-1];
          } //}}}
          BEAD *b_1 = &System->Bead[id1],
               *b_2 = &System->Bead[id2];
          // distance calculation
          vec3d r12 = Vector(b_1->Position, b_2->Position);
          double dist = VectLength(r12);
          AddArr2D(step_t, mol_i->Type, j / 2, dist);
          AddArr2D(step_t_c, mol_i->Type, j / 2, 1);
          // distance mins & maxes & averages //{{{
          if (dist < GetArr3D(bond_n_mma, mol_i->Type, j / n_per_set, 0)) {
            SetArr3D(bond_n_mma, mol_i->Type, j / n_per_set, 0, dist);
          }
          if (dist > GetArr3D(bond_n_mma, mol_i->Type,
                              j / n_per_set, 1)) {
            SetArr3D(bond_n_mma, mol_i->Type, j / n_per_set, 1, dist);
          }
          AddArr3D(bond_n_mma, mol_i->Type, j / n_per_set, 2, dist);
          //}}}
          int k = dist / width; // distribution 'bin'
          if (k < bins) {
            AddArr3D(bond_n, mol_i->Type, j / n_per_set, k, 1);
          }
        }
      }
    }
  } //}}}
  // write to per-timestep file? //{{{
  if (opt.t_file[0] != '\0') {
    FILE *fw = OpenFile(opt.t_file, "a");
    fprintf(fw, " %7d", step.coor);
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System->MoleculeType[i];
      if (opt.mt[i] && mt->nBonds > 0) {
        for (int j = 0; j < mt->nBTypes; j++) {
          for (int k = j; k < mt->nBTypes; k++) {
            int btype[2] = {mt->BType[j], mt->BType[k]};
            if (BondedPair(*mt, btype)) {
              double value = GetArr3D(per_step, i, btype[0], btype[1]);
              if (GetArr3D(per_step_c, i, btype[0], btype[1]) > 0) {
                value /= GetArr3D(per_step_c, i, btype[0], btype[1]);
              }
              fprintf(fw, "%10f", value);
            }
          }
        }
      }
    }
    if (opt.all) {
      for (int i = 0; i < Count->MoleculeType; i++) {
        MOLECULETYPE *mt = &System->MoleculeType[i];
        if (opt.mt[i] && mt->nBonds > 0) {
          for (int j = 0; j < mt->nBonds; j++) {
            double value = GetArr2D(step_bond_all, i, j);
            if (GetArr2D(step_bond_all_c, i, j) > 0) {
              value /= GetArr2D(step_bond_all_c, i, j);
            }
            fprintf(fw, "%10f", value);
          }
        }
      }
    }
    if (opt.n_file[0] != '\0') {
      for (int i = 0; i < Count->MoleculeType; i++) {
        if (opt.mt[i]) {
          for (int j = 0; j < opt.n_number; j += n_per_set) {
            // skip id pairs if both are too high for the molecule //{{{
            if (opt.n_list[j] >= System->MoleculeType[i].nBeads &&
                opt.n_list[j+1] >= System->MoleculeType[i].nBeads) {
              continue;
            } //}}}
            double value = GetArr2D(step_t, i, j / n_per_set) /
                           GetArr2D(step_t_c, i, j / n_per_set);
            fprintf(fw, " %10f", value);
          }
        }
      }
    }
    putc('\n', fw);
    fclose(fw);
  } //}}}
  // free per-timestep array
  FreeArrND(per_step);
  FreeArrND(per_step_c);
  if (opt.t_file[0] != '\0' || opt.n_file[0] != '\0') {
    FreeArrND(step_t);
    FreeArrND(step_t_c);
  }
  if (opt.t_file[0] != '\0' && opt.all) {
    FreeArrND(step_bond_all);
    FreeArrND(step_bond_all_c);
  }
}
// structure for the callback function
struct user_data {
  OPT opt;
  ArrNDd *bond_bt, *bond_bt_mma,
         *bond_all, *bond_all_mma,
         *bond_n, *bond_n_mma;
  double width;
  int bins, n_per_set;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, *step, p->opt, p->bond_bt, p->bond_bt_mma,
              p->bond_all, p->bond_all_mma,
              p->bond_n, p->bond_n_mma, p->width, p->bins, p->n_per_set);
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}
  // <width> - width of a single bin //{{{
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <output> - file name with bond length distribution
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option - if present, don't join molecules //{{{
  opt.join = !BoolOption(argc, argv, "--joined");
  // --all option - use all bonds in the molecules
  opt.all = BoolOption(argc, argv, "--all");
  // -n option - specify bead ids to calculate distance between //{{{
  FileNumbersOption(argc, argv, 0, 100, "-n", opt.n_list,
                    &opt.n_number, opt.n_file, 'i');
  // if '-n' is present without numbers, use first and last for each molecule
  int n_per_set = 2; // it's a bond, so there two beads in each
  if (opt.n_file[0] != '\0' && opt.n_number == 0) {
    opt.n_number = n_per_set;
    opt.n_list[0] = 1;
    opt.n_list[1] = HIGHNUM; // large number to specify last bead
  }
  int n_pair_num = opt.n_number / n_per_set;
  // Error: wrong number of integers //{{{
  if (opt.n_file[0] != '\0' && (opt.n_number % n_per_set) != 0) {
    err_msg("number of bead indices must be even");
    PrintErrorOption("-n");
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < opt.n_number; i += n_per_set) {
    if (opt.n_list[i] == opt.n_list[i+1] ||
        opt.n_list[i] == 0 || opt.n_list[i+1] == 0) {
      err_msg("each pair of bead ids must be non-zero and different");
      PrintErrorOption("-n");
      exit(1);
    }
  } //}}}
  //}}}
  // -w option - bond length warning
  if (!OneNumberOption(argc, argv, "-w", &opt.warn, 'd')) {
    opt.warn = HIGHNUM;
  }
  // -t option - save per-timestep data
  if (!FileOption(argc, argv, "-t", opt.t_file)) {
    opt.t_file[0] = '\0';
  }
  //}}}
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  vec3d box = System.Box.Length;

  // '-mt <name(s)>' option
  if (!(opt.mt = calloc(System.Count.MoleculeType, sizeof *opt.mt))) {
    ErrorAlloc("opt.mt");
  }
  if (!TypeOption(argc, argv, "-mt", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

  /*
   * number of bins: *10 because of the -n option; bondlength should be at most
   * half boxsize, but distance between any two beads in a molecule can be large
  */
  int bins = Max3(box.x, box.y, box.z) / width * 10;

  // arrays for BeadType-BeadType bonds //{{{
  size_t shape4D[4] = {Count->MoleculeType,
                       Count->BeadType,
                       Count->BeadType,
                       bins};
  ArrNDd *bond_bt = CreateArrNDd(4, shape4D);
  if (!bond_bt) {
    ErrorAlloc("bond_bt");
  }
  shape4D[3] = 3;
  ArrNDd *bond_bt_mma = CreateArrNDd(4, shape4D);
  if (!bond_bt_mma) {
    ErrorAlloc("bond_bt_mma");
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        // set high number as initial minium bond length
        size_t id[4] = {i, j, k, 0};
        SetArrND(bond_bt_mma, id, HIGHNUM);
      }
    }
  } //}}}
  // arrays for all bonds in molecules //{{{
  ArrNDd *bond_all = NULL;
  ArrNDd *bond_all_mma = NULL;
  // maximum number of bonds in all molecules
  int max_bonds = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (System.MoleculeType[i].nBonds > max_bonds) {
      max_bonds = System.MoleculeType[i].nBonds;
    }
  }
  if (opt.all) {
    if (!(bond_all = CreateArr3Dd(Count->MoleculeType, max_bonds, bins)) ||
        !(bond_all_mma = CreateArr3Dd(Count->MoleculeType, max_bonds, 3))) {
      ErrorAlloc("bond_all/bond_all_mma");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        // set high number as initial minium bond length
        SetArr3D(bond_all_mma, i, j, 0, HIGHNUM);
      }
    }
  } //}}}
  // extra arrays for -n option //{{{
  ArrNDd *bond_n = NULL;
  ArrNDd *bond_n_mma = NULL;
  if (opt.n_file[0] != '\0') {
    if (!(bond_n = CreateArr3Dd(Count->MoleculeType, n_pair_num, bins)) ||
        !(bond_n_mma = CreateArr3Dd(Count->MoleculeType, n_pair_num, 3))) {
      ErrorAlloc("bond_n/bond_n_mma");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < n_pair_num; j++) {
        // set high number as initial minimum
        SetArr3D(bond_n_mma, i, j, 0, HIGHNUM);
      }
    }
  } //}}}

  // write initial stuff for per-timestep averages? //{{{
  if (opt.t_file[0] != '\0') {
    FILE *fw = PrintBylineOpenFile(opt.t_file, argc, argv);
    fprintf(fw, "# (1) distance\n");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      if (opt.mt[i] && mt->nBonds > 0) {
        fprintf(fw, "# %s molecule:", mt->Name);
        for (int j = 0; j < mt->nBTypes; j++) {
          for (int k = j; k < mt->nBTypes; k++) {
            int btype[2] = {mt->BType[j], mt->BType[k]};
            if (BondedPair(*mt, btype)) {
              BEADTYPE *bt_lo = &System.BeadType[btype[0]];
              BEADTYPE *bt_hi = &System.BeadType[btype[1]];
              fprintf(fw, " (%d) %s-%s", ++count, bt_lo->Name, bt_hi->Name);
            }
          }
        }
        if (opt.all) {
          count++;
          fprintf(fw, " (%d)-(%d) individual bonds", count,
                  count + mt->nBonds - 1);
          count += mt->nBonds - 1;
        }
        putc('\n', fw);
      }
    }
    if (opt.n_file[0] != '\0') { //{{{
      // print the second line - molecule names and ids with column numbers
      fprintf(fw, "# from -n option:\n");
      for (int i = 0; i < Count->MoleculeType; i++) {
        MOLECULETYPE *mt_i = &System.MoleculeType[i];
        if (opt.mt[i]) {
          fprintf(fw, "# %s:", mt_i->Name);

          for (int j = 0; j < opt.n_number; j += n_per_set) {
            // skip id pairs if both are too high for the molecule //{{{
            if (opt.n_list[j] >= mt_i->nBeads &&
                opt.n_list[j+1] >= mt_i->nBeads) {
              continue;
            } //}}}
            // bead ids for the distance //{{{
            int id1, id2;
            // use first molecule bead if bead index too high or -1
            if (opt.n_list[j] >= mt_i->nBeads) {
              id1 = mt_i->nBeads;
            } else { // use specified index otherwise
              id1 = opt.n_list[j];
            }
            // use last molecule bead if bead index too high or -1
            if (opt.n_list[j+1] >= mt_i->nBeads) {
              id2 = mt_i->nBeads;
            } else { // use specified index otherwise
              id2 = opt.n_list[j+1];
            }
            // write the numbers so the first is higher
            if (id1 > id2) {
              SwapInt(&id1, &id2);
            } //}}}
            fprintf(fw, " (%d) %d-%d", ++count, id1, id2);
          }
        }
        putc('\n', fw);
      }
    } //}}}
    fclose(fw);
  } //}}}

  STEP step = InitStep;
  struct user_data ud = { opt, bond_bt, bond_bt_mma,
                          bond_all, bond_all_mma,
                          bond_n, bond_n_mma,
                          width, bins, n_per_set };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // sum up all bonds in molecules (normalization factor) //{{{
  // BeadType-BeadType bonds
  ArrNDi *bond_bt_norm = CreateArr3Di(Count->MoleculeType,
                                      Count->BeadType, Count->BeadType);
  if (!bond_bt_norm) {
    ErrorAlloc("bond_bt_norm");
  }
  int range[2] = {bins, 0}; // min/max distance
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = 0; l < bins; l++) {
          size_t id[4] = {i, j, k, l};
          if (GetArrND(bond_bt, id) > 0) {
            AddArr3D(bond_bt_norm, i, j, k, GetArrND(bond_bt, id));
            if (l < range[0]) {
              range[0] = l;
            }
            if (l > range[1]) {
              range[1] = l;
            }
          }
        }
      }
    }
  }
  // all molecules' bonds
  ArrNDi *bond_all_norm = NULL;
  if (opt.all) {
    if (!(bond_all_norm = CreateArr2Di(Count->MoleculeType, max_bonds))) {
      ErrorAlloc("bond_all_norm");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        for (int k = 0; k < bins; k++) {
          AddArr2D(bond_all_norm, i, j, GetArr3D(bond_all, i, j, k));
          if (GetArr3D(bond_all, i, j, k) && k < range[0]) {
            range[0] = k;
          }
          if (GetArr3D(bond_all, i, j, k) && k > range[1]) {
            range[1] = k;
          }
        }
      }
    }
  }
  // include nearest 0 values in the range of distances
  if (range[0] > 0) {
    range[0]--;
  }
  range[1] += 2; // +2 as for loop is range[0]...range[1]-1
  if (range[1] > bins) { // clamp range[1] if it's too high
    range[1] = bins;
  } //}}}

  if (opt.t_file[0] != '\0') {
    FILE *fw = OpenFile(opt.t_file, "a");
    WriteMinsMaxesAvgs(fw, System, opt, bond_bt_mma, bond_bt_norm,
                       bond_all_mma, bond_all_norm);
    fclose(fw);
  }

  // write distribution of bond lengths //{{{
  // print first lines of output file - molecule names and beadtype pairs //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# (1) distance\n");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (opt.mt[i] && mt->nBonds > 0) {
      fprintf(fw, "# %s molecule:", mt->Name);
      for (int j = 0; j < mt->nBTypes; j++) {
        for (int k = j; k < mt->nBTypes; k++) {
          int btype1 = mt->BType[j],
              btype2 = mt->BType[k];
          SortPairAsc(&btype1, &btype2);
          if (GetArr3D(bond_bt_norm, i, btype1, btype2) > 0) {
            count++;
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                                              System.BeadType[btype2].Name);
          }
        }
      }
      if (opt.all) {
        count++;
        fprintf(fw, " (%d)-(%d) individual bonds", count,
                count + mt->nBonds - 1);
        count += mt->nBonds - 1;
      }
      putc('\n', fw);
    }
  } //}}}
  int ncols = count;
  int nrows = range[1] - range[0];
  ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
  if (!data) {
    ErrorAlloc("data");
  }
  // collate data //{{{
  for (int i = range[0]; i < range[1]; i++) {
    int row = i - range[0];
    count = 0;
    SetArr2D(data, row, count++, width * (2 * i + 1) / 2);
    for (int j = 0; j < Count->MoleculeType; j++) {
      MOLECULETYPE *mt = &System.MoleculeType[j];
      if (opt.mt[j] && mt->nBonds > 0) {
        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < mt->nBTypes; k++) {
          for (int l = k; l < mt->nBTypes; l++) {
            int btype1 = mt->BType[k],
                btype2 = mt->BType[l];
            SortPairAsc(&btype1, &btype2);
            if (GetArr3D(bond_bt_norm, j, btype1, btype2) > 0) {
              size_t id[4] = {j, btype1, btype2, i};
              double value = GetArrND(bond_bt, id) /
                             GetArr3D(bond_bt_norm, j, btype1, btype2);
              SetArr2D(data, row, count++, value);
            }
          }
        }
        if (opt.all) {
          for (int k = 0; k < mt->nBonds; k++) {
            double value = GetArr3D(bond_all, j, k, i);
            if (GetArr2D(bond_all_norm, j, k) > 0) {
              value /= GetArr2D(bond_all_norm, j, k);
            }
            SetArr2D(data, row, count++, value);
          }
        }
      }
    }
  } //}}}
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(fw, nrows, ncols, data);
  FreeArrND(data);
  WriteMinsMaxesAvgs(fw, System, opt, bond_bt_mma, bond_bt_norm,
                         bond_all_mma, bond_all_norm);
  fclose(fw); //}}}

  // write distribution of distances from '-n' option //{{{
  if (opt.n_file[0] != '\0') {
    // sum up all calculated distances (normalization factors) //{{{
    ArrNDi *n_norm = CreateArr2Di(Count->MoleculeType, n_pair_num);
    if (!n_norm) {
      ErrorAlloc("n_norm");
    }
    range[0] = bins;
    range[1] = 0;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt.mt[i]) {
        for (int j = 0; j < n_pair_num; j++) {
          for (int k = 0; k < bins; k++) {
            if (GetArr3D(bond_n, i, j, k) > 0) {
              AddArr2D(n_norm, i, j, GetArr3D(bond_n, i, j, k));
              if (k < range[0]) {
                range[0] = k;
              }
              if (k > range[1]) {
                range[1] = k;
              }
            }
          }
        }
      }
    }
    // include nearest 0 values in the range of distances
    if (range[0] > 0) {
      range[0]--;
    }
    range[1] += 2; // +2 as for loop is range[0]...range[1]-1
    if (range[1] > bins) { // clamp range[1] if it's too high
      range[1] = bins;
    } //}}}
    FILE *fw = PrintBylineOpenFile(opt.n_file, argc, argv);
    // print the first line - molecule names with bead order //{{{
    fprintf(fw, "# bead order in molecule(s) -");
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt.mt[i]) {
        fprintf(fw, " %s:", mt_i->Name);
        for (int j = 0; j < mt_i->nBeads; j++) {
          int btype = mt_i->Bead[j];
          fprintf(fw, " %s", System.BeadType[btype].Name);
        }
        putc(';', fw);
      }
    }
    putc('\n', fw); //}}}
    // print the second line - molecule names and ids with column numbers //{{{
    fprintf(fw, "# (1) distance\n");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt.mt[i]) {
        fprintf(fw, "# %s:", mt_i->Name);

        for (int j = 0; j < opt.n_number; j += n_per_set) {
          // skip id pairs if both are too high for the molecule //{{{
          // TODO: condition used more times - change into inline function
          if (opt.n_list[j] >= mt_i->nBeads &&
              opt.n_list[j+1] >= mt_i->nBeads) {
            continue;
          } //}}}
          // bead ids for the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt.n_list[j] >= mt_i->nBeads) {
            id1 = mt_i->nBeads;
          } else { // use specified index otherwise
            id1 = opt.n_list[j];
          }
          // use last molecule bead if bead index too high or -1
          if (opt.n_list[j+1] >= mt_i->nBeads) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt.n_list[j+1];
          }
          // write the numbers so the first is higher
          if (id1 > id2) {
            SwapInt(&id1, &id2);
          } //}}}
          fprintf(fw, " (%d) %d-%d", ++count, id1, id2);
        }
      }
      putc('\n', fw);
    } //}}}
    // collate data //{{{
    int ncols = count;
    int nrows = range[1] - range[0];
    ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
    if (!data) {
      ErrorAlloc("data");
    }
    for (int i = range[0]; i < range[1]; i++) {
      int row = i - range[0];
      count = 0;
      SetArr2D(data, row, count++, width * (2 * i + 1) / 2);
      for (int j = 0; j < Count->MoleculeType; j++) {
        if (opt.mt[j]) {
          for (int k = 0; k < opt.n_number; k += n_per_set) {
            // skip id pairs if both are too high for the molecule //{{{
            if (opt.n_list[k] >= System.MoleculeType[j].nBeads &&
                opt.n_list[k+1] >= System.MoleculeType[j].nBeads) {
              continue;
            } //}}}
            double val = (double)(GetArr3D(bond_n, j, k/n_per_set, i)) /
                         GetArr2D(n_norm, j, k/n_per_set);
            SetArr2D(data, row, count++, val);
          }
        }
      }
    } //}}}
    ComputeColumnWidths(nrows, ncols, data, 6);
    PrintDataAll(fw, nrows, ncols, data);
    FreeArrND(data);
    // write mins and maxes
    // legend line //{{{
    fprintf(fw, "# min(1st columns)/max(2nd columns)/average(3rd columns)\n");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt.mt[i]) {
        fprintf(fw, "# %s:", mt_i->Name);
        for (int j = 0; j < opt.n_number; j += n_per_set) {
          // skip id pairs if both are too high for the molecule //{{{
          if (opt.n_list[j] >= mt_i->nBeads &&
              opt.n_list[j+1] >= mt_i->nBeads) {
            continue;
          } //}}}
          // bead ids the distance //{{{
          int id1, id2;
          // use last molecule bead if bead index is too high
          if (opt.n_list[j] >= mt_i->nBeads) {
            id1 = mt_i->nBeads;
          } else { // use specified index otherwise
            id1 = opt.n_list[j];
          }
          // use last molecule bead if bead index is too high
          if (opt.n_list[j+1] >= mt_i->nBeads) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt.n_list[j+1];
          }
          // write the numbers so the first is lower
          if (id1 > id2) {
            SwapInt(&id1, &id2);
          } //}}}
          fprintf(fw, " (%d) %d-%d", count, id1, id2);
          count += 3;
        }
      }
      putc('\n', fw);
    } //}}}
    // data line //{{{
    putc('#', fw);
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt.mt[i]) {
        for (int j = 0; j < opt.n_number; j += n_per_set) {
          // skip id pairs if both are too high for the molecule //{{{
          if (opt.n_list[j] >= System.MoleculeType[i].nBeads &&
              opt.n_list[j+1] >= System.MoleculeType[i].nBeads) {
            continue;
          } //}}}
          // if this bin is filled, its max must be larger than 0
          int bond_id = j / n_per_set;
          if (GetArr2D(n_norm, i, bond_id) > 0) {
            fprintf(fw, " %lf", GetArr3D(bond_n_mma, i, bond_id, 0));
            fprintf(fw, " %lf", GetArr3D(bond_n_mma, i, bond_id, 1));
            double val = GetArr3D(bond_n_mma, i, bond_id, 2) /
                         GetArr2D(n_norm, i, bond_id);
            fprintf(fw, " %lf", val);
          }
        }
      }
    }
    putc('\n', fw); //}}}
    fclose(fw);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(opt.mt);
  FreeArrND(bond_bt);
  FreeArrND(bond_bt_mma);
  FreeArrND(bond_bt_norm);
  if (opt.all) {
    FreeArrND(bond_all);
    FreeArrND(bond_all_mma);
    FreeArrND(bond_all_norm);
  }
  // free arrays for -n option
  if (opt.n_file[0] != '\0') {
    FreeArrND(bond_n);
    FreeArrND(bond_n_mma);
  }
  FreeSystem(&System);
  //}}}

  return 0;
}

// helper functions
// is given bead type pair in a bond in the given molecule? //{{{
static bool BondedPair(MOLECULETYPE mt, int btype[2]) {
  SortPairAsc(&btype[0], &btype[1]);
  for (int l = 0; l < mt.nBonds; l++) {
    int bead[2] = {mt.Bead[mt.Bond[l][0]],
                   mt.Bead[mt.Bond[l][1]]};
    if (bead[0] > bead[1]) {
      SwapInt(&bead[0], &bead[1]);
    }
    if (bead[0] == btype[0] && bead[1] == btype[1]) {
      return true;
    }
  }
  return false;
} //}}}
// write mins, maxes, and averages //{{{
static void WriteMinsMaxesAvgs(FILE *fw, SYSTEM System, OPT opt,
                               ArrNDd *bond_bt_mma, ArrNDi *bond_bt_norm,
                               ArrNDd *bond_all_mma, ArrNDi *bond_all_norm) {
  fprintf(fw, "# min(1st columns)/max(2nd columns)/average(3rd columns)\n");
  int count = 1;
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (opt.mt[i] && mt->nBonds > 0) {
      fprintf(fw, "# %s molecule:", mt->Name);
      for (int j = 0; j < mt->nBTypes; j++) {
        for (int k = j; k < mt->nBTypes; k++) {
          int btype1 = mt->BType[j],
              btype2 = mt->BType[k];
          SortPairAsc(&btype1, &btype2);
          if (GetArr3D(bond_bt_norm, i, btype1, btype2) > 0) {
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                                              System.BeadType[btype2].Name);
            count += 3;
          }
        }
      }
      if (opt.all) {
        fprintf(fw, " (%d)-(%d) individual bonds", count,
                count + 3 * mt->nBonds - 1);
        count += 3 * mt->nBonds;
      }
      putc('\n', fw);
    }
  }
  // data line
  putc('#', fw);
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    for (int j = 0; j < System.Count.BeadType; j++) {
      for (int k = j; k < System.Count.BeadType; k++) {
        // if this bin is filled, its max must be larger than 0
        size_t id_0[4] = {i, j, k, 0};
        size_t id_1[4] = {i, j, k, 1};
        size_t id_2[4] = {i, j, k, 2};
        if (GetArrND(bond_bt_mma, id_1) > 0) {
          fprintf(fw, " %lf", GetArrND(bond_bt_mma, id_0));
          fprintf(fw, " %lf", GetArrND(bond_bt_mma, id_1));
          double val = GetArrND(bond_bt_mma, id_2) /
                       GetArr3D(bond_bt_norm, i, j, k);
          fprintf(fw, " %lf", val);
        }
      }
    }
    if (opt.all) {
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        fprintf(fw, " %lf", GetArr3D(bond_all_mma, i, j, 0));
        fprintf(fw, " %lf", GetArr3D(bond_all_mma, i, j, 1));
        double val = GetArr3D(bond_all_mma, i, j, 2) /
                     GetArr2D(bond_all_norm, i, j);
        fprintf(fw, " %lf", val);
      }
    }
  }
  putc('\n', fw);
} //}}}
