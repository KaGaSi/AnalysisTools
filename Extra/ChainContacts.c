#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "TBW",

  "Usage: ChainContacts <input> <output> <bead(s)> [options]",
  .args = 2, // number of mandatory arguments
  .all = 15, // number of valid lines OptSpec (not counting last {NULL})
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
  {"<output>", NULL, "output file with pair correlation function(s)", OPT_ARG},
  {"<bead(s)>", NULL, "bead name(s) for calculation (optional and ignored if '--all' is used)", OPT_ARG},
  {"-d", "<dist>", "maximum distance for RDF calculation (default: 1/3 of the shortest box side length)", OPT_EXTRA},
  {"-skip", "<num>", "number of in-between beads to skip in the same molecule", OPT_EXTRA},
  {"--all", NULL, "use all bead types (overwrites <bead(s)>)", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  bool all;        // --all
  bool *bt;
  double max_dist; // -d
  int skip;        // -skip
}; //}}}

// get bead position in molecule //{{{
int PosInMol(const int b_id, const SYSTEM System) {
  int m_id = System.Bead[b_id].Molecule;
  int mtype = System.Molecule[m_id].Type;
  for (int i = 0; i < System.MoleculeType[mtype].nBeads; i++) {
    if (b_id == System.Molecule[m_id].Bead[i]) {
      return i;
    }
  }
  err_msg("bead not in the molecule!");
  PrintError();
  exit(1);
} //}}}

// TODO: move to some library file
static inline void CorrectBTypeOrder(int *btype_i, int *btype_j) {
  if (*btype_i > *btype_j) {
    SwapInt(btype_i, btype_j);
  }
}

// functions to plug into traversal functions
// calculate contacts between bead types //{{{
// the calculation itself
static void CalculateContacts(int id_i, int id_j, SYSTEM System, OPT opt,
                              ArrNDi *contacts_step, double max_dist) {
  int i = System.BeadCoor[id_i];
  int j = System.BeadCoor[id_j];
  BEAD *b_i = &System.Bead[i];
  BEAD *b_j = &System.Bead[j];
  // calculate distance between the two beads
  vec3d d = DistancePBC(b_i->Position, b_j->Position, &System.Box);
  double dist = VectLength(d);
  if (b_i->Molecule != b_j->Molecule ||
      abs(PosInMol(i, System) - PosInMol(j, System)) > opt.skip) {
    if (dist < max_dist) {
      int btype_i = b_i->Type;
      int btype_j = b_j->Type;
      CorrectBTypeOrder(&btype_i, &btype_j);
      AddArr2D(contacts_step, btype_i, btype_j, 1);
    }
  }
}
// structure for the callback function
struct contacts_args {
  ArrNDi *contacts_step;
  double max_dist;
  OPT opt;
};
// adaptor for the CalculatePCF() function
static void CalculateContacts_adaptor(int id_i, int id_j,
                                      const SYSTEM System, void *ud) {
  struct contacts_args *p = (struct contacts_args*)ud;
  CalculateContacts(id_i, id_j, System, p->opt, p->contacts_step, p->max_dist);
} //}}}
// condition for using specified beads //{{{
// check based on supplied type (needed for writing to file)
static bool CheckBeadType(int btype, OPT opt) {
  if (!opt.bt[btype]) {
    return false;
  } else {
    return true;
  }
}
// check based on supplied System.BeedCoor id
static bool CheckBead(int id_i, SYSTEM System, OPT opt) {
  int i = System.BeadCoor[id_i];
  return CheckBeadType(System.Bead[i].Type, opt);
}
// structure for the callback function
struct check_args {
  OPT opt;
};
// adaptor for the CalculatePCF() function
static bool CheckBeadType_adaptor(int id_i, SYSTEM System, void *ud) {
  struct check_args *p = (struct check_args*)ud;
  return CheckBead(id_i, System, p->opt);
}
//}}}

void Calculation(SYSTEM *System, STEP *step, OPT opt, ArrNDd *contacts,
                 double cell_size, char output[LINE]) {
  COUNT *Count = &System->Count;
  ArrNDi *contacts_step = CreateArr2Di(Count->BeadType, Count->BeadType);
  FillArrND(contacts_step, 0);
  WrapJoinCoordinates(System, true, false);
  struct contacts_args args = { contacts_step, opt.max_dist, opt };
  struct check_args check = { opt };
  TraversePairs(*System, cell_size, CalculateContacts_adaptor, &args,
                CheckBeadType_adaptor, &check);
  // calculate molecules for normalisation - all containing used beads
  int count_mols = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < System->MoleculeType[i].nBTypes; j++) {
      if (CheckBeadType(System->MoleculeType[i].BType[j], opt)) {
        count_mols += System->MoleculeType[i].Number;
        break;
      }
    }
  }
  // save per-step values
  FILE *fw = OpenFile(output, "a");
  fprintf(fw, "%6d", step->coor);
  for (int i = 0; i < Count->BeadType; i++) {
    if (CheckBeadType(i, opt)) {
      for (int j = i; j < Count->BeadType; j++) {
        if (CheckBeadType(j, opt)) {
          double val = (double)(GetArr2D(contacts_step, i, j)) / count_mols;
          fprintf(fw, " %lf", val);
          AddArr2D(contacts, i, j, val);
        }
      }
    }
  }
  fprintf(fw, "\n");
  fclose(fw);
  FreeArrND(contacts_step);
}
// structure for the callback function
struct user_data {
  OPT opt;
  ArrNDd *contacts;
  double cell_size;
  char *output;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, step, p->opt, p->contacts, p->cell_size, p->output);
};

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, false, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}
  // <output> - filename with pcf(s)
  char output[LINE] = "";
  s_strcpy(output, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }
  // use all bead types in the structure file?
  opt.all = BoolOption(argc, argv, "--all");
  //}}}

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  vec3d box = System.Box.Length;

  // <bead(s)> - names of bead types to use //{{{
  if (!(opt.bt = calloc(Count->BeadType, sizeof *opt.bt))) {
    ErrorAlloc("opt.bt");
  }
  if (opt.all) {
    for (int i = 0; i < Count->BeadType; i++) {
      opt.bt[i] = true;
    }
  } else {
    for (int i = 0; i < Count->BeadType; i++) {
      opt.bt[i] = false;
    }
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], System);
      if (type == -1) {
        ErrorBeadType(argv[count], System);
        exit(1);
      }
      if (opt.bt[type]) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
                 ErrYellow(), argv[count], ErrCyan());
        PrintWarning();
      }
      opt.bt[type] = true;
    }
    count--; // while always increments count at least once
    if (count < (HelpDesc.args + 1)) {
      err_msg("missing <bead(s)> or --all option");
      PrintError();
      PrintCommand(stderr, argc, argv);
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  opt.max_dist = Min3(box.x, box.y, box.z) / 3;
  if (OneNumberOption(argc, argv, "-d", &opt.max_dist, 'd') &&
      opt.max_dist <= 0) {
    err_msg("distance must be a positive number");
    ErrorOption("-d");
  }
  double cell_size = opt.max_dist;
  opt.skip = 0;
  if (OneNumberOption(argc, argv, "-skip", &opt.skip, 'i') &&
      opt.skip < 0) {
    err_msg("non-negative number of beads must be skipped");
    ErrorOption("-skip");
  }

  // pair correlation function
  ArrNDd *contacts = CreateArr2Dd(Count->BeadType, Count->BeadType);
  if (!contacts) {
    ErrorAlloc("contacts");
  }

  FILE *fw = PrintBylineOpenFile(output, argc, argv);
  count = 1;
  fprintf(fw, "# (%d) step", count);
  count++;
  for (int i = 0; i < Count->BeadType; i++) {
    if (opt.bt[i]) {
      for (int j = i; j < Count->BeadType; j++) {
        if (opt.bt[j]) {
          fprintf(fw, "; (%d) %s-%s", count, System.BeadType[i].Name,
                                             System.BeadType[j].Name);
          count++;
        }
      }
    }
  }
  putc('\n', fw);
  fclose(fw);

  STEP step = InitStep;
  struct user_data ud = { opt, contacts, cell_size, output };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // append overall averages //{{{
  FILE *fw_avg = OpenFile(output, "a");
  fprintf(fw_avg, "# overall averages (%d steps):\n", step.used);
  fprintf(fw_avg, "#");
  for (int i = 0; i < Count->BeadType; i++) {
    if (CheckBeadType(i, opt)) {
      for (int j = i; j < Count->BeadType; j++) {
        if (CheckBeadType(j, opt)) {
          fprintf(fw_avg, " <%s-%s>", System.BeadType[i].Name,
                                      System.BeadType[j].Name);
        }
      }
    }
  }
  fprintf(fw_avg, "\n#");
  for (int i = 0; i < Count->BeadType; i++) {
    if (CheckBeadType(i, opt)) {
      for (int j = i; j < Count->BeadType; j++) {
        if (CheckBeadType(j, opt)) {
          fprintf(fw_avg, " %lf", GetArr2D(contacts, i, j) / step.used);
        }
      }
    }
  }
  putc('\n', fw_avg);
  fclose(fw_avg); //}}}

  // free memory //{{{
  FreeArrND(contacts);
  free(opt.bt);
  FreeSystem(&System);
  //}}}

  return 0;
}

// backup - will the MainLoop() function work?
