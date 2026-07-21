#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "TBW",

  "Usage: ChainContacts <input> <output> <bead(s)> [options]",
  .args = 2, // number of mandatory arguments
  .all = 15, // number of valid lines OptSpec (not counting last {nullptr})
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
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<output>", nullptr, "output file with number of isolated beads per chain",
    OPT_ARG},
  {"<bead(s)>", nullptr, "bead name(s) for calculation "
    "(optional and ignored if '--all' is used)", OPT_ARG},
  {"-d", "<dist>", "maximum distance for contact calculation "
    "(default: 1/3 of the shortest box side length)", OPT_EXTRA},
  {"-skip", "<num>", "number of in-between beads to skip in the same molecule",
    OPT_EXTRA},
  {"--all", nullptr, "use all bead types (overwrites <bead(s)>)", OPT_EXTRA},
  {nullptr}
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

// functions to plug into traversal functions
// mark beads that have at least one contact //{{{
static void MarkContact(int id_i, int id_j, SYSTEM System, OPT opt,
                        bool *in_contact, double max_dist) {
  int i = System.BeadCoor[id_i];
  int j = System.BeadCoor[id_j];
  BEAD *b_i = &System.Bead[i];
  BEAD *b_j = &System.Bead[j];
  vec3d d = DistancePBC(b_i->Position, b_j->Position, &System.Box);
  double dist = VectLength(d);
  if (b_i->Molecule != b_j->Molecule ||
      abs(PosInMol(i, System) - PosInMol(j, System)) > opt.skip) {
    if (dist < max_dist) {
      in_contact[i] = true;
      in_contact[j] = true;
    }
  }
}
// structure for the callback function
struct contacts_args {
  bool *in_contact;
  double max_dist;
  OPT opt;
};
// adaptor for the TraversePairs() function
static void MarkContact_adaptor(int id_i, int id_j,
                                const SYSTEM System, void *ud) {
  struct contacts_args *p = (struct contacts_args*)ud;
  MarkContact(id_i, id_j, System, p->opt, p->in_contact, p->max_dist);
} //}}}
// condition for using specified beads //{{{
static bool CheckBeadType(int btype, OPT opt) {
  return opt.bt[btype];
}
static bool CheckBead(int id_i, SYSTEM System, OPT opt) {
  int i = System.BeadCoor[id_i];
  return CheckBeadType(System.Bead[i].Type, opt);
}
struct check_args {
  OPT opt;
};
static bool CheckBeadType_adaptor(int id_i, SYSTEM System, void *ud) {
  struct check_args *p = (struct check_args*)ud;
  return CheckBead(id_i, System, p->opt);
}
//}}}

void Calculation(SYSTEM *System, STEP *step, OPT opt, double *isolated_avg,
                 double cell_size, char output[LINE]) {
  COUNT *Count = &System->Count;
  bool *in_contact = calloc(Count->Bead, sizeof *in_contact);
  if (!in_contact) {
    ErrorAlloc("in_contact");
  }
  WrapJoinCoordinates(System, true, false);
  struct contacts_args args = { in_contact, opt.max_dist, opt };
  struct check_args check = { opt };
  TraversePairs(*System, cell_size, MarkContact_adaptor, &args,
                CheckBeadType_adaptor, &check);
  // count isolated selected beads per molecule type
  int *isolated_step = calloc(Count->MoleculeType, sizeof *isolated_step);
  if (!isolated_step) {
    ErrorAlloc("isolated_step");
  }
  for (int m = 0; m < Count->Molecule; m++) {
    int mtype = System->Molecule[m].Type;
    for (int k = 0; k < System->MoleculeType[mtype].nBeads; k++) {
      int bid = System->Molecule[m].Bead[k];
      if (opt.bt[System->Bead[bid].Type] && !in_contact[bid]) {
        isolated_step[mtype]++;
      }
    }
  }
  // save per-step values
  FILE *fw = OpenFile(output, "a");
  fprintf(fw, "%6d", step->coor);
  for (int i = 0; i < Count->MoleculeType; i++) {
    bool has_selected = false;
    for (int j = 0; j < System->MoleculeType[i].nBTypes; j++) {
      if (opt.bt[System->MoleculeType[i].BType[j]]) {
        has_selected = true;
        break;
      }
    }
    if (has_selected) {
      double val = (double)isolated_step[i] / System->MoleculeType[i].Number;
      fprintf(fw, " %lf", val);
      isolated_avg[i] += val;
    }
  }
  fprintf(fw, "\n");
  fclose(fw);
  free(in_contact);
  free(isolated_step);
}
// structure for the callback function
struct user_data {
  OPT opt;
  double *isolated_avg;
  double cell_size;
  char *output;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, step, p->opt, p->isolated_avg, p->cell_size, p->output);
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
  // <output> - filename with isolated bead counts
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

  // per-molecule-type accumulator for isolated bead counts
  double *isolated_avg = calloc(Count->MoleculeType, sizeof *isolated_avg);
  if (!isolated_avg) {
    ErrorAlloc("isolated_avg");
  }

  FILE *fw = PrintBylineOpenFile(output, argc, argv);
  count = 1;
  fprintf(fw, "# (%d) step", count);
  count++;
  for (int i = 0; i < Count->MoleculeType; i++) {
    bool has_selected = false;
    for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
      if (opt.bt[System.MoleculeType[i].BType[j]]) {
        has_selected = true;
        break;
      }
    }
    if (has_selected) {
      fprintf(fw, "; (%d) %s", count, System.MoleculeType[i].Name);
      count++;
    }
  }
  putc('\n', fw);
  fclose(fw);

  STEP step = InitStep;
  struct user_data ud = { opt, isolated_avg, cell_size, output };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // append overall averages //{{{
  FILE *fw_avg = OpenFile(output, "a");
  fprintf(fw_avg, "# overall averages (%d steps):\n", step.used);
  fprintf(fw_avg, "#");
  for (int i = 0; i < Count->MoleculeType; i++) {
    bool has_selected = false;
    for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
      if (opt.bt[System.MoleculeType[i].BType[j]]) {
        has_selected = true;
        break;
      }
    }
    if (has_selected) {
      fprintf(fw_avg, " <%s>", System.MoleculeType[i].Name);
    }
  }
  fprintf(fw_avg, "\n#");
  for (int i = 0; i < Count->MoleculeType; i++) {
    bool has_selected = false;
    for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
      if (opt.bt[System.MoleculeType[i].BType[j]]) {
        has_selected = true;
        break;
      }
    }
    if (has_selected) {
      fprintf(fw_avg, " %lf", isolated_avg[i] / step.used);
    }
  }
  putc('\n', fw_avg);
  fclose(fw_avg); //}}}

  // free memory //{{{
  free(isolated_avg);
  free(opt.bt);
  FreeSystem(&System);
  //}}}

  return 0;
}
