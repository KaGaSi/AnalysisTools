#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "GyrationMolecules calculates the gyration tensor for molecules and "
  "determines their shape descriptors like the radius of gyration, "
  "acylindricity, asphericity, or relative shape anisotropy. It writes "
  "per-timestep averages to output files (one per molecule type) and "
  "appends overall averages to those files.",

  "Usage: GyrationMolecules <input> <output> [options]",
  .args = 2,  // number of mandatory arguments
  .all = 14,  // number of valid lines in OptSpec (not counting last {nullptr})
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
  {"<input>",  nullptr, "input coordinate file", OPT_ARG},
  {"<output>", nullptr, "output base name (appends '-<molname>.txt')", OPT_ARG},
  {"--joined", nullptr, "<input> contains joined coordinates", OPT_EXTRA},
  {"-mt", "<name(s)>", "molecule types to use (default: all)", OPT_EXTRA},
  {"-bt", "<name(s)>", "bead types to use (default: all)", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  bool join,  // --joined
       *mt,   // -mt (per molecule type)
       *bt;   // -bt (per bead type)
}; //}}}

// column indices for per-type accumulator arrays; COL_RE has its own
// normalization count (COL_RE_N) as end beads may be absent from a frame
enum { COL_RG, COL_SQRRG, COL_ANIS, COL_ACYL, COL_ASPHER,
       COL_EIGEN0, COL_EIGEN1, COL_EIGEN2, COL_RE, COL_RE_N, N_COLS };

// all state shared between main() and the per-timestep callback //{{{
struct user_data {
  struct OPT opt;
  ArrNDd *sums;      // overall sums across all timesteps
  int *mol_count;    // total molecule count across all timesteps
  ArrNDd *step_vals; // per-timestep sums (zeroed each call)
  int *step_count;   // per-timestep count (zeroed each call)
  FILE **files;      // open output file handles
  int *list;         // reusable bead-id buffer
}; //}}}

// per-timestep calculation and output //{{{
static void Calculation(SYSTEM *System, STEP *step, struct user_data *ud) {
  const struct OPT opt = ud->opt;
  COUNT *Count = &System->Count;

  // zero per-step accumulators (reuse pre-allocated memory)
  FillArrND(ud->step_vals, 0.0);
  memset(ud->step_count, 0, Count->MoleculeType * sizeof *ud->step_count);

  // wrap and join only when coordinates are not already joined
  WrapJoinCoordinates(System, opt.join, false);
  if (opt.join) {
    for (int i = 0; i < Count->Molecule; i++) {
      if (opt.mt[System->Molecule[i].Type]) {
        RemovePBCMolecule(i, System);
      }
    }
  }

  // calculate shape descriptors for each molecule
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol = &System->Molecule[i];
    MOLECULETYPE *mtype = &System->MoleculeType[mol->Type];
    if (!opt.mt[mol->Type]) {
      continue;
    }

    // build list of bead ids present in this frame and matching -bt
    int n = 0;
    for (int j = 0; j < mtype->nBeads; j++) {
      int id = mol->Bead[j];
      if (opt.bt[System->Bead[id].Type] && System->Bead[id].InTimestep) {
        ud->list[n] = id;
        n++;
      }
    }
    if (n < 2) { // need at least 2 beads for gyration
      continue;
    }

    // end-to-end distance is computed before Gyration() translates bead
    // positions, so both end beads must to be present in the timestep
    BEAD *b_first = &System->Bead[mol->Bead[0]];
    BEAD *b_last = &System->Bead[mol->Bead[mtype->nBeads-1]];
    bool re_valid = false;
    if (b_first->InTimestep && b_last->InTimestep) {
      re_valid = true;
    }
    double Re = 0;
    if (re_valid) {
      Re = VectLength(Vector(b_first->Position, b_last->Position));
    }

    vec3d eigen = Gyration(n, ud->list, System);
    // skip degenerate (all-zero) or numerically corrupt (negative) eigenvalues;
    // Gyration() already warns about the latter
    if (eigen.x < 0 || (eigen.x == 0 && eigen.y == 0 && eigen.z == 0)) {
      continue;
    }

    int mt = mol->Type;
    // radius of gyration
    double Rgi = sqrt(eigen.x + eigen.y + eigen.z);
    AddArr2D(ud->step_vals, mt, COL_RG, Rgi);
    double val = Square(Rgi);
    AddArr2D(ud->step_vals, mt, COL_SQRRG, val);
    // relative shape anisotropy
    val = 1.5 * SqVectLength(eigen) / Square(eigen.x + eigen.y + eigen.z) - 0.5;
    AddArr2D(ud->step_vals, mt, COL_ANIS, val);
    // acylindricity
    val = eigen.y - eigen.x;
    AddArr2D(ud->step_vals, mt, COL_ACYL, val);
    // asphericity
    val = eigen.z - 0.5 * (eigen.x + eigen.y);
    AddArr2D(ud->step_vals, mt, COL_ASPHER, val);
    // eigenvalues
    for (int dd = 0; dd < 3; dd++) {
      AddArr2D(ud->step_vals, mt, COL_EIGEN0 + dd, eigen.v[dd]);
    }
    // end-to-end distance
    if (re_valid) {
      AddArr2D(ud->step_vals, mt, COL_RE, Re);
      AddArr2D(ud->step_vals, mt, COL_RE_N, 1);
    }
    ud->step_count[mt]++;
  }

  // add to overall sums and write per-timestep data
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i] || ud->step_count[i] == 0) {
      continue;
    }
    for (int c = 0; c < N_COLS; c++) {
      AddArr2D(ud->sums, i, c, GetArr2D(ud->step_vals, i, c));
    }
    ud->mol_count[i] += ud->step_count[i];

    int n = ud->step_count[i];
    FILE *f = ud->files[i];
    fprintf(f, "%5d", step->coor);
    fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_RG) / n);
    fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_SQRRG) / n);
    fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_ANIS) / n);
    fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_ACYL) / n);
    fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_ASPHER) / n);
    for (int dd = 0; dd < 3; dd++) {
      fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_EIGEN0 + dd) / n);
    }
    double re_n = GetArr2D(ud->step_vals, i, COL_RE_N);
    if (re_n > 0) {
      fprintf(f, " %8.5f", GetArr2D(ud->step_vals, i, COL_RE) / re_n);
    } else {
      fprintf(f, " %8.5f", NAN);
    }
    putc('\n', f);
  }
}
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  Calculation(System, step, (struct user_data *)userdata);
} //}}}

int main(int argc, char *argv[]) {

  // command-line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  struct OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}
  // <output> - output file base name
  char output[LINE];
  s_strcpy(output, argv[++count], LINE);
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option - if present, don't join molecules
  opt.join = !BoolOption(argc, argv, "--joined");
  //}}}

  // print command to stdout
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // -mt option - molecule types to use //{{{
  opt.mt = calloc(Count->MoleculeType, sizeof *opt.mt);
  if (!opt.mt) {
    ErrorAlloc("opt.mt");
  }
  if (!TypeOption(argc, argv, "-mt", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  } //}}}

  // -bt option - bead types to use //{{{
  opt.bt = calloc(Count->BeadType, sizeof *opt.bt);
  if (!opt.bt) {
    ErrorAlloc("opt.bt");
  }
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System)) {
    InitBoolArray(opt.bt, Count->BeadType, true);
  } //}}}

  // warn if any selected molecule type is not a linear chain //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i]) {
      continue;
    }
    MOLECULETYPE *mtype = &System.MoleculeType[i];
    // count bonds per bead
    int *degree = calloc(mtype->nBeads, sizeof *degree);
    if (!degree) {
      ErrorAlloc("degree");
    }
    for (int j = 0; j < mtype->nBonds; j++) {
      degree[mtype->Bond[j][0]]++;
      degree[mtype->Bond[j][1]]++;
    }
    int max_degree = 0;
    for (int j = 0; j < mtype->nBeads; j++) {
      if (degree[j] > max_degree) {
        max_degree = degree[j];
      }
    }
    free(degree);
    if (max_degree > 2) {
      snprintf(ERROR_MSG, LINE, "molecule type %s%s%s has branched topology "
               "(max bead degree: %s%d%s); Re is ill-defined",
               ErrYellow(), mtype->Name, ErrCyan(),
               ErrYellow(), max_degree, ErrCyan());
      PrintWarning();
    } else if (mtype->nBonds != mtype->nBeads - 1) {
      snprintf(ERROR_MSG, LINE, "molecule type %s%s%s is not a linear chain "
               "(nBonds=%s%d%s, nBeads=%s%d%s); Re is ill-defined",
               ErrYellow(), mtype->Name, ErrCyan(),
               ErrYellow(), mtype->nBonds, ErrCyan(),
               ErrYellow(), mtype->nBeads, ErrCyan());
      PrintWarning();
    }
  } //}}}

  // open output files and write headers //{{{
  FILE **file = calloc(Count->MoleculeType, sizeof *file);
  if (!file) {
    ErrorAlloc("files");
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i]) {
      continue;
    }
    char fname[LINE+MOL_NAME+5];
    snprintf(fname, sizeof fname, "%s-%s.txt", output,
             System.MoleculeType[i].Name);
    file[i] = PrintBylineOpenFile(fname, argc, argv);
    fprintf(file[i], "# %s\n", System.MoleculeType[i].Name);
    int col = 1;
    fprintf(file[i], "# column:");
    fprintf(file[i], " (%d) timestep", col++);
    fprintf(file[i], ", (%d) <Rg>", col++);
    fprintf(file[i], ", (%d) <Rg^2>", col++);
    fprintf(file[i], ", (%d) <Anis>", col++);
    fprintf(file[i], ", (%d) <Acyl>", col++);
    fprintf(file[i], ", (%d) <Aspher>", col++);
    fprintf(file[i], ", (%d) <eigen[0]>", col++);
    fprintf(file[i], ", (%d) <eigen[1]>", col++);
    fprintf(file[i], ", (%d) <eigen[2]>", col++);
    fprintf(file[i], ", (%d) <Re>", col++);
    putc('\n', file[i]);
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // allocate sum/step arrays and bead-list buffer //{{{
  ArrNDd *sums = CreateArr2Dd(Count->MoleculeType, N_COLS);
  int *mol_count = calloc(Count->MoleculeType, sizeof *mol_count);
  ArrNDd *step_vals = CreateArr2Dd(Count->MoleculeType, N_COLS);
  int *step_count = calloc(Count->MoleculeType, sizeof *step_count);
  if (!sums || !mol_count || !step_vals || !step_count) {
    ErrorAlloc("sum/step arrays");
  }
  int list_cap = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt.mt[i] && System.MoleculeType[i].nBeads > list_cap) {
      list_cap = System.MoleculeType[i].nBeads;
    }
  }
  int *list = nullptr;
  if (list_cap > 0) {
    if (!(list = malloc(list_cap * sizeof *list))) {
      ErrorAlloc("list");
    }
  } //}}}

  struct user_data ud = {
    .opt = opt,
    .sums = sums,
    .mol_count = mol_count,
    .step_vals = step_vals,
    .step_count = step_count,
    .files = file,
    .list = list,
  };

  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // append overall averages and close files //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i]) {
      continue;
    }
    if (mol_count[i] > 0) {
      int n = mol_count[i];
      FILE *f = file[i];
      fprintf(f, "# overall averages (%d molecules total):\n", n);
      fprintf(f, "# <Rg> <Rg^2> <Anis> <Acyl> <Aspher>"
                 " <eigen[0]> <eigen[1]> <eigen[2]> <Re>\n");
      fprintf(f, "#");
      fprintf(f, " %lf", GetArr2D(sums, i, COL_RG) / n);
      fprintf(f, " %lf", GetArr2D(sums, i, COL_SQRRG) / n);
      fprintf(f, " %lf", GetArr2D(sums, i, COL_ANIS) / n);
      fprintf(f, " %lf", GetArr2D(sums, i, COL_ACYL) / n);
      fprintf(f, " %lf", GetArr2D(sums, i, COL_ASPHER) / n);
      for (int dd = 0; dd < 3; dd++) {
        fprintf(f, " %lf", GetArr2D(sums, i, COL_EIGEN0 + dd) / n);
      }
      double re_n = GetArr2D(sums, i, COL_RE_N);
      if (re_n > 0) {
        fprintf(f, " %lf", GetArr2D(sums, i, COL_RE) / re_n);
      } else {
        fprintf(f, " %lf", NAN);
      }
      putc('\n', f);
    }
    fclose(file[i]);
  } //}}}

  // free memory //{{{
  FreeSystem(&System);
  free(opt.mt);
  free(opt.bt);
  FreeArrND(sums);
  free(mol_count);
  FreeArrND(step_vals);
  free(step_count);
  free(file);
  free(list);
  //}}}

  return 0;
}
