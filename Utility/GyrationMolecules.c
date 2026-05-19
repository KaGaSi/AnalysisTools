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
  .all = 14,  // number of valid lines in OptSpec (not counting last {NULL})
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
  {"<input>",  NULL,        "input coordinate file",                              OPT_ARG},
  {"<output>", NULL,        "output base name (gets '-<molname>.txt' appended)",  OPT_ARG},
  {"--joined", NULL,        "<input> contains joined coordinates",                OPT_EXTRA},
  {"-mt",      "<name(s)>", "molecule types to use (default: all)",               OPT_EXTRA},
  {"-bt",      "<name(s)>", "bead types used for calculation (default: all)",     OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  bool join,  // --joined
       *mt,   // -mt (per molecule type)
       *bt;   // -bt (per bead type)
}; //}}}

// per-timestep calculation and output //{{{
static void Calculation(SYSTEM *System, STEP *step,
                        const struct OPT opt, const char *output,
                        double *Rg_sum, double *sqrRg_sum,
                        double *Anis_sum, double *Acyl_sum, double *Aspher_sum,
                        double *Re_sum,
                        vec3d *eigen_sum,
                        int *total_mol_count) {
  COUNT *Count = &System->Count;
  // wrap all beads; join only selected molecule types
  WrapJoinCoordinates(System, true, false);
  if (opt.join) {
    for (int i = 0; i < Count->Molecule; i++) {
      if (opt.mt[System->Molecule[i].Type]) {
        RemovePBCMolecule(i, System);
      }
    }
  }

  // per-timestep accumulators
  double *Rg_step = calloc(Count->MoleculeType, sizeof *Rg_step);
  double *sqrRg_step = calloc(Count->MoleculeType, sizeof *sqrRg_step);
  double *Anis_step = calloc(Count->MoleculeType, sizeof *Anis_step);
  double *Acyl_step = calloc(Count->MoleculeType, sizeof *Acyl_step);
  double *Aspher_step = calloc(Count->MoleculeType, sizeof *Aspher_step);
  double *Re_step = calloc(Count->MoleculeType, sizeof *Re_step);
  vec3d *eigen_step = calloc(Count->MoleculeType, sizeof *eigen_step);
  int *mol_count_step = calloc(Count->MoleculeType, sizeof *mol_count_step);
  if (!Rg_step || !sqrRg_step || !Anis_step || !Acyl_step || !Aspher_step ||
      !Re_step || !eigen_step || !mol_count_step) {
    ErrorAlloc("step arrays");
  }

  // calculate shape descriptors for each molecule
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol = &System->Molecule[i];
    MOLECULETYPE *mtype = &System->MoleculeType[mol->Type];
    if (!opt.mt[mol->Type]) {
      continue;
    }

    // build list of bead ids filtered by -bt
    int *list = malloc(mtype->nBeads * sizeof *list);
    if (!list) {
      ErrorAlloc("list");
    }
    int n = 0;
    for (int j = 0; j < mtype->nBeads; j++) {
      int id = mol->Bead[j];
      if (opt.bt[System->Bead[id].Type]) {
        list[n++] = id;
      }
    }
    if (n < 2) { // need at least 2 beads for gyration
      free(list);
      continue;
    }

    vec3d eigen = Gyration(n, list, System);
    free(list);
    if (eigen.x == 0 && eigen.y == 0 && eigen.z == 0) {
      continue;
    }

    int mt = mol->Type;
    double Rgi = sqrt(eigen.x + eigen.y + eigen.z);
    Rg_step[mt] += Rgi;
    sqrRg_step[mt] += Square(Rgi);
    Anis_step[mt] += 1.5 * SqVectLength(eigen) /
                       Square(eigen.x + eigen.y + eigen.z) - 0.5;
    Acyl_step[mt] += eigen.y - eigen.x;
    Aspher_step[mt] += eigen.z - 0.5 * (eigen.x + eigen.y);
    for (int dd = 0; dd < 3; dd++) {
      eigen_step[mt].v[dd] += eigen.v[dd];
    }
    // end-to-end distance (always uses first and last bead of molecule)
    int first = mol->Bead[0];
    int last  = mol->Bead[mtype->nBeads - 1];
    vec3d dist = Vector(System->Bead[last].Position, System->Bead[first].Position);
    Re_step[mt] += VectLength(dist);
    mol_count_step[mt]++;
  }

  // add to overall sums and write per-timestep data to output files
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i] || mol_count_step[i] == 0) {
      continue;
    }
    Rg_sum[i] += Rg_step[i];
    sqrRg_sum[i] += sqrRg_step[i];
    Anis_sum[i] += Anis_step[i];
    Acyl_sum[i] += Acyl_step[i];
    Aspher_sum[i] += Aspher_step[i];
    Re_sum[i] += Re_step[i];
    for (int dd = 0; dd < 3; dd++) {
      eigen_sum[i].v[dd] += eigen_step[i].v[dd];
    }
    total_mol_count[i] += mol_count_step[i];

    char fname[LINE + MOL_NAME + 5];
    snprintf(fname, sizeof fname, "%s-%s.txt", output,
             System->MoleculeType[i].Name);
    FILE *out = OpenFile(fname, "a");
    int n = mol_count_step[i];
    fprintf(out, "%5d", step->coor);
    fprintf(out, " %8.5f", Rg_step[i] / n);
    fprintf(out, " %8.5f", sqrRg_step[i] / n);
    fprintf(out, " %8.5f", Anis_step[i] / n);
    fprintf(out, " %8.5f", Acyl_step[i] / n);
    fprintf(out, " %8.5f", Aspher_step[i] / n);
    for (int dd = 0; dd < 3; dd++) {
      fprintf(out, " %8.5f", eigen_step[i].v[dd] / n);
    }
    fprintf(out, " %8.5f", Re_step[i] / n);
    putc('\n', out);
    fclose(out);
  }

  free(Rg_step);
  free(sqrRg_step);
  free(Anis_step);
  free(Acyl_step);
  free(Aspher_step);
  free(Re_step);
  free(eigen_step);
  free(mol_count_step);
} //}}}

// userdata struct for the callback
struct user_data {
  struct OPT opt;
  char output[LINE];
  double *Rg_sum, *sqrRg_sum, *Anis_sum, *Acyl_sum, *Aspher_sum, *Re_sum;
  vec3d *eigen_sum;
  int *total_mol_count;
};

static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data *)userdata;
  Calculation(System, step, p->opt, p->output,
              p->Rg_sum, p->sqrRg_sum, p->Anis_sum, p->Acyl_sum, p->Aspher_sum,
              p->Re_sum, p->eigen_sum,
              p->total_mol_count);
}

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
  // --joined option (opt.join == true -> needs joining)
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

  // write headers to per-molecule-type output files //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i]) {
      continue;
    }
    char fname[LINE + MOL_NAME + 5];
    snprintf(fname, sizeof fname, "%s-%s.txt", output, System.MoleculeType[i].Name);
    FILE *out = PrintBylineOpenFile(fname, argc, argv);
    fprintf(out, "# %s\n", System.MoleculeType[i].Name);
    int col = 1;
    fprintf(out, "# column:");
    fprintf(out, " (%d) timestep", col++);
    fprintf(out, ", (%d) <Rg>",       col++);
    fprintf(out, ", (%d) <Rg^2>",     col++);
    fprintf(out, ", (%d) <Anis>",     col++);
    fprintf(out, ", (%d) <Acyl>",     col++);
    fprintf(out, ", (%d) <Aspher>",   col++);
    fprintf(out, ", (%d) <eigen[0]>", col++);
    fprintf(out, ", (%d) <eigen[1]>", col++);
    fprintf(out, ", (%d) <eigen[2]>", col++);
    fprintf(out, ", (%d) <Re>",       col++);
    putc('\n', out);
    fclose(out);
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // allocate overall sum arrays //{{{
  double *Rg_sum       = calloc(Count->MoleculeType, sizeof *Rg_sum);
  double *sqrRg_sum    = calloc(Count->MoleculeType, sizeof *sqrRg_sum);
  double *Anis_sum     = calloc(Count->MoleculeType, sizeof *Anis_sum);
  double *Acyl_sum     = calloc(Count->MoleculeType, sizeof *Acyl_sum);
  double *Aspher_sum   = calloc(Count->MoleculeType, sizeof *Aspher_sum);
  double *Re_sum       = calloc(Count->MoleculeType, sizeof *Re_sum);
  vec3d  *eigen_sum    = calloc(Count->MoleculeType, sizeof *eigen_sum);
  int    *total_mol_count = calloc(Count->MoleculeType, sizeof *total_mol_count);
  if (!Rg_sum || !sqrRg_sum || !Anis_sum || !Acyl_sum || !Aspher_sum ||
      !Re_sum || !eigen_sum || !total_mol_count) {
    ErrorAlloc("sum arrays");
  } //}}}

  struct user_data ud = {
    .opt = opt,
    .Rg_sum = Rg_sum,
    .sqrRg_sum = sqrRg_sum,
    .Anis_sum = Anis_sum,
    .Acyl_sum = Acyl_sum,
    .Aspher_sum = Aspher_sum,
    .Re_sum = Re_sum,
    .eigen_sum = eigen_sum,
    .total_mol_count = total_mol_count,
  };
  s_strcpy(ud.output, output, LINE);

  STEP step = InitStep;
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // append overall averages to output files //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!opt.mt[i] || total_mol_count[i] == 0) {
      continue;
    }
    char fname[LINE + MOL_NAME + 5];
    snprintf(fname, sizeof fname, "%s-%s.txt", output, System.MoleculeType[i].Name);
    FILE *out = OpenFile(fname, "a");
    int n = total_mol_count[i];
    fprintf(out, "# overall averages (%d molecules total):\n", n);
    fprintf(out, "# <Rg> <Rg^2> <Anis> <Acyl> <Aspher>"
                 " <eigen[0]> <eigen[1]> <eigen[2]> <Re>\n");
    fprintf(out, "#");
    fprintf(out, " %lf", Rg_sum[i] / n);
    fprintf(out, " %lf", sqrRg_sum[i] / n);
    fprintf(out, " %lf", Anis_sum[i] / n);
    fprintf(out, " %lf", Acyl_sum[i] / n);
    fprintf(out, " %lf", Aspher_sum[i] / n);
    for (int dd = 0; dd < 3; dd++) {
      fprintf(out, " %lf", eigen_sum[i].v[dd] / n);
    }
    fprintf(out, " %lf", Re_sum[i] / n);
    putc('\n', out);
    fclose(out);
  } //}}}

  // free memory //{{{
  FreeSystem(&System);
  free(opt.mt);
  free(opt.bt);
  free(Rg_sum);
  free(sqrRg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(Re_sum);
  free(eigen_sum);
  free(total_mol_count);
  //}}}

  return 0;
}
