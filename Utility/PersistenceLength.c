#include "../src/AnalysisTools.h"
// TODO: very messy!!!
// TODO: output handling - only after it's decided what to print!

// Observables as a function of bond separation s, per selected molecule type,
// averaged over molecules and timesteps; fit these to obtain l_p:
//   S1 - <cos> of each bond with the first bond (and, reversed, the last bond);
//        the running "sum S1" x <bond length> is Flory's persistence length
//   S2 - lag-averaged bond orientational autocorrelation <cos th(s)> ~ e^-s/lp
//   S3 - mean-square internal distance <R^2(s)> of a subchain of s+1 bonds

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "PersistenceLength calculates correlation of bond vectors and its standard "
  "deviation. To get the peristence length, the data must be fitted via, "
  "typically, an exponential function. Note the utility expects a linear chain "
  "with ordered beads ids (e.g., for a 4-bead chain, the order must be 1-2-3-4, "
  "leading to ordered bonds 1-2, 2-3, 3-4; connectivity like 1-4-2-3 with "
  "bonds 1-4, 2-4, 2-3 could give unexpected results).",

  "Usage: %s <input> <output> [options]",
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
  {"<output>", nullptr, "output file with the persistence length", OPT_ARG},
  {"-m", "<name(s)>", "molecule types to calculate bond lengths for "
    "(default: all molecule types)", OPT_EXTRA},
  {"--joined", nullptr, "specify that <input> contains joined coordinates",
    OPT_EXTRA},
  {"-ns", "<int>", "start with <int>-th bead in a molecule", OPT_EXTRA},
  {"-ne", "<int>", "end with <int>-th bead in a molecule", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  bool join,  // --joined
       *mt;   // -m
  int ns, ne; // -ns/-ne; first bead and bond and last bead and bond
}; //}}}

// per-molecule-type accumulators //{{{
// a correlation observable: running sum with its matched sample count
typedef struct {
  ArrNDd *sum;
  ArrNDi *count;
} CORR;
// bundled so Calculation() takes one accumulator instead of eight arrays
typedef struct {
  CORR s1, s2, s3;    // S1 (forward + reverse), S2, S3
  double *bondlength; // running sum of bond lengths (with count_bonds -> mean)
  int *count_bonds;
} ACC; //}}}

// go through all molecules and calcule l_p & Co. //{{{
// warn once (per run) that coincident beads produced a zero-length bond vector
static void WarnDegenerateBond(void) {
  static bool warned = false;
  if (!warned) {
    err_msg("zero-length bond vector (coincident beads); affected S1/S2 terms "
            "are skipped (the average is therefore taken over fewer samples)");
    PrintWarning();
    warned = true;
  }
}
void Calculation(SYSTEM *System, OPT opt, ACC *a) {
  WrapJoinCoordinates(System, false, opt.join);
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    // last bond id
    int last_bond = mt->nBonds;
    if (opt.ne != HIGHNUM) {
      last_bond = opt.ne;
    }
    int first_bond = 0;
    if (opt.ns > 0) {
      first_bond = opt.ns;
    }
    // use only specified molecule types that are long enough
    if (!opt.mt[i] || mt->nBonds < opt.ns) {
      continue;
    }

    for (int j = 0; j < mt->Number; j++) {
      MOLECULE *mol = &System->Molecule[mt->Index[j]];
      // S1 function
      // first bond vector (for S1)
      int b1 = mol->Bead[mt->Bond[first_bond][0]],
          b2 = mol->Bead[mt->Bond[first_bond][1]];
      vec3d bond1 = Vector(System->Bead[b1].Position,
                           System->Bead[b2].Position);
      // last bond vector (for reversed S1)
      b1 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][0]],
      b2 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][1]];
      vec3d bondN = Vector(System->Bead[b1].Position,
                           System->Bead[b2].Position);
      for (int k = first_bond; k < last_bond; k++) {
        // S1 function & bondlengths //{{{
        // 1->N S1
        b1 = mol->Bead[mt->Bond[k][0]];
        b2 = mol->Bead[mt->Bond[k][1]];
        vec3d bondj = Vector(System->Bead[b1].Position,
                             System->Bead[b2].Position);
        double s1_fwd = CosAngle(bondj, bond1);
        if (!isnan(s1_fwd)) {
          AddArr3D(a->s1.sum, mol->Type, k - first_bond, 0, s1_fwd);
          AddArr3D(a->s1.count, mol->Type, k - first_bond, 0, 1);
        } else { // degenerate (zero-length) bond vector
          WarnDegenerateBond();
        }
        // bondlength & count bonds
        a->bondlength[mol->Type] += VectLength(bondj);
        a->count_bonds[mol->Type]++;
        // reverse S1
        int bond_id = mt->nBonds - k - 1;
        int bin_id = k - first_bond;
        b1 = mol->Bead[mt->Bond[bond_id][0]];
        b2 = mol->Bead[mt->Bond[bond_id][1]];
        bondj = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
        double s1_rev = CosAngle(bondN, bondj);
        if (!isnan(s1_rev)) {
          AddArr3D(a->s1.sum, mol->Type, bin_id, 1, s1_rev);
          AddArr3D(a->s1.count, mol->Type, bin_id, 1, 1);
        } else { // degenerate (zero-length) bond vector
          WarnDegenerateBond();
        }
        //}}}
        for (int l = k; l < last_bond; l++) {
          int lag = l - k;
          // S2 function (classic bond correlation) //{{{
          // first bond vector
          int b1 = mol->Bead[mt->Bond[k][0]],
              b2 = mol->Bead[mt->Bond[k][1]];
          vec3d bondk = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
          // second bond vector
          b1 = mol->Bead[mt->Bond[l][0]];
          b2 = mol->Bead[mt->Bond[l][1]];
          vec3d bondl = Vector(System->Bead[b1].Position,
                               System->Bead[b2].Position);
          // autocorrelation
          double s2 = CosAngle(bondk, bondl);
          if (!isnan(s2)) {
            AddArr2D(a->s2.sum, mol->Type, lag, s2);
            AddArr2D(a->s2.count, mol->Type, lag, 1);
          } else { // degenerate (zero-length) bond vector
            WarnDegenerateBond();
          }
          //}}}
          // S3 function (end-to-end distances) //{{{
          b1 = mol->Bead[mt->Bond[k][0]];
          b2 = mol->Bead[mt->Bond[l][1]];
          vec3d Re = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
          AddArr2D(a->s3.sum, mol->Type, lag, SqVectLength(Re));
          AddArr2D(a->s3.count, mol->Type, lag, 1);
          //}}}
        }
      }
    }
  }
} //}}}
// structure for the callback function
struct user_data {
  OPT opt;
  ACC acc;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, p->opt, &p->acc);
};

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
  // <output> - file name with persistence lengths
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option (opt.join == true -> needs joining)
  opt.join = !BoolOption(argc, argv, "--joined");
  if (!OneNumberOption(argc, argv, "-ns", &opt.ns, 'i')) {
    opt.ns = 1;
  }
  opt.ns--; // indexing starts from 0
  if (!OneNumberOption(argc, argv, "-ne", &opt.ne, 'i')) {
    opt.ne = HIGHNUM;
  } else {
    opt.ne--; // indexing starts from 0
  }
  if (opt.ns != HIGHNUM && opt.ne != HIGHNUM &&
      (opt.ne - opt.ns) < 2) {
    err_msg("at least three beads are necessary, i.e., <-ns> - <-ne> > 1; "
            "(note that the calculation is meaningful only for longer chains)");
    PrintErrorOption("-ns/-ne");
    exit(1);
  }
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-m <name(s)>' option
  if (!(opt.mt = calloc(Count->MoleculeType, sizeof *opt.mt))) {
    ErrorAlloc("opt.mt");
  }
  if (!TypeOption(argc, argv, "-m", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // maximum number of bonds & beads //{{{
  int max_bonds = 0;
  if (opt.ne == HIGHNUM) { // no -ne option -> find longest molecule
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt.mt[i] && System.MoleculeType[i].nBonds > max_bonds) {
        max_bonds = System.MoleculeType[i].nBonds;
      }
    }
  } else { // -ne option -> cannot be longer than the specified length
    max_bonds = opt.ne;
  }
  if (max_bonds < opt.ns) {
    err_msg("starting bead is larger than the length of any molecule");
    ErrorOption("-ns");
    exit(1);
  } //}}}

  // accumulators for the observables (bundled into one struct)
  ACC acc = {0};
  acc.bondlength = calloc(Count->MoleculeType, sizeof *acc.bondlength);
  acc.count_bonds = calloc(Count->MoleculeType, sizeof *acc.count_bonds);
  acc.s1.sum   = CreateArr3Dd(Count->MoleculeType, max_bonds, 2);
  acc.s1.count = CreateArr3Di(Count->MoleculeType, max_bonds, 2);
  acc.s2.sum   = CreateArr2Dd(Count->MoleculeType, max_bonds);
  acc.s2.count = CreateArr2Di(Count->MoleculeType, max_bonds);
  acc.s3.sum   = CreateArr2Dd(Count->MoleculeType, max_bonds);
  acc.s3.count = CreateArr2Di(Count->MoleculeType, max_bonds);
  if (!acc.s1.sum || !acc.s2.sum || !acc.s3.sum ||
      !acc.s1.count || !acc.s2.count || !acc.s3.count ||
      !acc.bondlength || !acc.count_bonds) {
    ErrorAlloc("persistence-length accumulators");
  }

  STEP step = InitStep;
  struct user_data ud = { opt, acc };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // write to output file
  // determine width of each column & collate data //{{{
  int datalines = max_bonds - opt.ns;
  // count used molecule types
  count = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt.mt[i]) {
      count++;
    }
  }
  int data_per_mtype = 7;
  int columns = count * data_per_mtype + 1;
  int digits[columns][2];
  InitInt2DArray((int *)digits, columns, 2, 0);
  double *data[datalines];
  // arrays for integrated functions
  double sum_S1[Count->MoleculeType][2]; // [0] ... 1->N; [1] ... reverse
  double sum_S2[Count->MoleculeType];
  // average bond length
  for (int i = 0; i < Count->MoleculeType; i++) {
    sum_S1[i][0] = 0;
    sum_S1[i][1] = 0;
    sum_S2[i] = 0;
  }
  for (int lag = 0; lag < datalines; lag++) {
    data[lag] = calloc(columns, sizeof *data[lag]);
    if (!data[lag]) {
      ErrorAlloc("data[lag]");
    }
    count = -1;
    // bond lag for x-axis
    data[lag][++count] = lag;
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (opt.mt[j]) {
        // S1 function (from either end)
        for (int dd = 0; dd < 2; dd++) {
          double avg = GetArr3D(acc.s1.sum, j, lag, dd) /
                       GetArr3D(acc.s1.count, j, lag, dd);
          sum_S1[j][dd] += avg;
          data[lag][++count] = avg;
          data[lag][++count] = sum_S1[j][dd];
        }
        // S2 (autocorrelation function)
        double avg = GetArr2D(acc.s2.sum, j, lag) / GetArr2D(acc.s2.count, j, lag);
        sum_S2[j] += avg;
        data[lag][++count] = avg;
        data[lag][++count] = sum_S2[j];
        // S3 (end-to-end distances)
        data[lag][++count] = GetArr2D(acc.s3.sum, j, lag) /
                             GetArr2D(acc.s3.count, j, lag);
      }
    }
  }
  FillMaxDigits(columns, datalines, data, digits); //}}}
  // print the data //{{{
  // headers
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# for each molecule type: ");
  count = 1;
  fprintf(fw, "(%d) S1, ", count++);
  fprintf(fw, "(%d) sum S1, ", count++);
  fprintf(fw, "(%d) S1 (rev), ", count++);
  fprintf(fw, "(%d) sum S1 (rev), ", count++);
  fprintf(fw, "(%d) S2 (autocorr), ", count++);
  fprintf(fw, "(%d) sum S2 (autocorr), ", count++);
  fprintf(fw, "(%d) S3", count++);
  putc('\n', fw);
  fprintf(fw, "# ");
  count = 1;
  fprintf(fw, "(%d) lag; ", count++);
  fprintf(fw, "molecule types: ");
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (opt.mt[i]) {
      fprintf(fw, "(%d)-(%d) %s", count, count+data_per_mtype-1, mt->Name);
      count += data_per_mtype;
      if (i != (Count->MoleculeType - 1)) {
        fprintf(fw, ", ");
      }
    }
  }
  putc('\n', fw);
  // datalines
  for (int lag = 0; lag < datalines; lag++) {
    FprintfRow(fw, columns, data[lag], digits);
    free(data[lag]);
  } //}}}
  fclose(fw);

  // free memory - to make valgrind happy //{{{
  free(opt.mt);
  FreeArrND(acc.s1.sum);
  FreeArrND(acc.s2.sum);
  FreeArrND(acc.s3.sum);
  FreeArrND(acc.s1.count);
  FreeArrND(acc.s2.count);
  FreeArrND(acc.s3.count);
  free(acc.bondlength);
  free(acc.count_bonds);
  FreeSystem(&System);
  //}}}

  return 0;
}
