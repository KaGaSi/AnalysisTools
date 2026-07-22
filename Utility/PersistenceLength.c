#include "../src/AnalysisTools.h"
// TODO: very messy!!!
// TODO: explain S1 through S3
// TODO: output handling - only after it's decided what to print!

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

// go through all molecules and calcule l_p & Co. //{{{
// warn once (per run) that coincident beads produced a zero-length bond vector
static void WarnDegenerateBond(void) {
  static bool warned = false;
  if (!warned) {
    err_msg("zero-length bond vector (coincident beads); affected S1/S2 terms "
            "are skipped - note S1 is normalized by the step count, so its "
            "value for the affected chains is biased low");
    PrintWarning();
    warned = true;
  }
}
void Calculation(SYSTEM *System, OPT opt, double *bondlength, int *count_bonds,
                 ArrNDd *S1, ArrNDd *S2, ArrNDd *S3,
                 ArrNDi *count_S2, ArrNDi *count_S3) {
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
      vec3d bond1 = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
      // last bond vector (for reversed S1)
      b1 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][0]],
      b2 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][1]];
      vec3d bondN = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
      for (int k = first_bond; k < last_bond; k++) {
        // S1 function & bondlengths //{{{
        // 1->N S1
        b1 = mol->Bead[mt->Bond[k][0]];
        b2 = mol->Bead[mt->Bond[k][1]];
        vec3d bondj = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
        double s1_fwd = CosAngle(bondj, bond1);
        if (!isnan(s1_fwd)) {
          AddArr3D(S1, mol->Type, k - first_bond, 0, s1_fwd);
        } else { // degenerate (zero-length) bond vector
          WarnDegenerateBond();
        }
        // bondlength & count bonds
        bondlength[mol->Type] += VectLength(bondj);
        count_bonds[mol->Type]++;
        // reverse S1
        int bond_id = mt->nBonds - k - 1;
        int bin_id = k - first_bond;
        b1 = mol->Bead[mt->Bond[bond_id][0]];
        b2 = mol->Bead[mt->Bond[bond_id][1]];
        bondj = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
        double s1_rev = CosAngle(bondN, bondj);
        if (!isnan(s1_rev)) {
          AddArr3D(S1, mol->Type, bin_id, 1, s1_rev);
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
          vec3d bondl = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
          // autocorrelation
          double s2 = CosAngle(bondk, bondl);
          if (!isnan(s2)) {
            AddArr2D(S2, mol->Type, lag, s2);
            AddArr2D(count_S2, mol->Type, lag, 1);
          } else { // degenerate (zero-length) bond vector
            WarnDegenerateBond();
          }
          //}}}
          // S3 function (end-to-end distances) //{{{
          b1 = mol->Bead[mt->Bond[k][0]];
          b2 = mol->Bead[mt->Bond[l][1]];
          vec3d Re = Vector(System->Bead[b1].Position, System->Bead[b2].Position);
          AddArr2D(S3, mol->Type, lag, SqVectLength(Re));
          AddArr2D(count_S3, mol->Type, lag, 1);
          //}}}
        }
      }
    }
  }
} //}}}
// structure for the callback function
struct user_data {
  OPT opt;
  double *bondlength;
  int *count_bonds;
  ArrNDd *S1, *S2, *S3;
  ArrNDi *count_S2, *count_S3;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, p->opt, p->bondlength, p->count_bonds,
              p->S1, p->S2, p->S3, p->count_S2, p->count_S3);
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

  // arrays for the observables
  double *bondlength = calloc(Count->MoleculeType, sizeof *bondlength);
  int *count_bonds = calloc(Count->MoleculeType, sizeof *count_bonds);
  ArrNDd *S1 = CreateArr3Dd(Count->MoleculeType, max_bonds, 2);
  ArrNDd *S2 = CreateArr2Dd(Count->MoleculeType, max_bonds);
  ArrNDi *count_S2 = CreateArr2Di(Count->MoleculeType, max_bonds);
  ArrNDd *S3 = CreateArr2Dd(Count->MoleculeType, max_bonds);
  ArrNDi *count_S3 = CreateArr2Di(Count->MoleculeType, max_bonds);
  if (!S1 || !S2 || !S3 || !count_S2 || !count_S3 ||
      !bondlength || !count_bonds) {
    ErrorAlloc("S1/S2/S3/count_S2/count_S3/bondlength/count_bonds");
  }

  STEP step = InitStep;
  struct user_data ud = { opt, bondlength, count_bonds,
                          S1, S2, S3, count_S2, count_S3 };
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
          double avg = GetArr3D(S1, j, lag, dd) / step.used;
          sum_S1[j][dd] += avg;
          data[lag][++count] = avg;
          data[lag][++count] = sum_S1[j][dd];
        }
        // S2 (autocorrelation function)
        double avg = GetArr2D(S2, j, lag) / GetArr2D(count_S2, j, lag);
        sum_S2[j] += avg;
        data[lag][++count] = avg;
        data[lag][++count] = sum_S2[j];
        // S3 (end-to-end distances)
        data[lag][++count] = GetArr2D(S3, j, lag) / GetArr2D(count_S3, j, lag);
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
  FreeArrND(S1);
  FreeArrND(S2);
  FreeArrND(S3);
  FreeArrND(count_S2);
  FreeArrND(count_S3);
  free(bondlength);
  free(count_bonds);
  FreeSystem(&System);
  //}}}

  return 0;
}

// backup - will the MainLoop() function work?
// // main loop //{{{
// FILE *fr = OpenFile(in.coor.name, "r");
// int count_coor = 0, // count steps in the vcf file
//     count_used = 0, // count steps in output file
//     line_count = 0; // count lines in the vcf file
// while (true) {
//   PrintStep(&count_coor, commons.start, commons.silent);
//   // use every skip-th timestep between start and end
//   bool use = false;
//   if (UseStep(commons, count_coor)) {
//     use = true;
//   }
//   if (use) { //{{{
//     if (!ReadTimestep(in, fr, &System, &line_count)) {
//       count_coor--;
//       break;
//     }
//     count_used++;
//     WrapJoinCoordinates(&System, false, opt.join);
//     // go through all molecules //{{{
//     for (int i = 0; i < Count->MoleculeType; i++) {
//       MOLECULETYPE *mt = &System.MoleculeType[i];
//       // last bond id
//       int last_bond = mt->nBonds;
//       if (opt.ne != HIGHNUM) {
//         last_bond = opt.ne;
//       }
//       int first_bond = 0;
//       if (opt.ns > 0) {
//         first_bond = opt.ns;
//       }
//       // use only specified molecule types that are long enough
//       if (!opt.mt[i] || mt->nBonds < opt.ns) {
//         continue;
//       }
//
//       for (int j = 0; j < mt->Number; j++) {
//         MOLECULE *mol = &System.Molecule[mt->Index[j]];
//         // S1 function
//         // first bond vector (for S1)
//         int b1 = mol->Bead[mt->Bond[first_bond][0]],
//             b2 = mol->Bead[mt->Bond[first_bond][1]];
//         vec3d bond1 = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//         // last bond vector (for reversed S1)
//         b1 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][0]],
//         b2 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][1]];
//         vec3d bondN = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//         for (int k = first_bond; k < last_bond; k++) {
//           // S1 function & bondlengths //{{{
//           // 1->N S1
//           b1 = mol->Bead[mt->Bond[k][0]];
//           b2 = mol->Bead[mt->Bond[k][1]];
//           vec3d bondj = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//           AddArr3D(S1, mol->Type, k - first_bond, 0, CosAngle(bondj, bond1));
//           // bondlength & count bonds
//           bondlength[mol->Type] += VectLength(bondj);
//           count_bonds[mol->Type]++;
//           // reverse S1
//           int bond_id = mt->nBonds - k - 1;
//           int bin_id = k - first_bond;
//           b1 = mol->Bead[mt->Bond[bond_id][0]];
//           b2 = mol->Bead[mt->Bond[bond_id][1]];
//           bondj = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//           AddArr3D(S1, mol->Type, bin_id, 1, CosAngle(bondN, bondj));
//           //}}}
//           for (int l = k; l < last_bond; l++) {
//             int lag = l - k;
//             // S2 function (classic bond correlation) //{{{
//             // first bond vector
//             int b1 = mol->Bead[mt->Bond[k][0]],
//                 b2 = mol->Bead[mt->Bond[k][1]];
//             vec3d bondk = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//             // second bond vector
//             b1 = mol->Bead[mt->Bond[l][0]];
//             b2 = mol->Bead[mt->Bond[l][1]];
//             vec3d bondl = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//             // autocorrelation
//             AddArr2D(S2, mol->Type, lag, CosAngle(bondk, bondl));
//             AddArr2D(count_S2, mol->Type, lag, 1);
//             //}}}
//             // S3 function (end-to-end distances) //{{{
//             b1 = mol->Bead[mt->Bond[k][0]];
//             b2 = mol->Bead[mt->Bond[l][1]];
//             vec3d Re = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
//             AddArr2D(S3, mol->Type, lag, SqVectLength(Re));
//             AddArr2D(count_S3, mol->Type, lag, 1);
//             //}}}
//           }
//         }
//       }
//     } //}}}
//   //}}}
//   } else { //{{{
//     if (!SkipTimestep(in, fr, &line_count)) {
//       count_coor--;
//       break;
//     }
//   } //}}}
//   // exit the main loop if reached user-specied end timestep
//   if (count_coor == commons.end) {
//     break;
//   }
// }
// fclose(fr);
// PrintLastStep(count_coor, count_used, commons.silent); //}}}
