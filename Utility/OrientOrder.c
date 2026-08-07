#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "OrientOrder utility calculates the orientational order parameter "
  "S=0.5*(3*cos^2<angle>-1) for specified bead pairs in specified molecule "
  "type(s), where θ is the angle between the bead-pair vector and the "
  "bilayer normal axis. For each pair it outputs the distribution of S "
  "over [-0.5, 1] and the average value appended at the bottom of the file. "
  "With --bonds it also writes <output>-bonds.txt containing the average S_n "
  "for each consecutive bead-to-bead bond along the chain (bond n = beads "
  "n → n+1), one column per active molecule type.",

  "Usage: OrientOrder <input> <width> <output> [options]",
  .args = 3,
  .all = 17,
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
  {"<input>", nullptr, "input coordinate file",                                                        OPT_ARG},
  {"<width>", nullptr, "width of a distribution bin",                                                  OPT_ARG},
  {"<output>", nullptr, "output file with S distribution and averages",                                 OPT_ARG},
  {"-m", "<name(s)>", "molecule types to use (default: all)",                                         OPT_EXTRA},
  {"--joined", nullptr, "input contains joined coordinates",                                            OPT_EXTRA},
  {"-a", "<axis>", "bilayer normal axis: x, y, or z (default: z)",                                 OPT_EXTRA},
  {"-n", "<int int ...>","bead pairs (must be in pairs; default: first and last bead)",                   OPT_EXTRA},
  {"--bonds", nullptr, "also write per-bond S_n profile to <output>-bonds.txt",                        OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  bool join;
  bool bonds; // write per-bond S_n profile
  bool *mt;
  vec3d normal;
  int n_list[100], // bead pairs, 1-indexed (always even number of elements)
      n_number;    // number of elements in n_list
}; //}}}

// per-molecule-type accumulators //{{{
// a running sum with its matched sample count
typedef struct {
  ArrNDd *sum;
  ArrNDd *count;
} CORR;
// bundled so Calculation() takes one accumulator instead of five arrays
typedef struct {
  ArrNDd *dist; // S histogram [type][pair][bin]
  CORR pair;    // per-pair average S: sum + count [type][pair]
  CORR bond;    // per-bond average S_n: sum + count [type][bond] (--bonds only)
} ACC; //}}}

static void Calculation(SYSTEM *System, OPT opt, ACC *a,
                        double width, int bins, int n_pairs) {
  COUNT *Count = &System->Count;
  WrapJoinCoordinates(System, true, opt.join);

  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol = &System->Molecule[i];
    MOLECULETYPE *mt = &System->MoleculeType[mol->Type];
    if (!opt.mt[mol->Type]) {
      continue;
    }

    // pair-based S distribution (existing behaviour)
    for (int p = 0; p < n_pairs; p++) {
      int pos1, pos2;
      if (opt.n_number == 0) {
        pos1 = 0;
        pos2 = mt->nBeads - 1;
      } else {
        pos1 = opt.n_list[2*p] - 1;
        pos2 = opt.n_list[2*p+1] - 1;
        // clamp out-of-range ids to last bead (matches v3.5 behaviour)
        if (pos1 >= mt->nBeads) {
          pos1 = mt->nBeads - 1;
        }
        if (pos2 >= mt->nBeads) {
          pos2 = mt->nBeads - 1;
        }
      }
      if (pos1 == pos2) {
        continue; // degenerate pair: skip
      }

      vec3d bvec = Vector(System->Bead[mol->Bead[pos1]].Position,
                          System->Bead[mol->Bead[pos2]].Position);
      double len = VectLength(bvec);
      if (len == 0) {
        continue;
      }

      double cos_theta = Dot(bvec, opt.normal) / len;
      double S = 0.5 * (3 * cos_theta * cos_theta - 1);

      AddArr2D(a->pair.sum, mol->Type, p, S);
      AddArr2D(a->pair.count, mol->Type, p, 1);

      int k = (int)((S + 0.5) / width);
      if (k >= bins) {
        k = bins - 1;
      }
      if (k < 0) {
        k = 0;
      }
      AddArr3D(a->dist, mol->Type, p, k, 1);
    }

    // per-bond S_n profile
    if (opt.bonds) {
      for (int b = 0; b < mt->nBeads - 1; b++) {
        vec3d bvec = Vector(System->Bead[mol->Bead[b]].Position,
                            System->Bead[mol->Bead[b+1]].Position);
        double len = VectLength(bvec);
        if (len == 0) {
          continue;
        }
        double cos_theta = Dot(bvec, opt.normal) / len;
        double S = 0.5 * (3 * cos_theta * cos_theta - 1);
        AddArr2D(a->bond.sum, mol->Type, b, S);
        AddArr2D(a->bond.count, mol->Type, b, 1);
      }
    }
  }
}

struct user_data {
  OPT opt;
  ACC acc;
  double width;
  int bins, n_pairs;
};
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data *)userdata;
  Calculation(System, p->opt, &p->acc, p->width, p->bins, p->n_pairs);
}

int main(int argc, char *argv[]) {

  // command line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input>
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <width>
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  // <output>
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // common options
  COMMON_OPT commons = CommonOptions(argc, argv, in);

  // --joined
  opt.join = !BoolOption(argc, argv, "--joined");

  // --bonds
  opt.bonds = BoolOption(argc, argv, "--bonds");

  // -a: bilayer normal axis (default: z)
  opt.normal = (vec3d){.x = 0, .y = 0, .z = 1};
  for (int i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "-a") == 0) {
      switch (argv[i+1][0]) {
        case 'x':
          opt.normal = (vec3d){.x = 1, .y = 0, .z = 0};
          break;
        case 'y':
          opt.normal = (vec3d){.x = 0, .y = 1, .z = 0};
          break;
        case 'z':
          opt.normal = (vec3d){.x = 0, .y = 0, .z = 1};
          break;
        default:
          err_msg("argument must be x, y, or z");
          PrintErrorOption("-a");
          exit(1);
      }
      break;
    }
  }

  // -n: bead pairs (1-indexed, must come in pairs)
  opt.n_number = 0;
  NumbersOption(argc, argv, 100, "-n", &opt.n_number, opt.n_list, 'i');
  if (opt.n_number % 2 != 0) {
    err_msg("bead indices must be given in pairs");
    PrintErrorOption("-n");
    exit(1);
  }
  for (int i = 0; i < opt.n_number; i += 2) {
    if (opt.n_list[i] == opt.n_list[i+1] ||
        opt.n_list[i] <= 0 || opt.n_list[i+1] <= 0) {
      err_msg("each pair must contain two different positive bead ids");
      PrintErrorOption("-n");
      exit(1);
    }
  }
  int n_pairs = opt.n_number / 2;
  if (opt.n_number == 0) {
    n_pairs = 1;
  } //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // S \in [-0.5, 1.0], range = 1.5
  int bins = (int)(1.5 / width);
  if (bins < 1) {
    bins = 1;
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // -m: molecule types to use (default: all)
  if (!(opt.mt = calloc(Count->MoleculeType, sizeof *opt.mt))) {
    ErrorAlloc("opt.mt");
  }
  if (!TypeOption(argc, argv, "-m", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // allocate distribution, average, and count arrays //{{{
  ACC acc = {0};
  acc.dist = CreateArr3Dd(Count->MoleculeType, n_pairs, bins);
  acc.pair.sum = CreateArr2Dd(Count->MoleculeType, n_pairs);
  acc.pair.count = CreateArr2Dd(Count->MoleculeType, n_pairs);
  if (!acc.dist || !acc.pair.sum || !acc.pair.count) {
    ErrorAlloc("dist/avg/count");
  } //}}}

  // allocate per-bond arrays if --bonds //{{{
  int max_bonds = 0;
  if (opt.bonds) {
    for (int mt = 0; mt < Count->MoleculeType; mt++) {
      if (!opt.mt[mt]) {
        continue;
      }
      int nb = System.MoleculeType[mt].nBeads - 1;
      if (nb > max_bonds) {
        max_bonds = nb;
      }
    }
    if (max_bonds > 0) {
      acc.bond.sum   = CreateArr2Dd(Count->MoleculeType, max_bonds);
      acc.bond.count = CreateArr2Dd(Count->MoleculeType, max_bonds);
      if (!acc.bond.sum || !acc.bond.count) {
        ErrorAlloc("avg_bonds/count_bonds");
      }
    }
  } //}}}

  STEP step = InitStep;
  struct user_data ud = { opt, acc, width, bins, n_pairs };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // write S distribution output //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);

  // header //{{{
  fprintf(fw, "# (1) S");
  count = 1;
  for (int mt = 0; mt < Count->MoleculeType; mt++) {
    if (!opt.mt[mt]) {
      continue;
    }
    MOLECULETYPE *mtype = &System.MoleculeType[mt];
    for (int p = 0; p < n_pairs; p++) {
      int lo, hi;
      if (opt.n_number == 0) {
        lo = 1;
        hi = mtype->nBeads;
      } else {
        lo = opt.n_list[2*p];
        hi = opt.n_list[2*p+1];
        if (lo - 1 >= mtype->nBeads) {
          lo = mtype->nBeads;
        }
        if (hi - 1 >= mtype->nBeads) {
          hi = mtype->nBeads;
        }
      }
      if (lo == hi) {
        continue; // degenerate pair for this molecule type
      }
      fprintf(fw, "; (%d) %s:%d-%d", ++count, mtype->Name, lo, hi);
    }
  }
  putc('\n', fw); //}}}

  // distribution //{{{
  int ncols = count;
  int nrows = bins;
  ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
  if (!data) {
    ErrorAlloc("data");
  }
  for (int i = 0; i < nrows; i++) {
    count = 0;
    double s_center = -0.5 + width * (i + 0.5);
    SetArr2D(data, i, count++, s_center);
    for (int mt = 0; mt < Count->MoleculeType; mt++) {
      if (!opt.mt[mt]) {
        continue;
      }
      MOLECULETYPE *mtype = &System.MoleculeType[mt];
      for (int p = 0; p < n_pairs; p++) {
        int lo, hi;
        if (opt.n_number == 0) {
          lo = 1;
          hi = mtype->nBeads;
        } else {
          if ((opt.n_list[2*p] - 1) >= mtype->nBeads) {
            lo = mtype->nBeads;
          } else {
            lo = opt.n_list[2*p];
          }
          if ((opt.n_list[2*p+1] - 1) >= mtype->nBeads) {
            hi = mtype->nBeads;
          } else {
            hi = opt.n_list[2*p+1];
          }
        }
        if (lo == hi) {
          continue;
        }
        double total = GetArr2D(acc.pair.count, mt, p);
        double val = 0;
        if (total > 0) {
          val = GetArr3D(acc.dist, mt, p, i) / total;
        }
        SetArr2D(data, i, count++, val);
      }
    }
  }
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(fw, nrows, ncols, data);
  FreeArrND(data); //}}}

  // averages //{{{
  fprintf(fw, "# Average S:");
  for (int mt = 0; mt < Count->MoleculeType; mt++) {
    if (!opt.mt[mt]) {
      continue;
    }
    MOLECULETYPE *mtype = &System.MoleculeType[mt];
    for (int p = 0; p < n_pairs; p++) {
      int lo, hi;
      if (opt.n_number == 0) {
        lo = 1;
        hi = mtype->nBeads;
      } else {
        if (opt.n_list[2*p] - 1 >= mtype->nBeads) {
          lo = mtype->nBeads;
        } else {
          lo = opt.n_list[2*p];
        }
        if (opt.n_list[2*p+1] - 1 >= mtype->nBeads) {
          hi = mtype->nBeads;
        } else {
          hi = opt.n_list[2*p+1];
        }
      }
      if (lo == hi) {
        continue;
      }
      double total = GetArr2D(acc.pair.count, mt, p);
      double s_avg = 0;
      if (total > 0) {
        s_avg = GetArr2D(acc.pair.sum, mt, p) / total;
      }
      fprintf(fw, " %s:%d-%d=%.6f", mtype->Name, lo, hi, s_avg);
    }
  }
  putc('\n', fw); //}}}

  fclose(fw); //}}}

  // write per-bond S_n profile //{{{
  if (opt.bonds && max_bonds > 0) {
    char fbonds[LINE+16];
    snprintf(fbonds, sizeof fbonds, "%s-bonds.txt", fout);
    fw = PrintBylineOpenFile(fbonds, argc, argv);

    // header: one column per active mol type with at least 2 beads
    fprintf(fw, "# (1) bond");
    int col = 1;
    for (int mt = 0; mt < Count->MoleculeType; mt++) {
      if (!opt.mt[mt] || System.MoleculeType[mt].nBeads < 2) {
        continue;
      }
      fprintf(fw, "; (%d) %s", ++col, System.MoleculeType[mt].Name);
    }
    putc('\n', fw);

    // count active columns for ArrNDd
    int n_active = 0;
    for (int mt = 0; mt < Count->MoleculeType; mt++) {
      if (opt.mt[mt] && System.MoleculeType[mt].nBeads >= 2) {
        n_active++;
      }
    }
    int bncols = 1 + n_active;

    ArrNDd *bdata = CreateArr2Dd(max_bonds + 2, bncols);
    if (!bdata) {
      ErrorAlloc("bdata");
    }

    for (int b = 0; b < max_bonds; b++) {
      int c = 0;
      SetArr2D(bdata, b, c++, (double)(b + 1)); // 1-indexed bond number
      for (int mt = 0; mt < Count->MoleculeType; mt++) {
        if (!opt.mt[mt] || System.MoleculeType[mt].nBeads < 2) {
          continue;
        }
        double s_avg = 0.0;
        if (b < System.MoleculeType[mt].nBeads - 1) {
          double total = GetArr2D(acc.bond.count, mt, b);
          if (total > 0) {
            s_avg = GetArr2D(acc.bond.sum, mt, b) / total;
          }
        }
        SetArr2D(bdata, b, c++, s_avg);
      }
    }
    ComputeColumnWidths(max_bonds, bncols, bdata, 6);
    PrintDataAll(fw, max_bonds, bncols, bdata);
    FreeArrND(bdata);
    fclose(fw);
  } //}}}

  // free memory
  free(opt.mt);
  FreeArrND(acc.dist);
  FreeArrND(acc.pair.sum);
  FreeArrND(acc.pair.count);
  FreeArrND(acc.bond.sum);
  FreeArrND(acc.bond.count);
  FreeSystem(&System);

  return 0;
}
